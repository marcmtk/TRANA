#!/usr/bin/env python
"""Put every read of one sample in a bucket: the stage that dropped it, what emu made of it, and,
with --structure, how the read lies on 16S.

Fate: barbell's output for the sample is followed through host removal, the length and quality
filter, and emu. The filter rule is recomputed per read (length, and mean per-base accuracy as
filtlong averages it) and compared with what filtlong kept, so a change in filtlong's behaviour
shows up as a mismatch instead of silently moving reads between buckets.

Structure: every read is aligned to a few 16S references of each species emu called in the sample
(the references most of that species' reads hit). The walk of the read along the gene is named:
one copy, full length or truncated, a fold-back onto the other strand, a same-strand repeat, a
strand switch without overlap, or a chimera, where references of two species with one breakpoint
explain the read better than any one species does. Positions are on E. coli 16S (K-12 rrsH),
lifted through each reference's own alignment to it.

The chimera thresholds were calibrated on one run (2026-07-02): pure cultures as negatives and
reads of two cultures spliced together as positives. See docs/output.md.

Output:
  <prefix>.read_buckets.tsv.gz  one row per read
  <prefix>.read_buckets.json    sample-level facts for the run report
  <prefix>.read_buckets_refs.fasta  the references the structure step used
"""

import argparse
import collections
import csv
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

import numpy as np

VERSION = "0.1.0"

MIN_SEG = 100  # aligned query bases for a segment to count
MIN_ARM = 100  # bases each reference must explain on its side of a breakpoint
GRID = 5  # breakpoints tried every this many bases
GENE_START, GENE_END = 28, 1492  # the amplicon between the 16S_mix primers, E. coli positions
FULL_MARGIN = 100  # a read spanning the amplicon to within this on both ends is full length
CS = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[+\-][a-zA-Z]+)")
ACC = (1 - 10 ** (-(np.arange(256) - 33) / 10.0)) * 100
ACC[:33] = 0

COLUMNS = [
    "read_id", "length", "accuracy", "deacon", "filter", "filtlong", "emu", "emu_tax_id", "emu_p",
    "struct", "strand", "fstart", "fend", "head", "tail", "arm1", "arm2", "turn1", "turn2",
    "organism", "identity", "s_single", "s_intra", "s_inter", "gain", "bp_frame",
    "org_left", "org_right", "sup_left", "sup_right", "bucket",
]
STRUCT_COLUMNS = COLUMNS[9:30]
BUCKETS = ["full", "truncated", "tail", "colinear", "foldback", "repeat", "switch", "chimera", "poorfit", "unaligned"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--reads", required=True, help="the sample's reads as barbell (or the input) gave them")
    p.add_argument("--dehosted", default="", help="reads left after host removal")
    p.add_argument("--filtered", default="", help="reads filtlong kept")
    p.add_argument("--emu-reads", default="", help="reads emu got, when they differ from --filtered (downsampling)")
    p.add_argument("--abundance", default="", help="emu's rel-abundance table")
    p.add_argument("--assignments", default="", help="emu's read-assignment-distributions table")
    p.add_argument("--db", default="", help="emu database folder (species_taxid.fasta, taxonomy.tsv)")
    p.add_argument("--frame", default="", help="E. coli 16S reference that sets the coordinates")
    p.add_argument("--asvs", default="", help="savont's final ASVs for the sample, used as references too")
    p.add_argument("--asv-mappings", default="", help="savont's asv_mappings.tsv, which names each ASV's species")
    p.add_argument("--structure", action="store_true", help="align reads and name their structure")
    p.add_argument("--min-length", type=int, default=1200)
    p.add_argument("--max-length", type=int, default=1800)
    p.add_argument("--min-accuracy", type=float, default=94)
    p.add_argument("--min-gain", type=int, default=15, help="chimera: score gain of two species over the best one")
    p.add_argument("--min-side", type=int, default=5, help="chimera: score gain on each side of the breakpoint")
    p.add_argument("--min-fit", type=float, default=0.85, help="score per base below which no reference explains a read")
    p.add_argument("--refs-per-species", type=int, default=3, help="database references per called species; one if an ASV covers the species")
    p.add_argument("--max-species", type=int, default=50)
    p.add_argument("--select-reads", type=int, default=100, help="reads per species used to choose its references")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--prefix", default="")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p.parse_args()


# ---------------------------------------------------------------- reads

def fastq(path):
    """(id, sequence, quality) of each record; the id stops at the first whitespace."""
    with gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb") as h:
        while True:
            head = h.readline()
            if not head:
                return
            seq = h.readline().rstrip(b"\r\n")
            h.readline()
            qual = h.readline().rstrip(b"\r\n")
            yield head[1:].split()[0].decode(), seq, qual


def fastq_ids(path):
    if not path:
        return None
    return {rid for rid, _, _ in fastq(path)}


def mean_accuracy(qual):
    return float(ACC[np.frombuffer(qual, dtype=np.uint8)].mean()) if qual else 0.0


# ---------------------------------------------------------------- emu

def read_emu(abundance, assignments):
    calls = []
    if abundance and os.path.getsize(abundance):
        with open(abundance) as h:
            head = h.readline().rstrip("\n").split("\t")
            for line in h:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                if r.get("tax_id", "").isdigit():
                    calls.append(dict(tax_id=r["tax_id"], species=r.get("species") or r["tax_id"],
                                      counts=float(r.get("estimated counts") or 0), abundance=float(r.get("abundance") or 0)))
    calls.sort(key=lambda c: -c["counts"])
    assign = {}
    if assignments and os.path.getsize(assignments):
        with open(assignments) as h:
            taxa = h.readline().rstrip("\n").split("\t")[1:]
            for line in h:
                f = line.rstrip("\n").split("\t")
                probs = [(taxa[i], float(v)) for i, v in enumerate(f[1:]) if v not in ("", "NA") and float(v) > 0]
                if probs:
                    assign[f[0]] = probs
    return calls, assign


# ---------------------------------------------------------------- references

def minimap2(args, capture=True):
    return subprocess.run(["minimap2"] + args, capture_output=capture, text=True, check=True).stdout


def read_fasta(path):
    name, out = None, {}
    with open(path) as h:
        for line in h:
            if line[0] == ">":
                name = line[1:].split()[0]
                out[name] = []
            elif name:
                out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtRYKMBVDHNrykmbvdhn", "TGCAtgcaYRMKVBHDNyrmkvbhdn"))[::-1]


def choose_references(calls, assign, db, reads_path, tmp, a, asvs=()):
    """A few references per called species: those that most of the species' own reads hit best.

    The sample's ASVs go in as references too. An ASV that savont names after a species emu did not
    call is renamed after the closest species emu called, if one is within 1% of it, so that a strain's
    own copy and its species' database references count as one species; among equally close species,
    the one with most reads wins. A species with an ASV gets one database reference, others up to
    --refs-per-species.
    """
    species = calls[: a.max_species]
    wanted = {c["tax_id"] for c in species}
    cand = collections.defaultdict(dict)  # tax_id -> sequence -> database id; identical copies once

    def keep(name, seq):
        if name and seq:
            cand[name.split(":")[0]].setdefault("".join(seq).upper(), name)

    with open(os.path.join(db, "species_taxid.fasta")) as h:
        name, seq = None, []
        for line in h:
            if line[0] == ">":
                keep(name, seq)
                name, seq = line[1:].split()[0], []
                if name.split(":")[0] not in wanted:
                    name = None
            elif name:
                seq.append(line.strip())
        keep(name, seq)
    cfa = os.path.join(tmp, "candidates.fa")
    with open(cfa, "w") as o:
        for t, seqs in cand.items():
            for s, n in seqs.items():
                o.write(f">{n}\n{s}\n")
    called = {c["species"] for c in species}
    rename = [r for r in asvs if r["species"] not in called]
    if rename:
        afa = os.path.join(tmp, "asvs.fa")
        with open(afa, "w") as o:
            for i, r in enumerate(rename):
                o.write(f">{i}\n{r['seq']}\n")
        name_of = {c["tax_id"]: c["species"] for c in species}
        rank = {c["tax_id"]: i for i, c in enumerate(species)}
        best = {}
        for line in minimap2(["-x", "map-ont", "-c", "-N", "50", "-t", str(a.threads), cfa, afa]).splitlines():
            f = line.split("\t")
            if not int(f[10]) or int(f[9]) / int(f[10]) < 0.99 or (int(f[3]) - int(f[2])) < 0.9 * int(f[1]):
                continue
            t = f[5].split(":")[0]
            key = (-int(f[9]), rank[t])
            if f[0] not in best or key < best[f[0]][0]:
                best[f[0]] = (key, t)
        for i, (key, t) in best.items():
            rename[int(i)].update(tax_id=t, species=name_of[t], named_by="closest called species")
    covered = {r["species"] for r in asvs}
    # reads to choose with: those emu put on the species with P >= 0.5
    per = collections.Counter()
    pick = {}
    for rid, probs in assign.items():
        t, p = max(probs, key=lambda x: x[1])
        if p >= 0.5 and t in wanted and per[t] < a.select_reads:
            per[t] += 1
            pick[rid] = t
    hits = collections.Counter()
    if pick:
        sel = os.path.join(tmp, "select.fq")
        with open(sel, "wb") as o:
            for rid, s, q in fastq(reads_path):
                if rid in pick:
                    o.write(b"@" + rid.encode() + b"\n" + s + b"\n+\n" + q + b"\n")
        paf = minimap2(["-x", "map-ont", "-c", "--secondary=no", "-t", str(a.threads), cfa, sel])
        for line in paf.splitlines():
            f = line.split("\t")
            if f[5].split(":")[0] == pick.get(f[0]):
                hits[f[5]] += 1
    refs = []
    for c in species:
        t = c["tax_id"]
        names = {n: s for s, n in cand.get(t, {}).items()}
        ranked = sorted(names, key=lambda n: (-hits[n], n))[: 1 if c["species"] in covered else a.refs_per_species]
        for n in ranked:
            refs.append(dict(db_id=n, tax_id=t, species=c["species"], hits=hits[n], seq=names[n], source="database"))
    return refs + list(asvs)


def asv_references(fasta, mappings, calls):
    """The sample's own ASVs, named after the species of their closest database references.

    Where species tie, one emu called wins, the one with most reads first. These are the strain's own
    16S copies, which matter when two species in a mix differ by a few bases.
    """
    if not (fasta and os.path.exists(fasta) and os.path.getsize(fasta)):
        return []
    rank = {c["species"]: i for i, c in enumerate(calls)}
    best = {}
    if mappings and os.path.exists(mappings):
        with open(mappings) as h:
            for r in csv.DictReader(h, delimiter="\t"):
                try:
                    mm = int(r["number_mismatches"])
                except (KeyError, ValueError):
                    continue
                k = r["asv_header"]
                if k not in best or mm < best[k][0]:
                    best[k] = (mm, {r["species"]: r["tax_id"]})
                elif mm == best[k][0]:
                    best[k][1][r["species"]] = r["tax_id"]
    refs = []
    for name, seq in read_fasta(fasta).items():
        if name not in best:
            continue
        mm, species = best[name]
        sp = min(species, key=lambda x: (rank.get(x, len(rank)), x))
        depth = re.search(r"_depth_(\d+)", name)
        refs.append(dict(db_id=name, tax_id=species[sp], species=sp, hits=int(depth.group(1)) if depth else 0, seq=seq.upper(), source="asv", mismatches=mm, named_by="savont"))
    return refs


def orient_and_frame(refs, frame_path, tmp):
    """Turn each reference to the gene's sense strand and map its positions onto E. coli 16S."""
    rfa = os.path.join(tmp, "refs_raw.fa")
    with open(rfa, "w") as o:
        for i, r in enumerate(refs):
            o.write(f">R{i}\n{r['seq']}\n")
    paf = minimap2(["-x", "map-ont", "-k9", "-w1", "-c", "--secondary=no", "-t", "2", frame_path, rfa])
    best = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if f[0] not in best or int(f[9]) > best[f[0]][1]:
            best[f[0]] = (f[4], int(f[9]))
    for i, r in enumerate(refs):
        r["name"] = f"R{i}"
        r["frame_strand"] = best.get(r["name"], (".", 0))[0]
        if r["frame_strand"] == "-":
            r["seq"] = revcomp(r["seq"])
    ffa = os.path.join(tmp, "refs.fa")
    with open(ffa, "w") as o:
        for r in refs:
            o.write(f">{r['name']}\n{r['seq']}\n")
    paf = minimap2(["-x", "map-ont", "-k9", "-w1", "-c", "--cs=short", "--secondary=no", "-t", "2", frame_path, ffa])
    frame = {r["name"]: np.full(len(r["seq"]) + 1, -1, dtype=np.int32) for r in refs}
    done = set()
    for line in paf.splitlines():
        f = line.split("\t")
        q, qs, ts = f[0], int(f[2]), int(f[7])
        if q in done or f[4] != "+":
            continue
        done.add(q)
        cs = next(x for x in f[12:] if x.startswith("cs:Z:"))[5:]
        for t in CS.findall(cs):
            if t[0] == ":":
                n = int(t[1:])
                frame[q][qs:qs + n] = np.arange(ts, ts + n)
                qs += n
                ts += n
            elif t[0] == "*":
                frame[q][qs] = ts
                qs += 1
                ts += 1
            elif t[0] == "+":
                n = len(t) - 1
                frame[q][qs:qs + n] = ts
                qs += n
            else:
                ts += len(t) - 1
    for k, v in frame.items():
        idx = np.nonzero(v >= 0)[0]
        if len(idx) == 0:
            continue  # no alignment to E. coli 16S: positions stay -1
        # reference ends outside the alignment carry on at one base per base
        v[: idx[0]] = v[idx[0]] - (idx[0] - np.arange(idx[0]))
        v[idx[-1] + 1:] = v[idx[-1]] + np.arange(1, len(v) - idx[-1])
        v += 1  # 1-based E. coli positions
    return ffa, frame


# ---------------------------------------------------------------- structure, ported from the read bucket analysis

def parse_hit(line):
    f = line.rstrip("\n").split("\t")
    return dict(q=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], t=f[5], tl=int(f[6]), ts=int(f[7]), te=int(f[8]),
                m=int(f[9]), a=int(f[10]), cs=next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), ""))


def scorevec(h):
    """Per query base in [qs,qe): +1 match, -1 mismatch or inserted base; a deletion is charged at its position."""
    if "v" in h:
        return h["v"]
    lens, vals, dels = [], [], []
    q = 0
    for t in CS.findall(h["cs"]):
        c = t[0]
        if c == ":":
            n = int(t[1:])
            lens.append(n)
            vals.append(1)
            q += n
        elif c == "*":
            lens.append(1)
            vals.append(-1)
            q += 1
        elif c == "+":
            n = len(t) - 1
            lens.append(n)
            vals.append(-1)
            q += n
        else:
            dels.append((q, len(t) - 1))
    v = np.repeat(np.array(vals, dtype=np.int32), np.array(lens, dtype=np.int64))
    for pos, n in dels:
        if pos < len(v):
            v[pos] -= n
    if h["st"] == "-":
        v = v[::-1]
    h["v"] = v
    return v


def clip(h, a, b):
    qs, qe = max(h["qs"], a), min(h["qe"], b)
    if qe <= qs:
        return None
    ts, te = h["ts"], h["te"]
    if h["st"] == "+":
        ts += qs - h["qs"]
        te -= h["qe"] - qe
    else:
        te -= qs - h["qs"]
        ts += h["qe"] - qe
    return qs, qe, max(0, ts), max(ts + 1, te)


def segments(hits, a, b):
    """Greedy best-first alignments that don't overlap on the read: the walk of the read along 16S."""
    used, out = [], []
    for h in sorted(hits, key=lambda x: -x["m"]):
        c = clip(h, a, b)
        if not c or c[1] - c[0] < MIN_SEG:
            continue
        qs, qe, ts, te = c
        if sum(max(0, min(qe, y) - max(qs, x)) for x, y in used) > 0.5 * (qe - qs):
            continue
        used.append((qs, qe))
        out.append((qs, qe, h["st"], h["t"], ts, te, h))
    return sorted(out, key=lambda s: s[0])


def fspan(frame, s):
    f = frame[s[3]]
    i, j = min(s[4], len(f) - 1), min(s[5], len(f) - 1)
    return int(f[i]), int(f[j])


def structure(segs, frame, a, b):
    if not segs:
        return "unaligned", dict(strand=".", fstart=-1, fend=-1, head=b - a, tail=0, arm1=0, arm2=0, turn1=-1, turn2=-1)
    sp = [fspan(frame, s) for s in segs]
    first, last = segs[0], segs[-1]
    fs = sp[0][0] if first[2] == "+" else sp[0][1]
    fe = sp[-1][1] if last[2] == "+" else sp[-1][0]
    head, tail = first[0] - a, b - last[1]
    info = dict(strand=first[2], fstart=fs, fend=fe, head=head, tail=tail, arm1=first[1] - first[0], arm2=0, turn1=-1, turn2=-1)
    if len(segs) == 1:
        return ("single_tail" if tail >= MIN_SEG else "single_head" if head >= MIN_SEG else "single"), info
    kinds = []
    for k in range(len(segs) - 1):
        x, y = segs[k], segs[k + 1]
        o = max(0, min(sp[k][1], sp[k + 1][1]) - max(sp[k][0], sp[k + 1][0]))
        small = min(sp[k][1] - sp[k][0], sp[k + 1][1] - sp[k + 1][0])
        if x[2] != y[2]:
            kinds.append("foldback" if o >= 0.25 * small or o >= 50 else "switch")
        else:
            kinds.append("repeat" if o >= 50 else "colinear")
    for name in ("foldback", "repeat", "switch"):
        if name in kinds:
            k = kinds.index(name)
            x, y = segs[k], segs[k + 1]
            info.update(arm1=x[1] - x[0], arm2=y[1] - y[0],
                        turn1=sp[k][1] if x[2] == "+" else sp[k][0],
                        turn2=sp[k + 1][0] if y[2] == "+" else sp[k + 1][1])
            return name, info
    return "colinear", info


def twoparent(hits, a, b, org, strand):
    """Best single reference, best pair within one species, and best pair across species, one breakpoint."""
    L = max(h["ql"] for h in hits)
    per = {}
    for h in hits:
        c = clip(h, a, b)
        if h["st"] != strand or not c or c[1] - c[0] < MIN_SEG:
            continue
        per.setdefault(h["t"], []).append(h)
    if not per:
        return None
    rank = sorted(per, key=lambda t: -max(h["m"] for h in per[t]))
    n, keep = {}, []
    for t in rank:
        o = org[t]
        n[o] = n.get(o, 0) + 1
        if n[o] <= 2:
            keep.append(t)
    M = len(keep)
    K = np.arange(a + MIN_ARM, b - MIN_ARM + 1, GRID)
    P = np.zeros((M, b - a + 1), dtype=np.int32)
    for i, t in enumerate(keep):
        v = np.zeros(L + 1, dtype=np.int32)
        for h in sorted(per[t], key=lambda x: x["m"]):
            v[h["qs"]:h["qe"]] = scorevec(h)
        P[i, 1:] = np.cumsum(v[a:b])
    tot = P[:, -1]
    s1i = int(np.argmax(tot))
    single = int(tot[s1i])
    out = dict(s_single=single, s_intra=single, s_inter="", gain="", bp="", p_left="", org_left="", org_right="", sup_left="", sup_right="")
    if M < 2 or len(K) == 0:
        return out
    ko = K - a
    Pk = P[:, ko]
    S = Pk[:, None, :] - Pk[None, :, :] + tot[None, :, None]
    bk = S.argmax(axis=2)
    bv = np.take_along_axis(S, bk[:, :, None], axis=2)[:, :, 0]
    orgs = np.array([org[t] for t in keep])
    same = orgs[:, None] == orgs[None, :]
    out["s_intra"] = int(max(single, bv[same].max()))
    if (~same).any():
        vinter = np.where(same, -10 ** 9, bv)
        i, j = np.unravel_index(int(vinter.argmax()), vinter.shape)
        k = int(ko[bk[i, j]])
        out.update(s_inter=int(vinter[i, j]), bp=k + a, p_left=keep[i], org_left=org[keep[i]], org_right=org[keep[j]],
                   sup_left=int(P[i, k] - P[j, k]), sup_right=int((tot[j] - P[j, k]) - (tot[i] - P[i, k])))
        out["gain"] = out["s_inter"] - out["s_intra"]
    return out


def describe(hits, L, org, frame):
    segs = segments(hits, 0, L)
    name, info = structure(segs, frame, 0, L)
    row = dict(struct=name, **info)
    row.update(organism="", identity="", s_single="", s_intra="", s_inter="", gain="", bp_frame="", org_left="", org_right="", sup_left="", sup_right="")
    if segs:
        best = max(segs, key=lambda s: s[1] - s[0])
        row["organism"] = org[best[3]]
        row["identity"] = f"{sum(s[6]['m'] for s in segs) / sum(s[6]['a'] for s in segs):.4f}"
    # a fold-back or repeat already explains its own extra length: only single walks get the chimera test
    if name in ("single", "single_tail", "single_head", "colinear") and L >= 2 * MIN_ARM + GRID:
        tp = twoparent(hits, 0, L, org, segs[0][2])
        if tp:
            row.update({k: v for k, v in tp.items() if k in row})
            if tp["bp"] != "":
                h = max((x for x in hits if x["t"] == tp["p_left"] and x["st"] == segs[0][2]), key=lambda x: x["m"])
                frac = min(1, max(0, (tp["bp"] - h["qs"]) / max(1, h["qe"] - h["qs"])))
                tpos = h["ts"] + frac * (h["te"] - h["ts"]) if h["st"] == "+" else h["te"] - frac * (h["te"] - h["ts"])
                fr = frame[h["t"]]
                row["bp_frame"] = int(fr[min(int(tpos), len(fr) - 1)])
    return row


def bucket(row, L, a):
    k = row["struct"]
    if k == "unaligned":
        return "unaligned"
    if k in ("foldback", "repeat", "switch"):
        return k
    if row["gain"] != "" and row["gain"] >= a.min_gain and row["sup_left"] >= a.min_side and row["sup_right"] >= a.min_side \
            and row["s_inter"] >= a.min_fit * L:
        return "chimera"
    if row["s_single"] != "" and row["s_intra"] < a.min_fit * L:
        return "poorfit"
    if k in ("single_tail", "single_head"):
        return "tail"
    if k == "colinear":
        return "colinear"
    if abs(row["fend"] - row["fstart"]) >= (GENE_END - GENE_START) - 2 * FULL_MARGIN:
        return "full"
    return "truncated"


_W = {}


def _init(org, frame, a):
    _W.update(org=org, frame=frame, a=a)


def _work(batch):
    out = []
    for rid, lines in batch:
        hits = [parse_hit(x) for x in lines]
        L = hits[0]["ql"]
        row = describe(hits, L, _W["org"], _W["frame"])
        row["bucket"] = bucket(row, L, _W["a"])
        out.append((rid, row))
    return out


def map_structure(reads_path, refs_fa, org, frame, a):
    # minimap2 takes about 60% of the CPU time, the scoring the rest
    mm_threads = max(1, round(0.6 * a.threads))
    workers = max(1, a.threads - mm_threads)
    mm = subprocess.Popen(["minimap2", "-x", "map-ont", "-c", "--cs=short", "-P", "-t", str(mm_threads), "-K", "100M", refs_fa, reads_path],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)

    def batches(size=200):
        cur, lines, batch = None, [], []
        for line in mm.stdout:
            q = line.split("\t", 1)[0]
            if q != cur:
                if cur is not None:
                    batch.append((cur, lines))
                    if len(batch) >= size:
                        yield batch
                        batch = []
                cur, lines = q, []
            lines.append(line)
        if cur is not None:
            batch.append((cur, lines))
        if batch:
            yield batch

    res = {}
    with Pool(workers, initializer=_init, initargs=(org, frame, a)) as pool:
        for out in pool.imap_unordered(_work, batches(), chunksize=1):
            res.update(out)
    if mm.wait() != 0:
        sys.exit("minimap2 failed on the sample's reads")
    return res


# ---------------------------------------------------------------- main

def main():
    a = parse_args()
    prefix = a.prefix or a.sample
    t0 = time.time()
    timing = {}
    calls, assign = read_emu(a.abundance, a.assignments)
    dehosted = fastq_ids(a.dehosted)
    filtered = fastq_ids(a.filtered)
    emu_in = fastq_ids(a.emu_reads) if a.emu_reads else filtered if filtered is not None else dehosted
    timing["read_sets"] = time.time() - t0

    reads = []
    mismatch = collections.Counter()
    examples = []
    seen = set()
    for rid, seq, qual in fastq(a.reads):
        L, acc = len(seq), mean_accuracy(qual)
        seen.add(rid)
        d = "" if dehosted is None else "kept" if rid in dehosted else "removed"
        rule = fl = ""
        if filtered is not None and d != "removed":
            rule = "short" if L < a.min_length else "long" if L > a.max_length else "qfail" if acc < a.min_accuracy else "pass"
            fl = "kept" if rid in filtered else "dropped"
            if (rule == "pass") != (fl == "kept"):
                mismatch[f"{rule}|{fl}"] += 1
                if len(examples) < 5:
                    examples.append(dict(read_id=rid, length=L, accuracy=round(acc, 3), rule=rule, filtlong=fl))
        reached = d != "removed" and (filtered is None or fl == "kept")
        e, et, ep = "", "", ""
        if reached:
            if emu_in is not None and rid not in emu_in:
                e = "not_sampled"
            elif rid in assign:
                t, p = max(assign[rid], key=lambda x: x[1])
                e, et, ep = "assigned", t, f"{p:.4f}"
            else:
                e = "not_assigned"
        reads.append([rid, L, f"{acc:.3f}", d, rule, fl, e, et, ep])
    timing["fate"] = time.time() - t0
    unknown = {
        "dehosted": len(dehosted - seen) if dehosted else 0,
        "filtered": len(filtered - seen) if filtered else 0,
        "assigned": len(set(assign) - seen),
    }

    refs, struct, why = [], {}, ""
    if a.structure:
        if not calls:
            why = "emu called no species, so there were no references to align to"
        elif not reads:
            why = "no reads"
        else:
            with tempfile.TemporaryDirectory(dir=".") as tmp:
                asvs = asv_references(a.asvs, a.asv_mappings, calls)
                refs = choose_references(calls, assign, a.db, a.reads, tmp, a, asvs)
                timing["references"] = time.time() - t0
                refs_fa, frame = orient_and_frame(refs, a.frame, tmp)
                org = {r["name"]: r["species"] for r in refs}
                struct = map_structure(a.reads, refs_fa, org, frame, a)
                timing["structure"] = time.time() - t0
                with open(f"{prefix}.read_buckets_refs.fasta", "w") as o:
                    for r in refs:
                        o.write(f">{r['name']} {r['db_id']} {r['species']}\n{r['seq']}\n")
    if not os.path.exists(f"{prefix}.read_buckets_refs.fasta"):
        open(f"{prefix}.read_buckets_refs.fasta", "w").close()

    # rows, and what each emu call rests on
    sum_p = collections.defaultdict(collections.Counter)
    n_p50 = collections.defaultdict(collections.Counter)
    with gzip.open(f"{prefix}.read_buckets.tsv.gz", "wt", compresslevel=6) as o:
        o.write("\t".join(COLUMNS) + "\n")
        for r in reads:
            rid, L = r[0], r[1]
            if a.structure and not why:
                s = struct.get(rid)
                if s is None:
                    s = dict(struct="unaligned", strand=".", fstart=-1, fend=-1, head=L, tail=0, arm1=0, arm2=0, turn1=-1, turn2=-1, bucket="unaligned")
                vals = [s.get(k, "") for k in STRUCT_COLUMNS] + [s["bucket"]]
                b = s["bucket"]
            else:
                vals = [""] * (len(STRUCT_COLUMNS) + 1)
                b = ""
            o.write("\t".join(str(x) for x in r + vals) + "\n")
            for t, p in assign.get(rid, ()):
                sum_p[t][b] += p
                n_p50[t][b] += p >= 0.5
    for c in calls:
        c["sum_p"] = dict(sum_p.get(c["tax_id"], {}))
        c["reads_p50"] = dict(n_p50.get(c["tax_id"], {}))

    summary = dict(
        sample=a.sample,
        version=VERSION,
        reads=len(reads),
        structure=bool(a.structure and not why),
        structure_skipped=why,
        stages=dict(deacon=dehosted is not None, filter=filtered is not None, downsampled=bool(a.emu_reads), assignments=bool(assign)),
        filter_rule=dict(min_length=a.min_length, max_length=a.max_length, min_accuracy=a.min_accuracy),
        filter_mismatch=dict(mismatch), filter_mismatch_examples=examples,
        reads_not_in_input=unknown,
        thresholds=dict(min_gain=a.min_gain, min_side=a.min_side, min_fit=a.min_fit, gene_start=GENE_START, gene_end=GENE_END, full_margin=FULL_MARGIN),
        references=[{k: v for k, v in r.items() if k != "seq"} | dict(length=len(r["seq"])) for r in refs],
        calls=calls,
        seconds={k: round(v, 1) for k, v in timing.items()},
    )
    with open(f"{prefix}.read_buckets.json", "w") as o:
        json.dump(summary, o, indent=1)


if __name__ == "__main__":
    main()
