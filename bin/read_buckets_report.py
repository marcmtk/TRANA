#!/usr/bin/env python
"""Run-level read bucket report: where every read of the run went, and what the reads that reach
emu are made of.

Reads the per-sample tables of read_buckets.py, barbell's per-read tables and the run's raw reads,
and writes a self-contained HTML report, per-sample and per-call tables, and MultiQC custom content.
The counts are checked to add up at every stage; a stage that doesn't add up is reported, not hidden.
"""

import argparse
import collections
import csv
import gzip
import json
import math
import os
import re
import statistics
from html import escape as esc

VERSION = "0.2.0"

BUCKET_ORDER = ["full", "truncated", "foldback", "repeat", "tail", "chimera", "switch", "colinear", "poorfit", "unaligned"]
BUCKET_LABEL = {
    "full": "Full length, one species", "truncated": "Truncated", "foldback": "Fold-back", "repeat": "Same-strand repeat",
    "tail": "Unaligned tail or head", "chimera": "Chimera", "switch": "Strand switch, no overlap", "colinear": "Gapped, same strand",
    "poorfit": "Fits no reference", "unaligned": "Unaligned",
}
# chart series: the six categorical slots, full length in neutral, the rare kinds together
BUCKET_SERIES = [("truncated", "Truncated", "--s1"), ("foldback", "Fold-back", "--s2"), ("repeat", "Same-strand repeat", "--s3"),
                 ("tail", "Unaligned tail or head", "--s4"), ("chimera", "Chimera", "--s5"), ("other", "Fits no reference, or other", "--s6")]
OTHER = {"switch", "colinear", "poorfit", "unaligned"}
GROUP_ORDER = ["culture", "mix", "control", "clinical", ""]
GROUP_LABEL = {"culture": "Pure cultures", "mix": "Mixes", "control": "Controls", "clinical": "Clinical samples", "": "Samples"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--tables", nargs="*", default=[], help="<sample>.read_buckets.tsv.gz files")
    p.add_argument("--summaries", nargs="*", default=[], help="<sample>.read_buckets.json files")
    p.add_argument("--barbell-annotation", default="")
    p.add_argument("--barbell-patterns", default="")
    p.add_argument("--barbell-filtered", default="")
    p.add_argument("--barbell-counts", default="", help="barbell_read_counts.tsv, whose run_input row is the run's read count")
    p.add_argument("--sheet", default="", help="barcodes sample sheet, optionally with sample_type and expected columns")
    p.add_argument("--call-check", default="", help="emu_call_check.tsv")
    p.add_argument("--kit", default="")
    p.add_argument("--run-name", default="")
    p.add_argument("--run-path", default="", help="the run's fastq_pass folder; the run is named after the folder above it")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p.parse_args()


# ================================================================== figures (inline SVG, no dependencies)

FS = 13


def nice(v, n=4):
    if v <= 0:
        return 1, 0.25
    raw = v / n
    k = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * k >= raw:
            step = m * k
            break
    return step * math.ceil(v / step - 1e-9), step


def fmt(n):
    return f"{int(round(float(n))):,}"


def pct(a, b, d=1):
    return f"{100 * a / b:.{d}f}%" if b else "–"


CSS = """
:root{color-scheme:light dark;--ink:#0b0b0b;--ink2:#52514e;--ink3:#898781;--grid:#e1e0d9;--base:#c3c2b7;--surface:#fcfcfb;--panel:#f6f5f1;
--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--s4:#eda100;--s5:#e87ba4;--s6:#008300;--neutral:#b9b7ae;--seq0:#eef3f9;--seq1:#2a78d6;
--tagL:#cfccc4;--tagB:#f7dc8f;--tagR:#6e6c66;--ins:#8fb4e3;--ins2:#e7a9c0;--warn:#b54708;--tipbg:#1a1a19;--tipfg:#f0efec}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ink:#f0efec;--ink2:#c3c2b7;--ink3:#8f8d86;--grid:#2c2c2a;--base:#383835;--surface:#1a1a19;--panel:#232321;
--s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--s5:#d55181;--s6:#008300;--neutral:#6f6e68;--seq0:#20262e;--seq1:#3987e5;
--tagL:#77756e;--tagB:#c9a93f;--tagR:#a9a79f;--ins:#4f7fb8;--ins2:#b5587e;--warn:#f79009;--tipbg:#f0efec;--tipfg:#1a1a19}}
*{box-sizing:border-box}
html,body{margin:0;background:var(--surface);color:var(--ink)}
.rb{font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif;max-width:1000px;margin:0 auto;padding:28px 16px 60px}
.rb h1{font-size:26px;line-height:1.25;margin:0 0 6px;font-weight:650}
.rb h2{font-size:20px;margin:42px 0 8px;font-weight:650;padding-top:10px;border-top:1px solid var(--grid)}
.rb h3{font-size:16px;margin:24px 0 6px;font-weight:650}
.rb p{margin:8px 0}.rb ul{margin:8px 0;padding-left:20px}.rb li{margin:5px 0}
.rb code{font:12.5px/1.45 ui-monospace,Menlo,monospace;overflow-wrap:anywhere}
.rb .sub{color:var(--ink2)}.rb .cap{color:var(--ink2);font-size:13px;margin:4px 0 16px}
.rb .warn{border-left:3px solid var(--warn);padding:4px 12px;margin:12px 0;background:var(--panel)}
.rb .note{border-left:3px solid var(--base);padding:4px 12px;margin:12px 0;background:var(--panel)}
.rb table{border-collapse:collapse;font-size:13px;width:100%;margin:8px 0}
.rb th{font-weight:600;text-align:left;border-bottom:1px solid var(--base);padding:4px 8px 4px 0;vertical-align:bottom}
.rb td{border-bottom:1px solid var(--grid);padding:3px 8px 3px 0;vertical-align:top}
.rb td.n,.rb th.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.rb details{margin:6px 0 18px}.rb summary{cursor:pointer;color:var(--ink2);font-size:13px}
.rb .tw{overflow-x:auto}
.rb svg{display:block;width:100%;height:auto;margin:6px 0}
.rb svg text{font-family:system-ui,sans-serif;fill:var(--ink2)}
.rb svg text.mono{font-family:ui-monospace,Menlo,monospace}
.rb svg text.pri{fill:var(--ink)}
.rb svg [data-t]:hover{opacity:.78}
.rb .kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px;margin:12px 0 4px}
.rb .kv b{display:block;font-size:24px;font-weight:600;line-height:1.2}
.rb .kv span{font-size:12.5px;color:var(--ink2)}
.rb .swatch{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px}
.rb .legend{display:flex;flex-wrap:wrap;gap:4px 14px;font-size:12.5px;color:var(--ink2);margin:6px 0 2px}
#rbtip{position:fixed;pointer-events:none;display:none;background:var(--tipbg);color:var(--tipfg);padding:6px 9px;border-radius:6px;font:12.5px/1.4 system-ui,sans-serif;max-width:320px;z-index:9}
#rbtip b{font-weight:650}
"""

JS = """
(function(){var t=document.getElementById('rbtip');
document.addEventListener('pointermove',function(e){var n=e.target.closest&&e.target.closest('[data-t]');
if(!n){t.style.display='none';return}
t.textContent='';var p=n.getAttribute('data-t').split('\\n');var b=document.createElement('b');b.textContent=p[0];t.appendChild(b);
for(var i=1;i<p.length;i++){t.appendChild(document.createElement('br'));t.appendChild(document.createTextNode(p[i]))}
t.style.display='block';var w=t.offsetWidth,h=t.offsetHeight,x=e.clientX+14,y=e.clientY+14;
if(x+w>window.innerWidth-6)x=e.clientX-w-10;if(y+h>window.innerHeight-6)y=e.clientY-h-10;t.style.left=x+'px';t.style.top=y+'px'})})();
"""


def svg(w, h, body, label):
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{esc(label, quote=True)}">{body}</svg>'


def text(x, y, s, size=FS, anchor="start", cls="", weight=None):
    w = f' font-weight="{weight}"' if weight else ""
    c = f' class="{cls}"' if cls else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}"{c}{w}>{esc(str(s))}</text>'


def rect(x, y, w, h, fill, tip=None):
    t = f' data-t="{esc(tip, quote=True)}"' if tip else ""
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0):.1f}" height="{h:.1f}" fill="var({fill})"{t}/>'


def endrect(x, y, w, h, fill, tip=None, r=4):
    t = f' data-t="{esc(tip, quote=True)}"' if tip else ""
    r = min(r, w, h / 2)
    d = f"M{x:.1f},{y:.1f}H{x + w - r:.1f}Q{x + w:.1f},{y:.1f} {x + w:.1f},{y + r:.1f}V{y + h - r:.1f}Q{x + w:.1f},{y + h:.1f} {x + w - r:.1f},{y + h:.1f}H{x:.1f}Z"
    return f'<path d="{d}" fill="var({fill})"{t}/>'


def colrect(x, y, w, h, fill, tip=None, r=3):
    t = f' data-t="{esc(tip, quote=True)}"' if tip else ""
    r = min(r, w / 2, h)
    if h <= 0:
        return ""
    d = f"M{x:.1f},{y + h:.1f}V{y + r:.1f}Q{x:.1f},{y:.1f} {x + r:.1f},{y:.1f}H{x + w - r:.1f}Q{x + w:.1f},{y:.1f} {x + w:.1f},{y + r:.1f}V{y + h:.1f}Z"
    return f'<path d="{d}" fill="var({fill})"{t}/>'


def legend(items):
    return '<div class="legend">' + "".join(f'<span><span class="swatch" style="background:var({c})"></span>{esc(n)}</span>' for n, c in items) + "</div>"


class Raw(str):
    """Table cell holding markup the script wrote itself."""


def table(head, rows, numeric=()):
    h = "".join(f'<th{" class=n" if i in numeric else ""}>{esc(str(x))}</th>' for i, x in enumerate(head))
    b = ""
    for r in rows:
        b += "<tr>" + "".join(f'<td{" class=n" if i in numeric else ""}>{x if isinstance(x, Raw) else esc(str(x))}</td>' for i, x in enumerate(r)) + "</tr>"
    return f'<div class="tw"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def details(summary, inner):
    return f"<details><summary>{esc(summary)}</summary>{inner}</details>"


def stacked(rows, series, width=900, title="", value_fmt=None, xmax=None, row_h=20, label_w=110, right_w=80, tickfmt=None):
    """Horizontal stacked bars. rows = [(group, label, {key: value}, right label)]."""
    if not rows:
        return ""
    x0, x1 = label_w, width - right_w
    y = 26
    maxv, stp = nice(xmax or max(sum(r[2].get(k, 0) for k, _, _ in series) for r in rows) or 1)
    ticks = int(round(maxv / stp))
    top = y
    marks = ""
    prev = None
    for grp, lab, val, tot in rows:
        if grp != prev:
            if prev is not None:
                y += 8
            if grp:
                marks += text(0, y + 11, grp, 11.5, cls="pri", weight=600)
                y += 16
            prev = grp
        marks += text(label_w - 8, y + row_h / 2 + 4, lab, 12.5, "end")
        x = x0
        vals = [(k, n, c) for k, n, c in series if val.get(k, 0) > 0]
        for i, (k, n, c) in enumerate(vals):
            w = (x1 - x0) * val[k] / maxv
            tip = f"{lab}: {n}\n{(value_fmt or fmt)(val[k])}"
            gap = 2 if i < len(vals) - 1 else 0
            marks += (endrect if i == len(vals) - 1 else rect)(x, y, max(w - gap, 0.8), row_h - 4, c, tip)
            x += w
        if tot:
            marks += text(x1 + 6, y + row_h / 2 + 3, tot, 11.5)
        y += row_h
    grid = ""
    for i in range(ticks + 1):
        v = maxv * i / ticks
        xx = x0 + (x1 - x0) * i / ticks
        grid += f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="{top - 4}" y2="{y}" stroke="var(--grid)"/>'
        grid += text(xx, top - 9, (tickfmt or fmt)(v), 11, "middle")
    return svg(width, y + 6, grid + marks, title or "stacked bars")


def _compact(v):
    if v >= 1000:
        return f"{v / 1000:.1f}K".replace(".0K", "K")
    if v == 0:
        return "0"
    if v < 1:
        return f"{v:.1f}"
    return f"{int(round(v))}"


def columns(bins, series, width=900, height=170, title="", xlab="", xfmt=None, yfmt=None, left=52, ticks_every=None):
    """Stacked columns on a shared x. series = [(name, colour, {bin: y})]."""
    bottom, top = 32, 18
    x0, x1 = left, width - 8
    y0, y1 = height - bottom, top
    step = (x1 - x0) / max(1, len(bins))
    tot = {b: sum(s[2].get(b, 0) for s in series) for b in bins}
    my, stp = nice(max(tot.values(), default=1) or 1)
    nt = int(round(my / stp))
    body = ""
    for i in range(nt + 1):
        v = stp * i
        yy = y0 + (y1 - y0) * v / my
        body += f'<line x1="{x0}" x2="{x1}" y1="{yy:.1f}" y2="{yy:.1f}" stroke="var(--grid)"/>'
        body += text(x0 - 6, yy + 4, (yfmt or _compact)(v), 11, "end")
    bw = max(step - 1.5, 1)
    for k, b in enumerate(bins):
        acc = 0
        for name, c, d in series:
            v = d.get(b, 0)
            if not v:
                continue
            h = (y0 - y1) * v / my
            yy = y0 - (y0 - y1) * (acc + v) / my
            tip = f"{(xfmt or str)(b)}: {name}\n{fmt(v) if v >= 1 else round(v, 3)}"
            body += colrect(x0 + k * step + .75, yy, bw, max(h - 1, .6), c, tip)
            acc += v
    every = ticks_every or max(1, len(bins) // 8)
    for k, b in enumerate(bins):
        if k % every == 0:
            body += text(x0 + k * step + step / 2, y0 + 15, (xfmt or str)(b), 11, "middle")
    if xlab:
        body += text((x0 + x1) / 2, height - 3, xlab, 11.5, "middle")
    body += f'<line x1="{x0}" x2="{x1}" y1="{y0}" y2="{y0}" stroke="var(--base)"/>'
    if title:
        body += text(x0, 11, title, 12.5, cls="pri", weight=600)
    return svg(width, height, body, title or "columns")


# read architecture in the style of Beeloo et al. Fig. 1

def tag(x, y, h=14, flank_only=False):
    out, cx = "", x
    for w, c in ((12, "--tagL"), (20, "--tagB"), (22, "--tagR")):
        out += rect(cx, y, w, h, c)
        cx += w
    if flank_only:
        bx = x + 12
        out += f'<path d="M{bx + 10:.1f},{y}l-3,3.5l6,3.5l-6,3.5l3,3.5" stroke="#c0392b" fill="none" stroke-width="1.6"/>'
    return out, 54


def insert(x, y, w, h, fill="--ins", rc=False, label=None):
    a = min(7, w / 2)
    if not rc:
        d = f"M{x:.1f},{y}H{x + w - a:.1f}L{x + w:.1f},{y + h / 2}L{x + w - a:.1f},{y + h}H{x:.1f}Z"
    else:
        d = f"M{x + w:.1f},{y}H{x + a:.1f}L{x:.1f},{y + h / 2}L{x + a:.1f},{y + h}H{x + w:.1f}Z"
    s = f'<path d="{d}" fill="var({fill})"/>'
    if label:
        s += text(x + w / 2, y + h - 3.5, label, 10.5, "middle", cls="pri")
    return s


def arch_panel(rows, width=900, row_h=46, x0=130, scale=0.215, h=14):
    body, y = "", 8
    for name, parts, pattern in rows:
        body += text(0, y + h - 2, name, 12.5, cls="pri", weight=600)
        x, boxes, ks = x0, [], x0
        for p in parts:
            k = p[0]
            if k == "keep":
                ks = x
            elif k == "/keep":
                boxes.append((ks, x))
            elif k in ("tag", "flank"):
                s, w = tag(x, y, h, flank_only=k == "flank")
                body += s
                x += w
            elif k.startswith("ins"):
                w = p[1] * scale
                body += insert(x, y, w, h, "--ins2" if k.startswith("ins2") else "--ins", rc=k.endswith("rc"), label=p[2] if len(p) > 2 else None)
                x += w
            elif k == "tail":
                w = p[1] * scale
                body += f'<rect x="{x:.1f}" y="{y + 3}" width="{w:.1f}" height="{h - 6}" fill="none" stroke="var(--ink3)"/>'
                x += w
        for a_, b_ in boxes:
            body += f'<rect x="{a_ - 1:.1f}" y="{y - 3}" width="{b_ - a_ + 2:.1f}" height="{h + 6}" fill="none" stroke="var(--ink2)" stroke-dasharray="3 2" rx="2"/>'
        body += text(x0, y + h + 15, pattern, 11, cls="mono")
        y += row_h
    return svg(width, y, body, "read architectures")


# ================================================================== data

def read_sheet(path):
    """sample -> (type, [expected names]); both optional columns."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as h:
        for r in csv.DictReader(h):
            r = {(k or "").strip().lower(): (v or "").strip() for k, v in r.items()}
            s = r.get("sample_id") or r.get("sample")
            if not s:
                continue
            t = r.get("sample_type", "").lower()
            exp = [x.strip() for x in re.split(r"[;|]", r.get("expected", "")) if x.strip()]
            out[s] = (t if t in GROUP_ORDER else ("" if not t else t), exp)
    return out


def expected_status(species, expected):
    if not expected:
        return ""
    s = species.lower()
    for e in expected:
        e = e.lower()
        if s == e or (" " not in e and s.split(" ")[0] == e):
            return "expected"
    return "not expected"


def raw_id(rid):
    return rid.split("_")[0]


def run_input(path):
    """The run's read count, from the run_input row of barbell_read_counts.tsv."""
    if path and os.path.exists(path):
        with open(path) as h:
            for line in h:
                f = line.rstrip("\n").split("\t")
                if f[0] == "run_input" and f[-1].isdigit():
                    return int(f[-1])
    return None


def pattern_group(flanks):
    """Our grouping of barbell's patterns, from the flank positions of each match."""
    n = len(flanks)
    if n == 0:
        return "None"
    if n == 1:
        s, e, mt, lab = flanks[0]
        return ("Ftag_left" if mt == "Ftag" else "Fflank_left") if s <= 250 else "single_elsewhere"
    if n > 2:
        return "three_plus"
    a, b = flanks
    if a[0] > 250:
        return "two_elsewhere"
    t = ("T" if a[2] == "Ftag" else "F") + ("T" if b[2] == "Ftag" else "F")
    gap = b[0] - a[0]
    if gap <= 250:
        return "tandem_" + t
    if 1350 <= gap <= 1750:
        return "amplicon_" + t
    return "midread_" + t


PATTERN_GROUPS = [
    ("Ftag_left", "One Ftag at the start", "Ftag[fw,*,@left(0..250)]"),
    ("Fflank_left", "One flank at the start, no barcode", "Fflank[fw,*,@left(0..250)]"),
    ("None", "No tag found", "None"),
    ("amplicon_TF", "Ftag, then a flank one amplicon in", "Ftag…Fflank@prev_left(1350..1750)"),
    ("amplicon_TT", "Ftag, then another Ftag one amplicon in", "Ftag…Ftag@prev_left(1350..1750)"),
    ("amplicon_FT", "Flank, then an Ftag one amplicon in", "Fflank…Ftag@prev_left(1350..1750)"),
    ("amplicon_FF", "Two flanks one amplicon apart", "Fflank…Fflank@prev_left(1350..1750)"),
    ("tandem_TT", "Two Ftags within 250 bp", "Ftag…Ftag@prev_left(0..250)"),
    ("tandem_TF", "Ftag and flank within 250 bp", "Ftag…Fflank@prev_left(0..250)"),
    ("tandem_FT", "Flank and Ftag within 250 bp", "Fflank…Ftag@prev_left(0..250)"),
    ("tandem_FF", "Two flanks within 250 bp", "Fflank…Fflank@prev_left(0..250)"),
    ("midread_TT", "Two Ftags, other spacing", "Ftag…Ftag, other gap"),
    ("midread_TF", "Ftag and flank, other spacing", "Ftag…Fflank, other gap"),
    ("midread_FT", "Flank and Ftag, other spacing", "Fflank…Ftag, other gap"),
    ("midread_FF", "Two flanks, other spacing", "Fflank…Fflank, other gap"),
    ("single_elsewhere", "One match, not at the start", "…@left(250..)"),
    ("two_elsewhere", "Two matches, first not at the start", "…"),
    ("three_plus", "Three or more matches", "…"),
]


def load(a):
    sheet = read_sheet(a.sheet)
    S = {}
    for p in a.summaries:
        d = json.load(open(p))
        S[d["sample"]] = d
    order = [s for s in sheet if s in S] + sorted(s for s in S if s not in sheet)
    group = {s: sheet.get(s, ("", []))[0] for s in order}
    expected = {s: sheet.get(s, ("", []))[1] for s in order}

    check = {}
    if a.call_check and os.path.exists(a.call_check) and os.path.getsize(a.call_check):
        with open(a.call_check) as h:
            for r in csv.DictReader(h, delimiter="\t"):
                check[(r["sample"], r["tax_id"])] = r["status"]

    # per-read tables
    T = dict(fate=collections.defaultdict(collections.Counter), bucket=collections.defaultdict(collections.Counter),
             bucket_emu=collections.defaultdict(collections.Counter), bucket_long=collections.defaultdict(collections.Counter),
             strand=collections.defaultdict(collections.Counter), pairs=collections.defaultdict(collections.Counter),
             short_end=collections.defaultdict(collections.Counter), turn=collections.defaultdict(collections.Counter),
             bp=collections.Counter(), acc_sample=collections.defaultdict(collections.Counter), short_genus=collections.Counter())
    raw_ids = collections.defaultdict(set)
    for p in a.tables:
        with gzip.open(p, "rt") as h:
            head = h.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(head)}
            s = os.path.basename(p)[: -len(".read_buckets.tsv.gz")]
            for line in h:
                f = line.rstrip("\n").split("\t")
                rid, L = f[ix["read_id"]], int(f[ix["length"]])
                raw_ids[s].add(raw_id(rid))
                d, rule, fl, e = f[ix["deacon"]], f[ix["filter"]], f[ix["filtlong"]], f[ix["emu"]]
                if d == "removed":
                    k = "deacon"
                elif fl == "dropped":
                    k = rule if rule in ("short", "long", "qfail") else "filter_other"
                elif e in ("assigned", "not_assigned", "not_sampled"):
                    k = "emu_" + e
                else:
                    k = "kept"
                T["fate"][s][k] += 1
                lc = rule if rule in ("short", "long", "qfail") else "window" if rule == "pass" else ""
                if lc in ("window", "qfail"):
                    T["acc_sample"][s][int(float(f[ix["accuracy"]]))] += 1
                b = f[ix["bucket"]]
                if not b:
                    continue
                T["bucket"][s][b] += 1
                if e:
                    T["bucket_emu"][s][b] += 1
                if lc == "long":
                    T["bucket_long"][s][b] += 1
                st = f[ix["strand"]]
                if st in ("+", "-"):
                    T["strand"][s][st] += 1
                if b == "chimera":
                    x, y = sorted((f[ix["org_left"]], f[ix["org_right"]]))
                    T["pairs"][s][f"{x}|{y}"] += 1
                    if f[ix["bp_frame"]] not in ("", "-1"):
                        T["bp"][int(f[ix["bp_frame"]]) // 25 * 25] += 1
                if lc == "short" and f[ix["struct"]] == "single" and b in ("truncated", "full"):
                    T["short_end"][st][int(f[ix["fend"]]) // 25 * 25] += 1
                    T["short_genus"][f[ix["organism"]].split(" ")[0]] += 1
                if b == "foldback":
                    T["turn"][st][min(int(f[ix["turn1"]]) // 50 * 50, 1500)] += 1
    return S, order, group, expected, check, T, raw_ids


def barbell_tables(a, raw_ids):
    """Run-level flow and pattern groups from barbell's tables and the raw reads."""
    if not (a.barbell_annotation and os.path.exists(a.barbell_annotation)):
        return None
    flanks = collections.defaultdict(list)
    length = {}
    with gzip.open(a.barbell_annotation, "rt") as h:
        head = h.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        for line in h:
            f = line.rstrip("\n").split("\t")
            flanks[f[ix["read_id"]]].append((int(f[ix["read_start_flank"]]), int(f[ix["read_end_flank"]]), f[ix["match_type"]], f[ix["label"]]))
            length[f[ix["read_id"]]] = int(f[ix["read_len"]])
    pattern = {}
    if a.barbell_patterns and os.path.exists(a.barbell_patterns):
        with gzip.open(a.barbell_patterns, "rt") as h:
            for line in h:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[0] != "read_id":
                    pattern[f[0]] = f[1].replace(" ", "")
    written = set()
    if a.barbell_filtered and os.path.exists(a.barbell_filtered):
        with gzip.open(a.barbell_filtered, "rt") as h:
            h.readline()
            for line in h:
                written.add(line.split("\t", 1)[0])
    raw = run_input(a.barbell_counts)
    in_sample = set().union(*raw_ids.values()) if raw_ids else set()
    B = dict(raw=raw, annotated=len(flanks), written=len(written), in_samples=len(in_sample),
             written_not_in_sample=len(written - in_sample), in_sample_not_written=len(in_sample - written) if written else 0)
    grp = collections.Counter()
    grp_written = collections.Counter()
    pat = collections.Counter()
    pat_written = collections.Counter()
    pat_len = collections.defaultdict(list)
    same = collections.Counter()
    sample_grp = collections.defaultdict(collections.Counter)
    owner = {}
    for s, ids in raw_ids.items():
        for r in ids:
            owner[r] = s
    for rid, fl in flanks.items():
        fl.sort()
        g = pattern_group(fl)
        pt = pattern.get(rid, "None")
        grp[g] += 1
        pat[pt] += 1
        if len(pat_len[pt]) < 5000:
            pat_len[pt].append(length[rid])
        if rid in written:
            grp_written[g] += 1
            pat_written[pt] += 1
        if g in ("amplicon_TT", "tandem_TT", "midread_TT"):
            same[(g, "same" if fl[0][3] == fl[1][3] else "different")] += 1
        if rid in owner:
            sample_grp[owner[rid]][g] += 1
    # barbell's tables only list reads with a match; the rest of the run had no tag at all
    if raw is not None:
        grp["None"] += raw - len(flanks)
        pat["None"] += raw - len(flanks)
    B["no_tag"] = grp["None"] if raw is not None else None
    B["refused"] = B["annotated"] - B["written"]
    B["groups"] = grp
    B["groups_written"] = grp_written
    B["patterns"] = [(p, n, pat_written[p], int(statistics.median(pat_len[p])) if pat_len[p] else None) for p, n in pat.most_common(15)]
    B["n_patterns"] = len(pat)
    B["same"] = same
    B["sample_groups"] = sample_grp
    return B


# ================================================================== report

def main():
    a = parse_args()
    if a.run_path and not a.run_name:
        path = os.path.normpath(a.run_path)
        a.run_name = os.path.basename(os.path.dirname(path)) if os.path.basename(path) == "fastq_pass" else os.path.basename(path)
    S, order, group, expected, check, T, raw_ids = load(a)
    B = barbell_tables(a, raw_ids)
    structure = any(S[s]["structure"] for s in order)
    typed = any(group[s] for s in order)
    gl = lambda s: GROUP_LABEL.get(group[s], group[s].capitalize()) if typed else ""
    ordered = sorted(order, key=lambda s: (GROUP_ORDER.index(group[s]) if group[s] in GROUP_ORDER else len(GROUP_ORDER), order.index(s)))

    html = []
    w = html.append

    def para(s):
        w(f"<p>{s}</p>")

    title = f"Read buckets{(' — ' + esc(a.run_name)) if a.run_name else ''}"
    w(f'<h1>{title}</h1>')
    w(f'<p class="sub">Where every read of the run went, and what the reads behind each emu call are made of. '
      f'{len(order)} samples{", kit " + esc(a.kit) if a.kit else ""}. read_buckets {VERSION}.</p>')

    # ---------------------------------------------------------------- checks that the stages add up
    problems = []
    for s in order:
        d = S[s]
        n = sum(T["fate"][s].values())
        if n != d["reads"]:
            problems.append(f"{esc(s)}: the read table has {n:,} rows but the sample had {d['reads']:,} reads.")
        mm = d.get("filter_mismatch") or {}
        if mm:
            problems.append(f"{esc(s)}: the length and quality rule ({d['filter_rule']['min_length']}–{d['filter_rule']['max_length']} bp, mean accuracy ≥ {d['filter_rule']['min_accuracy']}%) "
                            f"and filtlong disagree on {sum(mm.values()):,} reads ({', '.join(f'{k}: {v}' for k, v in mm.items())}).")
        for k, v in (d.get("reads_not_in_input") or {}).items():
            if v:
                problems.append(f"{esc(s)}: {v:,} reads in the {k} set are not among the sample's reads.")
    if B and B["in_sample_not_written"]:
        problems.append(f"{B['in_sample_not_written']:,} reads in the samples are missing from barbell's filtered table.")
    if problems:
        w('<div class="warn"><b>The counts don\'t add up everywhere.</b><ul>' + "".join(f"<li>{x}</li>" for x in problems) + "</ul></div>")

    # ---------------------------------------------------------------- headline numbers
    tot = collections.Counter()
    for s in order:
        tot.update(T["fate"][s])
    n_all = sum(tot.values())
    kv = []
    if B and B["raw"] is not None:
        kv.append((fmt(B["raw"]), "reads in the run"))
    kv.append((fmt(n_all), "reads in the samples"))
    reach = tot["emu_assigned"] + tot["emu_not_assigned"]
    kv.append((pct(reach, n_all), "reach emu"))
    if structure:
        chim = sum(T["bucket_emu"][s]["chimera"] for s in order)
        n_emu = sum(sum(T["bucket_emu"][s].values()) for s in order)
        kv.append((pct(chim, n_emu), "of the reads emu sees are chimeric"))
    w('<div class="kv">' + "".join(f"<div><b>{x}</b><span>{y}</span></div>" for x, y in kv) + "</div>")
    if structure:
        w('<div class="note"><b>Uncalibrated.</b> The chimera call and the read structure use thresholds chosen on one run, 2026-07-02 (pure cultures as negatives, reads of two cultures spliced together as positives), '
          'and haven\'t been checked on another run with known composition. Treat the rates as descriptive. Calls between closely related species (two staphylococci, two klebsiellae) are under-counted.</div>')

    # ---------------------------------------------------------------- flow
    w("<h2>Where the reads go</h2>")
    if B:
        rows = []
        if B["raw"] is not None:
            rows.append(("Reads in the run (barbell's input)", B["raw"]))
            rows.append(("No tag found (pattern None)", B["no_tag"]))
        rows.append(("Tag found, pattern refused by barbell", B["refused"]))
        rows.append(("Kept by barbell", B["written"]))
        rows.append(("… in a barcode with no sample, or empty after trimming", B["written_not_in_sample"]))
        rows.append(("… in the samples", B["in_samples"]))
        rows.append(("Reads in the samples, after barbell split some in two", n_all))
        w(table(["Stage", "Reads"], [(x, fmt(y)) for x, y in rows if y is not None], numeric=(1,)))
    fate_series = [("deacon", "Removed as human", "--s4"), ("short", f"Under the length window", "--s1"), ("long", "Over the length window", "--s2"),
                   ("qfail", "Below the accuracy cut-off", "--s3"), ("filter_other", "Dropped by filtlong for another reason", "--s6"),
                   ("emu_not_sampled", "Not in the downsample", "--s6"), ("emu_not_assigned", "Reached emu, not assigned", "--s5"),
                   ("emu_assigned", "Assigned by emu", "--neutral"), ("kept", "Kept", "--neutral")]
    used = [x for x in fate_series if any(T["fate"][s][x[0]] for s in order)]
    rows = []
    for s in ordered:
        f = T["fate"][s]
        n = sum(f.values()) or 1
        rows.append((gl(s), s, {k: 100 * f[k] / n for k, _, _ in used}, fmt(n)))
    w(legend([(n, c) for _, n, c in used]) + stacked(rows, used, xmax=100, value_fmt=lambda v: f"{v:.2f}% of the sample", tickfmt=lambda v: f"{v:.0f}%", title="Fate of each sample's reads"))
    rule = S[order[0]]["filter_rule"] if order else {}
    w(f'<p class="cap">Share of each sample\'s reads by where they stop. The length window is {rule.get("min_length", "")}–{rule.get("max_length", "")} bp and the accuracy cut-off a mean per-base accuracy of '
      f'{rule.get("min_accuracy", "")}%, as filtlong applies them. Right: reads in the sample.</p>')
    w(details("Table: fate of each sample's reads", table(["Sample"] + [n for _, n, _ in used] + ["Reads"],
                                                          [[s] + [fmt(T["fate"][s][k]) for k, _, _ in used] + [fmt(sum(T["fate"][s].values()))] for s in ordered],
                                                          numeric=tuple(range(1, len(used) + 2)))))

    # ---------------------------------------------------------------- notation
    if B or structure:
        w("<h2>Notation</h2>")
        para("Barcode patterns follow Beeloo et al., <i>Barbell reveals and resolves demultiplexing and trimming issues in Nanopore data</i>: a read is a chain of matches joined by "
             "<code>__</code>, each <code>type[orientation, label, position]</code>. For the 16S insert this report adds elements in the same style, marked <i>ours</i>.")
        terms = [
            ("Ftag", "paper", "Left flank, barcode and right flank found; the barcode is assigned."),
            ("Fflank", "paper", "The flank is found but no barcode passes barbell's cut-offs."),
            ("@left(0..250)", "paper", "Within 0–250 bp of the read start. <code>@prev_left(i..j)</code> counts from the previous match."),
            ("None", "paper", "No tag found anywhere in the read."),
            ("Ins[fw|rc, A, s..e]", "ours", "An aligned stretch of 16S. <i>fw</i>: the read runs along the gene from the forward primer; <i>rc</i>: from the reverse primer. "
                                            "<i>A</i> is the species, <i>s..e</i> the E. coli 16S positions where the stretch starts and ends, in read order."),
            ("|", "ours", "A breakpoint where the species changes: <code>Ins[fw,A,28..k]|Ins[fw,B,k..1492]</code>."),
            ("Tail(n)", "ours", "n bases after the last aligned stretch that align to no reference."),
            ("Fold-back", "ours", "Two inserts in opposite orientations over the same part of the gene: the read turns and reads back along the other strand."),
            ("Same-strand repeat", "ours", "Two inserts in the same orientation over the same part of the gene."),
            ("Chimera", "ours", "One walk along the gene whose two ends fit two different species better than any one species does."),
        ]
        w(table(["Term", "Source", "Meaning"], [(Raw(f"<code>{esc(x)}</code>"), y, Raw(z)) for x, y, z in terms]))
        arch = [
            ("Expected", [("tag",), ("keep",), ("ins", 1465, "Ins[fw,A,28..1492]"), ("/keep",)], "Ftag[fw,*,@left(0..250)]__Ins[fw,A,28..1492]"),
            ("Truncated", [("tag",), ("keep",), ("insrc", 620, "Ins[rc,A,1492..870]"), ("/keep",)], "Ftag[fw,*,@left(0..250)]__Ins[rc,A,1492..870]"),
            ("Chimera", [("tag",), ("keep",), ("insrc", 700, "Ins[rc,A,1492..790]"), ("ins2rc", 765, "Ins[rc,B,790..28]"), ("/keep",)],
             "Ftag[fw,*,@left(0..250)]__Ins[rc,A,1492..790]|Ins[rc,B,790..28]"),
            ("Fold-back", [("tag",), ("keep",), ("insrc", 1465, "Ins[rc,A,1492..28]"), ("/keep",), ("flank",), ("ins", 1000, "Ins[fw,A,28..1030]")],
             "Ftag[fw,*,@left(0..250)]__Ins[rc,A,1492..28]__Fflank[fw,*,@prev_left(1350..1750)]__Ins[fw,A,28..1030]"),
            ("Two molecules", [("tag",), ("keep",), ("ins", 1465, "Ins[fw,A,28..1492]"), ("/keep",), ("tag",), ("keep",), ("insrc", 1465, "Ins[rc,B,1492..28]"), ("/keep",)],
             "Ftag[fw,AB01,@left(0..250)]__Ins[fw,A,28..1492]__Ftag[fw,AB05,@prev_left(1350..1750)]__Ins[rc,B,1492..28]"),
            ("Unaligned tail", [("tag",), ("keep",), ("ins", 1465, "Ins[fw,A,28..1492]"), ("tail", 300), ("/keep",)], "Ftag[fw,*,@left(0..250)]__Ins[fw,A,28..1492]__Tail(300)"),
        ]
        w(arch_panel(arch))
        w('<p class="cap">Read architectures. Barcode tags as in the paper\'s Fig. 1: grey flanks, yellow barcode, a red zigzag where the barcode is missing. '
          'The insert is an arrow towards the gene\'s 3′ end; a dashed box marks what barbell keeps after trimming. Tags are not to scale; inserts are.</p>')

    # ---------------------------------------------------------------- barcode patterns
    if B:
        w("<h2>Barcode patterns</h2>")
        g = B["groups"]
        ntot = sum(g.values()) or 1
        rows = []
        for k, lab, pat in PATTERN_GROUPS:
            if g[k]:
                rows.append((Raw(f"<code>{esc(pat)}</code>"), lab, fmt(g[k]), pct(g[k], ntot, 2), fmt(B["groups_written"][k])))
        w(table(["Pattern", "Meaning", "Reads", "Share", "Kept by barbell"], rows, numeric=(2, 3, 4)))
        same = B["same"]
        if same[("amplicon_TT", "same")] + same[("amplicon_TT", "different")]:
            sa, di = same[("amplicon_TT", "same")], same[("amplicon_TT", "different")]
            para(f"Of the {fmt(sa + di)} reads with a second Ftag one amplicon in, {fmt(di)} carry a different barcode the second time: two molecules read as one. "
                 f"The other {fmt(sa)} repeat the first barcode. A flank without a barcode one amplicon in (<code>Ftag…Fflank</code>, {fmt(g['amplicon_TF'])} reads) "
                 "was the signature of a fold-back in the 2026-07-02 run, where the read turns at the far end of the amplicon.")
        w(details(f"Table: the {len(B['patterns'])} most common of {B['n_patterns']:,} barbell patterns",
                  table(["Pattern", "Reads", "Kept by barbell", "Median raw length (bp)"],
                        [(Raw(f"<code>{esc(p)}</code>"), fmt(n), fmt(wr), fmt(m) if m else "–") for p, n, wr, m in B["patterns"]], numeric=(1, 2, 3))))
        ser = [("Fflank_left", "Flank only at the start", "--s6"), ("amplicon_TF", "Ftag…Fflank one amplicon apart", "--s2"),
               ("amplicon_TT", "Ftag…Ftag one amplicon apart", "--s1"), ("tandem_TT", "Two Ftags within 250 bp", "--s3"), ("rest", "Anything else", "--s4")]
        rows = []
        for s in ordered:
            v = B["sample_groups"].get(s, collections.Counter())
            n = sum(v.values())
            if not n:
                continue
            main_ = v["Ftag_left"]
            d = {k: 100 * v[k] / n for k, _, _ in ser if k != "rest"}
            d["rest"] = 100 * (n - main_ - sum(v[k] for k, _, _ in ser if k != "rest")) / n
            rows.append((gl(s), s, d, f"{pct(main_, n)} plain"))
        if rows:
            w(legend([(n, c) for _, n, c in ser]) + stacked(rows, ser, value_fmt=lambda v: f"{v:.2f}% of the sample's reads", tickfmt=lambda v: f"{v:g}%",
                                                            title="Barcode patterns other than one Ftag at the start, by sample"))
            w('<p class="cap">Reads of each sample whose raw read had a pattern other than <code>Ftag[fw,*,@left(0..250)]</code>, as a share of the sample\'s raw reads. Right: share with the plain pattern.</p>')

    # ---------------------------------------------------------------- structure
    if structure:
        w("<h2>What the reads are made of</h2>")
        para("Each read is aligned to its sample's 16S references: savont's ASVs for the sample, named after the species of their closest database references, "
             "and, for each species emu called, database references that most of the species' own reads hit best (one if an ASV already has the species, up to three if not). "
             "The read's path along the gene names its bucket. A read can be in only one bucket; fold-backs and repeats are named before the chimera test, "
             "because their extra length already explains them.")
        ref_rows = []
        for s in ordered:
            for r in S[s].get("references", []):
                ref_rows.append((s, Raw(f"<i>{esc(r['species'])}</i>"), "savont ASV" if r.get("source") == "asv" else "emu database", r["db_id"],
                                 fmt(r["hits"]) if r.get("source") == "asv" else f"{r['hits']} of the species' reads"))
        if ref_rows:
            w(details("Table: the references each sample was aligned to", table(["Sample", "Species", "Source", "Reference", "Depth or reads that hit it best"], ref_rows)))

        def brow(src, s):
            v = src[s]
            n = sum(v.values())
            if not n:
                return None
            d = collections.Counter()
            for k, x in v.items():
                d["other" if k in OTHER else k] += x
            return (gl(s), s, {k: 100 * d[k] / n for k, _, _ in BUCKET_SERIES}, fmt(n))

        for src, ttl, cap in ((T["bucket_emu"], "Reads that reach emu", "the reads emu got"), (T["bucket"], "All reads", "all of the sample's reads")):
            rows = [r for r in (brow(src, s) for s in ordered) if r]
            w(legend([(n, c) for _, n, c in BUCKET_SERIES]) + stacked(rows, BUCKET_SERIES, value_fmt=lambda v: f"{v:.2f}% of the reads", tickfmt=lambda v: f"{v:g}%",
                                                                      title=f"{ttl}: buckets other than full length"))
            w(f'<p class="cap">Share of {cap} in each bucket. Full-length reads of one species make up the rest. Right: reads.</p>')
        cols = [b for b in BUCKET_ORDER if any(T["bucket"][s][b] for s in order)]
        w(details("Table: reads per bucket, all reads", table(["Sample"] + [BUCKET_LABEL[b] for b in cols],
                                                             [[s] + [fmt(T["bucket"][s][b]) for b in cols] for s in ordered], numeric=tuple(range(1, len(cols) + 1)))))
        w(details("Table: reads per bucket, reads that reach emu", table(["Sample"] + [BUCKET_LABEL[b] for b in cols],
                                                                        [[s] + [fmt(T["bucket_emu"][s][b]) for b in cols] for s in ordered], numeric=tuple(range(1, len(cols) + 1)))))

        # primer
        w("<h3>Starting primer</h3>")
        rows = []
        for s in ordered:
            v = T["strand"][s]
            n = v["+"] + v["-"]
            if n:
                rows.append((gl(s), s, {"-": 100 * v["-"] / n, "+": 100 * v["+"] / n}, f"{pct(v['-'], n, 0)} reverse"))
        ser = [("-", "From the reverse primer (rc)", "--s1"), ("+", "From the forward primer (fw)", "--s2")]
        w(legend([(n, c) for _, n, c in ser]) + stacked(rows, ser, xmax=100, value_fmt=lambda v: f"{v:.1f}%", tickfmt=lambda v: f"{v:.0f}%", title="Starting primer by sample"))
        w('<p class="cap">The primer a read starts from, read off the strand its first stretch aligns to. Barbell trims the primer with the flank, so it is not read directly; '
          'in the 2026-07-02 run, where it was, strand and primer agreed for 99.95% of reads.</p>')

        # short reads
        if T["short_end"]:
            w("<h3>Short reads: where they stop</h3>")
            bins = list(range(0, 1525, 25))
            fw, rv = T["short_end"].get("+", {}), T["short_end"].get("-", {})
            nf, nr = sum(fw.values()), sum(rv.values())
            w(columns(bins, [("From the forward primer", "--s2", dict(fw))], height=150, title=f"Short reads from the forward primer ({fmt(nf)})",
                      xlab="E. coli 16S position where the read ends (25 bp bins)", ticks_every=8))
            w(columns(bins, [("From the reverse primer", "--s1", dict(rv))], height=150, title=f"Short reads from the reverse primer ({fmt(nr)})",
                      xlab="E. coli 16S position where the read ends (25 bp bins)", ticks_every=8))
            top = ", ".join(f"<i>{esc(k)}</i> {fmt(v)}" for k, v in T["short_genus"].most_common(6))
            w(f'<p class="cap">Reads under the length window that align as one stretch, by the E. coli position of their last aligned base. Reads from the forward primer start near position 28 and run up the gene; '
              f'reads from the reverse primer start near 1492 and run down it. A stop position shared across species points to the library or the chemistry rather than to one organism. '
              f'Short reads by genus: {top}.</p>')

        # long reads
        if any(T["bucket_long"][s] for s in order):
            w("<h3>Long reads</h3>")
            rows = []
            for s in ordered:
                v = T["bucket_long"][s]
                n = sum(v.values())
                if n:
                    d = collections.Counter()
                    for k, x in v.items():
                        d["other" if k in OTHER else k] += x
                    rows.append((gl(s), s, {k: d[k] for k, _, _ in BUCKET_SERIES}, fmt(n)))
            w(legend([(n, c) for _, n, c in BUCKET_SERIES]) + stacked(rows, BUCKET_SERIES, title="Reads over the length window by bucket"))
            w('<p class="cap">Reads over the length window that barbell kept, by bucket; reads, not shares. Right: all reads over the window, full-length ones included.</p>')
            if T["turn"]:
                for st, ttl, c in (("+", "Fold-backs that start at the forward primer", "--s2"), ("-", "Fold-backs that start at the reverse primer", "--s1")):
                    d = T["turn"].get(st, {})
                    if d:
                        w(columns(list(range(0, 1550, 50)), [("Turn position", c, dict(d))], height=150, title=f"{ttl} ({fmt(sum(d.values()))})",
                                  xlab="E. coli 16S position where the read turns (50 bp bins)", ticks_every=4))
                w('<p class="cap">Where fold-backs turn, all lengths. A turn at the far end of the amplicon reads one strand and then the other.</p>')

        # chimeras
        w("<h2>Chimeras (uncalibrated)</h2>")
        thr = S[order[0]]["thresholds"]
        para(f"A read is chimeric when two species with one breakpoint explain it by at least {thr['min_gain']} score units more than the best single species does, "
             f"each side gains at least {thr['min_side']}, and the pair fits at least {thr['min_fit']} per base. Score: +1 per matching base, −1 per mismatch, inserted or deleted base. "
             "A read is scored only if it walks along the gene once; fold-backs and repeats are not.")
        rows = []
        tab = []
        for s in ordered:
            v, ve = T["bucket"][s], T["bucket_emu"][s]
            n, ne = sum(v.values()), sum(ve.values())
            if not n:
                continue
            rows.append((gl(s), s, {"all": 100 * v["chimera"] / n}, pct(v["chimera"], n, 2)))
            top = ", ".join(f"{esc(p.replace('|', ' × '))} {fmt(c)}" for p, c in T["pairs"][s].most_common(3))
            tab.append((s, fmt(n), fmt(v["chimera"]), pct(v["chimera"], n, 2), fmt(ne), fmt(ve["chimera"]), pct(ve["chimera"], ne, 2), Raw(top)))
        w(stacked(rows, [("all", "Chimeric reads", "--s5")], value_fmt=lambda v: f"{v:.2f}% of the sample's reads", tickfmt=lambda v: f"{v:g}%", title="Chimeric reads by sample"))
        cult = [s for s in order if group[s] == "culture"]
        if cult:
            c = sum(T["bucket"][s]["chimera"] for s in cult)
            n = sum(sum(T["bucket"][s].values()) for s in cult)
            para(f"In the pure cultures {fmt(c)} of {fmt(n)} reads ({pct(c, n, 3)}) are called chimeric. This is a weak negative control: a culture's references are its own ASVs and the species emu called in it, "
                 "mostly one species and the database names that tie with it, so a false call between unrelated species can't happen there. "
                 "In the 2026-07-02 calibration, where every culture was scored against all of the run's species, 29 of about 427,000 pure-culture reads were called chimeric.")
        w(table(["Sample", "Reads", "Chimeric", "Share", "Reach emu", "Chimeric", "Share", "Most common pairs"], tab, numeric=(1, 2, 3, 4, 5, 6)))
        if T["bp"]:
            bins = list(range(0, 1525, 25))
            w(columns(bins, [("Breakpoint", "--s5", dict(T["bp"]))], height=160, title=f"Where the species change, all chimeric reads ({fmt(sum(T['bp'].values()))})",
                      xlab="E. coli 16S position of the breakpoint (25 bp bins)", ticks_every=8))

        # emu calls
        w("<h2>What each emu call rests on</h2>")
        para("For each species emu called, its reads by bucket, weighted by emu's assignment probability. "
             "A call whose reads are mostly chimeric is likely made by chimeras of the sample's main species; "
             + ("the label from emu's call check (near-exact reads) is shown next to it." if check else "turn on <code>--check_emu_calls</code> to see emu's own near-exact read label next to it."))
        rows, tab = [], []
        for s in ordered:
            for c in S[s]["calls"]:
                sp = c["sum_p"]
                t = sum(sp.values()) or 1
                d = collections.Counter()
                for k, x in sp.items():
                    d["other" if k in OTHER or not k else k] += x
                st = expected_status(c["species"], expected[s])
                lab = check.get((s, c["tax_id"]), "")
                share = 100 * sp.get("chimera", 0) / t
                grp = st.capitalize() if st else ("No species listed in the sample sheet" if any(expected.values()) else "")
                rows.append((grp, f"{s} {c['species']}", {k: 100 * d[k] / t for k, _, _ in BUCKET_SERIES}, f"{c['counts']:,.0f}" + (f"  {lab}" if lab else "")))
                tab.append((s, Raw(f"<i>{esc(c['species'])}</i>"), f"{c['counts']:,.0f}", f"{100 * c['abundance']:.2f}%", f"{share:.1f}%", st, lab))
        if any(r[0] for r in rows):
            rows.sort(key=lambda r: {"Expected": 0, "Not expected": 1}.get(r[0], 2))
        w(legend([(n, c) for _, n, c in BUCKET_SERIES]) + stacked(rows, BUCKET_SERIES, xmax=100, label_w=300, right_w=150 if check else 60, row_h=19,
                                                                  value_fmt=lambda v: f"{v:.1f}% of the call", tickfmt=lambda v: f"{v:.0f}%", title="Buckets behind each emu call"))
        w('<p class="cap">Each call\'s reads by bucket, as a share of the call\'s summed assignment probability. Right: emu\'s estimated read count'
          + (", and the call check\'s label." if check else ".")
          + (" Grouped by whether the sample sheet lists the species; a database name that ties with a listed species (such as an FDAARGOS strain name) shows as not expected." if any(expected.values()) else "")
          + "</p>")
        w(details("Table: emu calls", table(["Sample", "Species", "Estimated reads", "Abundance", "Chimeric", "Sample sheet", "Call check"], tab, numeric=(2, 3, 4))))

    # ---------------------------------------------------------------- quality
    if any(T["acc_sample"][s] for s in order):
        w("<h2>Quality</h2>")
        rows = []
        for s in ordered:
            h = T["acc_sample"][s]
            n = sum(h.values())
            if not n:
                continue
            q = sum(v for k, v in h.items() if k < rule.get("min_accuracy", 94))
            cum, med = 0, None
            for k in sorted(h):
                cum += h[k]
                if cum >= n / 2:
                    med = k
                    break
            rows.append((s, fmt(n), fmt(q), pct(q, n, 2), f"{med}%"))
        w(table(["Sample", "Reads in the length window", "Below the accuracy cut-off", "Share", "Median accuracy"], rows, numeric=(1, 2, 3, 4)))
        w(f'<p class="cap">Reads in the length window by mean per-base accuracy, as filtlong computes it; the cut-off is {rule.get("min_accuracy", "")}%. Median rounded down to a whole percent.</p>')

    # ---------------------------------------------------------------- methods
    w("<h2>Methods and limits</h2>")
    lim = [
        "Every number is counted from the run's own files: barbell's per-read tables, each sample's reads before and after host removal and filtlong, and emu's per-read assignments. Nothing is estimated.",
        "Barbell's patterns cover the raw reads; the read structure covers the reads barbell kept for a sample. Reads barbell refused or put in a barcode with no sample are counted, not aligned.",
    ]
    if structure:
        lim += [
            "References: the sample's savont ASVs, each named after the species of its closest database references (an emu-called species wins a tie), and for each species emu called (up to 50) "
            "database references chosen as those that most of the species' own reads (P ≥ 0.5) hit best: one if an ASV has the species, up to three if not. "
            "A species with neither can't be a parent, so its reads show up as poor fits or chimeras of its relatives.",
            "Positions are on E. coli 16S (K-12 MG1655 rrsH, from emu's database), lifted through each reference's alignment to it; they can be a few bases off in variable regions.",
            "Bucket order: unaligned, fold-back, same-strand repeat, strand switch, chimera, poor fit, unaligned tail or head, gapped, full length (spans the amplicon to within 100 bp at both ends), truncated.",
            "The chimera thresholds come from the 2026-07-02 run, where the consensuses of the run's pure cultures were the references for every sample. With these references instead, "
            "96.8–99.9% of each sample's reads in that run land in the same bucket, no pure-culture read is called chimeric (against 2–6 per culture), and four of the five mixes have between 11% fewer and 0.3% more chimeric reads. "
            "In AB13, most chimeras of its two klebsiellae are missed (379 against 1,198): a database variant of one species can fit a chimera of two close relatives almost as well as the pair does, so treat the rates as lower bounds. "
            "The same nine unexpected emu calls rest at least half on chimeric reads either way.",
            "The starting primer is inferred from the strand of the first aligned stretch, because barbell trims the primer.",
        ]
    lim.append("Reads MinKNOW failed (fastq_fail) are not part of the run's input and are not counted.")
    w("<ul>" + "".join(f"<li>{x}</li>" for x in lim) + "</ul>")
    w('<p class="cap">Tables: <code>read_buckets_summary.tsv</code> (per sample), <code>read_buckets_calls.tsv</code> (per emu call), and per-sample read tables if kept.</p>')

    page = ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>Read buckets</title><style>{CSS}</style></head><body><div class="rb">' + "".join(html) + f'</div><div id="rbtip"></div><script>{JS}</script></body></html>')
    with open("read_buckets.html", "w") as o:
        o.write(page)

    # ---------------------------------------------------------------- tables
    fate_keys = [k for k, _, _ in fate_series]
    with open("read_buckets_summary.tsv", "w") as o:
        cols = ["sample", "sample_type", "reads"] + fate_keys + [f"bucket_{b}" for b in BUCKET_ORDER] + [f"emu_bucket_{b}" for b in BUCKET_ORDER] + ["strand_fw", "strand_rc"]
        o.write("\t".join(cols) + "\n")
        for s in order:
            f = T["fate"][s]
            vals = [s, group[s], sum(f.values())] + [f[k] for k in fate_keys] + [T["bucket"][s][b] for b in BUCKET_ORDER] + [T["bucket_emu"][s][b] for b in BUCKET_ORDER] + \
                   [T["strand"][s]["+"], T["strand"][s]["-"]]
            o.write("\t".join(str(x) for x in vals) + "\n")
    with open("read_buckets_calls.tsv", "w") as o:
        o.write("\t".join(["sample", "tax_id", "species", "estimated_counts", "abundance", "sample_sheet", "call_check", "sum_p"] + [f"sum_p_{b}" for b in BUCKET_ORDER] +
                          [f"reads_p50_{b}" for b in BUCKET_ORDER]) + "\n")
        for s in order:
            for c in S[s]["calls"]:
                sp, rp = c["sum_p"], c["reads_p50"]
                o.write("\t".join(str(x) for x in [s, c["tax_id"], c["species"], f"{c['counts']:.2f}", f"{c['abundance']:.6f}", expected_status(c["species"], expected[s]),
                                                   check.get((s, c["tax_id"]), ""), f"{sum(sp.values()):.2f}"] + [f"{sp.get(b, 0):.2f}" for b in BUCKET_ORDER] +
                                [rp.get(b, 0) for b in BUCKET_ORDER]) + "\n")

    # ---------------------------------------------------------------- MultiQC custom content
    with open("read_buckets_fate_mqc.tsv", "w") as o:
        used_keys = [(k, n) for k, n, _ in fate_series if any(T["fate"][s][k] for s in order)]
        o.write("# id: 'read_buckets_fate'\n# section_name: 'Read fate'\n# description: 'Where the reads of each sample stop in the pipeline. From the read bucket report.'\n"
                "# plot_type: 'bargraph'\n# pconfig:\n#     id: 'read_buckets_fate_plot'\n#     title: 'Read buckets: fate of the reads of each sample'\n#     ylab: 'Reads'\n")
        o.write("Sample\t" + "\t".join(n for _, n in used_keys) + "\n")
        for s in order:
            o.write(s + "\t" + "\t".join(str(T["fate"][s][k]) for k, _ in used_keys) + "\n")
    if structure:
        with open("read_buckets_structure_mqc.tsv", "w") as o:
            cols = [b for b in BUCKET_ORDER if any(T["bucket_emu"][s][b] for s in order)]
            o.write("# id: 'read_buckets_structure'\n# section_name: 'Read structure'\n# description: 'What the reads that reach emu are made of. Chimera calls are uncalibrated. From the read bucket report.'\n"
                    "# plot_type: 'bargraph'\n# pconfig:\n#     id: 'read_buckets_structure_plot'\n#     title: 'Read buckets: structure of the reads emu sees'\n#     ylab: 'Reads'\n")
            o.write("Sample\t" + "\t".join(BUCKET_LABEL[b] for b in cols) + "\n")
            for s in order:
                o.write(s + "\t" + "\t".join(str(T["bucket_emu"][s][b]) for b in cols) + "\n")


if __name__ == "__main__":
    main()
