#!/usr/bin/env python3
"""PROTOTYPE, not wired into TRANA. Post-emu evidence check: does each emu species call have reads that match one of its own references
almost exactly? Chimera-driven calls do not: their reads are mosaics that fit no single reference well.

Inputs are emu's own outputs (run with --keep-files --keep-read-assignments):
  SAM, <sample>_read-assignment-distributions.tsv, <sample>_rel-abundance.tsv
A read counts for a species if emu gives it P >= --min-p for that species. It is near-exact if its best alignment to any reference of that species
scores >= (1 - --max-err) * read length, where score = matching bases - mismatches - inserted bases - deleted bases.
A call is 'confirmed' with >= --min-reads near-exact reads, otherwise 'unconfirmed'. Nothing is removed; the call is labelled.

usage: prototype_near_exact.py SAM ASSIGN ABUND [--max-err 0.01] [--min-reads 2] [--min-p 0.5] > table.tsv"""
import sys, re, argparse, csv
ap = argparse.ArgumentParser()
ap.add_argument('sam'); ap.add_argument('assign'); ap.add_argument('abund')
ap.add_argument('--max-err', type=float, default=0.01); ap.add_argument('--min-reads', type=int, default=2); ap.add_argument('--min-p', type=float, default=0.5)
a = ap.parse_args()
# read -> species it is assigned to
top = {}
with open(a.assign) as h:
    r = csv.reader(h, delimiter='\t'); head = next(r)
    tids = head[1:]
    for row in r:
        best = None; bp = 0.0
        for t, v in zip(tids, row[1:]):
            if v and v != 'NA' and float(v) > bp: bp = float(v); best = t
        if best is not None and bp >= a.min_p: top[row[0]] = best
cig = re.compile(r'(\d+)([MIDNSHP=X])')
best = {}; qlen = {}
with open(a.sam) as h:
    for l in h:
        if l[0] == '@': continue
        f = l.split('\t', 11)
        q = f[0]; t = f[2].split(':')[0]
        if f[2] == '*' or top.get(q) != t: continue
        c = {'M': 0, 'I': 0, 'D': 0, 'S': 0, 'H': 0}
        for n, o in cig.findall(f[5]):
            if o in c: c[o] += int(n)
        nm = int(re.search(r'NM:i:(\d+)', f[11]).group(1))
        x = nm - c['I'] - c['D']
        sc = c['M'] - 2 * x - c['I'] - c['D']
        ql = c['M'] + c['I'] + c['S'] + c['H']
        qlen[q] = max(qlen.get(q, 0), ql)
        if sc > best.get(q, -1e9): best[q] = sc
n_ass = {}; n_near = {}
for q, t in top.items():
    n_ass[t] = n_ass.get(t, 0) + 1
    if q in best and best[q] >= (1 - a.max_err) * qlen[q]: n_near[t] = n_near.get(t, 0) + 1
w = csv.writer(sys.stdout, delimiter='\t')
w.writerow(['tax_id', 'species', 'estimated_counts', 'reads_assigned', 'near_exact_reads', 'status'])
with open(a.abund) as h:
    for row in csv.DictReader(h, delimiter='\t'):
        t = row['tax_id']
        if not t.isdigit(): continue
        k = n_near.get(t, 0)
        w.writerow([t, row['species'], row.get('estimated counts', ''), n_ass.get(t, 0), k, 'confirmed' if k >= a.min_reads else 'unconfirmed'])
