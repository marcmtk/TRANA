#!/usr/bin/env python
"""Label each emu species call confirmed or unconfirmed by its near-exact reads.

A chimeric read is a mosaic of two species and fits no single reference well, so a
call made from chimeras has reads that each align only roughly. A real species has
reads that match one of its references almost exactly. A read counts for a species if
emu assigns it there with probability >= --min-probability, and it is near-exact if its
best alignment to a reference of that species scores >= (1 - --max-error) * read length,
where score = matches - mismatches - inserted bases - deleted bases. A call is
confirmed with >= --min-reads near-exact reads.

Nothing is removed and emu's abundances are untouched. The share of near-exact reads
falls steeply with read quality, so when fewer than --min-share of the main species'
reads are near-exact, the calls are labelled not_assessed instead.

The thresholds were chosen on one run and haven't been checked on another with known
composition, so every row says calibration = uncalibrated until that is done.

See docs/evaluation/emu_chimera_step.md for the evaluation behind the defaults.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import pysam

VERSION = "0.3.0"
COLUMNS = [
    "sample",
    "tax_id",
    "species",
    "estimated_counts",
    "reads_assigned",
    "near_exact_reads",
    "status",
    "calibration",
    "main_species_near_exact_share",
]
CALIBRATION = "uncalibrated"
# pysam's CIGAR operation codes; = and X are aligned bases like M, and NM counts their mismatches
ALIGNED, INSERTED, DELETED, SOFT_CLIP, HARD_CLIP = {0, 7, 8}, 1, 2, 4, 5


def assigned_reads(path, min_probability):
    """Return each read's most likely tax id, if emu gives it at least min_probability."""
    assigned = {}
    with open(path, newline="") as fin:
        reader = csv.reader(fin, delimiter="\t")
        tax_ids = next(reader)[1:]
        for read, *probabilities in reader:
            best, best_p = None, 0.0
            for tax_id, p in zip(tax_ids, probabilities):
                if p and p != "NA" and float(p) > best_p:
                    best, best_p = tax_id, float(p)
            if best is not None and best_p >= min_probability:
                assigned[read] = best
    return assigned


def best_scores(path, assigned):
    """Return each assigned read's best alignment score and length against its own species."""
    score, length = {}, {}
    with pysam.AlignmentFile(str(path), check_sq=False) as alignments:
        for aln in alignments:
            if aln.is_unmapped:
                continue
            read = aln.query_name
            # emu's reference names start with the species tax id, e.g. 1280:emu_db:1
            if assigned.get(read) != aln.reference_name.split(":")[0]:
                continue
            ops = Counter()
            for op, n in aln.cigartuples:
                ops["M" if op in ALIGNED else op] += n
            mismatches = aln.get_tag("NM") - ops[INSERTED] - ops[DELETED]
            read_score = ops["M"] - 2 * mismatches - ops[INSERTED] - ops[DELETED]
            # secondary alignments are soft- or hard-clipped, so both count towards length
            length[read] = max(length.get(read, 0), ops["M"] + ops[INSERTED] + ops[SOFT_CLIP] + ops[HARD_CLIP])
            score[read] = max(score.get(read, read_score), read_score)
    return score, length


def check(sample, report, assigned, score, length, max_error, min_reads, min_share):
    n_assigned, n_near = Counter(), Counter()
    for read, tax_id in assigned.items():
        n_assigned[tax_id] += 1
        # a read scoring exactly (1 - max_error) * length counts; the tolerance keeps float rounding from dropping it
        if read in score and score[read] >= (1 - max_error) * length[read] - 1e-9:
            n_near[tax_id] += 1

    # emu's unassigned and unclassified rows have no tax id
    calls = [row for row in report if row["tax_id"].isdigit()]
    share = ""
    if calls:
        main = max(calls, key=lambda row: n_assigned[row["tax_id"]])["tax_id"]
        if n_assigned[main]:
            share = n_near[main] / n_assigned[main]

    rows = []
    for row in calls:
        near = n_near[row["tax_id"]]
        if share == "" or share < min_share:
            status = "not_assessed"
        else:
            status = "confirmed" if near >= min_reads else "unconfirmed"
        rows.append({
            "sample": sample,
            "tax_id": row["tax_id"],
            "species": row["species"],
            "estimated_counts": row.get("estimated counts", ""),
            "reads_assigned": n_assigned[row["tax_id"]],
            "near_exact_reads": near,
            "status": status,
            "calibration": CALIBRATION,
            "main_species_near_exact_share": f"{share:.3f}" if share != "" else "",
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--report", type=Path, required=True, help="emu rel-abundance table")
    parser.add_argument("--assignments", type=Path, required=True, help="emu read-assignment-distributions table")
    parser.add_argument("--alignments", type=Path, required=True, help="emu's alignments as BAM or SAM, kept with --keep-files")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-error", type=float, default=0.01)
    parser.add_argument("--min-reads", type=int, default=2)
    parser.add_argument("--min-probability", type=float, default=0.5)
    parser.add_argument("--min-share", type=float, default=0.2)
    parser.add_argument("--version", action="version", version=f"check_emu_calls.py {VERSION}")
    args = parser.parse_args()

    with open(args.report, newline="") as fin:
        report = list(csv.DictReader(fin, delimiter="\t"))
    assigned = assigned_reads(args.assignments, args.min_probability)
    score, length = best_scores(args.alignments, assigned)
    rows = check(args.sample, report, assigned, score, length, args.max_error, args.min_reads, args.min_share)
    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
