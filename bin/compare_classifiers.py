#!/usr/bin/env python
"""Put the emu and savont calls for one sample side by side.

savont splits an ASV between every reference tied at the lowest edit distance, so a
mislabelled reference only has to tie with the right one to produce a wrong-genus
call. savont_tied_genera lists the genera behind each species so those calls stand
out. E. coli/Shigella and other pairs that 16S can't separate will also show up there.

With --emu-check, emu_call_check gives the label from check_emu_calls.py for each
emu call: confirmed, unconfirmed or not_assessed.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

VERSION = "0.2.0"
COLUMNS = [
    "sample",
    "species",
    "genus",
    "emu_abundance",
    "emu_estimated_counts",
    "emu_call_check",
    "savont_abundance",
    "detected_by",
    "savont_tied_genera",
]


def read_tsv(path):
    with open(path, newline="") as fin:
        return list(csv.DictReader(fin, delimiter="\t"))


def tied_genera_per_species(mappings):
    """Return the genera tied at the lowest edit distance in any ASV holding the species."""
    lowest = {}
    for row in mappings:
        edits = int(row["number_mismatches"])
        lowest[row["asv_header"]] = min(edits, lowest.get(row["asv_header"], edits))

    asvs = defaultdict(list)
    for row in mappings:
        if int(row["number_mismatches"]) == lowest[row["asv_header"]]:
            asvs[row["asv_header"]].append(row)

    tied = defaultdict(set)
    for hits in asvs.values():
        genera = {hit["genus"] for hit in hits if hit["genus"]}
        for hit in hits:
            tied[hit["species"]] |= genera
    return tied


def compare(sample, emu_rows, savont_rows, mappings, check_rows=()):
    calls = {}

    def call(species, genus):
        return calls.setdefault(species, {"sample": sample, "species": species, "genus": genus})

    for row in emu_rows:
        # emu's unassigned and unclassified rows have no species
        if row.get("species"):
            entry = call(row["species"], row.get("genus", ""))
            entry["emu_abundance"] = row["abundance"]
            entry["emu_estimated_counts"] = row.get("estimated counts", "")

    for row in check_rows:
        if row["species"] in calls:
            calls[row["species"]]["emu_call_check"] = row["status"]

    for row in savont_rows:
        entry = call(row["species"], row.get("genus", ""))
        entry["savont_abundance"] = row["abundance"]

    tied = tied_genera_per_species(mappings)
    for species, entry in calls.items():
        emu = "emu_abundance" in entry
        savont = "savont_abundance" in entry
        entry["detected_by"] = "both" if emu and savont else "emu" if emu else "savont"
        genera = tied.get(species, set())
        entry["savont_tied_genera"] = ",".join(sorted(genera)) if savont and len(genera) > 1 else ""

    order = {"both": 0, "savont": 1, "emu": 2}
    return sorted(
        calls.values(),
        key=lambda c: (order[c["detected_by"]], -float(c.get("savont_abundance") or c.get("emu_abundance") or 0)),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--emu", type=Path, required=True)
    parser.add_argument("--savont-species", type=Path, required=True)
    parser.add_argument("--savont-mappings", type=Path, required=True)
    parser.add_argument("--emu-check", type=Path, help="table from check_emu_calls.py")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--version", action="version", version=f"compare_classifiers.py {VERSION}")
    args = parser.parse_args()

    check = read_tsv(args.emu_check) if args.emu_check else []
    rows = compare(args.sample, read_tsv(args.emu), read_tsv(args.savont_species), read_tsv(args.savont_mappings), check)
    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, COLUMNS, delimiter="\t", lineterminator="\n", restval="")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
