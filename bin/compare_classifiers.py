#!/usr/bin/env python
"""Put the emu and savont calls for one sample side by side.

savont names each ASV after one reference, even when several references are equally
close. Its asv_tied_classification.tsv names all of them, so for every species in the
table the savont_tie_* columns say whether an ASV behind it is tied, and with what:
savont_tie_level is none, species (several species of one genus) or cross_genus,
savont_tied_species gives the tied names, and savont_tied_genera lists the genera when
there is more than one. A species emu calls and savont only has inside a tie, such as
S. agalactiae in an ASV savont names after Streptococcus sp. FDAARGOS_520, is detected_by
emu but has the tie columns filled. E. coli/Shigella and other pairs that 16S can't
separate show up as ties too, and so does a mislabelled reference that ties with the
right one.

With --emu-check, emu_call_check gives the label from check_emu_calls.py for each
emu call: confirmed, unconfirmed or not_assessed.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

VERSION = "0.3.0"
COLUMNS = [
    "sample",
    "species",
    "genus",
    "emu_abundance",
    "emu_estimated_counts",
    "emu_call_check",
    "savont_abundance",
    "detected_by",
    "savont_tie_level",
    "savont_tied_species",
    "savont_tied_genera",
]
TIE_SEPARATOR = " / "
TIE_LEVELS = ["none", "species", "cross_genus"]


def read_tsv(path):
    with open(path, newline="") as fin:
        return list(csv.DictReader(fin, delimiter="\t"))


def ties_per_species(mappings, tied):
    """Return, per species, the tie rows of the ASVs behind it.

    An ASV is behind a species if savont named it after the species, or if the species is
    one of the names tied in it.
    """
    by_asv = {row["asv_header"]: row for row in tied if row["species"] != "UNCLASSIFIED"}
    behind = defaultdict(dict)
    for row in mappings:
        if row["asv_header"] in by_asv:
            behind[row["species"]][row["asv_header"]] = by_asv[row["asv_header"]]
    for asv, row in by_asv.items():
        for species in row["species"].split(TIE_SEPARATOR):
            behind[species][asv] = row
    return behind


def tie_columns(rows):
    """Summarise the tie rows of the ASVs behind one species."""
    if not rows:
        return {}
    level = max((row["tie_level"] for row in rows), key=TIE_LEVELS.index)
    names = sorted({row["species"] for row in rows if row["tie_level"] != "none"})
    genera = sorted({genus for row in rows for genus in row["genus"].split(TIE_SEPARATOR) if genus})
    return {
        "savont_tie_level": level,
        "savont_tied_species": "; ".join(names),
        "savont_tied_genera": ",".join(genera) if len(genera) > 1 else "",
    }


def compare(sample, emu_rows, savont_rows, mappings, tied, check_rows=()):
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

    behind = ties_per_species(mappings, tied)
    for species, entry in calls.items():
        emu = "emu_abundance" in entry
        savont = "savont_abundance" in entry
        entry["detected_by"] = "both" if emu and savont else "emu" if emu else "savont"
        entry.update(tie_columns(list(behind.get(species, {}).values())))

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
    parser.add_argument("--savont-tied", type=Path, required=True, help="asv_tied_classification.tsv from savont classify")
    parser.add_argument("--emu-check", type=Path, help="table from check_emu_calls.py")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--version", action="version", version=f"compare_classifiers.py {VERSION}")
    args = parser.parse_args()

    check = read_tsv(args.emu_check) if args.emu_check else []
    rows = compare(args.sample, read_tsv(args.emu), read_tsv(args.savont_species), read_tsv(args.savont_mappings), read_tsv(args.savont_tied), check)
    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, COLUMNS, delimiter="\t", lineterminator="\n", restval="")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
