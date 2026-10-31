#!/usr/bin/env python
"""Turn an emu abundance report into the text format of ktImportText.

The lineage comes from the report, and so from the database the reads were classified
against, which means the plot can't lose taxa that a separate Krona taxonomy lacks.
"""

import argparse
import csv
from pathlib import Path

VERSION = "0.1.0"
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
COUNTS = "estimated counts"


def lineage(row):
    """Return the named ranks from top to bottom, or emu's label for reads without a taxon."""
    named = [row[rank] for rank in RANKS if row.get(rank)]
    return named or [row["tax_id"]]


def convert(report, output):
    with open(report, newline="") as fin, open(output, "w") as fout:
        for row in csv.DictReader(fin, delimiter="\t"):
            if float(row[COUNTS]) > 0:
                fout.write("\t".join([row[COUNTS], *lineage(row)]) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--version", action="version", version=f"emu_report_to_krona.py {VERSION}")
    args = parser.parse_args()
    convert(args.report, args.output)


if __name__ == "__main__":
    main()
