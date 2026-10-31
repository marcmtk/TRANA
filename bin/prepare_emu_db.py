#!/usr/bin/env python
"""Unpack an emu database into the one form that both emu and savont can read.

The 2026 emu database trims trailing empty fields from taxonomy.tsv. emu copes, but
savont skips any row with fewer than 12 fields, which would silently lose most of
the species, so the rows are padded here. savont also needs a marker file to know the
layout, since the folder isn't named after a database it recognises.
"""

import argparse
import gzip
import hashlib
import sys
import tarfile
import tempfile
from pathlib import Path

VERSION = "0.1.0"
FASTA_NAME = "species_taxid.fasta"
TAXONOMY_NAME = "taxonomy.tsv"
SAVONT_MARKER = ".savont_db"
SAVONT_DB_NAME = "emu-1"


def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def find_file(root, name):
    """Find name below root, preferring an unpacked copy over a gzipped one."""
    for candidate in (name, f"{name}.gz"):
        hits = sorted(root.rglob(candidate))
        if hits:
            return hits[0]
    sys.exit(f"ERROR: no {name} or {name}.gz found in the database at {root}")


def write_taxonomy(source, dest):
    """Copy taxonomy.tsv with every row padded to the width of the header."""
    rows = padded = 0
    with open_text(source) as fin, open(dest, "w") as fout:
        width = len(fin.readline().rstrip("\r\n").split("\t"))
        fin.seek(0)
        for line_no, line in enumerate(fin, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) > width:
                sys.exit(f"ERROR: {source} line {line_no} has {len(fields)} fields, the header has {width}")
            padded += len(fields) < width
            rows += 1
            fout.write("\t".join(fields + [""] * (width - len(fields))) + "\n")
    return rows - 1, padded


def write_fasta(source, dest, exclude):
    """Copy the FASTA without the listed records, and return how many were kept and dropped."""
    kept = 0
    dropped = set()
    keep = True
    with open_text(source) as fin, open(dest, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                record_id = line[1:].split(None, 1)[0]
                keep = record_id not in exclude
                if keep:
                    kept += 1
                else:
                    dropped.add(record_id)
            if keep:
                fout.write(line)
    return kept, dropped


def read_ids(path):
    with open(path) as fin:
        return {line.strip() for line in fin if line.strip() and not line.startswith("#")}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fin:
        for chunk in iter(lambda: fin.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare(db_root, outdir, exclude_path, summary_path):
    outdir.mkdir(parents=True, exist_ok=True)
    exclude = read_ids(exclude_path) if exclude_path else set()

    source_fasta = find_file(db_root, FASTA_NAME)
    n_taxa, n_padded = write_taxonomy(find_file(db_root, TAXONOMY_NAME), outdir / TAXONOMY_NAME)
    n_kept, dropped = write_fasta(source_fasta, outdir / FASTA_NAME, exclude)

    # A list made for another database release must not pass unnoticed
    missing = exclude - dropped
    if missing:
        sys.exit(
            f"ERROR: {len(missing)} of {len(exclude)} IDs in {exclude_path} are not in the database, "
            f"e.g. {sorted(missing)[0]}. Was the list made for another database release?"
        )

    (outdir / SAVONT_MARKER).write_text(SAVONT_DB_NAME)

    with open(summary_path, "w") as fout:
        fout.write(f"reference_sequences\t{n_kept}\n")
        fout.write(f"excluded_sequences\t{len(dropped)}\n")
        fout.write(f"exclusion_list\t{exclude_path.name if exclude_path else ''}\n")
        fout.write(f"taxa\t{n_taxa}\n")
        fout.write(f"taxonomy_rows_padded\t{n_padded}\n")
        fout.write(f"species_taxid_fasta_sha256\t{sha256(outdir / FASTA_NAME)}\n")
        fout.write(f"taxonomy_tsv_sha256\t{sha256(outdir / TAXONOMY_NAME)}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True, help="Folder or .tar.gz holding the database")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--exclude-ids", type=Path, help="FASTA header IDs to drop, one per line")
    parser.add_argument("--version", action="version", version=f"prepare_emu_db.py {VERSION}")
    args = parser.parse_args()

    if args.db.is_file():
        with tempfile.TemporaryDirectory(dir=".") as tmp, tarfile.open(args.db) as tar:
            tar.extractall(tmp, filter="data")
            prepare(Path(tmp), args.outdir, args.exclude_ids, args.summary)
    else:
        prepare(args.db, args.outdir, args.exclude_ids, args.summary)


if __name__ == "__main__":
    main()
