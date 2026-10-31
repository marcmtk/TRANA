#!/usr/bin/env bash
# Name barbell's barcode bins after their samples and count the reads in every bin.
#
# Usage: barbell_bins_to_samples.sh <input_file_list> <barbell_dir> <out_dir> [barcode,sample_id CSV]
#
# Barcodes are matched on their number, so a sample sheet written for MinKNOW's
# barcode01 also finds barbell's AB01. Bins not in the sheet aren't analysed, but
# stay in the counts: an unused barcode is a blank well.
set -euo pipefail

input_files=$1
barbell_dir=$2
out_dir=$3
sheet=${4:-}

mkdir -p "$out_dir"

barcode_number() {
    local digits
    digits=$(grep -oE '[0-9]+$' <<< "$1" || true)
    [[ -n "$digits" ]] && echo $((10#$digits))
}

count_reads() {
    echo $(( $(gzip -cd "$1" | wc -l) / 4 ))
}

declare -A sample_of
if [[ -n "$sheet" ]]; then
    while IFS=, read -r barcode sample _ || [[ -n "$barcode" ]]; do
        sample=${sample//$'\r'/}
        [[ -z "$barcode" || "$barcode" == "barcode" ]] && continue
        number=$(barcode_number "$barcode")
        if [[ -z "$number" ]]; then
            echo "Can't read a barcode number from '$barcode' in $sheet" >&2
            exit 1
        fi
        if [[ -n "${sample_of[$number]:-}" ]]; then
            echo "Barcode $barcode is listed twice in $sheet" >&2
            exit 1
        fi
        sample_of[$number]=$sample
    done < "$sheet"
fi

declare -A seen
printf 'barcode\tsample_id\treads\n' > barbell_read_counts.tsv
for bin in "$barbell_dir"/*.trimmed.fastq.gz; do
    [[ -e "$bin" ]] || continue
    barcode=$(basename "$bin" .trimmed.fastq.gz)
    number=$(barcode_number "$barcode")
    reads=$(count_reads "$bin")
    sample=
    if [[ -z "$sheet" ]]; then
        sample=$barcode
    elif [[ -n "$number" ]]; then
        sample=${sample_of[$number]:-}
    fi
    printf '%s\t%s\t%s\n' "$barcode" "$sample" "$reads" >> barbell_read_counts.tsv
    if [[ -n "$sample" ]]; then
        seen[$number]=1
        mv "$bin" "$out_dir/$sample.fastq.gz"
    fi
done

input_reads=0
while IFS= read -r file; do
    input_reads=$(( input_reads + $(count_reads "$file") ))
done < "$input_files"
printf 'run_input\t\t%s\n' "$input_reads" >> barbell_read_counts.tsv

: > samples_without_reads.txt
for number in "${!sample_of[@]}"; do
    if [[ -z "${seen[$number]:-}" ]]; then
        echo "${sample_of[$number]}" >> samples_without_reads.txt
    fi
done

if [[ -z "$(ls -A "$out_dir")" ]]; then
    echo "barbell assigned no reads to any sample. Were the barcodes already trimmed, or is this another kit?" >&2
    exit 1
fi
