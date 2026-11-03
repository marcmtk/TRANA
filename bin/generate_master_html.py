#!/usr/bin/env python

"""Generate a master html template."""

import os
import re
import argparse
import pandas as pd
from jinja2 import Template
from datetime import datetime

description = '''
------------------------
Title: generate_master_html.py
Date: 2024-12-16
Author(s): Ryan Kennedy
------------------------
Description:
    This script creates master html file that points to all html files that were outputted from EMU.

List of functions:
    get_sample_ids, generate_master_html.

List of standard modules:
    csv, os, argparse.

List of "non standard" modules:
    pandas, jinja2.

Procedure:
    1. Get sample IDs by parsing samplesheet csv.
    2. Render html using template.
    3. Write out master.html file.

-----------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------
Generates master html file that points to all html files.
Executed using: python3 ./generate_master_html.py -i <Input_Directory> -o <Output_Filepath>
-----------------------------------------------------------------------------------------------------------
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=usage
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.3.0'
    )
parser.add_argument(
    '-c', '--csv',
    help='input samplesheet csv filepath',
    metavar='SAMPLESHEET_CSV_FILEPATH',
    dest='csv',
    required=True
    )
parser.add_argument(
    '-m', '--html',
    help='input master html template filepath',
    metavar='MASTER_HTML_TEMPLATE_FILEPATH',
    dest='html',
    required=True
    )
parser.add_argument(
    '-t', '--timestamp',
    help='pipeline execution timestamp',
    metavar='PIPELINE_EXECUTION_TIMESTAMP',
    dest='timestamp',
    required=True
    )
parser.add_argument(
    '-s', '--sections',
    help='comma separated outputs this run produces, which the template links to',
    metavar='SECTIONS',
    dest='sections',
    default=''
    )
parser.add_argument(
    '-r', '--run-path',
    help='sequencing run path to read the run date from, instead of the samplesheet path',
    metavar='RUN_PATH',
    dest='run_path',
    default=''
    )
parser.add_argument(
    '--call-check',
    help='emu_call_check.tsv, whose labels are shown per sample',
    metavar='CALL_CHECK_TSV',
    dest='call_check',
    default=''
    )
parser.add_argument(
    '-o', '--output',
    help='output filepath',
    metavar='OUTPUT_FILEPATH',
    dest='output',
    required=True
    )

args = parser.parse_args()

def find_date_in_string(input_string, date_pattern):
    """Searches for a date within a given string."""
    date = "(No date found)"
    match = re.search(date_pattern, input_string)
    if match:
        date_matched = match.group(1)
        if len(date_matched) == 8:
            date = datetime.strptime(date_matched, "%Y%m%d").strftime("%d-%m-%Y")
        elif len(date_matched) > 8:
            date = date_matched
    return date

def get_sample_ids(samplesheet_csv):
    """Get sample id from csv."""
    df = pd.read_csv(samplesheet_csv)
    sample_ids = df['sample'].tolist()
    return sample_ids

def summarise_call_check(call_check_tsv, sample_ids):
    """Per sample, the emu calls by label from check_emu_calls.py, in order of estimated counts."""
    df = pd.read_csv(call_check_tsv, sep='\t', dtype={'sample': str})
    df = df.sort_values('estimated_counts', ascending=False)
    summary = []
    for sample_id in sample_ids:
        calls = df[df['sample'] == str(sample_id)]
        by_status = {
            status: [
                {'species': row.species, 'reads': int(row.reads_assigned), 'near_exact': int(row.near_exact_reads)}
                for row in calls[calls.status == status].itertuples()
            ]
            for status in ('confirmed', 'unconfirmed', 'not_assessed')
        }
        share = calls.main_species_near_exact_share.dropna()
        summary.append({
            'sample_id': sample_id,
            'checked': len(calls) > 0,
            'share': float(share.iloc[0]) if len(share) else None,
            **by_status,
        })
    return summary

def generate_master_html(template_html_fpath, sample_ids, seqrun_date, timestamp, sections, call_check):
    """Read the template from an HTML file."""
    with open(template_html_fpath, "r") as file:
        master_template = file.read()
    template = Template(master_template)
    rendered_html = template.render(sample_ids=sample_ids, seqrun_date=seqrun_date, timestamp=timestamp, sections=sections, call_check=call_check)
    return rendered_html

def main():
    sample_ids = get_sample_ids(args.csv)
    seqrun_date = find_date_in_string(args.run_path or args.csv, r'/(\d{8})_')
    sections = [section for section in args.sections.split(',') if section]
    call_check = summarise_call_check(args.call_check, sample_ids) if args.call_check else None
    rendered_html = generate_master_html(args.html, sample_ids, seqrun_date, args.timestamp, sections, call_check)
    with open(args.output, "w") as fout:
        fout.write(rendered_html)

if __name__ == "__main__":
    main()
