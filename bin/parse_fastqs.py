#!/usr/bin/env python3
import argparse, os, re, glob, csv
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a samplesheet CSV from FASTQ(s) in a dir or single files."
    )
    p.add_argument('--sample',      required=True)
    p.add_argument('--platform',    required=True)
    p.add_argument('--fastq1',      required=True,
                   help="Path to a FASTQ file or a directory to scan")
    p.add_argument('--fastq2',      default='',
                   help="Optional FASTQ2 file (only for single-sample mode)")
    p.add_argument('--seq-summary', default='',   dest='seq_summary')
    p.add_argument('--trim',        default='false')
    p.add_argument('--type',        default='UNKNOWN')
    p.add_argument('--output', '-o',required=True,
                   help="Path to write the CSV")
    return p.parse_args()

def scan_dir(dirpath):
    pats = ['*.fastq','*.fastq.gz','*.fq','*.fq.gz']
    files = []
    for pat in pats:
        files.extend(glob.glob(os.path.join(dirpath, pat)))
    return files

def group_fastqs(files):
    # capture: base name, optional _R1 or _2, then ext
    regex = re.compile(r'(.+?)(?:[_.]R?([12]))?\.(?:fastq|fq)(?:\.gz)?$', re.IGNORECASE)
    groups = defaultdict(dict)
    for fp in files:
        name = os.path.basename(fp)
        m = regex.match(name)
        if not m:
            continue
        base, rd = m.group(1), (m.group(2) or '1')
        groups[base][f'fastq_{rd}'] = os.path.abspath(fp)
    return groups

def main():
    args = parse_args()
    rows = []

    if os.path.isdir(args.fastq1):
        files  = scan_dir(args.fastq1)
        grouped = group_fastqs(files)
        for samp, mapping in grouped.items():
            rows.append({
                'sample'            : samp,
                'platform'          : args.platform,
                'fastq_1'           : mapping.get('fastq_1',''),
                'fastq_2'           : mapping.get('fastq_2',''),
                'sequencing_summary': args.seq_summary,
                'trim'              : args.trim,
                'type'              : args.type
            })
    else:
        # single-sample mode
        rows.append({
            'sample'            : args.sample,
            'platform'          : args.platform,
            'fastq_1'           : os.path.abspath(args.fastq1),
            'fastq_2'           : os.path.abspath(args.fastq2) if args.fastq2 else '',
            'sequencing_summary': args.seq_summary,
            'trim'              : args.trim,
            'type'              : args.type
        })

    # write out the CSV
    with open(args.output, 'w', newline='') as csvf:
        writer = csv.DictWriter(
            csvf,
            fieldnames=['sample','platform','fastq_1','fastq_2','sequencing_summary','trim','type']
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

if __name__ == '__main__':
    main()
