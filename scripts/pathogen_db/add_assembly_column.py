#!/usr/bin/env python3
"""Add an assembly_accession column to the pathogen sheet.

Maps each row's taxid to its downloaded assembly accession (GCF_*) using the
taxid->assembly mapping produced by append_taxid.py, then writes a new sheet
with the extra column. Rows whose taxid has no downloaded assembly get an empty
value in the new column.
"""
import argparse
import csv


def load_taxid_to_assembly(mapsheet):
    """Build {taxid: assembly_accession} from the mapping tsv.

    Expected columns: Acc, Assembly, Organism_Name, Description, Mapped_Value
    (Mapped_Value is the taxid).
    """
    mapping = {}
    with open(mapsheet, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            taxid = (row.get("Mapped_Value") or "").strip()
            assembly = (row.get("Assembly") or "").strip()
            if taxid and assembly and taxid not in mapping:
                mapping[taxid] = assembly
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pathogens_sheet",
        required=True,
        help="Path to the input pathogen sheet (CSV) with a 'taxid' column",
    )
    parser.add_argument(
        "--mapsheet",
        required=True,
        help="Path to the taxid->assembly mapping (tab-delimited) from append_taxid.py",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV path for the pathogen sheet with the assembly_accession column",
    )
    parser.add_argument(
        "--column",
        default="assembly_accession",
        help="Name of the new column (default: assembly_accession)",
    )
    parser.add_argument(
        "--position",
        choices=["after_taxid", "end"],
        default="after_taxid",
        help="Where to place the new column. Use 'end' for the shipped asset so "
             "scripts that read the sheet positionally (e.g. match_paths.py) are "
             "not shifted. Default: after_taxid.",
    )
    parser.add_argument(
        "--download_input",
        default=None,
        help="Optional path to also write a 2-column TSV (taxid<TAB>assembly_accession) "
             "for download_fastas.py --accession_col, kept in sync with the sheet.",
    )
    args = parser.parse_args()

    taxid_to_assembly = load_taxid_to_assembly(args.mapsheet)

    di = open(args.download_input, "w", newline="") if args.download_input else None

    matched = 0
    kept = 0
    total = 0
    with open(args.pathogens_sheet, newline="") as fin, \
            open(args.output, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        fieldnames = list(reader.fieldnames)
        if args.column not in fieldnames:
            if args.position == "after_taxid" and "taxid" in fieldnames:
                # Insert right after taxid for readability
                fieldnames.insert(fieldnames.index("taxid") + 1, args.column)
            else:
                # Append at the end (safest: no positional column shift)
                fieldnames.append(args.column)
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            total += 1
            taxid = (row.get("taxid") or "").strip()
            assembly = taxid_to_assembly.get(taxid, "")
            if assembly:
                matched += 1
            elif (row.get(args.column) or "").strip():
                # Not downloaded this run (NCBI hiccup, suppressed assembly...): keep
                # the curated accession so the next run still tries it first.
                assembly = row[args.column].strip()
                kept += 1
            row[args.column] = assembly
            writer.writerow(row)
            if di and taxid:
                di.write(f"{taxid}\t{assembly}\n")

    if di:
        di.close()

    print(f"Wrote '{args.output}' with column '{args.column}'.")
    if args.download_input:
        print(f"Wrote download input '{args.download_input}' (taxid<TAB>assembly_accession).")
    print(f"Rows: {total} | assembly from this build: {matched} | existing accession kept "
          f"(not downloaded this run): {kept} | without: {total - matched - kept}")


if __name__ == "__main__":
    main()
