#!/usr/bin/env python3
"""Check that every pathogen-sheet taxid downloaded an assembly with the expected name.

Left-joins the pathogen sheet (name, taxid) to the mapping written by
bin/append_taxid.py (Organism_Name, Mapped_Value = taxid) and writes, in sheet
order, one row per sheet entry:

    name    Organism_Name    ratio

ratio is Mapped_Value / taxid: 1.0 when the taxid matched, empty when the sheet
taxid has no downloaded assembly (Organism_Name then repeats the sheet name).
Rows where the two names differ are worth a look: the assembly picked for that
taxid may be a different strain/species label than the sheet expects.

Example:
    validate_names.py --mapsheet pathogens.mapping.taxid.tsv \\
        --pathogens_sheet assets/pathogen_sheet.csv -o pathogens.mapping.taxid.validated.tsv
"""
import argparse
import sys

import pandas as pd


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mapsheet", required=True,
                        help="Tab-delimited map with Organism_Name and Mapped_Value columns "
                             "(pathogens.mapping.taxid.tsv)")
    parser.add_argument("--pathogens_sheet", required=True,
                        help="Pathogen sheet CSV with name and taxid columns")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument("--mismatches", default=None,
                        help="Optional TSV of only the rows whose names differ or have no match")
    args = parser.parse_args(argv)

    try:
        maps = pd.read_csv(args.mapsheet, sep="\t")[["Organism_Name", "Mapped_Value"]]
    except Exception as e:  # noqa: BLE001
        sys.exit(f"ERROR reading mapsheet {args.mapsheet}: {e}")
    try:
        sheet = pd.read_csv(args.pathogens_sheet)[["name", "taxid"]]
    except Exception as e:  # noqa: BLE001
        sys.exit(f"ERROR reading pathogen sheet {args.pathogens_sheet}: {e}")

    maps["Mapped_Value"] = pd.to_numeric(maps["Mapped_Value"], errors="coerce")
    sheet["taxid"] = pd.to_numeric(sheet["taxid"], errors="coerce")
    sheet["ps_index"] = range(len(sheet))

    merged = sheet.merge(maps, left_on="taxid", right_on="Mapped_Value", how="left", sort=False)
    merged["ratio"] = merged["Mapped_Value"] / merged["taxid"]
    # several mapping rows (contigs) per taxid -> keep the first per sheet row
    out = (merged.sort_values("ps_index")
                 .groupby("ps_index", as_index=False).first()
                 [["name", "Organism_Name", "ratio"]])
    no_match = out["Organism_Name"].isna()
    out["Organism_Name"] = out["Organism_Name"].fillna(out["name"])
    out.to_csv(args.output, sep="\t", index=False)

    differ = out["name"] != out["Organism_Name"]
    if args.mismatches:
        out[differ | no_match].to_csv(args.mismatches, sep="\t", index=False)
    n = len(out)
    print(f"Wrote {args.output}: {n} sheet rows | no downloaded assembly: {int(no_match.sum())} | "
          f"name differs: {int(differ.sum())} ({differ.sum() / max(n, 1) * 100:.1f}%)")


if __name__ == "__main__":
    main()
