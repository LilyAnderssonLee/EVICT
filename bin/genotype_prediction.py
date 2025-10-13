#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import argparse
import os
import sys

def suggest_genotype(df, suggest_min_rows=20,
                     suggest_min_identity=90.0, suggest_min_bitscore=400.0):
    """
    Genotype suggestion logic (identical to html_report.py)
    """
    try:
        max_pident = df.loc[df['pident'].idxmax()]
        max_bitscore = df.loc[df['bitscore'].idxmax()]
        grouped_pident_medians = df.groupby('scomname')['pident'].median()
        highest_median_pident = grouped_pident_medians.idxmax()
        grouped_bitscore_medians = df.groupby('scomname')['bitscore'].median()
        highest_median_bitscore = grouped_bitscore_medians.idxmax()

        # Count number of hits for the top species
        species_counts = df['scomname'].value_counts()
        top_species = max_pident['scomname']
        top_species_count = species_counts.get(top_species, 0)

        # Apply suggestion criteria
        if (
            (max_pident['scomname'] == max_bitscore['scomname']) and
            (top_species_count >= suggest_min_rows) and
            (df['pident'].max() >= suggest_min_identity) and
            (df['bitscore'].max() >= suggest_min_bitscore) and
            (highest_median_pident == max_pident['scomname']) and
            (highest_median_bitscore == max_bitscore['scomname'])
        ):
            suggestion = str(max_pident['scomname'])
        else:
            suggestion = "Var god bedöm manuellt."
    except Exception as e:
        print(f"⚠️ Suggestion logic failed: {e}")
        suggestion = "Var god bedöm manuellt."

    return suggestion


def main():
    parser = argparse.ArgumentParser(description="Summarize genotype suggestion from BLAST file")
    parser.add_argument("--ticket", type=str, required=True, help="Ticket")
    parser.add_argument("--blast-file", type=str, required=True, help="Path to BLAST result file (.blast)")
    parser.add_argument("--output", type=str, default="genotype.csv", help="Output CSV path")
    parser.add_argument("--suggest-min-rows", type=int, default=20,
                        help="Minimum number of hits (rows) required to consider auto-suggestion (default: 20)")
    parser.add_argument("--suggest-min-identity", type=float, default=90.0,
                        help="Minimum max %% identity required to consider auto-suggestion (default: 90)")
    parser.add_argument("--suggest-min-bitscore", type=float, default=400.0,
                        help="Minimum max bitscore required to consider auto-suggestion (default: 400)")
    args = parser.parse_args()

    if not os.path.exists(args.blast_file):
        print(f"❌ Error: BLAST file not found: {args.blast_file}")
        sys.exit(1)

    file_name = os.path.basename(args.blast_file)
    sample_name = file_name.replace(".blast", "")

    # Read BLAST file
    df = pd.read_csv(args.blast_file)

    # --- Filtering (same as html_report.py) ---
    df[["contig", "temp1"]] = df["qseqid"].str.split("_length_", expand=True)
    df[["length", "coverage"]] = df["temp1"].str.split("_cov_", expand=True)
    df.drop(["temp1", "length"], axis=1, inplace=True)
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    df = df.loc[df["qlen"] > 200]
    df = df.loc[df["coverage"] > 50]
    # ------------------------------------------

    if df.empty or df["scomname"].nunique() == 0:
        genotype = "Ingen giltig contig (för kort/låg täckning)"
    else:
        genotype = suggest_genotype(
            df,
            suggest_min_rows=args.suggest_min_rows,
            suggest_min_identity=args.suggest_min_identity,
            suggest_min_bitscore=args.suggest_min_bitscore
        )

    out_df = pd.DataFrame([{
        "Ticket": args.ticket,
        "Sample": sample_name,
        "Genotype": genotype
    }])

    # Append or create CSV
    header = not os.path.exists(args.output)
    out_df.to_csv(args.output, index=False, mode="a", header=header, encoding="utf-8-sig")

    print(f"✅ Added: Ticket {args.ticket}, Sample {sample_name}, Genotype: {genotype}")
    print(f"📄 Saved/updated summary: {args.output}")


if __name__ == "__main__":
    main()
