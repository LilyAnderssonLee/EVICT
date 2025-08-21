import pandas as pd
import sys

def summarize_blast(file):
    columns = [
        "qseqid", "sseqid", "evalue","bitscore", "pident", "qlen",
        "qstart", "qend", "sstart", "send", "taxid",
        "scomname", "length"
    ]

    blast_df = pd.read_csv(file, sep=",", names=columns, header=None)

    # Convert numeric columns to proper data types
    numeric_columns = ['evalue', 'bitscore', 'pident', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'taxid', 'length']

    for col in numeric_columns:
        blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce')

    # Remove rows where conversion failed (NaN values)
    blast_df = blast_df.dropna(subset=numeric_columns)

    summary = (
        blast_df.groupby(['qseqid', 'taxid', 'scomname'])
        .agg(
            count=('taxid', 'count'),
            min_pident=('pident', 'min'),
            max_pident=('pident', 'max'),
            median_pident=('pident', 'median'),
            min_length=('length', 'min'),
            max_length=('length', 'max'),
            median_length=('length', 'median'),
            min_bitscore=('bitscore', 'min'),
            max_bitscore=('bitscore', 'max'),
            median_bitscore=('bitscore', 'median')
        )
        .reset_index()
    )

    output_file = f"{file}.summary.csv"
    summary.to_csv(output_file, index=False)
    print(f"Summary of blast saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python blast_summary.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    summarize_blast(input_file)