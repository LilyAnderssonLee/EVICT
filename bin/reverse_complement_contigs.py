#!/usr/bin/env python3

import argparse
import csv

def parse_args():
    parser = argparse.ArgumentParser(
        description="Reverse-complement contigs whose BLAST hits are in reverse orientation (sstart > send)."
    )
    parser.add_argument("--blast", required=True, help="Comma-delimited BLAST file.")
    parser.add_argument("--fasta", required=True, help="Input FASTA file.")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file with corrected orientations.")
    return parser.parse_args()


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def get_reverse_contigs(blast_file):
    contigs = set()
    with open(blast_file, newline="") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            try:
                sstart = int(row["sstart"])
                send = int(row["send"])
            except:
                continue
            if sstart > send:
                contigs.add(row["qseqid"])
    return contigs

def process_fasta(fasta_in, fasta_out, reverse_set):
    """Reverse-complement sequences whose ID is in reverse_set and update header."""
    def write_record(out, header, seq_lines):
        if not header:
            return
        seq = "".join(seq_lines)

        # Extract FASTA ID (everything after ">" up to first space)
        fasta_id = header.split()[0][1:]

        if fasta_id in reverse_set:
            seq = reverse_complement(seq)
            new_header = f">{fasta_id}_reverse_complement"
        else:
            new_header = header

        out.write(new_header + "\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        header = None
        seq_lines = []

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                write_record(fout, header, seq_lines)
                header = line
                seq_lines = []
            else:
                if line:
                    seq_lines.append(line.strip())

        write_record(fout, header, seq_lines)


def main():
    args = parse_args()

    reverse_set = get_reverse_contigs(args.blast)
    if reverse_set:
        print(f"Reverse-complementing {len(reverse_set)} contig(s): {', '.join(reverse_set)}")
    else:
        print("No reverse-oriented contigs detected.")

    process_fasta(args.fasta, args.output_fasta, reverse_set)

    print(f"Output FASTA written to: {args.output_fasta}")


if __name__ == "__main__":
    main()
