import os
import argparse
#initializes the parser so that you can start to add custom arguments
parser = argparse.ArgumentParser(description='Prepare a CSV file from a directory of sample data.')
parser.add_argument('directory', type=str, help='Path to directory containing sample data.')
args = parser.parse_args()

case_dir = args.directory

samples = []

for sample_dir in os.listdir(case_dir):
    sample_dir_path = os.path.join(case_dir, sample_dir)
    if not os.path.isdir(sample_dir_path):#check whether the specified path is an existing directory or not.
        continue

    fastq_files = sorted([f for f in os.listdir(sample_dir_path) if f.endswith('.fastq.gz')])
    print("",fastq_files)
    if len(fastq_files) == 0:
        continue

    sample = {
        'sample': sample_dir,
        'run_accession': sample_dir,
        'instrument_platform': 'ILLUMINA',
        'fastq_1': os.path.join(case_dir, sample_dir, fastq_files[0]),
        'fastq_2': os.path.join(case_dir, sample_dir, fastq_files[1]) if len(fastq_files) == 2 else '',
        'fasta': ''
    }

    samples.append(sample)

# Write samples to CSV file
csv_filename = os.path.join(case_dir, f"{os.path.basename(case_dir)}_samplesheet.csv")

with open(csv_filename, 'w') as f:
    f.write("sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta\n")
    for sample in samples:
        f.write(f"{sample['sample']},{sample['run_accession']},{sample['instrument_platform']},{sample['fastq_1']},{sample['fastq_2']},{sample['fasta']}\n")
print(f"CSV file '{csv_filename}' has been created.")
