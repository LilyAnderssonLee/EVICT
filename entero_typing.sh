#!/usr/bin/env bash
#SBATCH -A development
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH -J EV_genotyping
#SBATCH --qos normal
#SBATCH --mem 180G

# =====================
# Config & Environment
# =====================
set -eo pipefail

# ---- Usage & args ----
usage() {
    cat <<EOF
Usage: $(basename "$0") <ticket_id> <sample_id>

Runs EV genotyping on a given ticket and sample.

Arguments:
    ticket_id   Folder name under production inbox and used for all outputs.
    sample_id   Sample folder name inside the ticket (used for SPAdes/BLAST/report steps).

Environment:
    Requires: conda env D_EV_LAL, nextflow, seqkit, spades.py, blastn, python.

Outputs (under BASE_DIR):
    taxprofiler_results/<ticket>/
    results/<ticket>/{spades/<sample>/, blast/, ev_contig/}
    data/<ticket>/<ticket>_samplesheet.csv   (generated)
EOF
}

# Check if user asked for help (-h/--help) OR did not provide exactly 2 arguments
if [[ "${1:-}" =~ ^(-h|--help)$ || "$#" -ne 2 ]]; then
    usage; [[ "$#" -ne 2 ]] && exit 1 || exit 0
fi

ticket=$1
sample=$2

# ---- Conda env ----
source /home/proj/stage/bin/miniconda3/etc/profile.d/conda.sh
conda activate D_EV_LAL

set -u

# ---- Paths (single source of truth) ----
BASE_DIR="/home/proj/development/microbial/metagenomics/enterovirus"
ASSETS_DIR="$BASE_DIR/assets"
BIN_DIR="$BASE_DIR/bin"
DATA_DIR="$BASE_DIR/data"
RESULTS_DIR="$BASE_DIR/results"
TAXP_DIR="$BASE_DIR/taxprofiler_results"

CONFIG="$ASSETS_DIR/custom_taxprofiler.config"
PARAMS="$ASSETS_DIR/params.json"
PRIMERS="$ASSETS_DIR/primers_fw_rev_compliment.fasta"
DB_DIR="/home/proj/development/microbial/metagenomics/databases/archive_latest_databases/blastdb/nt_viruses"
BLAST_DB="$DB_DIR/nt_viruses"
REF="/home/proj/development/microbial/metagenomics/references/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
PROD_INBOX="/home/proj/production/customers/cust160/inbox"

# ---- Resources (fallback if not in Slurm) ----
THREADS="${SLURM_NTASKS:-4}"
MEM_MB="${SLURM_MEM_PER_NODE:-16384}"

# ---- Logging & error trap ----
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
trap 'rc=$?; echo "[ERROR] $(date "+%F %T") step failed (exit $rc)"; exit $rc' ERR

# ---- Dependency check ----
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 127; }; }
for c in nextflow python seqkit spades.py blastn awk sort cut grep find; do need_cmd "$c"; done

# Create core output dirs up front (idempotent)
mkdir -p "$DATA_DIR" "$RESULTS_DIR/$ticket"/{spades,blast,ev_contig} "$TAXP_DIR"

# =====================
# Functions
# =====================

prepare_input() {
    log "Preparing input for taxprofiler: ticket=$ticket"

    local ticket_src="$PROD_INBOX/$ticket"
    local ticket_dst="$DATA_DIR/$ticket"

    if [[ ! -d "$ticket_dst" ]]; then
        [[ -d "$ticket_src" ]] || { echo "Source ticket folder not found: $ticket_src" >&2; exit 2; }
        log "Copying $ticket_src -> $DATA_DIR"
        cp -r "$ticket_src" "$DATA_DIR"
    else
        log "Input already present at $ticket_dst"
    fi

    # Merge FASTQs per sample (idempotent)
    # Expecting: data/<ticket>/<sample>/*_R1_*fastq.gz and *_R2_*fastq.gz
    log "Merging FASTQs under $ticket_dst"
    shopt -s nullglob
    for sdir in "$ticket_dst"/*/; do
        sdir="${sdir%/}"
        local sname; sname="$(basename "$sdir")"

        # Skip if already merged
        if [[ -s "$sdir/${sname}_1.fastq.gz" && -s "$sdir/${sname}_2.fastq.gz" ]]; then
            log "  $sname: merged files already exist, skipping"
            continue
        fi

        # Find R1/R2 files robustly (group -name terms so -maxdepth/-type apply to both)
        mapfile -t r1s < <(find "$sdir" -maxdepth 1 -type f \
            \( -name '*_R1_*fastq.gz' -o -name '*_R1_*.fastq.gz' \) | sort)
        mapfile -t r2s < <(find "$sdir" -maxdepth 1 -type f \
            \( -name '*_R2_*fastq.gz' -o -name '*_R2_*.fastq.gz' \) | sort)

        if (( ${#r1s[@]} == 0 || ${#r2s[@]} == 0 )); then
            log "  $sname: no R1/R2 pairs found, skipping"
            continue
        fi
        log "  $sname: merging ${#r1s[@]} R1 files and ${#r2s[@]} R2 files"

        # Concatenate deterministically
        cat "${r1s[@]}" > "$sdir/${sname}_1.fastq.gz"
        cat "${r2s[@]}" > "$sdir/${sname}_2.fastq.gz"

        # Remove the original lane-split FASTQs (keep the merged ones)
        for f in "${r1s[@]}" "${r2s[@]}"; do
            [[ "$f" == "$sdir/${sname}_1.fastq.gz" || "$f" == "$sdir/${sname}_2.fastq.gz" ]] && continue
            rm -f -- "$f" || true
        done
    done
    shopt -u nullglob

    # Generate samplesheet for taxprofiler
    local samplesheet="$ticket_dst/${ticket}_samplesheet.csv"
    if [[ -s "$samplesheet" ]]; then
        log "Samplesheet already exists: $samplesheet"
    else
        log "Generating samplesheet: $samplesheet"
        python "$BIN_DIR/taxprofiler_samplesheet.py" "$ticket_dst"
    fi
}

run_taxprofiler() {
    log "Running nf-core/taxprofiler for ticket=$ticket"
    local tp_out="$TAXP_DIR/$ticket"

    if [[ -d "$tp_out" ]]; then
        log "Taxprofiler results already exist: $tp_out"
        return 0
    fi

    nextflow run nf-core/taxprofiler -r 1.2.3 -profile hasta,singularity \
        --input "$DATA_DIR/${ticket}/${ticket}_samplesheet.csv" \
        --databases "$ASSETS_DIR/databases.csv" \
        --outdir "$TAXP_DIR/${ticket}" \
        --save_preprocessed_reads --perform_shortread_qc \
        --perform_shortread_complexityfilter --save_complexityfiltered_reads \
        --perform_shortread_hostremoval --hostremoval_reference "$REF" \
        --save_hostremoval_index --save_hostremoval_bam --save_hostremoval_unmapped \
        --shortread_qc_adapterlist "$PRIMERS" \
        -params-file "$PARAMS" -c "$CONFIG" -resume
}

run_spades() {
    log "Running SPAdes for sample=$sample (ticket=$ticket)"
    local tp_unmap_dir="$TAXP_DIR/${ticket}/bowtie2/align"
    local read1="${tp_unmap_dir}/${sample}_${sample}.unmapped_1.fastq.gz"
    local read2="${tp_unmap_dir}/${sample}_${sample}.unmapped_2.fastq.gz"
    local outdir="$RESULTS_DIR/$ticket/spades/$sample"

    mkdir -p "$outdir"

    if [[ -s "$outdir/contigs.fasta" || -s "$outdir/scaffolds.fasta" ]]; then
        log "SPAdes outputs found in $outdir, skipping"
        return 0
    fi

    [[ -s "$read1" && -s "$read2" ]] || { echo "Missing unmapped reads: $read1 / $read2" >&2; exit 3; }

    spades.py -t "$THREADS" -m $(( MEM_MB / 1024 )) \
        -1 "$read1" -2 "$read2" \
        -o "$outdir" --rnaviral
}

run_blast() {
    log "Running BLAST for sample=$sample (ticket=$ticket)"
    local outdir="$RESULTS_DIR/$ticket/blast"
    mkdir -p "$outdir"

    local contigs_scaf="$RESULTS_DIR/$ticket/spades/$sample/scaffolds.fasta"
    local contigs_std="$RESULTS_DIR/$ticket/spades/$sample/contigs.fasta"
    local query
    if [[ -s "$contigs_scaf" ]]; then
        query="$contigs_scaf"
    elif [[ -s "$contigs_std" ]]; then
        query="$contigs_std"
    else
        echo "No contigs found for BLAST in $RESULTS_DIR/$ticket/spades/$sample" >&2
        return 0
    fi

    local final="$outdir/${sample}.blast"
    if [[ -s "$final" ]]; then
        log "BLAST already completed: $final"
        return 0
    fi

    local taxid_list="$ASSETS_DIR/taxid_EV.txt"
    local header="$ASSETS_DIR/blast_header.txt"
    local tmp
    tmp="$(mktemp "$outdir/${sample}.XXXXXX.tmp")"

    export BLASTDB="$DB_DIR/"
    blastn -db "$BLAST_DB" -taxidlist "$taxid_list" -query "$query" \
        -max_target_seqs 100 -word_size 28 -task megablast \
        -perc_identity 60 -evalue 0.05 -dust yes -qcov_hsp_perc 50 \
        -outfmt '10 qseqid sseqid evalue bitscore pident qlen qstart qend sstart send staxid scomname length' \
        -num_threads "$THREADS" \
        -out "$tmp"

    cat "$header" "$tmp" > "$final"
    rm -f "$tmp"

    python "$BIN_DIR/blast_summary.py" "$final"
}

extract_ev_contigs() {
    log "Extracting EV contigs for sample=$sample (ticket=$ticket)"
    local blast_csv="$RESULTS_DIR/$ticket/blast/${sample}.blast"
    [[ -s "$blast_csv" ]] || { log "No BLAST CSV for $sample; skipping contig extraction"; return 0; }

    local outdir="$RESULTS_DIR/$ticket/ev_contig"
    mkdir -p "$outdir"

    local contigs_scaf="$RESULTS_DIR/$ticket/spades/$sample/scaffolds.fasta"
    local contigs_std="$RESULTS_DIR/$ticket/spades/$sample/contigs.fasta"
    local source_fa
    if [[ -s "$contigs_scaf" ]]; then source_fa="$contigs_scaf"; else source_fa="$contigs_std"; fi
    [[ -s "$source_fa" ]] || { log "No contig fasta for $sample; skipping"; return 0; }

    local fasta_out="$outdir/${sample}.fasta"
    local fasta_200="$outdir/${sample}_200bp.fasta"
    local fasta_cov50="$outdir/${sample}_200bp_minCov50.fasta"

    # Unique contig IDs from BLAST (skip header)
    awk -F, 'NR>1{print $1}' "$blast_csv" | sort -u > "$outdir/${sample}.ids.txt"

    # Extract matching contigs
    if [[ -s "$outdir/${sample}.ids.txt" ]]; then
        seqkit grep -f "$outdir/${sample}.ids.txt" "$source_fa" > "$fasta_out"
    else
        log "No matching contig IDs in BLAST for $sample"
        return 0
    fi

    # Length >= 200, keep gaps (N) with -g
    seqkit seq -m 200 -g "$fasta_out" > "$fasta_200"

    # Coverage > 50 without bc (parse cov_XX.X from header)
    awk '
        /^>/{
            cov=0
            if (match($0,/cov_([0-9.]+)/,m)) cov=m[1]
            if (cov>50) {print substr($0,2)}
        }' "$fasta_200" \
        | seqkit grep -f - "$fasta_200" > "$fasta_cov50" || true

    rm -f "$outdir/${sample}.ids.txt"
}

generate_report() {
    log "Generating genotyping report for sample=$sample (ticket=$ticket)"
    local blast_csv="$RESULTS_DIR/$ticket/blast/${sample}.blast"
    if [[ -s "$blast_csv" ]]; then
        python "$BIN_DIR/html_report.py" --ticket-nr "$ticket" --blast-file "$blast_csv"
    else
        log "No BLAST hits for $sample; report skipped"
    fi
}

# =====================
# Main
# =====================
main() {
    log "START ticket=$ticket sample=$sample threads=$THREADS mem=${MEM_MB}MB"
    prepare_input
    run_taxprofiler
    run_spades
    run_blast
    extract_ev_contigs
    generate_report
    log "DONE ticket=$ticket sample=$sample"
}

main
