#!/usr/bin/env bash
set -euo pipefail

# build_context_genome_db_v4.sh
#
# Build one local nucleotide BLAST database from a small RefSeq genome panel.
#
# Difference from v3:
#   - First queries all taxa and creates one deduplicated accession list.
#   - Then downloads all accessions in a single datasets download call.
#   - Then unpacks all genomes once and builds one BLAST DB.
#
# Requirements:
#   conda install -c bioconda ncbi-datasets-cli blast
#
# Usage:
#   bash build_context_genome_db_v4.sh context_genomes 3
#
# Outputs:
#   context_genomes/context_genomes.fna
#   context_genomes/context_genomes.fna.* BLAST db files
#   context_genomes/all_accessions.txt
#   context_genomes/taxon_accessions.tsv
#   context_genomes/missing_taxa.tsv
#   context_genomes/query_logs/*.raw

#OUTDIR="${1:-context_genomes}"
OUTDIR="${1:-benchmark/context_genomes}"
N_PER_TAXON="${2:-3}"

mkdir -p "$OUTDIR" "$OUTDIR/query_logs" "$OUTDIR/unpacked"

TAXA_FILE="$OUTDIR/context_taxa.tsv"
ALL_ACC="$OUTDIR/all_accessions.txt"
TAXON_ACC="$OUTDIR/taxon_accessions.tsv"
MISSING_FILE="$OUTDIR/missing_taxa.tsv"
COMBINED="$OUTDIR/context_genomes.fna"
ZIP="$OUTDIR/context_genomes.zip"

cat > "$TAXA_FILE" <<'TAXA'
Acinetobacter_baumannii	Acinetobacter baumannii
Burkholderia_cepacia	Burkholderia cepacia
Burkholderia_pseudomallei	Burkholderia pseudomallei
Campylobacter	Campylobacter jejuni
Citrobacter_freundii	Citrobacter freundii
Clostridioides_difficile	Clostridioides difficile
Corynebacterium_diphtheriae	Corynebacterium diphtheriae
Enterobacter_asburiae	Enterobacter asburiae
Enterobacter_cloacae	Enterobacter cloacae
Enterococcus_faecalis	Enterococcus faecalis
Enterococcus_faecium	Enterococcus faecium
Escherichia	Escherichia coli
Haemophilus_influenzae	Haemophilus influenzae
Klebsiella_oxytoca	Klebsiella oxytoca
Klebsiella_pneumoniae	Klebsiella pneumoniae
Neisseria_gonorrhoeae	Neisseria gonorrhoeae
Neisseria_meningitidis	Neisseria meningitidis
Pseudomonas_aeruginosa	Pseudomonas aeruginosa
Salmonella	Salmonella enterica
Serratia_marcescens	Serratia marcescens
Staphylococcus_aureus	Staphylococcus aureus
Staphylococcus_pseudintermedius	Staphylococcus pseudintermedius
Streptococcus_agalactiae	Streptococcus agalactiae
Streptococcus_pneumoniae	Streptococcus pneumoniae
Streptococcus_pyogenes	Streptococcus pyogenes
Vibrio_cholerae	Vibrio cholerae
Vibrio_parahaemolyticus	Vibrio parahaemolyticus
Vibrio_vulnificus	Vibrio vulnificus
TAXA

: > "$ALL_ACC"
: > "$TAXON_ACC"
: > "$MISSING_FILE"
rm -f "$COMBINED" "$ZIP"
rm -rf "$OUTDIR/unpacked"
mkdir -p "$OUTDIR/unpacked"

check_tools() {
    command -v datasets >/dev/null 2>&1 || { echo "ERROR: datasets CLI not found" >&2; exit 1; }
    command -v unzip >/dev/null 2>&1 || { echo "ERROR: unzip not found" >&2; exit 1; }
    command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found" >&2; exit 1; }
}

extract_accessions() {
    local raw="$1"
    grep -oE 'GC[AF]_[0-9]+\.[0-9]+' "$raw" | sort -u | head -n "$N_PER_TAXON" || true
}

summary_refseq() {
    local taxon="$1"
    local assembly_level="$2"
    local raw_out="$3"
    datasets summary genome taxon "$taxon" \
        --assembly-source RefSeq \
        --assembly-level "$assembly_level" \
        --assembly-version latest \
        --mag exclude \
        --limit "$N_PER_TAXON" \
        --report ids_only \
        > "$raw_out" 2> "${raw_out}.stderr" || true
}

check_tools

echo "Querying taxa and building one accession list..." >&2
while IFS=$'\t' read -r label taxon; do
    [[ -z "${label}" ]] && continue
    [[ "${label}" == \#* ]] && continue

    echo "### ${label}: ${taxon}" >&2
    TAX_ACC="$OUTDIR/${label}.accessions.txt"
    RAW_COMPLETE="$OUTDIR/query_logs/${label}.complete_refseq.raw"
    RAW_CHROMOSOME="$OUTDIR/query_logs/${label}.chromosome_refseq.raw"
    : > "$TAX_ACC"

    summary_refseq "$taxon" complete "$RAW_COMPLETE"
    extract_accessions "$RAW_COMPLETE" > "$TAX_ACC"

    if [[ ! -s "$TAX_ACC" ]]; then
        echo "No complete RefSeq genome parsed for ${taxon}; trying chromosome-level RefSeq assemblies" >&2
        summary_refseq "$taxon" chromosome "$RAW_CHROMOSOME"
        extract_accessions "$RAW_CHROMOSOME" > "$TAX_ACC"
    fi

    if [[ ! -s "$TAX_ACC" ]]; then
        echo -e "${label}\t${taxon}\tNO_COMPLETE_OR_CHROMOSOME_REFSEQ_GENOME_FOUND" | tee -a "$MISSING_FILE" >&2
        continue
    fi

    while read -r acc; do
        [[ -z "$acc" ]] && continue
        echo "$acc" >> "$ALL_ACC"
        echo -e "${label}\t${taxon}\t${acc}" >> "$TAXON_ACC"
    done < "$TAX_ACC"

    echo "Selected accessions for ${taxon}:" >&2
    sed 's/^/  /' "$TAX_ACC" >&2

done < "$TAXA_FILE"

sort -u "$ALL_ACC" -o "$ALL_ACC"

if [[ ! -s "$ALL_ACC" ]]; then
    echo "ERROR: accession list is empty; cannot download genomes" >&2
    echo "Check raw query logs in: $OUTDIR/query_logs" >&2
    exit 1
fi

N_ACC=$(wc -l < "$ALL_ACC" | tr -d ' ')
echo "Downloading ${N_ACC} unique genome accessions in one datasets call..." >&2

datasets download genome accession \
    --inputfile "$ALL_ACC" \
    --include genome \
    --filename "$ZIP"

echo "Unpacking genome download..." >&2
unzip -q "$ZIP" -d "$OUTDIR/unpacked"

: > "$COMBINED"
found_fna=0
while IFS= read -r -d '' fna; do
    found_fna=1
    # Append all genomic FASTA sequences. Preserve original headers, but append a tag
    # so hits remain identifiable as coming from this context database.
    awk '
        /^>/ {print $0 " context_db=metapointfinder_context_genomes"; next}
        {print}
    ' "$fna" >> "$COMBINED"
done < <(find "$OUTDIR/unpacked" -type f \( -name '*_genomic.fna' -o -name '*.fna' \) -print0)

if [[ "$found_fna" -eq 0 || ! -s "$COMBINED" ]]; then
    echo "ERROR: downloaded archive did not contain genomic FASTA files" >&2
    echo "Archive: $ZIP" >&2
    echo "Unpacked directory: $OUTDIR/unpacked" >&2
    exit 1
fi

echo "Building local BLAST database..." >&2
makeblastdb \
    -in "$COMBINED" \
    -dbtype nucl \
    -parse_seqids \
    -title metapointfinder_context_genomes \
    -out "$COMBINED"

echo "Built local BLAST database:" >&2
echo "  FASTA:             $COMBINED" >&2
echo "  DB prefix:         $COMBINED" >&2
echo "  Unique accessions: $ALL_ACC" >&2
echo "  Taxon mapping:     $TAXON_ACC" >&2
if [[ -s "$MISSING_FILE" ]]; then
    echo "Some taxa were missing; see: $MISSING_FILE" >&2
fi
