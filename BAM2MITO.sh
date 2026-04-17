#!/usr/bin/env bash
# =============================================================================
# BAM2MITO.sh
# =============================================================================
# Extracts mitochondrial DNA reads from a recalibrated whole-genome or
# whole-exome BAM file, re-aligns them to the revised Cambridge Reference
# Sequence (rCRS, NC_012920.1), calls variants, and generates a consensus
# FASTA ready for upload to MITOMASTER (https://www.mitomap.org).
#
# Usage:
#   bash BAM2MITO.sh
#
# Dependencies (minimum versions):
#   samtools >= 1.12, bwa, bcftools
#
# Known limitations:
#   1. bcftools --ploidy 1 detects homoplasmic variants only. For low-level
#      heteroplasmy, use: gatk Mutect2 --mitochondria-mode
#   2. NuMT (nuclear mitochondrial segment) reads are not filtered. For
#      stricter analyses, re-align to a combined nuclear+rCRS reference and
#      retain only reads whose best alignment maps to the mt contig.
# =============================================================================

# Exit immediately if any command fails, treat unset variables as errors,
# and propagate errors through pipes (prevents silent failures).
set -euo pipefail

# Collect intermediate files; remove them on any exit (success or failure).
TMPFILES=()
trap 'rm -f "${TMPFILES[@]}" 2>/dev/null || true' EXIT

# =============================================================================
# CONFIGURATION
# =============================================================================

# Sample identifier — used as a prefix for every output file.
PREFIX="your_bam_file_prefix"

# Input BAM file: base-quality score recalibrated alignment (from GATK BQSR
# or equivalent) covering the full genome or exome including chrM.
BAM="${PREFIX}.recal.bam"

# Name of the mitochondrial contig as it appears in the BAM header.
# Typical values: "chrM" (GRCh38/hg38) or "MT" (GRCh37/hg19, Ensembl).
MTCONTIG="chrM"

# Reference FASTA: the revised Cambridge Reference Sequence (rCRS),
# accession NC_012920.1, downloaded from MITOMAP. Used for re-alignment
# and variant calling.
REF="sequence.fasta"

# Number of threads for parallelisable steps.
THREADS=16

# Maximum pileup depth for bcftools mpileup. The default (250×) is far too
# low for mitochondrial DNA, which routinely reaches thousands of fold
# coverage. Without raising this cap, high-depth positions are silently
# truncated and variants can be missed.
MAX_DEPTH=500000

# =============================================================================
# STEP 0 — DEPENDENCY AND INPUT CHECKS
# =============================================================================
# Fail early with a clear message rather than obscure mid-run errors.

for tool in samtools bwa bcftools; do
    command -v "$tool" >/dev/null 2>&1 \
        || { echo "ERROR: $tool not found in PATH" >&2; exit 1; }
done

# Confirm samtools is at least 1.12; the -N flag (filter by read-name file)
# used in Step 3 was introduced in that release.
SAMTOOLS_VER=$(samtools --version | awk 'NR==1 {print $2}')
SAMTOOLS_MAJOR=$(echo "$SAMTOOLS_VER" | cut -d. -f1)
SAMTOOLS_MINOR=$(echo "$SAMTOOLS_VER" | cut -d. -f2)
if [[ "$SAMTOOLS_MAJOR" -lt 1 ]] || \
   { [[ "$SAMTOOLS_MAJOR" -eq 1 ]] && [[ "$SAMTOOLS_MINOR" -lt 12 ]]; }; then
    echo "ERROR: samtools >= 1.12 required (found ${SAMTOOLS_VER})" >&2
    exit 1
fi

[[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[[ -f "$REF" ]] || { echo "ERROR: Reference not found: $REF" >&2; exit 1; }

# =============================================================================
# STEP 1 — INDEX THE INPUT BAM
# =============================================================================
# A BAI index is required for random-access retrieval of reads from a specific
# genomic region (here, the mitochondrial contig) and for idxstats in the
# contig-name check immediately below.

echo "[1/8] Indexing BAM..."
samtools index -@ "$THREADS" "$BAM"

# Verify the configured contig name matches what is actually in this BAM.
# Catches the common chrM vs MT vs NC_012920.1 naming mismatch, which would
# cause Step 2 to silently return zero reads.
CONTIG_READS=$(samtools idxstats "$BAM" \
    | awk -v c="$MTCONTIG" '$1==c {print $3+$4}') || true
if [[ -z "$CONTIG_READS" || "$CONTIG_READS" -eq 0 ]]; then
    echo "ERROR: No reads found for contig '${MTCONTIG}' in ${BAM}." >&2
    echo "       Contigs present in this BAM:" >&2
    samtools idxstats "$BAM" | awk '$3>0 {print "  "$1}' >&2
    echo "       Update MTCONTIG to the correct name and re-run." >&2
    exit 1
fi

# =============================================================================
# STEP 2 — COLLECT MITOCHONDRIAL READ NAMES
# =============================================================================
# Fetch all reads that aligned to the mitochondrial contig and extract their
# query names (column 1 of SAM). Sorting and deduplication (-u) ensure each
# read name appears only once even if it was reported multiple times.
#
# We collect names rather than reads at this stage because paired-end reads
# may have one mate on chrM and the other elsewhere in the genome. We need
# both mates for accurate re-alignment.

echo "[2/8] Collecting ${MTCONTIG} read names..."
# -F 2308 excludes UNMAP(4)+SECONDARY(256)+SUPPLEMENTARY(2048) so that
# only primary alignments seed the read-name list.
samtools view -F 2308 "$BAM" "$MTCONTIG" \
    | awk '{print $1}' | sort -u > "${PREFIX}.mt.readnames.txt"
TMPFILES+=("${PREFIX}.mt.readnames.txt")
echo "      Read names found: $(wc -l < "${PREFIX}.mt.readnames.txt")"

# =============================================================================
# STEP 3 — EXTRACT BOTH MATES OF EACH MITOCHONDRIAL READ PAIR
# =============================================================================
# Using the collected read names, pull every read that belongs to a chrM-linked
# pair from the full BAM — including mates that mapped elsewhere.
#
#   -b  : write output in BAM format
#   -N  : filter reads whose name is listed in the supplied file (requires
#         samtools >= 1.12; checked in Step 0)

echo "[3/8] Extracting both mates..."
samtools view -b -N "${PREFIX}.mt.readnames.txt" "$BAM" > "${PREFIX}.mt.anymate.bam"
TMPFILES+=("${PREFIX}.mt.anymate.bam")

# =============================================================================
# STEP 4 — NAME-SORT AND CONVERT TO FASTQ
# =============================================================================
# samtools fastq requires read pairs to be adjacent in the stream. Name-sorting
# guarantees that R1 and R2 of each pair appear consecutively regardless of
# their original genomic coordinates.
#
# FASTQ flags:
#   -1 / -2  : paired forward and reverse reads
#   -s        : singleton reads (mate is unmapped or absent)
#   -0        : reads that match none of the above → discarded via /dev/null
#   -n        : preserve original read names without appending /1 or /2,
#               required for correct downstream handling by BWA-MEM

echo "[4/8] Name-sorting and converting to FASTQ..."
samtools sort -n -@ "$THREADS" \
    -o "${PREFIX}.mt.anymate.namesort.bam" "${PREFIX}.mt.anymate.bam"
TMPFILES+=("${PREFIX}.mt.anymate.namesort.bam")

samtools fastq -@ "$THREADS" \
    -1 "${PREFIX}.mt_any_R1.fastq.gz" \
    -2 "${PREFIX}.mt_any_R2.fastq.gz" \
    -s "${PREFIX}.mt_any_singletons.fastq.gz" \
    -0 /dev/null -n \
    "${PREFIX}.mt.anymate.namesort.bam"
TMPFILES+=("${PREFIX}.mt_any_R1.fastq.gz" "${PREFIX}.mt_any_R2.fastq.gz" \
           "${PREFIX}.mt_any_singletons.fastq.gz")

# Each FASTQ record spans exactly 4 lines (header, sequence, '+', quality),
# so dividing the total line count by 4 gives the number of reads.
echo "R1 reads:        $(( $(gzip -cd "${PREFIX}.mt_any_R1.fastq.gz"         | wc -l) / 4 ))"
echo "R2 reads:        $(( $(gzip -cd "${PREFIX}.mt_any_R2.fastq.gz"         | wc -l) / 4 ))"
echo "Singleton reads: $(( $(gzip -cd "${PREFIX}.mt_any_singletons.fastq.gz" | wc -l) / 4 ))"
ls -lh "${PREFIX}.mt_any_R"*.fastq.gz "${PREFIX}.mt_any_singletons.fastq.gz"

# =============================================================================
# STEP 5 — INDEX THE rCRS REFERENCE
# =============================================================================
# Two separate indexes are required:
#   samtools faidx  → creates a .fai index for coordinate-based FASTA access
#   bwa index       → builds the BWT/suffix-array structures for alignment
# Each index is skipped if it already exists to avoid unnecessary re-work.

echo "[5/8] Indexing rCRS reference..."
samtools faidx "$REF"
if [[ ! -f "${REF}.amb" || ! -f "${REF}.ann" || ! -f "${REF}.bwt" || \
      ! -f "${REF}.pac" || ! -f "${REF}.sa" ]]; then
    bwa index "$REF"
fi

# =============================================================================
# STEP 6 — RE-ALIGN MITOCHONDRIAL READS TO rCRS
# =============================================================================
# Align paired reads and singleton reads separately to the rCRS reference,
# then merge the results before marking duplicates.
#
# Aligning paired and singleton reads in separate BWA-MEM calls is necessary
# because BWA-MEM paired-end mode expects matched R1/R2 pairs; mixing
# unpaired reads into a paired-end run would corrupt mate information.
#
# A read-group header (-R) is mandatory for downstream tools (GATK, bcftools)
# to correctly identify samples and libraries.
#
# Re-alignment to a mitochondria-only reference improves mapping quality
# compared to the original whole-genome alignment because reads no longer
# compete with NuMTs (nuclear copies of mitochondrial sequences) present
# in the whole-genome reference.

echo "[6/8] Aligning to rCRS..."
RG="@RG\tID:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}\tPU:${PREFIX}"

# Paired reads: align and name-sort (name-sort is required by fixmate in Step 7)
bwa mem -t "$THREADS" -R "$RG" "$REF" \
    "${PREFIX}.mt_any_R1.fastq.gz" \
    "${PREFIX}.mt_any_R2.fastq.gz" \
    | samtools sort -n -@ "$THREADS" -o "${PREFIX}.mt.paired.namesort.bam"

# Singleton reads: align only if any exist, then merge with the paired BAM.
SINGLETON_COUNT=$(( $(gzip -cd "${PREFIX}.mt_any_singletons.fastq.gz" | wc -l) / 4 ))
if [[ "$SINGLETON_COUNT" -gt 0 ]]; then
    echo "      Aligning ${SINGLETON_COUNT} singleton reads..."
    bwa mem -t "$THREADS" -R "$RG" "$REF" \
        "${PREFIX}.mt_any_singletons.fastq.gz" \
        | samtools sort -n -@ "$THREADS" -o "${PREFIX}.mt.singletons.namesort.bam"

    # -n : inputs are name-sorted; merge preserves that order
    samtools merge -f -n -@ "$THREADS" \
        "${PREFIX}.mt.merged.namesort.bam" \
        "${PREFIX}.mt.paired.namesort.bam" \
        "${PREFIX}.mt.singletons.namesort.bam"
    rm "${PREFIX}.mt.paired.namesort.bam" "${PREFIX}.mt.singletons.namesort.bam"
else
    echo "      No singleton reads to merge."
    mv "${PREFIX}.mt.paired.namesort.bam" "${PREFIX}.mt.merged.namesort.bam"
fi

# =============================================================================
# STEP 7 — FIXMATE, COORDINATE SORT, AND MARK DUPLICATES
# =============================================================================
# Three sub-steps are required for correct duplicate marking:
#
#   fixmate -m  : adds mate-score tags (ms) to each read. These tags record a
#                 quality-based score for the mate, allowing markdup to select
#                 the best-quality read in a duplicate cluster rather than an
#                 arbitrary one. Requires name-sorted input.
#
#   sort        : coordinate-sort the fixmate output, as markdup requires
#                 reads ordered by genomic position.
#
#   markdup     : identifies and flags duplicate read pairs. Duplicate reads
#                 remain in the BAM with the 0x400 (DUP) flag set; they are
#                 handled at the mpileup stage in Step 8.

echo "[7/8] fixmate → coordinate sort → markdup..."
samtools fixmate -m -@ "$THREADS" \
    "${PREFIX}.mt.merged.namesort.bam" - \
    | samtools sort -@ "$THREADS" -o "${PREFIX}.mt.prededup.bam"

rm "${PREFIX}.mt.merged.namesort.bam"

samtools markdup -@ "$THREADS" \
    "${PREFIX}.mt.prededup.bam" \
    "${PREFIX}.mt.bam"

samtools index -@ "$THREADS" "${PREFIX}.mt.bam"
rm "${PREFIX}.mt.prededup.bam"

echo "── Alignment QC ──"
samtools flagstat "${PREFIX}.mt.bam"
samtools coverage  "${PREFIX}.mt.bam"

# =============================================================================
# STEP 8 — VARIANT CALLING, FILTERING, AND CONSENSUS FASTA
# =============================================================================
# Three-step bcftools pipeline:
#
#   mpileup  : compute per-base pileup statistics from the re-aligned BAM,
#              using the rCRS as reference. Key flags:
#
#              -d MAX_DEPTH       : raise the per-position depth cap from the
#                                   default 250× to 500,000×. Without this,
#                                   high-coverage mtDNA positions are silently
#                                   truncated and variants are missed.
#
#              --ns               : Skip reads with any of these flags set.
#                                   DUP is deliberately excluded from this
#                                   list. For high-copy mitochondrial DNA,
#                                   excluding duplicate-flagged reads can
#                                   severely reduce effective depth and distort
#                                   allele frequencies. Unmapped, secondary,
#                                   and QC-failed reads are still excluded.
#                                   (bcftools 1.19+ uses --ns instead of the
#                                   older --excl-flags)
#
#              -a FORMAT/AD,DP    : emit allele depth (AD) and total depth (DP)
#                                   per sample; required for the quality filter.
#
#   call     : call SNPs and indels using the multiallelic model (-m) and
#              report only variant sites (-v). --ploidy 1 applies a haploid
#              model appropriate for mitochondria (homoplasmic variants). For
#              heteroplasmy, use GATK Mutect2 --mitochondria-mode instead.
#
#   norm     : left-align and normalise indel representations relative to
#              the reference (-f), then write a bgzipped VCF (-Oz).

echo "[8/8] Calling variants..."
bcftools mpileup \
    -f "$REF" \
    -d "$MAX_DEPTH" \
    --ns UNMAP,SECONDARY,QCFAIL \
    -a FORMAT/AD,FORMAT/DP \
    -Ou \
    "${PREFIX}.mt.bam" \
| bcftools call \
    --ploidy 1 \
    -mv \
    -Ou \
| bcftools norm \
    -f "$REF" \
    -Oz \
    -o "${PREFIX}.mt.raw.vcf.gz"

bcftools index "${PREFIX}.mt.raw.vcf.gz"

echo "── Raw variant stats ──"
bcftools stats "${PREFIX}.mt.raw.vcf.gz" | grep "^SN"

# Quality filter: mark sites that fail QUAL >= 20 or FORMAT/DP >= 100 as
# LOWQUAL. mtDNA routinely exceeds 1,000× coverage, so DP < 100 signals a
# genuine mapping or assembly gap rather than normal low-coverage variation.
# An additional allele-frequency gate (ALT AD / total AD >= 0.95) enforces
# that only near-homoplasmic variants — not NuMT noise or systematic errors —
# are applied to the consensus FASTA.
# The raw VCF is kept alongside the filtered VCF for traceability.
bcftools filter \
    -s LOWQUAL \
    -i 'QUAL>=20 && FORMAT/DP>=100' \
    "${PREFIX}.mt.raw.vcf.gz" \
| bcftools view \
    -i 'FILTER="PASS" && (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1])) >= 0.95' \
    -Oz -o "${PREFIX}.mt.vcf.gz"

bcftools index "${PREFIX}.mt.vcf.gz"

echo "── Filtered variant stats ──"
bcftools stats "${PREFIX}.mt.vcf.gz" | grep "^SN"

# Apply PASS variants (already the only entries in mt.vcf.gz) to the rCRS
# sequence to produce a sample-specific mitochondrial consensus.
bcftools consensus -f "$REF" "${PREFIX}.mt.vcf.gz" > "${PREFIX}.mt.consensus.fasta"

# Replace the default header line (inherited from the reference) with the
# sample ID. The resulting FASTA is ready for direct upload to MITOMASTER
# for haplogroup assignment and variant interpretation.
sed -i "1s/^>.*/>${PREFIX}/" "${PREFIX}.mt.consensus.fasta"

echo ""
echo "════════════════════════════════════════"
echo " Done. Output files:"
echo "   BAM (markdup):   ${PREFIX}.mt.bam"
echo "   VCF (raw):       ${PREFIX}.mt.raw.vcf.gz"
echo "   VCF (filtered):  ${PREFIX}.mt.vcf.gz"
echo "   Consensus FASTA: ${PREFIX}.mt.consensus.fasta"
echo "════════════════════════════════════════"
