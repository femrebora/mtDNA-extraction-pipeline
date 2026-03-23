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
# Dependencies:
#   samtools, bwa, bcftools
# =============================================================================

# Exit immediately if any command fails, treat unset variables as errors,
# and propagate errors through pipes (prevents silent failures).
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Sample identifier — used as a prefix for every output file.
PREFIX="ECE25-427-OHD"

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

# =============================================================================
# STEP 1 — INDEX THE INPUT BAM
# =============================================================================
# A BAI index is required for random-access retrieval of reads from a specific
# genomic region (here, the mitochondrial contig).

samtools index "$BAM"

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

samtools view "$BAM" "$MTCONTIG" \
  | awk '{print $1}' | sort -u > "${PREFIX}.mt.readnames.txt"

# =============================================================================
# STEP 3 — EXTRACT BOTH MATES OF EACH MITOCHONDRIAL READ PAIR
# =============================================================================
# Using the collected read names, pull every read that belongs to a chrM-linked
# pair from the full BAM — including mates that mapped elsewhere.
#
#   -b  : write output in BAM format
#   -N  : filter reads whose name is listed in the supplied file

samtools view -b -N "${PREFIX}.mt.readnames.txt" "$BAM" > "${PREFIX}.mt.anymate.bam"

# =============================================================================
# STEP 4 — NAME-SORT FOR FASTQ CONVERSION
# =============================================================================
# samtools fastq requires read pairs to be adjacent in the stream. Name-sorting
# guarantees that R1 and R2 of each pair appear consecutively regardless of
# their original genomic coordinates.
#
#   -n    : sort by read name
#   -@ 8  : use 8 threads

samtools sort -n -@ 8 -o "${PREFIX}.mt.anymate.namesort.bam" "${PREFIX}.mt.anymate.bam"

# =============================================================================
# STEP 5 — CONVERT TO PAIRED FASTQ
# =============================================================================
# Split the name-sorted BAM into separate gzipped FASTQ files for each read
# orientation and handle singletons separately.
#
#   -1 / -2  : paired forward and reverse reads
#   -s        : singleton reads (mate is unmapped or absent)
#   -0        : reads that don't match any of the above → discarded
#   -n        : preserve original read names without appending /1 or /2
#   -@ 8      : use 8 threads

samtools fastq -@ 8 \
  -1 "${PREFIX}.mt_any_R1.fastq.gz" \
  -2 "${PREFIX}.mt_any_R2.fastq.gz" \
  -s "${PREFIX}.mt_any_singletons.fastq.gz" \
  -0 /dev/null -n \
  "${PREFIX}.mt.anymate.namesort.bam"

# =============================================================================
# STEP 6 — QC: READ COUNT SUMMARY
# =============================================================================
# Each FASTQ record spans exactly 4 lines (header, sequence, '+', quality),
# so dividing the total line count by 4 gives the number of reads.

echo "R1 reads:"        $(( $(gzip -cd "${PREFIX}.mt_any_R1.fastq.gz"        | wc -l) / 4 ))
echo "R2 reads:"        $(( $(gzip -cd "${PREFIX}.mt_any_R2.fastq.gz"        | wc -l) / 4 ))
echo "Singleton reads:" $(( $(gzip -cd "${PREFIX}.mt_any_singletons.fastq.gz" | wc -l) / 4 ))
# List files with human-readable sizes (-h) to verify they were written correctly.
# Note: singletons are retained here for record-keeping but are NOT passed to
# BWA-MEM in Step 8. Only the properly paired R1/R2 files are re-aligned,
# because BWA-MEM paired-end mode expects matched read pairs.
ls -lh "${PREFIX}.mt_any_R"*.fastq.gz "${PREFIX}.mt_any_singletons.fastq.gz"

# =============================================================================
# STEP 7 — INDEX THE rCRS REFERENCE
# =============================================================================
# Two separate indexes are required:
#   samtools faidx  → creates a .fai index for coordinate-based FASTA access
#   bwa index       → builds the BWT/suffix-array structures for alignment

samtools faidx "$REF"
bwa index "$REF"

# =============================================================================
# STEP 8 — RE-ALIGN MITOCHONDRIAL READS TO rCRS
# =============================================================================
# Align the paired-end reads specifically to the rCRS reference using
# BWA-MEM with 16 threads. Piping directly into samtools sort avoids writing
# an unsorted intermediate SAM file to disk.
#
# Re-alignment to a mitochondria-only reference improves mapping quality
# compared to the original whole-genome alignment, as reads no longer compete
# with NUMTs (nuclear copies of mitochondrial sequences).

bwa mem -t 16 "$REF" \
  "${PREFIX}.mt_any_R1.fastq.gz" \
  "${PREFIX}.mt_any_R2.fastq.gz" \
| samtools sort -@ 8 -o "${PREFIX}.mt.bam"

# Index the coordinate-sorted mitochondrial BAM.
samtools index "${PREFIX}.mt.bam"

# =============================================================================
# STEP 9 — VARIANT CALLING
# =============================================================================
# Three-step bcftools pipeline:
#
#   mpileup  : compute per-base pileup statistics from the re-aligned BAM,
#              using the rCRS as reference. Outputs uncompressed BCF (-Ou)
#              to avoid redundant compression between pipe stages.
#
#   call     : call SNPs and small indels using the multiallelic model (-m),
#              reporting only variant sites (-v).
#
#   norm     : left-align and normalise indel representations relative to
#              the reference (-f), then write a bgzipped VCF (-Oz) that can
#              be indexed and used by downstream tools.

bcftools mpileup -f "$REF" -Ou "${PREFIX}.mt.bam" \
| bcftools call -mv -Ou \
| bcftools norm -f "$REF" -Oz -o "${PREFIX}.mt.vcf.gz"

# Index the compressed VCF for fast random access by bcftools consensus.
bcftools index "${PREFIX}.mt.vcf.gz"

# =============================================================================
# STEP 10 — BUILD CONSENSUS FASTA FOR MITOMASTER
# =============================================================================
# Apply the called variants to the rCRS sequence to produce a sample-specific
# mitochondrial consensus. The output FASTA header is replaced with the sample
# identifier so that MITOMASTER correctly labels the uploaded sequence.

bcftools consensus -f "$REF" "${PREFIX}.mt.vcf.gz" > "${PREFIX}.mt.consensus.fasta"

# Replace the default header line (inherited from the reference) with the
# sample ID. The resulting FASTA is ready for direct upload to MITOMASTER
# for haplogroup assignment and variant interpretation.
sed -i "1s/^>.*/>${PREFIX}/" "${PREFIX}.mt.consensus.fasta"
