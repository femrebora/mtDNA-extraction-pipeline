#!/usr/bin/env bash
# Exits immediately on error, treats unset variables as errors, and propagates pipe failures.
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Sample identifier used as a prefix for all output files.
PREFIX="ECE25-427-OHD"

# Input BAM file: base-quality score recalibrated whole-genome or whole-exome BAM.
BAM="${PREFIX}.recal.bam"

# Name of the mitochondrial contig as it appears in the BAM header (e.g., chrM or MT).
MTCONTIG="chrM"

# Reference FASTA: the revised Cambridge Reference Sequence (rCRS, NC_012920.1)
# downloaded from MITOMAP. Used for re-alignment and variant calling.
REF="sequence.fasta"

# ---------------------------------------------------------------------------
# Step 1 – Extract mitochondrial reads from the original BAM
# ---------------------------------------------------------------------------
# Index the input BAM so that samtools can perform random-access retrieval.
samtools index "$BAM"

# Collect the read names of every read that mapped to the mitochondrial contig.
# Both mates of a pair may not map to chrM, so we first gather names, then
# pull all reads (including the nuclear-mapping mate) in the next command.
samtools view "$BAM" "$MTCONTIG" \
  | awk '{print $1}' | sort -u > "${PREFIX}.mt.readnames.txt"

# Extract all reads whose name appears in the list (i.e., either mate maps to chrM).
# -b  : output in BAM format
# -N  : filter by read name file
samtools view -b -N "${PREFIX}.mt.readnames.txt" "$BAM" > "${PREFIX}.mt.anymate.bam"

# ---------------------------------------------------------------------------
# Step 2 – Sort by name and convert to FASTQ
# ---------------------------------------------------------------------------
# Name-sort is required so that paired reads are adjacent in the stream,
# which is a prerequisite for correct paired FASTQ interleaving by samtools fastq.
samtools sort -n -@ 8 -o "${PREFIX}.mt.anymate.namesort.bam" "${PREFIX}.mt.anymate.bam"

# Convert the name-sorted BAM to gzipped FASTQ files.
# -1 / -2 : forward and reverse reads of proper pairs
# -s       : singleton reads (mate unmapped or otherwise unpaired)
# -0       : reads that do not fit the above categories (discarded to /dev/null)
# -n       : use the original read name without appending /1 or /2
samtools fastq -@ 8 \
  -1 "${PREFIX}.mt_any_R1.fastq.gz" \
  -2 "${PREFIX}.mt_any_R2.fastq.gz" \
  -s "${PREFIX}.mt_any_singletons.fastq.gz" \
  -0 /dev/null -n \
  "${PREFIX}.mt.anymate.namesort.bam"

# ---------------------------------------------------------------------------
# Step 3 – Quick QC: print read counts for each output file
# ---------------------------------------------------------------------------
# FASTQ records are 4 lines each, so dividing the total line count by 4 gives
# the number of reads.
echo "R1 reads:" $(( $(gzip -cd "${PREFIX}.mt_any_R1.fastq.gz" | wc -l) / 4 ))
echo "R2 reads:" $(( $(gzip -cd "${PREFIX}.mt_any_R2.fastq.gz" | wc -l) / 4 ))
echo "Singleton reads:" $(( $(gzip -cd "${PREFIX}.mt_any_singletons.fastq.gz" | wc -l) / 4 ))
ls -lh "${PREFIX}.mt_any_R"*.fastq.gz "${PREFIX}.mt_any_singletons.fastq.gz"

# ---------------------------------------------------------------------------
# Step 4 – Re-align mitochondrial reads to the rCRS reference
# ---------------------------------------------------------------------------
# Index the reference FASTA for both samtools (FASTA index) and BWA (BWT index).
samtools faidx "$REF"
bwa index "$REF"

# Align paired-end reads to the rCRS reference using BWA-MEM with 16 threads,
# then coordinate-sort and write the output directly to a BAM file.
bwa mem -t 16 "$REF" \
  "${PREFIX}.mt_any_R1.fastq.gz" \
  "${PREFIX}.mt_any_R2.fastq.gz" \
| samtools sort -@ 8 -o "${PREFIX}.mt.bam"

# Index the final mitochondrial BAM for downstream tools.
samtools index "${PREFIX}.mt.bam"

# ---------------------------------------------------------------------------
# Step 5 – Variant calling and consensus FASTA generation
# ---------------------------------------------------------------------------
# mpileup : compute per-base pileup statistics from the BAM
# call    : call SNPs and indels (-m multiallelic caller, -v variants only)
# norm    : left-align and normalise indels relative to the reference,
#           then write a bgzipped VCF
bcftools mpileup -f "$REF" -Ou "${PREFIX}.mt.bam" \
| bcftools call -mv -Ou \
| bcftools norm -f "$REF" -Oz -o "${PREFIX}.mt.vcf.gz"

# Index the compressed VCF so that bcftools consensus can retrieve variants quickly.
bcftools index "${PREFIX}.mt.vcf.gz"

# Apply the called variants to the rCRS reference to produce a sample-specific
# mitochondrial consensus sequence in FASTA format.
bcftools consensus -f "$REF" "${PREFIX}.mt.vcf.gz" > "${PREFIX}.mt.consensus.fasta"

# Replace the default FASTA header with the sample identifier so that the
# sequence is correctly labelled when uploaded to MITOMASTER for haplogroup
# assignment and variant interpretation.
sed -i "1s/^>.*/>${PREFIX}/" "${PREFIX}.mt.consensus.fasta"
