# BAM2MITO

A Bash pipeline for extracting patient mitochondrial DNA from a whole-genome or whole-exome BAM file, re-aligning it to the revised Cambridge Reference Sequence (rCRS), calling variants, and producing a consensus FASTA ready for upload to [MITOMASTER](https://www.mitomap.org/MITOMAP/MitoMaster).

---

## Overview

Standard sequencing pipelines align reads to the full human reference genome (GRCh37/GRCh38), which includes the mitochondrial contig (`chrM`). However, dedicated analysis of mitochondrial variants requires:

1. Isolating reads that originated from the mitochondrial genome.
2. Re-aligning them specifically to the rCRS reference to maximize mapping accuracy.
3. Calling variants and building a sample consensus sequence.
4. Uploading the consensus to MITOMASTER for haplogroup assignment and pathogenicity assessment.

This pipeline automates all four steps.

---

## Reference Genome

The file `sequence.fasta` contains the **revised Cambridge Reference Sequence (rCRS)**, accession **NC_012920.1**, downloaded directly from [MITOMAP](https://www.mitomap.org). It is the standard human mitochondrial reference used by MITOMASTER for variant interpretation.

---

## Dependencies

| Tool | Purpose |
|------|---------|
| `samtools` | BAM indexing, read extraction, name-sorting, FASTQ conversion |
| `bwa` | Reference indexing and short-read alignment (BWA-MEM) |
| `bcftools` | Variant calling, VCF normalisation, consensus sequence generation |

---

## Input Requirements

The input BAM (`{PREFIX}.recal.bam`) must be a **base-quality score recalibrated (BQSR)** alignment file. BQSR is a GATK pre-processing step that corrects systematic errors in base quality scores reported by the sequencer, improving the accuracy of downstream variant calling. It is typically produced as the final step of the GATK Best Practices pre-processing workflow (MarkDuplicates → BQSR).

---

## Usage

1. Place your recalibrated BAM file in the same directory as the script.
2. Edit the `PREFIX` variable at the top of `BAM2MITO.sh` to match your sample ID.
3. Run the script:

```bash
bash BAM2MITO.sh
```

---

## Pipeline Steps

```
Input BAM (whole-genome/exome, BQSR recalibrated)
        │
        ▼
[Step 1]  Index BAM → samtools index
        │
        ▼
[Step 2]  Extract chrM read names → samtools view | awk | sort -u
        │
        ▼
[Step 3]  Pull both mates of each pair → samtools view -N
        │
        ▼
[Step 4]  Name-sort for FASTQ conversion → samtools sort -n
        │
        ▼
[Step 5]  Convert to paired FASTQ → samtools fastq
          (R1, R2, and singletons produced)
        │
        ▼
[Step 6]  QC: print read counts + file sizes
          (singletons kept for records but NOT re-aligned)
        │
        ▼
[Step 7]  Index rCRS reference → samtools faidx + bwa index
        │
        ▼
[Step 8]  Re-align paired reads (R1 + R2 only) to rCRS → bwa mem | samtools sort
        │
        ▼
[Step 9]  Call variants → bcftools mpileup | call | norm → compressed VCF
        │
        ▼
[Step 10] Build consensus FASTA → bcftools consensus + rename header
        │
        ▼
Output: {PREFIX}.mt.consensus.fasta  ← upload to MITOMASTER
```

---

## Output Files

| File | Description |
|------|-------------|
| `{PREFIX}.mt.readnames.txt` | Read names of all chrM-mapped reads |
| `{PREFIX}.mt.anymate.bam` | BAM containing both mates of chrM read pairs |
| `{PREFIX}.mt.anymate.namesort.bam` | Name-sorted version of the above |
| `{PREFIX}.mt_any_R1.fastq.gz` | Forward reads (gzipped FASTQ) |
| `{PREFIX}.mt_any_R2.fastq.gz` | Reverse reads (gzipped FASTQ) |
| `{PREFIX}.mt_any_singletons.fastq.gz` | Unpaired/singleton reads |
| `{PREFIX}.mt.bam` | Final re-aligned mitochondrial BAM |
| `{PREFIX}.mt.vcf.gz` | Called and normalised variants (bgzipped VCF) |
| `{PREFIX}.mt.consensus.fasta` | Sample consensus mtDNA sequence for MITOMASTER |

---

## MITOMASTER Upload

After the pipeline completes, upload `{PREFIX}.mt.consensus.fasta` to [MITOMASTER](https://www.mitomap.org/MITOMAP/MitoMaster) for:

- Haplogroup assignment
- Variant annotation against the MITOMAP database
- Pathogenicity and disease association reporting
