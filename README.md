# build_lof_matrix.R

Collapses Loss-of-Function (LoF) variants per gene from a SnpEff-annotated VCF and outputs a genotype matrix.

```
Rows    = samples
Columns = genes (or transcripts)
Values  = LoF dosage (integer, bounded by ploidy)
```

Works with any ploidy and handles both phased and unphased VCFs.

---

## Requirements

- R (≥ 4.0)
- Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
```

- `parallel` (base R, no install needed)

---

## Input requirements

1. **SnpEff-annotated VCF** — the script reads the `ANN` INFO field. If it's not present the script will stop with an error.

2. **Biallelic sites only** — multiallelic sites should be split before running. The script will warn you if it detects any.

```bash
# Annotate with SnpEff
java -jar snpEff.jar <genome> input.vcf.gz > annotated.vcf.gz

# Split multiallelic sites
bcftools norm -m -any annotated.vcf.gz \
    | bcftools norm -f reference.fa \
    -O z -o input_biallelic.vcf.gz

# Index
bcftools index -t input_biallelic.vcf.gz
```

---

## Usage

```bash
Rscript build_lof_matrix.R \
    --vcf    input_biallelic.vcf.gz \
    --out    LoF_matrix.tsv \
    --impact HIGH \
    --by     gene \
    --threads 4
```

### Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--vcf` | Yes | — | Input VCF/VCF.gz (SnpEff-annotated, biallelic) |
| `--out` | No | `LoF_genotype_matrix.tsv` | Output TSV file |
| `--impact` | No | `HIGH` | SnpEff impact level to include. Use `""` to keep all levels |
| `--by` | No | `gene` | Collapse by `gene` or `transcript` |
| `--threads` | No | `1` | Number of parallel threads. Linux/macOS only — ignored on Windows |

---

## Dosage encoding

### Unphased VCF

For each gene, the dosage is the **maximum ALT allele count** seen at any single LoF variant across that gene. This is a conservative per-site estimate since cis/trans phase cannot be determined without phasing.

| Genotype | Dosage |
|---|---|
| `0/0` | 0 |
| `0/1` | 1 |
| `1/1` | 2 |
| `0/1/1/1` (tetraploid) | 3 |
| `1/1/1/1` (tetraploid) | 4 |

### Phased VCF

For each gene, the dosage is the **number of distinct haplotypes carrying at least one LoF allele** across all LoF variants. This correctly distinguishes cis stacking (both LoF alleles on the same haplotype) from trans compound heterozygotes (LoF alleles on different haplotypes).

| Example | Dosage | Interpretation |
|---|---|---|
| `0\|1` only | 1 | one haplotype affected |
| `0\|1` + `0\|1` (cis) | 1 | both variants on same haplotype |
| `0\|1` + `1\|0` (trans) | 2 | compound heterozygote |
| `1\|1` | 2 | both haplotypes affected |

If any variant in a gene is unphased for a given sample, that sample falls back to the unphased (max dosage) calculation.

---

## Output

A tab-separated matrix with samples as rows and genes as columns. The first column is `Accession` (sample name).

```
Accession    GENE1    GENE2    GENE3    ...
SAMPLE_001   0        1        2
SAMPLE_002   1        0        0
...
```

A summary is printed to the console after the run:

```
Phasing path summary (gene × sample cells):
  Phased   : 8420
  Unphased : 312
  Missing  : 51

Dosage distribution:
   0      1      2
56893452 12750827 5628295
```

---

## Notes on phased vs unphased VCFs

If your VCF was phased with a tool like BEAGLE, be aware that BEAGLE re-estimates genotype probabilities during phasing. For low-coverage sequencing data in particular, low-confidence homozygous ALT calls (`1/1`) can be revised to `0|0` based on population haplotype context. This is expected behaviour, not a bug — the `1/1` call was uncertain to begin with. You can check the `GP` (genotype posterior) field in the original unphased VCF to assess confidence.

---

## Example SLURM script

```bash
#!/bin/bash
#SBATCH -J LoF_matrix
#SBATCH -t 24:00:00
#SBATCH --mem=128g
#SBATCH --cpus-per-task=8

module load R

Rscript build_lof_matrix.R \
    --vcf     input_biallelic.vcf.gz \
    --out     LoF_matrix.tsv \
    --threads 8
```

---

## Citation

If you use this script, please cite the tools it depends on:

- **SnpEff**: Cingolani et al. (2012) *Fly* 6(2):80-92
- **VariantAnnotation**: Morgan et al., Bioconductor
- **BEAGLE** (if used for phasing): Browning et al. (2021) *Am J Hum Genet* 108(10):1880-1890
