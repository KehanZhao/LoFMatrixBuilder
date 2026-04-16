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

The script works with both phased and unphased VCFs. If you intend to phase your VCF, it is recommended to do so **before** splitting multiallelic sites, as phasing tools generally expect the original multiallelic representation. Split multiallelic sites after phasing.

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

For each gene, the dosage is the **number of distinct haplotypes carrying at least one LoF allele** across all LoF variants. This correctly distinguishes cis stacking (both LoF alleles on the same haplotype) from trans compound heterozygotes (LoF alleles on different haplotypes). The maximum value is the sample's ploidy.

| Example | Ploidy | Dosage | Interpretation |
|---|---|---|---|
| `0\|1` only | 2 | 1 | one haplotype affected |
| `0\|1` + `0\|1` (cis) | 2 | 1 | both variants on same haplotype |
| `0\|1` + `1\|0` (trans) | 2 | 2 | compound heterozygote |
| `1\|1` | 2 | 2 | both haplotypes affected |
| `0\|1\|0\|1` | 4 | 2 | two haplotypes affected (tetraploid) |
| `0\|1\|1\|1` | 4 | 3 | three haplotypes affected (tetraploid) |
| `1\|1\|1\|1` | 4 | 4 | all haplotypes affected (tetraploid) |

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

- **SnpEff**: Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., … Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. *Fly*, 6(2), 80–92. https://doi.org/10.4161/fly.19695

- **BEAGLE** (if used for phasing): Browning, B. L., Tian, X., Zhou, Y., & Browning, S. R. (2021). Fast two-stage phasing of large-scale sequence data. *The American Journal of Human Genetics*, 108(10), 1880–1890. https://doi.org/10.1016/j.ajhg.2021.08.005
