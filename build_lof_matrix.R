#!/usr/bin/env Rscript
# build_lof_matrix.R
#
# Collapse LoF variants per gene from a SnpEff-annotated VCF and build
# a genotype matrix (rows = samples, columns = genes, values = LoF dosage).
#
# Dosage encoding:
#   Unphased: count of ALT alleles at the worst single variant in the gene
#             e.g. 0/0=0, 0/1=1, 1/1=2, tetraploid 0/1/1/1=3
#   Phased:   number of distinct haplotypes carrying at least one LoF allele
#             e.g. 0|1 + 0|1 (cis) = 1,  0|1 + 1|0 (trans) = 2
#
# Usage:
#   Rscript build_lof_matrix.R \
#       --vcf     input.vcf.gz \
#       --out     matrix.tsv \
#       --impact  HIGH        (default: HIGH; use "" for all impacts)
#       --by      gene        (default: gene; or "transcript")
#       --threads 8           (default: 1; Linux/macOS only)
#
# Requirements:
#   BiocManager::install("VariantAnnotation")

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(parallel)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  return(args[idx + 1])
}

vcf_file   <- get_arg("--vcf",     NULL)
out_file   <- get_arg("--out",     "LoF_genotype_matrix.tsv")
impact_lvl <- get_arg("--impact",  "HIGH")
by_field   <- get_arg("--by",      "gene")
n_threads  <- as.integer(get_arg("--threads", "1"))

if (is.null(vcf_file)) {
  stop(
    "ERROR: --vcf is required.\n",
    "Usage: Rscript build_lof_matrix.R --vcf input.vcf.gz --out matrix.tsv\n"
  )
}

if (.Platform$OS.type == "windows" && n_threads > 1) {
  message("NOTE: parallel forking not supported on Windows, running single-threaded.")
  n_threads <- 1L
}

cat("=== LoF Genotype Matrix Builder ===\n")
cat("VCF      :", vcf_file,   "\n")
cat("Output   :", out_file,   "\n")
cat("Impact   :", impact_lvl, "\n")
cat("Collapse :", by_field,   "\n")
cat("Threads  :", n_threads,  "\n\n")


# ── Progress bar ───────────────────────────────────────────────────────────────
make_progress_bar <- function(total, width = 40) {
  start_time <- proc.time()["elapsed"]
  last_pct   <- -1L

  function(i) {
    pct <- as.integer(100 * i / total)
    if (pct == last_pct && i != total) return(invisible(NULL))
    last_pct <<- pct

    filled  <- as.integer(width * i / total)
    bar     <- paste0(strrep("=", filled),
                      if (filled < width) ">" else "",
                      strrep(" ", max(0, width - filled - 1)))

    elapsed <- proc.time()["elapsed"] - start_time
    eta_str <- if (i == 0 || pct == 0) "  --:--" else {
      eta_sec <- elapsed / i * (total - i)
      sprintf("%3dm%02ds", as.integer(eta_sec %/% 60), as.integer(eta_sec %% 60))
    }
    elapsed_str <- sprintf("%3dm%02ds",
                           as.integer(elapsed %/% 60),
                           as.integer(elapsed %% 60))

    cat(sprintf("\r  [%-*s] %3d%%  %d/%d genes  elapsed %s  ETA %s   ",
                width, bar, pct, i, total, elapsed_str, eta_str))
    if (i == total) cat("\n")
    flush.console()
  }
}


# ── Step 0: Input checks ───────────────────────────────────────────────────────
# Check that the VCF has been SnpEff-annotated and has no multiallelic sites.
# These checks run on the header and a small sample of records — fast even for
# large VCFs since we don't load the whole file here.
cat("[0/5] Checking input VCF...\n")

hdr <- scanVcfHeader(vcf_file)

# Check 1: SnpEff annotation
if (!"ANN" %in% names(info(hdr))) {
  stop(
    "ERROR: No ANN field found in the VCF INFO header.\n",
    "       This script requires SnpEff annotation. Please run SnpEff first:\n",
    "       java -jar snpEff.jar <genome> input.vcf.gz > annotated.vcf\n"
  )
}
cat("      SnpEff ANN field: found\n")

# Check 2: multiallelic sites — scan the first 50000 records
cat("      Checking for multiallelic sites (scanning up to 50000 records)...\n")

vcf_check <- readVcf(vcf_file, genome = "",
                     param = ScanVcfParam(info = NA, geno = NA))
vcf_check <- vcf_check[seq_len(min(50000, nrow(vcf_check))), ]
alt_vals  <- as.character(unlist(alt(vcf_check)))
n_multi   <- sum(grepl(",", alt_vals))

if (n_multi > 0) {
  warning(
    "WARNING: ", n_multi, " multiallelic site(s) detected in the first ",
    nrow(vcf_check), " records.\n",
    "         The dosage calculation assumes biallelic sites (ALT index = 1).\n",
    "         Please split multiallelic sites first:\n",
    "         bcftools norm -m -any input.vcf.gz | bcftools norm -f ref.fa > split.vcf\n",
    call. = FALSE
  )
} else {
  cat("      Multiallelic sites: none found\n")
}
cat("\n")


# ── Step 1: Load VCF ───────────────────────────────────────────────────────────
cat("[1/5] Opening VCF file...\n")

vcf_header <- scanVcfHeader(vcf_file)
samples    <- samples(vcf_header)
cat("      Samples found:", length(samples), "\n")


# ── Step 2: Read genotypes and ANN field ───────────────────────────────────────
cat("[2/5] Reading genotypes and ANN annotations...\n")

svp <- ScanVcfParam(info = "ANN", geno = "GT")
vcf <- readVcf(vcf_file, genome = "", param = svp)

cat("      Variants loaded:", nrow(vcf), "\n")


# ── Step 3: Parse ANN field ────────────────────────────────────────────────────
cat("[3/5] Parsing ANN field (impact =", impact_lvl, ")...\n")

ann_raw <- info(vcf)$ANN

parse_ann <- function(ann_entries, impact_filter = "HIGH", by = "gene") {
  if (all(is.na(ann_entries))) return(character(0))
  genes <- c()
  for (entry in ann_entries) {
    records <- strsplit(entry, ",")[[1]]
    for (rec in records) {
      fields <- strsplit(rec, "\\|")[[1]]
      if (length(fields) < 5) next
      impact     <- trimws(fields[3])
      gene_name  <- trimws(fields[4])
      transcript <- trimws(fields[7])
      if (impact_filter != "" && impact != impact_filter) next
      if (nchar(gene_name) == 0) next
      # SnpEff sometimes joins multiple gene names with "&" when a variant
      # overlaps two genes — split these so each gene gets its own column.
      if (by == "transcript" && nchar(transcript) > 0) {
        genes <- c(genes, strsplit(transcript, "&", fixed = TRUE)[[1]])
      } else {
        genes <- c(genes, strsplit(gene_name,  "&", fixed = TRUE)[[1]])
      }
    }
  }
  return(unique(genes))
}

variant_genes <- lapply(ann_raw, parse_ann,
                        impact_filter = impact_lvl,
                        by = by_field)

has_lof       <- sapply(variant_genes, length) > 0
cat("      LoF variants passing filter:", sum(has_lof), "\n")

vcf_lof       <- vcf[has_lof, ]
variant_genes <- variant_genes[has_lof]


# ── Step 4: Build dosage and phase matrices ────────────────────────────────────
# Parse every GT string once into three flat structures used in Step 5:
#   dosage_mat  [variants × samples]  LoF allele count per call (NA if missing)
#   phased_mat  [variants × samples]  1 if phased ("|"), 0 if unphased ("/")
#   hap_array   [variants × samples × ploidy]  1 if that haplotype carries LoF
cat("[4/5] Pre-computing dosage and phase matrices...\n")

gt_matrix <- geno(vcf_lof)$GT
n_var     <- nrow(gt_matrix)
n_samp    <- ncol(gt_matrix)

first_gt <- gt_matrix[which(!is.na(gt_matrix))[1]]
has_pipe  <- grepl("|", first_gt, fixed = TRUE)
ploidy    <- length(strsplit(first_gt, if (has_pipe) "|" else "/", fixed = TRUE)[[1]])
cat("      Detected ploidy:", ploidy, "\n")

gt_flat    <- as.vector(gt_matrix)
is_phased  <- grepl("|", gt_flat, fixed = TRUE)
is_missing <- is.na(gt_flat) |
              gt_flat %in% c(".", "./.", ".|.",
                             paste(rep(".", ploidy), collapse = "/"))

gt_norm     <- gsub("|", "/", gt_flat, fixed = TRUE)
allele_list <- strsplit(gt_norm, "/", fixed = TRUE)

hap_mat <- matrix(0L, nrow = length(allele_list), ncol = ploidy)
for (h in seq_len(ploidy)) {
  allele_h <- sapply(allele_list, function(a)
    if (length(a) >= h) a[h] else NA_character_)
  hap_mat[, h] <- ifelse(is.na(allele_h) | allele_h %in% c(".", "0"), 0L, 1L)
}

dosage_flat             <- rowSums(hap_mat)
dosage_flat[is_missing] <- NA_integer_

dosage_mat <- matrix(dosage_flat,          nrow = n_var, ncol = n_samp,
                     dimnames = dimnames(gt_matrix))
phased_mat <- matrix(as.integer(is_phased), nrow = n_var, ncol = n_samp,
                     dimnames = dimnames(gt_matrix))
hap_array  <- array(hap_mat,
                    dim      = c(n_var, n_samp, ploidy),
                    dimnames = list(rownames(gt_matrix),
                                   colnames(gt_matrix),
                                   paste0("hap", seq_len(ploidy))))

cat("      Dosage matrix built:", n_var, "variants ×", n_samp, "samples\n")


# ── Step 5: Collapse per gene ──────────────────────────────────────────────────
cat("[5/5] Building gene-level LoF matrix...\n")

all_genes <- sort(unique(unlist(variant_genes)))
n_genes   <- length(all_genes)
cat("      Unique genes/transcripts:", n_genes, "\n\n")

gene_to_rows <- setNames(vector("list", n_genes), all_genes)
for (i in seq_along(variant_genes)) {
  for (g in variant_genes[[i]]) {
    gene_to_rows[[g]] <- c(gene_to_rows[[g]], i)
  }
}

process_gene <- function(gene) {

  rows  <- gene_to_rows[[gene]]
  d_sub <- dosage_mat[rows, , drop = FALSE]
  p_sub <- phased_mat[rows, , drop = FALSE]

  all_phased_samp  <- colSums(p_sub, na.rm = TRUE) == nrow(p_sub)
  all_missing_samp <- colSums(!is.na(d_sub)) == 0

  gene_dosage <- integer(n_samp)
  ph <- 0L; unph <- 0L; miss <- 0L

  # Unphased: take the max dosage observed across all LoF variants in the gene
  unph_idx <- which(!all_phased_samp & !all_missing_samp)
  if (length(unph_idx) > 0) {
    d_unph <- d_sub[, unph_idx, drop = FALSE]
    d_unph[is.na(d_unph)] <- 0L
    row_list <- lapply(seq_len(nrow(d_unph)), function(r) d_unph[r, ])
    gene_dosage[unph_idx] <- do.call(pmax, row_list)
    unph <- length(unph_idx)
  }

  # Phased: count how many haplotypes carry at least one LoF allele.
  # Reshape [variants × samples × ploidy] → 2D, use colSums to OR across
  # variants, then reshape back and rowSums to count affected haplotypes.
  ph_idx <- which(all_phased_samp & !all_missing_samp)
  if (length(ph_idx) > 0) {
    h_sub   <- hap_array[rows, ph_idx, , drop = FALSE]
    n_ph    <- length(ph_idx)
    h_2d    <- matrix(h_sub, nrow = length(rows), ncol = n_ph * ploidy)
    hap_any <- matrix(colSums(h_2d) > 0, nrow = n_ph, ncol = ploidy)
    gene_dosage[ph_idx] <- rowSums(hap_any)
    ph <- n_ph
  }

  miss <- sum(all_missing_samp)

  list(dosage = gene_dosage, ph = ph, unph = unph, miss = miss)
}


if (n_threads == 1L) {

  progress    <- make_progress_bar(n_genes)
  progress(0)
  lof_matrix  <- matrix(0L, nrow = n_samp, ncol = n_genes,
                        dimnames = list(samples, all_genes))
  path_counts <- c(phased = 0L, unphased = 0L, missing = 0L)

  for (gi in seq_along(all_genes)) {
    res <- process_gene(all_genes[gi])
    lof_matrix[, gi]        <- res$dosage
    path_counts["phased"]   <- path_counts["phased"]   + res$ph
    path_counts["unphased"] <- path_counts["unphased"] + res$unph
    path_counts["missing"]  <- path_counts["missing"]  + res$miss
    progress(gi)
  }

} else {

  cat(sprintf("  Distributing %d genes across %d threads...\n",
              n_genes, n_threads))

  chunk_size  <- ceiling(n_genes / n_threads)
  gene_chunks <- split(all_genes,
                       ceiling(seq_along(all_genes) / chunk_size))

  progress   <- make_progress_bar(n_genes)
  genes_done <- 0L
  progress(0)

  chunk_results <- mclapply(
    seq_along(gene_chunks),
    function(ci) lapply(gene_chunks[[ci]], process_gene),
    mc.cores = n_threads
  )

  lof_matrix  <- matrix(0L, nrow = n_samp, ncol = n_genes,
                        dimnames = list(samples, all_genes))
  path_counts <- c(phased = 0L, unphased = 0L, missing = 0L)

  gi <- 0L
  for (ci in seq_along(chunk_results)) {
    for (res in chunk_results[[ci]]) {
      gi <- gi + 1L
      lof_matrix[, gi]        <- res$dosage
      path_counts["phased"]   <- path_counts["phased"]   + res$ph
      path_counts["unphased"] <- path_counts["unphased"] + res$unph
      path_counts["missing"]  <- path_counts["missing"]  + res$miss
    }
    genes_done <- genes_done + length(chunk_results[[ci]])
    progress(genes_done)
  }
}

cat("      Matrix dimensions:", nrow(lof_matrix), "samples ×",
    ncol(lof_matrix), "genes\n")


# ── Write output ───────────────────────────────────────────────────────────────
cat("\n[Done] Writing matrix to:", out_file, "\n")

out_df <- as.data.frame(lof_matrix)
out_df <- cbind(Accession = rownames(out_df), out_df)
write.table(out_df, file = out_file, sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\nPhasing path summary (gene × sample cells):\n")
cat("  Phased   :", path_counts["phased"],   "\n")
cat("  Unphased :", path_counts["unphased"], "\n")
cat("  Missing  :", path_counts["missing"],  "\n")

cat("\nDosage distribution:\n")
print(table(as.vector(lof_matrix)))

cat("\nPreview (first 5 samples × first 5 genes):\n")
print(lof_matrix[1:min(5, nrow(lof_matrix)), 1:min(5, ncol(lof_matrix))])

cat("\nDone.\n")