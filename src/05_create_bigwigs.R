library(tidyverse)
setwd("/home/anthony/data/12Feb2024_Concordance_Analysis")
library(here)
library(foreach)
library(doParallel)
library(GenomicAlignments)
library(trackplot)
library(rtracklayer)
dir.create(here("results/05"))
dir.create(here("results/05/bigwig_files"))
registerDoParallel(cores=4)

bam_dirs <- list.dirs(here("data/samples")) %>% 
    grep("samples/HM[0-9]+(_2)?/whatshap$", value = TRUE, x = .)

bam_files <- list.files(
    bam_dirs, pattern = "haplotagged.bam$", full.names = TRUE
)

# Create list mapping patient ID to sample ID
samples <- list()
samples[["HG002"]]$sample_ids <- c("HM24385", "HM26105", "HM27730")
samples[["HG004"]]$sample_ids <- c("HM24143", "HM26077")
samples[["HG005"]]$sample_ids <- c("HM24631", "HM26107")
samples[["Fibroblast"]]$sample_ids <- c("HM23248", "HM20431", "HM23338")

# Map subject ID to files
foreach(subject = names(samples)) %do% {
    pattern = paste(samples[[subject]]$sample_ids, collapse = "|")
    samples[[subject]]$files <- grep(pattern, bam_files, value = TRUE)
}

# Convert bam files to big wigs
foreach(
  subject = names(samples)
) %dopar% {
  message("processing ", subject, "...")
  foreach(
    file = samples[[subject]]$files,
    i = seq_along(samples[[subject]]$files)
  ) %do% {
    bam_file <- file
    sample <- samples[[subject]]$sample_ids[i]
    message("reading ", sample, " alignments...")
    alignment <- readGAlignments(bam_file)
    reads_coverage <- coverage(alignment)
    message("exporting coverage...")
    export.bw(
        reads_coverage, 
        con = here("results/05/bigwig_files",paste0(sample, ".bw"))
    )
    message(sample, " done... \n")
    NULL
  }
}