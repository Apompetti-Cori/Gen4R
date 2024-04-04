library(tidyverse)
library(here)
library(foreach)
library(doParallel)
library(rtracklayer)
library(GenomicAlignments)
dir.create(here("results/05"), showWarnings=FALSE)
dir.create(here("results/05/bigwig_files"), showWarnings=FALSE)
dir.create(here("results/05/depth_files"), showWarnings=FALSE)
registerDoParallel(cores=4)

df <- data.frame(
  chr = "chr22",
  start = seq(from = 42126400, to = 42130800, by = 100),
  end = seq(from = 42126400, to = 42130800, by = 100) + 99 
)

data.table::fwrite(df, here("results/05/cyp2d6.bed"), sep = "\t", col.names = FALSE)

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

# Calculate coverage for sliding window using samtools depth
foreach(
  subject = names(samples)
) %dopar% {
  message("processing ", subject, "...")
  foreach(
    file = samples[[subject]]$files,
    i = seq_along(samples[[subject]]$files)
  ) %do% {
    bam_file <- file
    bed_file <- here("results/05/cyp2d6.bed")
    sample <- samples[[subject]]$sample_ids[i]
    output <- here("results/05/depth_files",paste0(sample,".depth"))
    message("reading ", sample, " alignments...")
    sprintf(
      "samtools depth -b %s -o %s %s",
      bed_file, output, bam_file
    ) %>% system()
  }
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
    region <- GRanges("chr22:42126499-42130791")
    param <- ScanBamParam(which = region, mapqFilter = 20)
    alignment <- readGAlignments(bam_file, param = param)
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