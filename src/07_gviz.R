library(shiny)
library(tidyverse)
library(here)
library(foreach)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

# Create list of bam and vcf files
bam_dirs <- list.dirs(here("data/samples")) %>% 
    grep("samples/HM[0-9]+(_2)?/whatshap$", value = TRUE, x = .)
bam_files <- list.files(
    bam_dirs, pattern = "haplotagged.bam$", full.names = TRUE
)
vcf_files <- list.files(
    here("results/04/filtered_vcf"), 
    pattern = "filteredx2\\.vcf\\.gz$", 
    full.names = TRUE,
    recursive = TRUE
)

# Create list mapping patient ID to sample ID
samples <- list()
samples$sample_ids <- c(
    c("HM24385", "HM26105", "HM27730"),
    c("HM24143", "HM26077"),
    c("HM24631", "HM26107"),
    c("HM23248", "HM20431", "HM23338")
)
samples$subject_ids <- c(
    rep("HG002",3), 
    rep("HG004", 2), 
    rep("HG005", 2), 
    rep("Fibroblast", 3)
)
names(samples$sample_ids) <- samples$subject_ids


# Map subject ID to bam files
foreach(sample = samples$sample_ids, i = seq_along(samples$sample_ids)) %do% {
    pattern = sample
    samples$bam_files[i] <- grep(pattern, bam_files, value = TRUE)
    samples$vcf_files[i] <- grep(pattern, vcf_files, value = TRUE)
}
df <- data.frame(
    subject_id = samples$subject_ids,
    sample_id = samples$sample_ids,
    bam_file = samples$bam_files,
    vcf_file = samples$vcf_files
)
rm(samples)

# Select which samples to analyze
sample_selection <- df %>% filter(subject_id %in% "HG002") %>% pull(sample_id)
bams <- df %>% filter(sample_id %in% sample_selection) %>% pull(bam_file)

# Select which region to analyze (CYP2D6)
region <- GRanges("chr22", IRanges(42126499, 42130810))

# Read in region of interest from bam files
param <- ScanBamParam(which = region, flag = scanBamFlag(isDuplicate = FALSE))
groups <- lapply(bams, function(x){readGAlignments(x, param = param)})
names(groups) <- sample_selection

# Calculate coverage
groupscov <- lapply(bams, function(x){GenomicAlignments::coverage(x, param = param)})
names(groupscov) <- sample_selection

# Convert to Granges
gr <- lapply(names(groupscov),
 function(x){
    name <- x
    x <- groupscov[[x]]
    x <- as(x, "GRanges")
    x <- subset(x, seqnames %in% "chr22")
    x <- as.data.frame(x)
    colnames(x)[colnames(x) == "score"] <- name
    return(x)
 }
)

gr <- Reduce(
    function(x, y){
        full_join(x, y, by = c("seqnames", "start", "end", "width", "strand"))
    }, 
    gr, 
    accumulate=FALSE
) %>%
makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

dTrack <- DataTrack(
    range = gr,
    genome = "hg38",
    chromosome = "chr22",
    name = "Coverage",
    type = "l",
    window = -1
)

# Set color for groups
col <- c("red", "blue", "yellow")
names(col) <- sample_selection



plotTracks(dTrack, groups = sample_selection, from = 42125498, to = 42131810)
