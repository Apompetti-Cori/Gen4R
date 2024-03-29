library(tidyverse)
library(here)
library(foreach)
library(igvR)
dir.create(here("results/06"))

bam_dirs <- list.dirs(here("data/samples")) %>% 
    grep("samples/HM[0-9]+(_2)?/whatshap$", value = TRUE, x = .)

bam_files <- list.files(
    bam_dirs, pattern = "haplotagged.bam$", full.names = TRUE
)

vcf_files <- list.files(
    here("results/01/filtered_vcfs"), pattern = "filtered.gz$", full.names = TRUE
)


# Create list mapping patient ID to sample ID
samples <- list()
samples[["HG002"]]$sample_ids <- c("HM24385", "HM26105", "HM27730")
samples[["HG004"]]$sample_ids <- c("HM24143", "HM26077")
samples[["HG005"]]$sample_ids <- c("HM24631", "HM26107")
samples[["Fibroblast"]]$sample_ids <- c("HM23248", "HM20431", "HM23338")

# Map subject ID to bam files
foreach(subject = names(samples)) %do% {
    pattern = paste(samples[[subject]]$sample_ids, collapse = "|")
    samples[[subject]]$bam_files <- grep(pattern, bam_files, value = TRUE)
    samples[[subject]]$vcf_files <- grep(pattern, vcf_files, value = TRUE)
}

igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "CYP2D6")
roi <- getGenomicRegion(igv)
gr.roi <- with(roi, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
param <- ScanBamParam(which = gr.roi, what = scanBamWhat())

foreach(
  subject = names(samples)
) %do% {
    foreach(
        bam_file = samples[[subject]]$bam_files,
        vcf_file = samples[[subject]]$vcf_files,
        sample = samples[[subject]]$sample_ids
    ) %do% {

        # Display BAM Track
        alignments <- readGAlignments(bam_file, use.names = TRUE, param = param)
        atrack <- GenomicAlignmentTrack(
            trackName=paste0(sample,".bam"), 
            alignments, 
            visibilityWindow=10000000, 
            trackHeight=260
        )
        displayTrack(igv, atrack)

        # Display VCF Track
        vcf.sub <- VariantAnnotation::readVcf(vcf_file, "hg38", param=gr.roi)
        vtrack <- VariantTrack(
            trackName = paste0(sample,".vcf"),
            vcf = vcf.sub,
            trackHeight=20
        )
        displayTrack(igv, vtrack)
    }
}
