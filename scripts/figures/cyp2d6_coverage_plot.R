library(here)
library(tidyverse)
library(foreach)
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
gr <- readRDS(here("results/07/cyp2d6.rds"))
gr_binned <- readRDS(here("results/07/cyp2d6_binned.rds"))
md <- readRDS(here("data/metadata/md.rds"))
rownames(md) <- md$sample_id
outdir <- here("output/figures",paste0(Sys.Date(),"_cyp2d6_plot"))
dir.create(outdir, showWarnings = FALSE)

from <- 42126500 - 1000
to <- 42130791 + 1000

# Create genome axis track and ideogram track
axTrack <- GenomeAxisTrack()
idTrack <- IdeogramTrack(genome = "hg38", "chr22")

# Create sequence track for hg38
sTrack <- SequenceTrack(Hsapiens, add53 = TRUE, chromosome = "chr22")

# Create gene region track for CYP2D6
grTrack <- GeneRegionTrack(
  here("data/metadata/gtf/gencode.v45.annotation.gtf"),
  collapseTranscripts = "meta", transcriptAnnotation = "symbol",
  genome = "hg38", chromosome = "chr22",
  from = from, to = to
)
grTrack <- grTrack[feature(grTrack) %in% c("CDS", "UTR") & symbol(grTrack) == "CYP2D6"]

# Create alignment tracks for all samples
alTracks <- foreach(sample = md$sample_id) %do% {
    bam <- md %>% filter(sample_id == sample) %>% pull(bam_file)
    subject <- md %>% filter(sample_id == sample) %>% pull(subject_id)
    track <- AlignmentsTrack(
        bam, name = paste(subject, sample, sep = "\n"),
        genome = "hg38", chromosome = "chr22",
        from = from, to = to,
        showIndels = TRUE, isPaired = TRUE,
        type = "coverage"
    )

    return(track)
}

# Create data tracks for all samples using the depth files
dTracks <- foreach(subject = unique(md$subject_id)) %do% {
  grp <- md %>% filter(subject_id == subject) %>% pull(sample_id)
  data <- gr_binned[,grp]
  track <- DataTrack(
    data, name = subject,
    genome = "hg38", chromosome = "chr22",
    from = from, to = to,
    size = 2, type = c("a"),
    groups = grp, ylim = c(0,35)
  )
  return(track)
}

# Combine all the tracks for plotting
tracks <- c(idTrack, axTrack, grTrack, alTracks, sTrack)

# Save plot to png
png(filename = here(outdir, "plot.png"), units = "in", width = 10, height = 15, res = 300)
plotTracks(
  tracks,
  from = from, to = to,
  showSampleNames = TRUE
)
dev.off()
