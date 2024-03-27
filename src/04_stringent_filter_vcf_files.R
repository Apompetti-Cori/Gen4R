library(tidyverse)
setwd("/home/anthony/data/12Feb2024_Concordance_Analysis")
library(here)
i_am(".here")
library(foreach)
dir.create(here("results/04"))
dir.create(here("results/04/filtered_vcf"))

#Get vcf files
vcf_files <- list.files(
    here("results/01/filtered_vcfs"), 
    recursive = TRUE, 
    full.names = TRUE,
    pattern = "GRCh38\\.deepvariant\\.vcf\\.filtered\\.gz$"
)

# Create list mapping patient ID to sample ID
samples <- list()
samples[["HG002"]]$sample_ids <- c("HM24385", "HM26105", "HM27730")
samples[["HG004"]]$sample_ids <- c("HM24143", "HM26077")
samples[["HG005"]]$sample_ids <- c("HM24631", "HM26107")
samples[["Fibroblast"]]$sample_ids <- c("HM23248", "HM20431", "HM23338")

# Map patient ID to files
foreach(subject = names(samples)) %do% {
    pattern = paste(samples[[subject]]$sample_ids, collapse = "|")

    samples[[subject]]$files <- grep(pattern, vcf_files, value = TRUE)
}

# Filter samples for GQ > 20
foreach(subject = names(samples)) %do% {
    message("processing", subject)
    dir.create(here("results/04/filtered_vcf",subject))
    foreach(file = samples[[subject]]$files, sample_id = samples[[subject]]$sample_ids) %do% {
        output <- here("results/04/filtered_vcf", subject, paste0(sample_id, ".GRCh38.filteredx2.vcf.gz"))
        sprintf(
            "bcftools view -i 'FMT/GQ>=20' %s > %s",
            file, output
        ) %>% 
        system()
    }
}
