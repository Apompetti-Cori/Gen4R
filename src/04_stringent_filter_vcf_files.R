library(tidyverse)
setwd("/home/anthony/data/12Feb2024_Concordance_Analysis")
library(here)
i_am(".here")
library(foreach)
dir.create(here("results/04"), showWarnings = FALSE)
dir.create(here("results/04/filtered_vcf"), showWarnings = FALSE)
dir.create(here("results/04/merged_vcf"), showWarnings = FALSE)

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
message("filtering vcfs...")
# Filter samples for GQ > 20
foreach(subject = names(samples)) %do% {
    message("processing ", subject, "...")
    dir.create(here("results/04/filtered_vcf",subject))
    foreach(file = samples[[subject]]$files, sample_id = samples[[subject]]$sample_ids) %do% {
        output <- here("results/04/filtered_vcf", subject, paste0(sample_id, ".GRCh38.filteredx2.vcf.gz"))
        sprintf(
            "bcftools view -i 'FMT/GQ>=20' -O z %s > %s && bcftools index -t %s",
            file, output, output
        ) %>%
        system()
    }
}

# Create list mapping patient ID to sample ID
samples <- list()
samples[["HG002"]] <- c("HM24385", "HM26105", "HM27730")
samples[["HG004"]] <- c("HM24143", "HM26077")
samples[["HG005"]] <- c("HM24631", "HM26107")
samples[["Fibroblast"]] <- c("HM23248", "HM20431", "HM23338")

message("merging vcfs...")
# Merge vcfs per patient
foreach(sample = names(samples)) %do%
    {
        pattern = paste(samples[[sample]], collapse = "|")
        files = list.files(here("results/04/filtered_vcf"), pattern = pattern, full.names = TRUE) %>% grep(., pattern = "\\.gz$", value = TRUE)
        input <- paste(files, collapse = " ")
        sprintf(
            "bcftools merge -O z %s > %s",
            input,
            here("results/04/merged_vcf",paste0(sample,".GRCh38.filteredx2.vcf.merged.gz"))
        ) %>%
        system()
    }