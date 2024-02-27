library(tidyverse)
setwd("/home/anthony/data/12Feb2024_Concordance_Analysis")
library(here)
i_am(".here")
library(foreach)

#Get vcf files
vcf_files <- list.files(
    here("data/samples"), 
    recursive = TRUE, 
    full.names = TRUE,
    pattern = "GRCh38\\.deepvariant\\.vcf\\.gz$"
)
names(vcf_files) <- str_extract(string = basename(vcf_files), pattern = "(.*)\\.GRCh38", group = 1)

message("filtering vcfs...")
# Filter vcf files for PASS filter
filtered_files <- foreach(sample = names(vcf_files)) %do%
    {
        file <- vcf_files[sample]
        sprintf(
            "bcftools view -f 'PASS,.' -O z %s > %s && bcftools index -t %s",
            file, 
            here("results/01/filtered_vcfs",paste0(sample,".GRCh38.deepvariant.vcf.filtered.gz")),
            here("results/01/filtered_vcfs",paste0(sample,".GRCh38.deepvariant.vcf.filtered.gz"))
        ) %>% 
        system()
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
        files = list.files(here("results/01/filtered_vcfs"), pattern = pattern, full.names = TRUE) %>% grep(., pattern = "\\.gz$", value = TRUE)
        input <- paste(files, collapse = " ")
        sprintf(
            "bcftools merge -O z %s > %s",
            input,
            here("results/01/merged_vcfs",paste0(sample,".GRCh38.deepvariant.vcf.merged.gz"))
        ) %>%
        system()
    }
