library(tidyverse)
library(vcfR)
setwd("/home/anthony/data/12Feb2024_Concordance_Analysis")
library(here)
i_am(".here")
library(foreach)
library(Hmisc)
library(data.table)
source(here("scripts/functions/functions.R"))
dir.create(here("results/02"))

# Get the merged files for input
vcf_files <- list.files(here("results/01/merged_vcfs"), pattern = "\\.merged\\.gz", full.names = TRUE)

# Name the list with the sample name
names(vcf_files) <- str_extract(string = vcf_files, pattern = (".*/(.*)\\.GRCh38"), group = 1)

foreach(sample = names(vcf_files)) %do% 
    {
        sprintf("reading in file for %s...", sample) %>% message()
        stats <- list("n_matching" = list())
        dir.create(here("results/02", sample))
        file =  vcf_files[[sample]]
        vcf <- read.vcfR(file = file)
        stats[["n_matching"]][["total"]] <- nrow(vcf)

        message("extracting GT...")
        gt <- extract.gt(vcf, element = "GT") %>% as.data.table()

        message("reordering GT...")
        # String split gt for each column and order by number. 
        # Reason: 1/4 and 4/1 GT are the same but will not be considered the same.
        gt <- gt %>% mutate_if(is.character, reorder_gt)

        sprintf("evaluating samples for %s...", sample) %>% message()
        # DO THIS IF THE VCF HAS 2 SAMPLES
        if(ncol(gt) == 2){
            message("vcf has 2 samples...")
            # Get observations where variant is present in all samples and all match
            idx <- (rowSums((is.na(gt))) == 0) & (gt[,1] == gt[,2])
            idx <- which(idx)
            sub_vcf <- vcf[idx,]

            stats[["n_matching"]][["two_present_two"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.two_present_two.gz")))

            # Get observations where variant is present in all samples and does not match
            idx <- (rowSums((is.na(gt))) == 0) & (gt[,1] != gt[,2])
            idx <- which(idx)
            sub_vcf <- vcf[idx,]
            stats[["n_matching"]][["two_present_zero"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.two_present_zero.gz")))

            # Get observations where variant is present in only 1 sample
            idx <- (rowSums((is.na(gt))) == 1)
            idx <- which(idx)
            sub_vcf <- vcf[idx,]
            stats[["n_matching"]][["one_present"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.one_present.gz")))

            # Save Summary Stats
            message("savings stats for: ", sample,"...\n")
            saveRDS(stats, here("results/02",sample,paste0(sample,"_stats.rds")))
            # Save Summary Stats to Tsv
            stats_df <- data.frame(unlist(stats$n_matching))
            colnames(stats_df) <- "n_matching"
            data.table::fwrite(stats_df, file = here("results/02",sample,paste0(sample,"_stats.csv")), bom = TRUE, row.names = TRUE)
        }

        # DO THIS IF THE VCF HAS 3 SAMPLES
        if(ncol(gt) == 3){
            message("vcf has 3 samples...")
            # Get observations where variant is present in all samples and all match
            idx <- (rowSums((is.na(gt))) == 0) & (gt[,1] == gt[,2]) & (gt[,1] ==  gt[,3])
            idx <- which(idx)
            sub_vcf <- vcf[idx,]
            stats[["n_matching"]][["three_present_three"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.three_present_three.gz")))

            # Get observations where variant is present in all samples and two match
            ## First get rows where all present and any match
            idx <- (rowSums((is.na(gt))) == 0) & ((gt[,1] == gt[,2]) | (gt[,1] ==  gt[,3]) | (gt[,2] ==  gt[,3]))
            idx <- which(idx)
            ## Subset vcf and gt
            sub_vcf <- vcf[idx,]
            sub_gt <- gt[idx,]

            ## Second filter out rows where all match
            idx <- !((sub_gt[,1] == sub_gt[,2]) & (sub_gt[,1] == sub_gt[,3]))
            idx <- which(idx)
            ## Subset again
            sub_vcf <- sub_vcf[idx,]
            sub_gt <- sub_gt[idx,]
            stats[["n_matching"]][["three_present_two"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.three_present_two.gz")))

            # Get observations where variant is present in all samples and none match
            idx <- (rowSums((is.na(gt))) == 0) & ((gt[,1] != gt[,2]) & (gt[,1] !=  gt[,3]) & (gt[,2] !=  gt[,3]))
            idx <- which(idx)
            ## Subset vcf and gt
            sub_vcf <- vcf[idx,]
            sub_gt <- gt[idx,]
            stats[["n_matching"]][["three_present_zero"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.three_present_zero.gz")))

            # Get observations where variant is present in two samples and match
            idx <- (rowSums((is.na(gt))) == 1) & ((gt[,1] == gt[,2]) | (gt[,1] ==  gt[,3]) | (gt[,2] ==  gt[,3]))
            idx <- which(idx)

            sub_vcf <- vcf[idx,]
            stats[["n_matching"]][["two_present_two"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.two_present_two.gz")))

            # Get observations where variant is present in two samples and don't match
            idx <- (rowSums((is.na(gt))) == 1) & ((gt[,1] != gt[,2]) | (gt[,1] !=  gt[,3]) | (gt[,2] !=  gt[,3]))
            idx <- which(idx)

            sub_vcf <- vcf[idx,]
            stats[["n_matching"]][["two_present_zero"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.two_present_zero.gz")))

            # Get observations where variant is present in one sample
            idx <- (rowSums((is.na(gt))) == 2)
            idx <- which(idx)

            sub_vcf <- vcf[idx,]
            stats[["n_matching"]][["one_present"]] <- nrow(sub_vcf)
            write.vcf(sub_vcf, file = here("results/02",sample,paste0(sample,".vcf.one_present.gz")))


            # Save Summary Stats
            message("savings stats for: ", sample,"...\n")
            saveRDS(stats, here("results/02",sample,paste0(sample,"_stats.rds")))

            # Save Summary Stats to Tsv
            stats_df <- data.frame(unlist(stats$n_matching))
            colnames(stats_df) <- "n_matching"
            data.table::fwrite(stats_df, file = here("results/02",sample,paste0(sample,"_stats.csv")), bom = TRUE, row.names = TRUE)
        }
    }