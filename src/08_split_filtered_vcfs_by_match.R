library(here)
library(tidyverse)
library(vcfR)
library(foreach)
library(matrixStats)
library(Hmisc)
library(data.table)
source(here("scripts/functions/functions.R"))
outdir <- here("results/08")
dir.create(outdir, showWarnings = FALSE)

# Get the merged files for input
vcf_files <- list.files(here("results/04/merged_vcf"), pattern = "\\.merged\\.gz", full.names = TRUE)

# Name the list with the sample name
names(vcf_files) <- str_extract(string = vcf_files, pattern = (".*/(.*)\\.GRCh38"), group = 1)

foreach(sample = names(vcf_files)) %do% 
    {
        sprintf("reading in file for %s...", sample) %>% message()
        stats <- list()
        dir.create(here(outdir, sample))
        file =  vcf_files[[sample]]
        vcf <- read.vcfR(file = file)
        stats[["n_matching"]][["total"]] <- nrow(vcf)

        message("extracting GT...")
        gt <- extract.gt(vcf, element = "GT")

        message("reordering GT...")
        # String split gt for each column and order by number. 
        # Reason: 1/4 and 4/1 GT are the same but will not be considered the same.
        key <- unique(na.omit(c(gt)))
        value <- unique(na.omit(c(gt)))
        value <- lapply(value, function(x){
            res <- as.numeric(strsplit(x, split = "/", fixed = TRUE)[[1]])
            res <- sort(res)
            res <- paste(res, collapse = "/")
        }) %>% unlist()
        names(value) <- key

        foreach(key = names(value), .final = message("done")) %do% {
            idx <- which(gt == key)
            gt[idx] <- value[[key]]
            return(NULL)
        }

        sprintf("evaluating samples for %s...", sample) %>% message()
        # DO THIS IF THE VCF HAS 2 SAMPLES
        if(ncol(gt) == 2){
            message("vcf has 2 samples...")
            # Get observations where variant is present in all samples and all match
            idx <- (rowSums((is.na(gt))) == 0) & (gt[,1] == gt[,2])
            idx <- which(idx)
            sub_vcf <- vcf[idx,]

            message("calculating variant count two present two...")
            stats[["n_matching"]][["two_present_two"]] <- nrow(sub_vcf)
            message("calculating variant stats two present two...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["two_present_two"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["two_present_two"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.two_present_two.gz")))

            # Get observations where variant is present in all samples and does not match
            idx <- (rowSums((is.na(gt))) == 0) & (gt[,1] != gt[,2])
            idx <- which(idx)
            sub_vcf <- vcf[idx,]
            message("calculating variant count two present zero...")
            stats[["n_matching"]][["two_present_zero"]] <- nrow(sub_vcf)
            message("calculating variant stats two present zero...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["two_present_zero"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["two_present_zero"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.two_present_zero.gz")))

            # Get observations where variant is present in only 1 sample
            idx <- (rowSums((is.na(gt))) == 1)
            idx <- which(idx)
            sub_vcf <- vcf[idx,]
            message("calculating variant count one present...")
            stats[["n_matching"]][["one_present"]] <- nrow(sub_vcf)
            message("calculating variant stats one present...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["one_present"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["one_present"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.one_present.gz")))

            # Save Summary Stats
            message("savings stats for: ", sample,"...\n")
            saveRDS(stats, here(outdir,sample,paste0(sample,"_stats.rds")))
        }

        # DO THIS IF THE VCF HAS 3 SAMPLES
        if(ncol(gt) == 3){
            message("vcf has 3 samples...")
            # Get observations where variant is present in all samples and all match
            idx <- (rowSums((is.na(gt))) == 0) & (gt[,1] == gt[,2]) & (gt[,1] ==  gt[,3])
            idx <- which(idx)
            sub_vcf <- vcf[idx,]
            message("calculating variant count thee present three...")
            stats[["n_matching"]][["three_present_three"]] <- nrow(sub_vcf)
            message("calculating variant stats three present three...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["three_present_three"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["three_present_three"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.three_present_three.gz")))

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
            message("calculating variant count thee present two...")
            stats[["n_matching"]][["three_present_two"]] <- nrow(sub_vcf)
            message("calculating variant stats three present two...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["three_present_two"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["three_present_two"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.three_present_two.gz")))

            # Get observations where variant is present in all samples and none match
            idx <- (rowSums((is.na(gt))) == 0) & ((gt[,1] != gt[,2]) & (gt[,1] !=  gt[,3]) & (gt[,2] !=  gt[,3]))
            idx <- which(idx)
            ## Subset vcf and gt
            sub_vcf <- vcf[idx,]
            sub_gt <- gt[idx,]
            message("calculating variant count thee present zero...")
            stats[["n_matching"]][["three_present_zero"]] <- nrow(sub_vcf)
            message("calculating variant stats three present zero...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["three_present_zero"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["three_present_zero"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.three_present_zero.gz")))

            # Get observations where variant is present in two samples and match
            idx <- (rowSums((is.na(gt))) == 1) & ((gt[,1] == gt[,2]) | (gt[,1] ==  gt[,3]) | (gt[,2] ==  gt[,3]))
            idx <- which(idx)

            sub_vcf <- vcf[idx,]
            message("calculating variant count two present two...")
            stats[["n_matching"]][["two_present_two"]] <- nrow(sub_vcf)
            message("calculating variant stats two present two...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["two_present_two"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["two_present_two"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.two_present_two.gz")))

            # Get observations where variant is present in two samples and don't match
            idx <- (rowSums((is.na(gt))) == 1) & ((gt[,1] != gt[,2]) | (gt[,1] !=  gt[,3]) | (gt[,2] !=  gt[,3]))
            idx <- which(idx)

            sub_vcf <- vcf[idx,]
            message("calculating variant count two present zero...")
            stats[["n_matching"]][["two_present_zero"]] <- nrow(sub_vcf)
            message("calculating variant stats two present zero...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["two_present_zero"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["two_present_zero"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.two_present_zero.gz")))

            # Get observations where variant is present in one sample
            idx <- (rowSums((is.na(gt))) == 2)
            idx <- which(idx)

            sub_vcf <- vcf[idx,]
            message("calculating variant count one present...")
            stats[["n_matching"]][["one_present"]] <- nrow(sub_vcf)
            message("calculating variant stats one present...")
            dp <- extract.gt(sub_vcf, element = "DP", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["dp_stats"]][["one_present"]] <- dp
            gq <- extract.gt(sub_vcf, element = "GQ", as.numeric = TRUE) %>% 
                data.frame() %>%
                mutate(
                    min = rowMins(as.matrix(.), na.rm = TRUE),
                    max = rowMaxs(as.matrix(.), na.rm = TRUE))
            stats[["gq_stats"]][["one_present"]] <- gq
            write.vcf(sub_vcf, file = here(outdir,sample,paste0(sample,".vcf.one_present.gz")))


            # Save Summary Stats
            message("savings stats for: ", sample,"...\n")
            saveRDS(stats, here(outdir,sample,paste0(sample,"_stats.rds")))
        }
    }