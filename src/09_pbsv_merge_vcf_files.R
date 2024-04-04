library(tidyverse)
library(here)
library(foreach)
library(doParallel)
md <- readRDS(here("data/metadata/md.rds")) %>% select(sample_id, subject_id, bam_file)
registerDoParallel(cores=4)
outdir <- here("results/09")
dir.create(outdir, showWarnings = FALSE)

sample_dirs <- list.dirs(here("data/samples"), recursive = TRUE) %>% 
    grep("HM[0-9]*(_2)?/pbsv/.*\\.chrom_vcfs$",.,value = TRUE)

# For each subject
foreach(subject = unique(md$subject_id)) %dopar% {
    message("processing ", subject, "...")
    dir.create(here("results/09",subject), showWarnings = FALSE)
    

    # For each sample
    samples <- md %>% filter(subject_id == subject) %>% pull(sample_id)
    foreach(sample = samples) %do% {
        dir.create(here("results/09", subject, sample), showWarnings = FALSE)
        dir.create(here("results/09", subject, sample, "fixed_vcf"), showWarnings = FALSE)
        sample_dir <- sample_dirs %>% grep(sample, ., value=TRUE)
        bam <- md %>% filter(sample_id == sample) %>% pull(bam_file)
        vcfs <- list.files(sample_dir, full.names = TRUE, pattern = ".*gz$")
        names(vcfs) <- str_extract(vcfs, pattern = "\\.chr([0-9]+|X|Y|M)\\.", group = 1)


        # Fix headers in genomic order for GatherVcfs to work
        message("fixing headers to genomic order ", sample, "...")
        foreach(vcf = vcfs, chr = names(vcfs)) %do% {
            output <- here("results/09", subject, sample, "fixed_vcf", paste0(sample, ".chr", chr,".fixed.vcf.gz"))
            sprintf(
                "gatk UpdateVCFSequenceDictionary -V %s --source-dictionary %s --output %s --replace true --verbosity ERROR",
                vcf, bam, output
            ) %>% system()
        }
    }
}

# For each subject
foreach(subject = unique(md$subject_id)) %dopar% {
    message("processing ", subject, "...")
    

    # For each sample
    samples <- md %>% filter(subject_id == subject) %>% pull(sample_id)
    foreach(sample = samples) %do% {
        # Gather chr vcfs into one vcf
        message("\n##### gathering vcfs for sample ", sample, "...")
        vcfs <- list.files(
            here("results/09", subject, sample, "fixed_vcf"), 
            pattern = ".*\\.gz$", full.names = TRUE
        )
        names(vcfs) <- str_extract(vcfs, pattern = "\\.chr([0-9]+|X|Y|M)\\.", group = 1)
        chr_order <- c(as.character(1:22), "X", "Y", "M")
        vcfs <- vcfs[chr_order] %>% na.omit()
        vcfs <- paste("-I",vcfs,collapse = " ")
        output <- here("results/09", subject, sample, paste0(sample, ".merged.vcf.gz"))
        sprintf(
            "gatk GatherVcfs %s -O %s --REORDER_INPUT_BY_FIRST_VARIANT true",
            vcfs, output
        ) %>% system()

        # Filter for PASS using bcftools and sort filtered vcf
        input <- output
        output <- here("results/09", subject, sample, paste0(sample, ".merged.filtered.vcf.gz"))
        output2 <- here("results/09", subject, sample, paste0(sample, ".merged.sorted.vcf.gz"))
        sprintf(
            "bcftools view -f 'PASS,.' -O z %s > %s && bcftools sort -O z %s > %s",
            input, output, output, output2
        ) %>% system()

        # Index sorted vcf using bcftools
        input <- output2
        sprintf(
            "bcftools index -t %s",
            input
        ) %>% system()
    }
}

# For each subject
foreach(subject = unique(md$subject_id)) %dopar% {
    message("processing ", subject, "...")
    
    output <- here("results/09", subject, paste0(subject,".merged.vcf.gz"))
    # Merge vcfs using bcftools merge
    vcfs <- list.files(
        here("results/09",subject), 
        full.names = TRUE, 
        recursive = TRUE, 
        pattern = ".*\\.merged\\.sorted\\.vcf\\.gz$"
    )
    vcfs <- paste(vcfs, collapse = " ")

    sprintf(
        "bcftools merge -O z %s > %s",
        vcfs, output
    ) %>% system()
}