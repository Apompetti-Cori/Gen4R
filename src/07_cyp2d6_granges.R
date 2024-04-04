library(shiny)
library(tidyverse)
library(here)
library(foreach)
library(Rsamtools)
library(rtracklayer)
library(Gviz)
library(Hmisc)

# Create list of depth files
dp_files <- list.files(
    here("results/05"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "\\.depth"
)


# Create list of bam files
bam_dirs <- list.dirs(here("data/samples")) %>% 
    grep("samples/HM[0-9]+(_2)?/whatshap$", value = TRUE, x = .)

bam_files <- list.files(
    bam_dirs, pattern = "haplotagged.bam$", full.names = TRUE
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
    samples$dp_files[i] <- grep(pattern, dp_files, value = TRUE)
    samples$bam_files[i] <- grep(pattern, bam_files, value = TRUE)
}

# Create DF mapping subject, sample, and files
df <- data.frame(
    subject_id = samples$subject_ids,
    sample_id = samples$sample_ids,
    dp_file = samples$dp_files,
    bam_file = samples$bam_files
)

df2 <- vroom::vroom(
    file = df$dp_file, 
    col_names = FALSE, 
    show_col_types = FALSE, 
    id = "file"
) %>%
mutate(
    file = str_extract(file, pattern = "/(HM.*)\\.depth", group = 1)
)
colnames(df2) <- c("id", "seqname", "start", "score")

df2 <- df2 %>%
mutate(end = start) %>%
select(id, seqname, start, end, score)

# Create granges with depth for every bp
gr <-  df2 %>% 
pivot_wider(names_from = "id", values_from = "score") %>%
makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Create granges averaging depth for every 100 bp
gr_binned <- df2 %>% 
pivot_wider(names_from = "id", values_from = "score")
cuts <- seq(min(gr_binned$start), max(gr_binned$start), 100)
gr_binned <- gr_binned %>%
mutate(
    grp = cut2(gr_binned$start, cuts) %>% gsub("\\[|\\]|\\(|\\)", "", .)
) %>%
select(-start, -end) %>%
group_by(grp, seqname) %>%
nest() %>%
mutate(
    data = map(
        .x = data, .f = function(x){
            x <- x %>%
                summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
            return(x)
        }
    )
) %>%
unnest(data) %>%
separate(., col = grp, into = c("start", "end"), sep = ",") %>%
mutate(
    start = as.numeric(start),
    end = as.numeric(end)
) %>%
makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(gr, here("results/07/cyp2d6.rds"))
saveRDS(gr_binned, here("results/07/cyp2d6_binned.rds"))
saveRDS(df, here("data/metadata/md.rds"))
