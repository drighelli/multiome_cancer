# prefetch GEX
# prefetch ATAC
sra <- list.files(recursive=FALSE)[-c(1:5)]
lapply(sra, function(x){ system(paste0("fastq-dump --split-files --gzip ", x, " -O fastqs"))})

fqs <- list.files(path="fastqs", pattern="*[1|2].fastq.gz$", full.names=TRUE)
srr <- unlist(lapply(strsplit(fqs,split="_"), function(x) {x[1]}))
## checking that all SRR have fastqs
length(sra)==length(unique(srr))
## checking how many fastqs each SRR has
table(srr) # 4 each

## workaround for local machine browsing
fn <- as.data.frame(readxl::read_excel("dataset_desc.xlsx"))
srr <- fn$SRR
##################
geotab <- c("mtd","da")
geotab="da"
lapply(seq_along(srr), function(i)
{
    browseURL(paste0("https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=",
                     srr[i],
                     ifelse(geotab=="mtd", "&display=metadata", "&display=data-access")))
})
## Based on the information retrieved from 
# 
# rep3 cancer RNAseq
# mv SRR21960678_1.fastq.gz  SRR21960678_S1_L001_I1_001.fastq.gz
# mv SRR21960678_2.fastq.gz  SRR21960678_S1_L001_R1_001.fastq.gz
# # mv SRR21960678_3.fastq.gz  SRR21960678_S1_L001_R2_001.fastq.gz
# # mv SRR21960678_4.fastq.gz  SRR21960678_S1_L001_R3_001.fastq.gz
# mv SRR21960678_S1_L001_R2_001.fastq.gz  SRR21960678_S1_L001_I2_001.fastq.gz
# mv SRR21960678_S1_L001_R3_001.fastq.gz  SRR21960678_S1_L001_R2_001.fastq.gz
# 
# mv SRR22482666_1.fastq.gz SRR22482666_S1_L001_I1_001.fastq.gz
# mv SRR22482666_2.fastq.gz  SRR22482666_S1_L001_R1_001.fastq.gz
# mv SRR22482666_3.fastq.gz  SRR22482666_S1_L001_R2_001.fastq.gz
# mv SRR22482666_4.fastq.gz  SRR22482666_S1_L001_R3_001.fastq.gz
# 
# mv SRR22482667_1.fastq.gz SRR22482667_S1_L001_I1_001.fastq.gz
# mv SRR22482667_2.fastq.gz  SRR22482667_S1_L001_R1_001.fastq.gz
# mv SRR22482667_3.fastq.gz  SRR22482667_S1_L001_R2_001.fastq.gz
# mv SRR22482667_4.fastq.gz  SRR22482667_S1_L001_R3_001.fastq.gz
# ################
# 
# ##Rep1 healthy
# mv SRR21960687_1.fastq.gz SRR21960687_S1_L001_I1_001.fastq.gz
# mv SRR21960687_2.fastq.gz  SRR21960687_S1_L001_R1_001.fastq.gz
# mv SRR21960687_3.fastq.gz  SRR21960687_S1_L001_I2_001.fastq.gz
# mv SRR21960687_4.fastq.gz  SRR21960687_S1_L001_R2_001.fastq.gz
# 
# mv SRR21960688_1.fastq.gz SRR21960688_S1_L001_I1_001.fastq.gz
# mv SRR21960688_2.fastq.gz  SRR21960688_S1_L001_R1_001.fastq.gz
# mv SRR21960688_3.fastq.gz  SRR21960688_S1_L001_I2_001.fastq.gz
# mv SRR21960688_4.fastq.gz  SRR21960688_S1_L001_R2_001.fastq.gz
# 
# 
# mv SRR22482677_1.fastq.gz SRR22482677_S1_L001_I1_001.fastq.gz
# mv SRR22482677_2.fastq.gz  SRR22482677_S1_L001_R1_001.fastq.gz
# mv SRR22482677_3.fastq.gz  SRR22482677_S1_L001_R2_001.fastq.gz
# mv SRR22482677_4.fastq.gz  SRR22482677_S1_L001_R3_001.fastq.gz
# 
# mv SRR22482678_1.fastq.gz SRR22482678_S1_L001_I1_001.fastq.gz
# mv SRR22482678_2.fastq.gz  SRR22482678_S1_L001_R1_001.fastq.gz
# mv SRR22482678_3.fastq.gz  SRR22482678_S1_L001_R2_001.fastq.gz
# mv SRR22482678_4.fastq.gz  SRR22482678_S1_L001_R3_001.fastq.gz

###############


# fn <- as.data.frame(readxl::read_excel("/Users/inzirio/Desktop/gDrive/works/analysis/multiome_cancer/dataset_desc.xlsx"))
# fn <- fn[c(1:23), c(1:9)]
# fn$rep<- gsub("rep","REP", fn$replicate)
# fn$cond <- gsub("healthy","HLTY", fn$condition)
# fn$cond <- gsub("cancer","CNCR", fn$cond)
# fn$sample <- paste(fn$rep, fn$cond, sep="_")
# writexl::write_xlsx(x=fn, path="dataset_desc.xlsx")
fn <- as.data.frame(readxl::read_excel("dataset_desc.xlsx"))
## [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
## Read Type using R1/R2/I1/I2 for GEX
srs <- unique(fn$SRR)
srs <- srs[-which(srs%in%c("SRR22482666", "SRR22482667", "SRR21960678",
        "SRR21960687", "SRR21960688", "SRR22482677", "SRR22482678"))] ## Already renamed runs

fn <- fn[ -which(fn$SRR %in%  c("SRR22482666", "SRR22482667", "SRR21960678",
                                "SRR21960687", "SRR21960688", "SRR22482677", "SRR22482678")),]
# rename_fastqs <- function(srs, fqs)
# {
#     for(i in seq_along(gexsrs))
#     {
#         fqsr <- fqs[grep(gexsrs[i], fqs)]
#         r12 <- NULL
#         for(fq in fqsr)
#         {
#             idx <- which(fn$SRR==gexsrs)
#             switch(strsplit(gsub(".fastq.gz","",basename(fq)),"_")[[1]][2],
#                    "1"={
#                        rt <- "I1"
#                    },
#                    "2"={rt <- "R1"},
#                    "3"={rt <- "I2"},
#                    "4"={rt <- "R2"},
#             )
#             nfn <- paste0("./", gexsrs[i],"_S1_L001_",rt,"_001.fastq.gz")
#             message("renaming ", fq, " in ", nfn)
#             # file.rename(from=fq, to=nfn)
#             fn[[rt]][idx] <- paste0(basename(fq), "->", nfn)
#         }
#     }
# }


## We refer to this reference for the renaming of the GEX/ATAC fastq files
## https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/fastq-input
## https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
## 

fn$F1 <- fn$F2 <- fn$F3 <- fn$F4 <- NA

gexsrs <- srs[fn$omic=="scRNAseq"]
for(i in seq_along(gexsrs))
{
    fqsr <- fqs[grep(gexsrs[i], fqs)]
    for(fq in fqsr)
    {
        idx <- which(fn$SRR %in% gexsrs[i])
        switch(strsplit(gsub(".fastq.gz","",basename(fq)),"_")[[1]][2],
               "1"={rt <- "I1"; ff <- "F1"},
               "2"={rt <- "R1"; ff <- "F2"},
               "3"={rt <- "I2"; ff <- "F3"},
               "4"={rt <- "R2"; ff <- "F4"}
        )
        nfn <- paste0(gexsrs[i],"_S1_L001_",rt,"_001.fastq.gz")
        message("renaming ", fq, " in ", nfn)
        file.rename(from=file.path("fastqs", fq), to=file.path("fastqs",nfn))
        fn[[ff]][idx] <- paste0(basename(fq), " -> ", basename(nfn))
    }
}


atacsrs <- srs[fn$omic=="scATACseq"]
for(i in seq_along(atacsrs))
{
    fqsr <- fqs[grep(atacsrs[i], fqs)]
    for(fq in fqsr)
    {
        idx <- which(fn$SRR %in% atacsrs[i])
        switch(strsplit(gsub(".fastq.gz","",basename(fq)),"_")[[1]][2],
               "1"={rt <- "I1"; ff <- "F1"},
               "2"={rt <- "R1"; ff <- "F2"},
               "3"={rt <- "R2"; ff <- "F3"},
               "4"={rt <- "R3"; ff <- "F4"}
        )
        nfn <- paste0("./", atacsrs[i],"_S1_L001_",rt,"_001.fastq.gz")
        message("renaming ", fq, " in ", nfn)
        fn[[ff]][idx] <- paste0(basename(fq), " -> ", basename(nfn))
        file.rename(from=file.path("fastqs", fq), to=file.path("fastqs",nfn))
    }
}




.importdesign <- function(exp_design)
{
    ext <- strsplit(basename(exp_design), split="[.]")[[1]][2]
    switch (ext,
        "xlsx" = {
            return(as.data.frame(readxl::read_excel("dataset_desc.xlsx")))
        }
    )
}


makemolibraries("dataset_desc.xlsx", proj_folder="/mnt/europa/drighelli/multiome_cancer")


makemolibraries <- function(exp_design, proj_folder=".", fqs_col_ids="SRR", 
    sample_col="sample", library="omic", outfolder="libraries", verbose=TRUE)
{

    des <- .importdesign(exp_design)
    stopifnot(all(c(sample_col, fqs_cols, library) %in% colnames(des)))
    ss <- unique(des[[sample_col]])
    for(s in ss)
    {
        idx <- grep(s, des[[sample_col]])
        # des[fqs_cols][idx,]
        libs <- des[idx, library]
        libs <- ifelse(libs=="scRNAseq", "Gene Expression", "Chromatin Accessibility")
        df <- data.frame("fastqs"=proj_folder,
                         "sample"=des[idx,]$SRR,
                         "library_type"=libs)
        filename <- file.path(outfolder, paste0(s, ".csv"))
        if(verbose) message("Writing ", filename)
        write.csv(df, file=filename, quote=FALSE, row.names=FALSE)
    }
}


## note that replicate 3 cancer isn't working so we make only the first two replicates
## make ALL_12





