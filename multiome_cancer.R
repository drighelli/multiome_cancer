# prefetch GEX
# prefetch ATAC
sra <- list.files(recursive=TRUE)[-c(1,2)]
lapply(sra, function(x){ system(paste0("fastq-dump --split-3 --gzip ", x))})

fqs <- list.files(pattern="*[1|2].fastq.gz$", full.names=TRUE)


# fn <- as.data.frame(readxl::read_excel("/Users/inzirio/Desktop/gDrive/works/analysis/multiome_cancer/dataset_desc.xlsx"))
# fn <- fn[c(1:23), c(1:9)]
fn$rep<- gsub("rep","REP", fn$replicate)
fn$cond <- gsub("healthy","HLTY", fn$condition)
fn$cond <- gsub("cancer","CNCR", fn$cond)
fn$sample <- paste(fn$rep, fn$cond, sep="_")
# writexl::write_xlsx(x=fn, path="dataset_desc.xlsx")
fn <- as.data.frame(readxl::read_excel("dataset_desc.xlsx"))
## [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
## Read Type using R1/R2
srs <- unique(fn$SRR)
fn$R1 <- NA
fn$R2 <- NA
for(i in seq_along(srs))
{
    fqsr <- fqs[grep(srs[i], fqs)]
    r12 <- NULL
    for(fq in fqsr)
    {
        rt <- paste0("R", ifelse(strsplit(gsub(".fastq.gz","",basename(fq)),"_")[[1]][2]=="1", "1", "2"))
        nfn <- paste0("./", fn$SRR[i],"_S1_L00",fn$lane[i],"_",rt,"_001.fastq.gz")
        print(nfn)
        if(rt=="R1")
        {
            fn$R1[i] <- nfn
        } else {
            fn$R2[i] <- nfn
        }
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









