{
    suppressMessages(library("stringr"))
    suppressMessages(library("DESeq2"))
    suppressMessages(library("tximport"))
    suppressMessages(library("jsonlite"))
    suppressMessages(library("glue"))  
    suppressMessages(library("DESeq2"))  
}
{
    WORKDIR <- dirname(dirname(dirname(getwd())))
    subject_info_table_1 <- glue("{WORKDIR}/Resources/Data/RNA/COVID19_master_table_20231007.Methyl_Overlap.20240402.Sex_M.txt")
    subject_info_table_2 <- glue("{WORKDIR}/Resources/Data/RNA/COVID19_master_table_20231007.Methyl_Overlap.20240402.Sex_F.txt")
    id_column <- "Project_ID_Alias"
    filein <- glue("{WORKDIR}/Results/3_rsem/ConfirmedRecovered/[id_column]/[id_column].genes.results")
    fileout_raw <- glue("{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.rawcount.tsv")
    fileout_norm <- glue("{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.normcount.tsv")
    fileout_vst <- glue("{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.normcount.vst.tsv")
    fileout_geoMeans <-  glue("{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.geoMeans.tsv")
    fileout_sizeFactors <-  glue("{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.sizeFactors.tsv")
    design <- "~Sample_Sex+Sample_Age"
    compare_var <- "Sample_Sex"
    meta_columns <- "Project_ID_Alias,Sample_Sex,Sample_Age,Sample_Trait,Sequencing_Type_Platform,Visit_order"
    control_group <- "M"
    case_group <- "F"
    mincount <- "1"
    num_cores <- "50"
    sep <- ","
}
{
    samplemetadf1 <- read.csv(subject_info_table_1, sep="\t")
    samplemetadf2 <- read.csv(subject_info_table_2, sep="\t")
    samplemetadfall <- rbind(samplemetadf1, samplemetadf2)

    file_exp <- filein
    filepathfun <- function(samplename){
        return(str_replace_all(file_exp, "\\[id_column\\]",samplename))
    }

    files_exp <- unlist(lapply(samplemetadfall[[id_column]], filepathfun))
    files_exp_exists <- files_exp[file.exists(files_exp)]
    files_exp_no <- files_exp[!file.exists(files_exp)]
    print(paste0("Samples not exists:",basename(files_exp_no)))
    txi_rsem <- tximport(files_exp_exists, type="rsem", txIn=FALSE, txOut=FALSE)
    counts <- as.data.frame(txi_rsem$counts)
    samples_files_exp_exists <- basename(dirname(files_exp_exists))[basename(dirname(files_exp_exists)) %in% samplemetadfall[[id_column]]]
    colnames(counts) <- samples_files_exp_exists
    samplemetadfallfilt <- samplemetadfall[samplemetadfall[[id_column]] %in% samples_files_exp_exists, ]
    metainfo <- samplemetadfallfilt[, unlist(str_split(meta_columns, pattern=sep))]

    for (col in colnames(metainfo)){
        unique_vals <- unique(metainfo[[col]])
        if(length(unique_vals) == 2){
            metainfo[[col]] <- factor(metainfo[[col]])
        }
    }
}
{
    dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=metainfo, design=as.formula(design))
}
{   
    ddssize <- estimateSizeFactors(dds)
    ddsdisp <- estimateDispersions(ddssize)
    ddsnorm <- nbinomWaldTest(ddsdisp, maxit=1000)
    vsd <- vst(ddsnorm, blind=FALSE)
}
{
    normcount <- counts(ddsnorm, normalized=TRUE)
    rawcount <- counts(ddsnorm, normalized=FALSE)
    geoMeans <- exp(rowMeans(log(rawcount + 1)))
}
{
    write.table(as.data.frame(round(rawcount)), file=fileout_raw, row.names=TRUE, sep="\t", quote=FALSE)
    write.table(as.data.frame(round(normcount)), file=fileout_norm, row.names=TRUE, sep="\t", quote=FALSE)
}
{
    write.table(as.data.frame(assay(vsd)), file=fileout_vst, row.names=TRUE, sep="\t", quote=FALSE)   
}
{
    write.table(as.data.frame(geoMeans), file=fileout_geoMeans, row.names=TRUE, sep="\t", quote=FALSE)
}
{
    write.table(as.data.frame(ddssize$sizeFactor), file=fileout_sizeFactors, row.names=TRUE, sep="\t", quote=FALSE)     
}
