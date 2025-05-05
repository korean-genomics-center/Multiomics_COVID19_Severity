{
    suppressMessages(library("stringr"))
    suppressMessages(library("DESeq2"))
    suppressMessages(library("tximport"))
    suppressMessages(library("jsonlite"))
    suppressMessages(library("glue"))  
    suppressMessages(library("DESeq2"))  
}
# {
#     expData <- "/BiO/Access/kyungwhan1998/Infectomics/Results/4_expmtx/HealthyPreVaccination/expression_matrix_genes.results_expected_count.tsv"
#     metaData <- "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/RNA/COVID19_master_table_added_CRF_20231007.txt"
#     design <- "~1"
#     id_column <- "Project_ID_Alias"
#     gene_column <- "ID"
#     geoMean_column <- "geoMeans"
#     filein_geoMeans <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.geoMeans.tsv"
#     fileout <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.HealthyControl.20250304.normcount.tsv"
#     fileout_vst <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.HealthyControl.20250304.normcount.vst.tsv"
#     fileout_geoMeans <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.HealthyControl.20250304.geoMeans.tsv"
#     fileout_sizeFactors <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.HealthyControl.20250304.sizeFactors.tsv"
# }
{
    expData <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.rawcount.tsv"
    metaData <- "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/RNA/COVID19_master_table_added_CRF_20250306.txt"
    design <- "~1"
    id_column <- "Project_ID_Alias"
    gene_column <- "ID"
    geoMean_column <- "geoMeans"
    filein_geoMeans <- NULL
    fileout <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.tsv"
    fileout_vst <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.tsv"
    fileout_geoMeans <- NULL
    fileout_sizeFactors <- NULL
}
{
    metainfo <- read.csv(metaData, sep="\t")
    metasamplenames <- metainfo[[id_column]]
    counts <- read.csv(expData, sep="\t")
    countsamplenamesraw <- names(counts)
    # DO NOT LET SAMPLE ID TO BE NUMERIC!!!
    countsamplenamesnew <- sapply(strsplit(countsamplenamesraw, "\\."), function(x) paste(x, collapse = "-"))
    names(counts) <- countsamplenamesnew
    countsfilt <- counts[, names(counts) %in% metasamplenames]
    roundcountsfilt <- as.data.frame(round(countsfilt))
    idcol <- as.data.frame(counts[[gene_column]])
    if(gene_column %in% names(counts)){
        names(idcol) <- gene_column
        countData <- cbind(idcol, roundcountsfilt)
    }else{
        countData <- roundcountsfilt
    }
    if (!is.null(filein_geoMeans)){
        genes <- rownames(read.csv(filein_geoMeans, sep="\t"))
        countData <- countData[genes%in%countData[["ID"]],]
    }
    if(gene_column %in% names(countData)){
        rownames(countData) <- countData[[gene_column]]
        countData <- countData[, -1]
    }
    colData <- metainfo[metainfo[[id_column]] %in% names(counts),]
}
{
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=as.formula(design))
}
{   
    if (!is.null(filein_geoMeans)){
        geoMeans <- read.csv(filein_geoMeans, sep="\t")[[geoMean_column]]
        ddssize <- estimateSizeFactors(dds, geoMeans=geoMeans)
    }
    ddssize <- estimateSizeFactors(dds)
    ddsdisp <- estimateDispersions(ddssize)
    ddsnorm <- nbinomWaldTest(ddsdisp, maxit=1000)
    vsd <- vst(ddsnorm, blind=FALSE)
    normcount <- counts(ddsnorm, normalized=TRUE)
    rawcount <- counts(ddsnorm, normalized=FALSE)
}
{
    write.table(as.data.frame(round(normcount)), file=fileout, row.names=TRUE, sep="\t", quote=FALSE)

    write.table(as.data.frame(assay(vsd)), file=fileout_vst, row.names=TRUE, sep="\t", quote=FALSE)   

    if (!is.null(fileout_geoMeans)){
        geoMeans <- exp(rowMeans(log(rawcount + 1)))
        write.table(as.data.frame(geoMeans), file=fileout_geoMeans, row.names=TRUE, sep="\t", quote=FALSE)
    }

    if (!is.null(fileout_sizeFactors)){
        write.table(as.data.frame(ddssize$sizeFactor), file=fileout_sizeFactors, row.names=TRUE, sep="\t", quote=FALSE)
    }
}
