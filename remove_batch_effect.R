{
    suppressMessages(library(sva))
}
{
   path_expData <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.rawcount.tsv"
   path_metaData <- "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/RNA/COVID19_master_table_added_CRF_20250306.txt"
   id_column <- "Project_ID_Alias"
   gene_column <- "ID"
   varBatch <- "Sequencing_Type_Platform"
   varGroup <- "Severity_Binary"
   file_out <- "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.rawcount.batch_adj.platform_corrected.severity_preserved.tsv"
}
{   
    counts <- read.csv(path_expData, sep="\t")
    colName <- names(counts)
    # DO NOT LET SAMPLE ID TO BE NUMERIC!!!
    colRename <- sapply(strsplit(colName, "\\."), function(x) paste(x, collapse = "-"))
    names(counts) <- colRename
    metaInfo <- read.csv(path_metaData, sep="\t")
    sampleData <- metaInfo[metaInfo[[id_column]] %in% colRename, ]
    geneName <- counts[[gene_column]]
    rownames(counts) <- geneName
    countData <- counts[, -1]
    sampleData[sampleData["Visit_order"] == "Recover", "Severity_Binary"] <- "Convalescent"
    sampleData[sampleData["Sample_Trait"] == "Healthy", "Severity_Binary"] <- "Healthy"
}
{   
    countRmvBatch <- ComBat_seq(counts = as.matrix(countData), batch=sampleData[[varBatch]], group=sampleData[[varGroup]])
}

{   # save counts matrix batch removed
    write.table(countRmvBatch, file=file_out, row.names=TRUE,sep="\t", quote=FALSE)
}
