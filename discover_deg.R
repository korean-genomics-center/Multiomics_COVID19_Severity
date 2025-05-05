library("glue")
library("stringr")
library("tximport")
library("jsonlite")
library("DESeq2")
library("BiocParallel")

WORKDIR = dirname(dirname(getwd()))
# dir_metatable = glue("{WORKDIR}/Resources/Data/Methylation")
# file_metatable_1 = "metatable_combined_all_firstVisit_WO_C19-C0{sample_loo}.cpg_table_file_path.Mild.tsv"
# file_metatable_2 = "metatable_combined_all_firstVisit_WO_C19-C0{sample_loo}.cpg_table_file_path.Severe.tsv"
# id_column = "Sample_ID"
# meta_columns = "Severity_group,Sample_ID,Sample_Sex"
# sep = ","
# full_model = "~Severity_group+Sample_Sex+hypertension+dyslipidemia+diabetes_mellitus+cancer+EverSmoked5PacksTobacco+MedicationUse"
# reduced_model = "~Sample_Sex+hypertension+dyslipidemia+diabetes_mellitus+cancer+EverSmoked5PacksTobacco+MedicationUse"
# design <- "~Severity_group+Sample_Sex"
# compare_var <- "Severity_group"
# mincount = 20
# case_group = "Severe"
# control_group = "Mild"
# num_cores = 10
id_column = "Project_ID_Alias"
meta_columns = "Severity_Binary,Project_ID_Alias,Sample_Sex"
sep = ","
design <- "~Severity_Binary+Sample_Sex"
compare_var <- "Severity_Binary"
mincount = 20
case_group = "Severe"
control_group = "Mild"
num_cores = 10


# list_sample_loo <- c("01", "08", "11", "12", "13", "18", "42", "46", "54", "58")
list_sample_loo <- c("58")
for (sample_loo in list_sample_loo){
    {
        dir_metatable = glue("{WORKDIR}/Results/5_deg")
        file_metatable_1 = glue("COVID19_Visit1_Mild_removed_loo_C19-C0{sample_loo}_20240327.tsv")
        file_metatable_2 = glue("COVID19_Visit1_Severe_removed_loo_C19-C0{sample_loo}_20240327.tsv")
        filein = glue("{WORKDIR}/Results/3_rsem/Confirmed_Matched_with_Methyl_FU_loss/[id_column]/[id_column].genes.results")
        dir_out = glue("{WORKDIR}/Results/5_deg")
        fileout = file.path(dir_out, glue("Visit1_Severe__Visit1_Mild_removed_loo_C19-C0{sample_loo}_20240327.tsv"))
        subject_info_table_1 = file.path(dir_metatable, file_metatable_1)
        print(subject_info_table_1)
        subject_info_table_2 = file.path(dir_metatable, file_metatable_2)
    }

    {   # get input data
        # get sample metadata
        samplemetadf1 <- read.csv(subject_info_table_1, sep="\t")
        samplemetadf2 <- read.csv(subject_info_table_2, sep="\t")
        samplemetadfall <- rbind(samplemetadf1, samplemetadf2)

        filepathfun <- function(samplename){
            return(str_replace_all(filein, "\\[id_column\\]",samplename))
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

    {   # make DESeq2 object
        # dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=metainfo, design=as.formula(full_model))
        # filter out by mincounts
        # keep <- rowSums(counts(dds)) >= mincount
        # ddsfiltered <- dds[keep,]
        # ddsnorm <- DESeq(ddsfiltered, test="LRT", reduced = as.formula(reduced_model))
        dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=metainfo, design=as.formula(design))
        # filter out by mincounts
        keep <- rowSums(counts(dds)) >= mincount
        ddsfiltered <- dds[keep,]
        condition <- compare_var
        ddsrelevel <- ddsfiltered
        ddsrelevel[[condition]] <- relevel(factor(ddsrelevel[[condition]]), ref=control_group)
        ddssize <- estimateSizeFactors(ddsrelevel)
        ddsdisp <- estimateDispersions(ddssize)
        ddsnorm <- nbinomWaldTest(ddsdisp, maxit=1000)
    }
    {   # save normalized count data (i.e. normalized = raw_count/sizefactor) N.B. dds object only has the raw count in dds$count throughout the analysis thus far 
        normcount <- counts(ddsnorm, normalized=TRUE)
        write.table(as.data.frame(normcount), file=paste0(fileout,".normcount"), row.names=TRUE, sep="\t", quote=FALSE)
         # save dds (deseq2 dataset) object used as an input for downstream visualization duirng EDA (QC)
        saveRDS(ddsnorm, paste0(fileout,".ddsobj"))
        # save vst (variance-stabilization transformation) object used as an input for downstream visualization during EDA (QC)
        vsd <- vst(ddsnorm, blind=FALSE)
        saveRDS(vsd, paste0(fileout,".vstobj"))
    }

    {   
        register(MulticoreParam(num_cores))
        resultDEG <- results(ddsnorm, contrast=c(compare_var, case_group, control_group), pAdjustMethod="BH", parallel=TRUE)
        # resultDEG <- results(ddsnorm, pAdjustMethod="BH")
        resultDEG <- resultDEG[order(resultDEG$log2FoldChange), ]
        logDEG <- resultsNames(resultDEG)
        summaryDEG <- summary(resultDEG)
        write.table(logDEG, file=paste0(fileout,".log"), row.names=TRUE, sep="\t", quote=FALSE, col.names=NA)
        write.table(summaryDEG, file=paste0(fileout,".summary"), row.names=TRUE, sep="\t", quote=FALSE, col.names=NA)
        write.table(data.frame("ID"=rownames(resultDEG),resultDEG), file=fileout, row.names=FALSE,sep="\t", quote=FALSE)
    }
}
