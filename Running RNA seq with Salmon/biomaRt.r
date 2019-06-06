source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

setwd('C:\\Users\\Justi\\Desktop\\Classes Fall 2018\\EQTL ANALYSIS ALL\\EQTL_SPENCER_Analysis\\matrixEQTL_redo')

#read in data frame with transcript names to annotate
df = read.csv('liang_ncrna_061118_TPM_Transcripts.csv')

#extract column with transcript names and rename for appropriate biomart search
df_transcript <- data.frame(unique(df$Name))
colnames(df_transcript)[1] <- 'ensembl_transcript_id_version'
colnames(df)[1] <- 'ensembl_transcript_id_version'



#Prepare Biomart
#use homosapien mart
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")


#place search ids in list in this variable.  in this case the transcript ids
affyids <-as.list(as.data.frame(t(df_transcript$ensembl_transcript_id_version)))

#query biomaRt where attributes are the values we want
results <- getBM(attributes = c("ensembl_transcript_id_version","chromosome_name", "transcript_start", "transcript_end"),
                 filters = "ensembl_transcript_id_version", 
                 values = affyids,
                 mart = mart)


#merge results to data original data frame here if desired
merged<-merge(results,df, by='ensembl_transcript_id_version', all=TRUE)


write.csv(merged, file='C:\\Users\\Justi\\Desktop\\MatrixEQTL_redo\\liang_ncrna_061118_TPM_Transcripts_annotated.csv')
write.csv(results, file='C:\\Users\\Justi\\Desktop\\MatrixEQTL_redo\\Transcript_location.csv')
