library(tidyverse)
library(tximport)
#library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(biomaRt)

#### FILES ####
study.design <- read_tsv('studydesign.txt')
path <- paste('mappedReads/', study.design$sample, '/abundance.h5', sep='')
#path <- file.path(targets$sample, "abundance.tsv")


#### BIOMART ENSEMBL ####
#listMarts()
#myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
#available.datasets <- listDatasets(myMart)

# 277081 obs. najwięcej, ale też dużo tx jest bez nazw genów

hsa.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(hsa.anno)
Tx <- getBM(attributes=c('ensembl_transcript_id_version',
                         'external_gene_name'), mart = hsa.anno)

Tx <- as_tibble(Tx)


#### ENSEMBLDB Hsa v 86 ####
# 216741 obs. wszystkie tx mają nazwy genów

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
#transcrip ID needs to be the first column in the dataframe
Tx <- dplyr::select(Tx, "target_id", "gene_name")

#### GENCODE DB ####
# 117952 obs. i w nazwach tx i genów są ensemblowe kody
# jak poznac nazwy genow?

txdb.gencode <- makeTxDbFromGFF('gencode.v45.basic.annotation.gff3',
                                organism = 'Homo sapiens')
k <- keys(txdb.gencode, keytype = "TXNAME")
tx2gene <- select(txdb.gencode, k, "GENEID", "TXNAME")

columns(txdb.gencode)



#### TXIMPORT ####
tx2gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

saveRDS(tx2gene, 'tx2gene.EnsDb.Hsa.v86.rds')



counts <- tximport$counts
colnames(counts) <- study.design$sample
genes <- rownames(counts)

counts_tb <- as_tibble(counts)
counts_tb[['gene']] <- genes

#######