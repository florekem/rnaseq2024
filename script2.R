library(edgeR)
library(tidyverse)
# ?cpm # CPM = TPM (CPM for genes, TPM for transcripts)
# ?calcNormFactors # TMM


tx2gene <- readRDS('tx2gene.EnsDb.Hsa.v86.rds')

study.design <- read_tsv('studydesign.txt')

counts <- tx2gene$counts
abundance <- tx2gene$abundance

myDGEList <- DGEList(counts)
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5
myDGEList.filtered.by.min.1.cpm <- myDGEList[keepers,] # ale ciągle DGEList nie jest CPM!!!
myDGEList.filtered.by.min.1.cpm.TMM.norm <- calcNormFactors(myDGEList.filtered.by.min.1.cpm,
                                                            method='TMM') # nadal nie!

saveRDS(myDGEList.filtered.by.min.1.cpm.TMM.norm, 'myDGEList.filtered.norm.rds')


# po tych krokach mam DGEList filtrowane i normalizowane, ale ciągle
# jako zwykłe countsy, nie cmp. Jeśli chcę to narysować, to muszę mięc cpm(log=TRUE)

# Do limmy TMM normalizowane będą jeszcze log2cpm przez voom!
# wszystkie log2 cmp tylko do rysunków i analizy


#log2.cpm.filtered.norm <- cpm(myDGEList.filtered.by.min.1.cpm.TMM.norm, log=TRUE) 
#log2.cpm.filtered.norm.tb <- as_tibble(log2.cpm.filtered.norm, rownames='geneID')
#colnames(log2.cpm.filtered.norm.tb) <- c('geneID', study.design$sample)


#saveRDS(log2.cpm.filtered.norm, 'log2.cpm.filtered.norm.rds')
#saveRDS(log2.cpm.filtered.norm.tb, 'log2.cpm.filtered.norm.tb.rds')




