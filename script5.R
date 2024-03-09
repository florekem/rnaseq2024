library(tidyverse)
library(limma)
library(edgeR)
#library(plotly)
#library(gt)
library(DT)

study.design <- read_tsv('studydesign.txt')
myDGEList.filtered.norm <- readRDS('myDGEList.filtered.norm.rds')

group <- factor(study.design$group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# voom > transform to log2-counts per million (logCPM),
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
# top, ale raczej wszystkie. Końcowy wynik DE limma
#myTopHits <- topTable(ebFit,
 #                     adjust ="BH", 
  #                    coef=1, 
   #                   number=40000,
    #                  p.value = 0.05,
     #                 confint = TRUE,
      #                sort.by="logFC")
#myTopHits.df <- myTopHits %>%
 # as_tibble(rownames = "geneID")

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- study.design$sample
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
allGenes <- v.DEGList.filtered.norm$E
allGenes.df <- as_tibble(allGenes, rownames='geneID')
#datatable(diffGenes.df,
 #         extensions = c('KeyTable', "FixedHeader"),
  #        caption = 'Table 1: DEGs in cutaneous leishmaniasis',
   #       options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  #formatRound(columns=c(2:11), digits=2)


results[,1]

?as.dist
?voom
