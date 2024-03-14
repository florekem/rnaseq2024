library(tidyverse)
library(limma)
library(edgeR)
#library(plotly)
#library(gt)
library(DT)
library(RColorBrewer)
library(gplots)
library(gprofiler2)
library(fgsea)

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
myTopHits <- topTable(ebFit,
                      adjust ="BH",
                      coef=1,
                      number=40000,
                      p.value = 0.001,
                      lfc = 4,
                      confint = TRUE,
                      sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

#results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=2)
# Although this function enables users to set 
# p-value and lfc cutoffs simultaneously, this combination criterion is not recommended.
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.001, lfc=4)
colnames(v.DEGList.filtered.norm$E) <- study.design$sample
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# heatmap
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2) # up/down regulated
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))
# select down and up regulated modules
diffGenes_upreg <- v.DEGList.filtered.norm$E[results[,1] ==1,]
diffGenes_downreg <- v.DEGList.filtered.norm$E[results[,1] ==-1,]
#
# gprofiler
gost.res <- gost(rownames(diffsenes_upreg), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res, interactive = T, capped = T)

# GESA
# krok 1 to uzyskanie listy myTopHits z wynikami testów limmy
# gdzie mam log2fc i nazwy genów, co jest inputem dla fgsea.
# pytanie czy powinny to byc wszystkie DGE, czy tylko te, które mają 
# interesujący mnie FC i pval.
# taraz idę tylko z tymi, które mają l2FC > 4 i pval 0.001, ale raczej
# w tym przypadku lepiej, jak tych genów jest chyba więcej, co pozwoli
# na więcej danych wejściowych do analizy gsea.
# zamiast lfc można brać wynik samej statystyki (https://stephenturner.github.io/deseq-to-fgsea/)
# a ta pani sama wylicza sobie ranking: (https://biostatsquid.com/fgsea-tutorial-gsea/)
# co jest mega dziwne, ale twierdzi że tak trzeba powołując się na chińczyków.
# anyway:
# select columns from myTopHits:
myTopHits_tb <- as_tibble(myTopHits, rownames = "Symbol")
myTopHits_tb <- myTopHits_tb |> 
  select("Symbol", "logFC")
  
ranks <- deframe(myTopHits_tb)
head(ranks)
pathways.hallmark <- gmtPathways("~/Pobrane/c2.all.v2023.2.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(padj)

# plot the most significantly enriched pathway
# plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
#                rankings) + 
#   labs(title = head(GSEAres[order(padj), ], 1)$pathway)

my_pathway <- fgseaResTidy[[1,1]]

plotEnrichment(pathways.hallmark[[my_pathway]],
                ranks) +
  labs(title = my_pathway)
#
#
# Plot the normalized enrichment scores. Color the bar indicating 
# whether or not the pathway was significant:
fgseaResTidyTop <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) |> 
  head(50)

ggplot(fgseaResTidyTop, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
#
#
#
#


#datatable(diffGenes.df,
 #         extensions = c('KeyTable', "FixedHeader"),
  #        caption = 'Table 1: DEGs in cutaneous leishmaniasis',
   #       options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  #formatRound(columns=c(2:11), digits=2)



