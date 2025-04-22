#Author: Bilge Elitok, 22.04.2025

#Dependencies:
library(clusterProfiler)
library(dplyr)
library(GSEABase)
library(DOSE)
library(ggplot)

#Data will be inputted.
PROJECT_FOLDER <- "C:/Users/../OneDrive/..../" #to be completed

#Hallmark gene  set collection from MSigDB:
hallmark.gmt <- file.path(PROJECT_FOLDER, "h.all.v7.2.symbols.gmt")
#term2gene conversion with read.gmt() function of GSEAbase:
pathways <- read.gmt(hallmark.gmt)

#check the unique pathways:
cat("Number of unique pathways found: ", length(unique(pathways$term)))


down <- read.table(file_path(PROJECT_FOLDER, "expression_matrix_for_down"),
                   sep = '\t',
                   header = TRUE,
                   stringsAsFactors = FALSE)

up <- read.table(file_path(PROJECT_FOLDER, "expression_matrix_for_up"),
                   sep = '\t',
                   header = TRUE,
                   stringsAsFactors = FALSE)

#GSEA is directional, we need ranked gene list (list of character vectors).
down_genes <- down$"logFC"
names(down_genes) <- down$"genes"

up_genes <- up$"logFC"
names(up_genes) <- up$"genes"

#Addition from BE: append rows and use both UP+DOWN for GSEA.
all_genes <- rbind.data.frame(up_genes, down_genes[2:nrow(down_genes),])
degs <- all_genes$logFC
names(degs) <- all_genes$Genes

degs <- sort(degs, decreasing = TRUE)


#GSEA main:
gsea_down <- clusterProfiler::GSEA(gene_list = down_genes,
                                   TERM2GENE = pathways,
                                   pvalueCutoff = 0.05)

#GSEA main:
gsea_up <- clusterProfiler::GSEA(gene_list = up_genes,
                                   TERM2GENE = pathways,
                                   pvalueCutoff = 0.05)

gsea_all  <- clusterProfiler::GSEA(gene_list = degs,
                                   TERM2GENE = pathways,
                                   pvalueCutoff = 0.05)


#Make figures, save figures.

p1 <- dotplot(gsea_down, showCategory=35)
ggsave(file.path(PROJECT_FOLDER, "downreg_BE_22042025.png"),
                 plot = p1, width = 25, height = 35)

p2 <- dotplot(gsea_up, showCategory=35)
ggsave(file.path(PROJECT_FOLDER, "upreg_BE_22042025.png"),
        plot = p2, width = 25, height = 35)

p2 <- dotplot(gsea_all, showCategory=35)
ggsave(file.path(PROJECT_FOLDER, "alldegs_BE_22042025.png"),
       plot = p2, width = 25, height = 35)


#Check project file.
