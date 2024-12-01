# Functional-Analysis-of-Hub1
Downloading and Installation of Miniconda 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local
sh ./Miniconda3-latest-Linux-x86_64.sh

Installation of SRA-Toolkit
conda install -c bioconda sra-tools

Retrieval of SRR Files as zip Files
fastq-dump --gzip 3822520

Data Quality Control
fastqc SRR3822520.fastq.gz 

Data Trimming
cutadapt -u 15 -U 15 -o SRR3822520_1_trim.fastq.gz -p SRR3822520_2_trim.fastq.gz

Data Quantification
conda install kallisto
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz
kallisto index -i transcripts.idx rna_coding.fasta

Quantification
kallisto quant -i transcripts.idx -o SRR3820520 -b 100 SRR3820520_1_trim.fastq.gz SRR3820520_2_trim.fastq.gz

Differential Expression Analysis
library("tximport")
library("readr")
library("htmltools")
library("ggplot2")
library("DESeq2")
library("ggrepel")
library("plotly")
library("dplyr")
library("tidyverse")

#read_samples_file
samples <- read.table("~/raw_file/Xrn1- hub/file_xrn1.csv",sep = ",", header = TRUE)
#Import_SRR_files
rownames(samples) <- samples$SRR_Number
#Read_kallisto_abundances_tsv_files
files <- file.path("~/raw_file/Xrn1-hub/trimmed", samples$SRR_Number, "abundance.tsv")
txi <- tximport(files, type = "kallisto", txOut = TRUE)
#analyze_with_DEseq2
ddsTxi <- DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ Treatment)
dds <- DESeq(ddsTxi)
#write_results
res <- results(dds, contrast= c("Treatment","wildtype","hub1"))
res
write.csv(res,file = "~/rawdata/wildtype_hub1.csv")
#MA_plot
ggmaplot(res, main = expression("wild_type" %->% "Hub1"),
         fdr = 0.05, fc=1,size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res$name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         xlab = "mean of normalized expression",
         ylab = "Log2 fold change",
         ggtheme = ggplot2::theme_minimal())

#volcano_plot
png(filename="Vol_wild_type_vs_Hub1.png", width = 480, height = 720, bg="white")
par(mar=c(5,6,4,1)+.1)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlim=c(-5,4), cex = 1.0, cex.axis=2.0, cex.lab=2.0))
with(subset(res, pvalue < 0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex = 1.0))
with(subset(res, pvalue < 0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex = 1.0))
dev.off()

Gene Ontology (GO)
library("dplyr")
library("tidyverse")
library("ggplot2")

#Filtration_of_significant_DE 
filter_res <- (res[ which(res$padj < .05 & abs(res$log2FoldChange) >1), ])
filter_res <- as.data.frame(filter_res)
#Creating_column_of_up_&_down_regulated_genes
filter_res <- filter_res %>% mutate( de_status = case_when(
  log2FoldChange > 0 ~ "up",
  log2FoldChange < 0 ~ "down"
))
view(filter_res)
write.csv(filter_res,file = "~/raw_file/Xrn1-hub/trimmed/
#Read_GO_file
go <- read.csv("~/raw_file/Xrn1-hub/trimmed/gProfiler_wild_type_vs_Hub1.csv")
#Read_significant_columns_of_GO_file
go <- go %>% select("term_name", "adjusted_p_value", "term_size", "intersection_size")
go <- go %>% mutate(Gene_ratio=intersection_size/term_size)
#Mutate function is used to create a new variable from the given data file
ggplot(go, aes(x=Gene_ratio, y=term_name, size=adjusted_p_value, colour=intersection_size)) + geom_point()+
  theme(axis.text = element_text(size = 5), axis.title = element_text(size = 5), text = element_text(size = 5))

Gene Set Enrichment Analysis 
library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)
library(pathview)
library(ggplot2)

#reading in data from deseq2
df = read.csv("~/raw_file/Xrn1-hub/trimmed/wildtype_vs_Hub1.csv", header=TRUE)
#we want the log2 fold change
original_gene_list <- df$log2FoldChange
#name the vector
names(original_gene_list) <- df$X
#omit any NA values
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
library("org.Sc.sgd.db")
organism = "org.Sc.sgd.db"
gse <- gseGO(geneList=gene_list,
                ont ="ALL",
                keyType = "ENSEMBL",
                nPerm = 10000,
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                verbose = TRUE,
                OrgDb = organism,
                pAdjustMethod = "none")
#Use the “Gene Set” parameter for the index in the title, and as the value for geneSetId
gseaplot(gse_3, by = "all", title = gse$Description[1], geneSetID = 1)
gseaplot(gse, geneSetID = 'GO:0071944')
gseaplot(gse, geneSetID = 'GO:0006811')

Pathway Analysis
ids <-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
#remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
#Create a new data frame df2 which has only the genes which were successfully mapped using the bitr function above
df_h11 = df[df$X %in% dedup_ids$ENSEMBL,]
#Create a new column in df2 with the corresponding ENTREZ IDs
df_h11$Y = dedup_ids$ENTREZID
#Create a vector of the gene unuiverse
kegg_gene_list <- df_h11$log2FoldChange
#Name vector with ENTREZ ids
names(kegg_gene_list) <- df_h11$Y
#omit any NA values
kegg_gene_list<-na.omit(kegg_gene_list)
#sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "sce"
kk <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")

#Produce_the_native_KEGG_plot 
dotplot(kk, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
