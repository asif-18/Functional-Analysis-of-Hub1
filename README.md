# Functional-Analysis-of-Hub1
3.1.1   Downloading and Installation of Miniconda 
The following command lines were used in the terminal interface of the mac operating system to download and install Miniconda. Miniconda is a well-designed, free, smaller version of Anaconda that installs Conda (a package manager feature). Miniconda is used to simplify package management and deployment.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local
sh ./Miniconda3-latest-Linux-x86_64.sh
3.1.2   Installation of SRA-Toolkit
        SRA-Tools was installed by running the following code:
conda install -c bioconda sra-tools
3.1.3   Retrieval of SRR Files as zip Files
         For downloading SRR files which is enlisted in the metadata table as Zip file the following command-lines was used:
fastq-dump --gzip 3822520
Here, --gzip = “--gzip” flag is used to download the sequence file as a compressed file.
SRR3822520 = Sequence Read Archive (SRA) data, available through multiple cloud providers and NCBI servers through SRA Results are called runs (SRR). Runs comprise the data gathered for a sample or sample bundle and refer to a defining experiment.
3.2   Data Quality Control

          The quality of the gene sequence data in the downloaded fastq files has been assessed by a bioinformatics tool FastQC v0.11.9 and using the following command lines.
fastqc SRR3822520.fastq.gz (for single SRR file)
 fastqc *fastq.gz (for all SRR files)
The output from FastQC, after analyzing a FASTQ file of sequence reads, is an html file. Several parameters have been assessed in the FastQC report to check the quality of the sequence data such as basic statistics, per base sequence quality, per sequence quality scores, per base sequence content, per base GC content, per sequence GC content, per base N content, sequence length distribution, sequence distribution levels, overrepresented sequences, and k-mer content (Brown et al., 2017).
3.3   Data Trimming

          The sequence data has been trimmed either at 15-20 base pair upstream or downstream region using a bioinformatics tool cutadapt 3.4 and using the following command.

          cutadapt -u 15 -U 15 -o SRR3822520_1_trim.fastq.gz -p SRR3822520_2_trim.fastq.gz

Here,
	cutadapt = Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
	-u = Remove bases from each read (only if single). If length is positive, remove bases from the beginning. To remove bases from the end “u” is negative. 
	-o = Write trimmed reads to file. FASTQ or FASTA format is chosen depending on input. 
	SRR3822520 trim.fastq = Output file name
	SRR3822520.fastq.gz = Read the sequence file
3.4   Data Quantification

         The sequence data has been quantified by indexing cDNAs of Saccharomyces (Saccharomyces_cerevisiae.R64-1-1.cdna.all) as a parameter in 100 bootstrips using a bioinformatics tool kallisto 0.46.2 (Pimentel et al., 2017). Kallisto was installed by running the following code:
conda install kallisto
The reference fasta file of cDNAs of Saccharomyces was downloaded by following the command line.
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz

Indexing was performed to create an index file using the above reference fasta file. The following command was used for indexing.
                                           kallisto index -i transcripts.idx rna_coding.fasta
Quantification of data was performed using SRR trimmed files with the following command.

kallisto quant -i transcripts.idx -o SRR3820520 -b 100 SRR3820520_1_trim.fastq.gz SRR3820520_2_trim.fastq.gz
3.5   Differential Expression Analysis

            Differential expression of genes has been analyzed by the DESeq2 of Bioconductor library of R version 4.1.1. The data has been imported in R in a csv file format. MA plot and volcano have been generated for differentially expressed genes. Log2 fold change, p-value, and adjusted p-value have been calculated according to Love et al., 2014.
          At first, Bioconductor has been installed in Rstudio and the following commands were used to generate MA plot and Volcano plot,

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
3.6   Gene Ontology (GO)

             Gene ontology enrichment analysis for biological processes, molecular function, and cellular components has been performed using g:Profiler (http://biit.cs.ut.ee/gprofiler/) is a public web server for characterizing and manipulating gene lists resulting from mining high-throughput genomic data. 
         The ClusterProfiler package has been used to visualize GO enrichment results according to (Yu et al., 2012).
In Rstudio Bioconductor installed integrated with the following packages:
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

3.7   Gene Set Enrichment Analysis 

         Gene set enrichment analysis (GSEA) has been executed by an R package clusterProfiler library version:4 according to (R. Z. Wu et al., 2007). 
In Rstudio Bioconductor installed integrated with the following packages:

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
3.8   Pathway Analysis
           The metabolic or functional pathway has been analyzed by KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway according to Zhang & Wiemann, 2009. 
#Convert gene IDs for gseKEGG function
#We will lose some genes here because not all IDs will be converted
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
