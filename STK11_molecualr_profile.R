# Load the GenVisR package
# BiocManager::install("GenVisR")
# library(GenVisR)
# # Plot with the MAF file type specified (default) The mainRecurCutoff
# # parameter is described in the next section
# 
# 
# mutation_data=read.table(file="~/PDXC/240402_Drug_Resistance/aa.txt", header=TRUE, sep="\t")
# colnames(mutation_data)[c(1,8,12)]=c("sample","gene_name","trv_type")
# waterfall(mutation_data, fileType="MGI", mainXlabel=TRUE, mainLabelCol="amino.acid.change", mainLabelSize=2)
# mutation_data

###########################################################
setwd("/DATA07/home/gangtl22/PDXC/240402_Drug_Resistance/")

library(readxl)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)

################################################################################
#1 PC9AR 을 control 로 PC9, PC9GR, HCC827AR, H1975AR, YU1150YH 세포주들과 비교 
################################################################################
total_genes = c("NDUFB4","NDUFA7","NDUFA11","NDUFB11","NDUFA13","NDUFB9","NDUFB7","NDUFB10","NDUFB8","NDUFS8","NDUFB6","NDUFS2",
                "NDUFA2","NDUFA10","NDUFA6","NDUFS7","NDUFV1","NDUFA12","NDUFS3","NDUFV3","NDUFS6","NDUFA3","NDUFA9","NDUFA8",
                "NDUFS1","SDHB","SDHC","SDHD","UQCR11","UQCRQ","UQCRH","CYC1","UQCR10","UQCRB","COX6A2","COX7A1","COX6B2",
                "COX7A2L","COX4I1","COX6A1","COX10","COX6B1","COX5B","COX5A","COX7B","COX17","ATP5C1","ATP5J2","ATP5O","ATP5D",
                "ATP5E","ATP6V1E1","ATP5H","ATP5I","ATP5G2","ATP5L","ATP5B","ATP5G1","ATP5A1","ATP5G3") # COX8BP


tpm.df= read.table("./RNAseq/STK11_TPM_240508.tsv",sep="\t",header = TRUE,
                    stringsAsFactors = FALSE, row.names = 1, check.names=FALSE, quote = "");dim(tpm.df)
tpm.df = as.data.frame(tpm.df)
colnames(tpm.df)
tpm.df = tpm.df[,c("Gene Name","PC9AR","PC9","PC9GR-R","PC9GR3-R","HCC827AR","H1975AR","YU1150YH","MUORNAHCL4006AR")]
colnames(tpm.df)=c("Gene Name","PC9AR","PC9","PC9GR-R","PC9GR3-R","HCC827AR","H1975AR","YU1150YH","H4006AR")

tpm.df=as.data.frame(tpm.df)
sel.tpm=tpm.df[tpm.df$`Gene Name` %in% total_genes,]
dim(sel.tpm)
sel.tpm[duplicated(sel.tpm$`Gene Name`),]

sel.tpm2=sel.tpm[ rownames(sel.tpm) != "ENSG00000255339.2",]
sel.tpm2=sel.tpm2[rownames(sel.tpm2) != "ENSG00000255292.2",]
sel.tpm2=sel.tpm2[rownames(sel.tpm2) != "ENSG00000267059.2",]
sel.tpm2=sel.tpm2[rownames(sel.tpm2) != "ENSG00000167774.2",]
rownames(sel.tpm2)=sel.tpm2$`Gene Name`
sel.tpm3=subset(sel.tpm2,select=-`Gene Name`)
sel.mtx2=as.matrix(sel.tpm3)
sel.mtx=log2(sel.mtx2+1)
sel.mtx
sel.scal=scal(sel.mtx,"z","row")


# Mito-related gene set

human_comp1 = c("NDUFB4","NDUFA7","NDUFA11","NDUFB11","NDUFA13","NDUFB9","NDUFB7","NDUFB10","NDUFB8","NDUFS8","NDUFB6","NDUFS2",
                 "NDUFA2","NDUFA10","NDUFA6","NDUFS7","NDUFV1","NDUFA12","NDUFS3","NDUFV3","NDUFS6","NDUFA3","NDUFA9","NDUFA8",
                 "NDUFS1")
human_comp2 = c("SDHB","SDHC","SDHD")
human_comp3 = c("UQCR11","UQCRQ","UQCRH","CYC1","UQCR10","UQCRB")
human_comp4 = c("COX6A2","COX7A1","COX6B2","COX7A2L","COX4I1","COX6A1","COX10","COX6B1","COX5B","COX5A","COX7B","COX17")
human_comp5 = c("ATP5C1","ATP5J2","ATP5O","ATP5D", "ATP5E","ATP6V1E1","ATP5H","ATP5I","ATP5G2","ATP5L","ATP5B","ATP5G1","ATP5A1","ATP5G3")
# ATP5C1","ATP5J2","ATP5O","ATP5D","ATP5E","ATP6V1E1","ATP5H","ATP5K","ATP5G2","ATP5L","ATP5B","ATP5G1","ATP5A1","ATP5G3")

hcomp1=data.frame("Gene"=human_comp1,"Gene set"="Complex 1")
hcomp2=data.frame("Gene"=human_comp2,"Gene set"="Complex 2")
hcomp3=data.frame("Gene"=human_comp3,"Gene set"="Complex 3")
hcomp4=data.frame("Gene"=human_comp4,"Gene set"="Complex 4")
hcomp5=data.frame("Gene"=human_comp5,"Gene set"="Complex 5")
ann_col=rbind(hcomp1,hcomp2)
ann_col2=rbind(hcomp3,hcomp4)
ann_col3=rbind(ann_col,ann_col2)
ann_tot = rbind(ann_col3,hcomp5)
rownames(ann_tot)=ann_tot$Gene
mtx=sel.scal
mtx = mtx[!rowSums(is.na(mtx)),]

d <- dist(mtx, method = "euclidean")
hc1 <- hclust(d, method = "ward.D2")
d2 <- dist(t(mtx), method = "euclidean")
hc2 <- hclust(d2, method = "ward.D2")
blre <- colorRampPalette(c("#6252ff","white","#ff1a09"), space = "rgb")(100)
colnames(ann_tot)=c("Gene","Gene set")
ann_color = list("Gene set"=c("Complex 1"="#1b9e78",
                              "Complex 2"="#d95f01",
                              "Complex 3"="#7570aa",
                              "Complex 4"="#e7298a",
                              "Complex 5"="#67a61f"))

ann_tot$`Gene set` =factor(ann_tot$`Gene set`,levels = c("Complex 1","Complex 2","Complex 3","Complex 4","Complex 5"))
ann_tot=ann_tot[rownames(ann_tot) %in% rownames(mtx),]
ann_tot
ann_tot2=subset(ann_tot,select="Gene set")
ann_tot2
mtx=mtx[rownames(ann_tot2),]
tiff("./heatmap_24.tiff",units="px", width=950, height=1050, res=120,compression = "none")

#ann_tot2
aa=ComplexHeatmap::pheatmap(mtx,
                         name = "Expression",
                         color = blre,
                         #cluster_cols=hc2,
                         cluster_rows=FALSE,
                         annotation_row=ann_tot2,
                         annotation_colors = ann_color, 
                         row_split = ann_tot2$`Gene set`,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         annotation_legend = TRUE,
                         annotation_names_row = F,
                         annotation_names_col = F,
                         treeheight_row = 30,
                         treeheight_col = 15,
                         fontsize_row = 8.5,
                         fontsize_col = 10,
                         cellwidth = 35,
                         cellheight = 9,
                         legend = TRUE,
                         main=" ", 
                         border_color = NA,
                         row_title_gp = gpar(fontsize = 8))


print(aa)
dev.off()
rownames(mtx)
########################## ClusterProfile :: GO ##########################
#BiocManager::install("goseq")
#https://bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/Pathway-Analysis.html#4_Enrichment
if(!"rWikiPathways" %in% installed.packages()){if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")} 
  BiocManager::install("rWikiPathways", update = FALSE)
};library(rWikiPathways)
load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

lung.expr <- read.csv(system.file("extdata","data-lung-cancer.csv", package="rWikiPathways"),stringsAsFactors = FALSE)
nrow(lung.expr)
head(lung.expr)
resSig
gnid = lung.expr[lung.expr$GeneName %in% resSig$Gene.Name,c("Gene ID", "Gene Name")]
gnid = data.frame(gnid)
rownames(gnid) = gnid$Gene.Name
gnid



################################################################################
# 3 LKB1, ccnd1 의 mRNA level 비교 (PC9AR, PC9, PC9GR, HCC827AR, H1975AR, YU1150YH) - pathway analysis 는 아니고, mRNA level 만 비교 부탁드립니다.
################################################################################

sim.tpm= read.table("./Drug_Resistance_TPM_240416.tsv",sep="\t",header = TRUE,
                    stringsAsFactors = FALSE, row.names = 1, check.names=FALSE, quote = "");dim(sim.tpm)
sim3 = sim.tpm[c(sim.tpm$`Gene Name` == "STK11" |  sim.tpm$`Gene Name` == "CCND1"),]
sim3 = data.frame(sim3)
head(sim3)
sim3.log = log2(subset(sim3,select=-Gene.Name)+1)

ggplot(data=sim3.log, aes(x=Response, y=value,color=Response))+
  geom_boxplot(color=c("Non-pCR"="gray50","pCR" = "#c41e0c")) +
  ylim(0,hgt) +
  xlab("Clinical Response")+
  ylab(ttl)+
  facet_grid(~Group, labeller = labeller(Group=ckd_label)) +
  stat_compare_means(label = "p.format", paired = FALSE, label.x.npc = 0.3,label.y.npc = 0.95) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black",size = 10),
        axis.text.y = element_text(color="black",size = 10),
        legend.position =C(.5,.5)) +
  theme(legend.title = element_text(color = "black", size = 15, face = "bold"))+
  #범레 텍스트 색,크기,진하게 설정
  theme(legend.text = element_text(color = "black", size = 15, face = "bold")) +
  theme(strip.text.x = element_text(size = 13))








tpm_mt = read.table("~/ACTS29/RNAseq/ACTS29_PrePost_tpm_210315.txt",sep="\t",header = TRUE,
                    stringsAsFactors = FALSE, row.names = 1, check.names=FALSE, quote = "");dim(pp_mt)
# half_mt = read.table("~/ACTS29/RNAseq/ACTS29_D15_rawcount_220207.txt",sep="\t",header = TRUE,
#                      stringsAsFactors = FALSE, row.names = 1, check.names=FALSE, quote = "");dim(half_mt)
# mon_mt = read.table("~/ACTS29/RNAseq/ACTS29_D29_rawcount_220207.txt",sep="\t",header = TRUE,
#                     stringsAsFactors = FALSE, row.names = 1, check.names=FALSE, quote = "");dim(mon_mt)
tpm_pre=tpm_mt[,c(rs1,nr1)]
writeDF(tpm_pre,"~/ACTS29/RNAseq/ACTS_TPM_PRE.txt")
