if(!require("maftools")){BiocManager::install("maftools")};library(maftools)
if(!require("dplyr")){install.packges("dplyr")};library(dplyr)

load.packg=c("maftools","dplyr")

###reading mutiple .maf files as a large list
setwd("/DATA07/home/gangtl22/CJ/Inhouse/WES/Analysis/MAF_files/")
maf.filenames = list.files(full.names=TRUE,pattern = "funcotator.maf")


clinc.info=read.table("/DATA07/home/gangtl22/CJ/Inhouse/RNAseq/Analysis/CibX/CJRNA92_ecotyper_annote.tsv",sep="\t",
                      header = TRUE,stringsAsFactors = TRUE)
clinc.info$Tumor_Sample_Barcode=paste0("CJWES",clinc.info$ID)
rownames(clinc.info)=clinc.info$Tumor_Sample_Barcode


if(".//MAFTOOLS" %in% list.dirs("./")){
  print("Bye")
}else{
  dir.create("./MAFTOOLS")
}
clinc.info$Tumor_Sample_Barcode

maf.filenames %in% clinc.info$Tumor_Sample_Barcode

for ( i in clinc.info$Tumor_Sample_Barcode){
  i="CJWESMB089"
  maf.fl = paste0("./",i,".funcotator.maf")
  if ( maf.fl  %in% maf.filenames ){
    d=read.maf(maf.fl)
    dim(d@data)
    d@data$Tumor_Sample_Barcode = paste0(i,"T")
    d@data$Matched_Norm_Sample_Barcode = paste0(i,"N")
    d@clinical.data = data.table(clinc.info[i,])
    writeDF(d@data,paste0("./MAFTOOLS/",i,"_nonSyn.maf"))
    #write.mafSummary(d,basename = paste0("./MAFTOOLS/",i))
  }else{
    print(paste0(maf.fl," does not exist"))
  }
}


mer.mf <- merge_mafs(lapply(Sys.glob("./MAFTOOLS/*_nonSyn.maf"), read.maf))
write.mafSummary(mer.mf,basename = "./MAFTOOLS/CJ70_nonSyn.maf")

################################################################################
######## Sample Number Discrepancy Problem 
######## After Merging, remove Tumor_Sample_Barcode == "__UNKNOWN__"
######## Maybe those "__UNKNOWN__" because 
######## "removeSilent" in read.maf option
######## logical. Whether to discard silent (with no functional impact) mutations 
######## ("Silent","Intron","RNA","3'UTR"). Default is TRUE.
######## so that they were filtered by maftools ########
################################################################################

mer.mf=read.maf("./MAFTOOLS/CJ70_nonSyn.maf",)
mer.mf@clinical.data = data.table(clinc.info)

######## Changing Colors for Variants Classification #######
vc_cols = RColorBrewer::brewer.pal(n = 7, name = 'Paired')
table(mer.mf@data$Tumor_Sample_Barcode)

names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'In_Frame_Del'
)
print(vc_cols)
ck_genes[1:10]
oncoplot(maf = mer.mf, colors = vc_cols, genes=ck_genes[1:10])
######## Bar plots #######
######## leftBarData, rightBarData and topBarData arguments can be used to #######
######## display additional values as barplots. Below example demonstrates #######
####### adding gene expression values and mutsig q-values as left and right #######
####### side bars respectively. #######

#Selected genes
ck_genes = c("EGFR","TP53","KRAS","BRAF","ALK")

#Variant allele frequcnies (Right bar plot)
ck_genes_vaf = subsetMaf(maf = mer.mf, genes = ck_genes, 
                         fields = "i_TumorVAF_WU", 
                         mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]



colnames(aml_genes_vaf)[2] = "VAF"
head(aml_genes_vaf)
#>    Hugo_Symbol      VAF
#>         <char>    <num>
#> 1:       ASXL1 37.11250
#> 2:       CEBPA 22.00235
#> 3:      DNMT3A 43.51556
#> 4:      DNMT3B 37.14000
#> 5:        EZH2 68.88500
#> 6:        FLT3 34.60294

#MutSig results (Right bar plot)
laml.mutsig = system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.mutsig = data.table::fread(input = laml.mutsig)[,.(gene, q)]
laml.mutsig[,q := -log10(q)] #transoform to log10
head(laml.mutsig)
#>      gene        q
#>    <char>    <num>
#> 1:   FLT3 12.64176
#> 2: DNMT3A 12.64176
#> 3:   NPM1 12.64176
#> 4:   IDH2 12.64176
#> 5:   IDH1 12.64176
#> 6:   TET2 12.64176
oncoplot(
  maf = laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  rightBarLims = c(0, 20)
)
  
  
  ########  Including annotations ########
  getClinicalData(x = mer.mf)
  typecolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
  names(typecolors) = c("E19del","L858R","WT","SQCC","E20insMutant","E19del&T790M","G719AMutant","SCLC")
  typecolors = list(Lung_Cancer_Type = typecolors)
  
#1140 500
tiff("/DATA07/home/gangtl22/GENESIS/Onco_1.tiff",
       units="px", width=1140, height= 500, res=120,compression = "none")
o1=oncoplot(maf = mer.mf, genes = ck_genes, clinicalFeatures = 'Lung_Cancer_Type',
            sortByAnnotation = TRUE,annotationColor = typecolors)
oncoplot(maf = mer.mf, genes = ck_genes, clinicalFeatures = 'Lung_Cancer_Type',
         sortByAnnotation = TRUE,annotationColor = typecolors)
print(o1)
dev.off()
