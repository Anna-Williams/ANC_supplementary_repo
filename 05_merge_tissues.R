---
  title: "2021_CombinedAnAllCells"
author: "Luise A. Seeker"
date: "19/01/2021"
output: html_document
---
  
  
  
  ```{r}
setwd("/exports/eddie/scratch/s1359339/HD/04_scran_normalised/")
library(Seurat)
library(ggplot2)
library(ggsci)
library("scales")
library(SingleCellExperiment)
library(BiocSingular)
library(scDblFinder)
library(scran)
library(scater)
library(RCurl)
library(AnnotationHub)
library(here)
```


Pick colour paletts

```{r}
mypal <- pal_npg("nrc", alpha = 0.7)(10)
mypal2 <-pal_tron("legacy", alpha = 0.7)(7)
mypal3 <- pal_lancet("lanonc", alpha = 0.7)(9)
mypal4 <- pal_simpsons(palette = c("springfield"), alpha = 0.7)(16)
mypal5 <- pal_rickandmorty(palette = c("schwifty"), alpha = 0.7)(6)
mypal6 <- pal_futurama(palette = c("planetexpress"), alpha = 0.7)(5)
mypal7 <- pal_startrek(palette = c("uniform"), alpha = 0.7)(5)
mycoloursP<- c(mypal, mypal2, mypal3, mypal4, mypal5, mypal6, mypal7)
show_col(mycoloursP, labels =F)
```


Set working directory and read in file

```{r}
norm_FrCx<-readRDS(here("FrCx/FrCx_sce_fil_norm.RDS"))
norm_CB<-readRDS(here("CB/CB_sce_fil_norm.RDS"))
norm_HC<-readRDS(here("HC/HC_sce_fil_norm.RDS"))
norm_CN<-readRDS(here("CN/CN_sce_fil_norm.RDS"))
combined_sce <- cbind(norm_FrCx, norm_CB, norm_HC, norm_CN)
```


Bioconductor offers packages that help identifying doublets in a dataset. 
We use scDblFinder which simuated doublets by randomly mixing the transcription
profile of cells in a dataset and then compares real samples to those simulated
ones. I like this method, because it allocated doublet scores to cells not 
clusters and therefore those scored remain if the clustering is changed. 

For more information look here:
  https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html
https://github.com/plger/scDblFinder



I tested the scDblFinder algorithm with and without running a dimensional 
reduction first and the results look very similar to me. Here the SCE dataset
is processed to add the PCA dimensional reduction which will be used if present.
Because this is so close to visualising the dataset, I do this here, too. 

Ideally, scDblFinder should be performed after empty droplets are removed and
before further filtering is performed. But on the other hand the authors say 
that it works with normalised data and with dimensional reduced data which 
are steps that are usually performed after QC. 

It does not make sense that doublets form between samples that are run on 
different chip. Therefore we specify a sample identifier and the algorithm 
first runs within sample. 

```{r}
set.seed(100)
combined_sce <- scDblFinder(combined_sce, samples="ProcessNumber")
table(combined_sce$scDblFinder.class)
```


```{r}
saveRDS(combined_sce,here("combined_SCE/20210309combined_SCE_normandraw.RDS"))
dir.create("here(doublet_score_file"))
write.csv(as.data.frame(colData(combined_sce)),here("doublet_score_file/20210309__doubletscore.csv"))
# the object above is still a singleCellExperiment which contains two 
# data slots: counts and logcounts. 
# counts is the  raw data, logcounts the normalised data. For now I am only 
# interested in the normalised data
```

Then combine all three datasets to a single one.

check that raw data is saved in seurat@assays$RNA@counts an log normalised
data is saved in seurat@assays$RNA@data. Raw data can be recognised as 
integers, whereas log normalised counts show decimal numbers.

```{r}
FrCx_seurat <- as.Seurat(norm_FrCx)
CB_seurat <- as.Seurat(norm_CB)
CN_seurat <- as.Seurat(norm_CN)
HC_seurat <- as.Seurat(norm_HC)

FrCx_seurat@meta.data$Tissue <- "FrCx"
CB_seurat@meta.data$Tissue <- "CB"
CN_seurat@meta.data$Tissue <- "CN"
HC_seurat@meta.data$Tissue <- "HC"


seur_comb <- merge(FrCx_seurat, y = c(CB_seurat, HC_seurat, CN_seurat), 
                   add.cell.ids = c("FrCx", "CB", "HC", "CN"), 
                   project = "HD_all_celltypes")
# check if cell order of seurat and sce objects are the same
summary(seur_comb$Barcode == combined_sce$Barcode)
##all true!


seur_comb$scDblFinder.weighted <- combined_sce$scDblFinder.weighted
seur_comb$scDblFinder.ratio <- combined_sce$scDblFinder.ratio
seur_comb$scDblFinder.class <- combined_sce$scDblFinder.class
seur_comb$scDblFinder.score <- combined_sce$scDblFinder.score
remove(combined_sce)
dir.create(here("combined_Seurat"))

saveRDS(seur_comb,here("20210309_seur_combined_rawandnorm.RDS"))
```

```{r}
sessionInfo()
```



