#Script for converting loom files to single cell experiments for all spliced and all unspliced matrices respectively
#LASeeker
#20200418
#LOAD LIBRARIES
library(lattice)
library(Matrix)
library(KernSmooth)
library(MASS)
library(nlme)
library(cluster)
library(survival)
library(devtools)


#load Seurat
#print(find.package("Seurat"))
#devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
#install.packages('Seurat', "/exports/cmvm/eddie/scs/groups/Williamsdata/Rlibraries", repos = "https://cran.ma.imperial.ac.uk/")
library(Seurat)
library(sctransform)
library(hdf5r)
library(pcaMethods)

#load velocyto packages
library(loomR)
#for macOS to install velocyto R do devtools::install_github('tractatus/velocyto.R') otherwise there might be compiler issues
library(velocyto.R)

#load other packages
#library(scater)
library(ggplot2)

library(cowplot)
library(dplyr)
library(future)



e<- 40000*1024^2
options(future.globals.maxSize= e)



######################Set paths to files###################################################

PathToAllLoomTxt<-"/media/annawilliams/Data/Sunniva/Velocyto/Loom.txt"
PathToMetadata<- "/media/annawilliams/Data/Sunniva/Velocyto/metadata.csv"
PathToLoomFiles<-"/media/annawilliams/Data/Sunniva/Velocyto/"


######################Set paths to output folders##########################################

OutputSCE_unspliced<-"/exports/eddie/scratch/s1359339/Mouse_Extraction/datasets/CombinedSingleCellObjects_unspliced"
OutputSCE_spliced<-"/exports/eddie/scratch/s1359339/Mouse_Extraction/datasets/CombinedSingleCellObjects_spliced"




folderNames<-readLines(PathToAllLoomTxt) 

## read in raw metadata
#metadataRaw<-read.csv(PathToMetadata)

########################Set variables fro creating Seurat object###########################
minCells<-0 #you can set QC thresholds here if you like (for example enter a 3 here and remove all genes that were expressed in less cells)
minFeatures<-0 # if you want to filter some more, you can remover here genes that are not frequently expressed (for example all genes below 200 copies detected across all cells)

#initialise iteration counters and empty lists:


endIt<-length(folderNames)
#i<-1

for(i in 1: endIt){
  fileNames<- list.files(path = paste(PathToLoomFiles, folderNames[i], sep="/"), 
                         pattern = "*.loom", all.files = FALSE,
                         full.names = FALSE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE) 
  
  #extract the unique identifier "processNumber"    from the folder names
  processNumber<-sapply(strsplit(folderNames[i],"_"), `[`, 2)
  processNumber<-sapply(strsplit(processNumber,"_"), `[`, 1) 
  
  #generate a project  and sample name that will be used for saving data and naming Seurat objects
  sampleName<-paste("10X", folderNames[i], sep="_") 
  
  
  
  #read in .loom file
  ldat2 <- read.loom.matrices(paste(PathToLoomFiles, "/" , folderNames[i], "/" , fileNames, sep="")) 
  
  #extract a matrix with counts for spliced mRNAs
  emat2 <- ldat2$spliced 
  #extract a matrix with counts for unspliced mRNAs
  nmat2 <- ldat2$unspliced 
  
  
  
  
  #The following command makes the rownames of the matrix (= feature names= gene names) unique. 
  #This is what the bioinformaticins in Stockholm seem to do. 
  #I did not like it, because it simply adds an ".1" to a feature that #has been detected above so you can potentially end up with count #data for the same gene in two different rows. 
  #I had a look at which genes are affected and it is usually the ones that are not highly expressed and are unlikely to be important for the downstream analysis.
  #However, I think it is a good idea to save data about which genes are affected and how much they are expressed.
  #So we will save a .csv file for each sample containing this information. 
  
  rownames(emat2)<-make.unique(rownames(emat2))
  rownames(nmat2)<-make.unique(rownames(nmat2))
  
  
  #prepare to link count matrix to metadata:
  
  df<-data.frame(Barcode= colnames(nmat2) , ProcessNumber= as.integer(paste(processNumber)))
  
  #Link barcodes to metadata
  metadata<-df
  
  #rownames must be the barcodes to allow merging
  
  rownames(metadata) <- metadata$Barcode
  
  
  y_unspliced<- CreateSeuratObject(counts =nmat2, min.cells= minCells, min.features = minFeatures,  project = sampleName, assay = "RNA", names.field = 1, names.delim = ".", meta.data = metadata ) #converts the combined matrix to a Seurat object
  y_unspliced[["percent.mt"]] <- PercentageFeatureSet(y_unspliced, pattern = "mt-")
  
  
  y_spliced<- CreateSeuratObject(counts =emat2, min.cells= minCells, min.features = minFeatures,  project = sampleName, assay = "RNA", names.field = 1, names.delim = ".", meta.data = metadata) #converts the combined matrix to a Seurat object
  y_spliced[["percent.mt"]] <- PercentageFeatureSet(y_spliced, pattern = "mt-")
  
  if(i==1){
    x_unspliced<-y_unspliced
    x_spliced<-y_spliced
  }else{
    x_unspliced<- merge(x = x_unspliced, y =  y_unspliced, project="HCA_raw_unspliced") 
    x_spliced<- merge(x = x_spliced, y =  y_spliced, project="HCA_raw_spliced") 
  }
  
}



NumbOfSamples<-i




singleCell_unspliced<-as.SingleCellExperiment(x_unspliced)
singleCell_spliced<-as.SingleCellExperiment(x_spliced)

saveRDS(singleCell_unspliced,  paste(OutputSCE_unspliced, "/", Sys.Date(), "_", "CombinedSCE_", NumbOfSamples, "unspliced_samples.RDS", sep=""))  
saveRDS(singleCell_spliced,  paste(OutputSCE_spliced, "/", Sys.Date(), "_", "CombinedSCE_", NumbOfSamples, "spliced_samples.RDS", sep=""))  


#End of Script