# script dividing SCE into tissues
# Luise A. Seeker
# 20210118


###
# 



# load libraries

library(scater)
library(scran)
library(SCE_comb)
library(devtools)
library(dplyr)

###### PART 4
###### SPLIT DATASET IN TISSUES AND SCRAN NORMALISATION

SCE_comb <- readRDS(file.path(getwd(), 
                              "splice_control_out",
                              "datasets",
                              "03_combined_matrices",
                              "combined_QC_filtered.RDS"))


# I noticed that there are 5 levels in the column "Tissue" instead of three, 
# because CB and BA4 may have an empty space behind their last letter or not,
# I am correcting this here
SCE_comb$donor_id <- ifelse(SCE_comb$donor_id == "NBB95-310 193A ", "NBB95-310",
                          ifelse(SCE_comb$donor_id == "EBBSD025/13 ", "EBBSD025/13",
                                 ifelse(SCE_comb$donor_id == "EBBSD031/14 ", "EBBSD031/14",
                                        ifelse(SCE_comb$donor_id == "NBB16-056 ", "NBB16-056",
                                               SCE_comb$donor_id))))
SCE_comb$tissue_id <- ifelse(SCE_comb$tissue_id == "IFG", "FrCx",
                           SCE_comb$tissue_id)

SCE_comb@colData<- within(SCE_comb@colData, donor_status[donor_status == 'HD' & donor_id == 'L-E14-13'] <- 'Ctrl')





normaliseSCE <- function(sce_dataset, split_id_column){
  for(i in 1: length(levels(as.factor(sce_dataset[[split_id_column]])))){
    tissue <- levels(as.factor(sce_dataset[[split_id_column]]))[i]
    split_boul <- sce_dataset[[split_id_column]] == tissue
    split_sce <- sce_dataset[, split_boul]
    
    # scran normalisation
    ## calculation of deconvolution factors
    set.seed(100)
    clust_sce <- quickCluster(split_sce) 
    #deconv_sf <- calculateSumFactors(split_sce, cluster=clust_sce)
    split_sce <- computeSumFactors(split_sce, 
                                   cluster=clust_sce, 
                                   min.mean=0.1)
    
    split_sce <- logNormCounts(split_sce)
    
    dir.create(file.path(getwd(), 
                         "splice_control_out",
                         "datasets",
                         "04_scran_normalised",
                         paste(tissue)))
    
    
    saveRDS(split_sce, file.path(getwd(), 
                                 "splice_control_out",
                                 "datasets",
                                 "04_scran_normalised",
                                 paste(tissue),
                                 paste(tissue,
                                       "sce_fil_norm.RDS",
                                       sep = "_")))
  }
  
}


normaliseSCE(SCE_comb, "donor_id")



print("done")