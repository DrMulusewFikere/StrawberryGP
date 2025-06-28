#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Population structure analysis using a raw Geno data #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # MFNote: Here the strategy is after population cluster analysis, then do sub-population based imputation #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. Define path and library #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(ggdendro)
library(dendextend)
library(RColorBrewer)
library(ggplot2)
library(cluster)
library(factoextra)

# setwd("/Users/uqmfiker/Dropbox/Post_201203/Data_post220524")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Read  data #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.1. Read phenotype data #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#CH    refinedPheno_2064 <- read.csv("adata_S01.csv", h=T, stringsAsFactors = F) 
    refinedPheno_2064 <- read.csv("../refineData/adata_S01.csv", h=T, stringsAsFactors = F) 
    length(unique(refinedPheno_2064$ID))
        
        #>>>>>>>>>>>>>>>>>>>>
        # 3.1.1. Congreunce of entries ###
        #>>>>>>>>>>>>>>>>>>>>
        phdata_congrev <- refinedPheno_2064[!is.na(refinedPheno_2064$SSC),]
        uTID <- unique(phdata_congrev[,c("T","ID")])
        uTIxG <- table(uTID$ID,uTID$T)
        G_Tcong <- t(uTIxG)%*%uTIxG
        G_Tcong
        
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 3.2. Read RAW genotype data #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
#CH    rawGeno_wzNA <- read.csv("sdata_withMissing3407x12951.csv", h=T, stringsAsFactors = F, na.strings = "NA") 
    rawGeno_wzNA <- read.csv("../rawData/sdata_withMissing3407x12951.csv", h=T, stringsAsFactors = F, na.strings = "NA") 
    dim(rawGeno_wzNA)
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Re-structure genotype data and keep the 2064 IDs #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    uID_2064 <- unique(refinedPheno_2064$ID)
    rawGeno_wzNA_2064 <- rawGeno_wzNA[match(uID_2064, rawGeno_wzNA$ID),]
    rawGeno_wzNA_2064_tmp <- rawGeno_wzNA_2064[,2:ncol(rawGeno_wzNA_2064),]
    rownames(rawGeno_wzNA_2064_tmp) <- rawGeno_wzNA_2064$ID
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 5. Population cluster analysis #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.1. Define number of optimal K
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    rawGeno_wzNA_012 <- as.matrix(rawGeno_wzNA_2064_tmp)
    rawGeno_wzNA_012_mat <-  rawGeno_wzNA_012/2
    rawGeno_wzNA_dist = as.matrix(dist(rawGeno_wzNA_012_mat))
    rawGeno_wzNA_Mds = cmdscale(rawGeno_wzNA_dist, eig=TRUE, x.ret=TRUE)
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.2. Visualize and save high resolution image
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
      # Elbow method
      optimal_nK_elbow <- fviz_nbclust(rawGeno_wzNA_dist, kmeans, method = "wss") +
        geom_vline(xintercept = 2, linetype = 2)+
        labs(subtitle = "Elbow method")
      
#CH      jpeg('../Analysis_post220605/1. Stucture/optimal_nK_elbow.jpeg', width = 566,height = 357, quality = 1100)
      jpeg('optimal_nK_elbow.jpeg', width = 566,height = 357, quality = 1100)
      print(optimal_nK_elbow)
      dev.off()
      
      # Silhouette method
      optimal_nK_silh <- fviz_nbclust(rawGeno_wzNA_dist, kmeans, method = "silhouette")+
        labs(subtitle = "Silhouette method")
      
#CH      jpeg('../Analysis_post220605/1. Stucture/optimal_nK_silh.jpeg', width = 566,height = 357, quality = 1100)
      jpeg('optimal_nK_silh.jpeg', width = 566,height = 357, quality = 1100)
      print(optimal_nK_silh)
      dev.off()

      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 5.3. Visualize individuals grouped in PCA dimensional reduction plot
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      rawGeno_wzNA_mds.var.per = round(rawGeno_wzNA_Mds$eig/sum(rawGeno_wzNA_Mds$eig)*100, 1)
      rawGeno_wzNA_mds.values = rawGeno_wzNA_Mds$points
      rawGeno_wzNA_mds.data = data.frame(Sample=rownames(rawGeno_wzNA_mds.values), X=rawGeno_wzNA_mds.values[,1], Y=rawGeno_wzNA_mds.values[,2])
      rawGeno_wzNA_scree_plot_dat = data.frame(PCoord = seq(1, length(rawGeno_wzNA_mds.var.per)), per.var.explained = rawGeno_wzNA_mds.var.per)
    
      rawGeno_wzNA_mds.data$Sample<- NULL
      colnames(rawGeno_wzNA_mds.data) [1] <- "PC1"
      colnames(rawGeno_wzNA_mds.data)[2] <- "PC2"
      rawGeno_wzNA_mds.data$Sample<- NULL
      rawGeno_wzNA_k2 <- kmeans(rawGeno_wzNA_mds.data, centers = 2, nstart = 25)
   
      xlab_name = paste("Principal Coord1 (", rawGeno_wzNA_mds.var.per[1], "%)", sep = "")
      ylab_name = paste("Principal Coord2 (", rawGeno_wzNA_mds.var.per[2], "%)", sep = "")
      
      PCA_dimenssionalRed <- fviz_cluster(rawGeno_wzNA_k2, geom = "point", data = rawGeno_wzNA_mds.data, ellipse.type = "convex",palette = "jco") + ggtitle("K = 2") + xlab(xlab_name) + ylab(ylab_name)+ theme_bw() 
      
#CH      jpeg('../Analysis_post220605/1. Stucture/PCA_dimenssionalRed.jpeg', width = 566,height = 357, quality = 1100)
      print(PCA_dimenssionalRed)
      dev.off()
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 6. Subset Population 1 and 2 from the "rawGeno_wzNA_2064" #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.1. Visualize individuals in PCA dimensional reduction plot
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
    ID_clust <- as.data.frame(rawGeno_wzNA_k2$cluster)
    colnames(ID_clust)[1] <- "Cluster"
    ID_clust$ID <- rownames(ID_clust)
    rownames(ID_clust) <- NULL
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.2. sub-set phenotypic data by population cluster index #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    refinedPheno_2064_temp <- merge (refinedPheno_2064, ID_clust, by= "ID")
    refinedPheno_pop1 <- refinedPheno_2064_temp[refinedPheno_2064_temp$Cluster == "1",]
    refinedPheno_pop2 <- refinedPheno_2064_temp[refinedPheno_2064_temp$Cluster == "2",]
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.3. sub-set genotype data by population cluster index #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    rawGeno_wzNA_2064_wizCluster <- merge(rawGeno_wzNA_2064, ID_clust, by= "ID")
    rawGeno_wzNA_pop1 <- rawGeno_wzNA_2064_wizCluster[rawGeno_wzNA_2064_wizCluster$Cluster == "1",]
    rawGeno_wzNA_pop2 <- rawGeno_wzNA_2064_wizCluster[rawGeno_wzNA_2064_wizCluster$Cluster == "2",]
    
    ## Remove the Cluster column from the genotype data
    rawGeno_wzNA_pop1$Cluster <- NULL
    rawGeno_wzNA_pop2$Cluster <- NULL

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 7. Save workspace and read-in this workspace for IMPUTATION analysis (refer section 2. Imputation) #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#CH save.image("../Analysis_post220605/1. Stucture/strawberryStructure.RData")
save.image("strawberryStructure.RData")
