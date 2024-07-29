#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Population structure based on IMPUTATION #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# MFNote: Here the strategy is after we defined population cluster analysis, followed by sub-population based imputation #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. Define path and library #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
options(rgl.useNULL = TRUE)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("impute")

library(rgl)
library(snpReady)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)

# setwd("/Users/uqmfiker/Dropbox/Post_201203/Data_post220524")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Load workspace from structure section and read mapFile #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#CH load("../Analysis_post220605/1. Stucture/strawberryStructure.RData")
load("strawberryStructure.RData")

mapData <- read.csv("../rawData/strawberry_mapFile.csv", h=T, stringsAsFactors = F)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Re-structure population 1 and 2 genotypype data #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 4.1. rawGenotypes of population 1 #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    rawGeno_wzNA_pop1_tmp <- rawGeno_wzNA_pop1[,2:ncol(rawGeno_wzNA_pop1)]
    rownames(rawGeno_wzNA_pop1_tmp) <- rawGeno_wzNA_pop1[,1]

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 4.2. rawGenotypes of population 2 #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    rawGeno_wzNA_pop2_tmp <- rawGeno_wzNA_pop2[,2:ncol(rawGeno_wzNA_pop2)]
    rownames(rawGeno_wzNA_pop2_tmp) <- rawGeno_wzNA_pop2[,1]

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 5. Visualize missing SNPs distribution across the genome #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.1. Convert the data into long format #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    rawGeno_wzNA_pop1_long <- reshape2::melt(rawGeno_wzNA_pop1, id.vars=c("ID")) 
    rawGeno_wzNA_pop2_long <- reshape2::melt(rawGeno_wzNA_pop2, id.vars=c("ID")) 
    
    colnames(rawGeno_wzNA_pop1_long) <- c("ID", "SNPID", "Call")
    colnames(rawGeno_wzNA_pop2_long) <- c("ID", "SNPID", "Call")
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.2. Combine the long format with the genetic mapData #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    rawGeno_wzNA_pop1_long_map <- merge(mapData, rawGeno_wzNA_pop1_long, "SNPID")
    rawGeno_wzNA_pop2_long_map <- merge(mapData, rawGeno_wzNA_pop2_long, "SNPID")
  
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 5.3. Plot the missing distribution #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ## Population 1 Missing SNPs distribution
    pop1_NAs_in_genome_tmp <- rawGeno_wzNA_pop1_long_map %>%
                              mutate(countfactor=cut(pos, breaks=c(0, 20000,500000, 1000000, 10000000, 20000000, max(pos, na.rm=T)),
                              labels=c("0","0-20k","20k-500k","500k-1Mb", "1Mb-10Mb", ">20Mb")))
    
    ## Population 2 Missing SNPs distribution
    pop2_NAs_in_genome_tmp <- rawGeno_wzNA_pop2_long_map %>%
                              mutate(countfactor=cut(pos, breaks=c(0, 20000,500000, 1000000, 10000000, 20000000, max(pos, na.rm=T)),
                              labels=c("0","0-20k","20k-500k","500k-1Mb", "1Mb-10Mb", ">20Mb")))
    
        #<<<<<<<<<<<<<<<<<<<<<
        # 5.3.1. Plot the missing rate per sample population 1#####
        #<<<<<<<<<<<<<<<<<<<<<

       missingDist_pop1 <- ggplot(pop1_NAs_in_genome_tmp, aes(x=countfactor, y=ID, fill=as.factor(Call))) +
          geom_tile() + facet_grid(~ chr) + 
          scale_fill_manual(values = c('white', 'white', 'white', 'red'), na.value = "red") +
          theme(strip.text = element_text(face="bold", size=12),
          strip.background = element_rect(fill="lightblue"))+
          ggtitle(paste("Chromosome", sep="")) +
          labs(x = "Genetic position interval", y = "Strawberry accessions - Population 1", color = "Call") +
          theme(legend.position = "bottom") + guides(fill=guide_legend(title="Call")) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                               axis.text.y = element_blank(),
                               panel.border = element_rect(colour = "black", fill=NA, size=1))
    
#CH        ggsave(missingDist_pop1, filename="../Analysis_post220605/2. Imputation/missingDist_population1.png", height=6.5, width=12.8, units="in", dpi=1000)
        ggsave(missingDist_pop1, filename="missingDist_population1.png", height=6.5, width=12.8, units="in", dpi=1000)
 
      #<<<<<<<<<<<<<<<<<<<<<
      # 5.3.2. Plot the missing rate per sample - population 2 #####
      #<<<<<<<<<<<<<<<<<<<<<
      
        missingDist_pop2 <- ggplot(pop2_NAs_in_genome_tmp, aes(x=countfactor, y=ID, fill=as.factor(Call))) +
                                  geom_tile() + facet_grid(~ chr) + 
                                  scale_fill_manual(values = c('white', 'white', 'white', 'red'), na.value = "red") +
                                  theme(strip.text = element_text(face="bold", size=12),
                                        strip.background = element_rect(fill="lightblue"))+
                                  ggtitle(paste("Chromosome", sep="")) +
                                  labs(x = "Genetic position interval", y = "Strawberry accessions - Population 2", color = "Call") +
                                  theme(legend.position = "bottom") + guides(fill=guide_legend(title="Call")) + 
                                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                        axis.text.y = element_blank(),
                                        panel.border = element_rect(colour = "black", fill=NA, size=1))
                                
          ggsave(missingDist_pop2, filename="missingDist_population2.png", height=469, width=892, units="mm", dpi=1000)
      
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 6. Do the imputation per population and extract the imputed 012 matrix #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.1. Re-structure rawSNP #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       
        ## Population 1
        rawGeno_wzNA_pop1_snpReady <- rawGeno_wzNA_pop1[,2:ncol(rawGeno_wzNA_pop1),]
        rownames(rawGeno_wzNA_pop1_snpReady) <- rawGeno_wzNA_pop1$ID
        
        ## Population 2
        rawGeno_wzNA_pop2_snpReady <- rawGeno_wzNA_pop2[,2:ncol(rawGeno_wzNA_pop2),]
        rownames(rawGeno_wzNA_pop2_snpReady) <- rawGeno_wzNA_pop2$ID

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # 6.2. Impute rawSNP population 1 #####
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    rawGeno_imputed_pop1_tmp <- raw.data(data = as.matrix(rawGeno_wzNA_pop1_snpReady), 
                                     frame = "wide", base = FALSE, maf = 0.05,
                                     call.rate = 0.95, imput = TRUE, imput.type = "mean")
        
    rawGeno_imputed_pop1 <- rawGeno_imputed_pop1_tmp$M.clean
    dim(rawGeno_imputed_pop1)
        
   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   # 6.3. Impute rawSNP population 2 #####
   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
   rawGeno_imputed_pop2_tmp <- raw.data(data = as.matrix(rawGeno_wzNA_pop2_snpReady), 
                                    frame = "wide", base = FALSE, maf = 0.05,
                                    call.rate = 0.95, imput = TRUE, imput.type = "mean")
        
   rawGeno_imputed_pop2 <- rawGeno_imputed_pop2_tmp$M.clean
   dim(rawGeno_imputed_pop2)
        
   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
   # 6.4. Impute rawSNP ALL population using "knni" imputation approach #####
   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
   ## Population ALL
   rawGeno_wzNA_2064_popALL_snpReady <- rawGeno_wzNA_2064[,2:ncol(rawGeno_wzNA_2064),]
   rownames(rawGeno_wzNA_2064_popALL_snpReady) <- rawGeno_wzNA_2064$ID
        
   rawGeno_imputed_popALL_tmp <- raw.data(data = as.matrix(rawGeno_wzNA_2064_popALL_snpReady), 
                                     frame = "wide", base = FALSE, maf = 0.01,
                                     call.rate = 0.8, imput = TRUE, imput.type = "knni")

   rawGeno_imputed_popALL <- rawGeno_imputed_popALL_tmp$M.clean
   dim(rawGeno_imputed_popALL)

   ## calculate concordance percentage
   ## extract column before and after imputation
   ## popall with pop1 common rows and column
    popallpop1_col <- rawGeno_imputed_popALL[, colnames(rawGeno_imputed_popALL) %in% colnames(rawGeno_imputed_pop1)]
    popallpop1_col_tmp <- popallpop1_col[, colnames(popallpop1_col) %in% colnames(rawGeno_imputed_pop1)]
    popallpop1_col_tmp2 <- popallpop1_col_tmp[, colnames(popallpop1_col_tmp) %in% colnames(rawGeno_imputed_pop1)] 
    popallpop1_col_row <- popallpop1_col_tmp2[rownames(rawGeno_imputed_pop1),] 

    ## pop1 with popall common rows and column
    pop1popall <- rawGeno_imputed_pop1[, colnames(rawGeno_imputed_pop1) %in% colnames(popallpop1_col_row)]
    
    ## concordance rate allpopImputation vs imputed pop1
    concordancePercentageallPop_vs_pop1 <- as.data.frame(colSums(popallpop1_col_row[,seq(1,ncol(popallpop1_col_row),1)]==pop1popall[,seq(1,ncol(pop1popall),1)])/nrow(pop1popall)*100)
    concordancePercentageallPop_vs_pop1$SNPID <- rownames(concordancePercentageallPop_vs_pop1)
    colnames(concordancePercentageallPop_vs_pop1)[1] <- "conc.popallvspop1"
    rownames(concordancePercentageallPop_vs_pop1) <-NULL
    
    ## popall with pop2 common rows and column
    popallpop2_col <- rawGeno_imputed_popALL[, colnames(rawGeno_imputed_popALL) %in% colnames(rawGeno_imputed_pop2)]
    popallpop2_col_tmp <- popallpop2_col[, colnames(popallpop2_col) %in% colnames(rawGeno_imputed_pop2)]
    popallpop2_col_tmp2 <- popallpop2_col_tmp[, colnames(popallpop2_col_tmp) %in% colnames(rawGeno_imputed_pop2)] 
    popallpop2_col_row <- popallpop2_col_tmp2[rownames(rawGeno_imputed_pop2),] 
    
    ## pop2 with popall common rows and column
    pop2popall <- rawGeno_imputed_pop2[, colnames(rawGeno_imputed_pop2) %in% colnames(popallpop2_col_row)]
    
    ## concordance rate allpopImputation vs imputed pop2
    concordancePercentageallPop_vs_pop2 <- as.data.frame(colSums(popallpop2_col_row[,seq(1,ncol(popallpop2_col_row),1)]==pop2popall[,seq(1,ncol(pop2popall),1)])/nrow(pop2popall)*100)
    concordancePercentageallPop_vs_pop2$SNPID <- rownames(concordancePercentageallPop_vs_pop2)
    colnames(concordancePercentageallPop_vs_pop2)[1] <- "conc.popallvspop2"
    rownames(concordancePercentageallPop_vs_pop2) <-NULL
    
    ## plot popallpop1
    
    conc.pop1_map <- merge(concordancePercentageallPop_vs_pop1, mapData, by="SNPID")
    plot(density(conc.pop1_map$conc.popallvspop1))
    
    ## plot popallpop2
    conc.pop2_map <- merge(concordancePercentageallPop_vs_pop2, mapData, by="SNPID")
    plot(density(conc.pop2_map$conc.popallvspop2))

    
    ##
    rawGeno_imputed_popALL_tmp_Mnimp <- raw.data(data = as.matrix(rawGeno_wzNA_2064_popALL_snpReady),
                                           frame = "wide", base = FALSE, maf = 0.01,
                                           call.rate = 0.8, imput = TRUE, imput.type = "mean", outfile = "012")

    rawGeno_imputed_popALL_meanImp <- rawGeno_imputed_popALL_tmp_Mnimp$M.clean
    
    
    concordancePercentage <- colSums(rawGeno_imputed_popALL[,seq(1,ncol(rawGeno_imputed_popALL),1)]==rawGeno_imputed_popALL_meanImp[,seq(1,ncol(rawGeno_imputed_popALL_meanImp),1)])/nrow(rawGeno_imputed_popALL_meanImp)*100
    range(concordancePercentage)
    plot(concordancePercentage, pch=16, cex=1)
    
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 7. Save work space #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
      
save.image("strawberryImputation.RData")
