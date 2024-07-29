#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Genomic relationship using IMPUTED DATA #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. Define path and library #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(asreml)
library(AGHmatrix)
library(ASRgenomics)
library(snpReady)

# setwd("/Users/uqmfiker/Dropbox/Post_201203/Data_post220524")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Load workspace from section 2 "strawberyImputation" #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load("strawberryStructure.RData")
load("strawberryImputation.RData")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Re-structure data for estimating GRM #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 4.1. ALL population GRM #####
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      rawGeno_imputed_popALL_tmp <- rawGeno_imputed_popALL
      
          #<<<<<<<<<<<<<<<<<<<
          # 4.1.1 Generate GRM for all population
          #<<<<<<<<<<<<<<<<<<<
          
          colnames(rawGeno_imputed_popALL_tmp) <- NULL
          rownames(rawGeno_imputed_popALL_tmp) <-NULL
          P_all_tmp <- matrix(as.numeric(unlist(rawGeno_imputed_popALL_tmp)),nrow=nrow(rawGeno_imputed_popALL_tmp))
          P_all_tmp <- as.matrix(P_all_tmp) 
          
          #<<<<<<<<<<<<<<<<<<
          # 4.1.2 Build grm_all for population
          #<<<<<<<<<<<<<<<<<<
    
          P_all<-matrix(colMeans(0.5*P_all_tmp),1,ncol(P_all_tmp))
          T_all<-matrix(1,nrow(P_all_tmp),1)
          S_all<-P_all_tmp-2*(T_all%*%P_all)
          ZalltZall<-S_all%*%t(S_all)
          Q_all<-1-P_all
          f_all<-2*sum(P_all*Q_all)
          grm_all <-ZalltZall*(1/f_all)
          
          rownames(grm_all) <- rownames(rawGeno_imputed_popALL)
          colnames(grm_all) <- rownames(rawGeno_imputed_popALL)
          
          ## Summary stat- below diag
          summary(c(grm_all[upper.tri(grm_all)]))
          summary(diag(grm_all))
          
          ## check GRM with R software
          G <- G.matrix(M = rawGeno_imputed_popALL, method = "VanRaden", format = "wide") 
          G$Ga [1:5, 1:5]   
          
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 4.2. Population specific GRM #####
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      uID_pop1 <- rownames(rawGeno_imputed_pop1)
      rawGeno_imputed_pop1 <- rawGeno_imputed_popALL[match(uID_pop1, rownames(rawGeno_imputed_popALL)),]

          #<<<<<<<<<<<<<<<<<<<<<
          # 4.2.1. Population - pop11 GRM #####
          #<<<<<<<<<<<<<<<<<<<<<
          
          rawGeno_imputed_pop1_tmp <- rawGeno_imputed_pop1

          ## Re-structure input matrix
      
          colnames(rawGeno_imputed_pop1_tmp) <- NULL
          rownames(rawGeno_imputed_pop1_tmp) <-NULL
          P_11_tmp <- matrix(as.numeric(unlist(rawGeno_imputed_pop1_tmp)),nrow=nrow(rawGeno_imputed_pop1_tmp))
          P_11_tmp <- as.matrix(P_11_tmp) 
      
              #<<<<<<<<<<<<<<
              # 4.3.1.1. Build grm for pop 11
              #<<<<<<<<<<<<<<
              
              P_11<-matrix(colMeans(0.5*P_11_tmp),1,ncol(P_11_tmp))
              T_11<-matrix(1,nrow(P_11_tmp),1)
              S_11<-P_11_tmp-2*(T_11%*%P_11)
              Z11tZ11<-S_11%*%t(S_11)
              Q_11<-1-P_11
              f_11<-2*sum(P_11*Q_11)
              grm_11 <-Z11tZ11*(1/f_11)

              ## Summary stat- below diag
              summary(c(grm_11[upper.tri(grm_11)]))
              summary(diag(grm_11))
          
          #<<<<<<<<<<<<<<<<<<<<<
          # 4.2.2. Population - pop22 GRM #####
          #<<<<<<<<<<<<<<<<<<<<<

          uID_pop2 <- rownames(rawGeno_imputed_pop2)
          rawGeno_imputed_pop2 <- rawGeno_imputed_popALL[match(uID_pop2, rownames(rawGeno_imputed_popALL)),]
          
          rawGeno_imputed_pop2_tmp <- rawGeno_imputed_pop2

          ## Re-structure input matrix
      
          colnames(rawGeno_imputed_pop2_tmp) <- NULL
          rownames(rawGeno_imputed_pop2_tmp) <-NULL
          P_22_tmp <- matrix(as.numeric(unlist(rawGeno_imputed_pop2_tmp)),nrow=nrow(rawGeno_imputed_pop2_tmp))
          P_22_tmp <- as.matrix(P_22_tmp) 
              
              #<<<<<<<<<<<<<<
              # 4.2.2.1. Build grm for pop 22
              #<<<<<<<<<<<<<<
              
              P_22<-matrix(colMeans(0.5*P_22_tmp),1,ncol(P_22_tmp))
              T_22<-matrix(1,nrow(P_22_tmp),1)
              S_22<-P_22_tmp-2*(T_22%*%P_22)
              Z22tZ22<-S_22%*%t(S_22)
              Q_22<-1-P_22
              f_22<-2*sum(P_22*Q_22)
              grm_22 <-Z22tZ22*(1/f_22)
              
              ## Summary stat- below diag
              summary(c(grm_22[upper.tri(grm_22)]))
              summary(diag(grm_22))
      
          #<<<<<<<<<<<<<<<<<<<<<
          #  4.2.3. Population - pop12 and pop21 GRM #####
          #<<<<<<<<<<<<<<<<<<<<<
              
          grm12 <- (S_11%*%t(S_22))*(1/(sqrt(f_11)*sqrt(f_22)))
          grm21 <- (S_22%*%t(S_11))*(1/(sqrt(f_11)*sqrt(f_22)))
          
          com12 <- cbind(grm_11,grm12)
          com21 <- cbind(grm21, grm_22)
          grm_popSpecific <- rbind(com12, com21)
          
          summary(c(grm_popSpecific[upper.tri(grm_popSpecific)]))
          summary(diag(grm_popSpecific))
              
          #<<<<<<<<<<<<<<<<<<<<<
          #  4.2.4. final step of population specific GRM #####
          #<<<<<<<<<<<<<<<<<<<<<
          
          rownames(grm_popSpecific) <- rownames(rawGeno_imputed_popALL)
          colnames(grm_popSpecific) <- rownames(rawGeno_imputed_popALL)

          #<<<<<<<<<<<<<<<<<<<<<
          # 4.3.5. Plot allele frequency between pop1&2 #####
          #<<<<<<<<<<<<<<<<<<<<<
          
          ## Global vs pop1
          plot(P_all,P_11, pch=16, cex=0.6, xlab="GlobalPop", ylab="pop1")
          ## Global vs pop2
          plot(P_all,P_22, pch=16, cex=0.6, xlab="GlobalPop", ylab="pop2")
          ## pop1 vs pop2
          plot(P_11,P_22, pch=16, cex=0.6, xlab="pop1", ylab="pop2")
 
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 4.3. PCA approach (eigen vector decomposition) GRM #####
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<         
          
      ## Gou et al., 2014 GRM = UDV, eigenvalue decomposition
     
      G_all_blend <- G.tuneup(G=grm_all, pblend=0.05, blend=TRUE)$Gb
      peig <- eigen(G_all_blend)
      ls(peig)
      eigenv <- peig$values # eigenvalues (vector of diag(D))
      loads <- peig$vectors # loadings (matrix U)
      Guo_mat <- loads %*% diag(sqrt(eigenv))

      # Preparing for ASReml-R
      colnames(Guo_mat) <- paste0('V', c(1:2064))
      rownames(Guo_mat) <- rownames(G_all_blend)
      Guo_mat <- data.frame(GAID = row.names(Guo_mat), Guo_mat)
      rownames(Guo_mat) <- NULL
      
      Guo_mat_pheno <- merge(Guo_mat, refinedPheno_2064, by="GAID")
      Guo_mat_pheno$PAID <- Guo_mat_pheno$GAID
      Guo_mat_pheno$GAID <- NULL

      # MFNote: Now Fit "Guo_mat_pheno" instead of GIV into ASReml as groups

      
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 5. See predPCABLUP (section 4) to see how PCA is fitted #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
      

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 6. Save work space #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
          
save.image("strawberryEstimateGRM.RData")

