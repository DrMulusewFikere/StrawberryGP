#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 1. Genomic prediction using PCA approach (eigenvalue decomposition) #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#### MFNote: Not ready for Craig's review (16/06/2022)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 2. Define path and library  #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(asreml)
library(reshape2)
library(tidyr)
library(dplyr)

# setwd("/Users/uqmfiker/Dropbox/Post_201203/Data_post220524")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 3. Load workspace from section 1 ("Structure") ("and 2 ("EstimateGRM") #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#CH load("../Analysis_post220605/3. Estimate GRMs/strawberryEstimateGRM.RData")
load("strawberryEstimateGRM.RData")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 4. Fit model #####
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 4.1. Test PCA approach (it works!) #####
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      ## Test PCA APPROACH (WORKS!)
      
    Guo_mat_pheno$L <- factor(Guo_mat_pheno$L)
    Guo_mat_pheno$V <- factor(Guo_mat_pheno$V)

      testpca.asr <- asreml(lnSSC ~ L,
                            random=~grp(V),
                            residual=~idv(units),
                            group=list(V=1:2064),
                            na.action=na.method(y='include'),
                            workspace=1e08,
                            data=Guo_mat_pheno)

      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 4.2a.  PCA full model #####
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      Guo_mat_pheno$L <- factor(Guo_mat_pheno$L )
      Guo_mat_pheno$T <- factor(Guo_mat_pheno$T)
      Guo_mat_pheno$S <- factor(Guo_mat_pheno$S)
      Guo_mat_pheno$PY <- factor(Guo_mat_pheno$PY)
      Guo_mat_pheno$B <- factor(Guo_mat_pheno$B)
      Guo_mat_pheno$D <- factor(Guo_mat_pheno$D)
      Guo_mat_pheno$X <- factor(Guo_mat_pheno$X)
      Guo_mat_pheno$SP <- factor(Guo_mat_pheno$SP)
      Guo_mat_pheno$P <- factor(Guo_mat_pheno$P)
      Guo_mat_pheno$U <- factor(Guo_mat_pheno$U)
      Guo_mat_pheno$ID <- factor(Guo_mat_pheno$ID)
      Guo_mat_pheno$PAID <- factor(Guo_mat_pheno$PAID)
 
      Guo_mat_pheno$AE23 <- factor(Guo_mat_pheno$AE23)
      Guo_mat_pheno$CE23 <- factor(Guo_mat_pheno$CE23)
      
      GAID <- Guo_mat_pheno$PAID
      
      ## MFNote: not tested yet
      pcaFull.asr <- asreml(fixed = lnSSC ~ 
                           -1 + T + 
                           at(T,c('E8','P4','P5')):B + 
                           at(T,c('B1','C1')):Y +
                           at(T,c('N8')):H,
                         random = ~grp(V),
                         ~ at(T,c('P4')):B:D:X +
                           fa(AE23,1):vm(GAID, list(V=1:2064)) +
                           at(CE23,c('E8','M4','M5','N8','P','W8')):idv(C),
                         residual = ~ 
                           dsum(~idv(units)| T),
                         group=list(V=1:2064),
                         na.action = na.method(y = "include", x = "include"),
                         maxiter = 20,
                         data = Guo_mat_pheno,
                         ai.loadings = TRUE,
                         workspace = "16gb")

      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # 4.2b.  PCA model with Gfa1 update #####
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
      ## MFNote: not tested yet
      
      MT_Pfa1.asr <- update(MT_Gfa1.asr,
                            random = ~ 
                              at(T,c('P4')):B:D:X +
                              fa(AE23,1):vm(GAID, GuoGRM) +
                              at(CE23,c('E8','M4','M5','N8','P','W8')):idv(C),
                            maxiter = 50,
                            workspace="16gb",
                            data=Guo_mat_pheno)
    


