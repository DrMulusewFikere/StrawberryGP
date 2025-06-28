###############################################################################
# Multi-trail Wfa model
###############################################################################

###############################################################################
# 1. Set up R
###############################################################################

    ############################################################
    # 1.1. Load packages
    ############################################################

    library(asreml, lib.loc = "C:/Users/uqchardn/Dropbox/Documents/R/R-3.6.0/Asreml_versions/Asreml4")
    library(asreml)
    library(snpReady)
 
  
###############################################################################
# 2. Read in Gfa1_RDATA
###############################################################################

load("../../../../Post_201203/Analysis_post210218/4. Fit models/4.2. MT analy/Gfa1RDATA.RData")

###############################################################################
    # 2.1. Prepare GRM 
###############################################################################

    sdata.all_012 <- sdata_snpReady_filt
    
    ###########################
    # 2.2. cluster analysis
    ###########################
    
    sdata.all_012 <- as.matrix( sdata.all_012)
    sdata.all_012_mat <-  sdata.all_012/2
    Mfdist = as.matrix(dist(sdata.all_012_mat))
    Mds = cmdscale(Mfdist, eig=TRUE, x.ret=TRUE)
    
    mds.var.per = round(Mds$eig/sum(Mds$eig)*100, 1)
    mds.values = Mds$points
    mds.data = data.frame(Sample=rownames(mds.values), X=mds.values[,1], Y=mds.values[,2])
    scree_plot_dat = data.frame(PCoord = seq(1, length(mds.var.per)), per.var.explained = mds.var.per)
    
    mds.data$Sample<- NULL
    colnames(mds.data) [1] <- "PC1"
    colnames(mds.data)[2] <- "PC2"
    mds.data$Sample<- NULL
    k2 <- kmeans(mds.data, centers = 2, nstart = 25)

    ###########################
    # 2.3. Combine the cluster ID list with the pdata
    ###########################
    
    ID_cluster_temp <-as.data.frame(k2$cluster)
    sdata_select_S01_ID <- data.frame(ID = row.names(sdata_snpReady_filt), sdata_snpReady_filt)
    rownames(ID_cluster_temp) <- sdata_select_S01_ID[,1]
    
    ID_cluster <- data.frame(ID = row.names(ID_cluster_temp), ID_cluster_temp)
    colnames(ID_cluster) [2] <- "cluster"
    rownames(ID_cluster) <-NULL
    
    adata_S01 <- as.data.frame(unclass(adata_S01))
    
    pdata_in_temp <- merge (adata_S01, ID_cluster, by= "ID")
    
    ###########################
    # 2.4. Extract IDs in cluster I & II
    ###########################
    
    Z11_pdata <- pdata_in_temp[pdata_in_temp$cluster == "1",]
    Z22_pdata <- pdata_in_temp[pdata_in_temp$cluster == "2",]
    
    ###########################
    # 2.5. Subset sdata by cluster (I & II)
    ###########################
    sdata_select_S01_ID_mod <- sdata_select_S01_ID[,-1]
    sdata_S01_mod <- data.frame(ID = row.names(sdata_select_S01_ID_mod), sdata_select_S01_ID_mod)
    rownames(sdata_S01_mod) <-NULL
 
    sdata_S01_temp <- merge (sdata_S01_mod, ID_cluster, by= "ID")
  
    Z11_sdata_temp <- sdata_S01_temp[sdata_S01_temp$cluster == "1",]
    Z11_sdata_temp$cluster <-NULL
    Z11_sdata <- Z11_sdata_temp[,2:ncol(Z11_sdata_temp)]
    rownames(Z11_sdata) <- Z11_sdata_temp[,1]
    
    Z22_sdata_temp <- sdata_S01_temp[sdata_S01_temp$cluster == "2",]
    Z22_sdata_temp$cluster <-NULL
    Z22_sdata <- Z22_sdata_temp[,2:ncol(Z22_sdata_temp)]
    rownames(Z22_sdata) <- Z22_sdata_temp[,1]
    
    Z11_sdata$cluster <- NULL
    Z22_sdata$cluster <- NULL

    ###########################
    # 2.6. Construct the GRM
    ###########################
    
        ########################################
        # 2.6.1. GRM for cluster I
        ########################################

        sdata.Z11 <- Z11_sdata
        colnames(sdata.Z11) <- NULL
        rownames(sdata.Z11) <-NULL
        M11 <- matrix(as.numeric(unlist(sdata.Z11)),nrow=nrow(sdata.Z11))
        
            #############################
            # 2.6.1.1. Calc. allele Freq and GRM for Z11 population
            #############################
            
            P11<-matrix(colMeans(0.5*M11),1,ncol(M11))
            T11<-matrix(1,nrow(M11),1)
            S11<-M11-2*(T11%*%P11)
            Z1tZ1<-S11%*%t(S11)
            Q11<-1-P11
            f11<-2*sum(P11*Q11)
            Z11<-Z1tZ1*(1/f11)
    
        ########################################
        # 2.6.2. GRM for cluster II 
        ########################################
        
        sdata.Z22 <- Z22_sdata 
        colnames(sdata.Z22) <- NULL
        rownames(sdata.Z22) <-NULL
        M22 <-matrix(as.numeric(unlist(sdata.Z22)),nrow=nrow(sdata.Z22))
        
            #############################
            # 2.6.2.1. Calc. allele Freq and GRM for Z22 population
            #############################
    
            P22<-matrix(colMeans(0.5*M22),1,ncol(M22))
            T22<-matrix(1,nrow(M22),1)
            S22<-M22-2*(T22%*%P22)
            Z2tZ2<-S22%*%t(S22)
            Q11<-1-P22
            f22<-2*sum(P22*Q11)
            Z22<-Z2tZ2*(1/f22)
            
            summary(c(Z22[upper.tri(Z22)]))
            summary(diag(Z22))
    
        ########################################
        # 2.6.3. GRM for Z12 and Z21 Matrix
        ########################################
    
        subZ12 <- (S11%*%t(S22))*(1/(sqrt(f11))*sqrt(f22))
        subZ21 <- (S22%*%t(S11))*(1/(sqrt(f11)*sqrt(f22)))
        
        ########################################
        # 2.6.4. Make a complete GRM matrix
        ########################################
        
        com1 <- cbind(Z11,subZ12)
        com2 <- cbind(subZ21, Z22)
        WGRM <- rbind(com1, com2)
        
        #summary(c(GRMnew[upper.tri(WGRM)]))
        summary(diag(WGRM))
        rownames(WGRM) <- sdata_S01_temp[,1]
        colnames(WGRM) <- sdata_S01_temp[,1]
        
    ###########################
    # 2.7. Build the g inverse
    ###########################
        
    WGRM_pd <- make.positive.definite.f(WGRM)
    WGRM.giv <- solve(WGRM_pd)
    hist(diag(WGRM.giv))

    ###########################
    # 2.8. Factorize
    ###########################

    adata_S01$WAID <- factor(as.character(adata_S01$ID, levels = rownames(WGRM.giv))) 
        
###############################################################################
# 3. Fit Wfa2
###############################################################################
        

    ############################################################
    # 4.1. Fit model
    ############################################################
        
    asreml.options(dense = ~ vm(WAID, WGRM.giv))
        
    MT_Wfa2.asr <- update(MT_Gfa1.asr,
                              random = ~ 
                                at(T,c('P4')):B:D:X +
                                fa(AE23,2):vm(WAID, WGRM.giv) +
                                at(CE23,c('E8','M4','M5','N8','P','W8')):idv(C),
                          maxiter = 100,
                          workspace="16gb",
                          data= adata_S01)

    ## Finally this converge

    MT_Wfa2_up8.asr <- update( MT_Wfa2_up7.asr, maxiter=100, extra = 50)
    
###############################################################################
# 5. Save Wfa2.RDATA
###############################################################################    
    
save.image("Wfa2RDATA.RData")
 
        