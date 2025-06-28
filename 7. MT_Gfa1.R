###############################################################################
# Multi-trail single pop GRM 
###############################################################################

###############################################################################
# 1. Set up R
###############################################################################

    ############################################################
    # 1.1. Load packages
    ############################################################

    library(asreml, lib.loc = "C:/Users/uqchardn/Dropbox/Documents/R/R-3.6.0/Asreml_versions/Asreml4")
    library(asreml)

###############################################################################
# 2. Input MT_Prelim_RDATA
###############################################################################

load("MT_Prelim_RDATA_CH210609.RData")

###############################################################################
# 4. Fit Gfa1
###############################################################################
        
    ############################################################
    # 4.1. Prepare GRM 
    ############################################################
    
    ### The refined data "sdata_snpReady_filt"

    sdata.Z <- sdata_snpReady_filt
    colnames(sdata.Z) <- NULL
    rownames(sdata.Z) <-NULL
    M <-matrix(as.numeric(unlist(sdata.Z)),nrow=nrow(sdata.Z))
    M <- as.matrix(M) 
    
    #############################
    # 5.6.1.1. Calc. allele Freq and GRM for Z population
    #############################
    
    P<-matrix(colMeans(0.5*M),1,ncol(M))
    T<-matrix(1,nrow(M),1)
    S<-M-2*(T%*%P)
    Z1tZ1<-S%*%t(S)
    Q<-1-P
    f<-2*sum(P*Q)
    Z<-Z1tZ1*(1/f)

    GRM <- Z
    colnames(GRM) <- rownames(sdata_snpReady_filt)
    rownames(GRM) <- rownames(sdata_snpReady_filt)
    
    #############################
    # 5.6.1.1. Make a positive definite and estimate GIV
    #############################
    
    GRM_pos <- make.positive.definite.f(GRM)
    giv_Gfa <- solve(GRM_pos)
    
    hist(diag(giv_Gfa))
    
    attr(giv_Gfa, "INVERSE") = TRUE
    
    adata_S01$GAID <- factor(adata_S01$ID, levels = rownames(giv_Gfa))
    
    ############################################################
    # 4.2. Fit model
    ############################################################
        
    #asreml.options(dense = ~ vm(GAID, giv_Gfa))
        
    MT_Gfa1.asr <- asreml(fixed = lnSSC ~ 
                                           -1 + T + 
                                           at(T,c('E8','P4','P5')):B + 
                                           at(T,c('B1','C1')):Y +
                                           at(T,c('N8')):H,
                           random = ~ 
                                       at(T,c('P4')):B:D:X +
                                       fa(AE23,1):vm(GAID, giv_Gfa) +
                                       at(CE23,c('E8','M4','M5','N8','P','W8')):idv(C),
                           residual = ~ 
                                        dsum(~idv(units)| T),
                           na.action = na.method(y = "include", x = "include"),
                           maxiter = 20,
                           data = adata_S01,
                           ai.loadings = TRUE,
                           workspace = "16gb")
    
    
###############################################################################
# 5. Save Gfa1.RDATA
###############################################################################    
    
save.image("Gfa1RDATA_210701.RData")
    
###############################################################################
# 6. Review of fit
###############################################################################
    
options(scipen = 999)
        
summary(MT_Gfa1.asr)$var

blups <- MT_Gfa1.asr$coefficients$random

bGAID <- data.frame(blups[grepl('GAID',rownames(blups)) & !grepl('Comp',rownames(blups)),])
colnames(bGAID) <- 'bGAID'
bGAID$AE23 <- sub('fa(AE23, 1)_','',sapply(strsplit(rownames(bGAID),':'),'[',1),fixed = T)
bGAID$GAID <- sub('vm(GAID, giv_Gfa)_','',sapply(strsplit(rownames(bGAID),':'),'[',2),fixed = T)
rownames(bGAID) <- NULL 
hist(bGAID$bGAID[bGAID$AE23 =='B1'])

GAID_at_B1 <- as.character(adata_S01$GAID[adata_S01$L == 'B'])

hist(bGAID$bGAID[bGAID$AE23 =='B1'])
hist(bGAID$bGAID[bGAID$AE23 =='B1' & bGAID$GAID%in%GAID_at_B1])
hist(bGAID$bGAID[bGAID$AE23 =='B1' & !bGAID$GAID%in%GAID_at_B1])
