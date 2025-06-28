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

load("Gfa1RDATA_210701.RData")

###############################################################################
# 4. Fit Gfa2
###############################################################################
    ############################################################
    # 4.2. Fit model
    ############################################################
        
    #asreml.options(dense = ~ vm(GAID, giv_Gfa))
        
    MT_Gfa2.asr <- update(MT_Gfa1.asr,
                           random = ~ 
                                       at(T,c('P4')):B:D:X +
                                       fa(AE23,1):vm(GAID, giv_Gfa) +
                                       at(CE23,c('E8','M4','M5','N8','P','W8')):idv(C),
                           )
    
###############################################################################
# 5. Save Gfa1.RDATA
###############################################################################    
    
save.image("Gfa2RDATA_210701.RData")
    
# ###############################################################################
# # 6. Review of fit
# ###############################################################################
#     
# options(scipen = 999)
#         
# summary(MT_Gfa1.asr)$var
# 
# blups <- MT_Gfa1.asr$coefficients$random
# 
# bGAID <- data.frame(blups[grepl('GAID',rownames(blups)) & !grepl('Comp',rownames(blups)),])
# colnames(bGAID) <- 'bGAID'
# bGAID$AE23 <- sub('fa(AE23, 1)_','',sapply(strsplit(rownames(bGAID),':'),'[',1),fixed = T)
# bGAID$GAID <- sub('vm(GAID, giv_Gfa)_','',sapply(strsplit(rownames(bGAID),':'),'[',2),fixed = T)
# rownames(bGAID) <- NULL 
# hist(bGAID$bGAID[bGAID$AE23 =='B1'])
# 
# GAID_at_B1 <- as.character(adata_S01$GAID[adata_S01$L == 'B'])
# 
# hist(bGAID$bGAID[bGAID$AE23 =='B1'])
# hist(bGAID$bGAID[bGAID$AE23 =='B1' & bGAID$GAID%in%GAID_at_B1])
# hist(bGAID$bGAID[bGAID$AE23 =='B1' & !bGAID$GAID%in%GAID_at_B1])
