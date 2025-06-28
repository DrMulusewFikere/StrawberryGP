###############################################################################
# Multi-trail Pfa model
###############################################################################

###############################################################################
# 1. Set up R
###############################################################################

    ############################################################
    # 1.1. Load packages
    ############################################################

    #library(asreml, lib.loc = "C:/Users/uqchardn/Dropbox/Documents/R/R-3.6.0/Asreml_versions/Asreml4")
    library(asreml)
    library(snpReady)


    ############################################################
    # 1.2. Define custom functions
    ############################################################

        ########################################
        # 1.2.1. Estimate vcov from fa model
        ########################################
        
        fatocov.f <- function(fatocov_in)
        {
            fc.Mpsi <- diag(fatocov_in$gammas[grepl(fatocov_in$fac,names(fatocov_in$gammas),fixed=T) &
                                                  grepl(fatocov_in$E,names(fatocov_in$gammas),fixed=T) &
                                                  grepl("var",names(fatocov_in$gammas))])
            rownames(fc.Mpsi) <- colnames(fc.Mpsi) <- sub(paste(fatocov_in$E,".",sep=""),"",
                                                          sub(".var","",
                                                              sapply(
                                                                  strsplit(
                                                                      names(fatocov_in$gammas[
                                                                          grepl(fatocov_in$fac,names(fatocov_in$gammas),fixed=T) &
                                                                              grepl(fatocov_in$E,names(fatocov_in$gammas),fixed=T) &
                                                                              grepl("var",names(fatocov_in$gammas))]),
                                                                      split = "!",fixed = T),
                                                                  "[[",2)
                                                          )
            )
            
            fc.Mlam <- matrix(fatocov_in$gammas[grepl(fatocov_in$fac,names(fatocov_in$gammas),fixed=T) &
                                                    grepl(fatocov_in$E,names(fatocov_in$gammas),fixed=T) &
                                                    !grepl("var",names(fatocov_in$gammas))],
                              nrow = nrow(fc.Mpsi))
            
            rownames(fc.Mlam) <- rownames(fc.Mpsi)
            colnames(fc.Mlam) <- paste("fa",1:ncol(fc.Mlam),sep="")
            
            fc.Mvcov <- as.matrix(fc.Mlam%*%t(fc.Mlam) + fc.Mpsi)
            
            fc_out <- NULL
            fc_out$Mpsi <- as.matrix(fc.Mpsi)
            fc_out$Mlam <- as.matrix(fc.Mlam)
            fc_out$Mvcov <- as.matrix(fc.Mvcov)
            
            fc_out$Mr <- cov2cor(fc_out$Mvcov)
            
            return(fc_out)
            
        }
        
        ########################################
        # 1.2.2. Test model differences
        ########################################
        
        tmd.f <- function(S01b.v.S01_tmd_in)
        {
            ####################
            # 1.3.3.1. Input tmd format
            #
            # S01b.v.S01_tmd_in <- NULL
            # S01b.v.S01_tmd_in$analy_name
            # S01b.v.S01_tmd_in$fullmod_name
            # S01b.v.S01_tmd_in$redmod_name
            # S01b.v.S01_tmd_in$fullmod.asr
            # S01b.v.S01_tmd_in$redmod.asr
            # S01b.v.S01_tmd_in$Lest ("boundary"/"other")
            #
            ####################
            
            analysum_temp <- NULL
            analysum_temp$analy_name <- S01b.v.S01_tmd_in$analy_name
            analysum_temp$full_name <- S01b.v.S01_tmd_in$fullmod_name
            analysum_temp$red_name <- S01b.v.S01_tmd_in$redmod_name
            
            analysum_temp$fullmod.logl <- S01b.v.S01_tmd_in$fullmod.asr$loglik
            analysum_temp$redmod.logl <- S01b.v.S01_tmd_in$redmod.asr$loglik
            
            analysum_temp$fullmod.df <- sum(summary(S01b.v.S01_tmd_in$fullmod.asr)$var$bound%in%c("P","U"))
            analysum_temp$redmod.df <- sum(summary(S01b.v.S01_tmd_in$redmod.asr)$var$bound%in%c("P","U"))
            
            analysum_temp$Dlogl <- round(2*(analysum_temp$fullmod.logl - analysum_temp$redmod.logl),3)
            analysum_temp$Ddf <- analysum_temp$fullmod.df - analysum_temp$redmod.df
            
            analysum_temp$mTrialest <- S01b.v.S01_tmd_in$mTrialest
            
            if(S01b.v.S01_tmd_in$mTrialest == 'boundary')
            {
                if(analysum_temp$Ddf <1)
                {
                    analysum_temp$pDlogl = (1 - pchisq(analysum_temp$Dlogl,1))/2
                } else
                {
                    analysum_temp$pDlogl <- (1 - pchisq(analysum_temp$Dlogl,analysum_temp$Ddf))/2
                }
            } else
            {
                if(analysum_temp$Ddf <1)
                {
                    analysum_temp$pDlogl = (1 - pchisq(analysum_temp$Dlogl,1))
                } else
                {
                    analysum_temp$pDlogl <- (1 - pchisq(analysum_temp$Dlogl,analysum_temp$Ddf))
                }
                
            }
            
            return(as.data.frame(analysum_temp))
        }
    
        ########################################
        # 1.2.3. Make a positive definite
        ########################################

        make.positive.definite.f <- function(M)
        {
          tol=1e-3
          eig <- eigen(M, symmetric=TRUE)
          rtol <- tol * eig$values[1]
          if(min(eig$values) < rtol)
          {
            vals <- eig$values
            vals[vals < rtol] <- rtol
            srev <- eig$vectors %*% (vals * t(eig$vectors))
            dimnames(srev) <- dimnames(M)
            return(srev)
          } else
          {
            return(M)
          }
        }
        
###############################################################################
# 2. Read in Gfa1_RDATA
###############################################################################

#load("../../../../Post_201203/Analysis_post210218/4. Fit models/4.2. MT analy/Gfa1RDATA.RData")
load("Gfa1RDATA.RData")

###############################################################################
# 3. Prepare GRM 
###############################################################################

# grmpd = GRM + min(GRM)
# grmcor = cov2cor(grmpd)
# grmdist = cor2dist(grmcor)
# 
# mds.stuff = cmdscale(grmdist,k=2063, eig=TRUE, x.ret=TRUE)

# U_S01 <- mds.stuff$points
# D_S01 <- mds.stuff$eig
#         
# Ut_S01 <- t(U_S01)
#        
# PCAgrm_S01 <- t(U_S01) %*% D_S01 %*% Ut_S01

        
  
#### G_S01 = UDUt (refer: DOI 10.1007/s00122-013-2255-x)

pca_S01 <- prcomp(GRM, center = TRUE,scale. = FALSE)
    
### U_S01 = nx(n-1): where U is an n ? (n  - 1) matrix of the eigen vectors of G with Ui the column i (i = 1, 2, ., n-1)
U_S01 = pca_S01$rotation
     
MD <- matrix(0,nrow = length(pca_S01$sdev),ncol = length(pca_S01$sdev))
diag(MD) <- pca_S01$sdev^2

### D_S01 = (n-1)x(n-1) is diagonal matrix with each diagonal element representing eigenvalues
    
D_S01 <- eigen(cor(scale(U_S01)))
Ut_S01 <- t(U_S01)
       
### After eigenvalue decomposition of PCAgrm_S01 as

PCAgrm_S01 <- U_S01 %*% diag(D_S01$values) %*% Ut_S01
hist(diag(PCAgrm_S01))

test <- solve(PCAgrm_S01)

grm_Pfa_pos <-  make.positive.definite.f(PCAgrm_S01)
    
giv_Pfa <- solve(grm_Pfa_pos)

hist(giv_Pfa)

    ############################################################
    # 3.1. Factorize
    ############################################################

    adata_S01$PAID <- factor(as.character(adata_S01$ID, levels = rownames(giv_Pfa)))
     
        
###############################################################################
# 4. Fit Pfa1
###############################################################################
        

    ############################################################
    # 4.1. Fit model
    ############################################################
        
    asreml.options(dense = ~ vm(PAID, giv_Pfa), ai.sing=TRUE)
        
    MT_Pfa1.asr <- update(MT_Gfa1.asr,
                              random = ~ 
                                at(T,c('P4')):B:D:X +
                                fa(AE23,1):vm(PAID, giv_Pfa) +
                                at(CE23,c('E8','M4','M5','N8','P','W8')):idv(C),
                          maxiter = 50,
                          workspace="16gb",
                          data=adata_S01)

###############################################################################
# 5. Save Pfa1.RDATA
###############################################################################    
    
save.image("Pfa1RDATA.RData")
 
        