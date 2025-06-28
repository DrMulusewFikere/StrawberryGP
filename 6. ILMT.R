###############################################################################
# Individual trial analyses ln transformed obs
###############################################################################

###############################################################################
# 1. Set up R
###############################################################################

    ############################################################
    # 1.1. Load packages
    ############################################################

    library(asreml, lib.loc = "C:/Users/uqchardn/Dropbox/Documents/R/R-3.6.0/Asreml_versions/Asreml4")
    library(asreml)
    library(sommer)
    library(asremlPlus)
    library(gtools)

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
        
        tmd.f <- function(tmd_in)
        {
            ####################
            # 1.3.3.1. Input tmd format
            #
            # tmd_in <- NULL
            # tmd_in$analy_name
            # tmd_in$fullmod_name
            # tmd_in$redmod_name
            # tmd_in$fullmod.asr
            # tmd_in$redmod.asr
            # tmd_in$Lest ("boundary"/"other")
            #
            ####################
            
            analysum_temp <- NULL
            analysum_temp$analy_name <- tmd_in$analy_name
            analysum_temp$full_name <- tmd_in$fullmod_name
            analysum_temp$red_name <- tmd_in$redmod_name
            
            analysum_temp$fullmod.logl <- tmd_in$fullmod.asr$loglik
            analysum_temp$redmod.logl <- tmd_in$redmod.asr$loglik
            
            analysum_temp$fullmod.df <- sum(summary(tmd_in$fullmod.asr)$var$bound%in%c("P","U"))
            analysum_temp$redmod.df <- sum(summary(tmd_in$redmod.asr)$var$bound%in%c("P","U"))
            
            analysum_temp$Dlogl <- round(2*(analysum_temp$fullmod.logl - analysum_temp$redmod.logl),3)
            analysum_temp$Ddf <- analysum_temp$fullmod.df - analysum_temp$redmod.df
            
            analysum_temp$test <- tmd_in$test
            
            if(tmd_in$test == 'boundary')
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
    
###############################################################################
# 2. Prepare data
###############################################################################

    ############################################################
    # 2.1. input data
    ############################################################
        
    ##### pdata
    pdata_in_temp <- read.csv("../../../0. Data/data_Ubal_post210120.csv", header = TRUE, stringsAsFactors = FALSE)
    
    ##### remove very low values
    pdata_in_temp[pdata_in_temp$SSC < 4 & !is.na(pdata_in_temp$SSC),]
    pdata_in <- pdata_in_temp
    pdata_in$SSC[pdata_in$SSC < 4] <- NA
    
    pdata_in$lnSSC <- log(pdata_in$SSC)*100
    
    ##### sdata
    sdata_in <- read.csv("../../../0. Data/sdata_imputed_2064.csv",header = T,stringsAsFactors = F)
    
    sdata <- sdata_in[,2:ncol(sdata_in)]
    rownames(sdata) <- sdata_in[,1]    
    
###############################################################################
# 3. M
###############################################################################
    
    ###########################################################
    # 3.1. Prepare data
    ###########################################################
    
        ########################################
        # 3.1.1. Select pdata
        ########################################
        
        sort(unique(pdata_in$L))
    
        pdata_M <- pdata_in[pdata_in$L%in%c("M"),]
        opdata_M <- pdata_M[order(pdata_M$U,pdata_M$Y,pdata_M$H),]
        
        ##### Make list of accession in pdata
        ID_M <- unique(opdata_M$ID)
        
        ########################################
        # 3.1.2. Select sdata and construct giv
        ########################################
        
        ##### Select sdata_B and construct giv
        sdata_M <- sdata[ID_M,]
        
        grm_M <-A.mat(sdata_M,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                      n.core=1,shrink=T,return.imputed=FALSE, ploidy=2)
        
        giv_M <- solve(grm_M)
        hist(diag(giv_M))
        attr(giv_M, "INVERSE") = TRUE
        
        ########################################
        # 3.1.3. Set initial factors
        ########################################
        
        adata_M <- opdata_M
        ##### Block
        adata_M$B <- factor(adata_M$B)
        
        ##### Row
        adata_M$X <- factor(adata_M$X)

        ##### Plot
        adata_M$P <- factor(adata_M$P)
        
        ##### Genetic factors
        adata_M$AID <- factor(adata_M$ID, levels = rownames(giv_M))
        
        ##### Genetic factors
        adata_M$C <- factor(adata_M$ID)
        
        ##### U = L,P,
        adata_M$U <- factor(adata_M$U)
        
        ##### Assessment Year
        adata_M$Y <- factor(adata_M$Y)
        
        ##### H = harvest
        adata_M$H <- factor(adata_M$H)
        
        ##### L = location
        adata_M$L <- factor(adata_M$L)
        
        ##### T = TRIAL
        adata_M$T <- factor(adata_M$T)
        
        ############################################################
        # 3.1.5. Review data
        ############################################################
        
        str(adata_M)
        
        tapply(adata_M$lnSSC,adata_M$Y,FUN = mean,na.rm = T)
        tapply(adata_M$lnSSC,adata_M$Y,FUN = var,na.rm = T)
    
        table(adata_M$ID,adata_M$T)
        
    ###########################################################
    # 3.2. Fit Full model (Ma)
    ###########################################################
    
        ########################################
        # 3.2.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_M))
     
        Ma.asr <- asreml(lnSSC ~ 1 + T, 
                         random = ~ vm(AID, giv_M) + vm(AID, giv_M):idv(T) + idv(C) + id(C):idv(T),
                         residual = ~ dsum(~idv(U)| T, levels= c("M4", "M5")),
                         na.action = na.method(y = "include", x = "include"),
                         data = adata_M,
                         ai.loadings = TRUE,
                         workspace = "16gb")

        ########################################
        # 3.2.2. Review Fit
        ########################################
        
            ####################
            # 3.2.2.2. Model parameters
            ####################
            
            wald(Ma.asr)
            summary(Ma.asr)$var
            
            ####################
            # 3.2.2.3. Collate results
            ####################
            
            Ma_analysum <- data.frame("Model"= "Ma")
            Ma_analysum$Scale <- 'LN'
            Ma_analysum$Status <- 'FULL'
            Ma_analysum$conv <- Ma.asr$converge
            
            Ma_analysum$prT <- round(wald(Ma.asr)[rownames(wald(Ma.asr))=='T',4],3) 
            
            Ma_analysum$vA <- round(summary(Ma.asr)$var[!grepl("vm(AID, giv_M):T",names(Ma.asr$vparameters),fixed=T) & grepl("vm(AID, giv_M)", names(Ma.asr$vparameters),fixed=T),"component"],3)
            Ma_analysum$vAxT <- round(summary(Ma.asr)$var[grepl("vm(AID, giv_M):T",names(Ma.asr$vparameters),fixed=T),"component"],3)
            
            Ma_analysum$vC <- round(summary(Ma.asr)$var[grepl("C!C",names(Ma.asr$vparameters),fixed=T),"component"],3)
            Ma_analysum$vCxT <- round(summary(Ma.asr)$var[grepl("C:T",names(Ma.asr$vparameters),fixed=T),"component"],3)
            
            Ma_analysum$vM4U <- round(summary( Ma.asr)$var[grepl("M4!U",names( Ma.asr$vparameters),fixed=T),"component"],3)
            Ma_analysum$vM5U <- round(summary( Ma.asr)$var[grepl("M5!U",names( Ma.asr$vparameters),fixed=T),"component"],3)
            
            Ma_analysum$logl <- Ma.asr$loglik
            Ma_analysum$df <- sum(summary(Ma.asr)$var$bound%in%c("P","U"))
            
            Ma_analysum$AIC <- 2*(Ma_analysum$df - Ma_analysum$logl)
        
    ###########################################################
    # 3.3. Fit reduced model (Ma - A:T)
    ###########################################################
        
        ########################################
        # 3.3.1. Fit model
        ########################################
            
        asreml.options(dense = ~ vm(AID, giv_M))
        
        Mb.asr <- update(Ma.asr,
                         random = ~ vm(AID, giv_M) + idv(C) + id(C):idv(T),)
        
        ########################################
        # 3.3.2. Review Fit
        ########################################
        
            ####################
            # 3.3.2.2. Model parameters
            ####################
            
            wald(Mb.asr)
            summary(Mb.asr)$var
            
            ####################
            # 3.3.2.3. Collate results
            ####################
            
            Mb_analysum <- data.frame("Model"= "Mb")
            Mb_analysum$Scale <- 'LN'
            Mb_analysum$Status <- 'RED'
            Mb_analysum$conv <- Mb.asr$converge
            
            Mb_analysum$prT <- round(wald(Mb.asr)[rownames(wald(Mb.asr))=='T',4],3) 
            
            Mb_analysum$vA <- round(summary(Mb.asr)$var[!grepl("vm(AID, giv_M):T",names(Mb.asr$vparameters),fixed=T) & grepl("vm(AID, giv_M)", names(Mb.asr$vparameters),fixed=T),"component"],3)
            
            Mb_analysum$vC <- round(summary(Mb.asr)$var[grepl("C!C",names(Mb.asr$vparameters),fixed=T),"component"],3)
            Mb_analysum$vCxT <- round(summary(Mb.asr)$var[grepl("C:T",names(Mb.asr$vparameters),fixed=T),"component"],3)
            
            Mb_analysum$vM4U <- round(summary( Mb.asr)$var[grepl("M4!U",names( Mb.asr$vparameters),fixed=T),"component"],3)
            Mb_analysum$vM5U <- round(summary( Mb.asr)$var[grepl("M5!U",names( Mb.asr$vparameters),fixed=T),"component"],3)
            
            Mb_analysum$logl <- Mb.asr$loglik
            Mb_analysum$df <- sum(summary(Mb.asr)$var$bound%in%c("P","U"))
            
            Mb_analysum$AIC <- 2*(Mb_analysum$df - Mb_analysum$logl)
        
        ########################################
        # 3.3.3. Test difference bn models
        ########################################
        
        Ma_v_Mb_tmd_in <- NULL
        Ma_v_Mb_tmd_in$analy_name <- 'Ma_v_Mb' 
        Ma_v_Mb_tmd_in$fullmod_name <- 'Ma' 
        Ma_v_Mb_tmd_in$redmod_name <- 'Mb'
        Ma_v_Mb_tmd_in$fullmod.asr <- Ma.asr
        Ma_v_Mb_tmd_in$redmod.asr <- Mb.asr
        Ma_v_Mb_tmd_in$test <- "boundary" #("boundary"/"other")
        
        Ma_v_Mb_tmd_out <- tmd.f(tmd_in = Ma_v_Mb_tmd_in)
        Ma_v_Mb_tmd_out
        
        Mb_analysum$Status <- 'BEST'
        
    ###########################################################
    # 3.4. Fit reduced model (Mb - C:T)
    ###########################################################
        
        ########################################
        # 3.4.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_M))
        
        Mc.asr <- update(Mb.asr,
                         random = ~ vm(AID, giv_M) + idv(C),)
        
        ########################################
        # 3.4.2. Review Fit
        ########################################
        
            ####################
            # 3.4.2.2. Model parameters
            ####################
            
            wald(Mc.asr)
            summary(Mc.asr)$var
            
            ####################
            # 3.4.2.3. Collate results
            ####################
            
            Mc_analysum <- data.frame("Model"= "Mc")
            Mc_analysum$Scale <- 'LN'
            Mc_analysum$Status <- 'RED'
            
            Mc_analysum$conv <- Mc.asr$converge
            
            Mc_analysum$prT <- round(wald(Mc.asr)[rownames(wald(Mc.asr))=='T',4],3) 
            
            Mc_analysum$vA <- round(summary(Mc.asr)$var[!grepl("vm(AID, giv_M):T",names(Mc.asr$vparameters),fixed=T) & grepl("vm(AID, giv_M)", names(Mc.asr$vparameters),fixed=T),"component"],3)
            
            Mc_analysum$vC <- round(summary(Mc.asr)$var[grepl("C!C",names(Mc.asr$vparameters),fixed=T),"component"],3)
            
            Mc_analysum$vM4U <- round(summary( Mc.asr)$var[grepl("M4!U",names( Mc.asr$vparameters),fixed=T),"component"],3)
            Mc_analysum$vM5U <- round(summary( Mc.asr)$var[grepl("M5!U",names( Mc.asr$vparameters),fixed=T),"component"],3)
            
            Mc_analysum$logl <- Mc.asr$loglik
            Mc_analysum$df <- sum(summary(Mc.asr)$var$bound%in%c("P","U"))
            
            Mc_analysum$AIC <- 2*(Mc_analysum$df - Mc_analysum$logl)
        
        ########################################
        # 3.4.3. Test difference bn models
        ########################################
        
        Mb_v_Mc_tmd_in <- NULL
        Mb_v_Mc_tmd_in$analy_name <- 'Mb_v_Mc' 
        Mb_v_Mc_tmd_in$fullmod_name <- 'Mb' 
        Mb_v_Mc_tmd_in$redmod_name <- 'Mc'
        Mb_v_Mc_tmd_in$fullmod.asr <- Mb.asr
        Mb_v_Mc_tmd_in$redmod.asr <- Mc.asr
        Mb_v_Mc_tmd_in$test <- "boundary" #("boundary"/"other")
        
        Mb_v_Mc_tmd_out <- tmd.f(tmd_in = Mb_v_Mc_tmd_in)
        Mb_v_Mc_tmd_out
        
#CHNOTE: Mc < Mb      
#CHNOTE: BEST MODEL = Mb        
        
    ###########################################################
    # 3.5. Fit reduced multivariate model (Md)
    ###########################################################
        
        ########################################
        # 3.5.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_M))
        
        Md.asr <- update(Mb.asr,
                         random = ~ vm(AID, giv_M) + corgh(T):id(C),)
        
        ########################################
        # 3.5.2. Review Fit
        ########################################
        
            ####################
            # 3.5.2.2. Model parameters
            ####################
            
            wald(Md.asr)
            summary(Md.asr)$var
            
            ####################
            # 3.5.2.3. Collate results
            ####################
            
            Md_analysum <- data.frame("Model"= "Md")
            Md_analysum$Status <- 'REDMV'
            Md_analysum$Scale <- 'LN'
            
            Md_analysum$conv <- Md.asr$converge
            
            Md_analysum$logl <- Md.asr$loglik
            Md_analysum$df <- sum(summary(Md.asr)$var$bound%in%c("P","U"))
            
            Md_analysum$AIC <- 2*(Md_analysum$df - Md_analysum$logl)
        
        ########################################
        # 3.5.3. Test difference bn models
        ########################################
        
        Md_v_Mb_tmd_in <- NULL
        Md_v_Mb_tmd_in$analy_name <- 'Md_v_Mb' 
        Md_v_Mb_tmd_in$fullmod_name <- 'Md' 
        Md_v_Mb_tmd_in$redmod_name <- 'Mb'
        Md_v_Mb_tmd_in$fullmod.asr <- Md.asr
        Md_v_Mb_tmd_in$redmod.asr <- Mb.asr
        Md_v_Mb_tmd_in$test <- "conventional" #("boundary"/"other")
        
        Md_v_Mb_tmd_out <- tmd.f(tmd_in = Md_v_Mb_tmd_in)
        Md_v_Mb_tmd_out
        
#CHNOTE: Mb = Md
#CHNOTE: BEST MODEL = Mb       
        
    ###########################################################
    # 3.6. Fit reduced (diagT) multivariate model (Me)
    ###########################################################
        
        ########################################
        # 3.6.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_M))
        
        Me.asr <- update(Mb.asr,
                         random = ~ vm(AID, giv_M) + diag(T):id(C),)
        
        ########################################
        # 3.6.2. Review Fit
        ########################################
        
            ####################
            # 3.6.2.2. Model parameters
            ####################
            
            wald(Me.asr)
            summary(Me.asr)$var
            
            ####################
            # 3.6.2.3. Collate results
            ####################
            
            Me_analysum <- data.frame("Model"= "Me")
            Me_analysum$Status <- 'MVRED'
            Me_analysum$Scale <- 'LN'
            Me_analysum$conv <- Me.asr$converge
           
            Me_analysum$logl <- Me.asr$loglik
            Me_analysum$df <- sum(summary(Me.asr)$var$bound%in%c("P","U"))
            
            Me_analysum$AIC <- 2*(Me_analysum$df - Me_analysum$logl)
        
        ####################
        # 3.6.3. Test difference bn models
        ####################
        
        Md_v_Me_tmd_in <- NULL
        Md_v_Me_tmd_in$analy_name <- 'Md_v_Me' 
        Md_v_Me_tmd_in$fullmod_name <- 'Md' 
        Md_v_Me_tmd_in$redmod_name <- 'Me'
        Md_v_Me_tmd_in$fullmod.asr <- Md.asr
        Md_v_Me_tmd_in$redmod.asr <- Me.asr
        Md_v_Me_tmd_in$test <- "conventional" #("boundary"/"other")
        
        Md_v_Me_tmd_out <- tmd.f(tmd_in = Md_v_Me_tmd_in)
        Md_v_Me_tmd_out
        
#CHNOTE: Md = Me
        Me_analysum$Status <- 'BEST'       
        
    ###########################################################
    # 3.6. Fit reduced multivariate model (Mf)
    ###########################################################
        
        ########################################
        # 3.6.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_M))
        
        Mf.asr <- update(Me.asr,
                         random = ~ vm(AID, giv_M) + corv(T):id(C),)
        
        summary(Mf.asr)$var
        
  
    ###########################################################
    # 3.99 Combine Model summary
    ###########################################################
        
    M_analysum <- smartbind(Ma_analysum,
                             Mb_analysum,
                             Mc_analysum,
                             Md_analysum,
                             Me_analysum)    
        
    write.csv(M_analysum)
        
###############################################################################
# 4. P
###############################################################################
    
    ###########################################################
    # 4.1. Prepare data
    ###########################################################
    
        ########################################
        # 4.1.1. Select pdata
        ########################################
        
        sort(unique(pdata_in$L))
    
        pdata_P <- pdata_in[pdata_in$T%in%c("P4", "P5"),]
        opdata_P <- pdata_P[order(pdata_P$U,pdata_P$Y,pdata_P$H),]
        
        ##### Pake list of accession in pdata
        ID_P <- unique(opdata_P$ID)
        
        ########################################
        # 4.1.2. Select sdata and construct giv
        ########################################
        
        ##### Select sdata_B and construct giv
        sdata_P <- sdata[ID_P,]
        
        grm_P <-A.mat(sdata_P,min.MAF = 0.05,max.missing = NULL,impute.method="mean",tol=0.02,
                      n.core=1,shrink=T,return.imputed=FALSE, ploidy=2)
        
        giv_P <- solve(grm_P)
        hist(diag(giv_P))
        attr(giv_P, "INVERSE") = TRUE
        
        ########################################
        # 4.1.3. Set initial factors
        ########################################
        
        adata_P <- opdata_P
        ##### Block
        adata_P$B <- factor(adata_P$B)
        
        ##### Row
        adata_P$X <- factor(adata_P$X)

        ##### Plot
        adata_P$P <- factor(adata_P$P)
        
        ##### Genetic factors
        adata_P$AID <- factor(adata_P$ID, levels = rownames(giv_P))
        
        ##### U = L,P,
        adata_P$U <- factor(adata_P$U)
        
        ##### Assessment Year
        adata_P$Y <- factor(adata_P$Y)
        
        ##### H = harvest
        adata_P$H <- factor(adata_P$H)
        
        ##### L = location
        adata_P$L <- factor(adata_P$L)
        
        ##### T = TRIAL
        adata_P$T <- factor(adata_P$T)
        
         ##### C = clonal
        adata_P$C <- factor(adata_P$AID)
        
        ############################################################
        # 4.1.4. Review data
        ############################################################
        
        str(adata_P)
        
        tapply(adata_P$SSC,adata_P$Y,FUN = mean,na.rm = T)
        tapply(adata_P$lnSSC,adata_P$Y,FUN = mean,na.rm = T)
        tapply(adata_P$SSC,adata_P$Y,FUN = var,na.rm = T)
        tapply(adata_P$lnSSC,adata_P$Y,FUN = var,na.rm = T)
    
    ###########################################################
    # 4.2. Fit Full model (Pa)
    ###########################################################
    
        ########################################
        # 4.2.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_P))
        
        Pa.asr = asreml(lnSSC ~ -1 + T + at(T):B + at(T):B:D, 
                            random = ~ at(T,c('P4')):B:D:X +vm(AID, giv_P) + vm(AID, giv_P):idv(T) + idv(C) + id(C):idv(T),
                            residual = ~ dsum(~idv(U)| T),
                          na.action = na.method(y = "include", x = "include"),
                          data = adata_P,
                          ai.loadings = TRUE,
                          workspace = "16gb")
 
        ########################################
        # 4.2.2. Review Fit
        ########################################
        
            ####################
            # 4.2.2.2. Model parameters
            ####################
            
            wald(Pa.asr)
            summary(Pa.asr)$var
            
            ####################
            # 4.2.2.3. Collate results
            ####################
            
            Pa_analysum <- data.frame("Model"= "Pa")
            Pa_analysum$Scale <- 'LN'
            Pa_analysum$Status <- 'FULL'
            
            Pa_analysum$conv = Pa.asr$converge
            
            Pa_analysum$Trial <- round(wald(Pa.asr)[1,4],3)
            
            Pa_analysum$prBP4T <- round(wald(Pa.asr)[rownames(wald(Pa.asr))=='at(T, P4):B',4],3)
            Pa_analysum$prBP5T <- round(wald(Pa.asr)[rownames(wald(Pa.asr))=='at(T, P5):B',4],3)
            
            Pa_analysum$prBDP4T <- round(wald(Pa.asr)[rownames(wald(Pa.asr))=='at(T, P4):B:D',4],3)
            Pa_analysum$prBDP5T <- round(wald(Pa.asr)[rownames(wald(Pa.asr))=='at(T, P5):B:D',4],3)
            
            Pa_analysum$vTP4XDB <- round(summary(Pa.asr)$var[grepl("at(T, P4):X:D:B",names(Pa.asr$vparameters),fixed=T),"component"],3)
            Pa_analysum$vTP5XDB <- round(summary(Pa.asr)$var[grepl("at(T, P5):X:D:B",names(Pa.asr$vparameters),fixed=T),"component"],3)
           
            Pa_analysum$vA <- round(summary(Pa.asr)$var[!grepl("vm(AID, giv_P):T",names(Pa.asr$vparameters),fixed=T) & 
                                                            grepl("vm(AID, giv_P)", names(Pa.asr$vparameters),fixed=T),"component"],3)
            Pa_analysum$vAxT <- round(summary(Pa.asr)$var[grepl("vm(AID, giv_P):T",names(Pa.asr$vparameters),fixed=T),"component"],3)
            
            Pa_analysum$vC <- round(summary(Pa.asr)$var[grepl("C!C",names(Pa.asr$vparameters),fixed=T),"component"],3)
            Pa_analysum$vCxT <- round(summary(Pa.asr)$var[grepl("C:T!T",names(Pa.asr$vparameters),fixed=T),"component"],3)
          
            Pa_analysum$vP4U <- round(summary( Pa.asr)$var[grepl("P4!U",names( Pa.asr$vparameters),fixed=T),"component"],3)
            Pa_analysum$vP5U <- round(summary( Pa.asr)$var[grepl("P5!U",names( Pa.asr$vparameters),fixed=T),"component"],3)
            
            Pa_analysum$logl <- Pa.asr$loglik
            Pa_analysum$df <- sum(summary(Pa.asr)$var$bound%in%c("P","U"))
            
            Pa_analysum$AIC <- 2*(Pa_analysum$df - Pa_analysum$logl)
            
    ###########################################################
    # 4.3. Fit reduced (Pa - C:T)
    ###########################################################
        
        ########################################
        # 4.3.1. Fit model
        ########################################

        asreml.options(dense = ~ vm(AID, giv_P))
            
        Pb.asr = update(Pa.asr,
                        random = ~ at(T,c('P4')):B:D:X +vm(AID, giv_P) + vm(AID, giv_P):idv(T) + idv(C),)
        
        ########################################
        # 4.3.2. Review Fit
        ########################################
        
            ####################
            # 4.3.2.2. Model parameters
            ####################
            
            wald(Pb.asr)
            summary(Pb.asr)$var
        
            ####################
            # 4.3.2.3. Collate results
            ####################
            
            Pb_analysum <- data.frame("Model"= "Pb")
            Pb_analysum$Scale <- 'LN'
            Pb_analysum$Status <- 'RED'
            
            Pb_analysum$conv = Pb.asr$converge
            
            Pb_analysum$Trial <- round(wald(Pb.asr)[1,4],3)
            
            Pb_analysum$prBP4T <- round(wald(Pb.asr)[rownames(wald(Pb.asr))=='at(T, P4):B',4],3)
            Pb_analysum$prBP5T <- round(wald(Pb.asr)[rownames(wald(Pb.asr))=='at(T, P5):B',4],3)
            
            Pb_analysum$prBDP4T <- round(wald(Pb.asr)[rownames(wald(Pb.asr))=='at(T, P4):B:D',4],3)
            Pb_analysum$prBDP5T <- round(wald(Pb.asr)[rownames(wald(Pb.asr))=='at(T, P5):B:D',4],3)
            
            Pb_analysum$vTP4XDB <- round(summary(Pb.asr)$var[grepl("at(T, P4):X:D:B",names(Pb.asr$vparameters),fixed=T),"component"],3)
            Pb_analysum$vTP5XDB <- round(summary(Pb.asr)$var[grepl("at(T, P5):X:D:B",names(Pb.asr$vparameters),fixed=T),"component"],3)
            
            Pb_analysum$vA <- round(summary(Pb.asr)$var[!grepl("vm(AID, giv_P):T",names(Pb.asr$vparameters),fixed=T) & 
                                                            grepl("vm(AID, giv_P)", names(Pb.asr$vparameters),fixed=T),"component"],3)
            Pb_analysum$vAxT <- round(summary(Pb.asr)$var[grepl("vm(AID, giv_P):T",names(Pb.asr$vparameters),fixed=T),"component"],3)
            
            Pb_analysum$vC <- round(summary(Pb.asr)$var[grepl("C!C",names(Pb.asr$vparameters),fixed=T),"component"],3)
            Pb_analysum$vP4U <- round(summary( Pb.asr)$var[grepl("P4!U",names( Pb.asr$vparameters),fixed=T),"component"],3)
            Pb_analysum$vP5U <- round(summary( Pb.asr)$var[grepl("P5!U",names( Pb.asr$vparameters),fixed=T),"component"],3)
            
            Pb_analysum$logl <- Pb.asr$loglik
            Pb_analysum$df <- sum(summary(Pb.asr)$var$bound%in%c("P","U"))
            
            Pb_analysum$AIC <- 2*(Pb_analysum$df - Pb_analysum$logl)
            
        ########################################
        # 3.5.3. Test difference bn models
        ########################################
        
        Pa_v_Pb_tmd_in <- NULL
        Pa_v_Pb_tmd_in$analy_name <- 'Pa_v_Pb' 
        Pa_v_Pb_tmd_in$fullmod_name <- 'Pa' 
        Pa_v_Pb_tmd_in$redmod_name <- 'Pb'
        Pa_v_Pb_tmd_in$fullmod.asr <- Pa.asr
        Pa_v_Pb_tmd_in$redmod.asr <- Pb.asr
        Pa_v_Pb_tmd_in$test <- "boundary" #("boundary"/"other")
        
        Pa_v_Pb_tmd_out <- tmd.f(tmd_in = Pa_v_Pb_tmd_in)
        Pa_v_Pb_tmd_out
        
#CHNOTE: Pb = Pa
#CHNOTE: BEST MODEL = Pb
        
        Pb_analysum$Status <- 'BEST'
        
    ###########################################################
    # 4.4. Fit reduced (Pc - A:T) model 
    ###########################################################
        
        ########################################
        # 4.4.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_P))
        
        Pc.asr = update(Pb.asr,
                        random = ~ at(T,c('P4')):B:D:X +vm(AID, giv_P) + idv(C),)
        
        ########################################
        # 4.4.2. Review Fit
        ########################################
        
            ####################
            # 4.4.2.2. Model parameters
            ####################
            
            wald(Pc.asr)
            summary(Pc.asr)$var
            
            ####################
            # 4.4.2.3. Collate results
            ####################
            
            Pc_analysum <- data.frame("Model"= "Pc")
            Pc_analysum$Status <- 'RED'
            Pc_analysum$Scale <- 'LN'
            Pc_analysum$conv = Pc.asr$converge
            
            Pc_analysum$Trial <- round(wald(Pc.asr)[1,4],3)
            
            Pc_analysum$prBP4T <- round(wald(Pc.asr)[rownames(wald(Pc.asr))=='at(T, P4):B',4],3)
            Pc_analysum$prBP5T <- round(wald(Pc.asr)[rownames(wald(Pc.asr))=='at(T, P5):B',4],3)
            
            Pc_analysum$prBDP4T <- round(wald(Pc.asr)[rownames(wald(Pc.asr))=='at(T, P4):B:D',4],3)
            Pc_analysum$prBDP5T <- round(wald(Pc.asr)[rownames(wald(Pc.asr))=='at(T, P5):B:D',4],3)
            
            Pc_analysum$vTP4XDB <- round(summary(Pc.asr)$var[grepl("at(T, P4):X:D:B",names(Pc.asr$vparameters),fixed=T),"component"],3)
            Pc_analysum$vTP5XDB <- round(summary(Pc.asr)$var[grepl("at(T, P5):X:D:B",names(Pc.asr$vparameters),fixed=T),"component"],3)
            
            Pc_analysum$vA <- round(summary(Pc.asr)$var[!grepl("vm(AID, giv_P):T",names(Pc.asr$vparameters),fixed=T) & 
                                                            grepl("vm(AID, giv_P)", names(Pc.asr$vparameters),fixed=T),"component"],3)
            
            Pc_analysum$vC <- round(summary(Pc.asr)$var[grepl("C!C",names(Pc.asr$vparameters),fixed=T),"component"],3)
            
            Pc_analysum$vP4U <- round(summary( Pc.asr)$var[grepl("P4!U",names( Pc.asr$vparameters),fixed=T),"component"],3)
            Pc_analysum$vP5U <- round(summary( Pc.asr)$var[grepl("P5!U",names( Pc.asr$vparameters),fixed=T),"component"],3)
            
            Pc_analysum$logl <- Pc.asr$loglik
            Pc_analysum$df <- sum(summary(Pc.asr)$var$bound%in%c("P","U"))
            
            Pc_analysum$AIC <- 2*(Pc_analysum$df - Pc_analysum$logl)
        
        ########################################
        # 4.4.3. Test difference bn models
        ########################################
        
        Pb_v_Pc_tmd_in <- NULL
        Pb_v_Pc_tmd_in$analy_name <- 'Pb_v_Pc' 
        Pb_v_Pc_tmd_in$fullmod_name <- 'Pb' 
        Pb_v_Pc_tmd_in$redmod_name <- 'Pc'
        Pb_v_Pc_tmd_in$fullmod.asr <- Pb.asr
        Pb_v_Pc_tmd_in$redmod.asr <- Pc.asr
        Pb_v_Pc_tmd_in$test <- "boundary" #("boundary"/"other")
        
        Pb_v_Pc_tmd_out <- tmd.f(tmd_in = Pb_v_Pc_tmd_in)
        Pb_v_Pc_tmd_out
        
#CHNOTE: Pc = Pb
#CHNOTE: BEST MODEL = Pc            
            
        Pb_analysum$Status <- 'RED'
        Pc_analysum$Status <- 'BEST'
        
    ###########################################################
    # 4.5. Fit reduced (Pc - X) model 
    ###########################################################
        
        ########################################
        # 4.5.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_P))
        
        Pd.asr = update(Pb.asr,
                        random = ~ vm(AID, giv_P) + idv(C),)
        
        ########################################
        # 4.5.2. Review Fit
        ########################################
        
            ####################
            # 4.5.2.2. Model parameters
            ####################
            
            wald(Pd.asr)
            summary(Pd.asr)$var
            
            ####################
            # 4.5.2.3. Collate results
            ####################
            
            Pd_analysum <- data.frame("Model"= "Pd")
            Pd_analysum$Status <- 'RED'
            Pd_analysum$Scale <- 'LN'
            Pd_analysum$conv = Pd.asr$converge
            
            Pd_analysum$Trial <- round(wald(Pd.asr)[1,4],3)
            
            Pd_analysum$prBP4T <- round(wald(Pd.asr)[rownames(wald(Pd.asr))=='at(T, P4):B',4],3)
            Pd_analysum$prBP5T <- round(wald(Pd.asr)[rownames(wald(Pd.asr))=='at(T, P5):B',4],3)
            
            Pd_analysum$prBDP4T <- round(wald(Pd.asr)[rownames(wald(Pd.asr))=='at(T, P4):B:D',4],3)
            Pd_analysum$prBDP5T <- round(wald(Pd.asr)[rownames(wald(Pd.asr))=='at(T, P5):B:D',4],3)
            
            Pd_analysum$vTP4XDB <- round(summary(Pd.asr)$var[grepl("at(T, P4):X:D:B",names(Pd.asr$vparameters),fixed=T),"component"],3)
            Pd_analysum$vTP5XDB <- round(summary(Pd.asr)$var[grepl("at(T, P5):X:D:B",names(Pd.asr$vparameters),fixed=T),"component"],3)
            
            Pd_analysum$vA <- round(summary(Pd.asr)$var[!grepl("vm(AID, giv_P):T",names(Pd.asr$vparameters),fixed=T) & 
                                                            grepl("vm(AID, giv_P)", names(Pd.asr$vparameters),fixed=T),"component"],3)
            
            Pd_analysum$vC <- round(summary(Pd.asr)$var[grepl("C!C",names(Pd.asr$vparameters),fixed=T),"component"],3)
            
            Pd_analysum$vP4U <- round(summary( Pd.asr)$var[grepl("P4!U",names( Pd.asr$vparameters),fixed=T),"component"],3)
            Pd_analysum$vP5U <- round(summary( Pd.asr)$var[grepl("P5!U",names( Pd.asr$vparameters),fixed=T),"component"],3)
            
            Pd_analysum$logl <- Pd.asr$loglik
            Pd_analysum$df <- sum(summary(Pd.asr)$var$bound%in%c("P","U"))
            
            Pd_analysum$AIC <- 2*(Pd_analysum$df - Pd_analysum$logl)
            
        ########################################
        # 4.5.3. Test difference bn models
        ########################################
        
        Pc_v_Pd_tmd_in <- NULL
        Pc_v_Pd_tmd_in$analy_name <- 'Pc_v_Pd' 
        Pc_v_Pd_tmd_in$fullmod_name <- 'Pc' 
        Pc_v_Pd_tmd_in$redmod_name <- 'Pd'
        Pc_v_Pd_tmd_in$fullmod.asr <- Pc.asr
        Pc_v_Pd_tmd_in$redmod.asr <- Pd.asr
        Pc_v_Pd_tmd_in$test <- "boundary" #("boundary"/"other")
        
        Pc_v_Pd_tmd_out <- tmd.f(tmd_in = Pc_v_Pd_tmd_in)
        Pc_v_Pd_tmd_out
        
        #CHNOTE: Pd < Pc
        #CHNOTE: BEST MODEL = Pd              
        
        
    ###########################################################
    # 4.6. Fit reduced (Pc - C) model 
    ###########################################################
        
        ########################################
        # 4.6.1. Fit model
        ########################################
        
        asreml.options(dense = ~ vm(AID, giv_P))
        
        Pe.asr = update(Pb.asr,
                        random = ~ at(T,c('P4')):B:D:X +vm(AID, giv_P),)
        
        ########################################
        # 4.6.2. Review Fit
        ########################################
        
            ####################
            # 4.6.2.2. Model parameters
            ####################
            
            wald(Pe.asr)
            summary(Pe.asr)$var
            
            ####################
            # 4.6.2.3. Collate results
            ####################
            
            Pe_analysum <- data.frame("Model"= "Pe")
            Pe_analysum$Status <- 'RED'
            Pe_analysum$Scale <- 'LN'
            Pe_analysum$conv = Pe.asr$converge
            
            Pe_analysum$Trial <- round(wald(Pe.asr)[1,4],3)
            
            Pe_analysum$prBP4T <- round(wald(Pe.asr)[rownames(wald(Pe.asr))=='at(T, P4):B',4],3)
            Pe_analysum$prBP5T <- round(wald(Pe.asr)[rownames(wald(Pe.asr))=='at(T, P5):B',4],3)
            
            Pe_analysum$prBDP4T <- round(wald(Pe.asr)[rownames(wald(Pe.asr))=='at(T, P4):B:D',4],3)
            Pe_analysum$prBDP5T <- round(wald(Pe.asr)[rownames(wald(Pe.asr))=='at(T, P5):B:D',4],3)
            
            Pe_analysum$vTP4XDB <- round(summary(Pe.asr)$var[grepl("at(T, P4):B:D:X",names(Pe.asr$vparameters),fixed=T),"component"],3)
            #Pe_analysum$vTP5XDB <- round(summary(Pe.asr)$var[grepl("at(T, P5):X:D:B",names(Pe.asr$vparameters),fixed=T),"component"],3)
            
            Pe_analysum$vA <- round(summary(Pe.asr)$var[grepl("vm(AID, giv_P)",names(Pe.asr$vparameters),fixed=T),"component"],3) 
                                                            # & 
                                                            #grepl("vm(AID, giv_P)", names(Pe.asr$vparameters),fixed=T)
            
           # Pe_analysum$vC <- round(summary(Pe.asr)$var[grepl("C!C",names(Pe.asr$vparameters),fixed=T),"component"],3)
            
            Pe_analysum$vP4U <- round(summary( Pe.asr)$var[grepl("P4!U",names( Pe.asr$vparameters),fixed=T),"component"],3)
            Pe_analysum$vP5U <- round(summary( Pe.asr)$var[grepl("P5!U",names( Pe.asr$vparameters),fixed=T),"component"],3)
            
            Pe_analysum$logl <- Pe.asr$loglik
            Pe_analysum$df <- sum(summary(Pe.asr)$var$bound%in%c("P","U"))
            
            Pe_analysum$AIC <- 2*(Pe_analysum$df - Pe_analysum$logl)
            
            ########################################
            # 4.6.3. Test difference bn models
            ########################################
            
            Pc_v_Pe_tmd_in <- NULL
            Pc_v_Pe_tmd_in$analy_name <- 'Pc_v_Pe' 
            Pc_v_Pe_tmd_in$fullmod_name <- 'Pc' 
            Pc_v_Pe_tmd_in$redmod_name <- 'Pe'
            Pc_v_Pe_tmd_in$fullmod.asr <- Pc.asr
            Pc_v_Pe_tmd_in$redmod.asr <- Pe.asr
            Pc_v_Pe_tmd_in$test <- "boundary" #("boundary"/"other")
            
            Pc_v_Pe_tmd_out <- tmd.f(tmd_in = Pc_v_Pe_tmd_in)
            Pc_v_Pe_tmd_out
            
            #CHNOTE: Pe < Pc
            #CHNOTE: BEST MODEL = Pe                 
        
    ############################################################
    # 4.6. Combine analysum
    ############################################################
        
    P_analysum <- smartbind(Pa_analysum,
                            Pb_analysum,
                            Pc_analysum,
                            Pd_analysum,
                            Pe_analysum)