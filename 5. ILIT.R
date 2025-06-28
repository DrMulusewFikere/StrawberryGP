#####################################################################################################
#                                                                                                   #
#                       R-Script for Genomic Prediction in Strawberry                               #
#                       Date:  2021-11-11                                                           #                                                                        #
#########################################################################################################
#                                                                                                       #
# This routine curates SNP data and Phenotypic data using the R Package ASRgenomics. Individuals that   #
#are duplicated in the SNP data are removed, Duplicated individuals in the phenotypic data are assigned #
# the same name to avoid maintain the size of the training population                                   #
# Markers with MAF < 0.05 are removed                                                                   #
#                                                                                                       #
########################   LIST OF ABBREVIATIONS                            #############################                       
#                                                                                                   #                                                                                    
#                                                                                                   #
#                                                                                                   #
#                                                                                                   #
#                                                                                                   #
#                                                                                                   #
#   pdata  - Phenotypic data                                                                        #
#   sdata  - SNP genotypic data                                                                     #
#   GRM    - Genomic relationship Matrix                                                            #
#                                                                                                   #                                                          #
#                                                                                                   #
#                                                                                                   #
#####################################################################################################
# Individual trial analyses
###############################################################################

###############################################################################
# 1. Set up R ----
###############################################################################

############################################################
# 1.1. Load packages ----
############################################################

#library(asreml, lib.loc = "C:/Users/uqchardn/Dropbox/Documents/R/R-3.6.0/Asreml_versions/Asreml4")
library(asreml)
#library(sommer)
library(asremlPlus)
library(gtools)
library(ASRgenomics)
library(rrBLUP)
library(AGHmatrix)
library(dplyr)

#### Set the working directory for the sdata and pdata ----
setwd("~/Dropbox/Post_201203/Analysis_post211021/4. Fit models/4.0. Xtra data curation")
#load("Strawberry.RData")
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
  # tmd_in$test ("boundary"/"other")
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
# 2. Prepare data ----
###############################################################################

############################################################
# 2.1. input data ----
############################################################

##### pdata ## NM ---- 
#pdata_in_temp <- read.csv("../../../0. Data/data_Ubal_post210120.csv", header = TRUE, stringsAsFactors = FALSE)
  pdata_in_temp <- read.csv("data_Ubal_QC_post211111.csv", header = TRUE, stringsAsFactors = FALSE)

##### remove very low values 
pdata_in_temp[pdata_in_temp$SSC < 4 & !is.na(pdata_in_temp$SSC),] #checks the number of rows with SSC < 4
pdata_in <- pdata_in_temp
pdata_in$SSC[pdata_in$SSC < 4] <- NA # set SSC values < 4 as missing (NA)

pdata_in$lnSSC <- log(pdata_in$SSC)*100

##### sdata # NM ----
#sdata_in <- read.csv("../../../0. Data/sdata_imputed_2064.csv",header = T,stringsAsFactors = F)
sdata_in <- read.csv("sdata_imputed_2064_qc_NM211111.csv",header = T,stringsAsFactors = F)

sdata <- sdata_in[,2:ncol(sdata_in)]
rownames(sdata) <- sdata_in[,1]

###############################################################################
# 3. B1 ----
###############################################################################

###########################################################
# 3.1. Prepare data 
###########################################################

########################################
# 3.1.1. Select pdata
########################################

sort(unique(pdata_in$T)) #gives the unique Trials in pdata_in

pdata_B1 <- pdata_in[pdata_in$T%in%c("B1"),]
opdata_B1 <- pdata_B1[order(pdata_B1$U,pdata_B1$Y,pdata_B1$H),] #orders the pdata_B1 by column U followed by Y & H

##### Make list of accession in pdata
ID_B1 <- unique(opdata_B1$ID)



########################################
# 3.1.2. Select sdata and construct giv ----
########################################

##### Create sdata_B1
sdata_B1 <- sdata[rownames(sdata)%in%ID_B1,]
sdata_B1 <- as.matrix(sdata_B1)

#create pheno_data containing individuals that are present in sdata_B1
pheno_data <- opdata_B1[opdata_B1$ID%in%rownames(sdata_B1),]
###FROM FM Make GRM positive definite

   
 
 
 # Construct the additive GRM using the method of VanRaden 
 G_matrix <- G.matrix(M = sdata_B1, method = "VanRaden", na.string = NA)$G
 
 #Checking if the G_matrix is positive definite 
 is.positive.definite(G_matrix) 
 
 ### tuneup G by blending and blending since its ill conditioned
 G_tuneup <- G.tuneup(G = G_matrix, bend = TRUE, message = T)$Gb
 
 #Checking if the G_tuneup is positive definite 
 is.positive.definite(G_tuneup) 
 
 ## Invert the GRM (The tuned-up G_tuneup)
 G_check <- G.inverse(G= G_tuneup, sparseform = FALSE)$Ginv
 
## Making sure that the inverted GRM is inverted
 attr(G_check, "INVERSE") = TRUE  

#checking if all the individuals in the phenotypic data file were genotyped and if all
     #individuals in the GRM are in the phenotypic datafile
 pheno.G <- match.kinship2pheno(K = G_check, pheno.data = pheno_data,
                                indiv = "ID", clean = FALSE, mism = TRUE)

#########################################################################################
giv_B1 <- G_check
########################################
# 3.1.3. Set initial factors
########################################
adata_B1 <- pheno_data #### using the phenotypic data from which we removed duplicates
#adata_B1 <- opdata_B1
##### Block
adata_B1$B <- factor(adata_B1$B)

##### Row
adata_B1$X <- factor(adata_B1$X)

##### Plot
adata_B1$P <- factor(adata_B1$P)

##### Genetic factors
adata_B1$AID <- factor(adata_B1$ID, levels = rownames(giv_B1))

##### U = L,P,
adata_B1$U <- factor(adata_B1$U)

##### Assessment Year
adata_B1$Y <- factor(adata_B1$Y)

##### H = harvest
adata_B1$H <- factor(adata_B1$H)

############################################################
# 3.1.5. Review data
############################################################

str(adata_B1)

tapply(adata_B1$SSC,adata_B1$Y,FUN = mean,na.rm = T)
tapply(adata_B1$SSC,adata_B1$Y,FUN = var,na.rm = T)

##### Number of observations by Year
adata_B1_nna <- adata_B1[!is.na(adata_B1$SSC),]

table(adata_B1_nna$H)

###########################################################
# 3.2. Fit Full model (B1a) - untransformed data
###########################################################

########################################

# 3.2.1. Fit model
########################################
asreml.options(dense = ~ vm(AID, giv_B1))

B1a.asr = asreml(SSC ~ 1 + Y, 
                 random = ~ vm(AID, giv_B1)+ vm(AID, giv_B1):idv(Y) + idv(U),
                 na.action = na.method(y = "include", x = "include"),
                 data = adata_B1,
                 ai.loadings = TRUE,
                 workspace = "10gb")
###########################
update(B1a.asr)
plot(B1a.asr)



## Obtaining Predictions - BLUP - 
BLUP<-summary(B1a.asr,coef=TRUE)$coef.random
wald(B1a.asr)
summary(B1a.asr)$var

#cor(pheno_data$SSC,BLUP,method='pearson',use="complete.obs")



########################################
# 3.2.2. Review Fit
########################################

options(scipen = 999)

####################
# 3.2.2.1. Normal distribution of residuals
####################

plot(B1a.asr)

B1a.sR <- as.numeric(sqrt(var(B1a.asr$residuals,na.rm=T)))

adata_B1_stdres <- cbind(adata_B1,(B1a.asr$residuals/B1a.sR))

B1a_shapiro <- shapiro.test(adata_B1_stdres$e)

#CHNOTE: OK            

##### Identify std res > 3.5
B1a_nstdres <- adata_B1_stdres[abs(adata_B1_stdres$e)>3.5 & !is.na(adata_B1_stdres$e),]
#CHNOTE: None

####################
# 3.2.2.2. Model parameters
####################

wald(B1a.asr)
summary(B1a.asr)$var

####################
# 3.2.2.3. Collate results
####################

B1a_analysum <- data.frame("Model"= "B1a")
B1a_analysum$Scale <- "ID"
B1a_analysum$Status = "FULL"

B1a_analysum$conv <- B1a.asr$converge

B1a_analysum$prY <- round(wald(B1a.asr)[rownames(wald(B1a.asr))=='Y',4],3)
B1a_analysum$vA <- round(summary(B1a.asr)$var[!grepl("vm(AID, giv_B1):Y",names(B1a.asr$vparameters),fixed=T) & grepl("vm(AID, giv_B1)", names(B1a.asr$vparameters),fixed=T),"component"],3)
B1a_analysum$vAxY <- round(summary(B1a.asr)$var[grepl("vm(AID, giv_B1):Y",names(B1a.asr$vparameters),fixed=T),"component"],3)
B1a_analysum$vU <- round(summary( B1a.asr)$var[grepl("!U",names( B1a.asr$vparameters),fixed=T),"component"],3)
B1a_analysum$vR <- round(summary( B1a.asr)$var[grepl("!R",names( B1a.asr$vparameters),fixed=T),"component"],3)

B1a_analysum$logl <- B1a.asr$loglik
B1a_analysum$df <- sum(summary(B1a.asr)$var$bound%in%c("P","U"))
B1a_analysum$AIC <- 2*(B1a_analysum$df - B1a_analysum$logl)

B1a_analysum$W <- B1a_shapiro$statistic
B1a_analysum$pW <- B1a_shapiro$p.value
B1a_analysum$nstdres <- nrow(B1a_nstdres)

###########################################################
# 3.3. Fit Full model (B1b) -  log transofrmed
###########################################################

########################################
# 3.3.1. Prepare data
########################################

adata_B1$lnSSC <- log(adata_B1$SSC)*100

########################################
# 3.3.1. Fit model
########################################

asreml.options(dense = ~ vm(AID, giv_B1))

B1b.asr = asreml(lnSSC ~ 1 + Y, 
                 random = ~ vm(AID, giv_B1) + vm(AID, giv_B1):idv(Y) + idv(U),
                 na.action = na.method(y = "include", x = "include"),
                 data = adata_B1,
                 ai.loadings = TRUE,
                 workspace = "12gb")
update(B1b.asr)
########################################
# 3.3.2. Review Fit
########################################

options(scipen = 999)
#plot(B1b.asr)
####################
# 3.3.2.1. Normal distribution of residuals
####################

plot(B1b.asr)

B1b.sR <- as.numeric(sqrt(var(B1b.asr$residuals,na.rm=T)))

adata_B1_stdres <- cbind(adata_B1,(B1b.asr$residuals/B1b.sR))

B1b_shapiro <- shapiro.test(adata_B1_stdres$e)

#CHNOTE: OK            

##### Identify std res > 3.5
B1b_nstdres <- adata_B1_stdres[abs(adata_B1_stdres$e)>3.5 & !is.na(adata_B1_stdres$e),]
#CHNOTE: None

####################
# 3.3.2.2. Model parameters
####################

wald(B1b.asr)
summary(B1b.asr)$var
update(B1b.asr)
####################
# 3.3.2.3. Collate results
####################

B1b_analysum <- data.frame("Model"= "B1b")
B1b_analysum$Scale <- "LN"
B1b_analysum$Status = "FULL"

B1b_analysum$conv <- B1b.asr$converge

B1b_analysum$prY <- round(wald(B1b.asr)[rownames(wald(B1b.asr))=='Y',4],3)
B1b_analysum$vA <- round(summary(B1b.asr)$var[!grepl("vm(AID, giv_B1):Y",names(B1b.asr$vparameters),fixed=T) & grepl("vm(AID, giv_B1)", names(B1b.asr$vparameters),fixed=T),"component"],3)
B1b_analysum$vAxY <- round(summary(B1b.asr)$var[grepl("vm(AID, giv_B1):Y",names(B1b.asr$vparameters),fixed=T),"component"],3)
B1b_analysum$vU <- round(summary( B1b.asr)$var[grepl("!U",names( B1b.asr$vparameters),fixed=T),"component"],3)
B1b_analysum$vR <- round(summary( B1b.asr)$var[grepl("!R",names( B1b.asr$vparameters),fixed=T),"component"],3)

B1b_analysum$logl <- B1b.asr$loglik
B1b_analysum$df <- sum(summary(B1b.asr)$var$bound%in%c("P","U"))
B1b_analysum$AIC <- 2*(B1b_analysum$df - B1b_analysum$logl)

B1b_analysum$W <- B1b_shapiro$statistic
B1b_analysum$pW <- B1b_shapiro$p.value
B1b_analysum$nstdres <- nrow(B1b_nstdres)

###########################################################
# 3.4. Fit reduced model (B1c) -  log transformed
###########################################################

########################################
# 3.4.2. Fit model
########################################

asreml.options(dense = ~ vm(AID, giv_B1))

B1c.asr = asreml(lnSSC ~ 1 + Y, 
                 random = ~ vm(AID, giv_B1)
                 + idv(U),
                 na.action = na.method(y = "include", x = "include"),
                 data = adata_B1,
                 ai.loadings = TRUE,
                 workspace = "8gb")

########################################
# 3.4.3. Review Fit
########################################

options(scipen = 999)
update(B1c.asr)
####################
# 3.4.3.2. Model parameters
####################

wald(B1c.asr)
summary(B1c.asr)$var
plot(B1c.asr)
####################
# 3.4.3.3. Collate results
####################

B1c_analysum <- data.frame("Model"= "B1c")
B1c_analysum$Scale <- "LN"
B1c_analysum$Status = "RED"
B1c_analysum$conv <- B1c.asr$converge

B1c_analysum$prY <- round(wald(B1c.asr)[rownames(wald(B1c.asr))=='Y',4],3)
B1c_analysum$vA <- round(summary(B1c.asr)$var[!grepl("vm(AID, giv_B1):Y",names(B1c.asr$vparameters),fixed=T) & grepl("vm(AID, giv_B1)", names(B1c.asr$vparameters),fixed=T),"component"],3)
B1c_analysum$vAxY <- round(summary(B1c.asr)$var[grepl("vm(AID, giv_B1):Y",names(B1c.asr$vparameters),fixed=T),"component"],3)
B1c_analysum$vU <- round(summary( B1c.asr)$var[grepl("!U",names( B1c.asr$vparameters),fixed=T),"component"],3)
B1c_analysum$vR <- round(summary( B1c.asr)$var[grepl("!R",names( B1c.asr$vparameters),fixed=T),"component"],3)

B1c_analysum$logl <- B1c.asr$loglik
B1c_analysum$df <- sum(summary(B1c.asr)$var$bound%in%c("P","U"))
B1c_analysum$AIC <- 2*(B1c_analysum$df - B1c_analysum$logl)


B1c_analysum$W <- B1c_shapiro$statistic
B1c_analysum$pW <- B1c_shapiro$p.value
B1c_analysum$nstdres <- nrow(B1c_nstdres)

########################################
# 3.4.4. Test difference bn models
########################################

B1b_v_B1c_tmd_in <- NULL
B1b_v_B1c_tmd_in$analy_name <- 'B1b_v_B1c' 
B1b_v_B1c_tmd_in$fullmod_name <- 'B1b' 
B1b_v_B1c_tmd_in$redmod_name <- 'B1c'
B1b_v_B1c_tmd_in$fullmod.asr <- B1b.asr
B1b_v_B1c_tmd_in$redmod.asr <- B1c.asr
B1b_v_B1c_tmd_in$test <- "boundary" #("boundary"/"other")

B1b_v_B1c_tmd_out <- tmd.f(tmd_in = B1b_v_B1c_tmd_in)
B1b_v_B1c_tmd_out

B1c_analysum$Status <- "BEST"

###########################################################
# 3.5. Fit reduced model (B1d) -  log transofrmed
###########################################################

########################################
# 3.5.2. Fit model
########################################

asreml.options(dense = ~ vm(AID, giv_B1))

B1d.asr = asreml(lnSSC ~ 1 + Y, 
                 random = ~ vm(AID, giv_B1),
                 na.action = na.method(y = "include", x = "include"),
                 data = adata_B1,
                 ai.loadings = TRUE,
                 workspace = "8gb")

########################################
# 3.5.3. Review Fit
########################################

options(scipen = 999)

####################
# 3.5.3.2. Model parameters
####################

wald(B1d.asr)
summary(B1d.asr)$var

####################
# 3.5.3.3. Collate results
####################

B1d_analysum <- data.frame("Model"= "B1d")
B1d_analysum$Scale <- "LN"
B1d_analysum$Status = "FULL"

B1d_analysum$conv <- B1d.asr$converge

B1d_analysum$prY <- round(wald(B1d.asr)[rownames(wald(B1d.asr))=='Y',4],3)
B1d_analysum$vA <- round(summary(B1d.asr)$var[!grepl("vm(AID, giv_B1):Y",names(B1d.asr$vparameters),fixed=T) & grepl("vm(AID, giv_B1)", names(B1d.asr$vparameters),fixed=T),"component"],3)
B1d_analysum$vAxY <- round(summary(B1d.asr)$var[grepl("vm(AID, giv_B1):Y",names(B1d.asr$vparameters),fixed=T),"component"],3)
B1d_analysum$vU <- round(summary( B1d.asr)$var[grepl("!U",names( B1d.asr$vparameters),fixed=T),"component"],3)
B1d_analysum$vR <- round(summary( B1d.asr)$var[grepl("!R",names( B1d.asr$vparameters),fixed=T),"component"],3)

B1d_analysum$logl <- B1d.asr$loglik
B1d_analysum$df <- sum(summary(B1d.asr)$var$bound%in%c("P","U"))
B1d_analysum$AIC <- 2*(B1d_analysum$df - B1d_analysum$logl)

B1d_analysum$W <- B1d_shapiro$statistic
B1d_analysum$pW <- B1d_shapiro$p.value
B1d_analysum$nstdres <- nrow(B1d_nstdres)

########################################
# 3.5.4. Test difference bn models
########################################

B1b_v_B1c_tmd_in <- NULL
B1b_v_B1c_tmd_in$analy_name <- 'B1b_v_B1c' 
B1b_v_B1c_tmd_in$fullmod_name <- 'B1b' 
B1b_v_B1c_tmd_in$redmod_name <- 'B1c'
B1b_v_B1c_tmd_in$fullmod.asr <- B1b.asr
B1b_v_B1c_tmd_in$redmod.asr <- B1c.asr
B1b_v_B1c_tmd_in$test <- "boundary" #("boundary"/"other")

B1b_v_B1c_tmd_out <- tmd.f(tmd_in = B1b_v_B1c_tmd_in)
B1b_v_B1c_tmd_out

B1c_analysum$Status <- "RED"        
B1d_analysum$Status <- "BEST"    

########################################
#3.5.5. Trait heritability and PA
#######################################

summary(B1d.asr)$var

### Trait heritability
h2 <-vpredict(B1d.asr,h2~V1/(V1+V2))
h2

### Extract BLUPs
BLUPs<-summary(B1d.asr,coef=T)$coef.random
BLUP_accessions <-as.matrix(BLUPs)
dim(BLUP_accessions)

### Predicted error variance
mean.PEV<-mean(BLUP_accessions[,2]^2)
mean.PEV

### Calculate sigma2G
sigma2G <-summary(B1d.asr)$varcomp[1,1]
sigma2G

### Calculate accuracy
rB1d <-1-mean.PEV/sigma2G
rB1d



###############################################################################
# 18. Combine IT results and save
###############################################################################

IT_results <- smartbind(B1a_analysum,
                        B1b_analysum,
                        B1c_analysum,
                        B1d_analysum)

## Generate object
currentDate <- Sys.Date()

write.csv(IT_results,paste("IT_results_",currentDate,".csv",sep=""), row.names = F)
