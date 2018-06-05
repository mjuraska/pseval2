# R packages and other R functions used in the CYD14+CYD15 analysis reported in the application section

#' @import graphics
NULL
#' @import stats
NULL
#' @import utils
NULL
#' @import osDesign
NULL
#' @import np
NULL
#' @import chngpt
NULL

library(osDesign)
library(np)
library(chngpt)

logit <- function(p){
  return(log(p/(1-p)))
}

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# 'tpsPredict' returns predicted values from a model fitted by tps
# columns of newMatrix in the same order as the coefficient vector from tps
tpsPredict <- function(fit, newMatrix){
  linPred <- newMatrix %*% fit$coef
  return(drop(1/(1+exp(-linPred))))
}

# 'hNum' returns function values at s0 of the integrand in the numerator of P{Y(0)=1 | S(1)=s1}
# s0 is a numeric vector, whereas s1 is a scalar
hNum <- function(s0, s1, lev, vars, data, tpsFit, changePoint=NULL, npcdensFit1, npcdensFit2){
  ss0 <- s0
  if (!is.null(changePoint)){ ss0 <- ifelse(s0>changePoint, s0-changePoint, 0) }
  phat.s0 <- tpsPredict(tpsFit, cbind(1, ss0, matrix(lev[-1], nrow=length(ss0), ncol=length(lev[-1]), byrow=TRUE)))

  lev1names <- names(which(lev==1))
  predictLevels <- NULL
  for (vi in 1:length(vars)){
    varLevels <- levels(data[,vars[vi]])
    predictLevelVar <- NULL
    for (li in 2:length(varLevels)){
      if (any(grepl(vars[vi], lev1names) & grepl(varLevels[li], lev1names))){
        predictLevelVar <- varLevels[li]
        break
      }
    }
    if (is.null(predictLevelVar)){ predictLevelVar <- varLevels[1] }
    predictLevels <- c(predictLevels, predictLevelVar)
  }

  dfVars <- as.data.frame(matrix(predictLevels, nrow=length(s0), ncol=length(predictLevels), byrow=TRUE))
  colnames(dfVars) <- vars

  fhat.s0 <- predict(npcdensFit1, newdata=data.frame(Sb=s0, dfVars, S=s1))
  ghat.s0 <- predict(npcdensFit2, newdata=data.frame(dfVars, S=s0))
  return(phat.s0*fhat.s0*ghat.s0)
}

# 'hDen' returns function values at s0 of the integrand in the denominator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1 is a scalar
hDen <- function(s0, s1, lev, vars, data, npcdensFit1, npcdensFit2){
  lev1names <- names(which(lev==1))
  predictLevels <- NULL
  for (vi in 1:length(vars)){
    varLevels <- levels(data[,vars[vi]])
    predictLevelVar <- NULL
    for (li in 2:length(varLevels)){
      if (any(grepl(vars[vi], lev1names) & grepl(varLevels[li], lev1names))){
        predictLevelVar <- varLevels[li]
        break
      }
    }
    if (is.null(predictLevelVar)){ predictLevelVar <- varLevels[1] }
    predictLevels <- c(predictLevels, predictLevelVar)
  }

  dfVars <- as.data.frame(matrix(predictLevels, nrow=length(s0), ncol=length(predictLevels), byrow=TRUE))
  colnames(dfVars) <- vars

  fhat.s0 <- predict(npcdensFit1, newdata=data.frame(Sb=s0, dfVars, S=s1))
  ghat.s0 <- predict(npcdensFit2, newdata=data.frame(dfVars, S=s0))
  return(fhat.s0*ghat.s0)
}

propX <- function(X, lev1names){
  # all levels are zero except the intercept
  if (length(lev1names)==1){
    # if no covariates
    if (NCOL(X)==1){
      return(1)
    } else {
      Xsum <- apply(as.matrix(X[,-1]), 1, sum)
      return(mean(Xsum==0))
    }
  } else {
    X1 <- as.matrix(X[,lev1names[-1]])
    Xsum <- apply(X1, 1, sum)
    return(mean(Xsum==NCOL(X1)))
  }
}

# 'riskP' returns the value of P{Y(0)=1 | S(1)=s1}
# s1 is a scalar
# all baseline covariates are assumed to be discrete variables with a finite number of categories
# 'formula' is a one-sided formula of baseline covariates
riskP <- function(s1, formula, data, tpsFit, npcdensFit1, npcdensFit2, changePoint=NULL){
  den <- num <- 0
  UL <- 1.1 * max(data$S, na.rm=TRUE) # if integration over (0,Inf) fails, use (0,UL)

  mf <- model.frame(formula, data)
  X <- model.matrix(terms(formula), mf)
  Xu <- unique(X)
  vars <- all.vars(formula)
  tmp <- NULL
  for (i in 1:NROW(Xu)){
    lev <- drop(Xu[i,])
    names(lev) <- colnames(Xu)
    lev1names <- names(which(lev==1))
    pX <- propX(X, lev1names)
    tmp <- c(tmp, pX)
    hNumInt <- try(integrate(hNum, 0, Inf, s1=s1, lev=lev, vars=vars, data=data, tpsFit=tpsFit, changePoint=changePoint, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value, silent=TRUE)
    if (inherits(hNumInt, 'try-error')){
      num <- num + pX*integrate(hNum, 0, UL, s1=s1, lev=lev, vars=vars, data=data, tpsFit=tpsFit, changePoint=changePoint, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value
    } else {
      num <- num + pX*hNumInt
    }
    den <- den + pX*integrate(hDen, 0, UL, s1=s1, lev=lev, vars=vars, data=data, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000, rel.tol=30*.Machine$double.eps^0.25)$value
  }

  return(list(num/den,tmp))
}

# 'riskV' returns the value of P{Y(1)=1 | S(1)=s1}
# s1 is a scalar
riskV <- function(s1, data, data2, changePoint=NULL){
  nTControls <- NROW(subset(data, Z==1 & Y==0))
  nTCases <- NROW(subset(data, Z==1 & Y==1))
  tpsFormula <- paste0("Y~",ifelse(is.null(changePoint),"S","Sc"))
  fit <- tps(tpsFormula, data=subset(data2, Z==1 & !is.na(Y)), nn0=nTControls, nn1=nTCases, group=rep(1, NROW(subset(data2, Z==1 & !is.na(Y)))), method="PL", cohort=TRUE)
  if (!is.null(changePoint)){ s1 <- ifelse(s1>changePoint,s1-changePoint,0) }
  return(tpsPredict(fit, cbind(1, s1)))
}

# 'risk' returns the estimates of risk in each study group for a given s1
# s1 is a scalar
risk <- function(s1, formula, data, data2, tpsFit, npcdensFit1, npcdensFit2, changePoint=NULL){
  risk1 <- riskV(s1, data, data2, changePoint)
  risk0 <- riskP(s1, formula, data, tpsFit, npcdensFit1, npcdensFit2, changePoint)
  return(list(plaRisk=risk0, txRisk=risk1))
}

# 'bRiskCurve' returns a list; one list component is a matrix (or vector in case iter=1) with rows
# being the estimated VE(s1) curves based on bootstrap samples
# 'data' is assumed to be the ITT set at-risk at month 13 with no prior infection
# 'markerName' is one of "AUC", "S1", "S2", "S3", "S4"
# 'iter' is the number of bootstrap iterations
# 'saveFile' is the name of the .RData file where the output list will be stored
# 'saveDir' specifies the output directory
# 'seed' is an integer for set.seed()
bRiskCurve <- function(data, markerName, iter, saveFile=NULL, saveDir=NULL, seed=NULL){
  # so that the below generic code can be used for each marker
  if (markerName!="AUC"){
    if (markerName=="Min"){  # assumes that there exist variables 'bMin' and 'IMPSTLOG.Min'
      data$bAUC <- data$bMin
      data$IMPSTLOG.AUCMB <- data$IMPSTLOG.Min
    } else {
      data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
      data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
      data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
      data$bAUC <- data[,paste0("b",markerName)]
      data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))]
    }
  }
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))

  # compute the change point for the association of the marker with the dengue disease risk in both the placebo and vaccine groups and consider their minimum
  # as the change point in subsequent logistic regression models
  # use the same change point for each bootstrap sample
  cpointP <- chngptm(formula.1=ofstatus_m13 ~ factor(AGE) + factor(COUNTRY), formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpointV <- chngptm(formula.1=ofstatus_m13 ~ 1, formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==1 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpoint <- min(cpointP, cpointV)

  # the truncated version of IMPSTLOG.AUCMB based on the change point 'cpoint' is used in all subsequent logistic regression models of the dengue disease risk
  data$IMPSTLOG.AUCMB.trunc <- with(data, ifelse(IMPSTLOG.AUCMB>cpoint, IMPSTLOG.AUCMB-cpoint, 0))

  # re-extract the immunogenicity set to include IMPSTLOG.AUCMB.trunc
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))

  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  # in addition, within each treatment group in the immunogenicity set, delete cases to recover
  # the case:control ratio in the ITT set at-risk at month 13 with no prior infection
  dataControls <- subset(data, ofstatus_m0==0)
  nPControlsI <- NROW(dataPControlsI <- subset(dataI, VACC==0 & ofstatus_m0==0))
  nVControlsI <- NROW(dataVControlsI <- subset(dataI, VACC==1 & ofstatus_m0==0))

  dataCases <- subset(data, ofstatus_m13==1)
  nPCasesI <- NROW(dataPCasesI <- subset(dataI, VACC==0 & ofstatus_m13==1))
  nVCasesI <- NROW(dataVCasesI <- subset(dataI, VACC==1 & ofstatus_m13==1))

  nPControls <- NROW(subset(dataControls, VACC==0))
  nVControls <- NROW(subset(dataControls, VACC==1))
  nPCases <- NROW(subset(dataCases, VACC==0))
  nVCases <- NROW(subset(dataCases, VACC==1))

  # within each treatment group, calculate the number of cases in the immunogenicity set needed to achieve
  # the correct case:control ratio
  nPCasesInew <- nPCases * nPControlsI / nPControls
  nVCasesInew <- nVCases * nVControlsI / nVControls

  # within each treatment group, sample as many cases in the immunogenicity set as needed to achieve
  # the correct case:control ratio
  dataPIcorrectRatio <- rbind(dataPControlsI, dataPCasesI[sample(1:nPCasesI, nPCasesInew),])
  dataVIcorrectRatio <- rbind(dataVControlsI, dataVCasesI[sample(1:nVCasesI, nVCasesInew),])
  rm(dataPControlsI); rm(dataPCasesI); rm(dataVControlsI); rm(dataVCasesI)

  # the overall numbers of controls and cases for resampling
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)

  # estimate the optimal bandwidths for kernel density estimation and use these in each bootstrap run
  # RUNS SLOWLY!
  fbw <- npcdensbw(IMPSTLOG.AUCMB ~ bAUC + factor(AGE) + factor(COUNTRY), data=dataVIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  gbw <- npcdensbw(IMPSTLOG.AUCMB ~ factor(AGE) + factor(COUNTRY), data=dataPIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  rm(dataPIcorrectRatio); rm(dataVIcorrectRatio)

  if(!is.null(seed)){ set.seed(seed) }
  bSampleControls <- matrix(sample(1:nControls, nControls*iter, replace=TRUE), nrow=nControls, ncol=iter)
  bSampleCases <- matrix(sample(1:nCases, nCases*iter, replace=TRUE), nrow=nCases, ncol=iter)

  markerVals <- seq(min(dataI$IMPSTLOG.AUCMB), max(dataI$IMPSTLOG.AUCMB), by=0.05)

  # 'bVE' is a list each of whose components is also a list with components:
  # 'VEcurve' - a vector with VE(s) estimates (a single curve)
  # 'bnI'     -  the size of the bootstrapped immunogenicity set
  bRiskCurveList <- lapply(1:iter, function(i){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the bootstrapped immunogenicity set
    bdataI <- subset(bdata, !is.na(IMPSTLOG.AUCMB))

    # compute the change point for the association of the marker with the dengue disease risk in both the placebo and vaccine groups and consider their minimum
    # as the change point in subsequent logistic regression models
    # use the same change point for each bootstrap sample
    bcpointP <- chngptm(formula.1=ofstatus_m13 ~ factor(AGE) + factor(COUNTRY), formula.2=~IMPSTLOG.AUCMB, data=subset(bdataI, VACC==0 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
    bcpointV <- chngptm(formula.1=ofstatus_m13 ~ 1, formula.2=~IMPSTLOG.AUCMB, data=subset(bdataI, VACC==1 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
    bcpoint <- min(bcpointP, bcpointV)

    # the truncated version of IMPSTLOG.AUCMB based on the change point 'cpoint' is used in all subsequent logistic regression models of the dengue disease risk
    bdata$IMPSTLOG.AUCMB.trunc <- with(bdata, ifelse(IMPSTLOG.AUCMB>bcpoint, IMPSTLOG.AUCMB-bcpoint, 0))

    # re-extract the immunogenicity set to include IMPSTLOG.AUCMB.trunc
    bdataI <- subset(bdata, !is.na(IMPSTLOG.AUCMB))

    bdataControls <- subset(bdata, ofstatus_m0==0)
    nPControlsI <- NROW(bdataPControlsI <- subset(bdataI, VACC==0 & ofstatus_m0==0))
    nVControlsI <- NROW(bdataVControlsI <- subset(bdataI, VACC==1 & ofstatus_m0==0))

    bdataCases <- subset(bdata, ofstatus_m13==1)
    nPCasesI <- NROW(bdataPCasesI <- subset(bdataI, VACC==0 & ofstatus_m13==1))
    nVCasesI <- NROW(bdataVCasesI <- subset(bdataI, VACC==1 & ofstatus_m13==1))

    nPControls <- NROW(subset(bdataControls, VACC==0))
    nVControls <- NROW(subset(bdataControls, VACC==1))
    nPCases <- NROW(subset(bdataCases, VACC==0))
    nVCases <- NROW(subset(bdataCases, VACC==1))

    # within each treatment group, calculate the number of cases in the bootstrapped immunogenicity set
    # needed to achieve the correct case:control ratio
    nPCasesInew <- nPCases * nPControlsI / nPControls
    nVCasesInew <- nVCases * nVControlsI / nVControls

    # within each treatment group, sample as many cases in the bootstrapped immunogenicity set as needed
    # to achieve the correct case:control ratio
    bdataPIcorrectRatio <- rbind(bdataPControlsI, bdataPCasesI[sample(1:nPCasesI, nPCasesInew),])
    bdataVIcorrectRatio <- rbind(bdataVControlsI, bdataVCasesI[sample(1:nVCasesI, nVCasesInew),])
    rm(bdataPControlsI); rm(bdataPCasesI); rm(bdataVControlsI); rm(bdataVCasesI)

    group <- rep(1, NROW(subset(bdataI, VACC==0)))

    # weighted logistic regression model using the placebo group in the bootstrapped immunogenicity set
    fit1 <- tps(ofstatus_m13 ~ IMPSTLOG.AUCMB.trunc + factor(AGE) + factor(COUNTRY), data=subset(bdataI, VACC==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)

    # kernel density estimator for f(s1|S_base=sbase, X=x) using the vaccine group in the bootstrapped immunogenicity set
    bfbw <- npcdensbw(IMPSTLOG.AUCMB ~ bAUC + factor(AGE) + factor(COUNTRY), data=bdataVIcorrectRatio, bws=fbw, bandwidth.compute=FALSE)
    fhat <- npcdens(bfbw)

    # kernel density estimator for g(s0|X=x) using the placebo group in the bootstrapped immunogenicity set
    bgbw <- npcdensbw(IMPSTLOG.AUCMB ~ factor(AGE) + factor(COUNTRY), data=bdataPIcorrectRatio, bws=gbw, bandwidth.compute=FALSE)
    ghat <- npcdens(bgbw)

    # a single bootstrap VE(s) curve, risk1(s) curve, and risk0(s) curve
    curves <- lapply(markerVals, function(s){ risk(s, data, dataI, fit1, fhat, ghat, markerName, cpoint) })
    plaRiskCurve <- sapply(curves, "[[", "plaRisk")
    txRiskCurve <- sapply(curves, "[[", "txRisk")

    return(list(plaRiskCurve=plaRiskCurve, txRiskCurve=txRiskCurve))
  })

  # cbind all bootstrap risk curves
  plaRiskCurveBootEst <- drop(do.call(cbind, lapply(bRiskCurveList,"[[","plaRiskCurve")))
  txRiskCurveBootEst <- drop(do.call(cbind, lapply(bRiskCurveList,"[[","txRiskCurve")))
  # bList <- list(markerVals=markerVals, plaRiskCurveBootEst=plaRiskCurveBootEst, txRiskCurveBootEst=txRiskCurveBootEst, bIdxControls=bSampleControls, bIdxCases=bSampleCases,
  #               cpointP=cpointP, cpointV=cpointV, fOptBandwidths=fbw, gOptBandwidths=gbw, seed=seed)
  bList <- list(markerVals=markerVals, plaRiskCurveBootEst=plaRiskCurveBootEst, txRiskCurveBootEst=txRiskCurveBootEst)

  if (!is.null(saveFile)){
    save(bList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n")
  }

  return(invisible(bList))
}

# 'riskCurve' returns the estimated P(Y(0)=1|S(1)=s1) and P(Y(1)=1|S(1)=s1) on a grid of s1 values
# 'data' is assumed to be the ITT set at-risk at month 13 with no prior infection (the target population)
#' @param formula a two-sided formula for the GLM; the LHS is the binary clinical endpoint; the first variable listed on the RHS is the biomarker response at t0 and any variables that follow are baseline covariates
#' @param bsm a character string specifying the variable name representing the baseline surrogate measure
#' @param tx a character string specifying the variable name representing the treatment group indicator
#' @param hinge shall a hinge model be used?
#' @param weights either a numeric vector of weights with values corresponding to rows in \code{data} or a character string specifying the variable name in \code{data}. The weights are passed on to GLMs.
#' @param biomarkerGrid a numeric vector of biomarker values at which the risk curves are evaluated and returned
riskCurve <- function(formula, bsm, tx, data, hinge=FALSE, weights=NULL, biomarkerGrid=NULL, saveFile=NULL, saveDir=NULL){
  if (missing(bsm)){ stop("The variable name in argument 'bsm' for the baseline surrogate measure is missing.") }
  if (missing(tx)){ stop("The variable name in argument 'tx' for the treatment group indicator is missing.") }
  if (missing(data)){ stop("The data frame 'data' for interpreting the variables in 'formula' is missing.") }

  if (!is.null(weights)){
    if (is.character(weights)){
      colnames(data)[colnames(data)==weights] <- "weights"
    } else {
      # assuming 'weights' is a numeric vector
      data$weights <- weights
    }
  }

  formulaDecomp <- strsplit(strsplit(paste(deparse(formula), collapse = ""), " *[~] *")[[1]], " *[+] *")

  # standardize the variable name for the treatment group indicator
  colnames(data)[colnames(data)==tx] <- "Z"

  # standardize the variable name for the biomarker measurement at fixed time t0 post-randomization
  # this line assumes that it is the first listed variable on the RHS of 'formula'
  colnames(data)[colnames(data)==formulaDecomp[[2]][1]] <- "S"

  # standardize the variable name for the biomarker's baseline measurement
  colnames(data)[colnames(data)==bsm] <- "Sb"

  # standardize the variable name for the binary clinical endpoint measured after t0
  colnames(data)[colnames(data)==formulaDecomp[[1]]] <- "Y"

  # extract the subset with available phase 2 data (i.e., with measured S)
  data2 <- subset(data, !is.na(S))

  # case deletion method applied to 'data2' to match the case:control ratio in 'data', separately in each study group,
  # in order to be able to use the standard kernel density estimation procedure
  dataControls <- subset(data, Y==0)
  nPControls <- NROW(subset(dataControls, Z==0))
  nTControls <- NROW(subset(dataControls, Z==1))

  nPControls2 <- NROW(dataPControls2 <- subset(data2, Z==0 & Y==0))
  nTControls2 <- NROW(dataTControls2 <- subset(data2, Z==1 & Y==0))

  dataCases <- subset(data, Y==1)
  nPCases <- NROW(subset(dataCases, Z==0))
  nTCases <- NROW(subset(dataCases, Z==1))

  nPCases2 <- NROW(dataPCases2 <- subset(data2, Z==0 & Y==1))
  nTCases2 <- NROW(dataTCases2 <- subset(data2, Z==1 & Y==1))

  # in each study group separately, calculate the required number of cases in 'data2' to match the case:control ratio in 'data'
  nPCases2new <- nPCases * nPControls2 / nPControls
  nTCases2new <- nTCases * nTControls2 / nTControls

  # in each study group separately, randomly sample cases to match the case:control ratio in 'data'
  dataP2correctRatio <- rbind(dataPControls2, dataPCases2[sample(1:nPCases2, max(1,round(nPCases2new,0))),])
  dataT2correctRatio <- rbind(dataTControls2, dataTCases2[sample(1:nTCases2, max(1,round(nTCases2new,0))),])
  rm(dataPControls2); rm(dataPCases2); rm(dataTControls2); rm(dataTCases2)

  # estimate optimal bandwidths for kernel estimates of the conditional densities
  # fbw <- npcdensbw(as.formula(paste0("S ~ ",paste(c("Sb",formulaDecomp[[2]][-1]),collapse="+"))), data=dataT2correctRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  # save(fbw, file="h:/SCHARP/Sanofi/manuscript_VEcurveMethod/code/test_pseval2/fbw.RData")
  # gbw <- npcdensbw(as.formula(paste0("S ~ ",paste(formulaDecomp[[2]][-1],collapse="+"))), data=dataP2correctRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  # save(gbw, file="h:/SCHARP/Sanofi/manuscript_VEcurveMethod/code/test_pseval2/gbw.RData")

  load("h:/SCHARP/Sanofi/manuscript_VEcurveMethod/code/test_pseval2/fbw.RData")
  load("h:/SCHARP/Sanofi/manuscript_VEcurveMethod/code/test_pseval2/gbw.RData")

  if (hinge){
    # calculate weights for passing on to 'chngptm'
    if (is.null(weights) & !("weights" %in% colnames(data))){
      nControls <- NROW(dataControls)
      nControls2 <- nPControls2 + nTControls2

      nCases <- NROW(dataCases)
      nCases2 <- nPCases2 + nTCases2

      data2$weights <- ifelse(data2$Y==1, nCases/nCases2, nControls/nControls2)
    }

    # in each study group separately, estimate the hinge point for the association of S with Y
    cpointP <- chngptm(formula.1=as.formula(paste0("Y ~ ",paste(formulaDecomp[[2]][-1],collapse="+"))), formula.2=~S, data=subset(data2, Z==0 & !is.na(Y)), family="binomial", type="hinge", prob.weights=weights)$coefficients["chngpt"]
    cpointT <- chngptm(formula.1=Y ~ 1, formula.2=~S, data=subset(data2, Z==1 & !is.na(Y)), family="binomial", type="hinge", prob.weights=weights)$coefficients["chngpt"]
    # use their minimum as the hinge point in the below specified GLMs
    cpoint <- min(cpointP, cpointT)

    # biomarker S left-censored at 'cpoint' for use in the below specified GLMs
    data$Sc <- with(data, ifelse(S>cpoint, S-cpoint, 0))

    # re-extract the subset with phase 2 data
    data2 <- subset(data, !is.na(S))

    # IPW logistic regression model fitted to placebo recipients in the phase 2 subset accounting for two-phase sampling of S
    fit1 <- tps(as.formula(paste0("Y ~ ",paste(c("Sc",formulaDecomp[[2]][-1]),collapse="+"))), data=subset(data2, Z==0 & !is.na(Y)), nn0=nPControls, nn1=nPCases, group=rep(1, NROW(subset(data2, Z==0 & !is.na(Y)))), method="PL", cohort=TRUE)
  } else {
    # for passing on to the function 'risk'
    cpoint <- NULL

    # IPW logistic regression model fitted to placebo recipients in the phase 2 subset accounting for two-phase sampling of S
    fit1 <- tps(as.formula(paste0("Y ~ ",paste(c("S",formulaDecomp[[2]][-1]),collapse="+"))), data=subset(data2, Z==0 & !is.na(Y)), nn0=nPControls, nn1=nPCases, group=rep(1, NROW(subset(data2, Z==0 & !is.na(Y)))), method="PL", cohort=TRUE)
  }

  # kernel density estimator for f(s1|Sb=s0, X=x) using the treatment group in the phase 2 subset
  fhat <- npcdens(fbw)

  # kernel density estimator for g(s0|X=x) using the placebo group in the phase 2 subset
  ghat <- npcdens(gbw)

  # a grid of values of S on which the estimated risk curves are returned
  if (is.null(biomarkerGrid)){
    biomarkerGrid <- seq(min(data2$S), max(data2$S), length.out==200)
  }

  # the first argument of 'risk' is a scalar
  curves <- lapply(biomarkerGrid, function(s){ risk(s, as.formula(paste0("~",paste(formulaDecomp[[2]][-1],collapse="+"))), data, data2, fit1, fhat, ghat, cpoint) })
  plaRiskCurve <- sapply(curves, "[[", "plaRisk")
  txRiskCurve <- sapply(curves, "[[", "txRisk")

  # the output list
  oList <- list(biomarkerGrid=biomarkerGrid, plaRiskCurve=plaRiskCurve, txRiskCurve=txRiskCurve, fOptBandwidths=fbw, gOptBandwidths=gbw)
  if (hinge){
    oList$cpointP <- cpointP
    oList$cpointT <- cpointT
  }

  if (!is.null(saveFile)){
    save(oList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n",file.path(saveDir, saveFile),"\n")
  }

  return(invisible(oList))
}

# returns a two-sided p-value from the test of {H01: mCEP(s1)=CE for all s1} or from the test of {H02: mCEP(s1)=a known 'nullConstant' for s1 in S}
# 'plaRiskCurvePointEst' and 'txRiskCurvePointEst' are numeric vectors
# 'plaRiskCurveBootEst' and 'txRiskCurveBootEst' are length(markerVals) x nBoot matrices
# 'plaRiskOverall' and 'txRiskOverall' are numeric values and need to be specified if null=="H01"
# 'tMCEPconstantNull' is a numeric value and needs to be specified if null=="H02"
# 'MCEPcontrast' must be one of "multiplicativeTE" or "additiveTE"
# 'null' must be one of "H01" and "H02"
# 'S1' is a numeric vector specifying a grid of marker values on which the risk curve estimates are calculated; it needs to be specified if 'limS' is specified
# 'limS1' is a numeric vector specifying the minimum and maximum S(1) value in H02 and must be specified if null=="H02"
testConstancy2 <- function(plaRiskCurvePointEst, txRiskCurvePointEst, plaRiskCurveBootEst, txRiskCurveBootEst, plaRiskOverall=NULL, txRiskOverall=NULL, tMCEPconstantNull=NULL,
                           MCEPcontrast, null, S1=NULL, limS1=NULL){
  # trim the risk curves if 'limS' is specified
  if (!is.null(limS1)){
    if (is.null(S1)){ stop("'S1' and 'limS1' must be specified concurrently.") }

    plaRiskCurvePointEst <- plaRiskCurvePointEst[S1>=limS1[1] & S1<=limS1[2]]
    txRiskCurvePointEst <- txRiskCurvePointEst[S1>=limS1[1] & S1<=limS1[2]]
    plaRiskCurveBootEst <- plaRiskCurveBootEst[S1>=limS1[1] & S1<=limS1[2],]
    txRiskCurveBootEst <- txRiskCurveBootEst[S1>=limS1[1] & S1<=limS1[2],]
  }

  if (MCEPcontrast=="multiplicativeTE"){
    # transformed estimated MCEP curve
    tMCEP <- log(txRiskCurvePointEst/plaRiskCurvePointEst)
    # transformed bootstrapped MCEP curves
    tbMCEP <- log(txRiskCurveBootEst/plaRiskCurveBootEst)
  }

  if (MCEPcontrast=="additiveTE"){
    # transformed MCEP curve
    tMCEP <- logit(((plaRiskCurvePointEst - txRiskCurvePointEst) + 1)/2)
    # transformed bootstrapped MCEP curves
    tbMCEP <- logit(((plaRiskCurveBootEst - txRiskCurveBootEst) + 1)/2)
  }

  # bootstrap SE of tMCEP estimates
  bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)

  # calculate the supremum statistic for each bootstrap sample
  supAbsZ <- NULL
  for (j in 1:NCOL(tbMCEP)){
    Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
  }

  if (null=="H01"){
    if (is.null(plaRiskOverall) | is.null(txRiskOverall)){ stop("'plaRiskOverall' and 'txRiskOverall' must be specified for the test of H01.") }

    if (MCEPcontrast=="multiplicativeTE"){
      tCEest <- log(txRiskOverall/plaRiskOverall)
    }

    if (MCEPcontrast=="additiveTE"){
      tCEest <- logit(((plaRiskOverall - txRiskOverall) + 1)/2)
    }

    testStat <- max(abs(tMCEP-tCEest)/bSE, na.rm=TRUE)

    return(mean(supAbsZ > testStat))
  }

  if (null=="H02"){
    if (is.null(tMCEPconstantNull)){ stop("'tMCEPconstantNull' must be specified for the test of H02.") }

    testStat <- max(abs(tMCEP-tMCEPconstantNull)/bSE, na.rm=TRUE)

    return(mean(supAbsZ > testStat))
  }
}

plotBvsM13inP <- function(data, markerName){
  require(robustbase)
  markerName2 <- switch(markerName, AUC="AUCMB", Min="Min", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, AUC="AUC-MB", Min="minimum titer", S1="$\\log_{10}$ serotype 1 titer", S2="$\\log_{10}$ serotype 2 titer", S3="$\\log_{10}$ serotype 3 titer", S4="$\\log_{10}$ serotype 4 titer")
  markerX <- switch(markerName, AUC="AUC-MB", Min="Minimum Titer", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  if (markerName!="AUC"){
    if (markerName=="Min"){  # assumes that there exist variables 'bMin' and 'IMPSTLOG.Min'
      data$bAUC <- data$bMin
      data$IMPSTLOG.AUCMB <- data$IMPSTLOG.Min
    } else {
      data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
      data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
      data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
      data$bAUC <- data[,paste0("b",markerName)]
      data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))]
    }
  }

  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  dataP <- subset(dataI, VACC==0 & !is.na(bAUC)) # subset of the placebo group in the imm set with baseline markers
  par(mar=c(4,5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  with(dataP, plot(bAUC, IMPSTLOG.AUCMB, xlab=paste0("Baseline ",markerX," in Placebo Recipients"),
                   ylab=paste0("Month 13 ",markerX,ifelse(markerName=="AUC","","\n")," in Placebo Recipients"),
                   cex=0.7))
  abline(0,1, lty="dotted", lwd=2)
  abline(lmrob(IMPSTLOG.AUCMB ~ bAUC, data=dataP, maxit.scale=500), lwd=2.5, col="darkgoldenrod2")
  with(dataP, lines(lowess(bAUC, IMPSTLOG.AUCMB), col="red", lwd=2, lty="dashed"))
  legend("topleft", lwd=2, lty=c("solid","dashed","dotted"), col=c("darkgoldenrod2", "red", "black"),
         legend=c("Robust linear reg (Yohai, 1987)", "LOWESS", "y=x"), bty="n", cex=0.7)
  legend("bottomright", legend=paste0("Spearman's r=",with(dataP, round(cor(bAUC, IMPSTLOG.AUCMB, method="spearman"), 3))),
         cex=0.8, bty="n")
  #cat("\\caption{Association between $S(0)$ and $S_b$ measurements of ",markerName3," in arm $Z=0$ in the CYD14 trial}\\label{Fig: Sbase vs S0 ",markerName2,"}", sep="")
}

# 'plaRiskCurvePointEst' and 'txRiskCurvePointEst' are vectors
# 'plaRiskCurveBootEst' and 'txRiskCurveBootEst' are length(markerVals)-by-1000 matrices
plotLogRRcurve <- function(markerVals, plaRiskCurvePointEst, txRiskCurvePointEst, plaRiskCurveBootEst=NULL, txRiskCurveBootEst=NULL, title=NULL,
                           hingePoint=NULL, plotLegend=TRUE, yLim=NULL){
  MCEPcurvePointEst <- log(txRiskCurvePointEst/plaRiskCurvePointEst)

  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    # transformed MCEP curve (identity)
    tMCEP <- MCEPcurvePointEst
    # transformed bootstrapped MCEP curves (identity)
    # assuming the matrices have the same dimensions
    tbMCEP <- log(txRiskCurveBootEst/plaRiskCurveBootEst)

    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)

    # pointwise confidence bounds for MCEP(s1)
    ptLB.MCEP <- tMCEP - qnorm(0.975) * bSE
    ptUB.MCEP <- tMCEP + qnorm(0.975) * bSE

    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)

    smLB.MCEP <- tMCEP - qSupAbsZ * bSE
    smUB.MCEP <- tMCEP + qSupAbsZ * bSE
  } else {
    ptLB.MCEP <- ptUB.MCEP <- smLB.MCEP <- smUB.MCEP <- NULL
  }

  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2
  if (is.null(yLim)){ yLim <- range(c(MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP), na.rm=TRUE) }

  par(mar=c(5,5,1.5,5), cex.lab=cexLab, cex.axis=cexAxis, las=1)

  plot(markerVals, MCEPcurvePointEst, type="n", xlab="Month 13 Average Titer of Vaccinees", ylab="", xlim=range(markerVals),
       ylim=yLim, lwd=3.5, xaxt="n", yaxt="n")
  axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
  #axis(side=2, at=seq(-5,0,by=0.5), cex.axis=cexAxis)
  axis(side=2, at=seq(-15,3,by=1), cex.axis=cexAxis)
  mtext("Log Relative Risk", side=2, las=0, line=3.4, cex=cexLab)
  #axis(side=4, at=log(1-c(0,0.3,0.5,0.7,0.8,0.9,0.95,0.97,0.99)), labels=c(0,0.3,0.5,0.7,0.8,0.9,0.95,0.97,0.99)*100, cex.axis=cexAxis)
  axis(side=4, at=log(1-c(-1,0,0.5,0.8,0.9,0.95,0.99)), labels=c(-1,0,0.5,0.8,0.9,0.95,0.99)*100, cex.axis=cexAxis)
  mtext("Vaccine Efficacy (%)", side=4, las=0, line=3.2, cex=cexLab)
  if (!is.null(title)){ mtext(title, side=3, cex=cexTitle, line=0) }

  abline(h=0, lty="dotted", lwd=2, col="gray50")

  lines(markerVals, MCEPcurvePointEst, lwd=3.5)

  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    lines(markerVals, ptLB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, ptUB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, smLB.MCEP, lty="dotdash", lwd=3)
    lines(markerVals, smUB.MCEP, lty="dotdash", lwd=3)
  }

  if (plotLegend){ legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c("Pointwise 95% CI","Simultaneous 95% CI"), cex=cexLegend, bty="n") }
  if (!is.null(hingePoint)){ legend("topright", paste0("Hinge Point = ", hingePoint,"  "), cex=cexLegend, bty="n") }

  text(3, -0.7, expression(H[0]^1: p < 0.001), cex=cexLegend, pos=4)
  text((log10(5)+log10(57))/2, -2, expression(H[0]^2: p < 0.001), cex=cexLegend)
  segments(x0=log10(5), x1=log10(57), y0=-1.4, lwd=2)
  segments(x0=log10(5), y0=-1, y1=-1.4, lwd=2)
  segments(x0=log10(57), y0=-1, y1=-1.4, lwd=2)
}

# 'plaRiskCurvePointEst' and 'txRiskCurvePointEst' are vectors
# 'plaRiskCurveBootEst' and 'txRiskCurveBootEst' are length(markerVals)-by-1000 matrices
plotRiskDiffCurve <- function(markerVals, plaRiskCurvePointEst, txRiskCurvePointEst, plaRiskCurveBootEst=NULL, txRiskCurveBootEst=NULL, title=NULL,
                              hingePoint=NULL, plotLegend=TRUE, yLim=NULL){
  MCEPcurvePointEst <- plaRiskCurvePointEst - txRiskCurvePointEst

  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    # transformed MCEP curve
    tMCEP <- logit((MCEPcurvePointEst + 1)/2)
    # transformed bootstrapped MCEP curves
    # assuming the matrices have the same dimensions
    tbMCEP <- logit((plaRiskCurveBootEst - txRiskCurveBootEst + 1)/2)

    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)

    # pointwise confidence bounds for MCEP(s1)
    ptLB.MCEP <- 2*expit(tMCEP - qnorm(0.975) * bSE) - 1
    ptUB.MCEP <- 2*expit(tMCEP + qnorm(0.975) * bSE) - 1

    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)

    smLB.MCEP <- 2*expit(tMCEP - qSupAbsZ * bSE) - 1
    smUB.MCEP <- 2*expit(tMCEP + qSupAbsZ * bSE) - 1
  } else {
    ptLB.MCEP <- ptUB.MCEP <- smLB.MCEP <- smUB.MCEP <- NULL
  }

  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2
  if (is.null(yLim)){ yLim <- range(c(MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP), na.rm=TRUE) }

  par(mar=c(5,5,1.5,5), cex.lab=cexLab, cex.axis=cexAxis, las=1)

  plot(markerVals, MCEPcurvePointEst, type="n", xlab="Month 13 Average Titer of Vaccinees", ylab="", xlim=range(markerVals),
       ylim=yLim, lwd=3.5, xaxt="n", yaxt="n")
  axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
  # axis(side=2, at=seq(0,0.05,by=0.005), labels=FALSE, cex.axis=cexAxis)
  # axis(side=2, at=seq(0,0.05,by=0.005), cex.axis=cexAxis, line=-0.3, tick=FALSE)
  axis(side=2, at=seq(-0.25,0.3,by=0.05), labels=FALSE, cex.axis=cexAxis)
  axis(side=2, at=seq(-0.25,0.3,by=0.05), cex.axis=cexAxis, line=-0.3, tick=FALSE)
  mtext("Risk Difference (Placebo - Vaccine)", side=2, las=0, line=3.8, cex=cexLab)
  if (!is.null(title)){ mtext(title, side=3, cex=cexTitle, line=0) }

  abline(h=0, lty="dotted", lwd=2, col="gray50")

  lines(markerVals, MCEPcurvePointEst, lwd=3.5)

  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    lines(markerVals, ptLB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, ptUB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, smLB.MCEP, lty="dotdash", lwd=3)
    lines(markerVals, smUB.MCEP, lty="dotdash", lwd=3)
  }

  if (plotLegend){ legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c("Pointwise 95% CI","Simultaneous 95% CI"), cex=cexLegend, bty="n") }
  if (!is.null(hingePoint)){ legend("topright", paste0("Hinge Point = ", hingePoint,"  "), cex=cexLegend, bty="n") }

  text(2.5, 0.15, expression(H[0]^1: p == 0.16), cex=cexLegend, pos=4)
  text((log10(5)+log10(57))/2, 0.1, expression(H[0]^2: p < 0.001), cex=cexLegend)
  segments(x0=log10(5), x1=log10(57), y0=0.07, lwd=2)
  segments(x0=log10(5), y0=0.05, y1=0.07, lwd=2)
  segments(x0=log10(57), y0=0.05, y1=0.07, lwd=2)
}

plotLogRRcurvePSN <- function(markerVals, MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP, title=NULL, hingePoint=NULL, plotLegend=TRUE){
  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2

  par(mar=c(5,5,1.5,5), cex.lab=cexLab, cex.axis=cexAxis, las=1)

  plot(markerVals, MCEPcurvePointEst, type="n", xlab="Month 13 Average Titer of Vaccinees", ylab="", xlim=range(markerVals),
       ylim=range(c(MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP), na.rm=TRUE), lwd=3.5, xaxt="n", yaxt="n")
  axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
  axis(side=2, at=seq(-15,3,by=1), cex.axis=cexAxis)
  mtext("Log Relative Risk", side=2, las=0, line=3.4, cex=cexLab)
  axis(side=4, at=log(1-c(-1,0,0.5,0.8,0.9,0.95,0.99)), labels=c(-1,0,0.5,0.8,0.9,0.95,0.99)*100, cex.axis=cexAxis)
  mtext("Vaccine Efficacy (%)", side=4, las=0, line=3.2, cex=cexLab)
  if (!is.null(title)){ mtext(title, side=3, cex=cexTitle, line=0) }

  abline(h=0, lty="dotted", lwd=2, col="gray50")

  lines(markerVals, MCEPcurvePointEst, lwd=3.5)

  lines(markerVals, ptLB.MCEP, lty="dashed", lwd=3)
  lines(markerVals, ptUB.MCEP, lty="dashed", lwd=3)
  lines(markerVals, smLB.MCEP, lty="dotdash", lwd=3)
  lines(markerVals, smUB.MCEP, lty="dotdash", lwd=3)

  if (plotLegend){ legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c("Pointwise 95% CI","Simultaneous 95% CI"), cex=cexLegend, bty="n") }
  if (!is.null(hingePoint)){ legend("topright", paste0("Hinge Point = ", hingePoint,"  "), cex=cexLegend, bty="n") }
}

plotRiskDiffCurvePSN <- function(markerVals, MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP, title=NULL, hingePoint=NULL, plotLegend=TRUE){
  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2

  par(mar=c(5,5,1.5,5), cex.lab=cexLab, cex.axis=cexAxis, las=1)

  plot(markerVals, MCEPcurvePointEst, type="n", xlab="Month 13 Average Titer of Vaccinees", ylab="", xlim=range(markerVals),
       ylim=range(c(MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP), na.rm=TRUE), lwd=3.5, xaxt="n", yaxt="n")
  axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
  axis(side=2, at=seq(-0.25,0.3,by=0.05), labels=FALSE, cex.axis=cexAxis)
  axis(side=2, at=seq(-0.25,0.3,by=0.05), cex.axis=cexAxis, line=-0.3, tick=FALSE)
  mtext("Risk Difference (Placebo - Vaccine)", side=2, las=0, line=3.8, cex=cexLab)
  if (!is.null(title)){ mtext(title, side=3, cex=cexTitle, line=0) }

  abline(h=0, lty="dotted", lwd=2, col="gray50")

  lines(markerVals, MCEPcurvePointEst, lwd=3.5)

  lines(markerVals, ptLB.MCEP, lty="dashed", lwd=3)
  lines(markerVals, ptUB.MCEP, lty="dashed", lwd=3)
  lines(markerVals, smLB.MCEP, lty="dotdash", lwd=3)
  lines(markerVals, smUB.MCEP, lty="dotdash", lwd=3)

  if (plotLegend){ legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c("Pointwise 95% CI","Simultaneous 95% CI"), cex=cexLegend, bty="n") }
  if (!is.null(hingePoint)){ legend("topright", paste0("Hinge Point = ", hingePoint,"  "), cex=cexLegend, bty="n") }
}

# 'applyBoot' calculates the bootstrap confidence bands as pointwise quantiles from the large number of VE curves
# the purpose of this function is to speed up plotting and report generation and to avoid running out of memory
# this function creates a new .RData file that can be loaded in a plotting function
applyBoot <- function(loadFile, saveFile, saveDir){
  load(file.path(saveDir, loadFile))
  bVE <- NULL
  for (i in 1:length(results)){
    if (!is.null(names(results[[i]]))){
      bVE <- rbind(bVE, results[[i]]$bVE)
    }
  }
  #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))

  s <- results[[1]]$markerVals
  UB <- apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)
  LB <- apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)
  bList <- list(s=s, UB=UB, LB=LB)
  save(bList, file=file.path(saveDir,saveFile))
}

npcdensbw.formula <- function (formula, data, subset, na.action, call, ...){
  orig.class <- if (missing(data))
    sapply(eval(attr(terms(formula), "variables"), environment(formula)),
           class)
  else sapply(eval(attr(terms(formula), "variables"), data,
                   environment(formula)), class)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf),
             nomatch = 0)
  mf <- mf[c(1, m)]
  if (!missing(call) && is.call(call)) {
    for (i in 1:length(call)) {
      if (tryCatch(class(eval(call[[i]])) == "formula",
                   error = function(e) FALSE))
        break
    }
    mf[[2]] <- call[[i]]
  }
  mf[[1]] <- as.name("model.frame")
  #if (m[2] > 0) { # use data as environment
  #  mf[["formula"]] = eval(mf[[m[1]]], environment(mf[[m[2]]]))
  #} else { # use parent frame
  mf[["formula"]] = eval(mf[[m[1]]], parent.frame())
  #}
  variableNames <- np:::explodeFormula(mf[["formula"]])
  varsPlus <- lapply(variableNames, paste, collapse = " + ")
  mf[["formula"]] <- as.formula(paste(" ~ ", varsPlus[[1]],
                                      " + ", varsPlus[[2]]), env = environment(formula))
  mf[["formula"]] <- terms(mf[["formula"]])
  if (all(orig.class == "ts")) {
    args <- (as.list(attr(mf[["formula"]], "variables"))[-1])
    attr(mf[["formula"]], "predvars") <- as.call(c(quote(as.data.frame),
                                                   as.call(c(quote(ts.intersect), args))))
  }
  else if (any(orig.class == "ts")) {
    arguments <- (as.list(attr(mf[["formula"]], "variables"))[-1])
    arguments.normal <- arguments[which(orig.class != "ts")]
    arguments.timeseries <- arguments[which(orig.class ==
                                              "ts")]
    ix <- sort(c(which(orig.class == "ts"), which(orig.class !=
                                                    "ts")), index.return = TRUE)$ix
    attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(cbind),
                                                            as.call(c(quote(as.data.frame), as.call(c(quote(ts.intersect),
                                                                                                      arguments.timeseries)))), arguments.normal, check.rows = TRUE)))[,
                                                                                                                                                                       .(ix)])
  }
  mf <- eval(mf, parent.frame())
  ydat <- mf[, variableNames[[1]], drop = FALSE]
  xdat <- mf[, variableNames[[2]], drop = FALSE]
  tbw = npcdensbw(xdat = xdat, ydat = ydat, ...)
  tbw$call <- match.call(expand.dots = FALSE)
  environment(tbw$call) <- parent.frame()
  tbw$formula <- formula
  tbw$rows.omit <- as.vector(attr(mf, "na.action"))
  tbw$nobs.omit <- length(tbw$rows.omit)
  tbw$terms <- attr(mf, "terms")
  tbw$variableNames <- variableNames
  tbw
}
