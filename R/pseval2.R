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
#' @import MASS
NULL

library(osDesign)
library(np)
library(chngpt)
library(MASS)

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

  # if any baseline covariates are specified
  if (length(vars)>=1){
    phat.s0 <- tpsPredict(tpsFit, cbind(1, ss0, matrix(lev[-1], nrow=length(ss0), ncol=length(lev[-1]), byrow=TRUE)))

    lev1names <- names(which(lev==1))
    predictLevels <- NULL
    for (vi in 1:length(vars)){
      varLevels <- levels(as.factor(data[,vars[vi]]))
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
  } else {
    phat.s0 <- tpsPredict(tpsFit, cbind(1, ss0))
    fhat.s0 <- predict(npcdensFit1, newdata=data.frame(Sb=s0, S=s1))
    ghat.s0 <- predict(npcdensFit2, newdata=data.frame(S=s0))
  }

  return(phat.s0*fhat.s0*ghat.s0)
}

# 'hDen' returns function values at s0 of the integrand in the denominator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1 is a scalar
hDen <- function(s0, s1, lev, vars, data, npcdensFit1, npcdensFit2){
  # if any baseline covariates are specified
  if (length(vars)>=1){
    lev1names <- names(which(lev==1))
    predictLevels <- NULL
    for (vi in 1:length(vars)){
      varLevels <- levels(as.factor(data[,vars[vi]]))
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
  } else {
    fhat.s0 <- predict(npcdensFit1, newdata=data.frame(Sb=s0, S=s1))
    ghat.s0 <- predict(npcdensFit2, newdata=data.frame(S=s0))
  }
  return(fhat.s0*ghat.s0)
}

rowMatch <- function(x, y){
  return(sum(abs(x-y))==0)
}

propX <- function(X, lev){
  Xmatch <- apply(as.matrix(X), 1, rowMatch, y=lev)
  return(mean(Xmatch))
}

# 'riskP' returns the value of P{Y(0)=1 | S(1)=s1}
# s1 is a scalar
# all baseline covariates are assumed to be discrete variables with a finite number of categories
# 'formula' is a one-sided formula of baseline covariates
riskP <- function(s1, formula, data, tpsFit, npcdensFit1, npcdensFit2, changePoint=NULL){
  den <- num <- 0
  UL <- 1.05 * max(data$S, na.rm=TRUE) # if integration over (0,Inf) fails, use (0,UL)

  mf <- model.frame(formula, data)
  X <- model.matrix(terms(formula), mf)
  Xu <- unique(X)
  vars <- all.vars(formula)
  for (i in 1:NROW(Xu)){
    lev <- drop(Xu[i,])
    names(lev) <- colnames(Xu)
    pX <- propX(X, lev)
    hNumInt <- try(integrate(hNum, 0, Inf, s1=s1, lev=lev, vars=vars, data=data, tpsFit=tpsFit, changePoint=changePoint, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value, silent=TRUE)
    if (inherits(hNumInt, 'try-error')){
      num <- num + pX*integrate(hNum, 0, UL, s1=s1, lev=lev, vars=vars, data=data, tpsFit=tpsFit, changePoint=changePoint, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value
    } else {
      num <- num + pX*hNumInt
    }
    den <- den + pX*integrate(hDen, 0, UL, s1=s1, lev=lev, vars=vars, data=data, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000, rel.tol=30*.Machine$double.eps^0.25)$value
  }

  return(num/den)
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
bootRiskCurve <- function(formula, bsm, tx, data, hinge=FALSE, weights=NULL, biomarkerGrid=NULL, iter, saveFile=NULL, saveDir=NULL, seed=NULL){
  if (missing(bsm)){ stop("The variable name in argument 'bsm' for the baseline surrogate measure is missing.") }
  if (missing(tx)){ stop("The variable name in argument 'tx' for the treatment group indicator is missing.") }
  if (missing(data)){ stop("The data frame 'data' for interpreting the variables in 'formula' is missing.") }
  if (missing(iter)){ stop("The number of bootstrap iterations in argument 'iter' is missing.") }

  if (!is.null(weights)){
    if (is.character(weights)){
      colnames(data)[colnames(data)==weights] <- "weights"
    } else {
      # assuming 'weights' is a numeric vector
      data$weights <- weights
    }
  }

  formulaDecomp <- strsplit(strsplit(paste(deparse(formula), collapse = ""), " *[~] *")[[1]], " *[+] *")
  anyBaselineCovar <- length(formulaDecomp[[2]])>1

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
  nControls <- NROW(dataControls)
  nPControls <- NROW(subset(dataControls, Z==0))
  nTControls <- NROW(subset(dataControls, Z==1))

  nPControls2 <- NROW(dataPControls2 <- subset(data2, Z==0 & Y==0))
  nTControls2 <- NROW(dataTControls2 <- subset(data2, Z==1 & Y==0))

  dataCases <- subset(data, Y==1)
  nCases <- NROW(dataCases)
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
  if (anyBaselineCovar){
    fm.fbw <- as.formula(paste0("S ~ ",paste(c("Sb",formulaDecomp[[2]][-1]),collapse="+")))
    fm.gbw <- as.formula(paste0("S ~ ",paste(formulaDecomp[[2]][-1],collapse="+")))
    gbw <- npcdensbw(fm.gbw, data=dataP2correctRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  } else {
    fm.fbw <- S ~ Sb
    fm.gbw <- ~ S
    # marginal density
    gbw <- npudensbw(fm.gbw, data=dataP2correctRatio, ckertype="epanechnikov")
  }
  fbw <- npcdensbw(fm.fbw, data=dataT2correctRatio, cxkertype="epanechnikov", cykertype="epanechnikov")

  if(!is.null(seed)){ set.seed(seed) }
  bSampleControls <- matrix(sample(1:nControls, nControls*iter, replace=TRUE), nrow=nControls, ncol=iter)
  bSampleCases <- matrix(sample(1:nCases, nCases*iter, replace=TRUE), nrow=nCases, ncol=iter)

  # a grid of values of S on which the bootstrapped risk curves are returned
  if (is.null(biomarkerGrid)){
    biomarkerGrid <- seq(min(data2$S), max(data2$S), length.out=200)
  }
  rm(data2)

  bRiskCurveList <- lapply(1:iter, function(i, formulaDecomp){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the bootstrap subset with phase 2 data
    bdata2 <- subset(bdata, !is.na(S))

    # estimate the hinge point in each bootstrap iteration
    if (hinge){
      # recalculate weights in each bootstrap iteration for passing on to 'chngptm'
      if (is.null(weights)){
        nControls2 <- NROW(subset(bdata2, Y==0))
        nCases2 <- NROW(subset(bdata2, Y==1))

        bdata2$weights <- ifelse(bdata2$Y==1, nCases/nCases2, nControls/nControls2)
      }

      # in each study group separately, estimate the hinge point for the association of S with Y
      if (anyBaselineCovar){
        fm <- as.formula(paste0("Y ~ ",paste(formulaDecomp[[2]][-1],collapse="+")))
      } else {
        fm <- Y ~ 1
      }
      cpointP <- chngptm(formula.1=fm, formula.2=~S, data=subset(bdata2, Z==0 & !is.na(Y)), family="binomial", type="hinge", prob.weights=weights)$coefficients["chngpt"]
      cpointT <- chngptm(formula.1=Y ~ 1, formula.2=~S, data=subset(bdata2, Z==1 & !is.na(Y)), family="binomial", type="hinge", prob.weights=weights)$coefficients["chngpt"]
      # use their minimum as the hinge point in the below specified GLMs
      cpoint <- min(cpointP, cpointT)

      # biomarker S left-censored at 'cpoint' for use in the below specified GLMs
      bdata$Sc <- with(bdata, ifelse(S>cpoint, S-cpoint, 0))

      # re-extract the subset with phase 2 data
      bdata2 <- subset(bdata, !is.na(S))

      # IPW logistic regression model fitted to placebo recipients in the phase 2 subset accounting for two-phase sampling of S
      if (anyBaselineCovar){
        fm <- as.formula(paste0("Y ~ ",paste(c("Sc",formulaDecomp[[2]][-1]),collapse="+")))
      } else {
        fm <- Y ~ Sc
      }
      fit1 <- tps(fm, data=subset(bdata2, Z==0 & !is.na(Y)), nn0=NROW(subset(bdata, Z==0 & Y==0)), nn1=NROW(subset(bdata, Z==0 & Y==1)), group=rep(1, NROW(subset(bdata2, Z==0 & !is.na(Y)))), method="PL", cohort=TRUE)
    } else {
      # for passing on to the function 'risk'
      cpoint <- NULL

      # IPW logistic regression model fitted to placebo recipients in the phase 2 subset accounting for two-phase sampling of S
      if (anyBaselineCovar){
        fm <- as.formula(paste0("Y ~ ",paste(c("S",formulaDecomp[[2]][-1]),collapse="+")))
      } else {
        fm <- Y ~ S
      }
      fit1 <- tps(fm, data=subset(bdata2, Z==0 & !is.na(Y)), nn0=NROW(subset(bdata, Z==0 & Y==0)), nn1=NROW(subset(bdata, Z==0 & Y==1)), group=rep(1, NROW(subset(bdata2, Z==0 & !is.na(Y)))), method="PL", cohort=TRUE)
    }

    bdataControls <- subset(bdata, Y==0)
    nPControls2 <- NROW(bdataPControls2 <- subset(bdata2, Z==0 & Y==0))
    nTControls2 <- NROW(bdataTControls2 <- subset(bdata2, Z==1 & Y==0))

    bdataCases <- subset(bdata, Y==1)
    nPCases2 <- NROW(bdataPCases2 <- subset(bdata2, Z==0 & Y==1))
    nTCases2 <- NROW(bdataTCases2 <- subset(bdata2, Z==1 & Y==1))

    nPControls <- NROW(subset(bdataControls, Z==0))
    nTControls <- NROW(subset(bdataControls, Z==1))
    nPCases <- NROW(subset(bdataCases, Z==0))
    nTCases <- NROW(subset(bdataCases, Z==1))

    # in each study group separately, recalculate the required number of cases in 'bdata2' to match the case:control ratio in 'bdata'
    nPCases2new <- nPCases * nPControls2 / nPControls
    nTCases2new <- nTCases * nTControls2 / nTControls

    # in each study group separately, randomly sample cases to match the case:control ratio in 'bdata'
    bdataP2correctRatio <- rbind(bdataPControls2, bdataPCases2[sample(1:nPCases2, max(1,round(nPCases2new,0))),])
    bdataT2correctRatio <- rbind(bdataTControls2, bdataTCases2[sample(1:nTCases2, max(1,round(nTCases2new,0))),])
    rm(bdataPControls2); rm(bdataPCases2); rm(bdataTControls2); rm(bdataTCases2)

    # kernel density estimator for g(s0|X=x) using the placebo group in the phase 2 subset
    if (anyBaselineCovar){
      # the formulas are reintroduced in the 'lapply' so that 'npcdensbw' and 'npudensbw' are able to evaluate them in the parent.frame() environment
      fm.fbw <- as.formula(paste0("S ~ ",paste(c("Sb",formulaDecomp[[2]][-1]),collapse="+")))
      fm.gbw <- as.formula(paste0("S ~ ",paste(formulaDecomp[[2]][-1],collapse="+")))
      bgbw <- npcdensbw(fm.gbw, data=bdataP2correctRatio, bws=gbw, bandwidth.compute=FALSE)
      ghat <- npcdens(bgbw)
    } else {
      fm.fbw <- S ~ Sb
      fm.gbw <- ~ S
      # marginal density
      bgbw <- npudensbw(fm.gbw, data=bdataP2correctRatio, bws=gbw, bandwidth.compute=FALSE)
      ghat <- npudens(bgbw)
    }

    # kernel density estimator for f(s1|Sb=s0, X=x) using the treatment group in the phase 2 subset
    bfbw <- npcdensbw(fm.fbw, data=bdataT2correctRatio, bws=fbw, bandwidth.compute=FALSE)
    fhat <- npcdens(bfbw)

    # the first argument of 'risk' is a scalar
    if (anyBaselineCovar){
      fm <- as.formula(paste0("~",paste(formulaDecomp[[2]][-1],collapse="+")))
    } else {
      fm <- ~ 1
    }
    curves <- lapply(biomarkerGrid, function(s){ risk(s, fm, bdata, bdata2, fit1, fhat, ghat, cpoint) })
    plaRiskCurve <- sapply(curves, "[[", "plaRisk")
    txRiskCurve <- sapply(curves, "[[", "txRisk")

    # the output list
    out <- list(plaRiskCurve=plaRiskCurve, txRiskCurve=txRiskCurve)
    if (hinge){
      out$cpointP <- cpointP
      out$cpointT <- cpointT
    }

    return(out)
  }, formulaDecomp=formulaDecomp)

  # cbind all bootstrap risk curves
  plaRiskCurveBootEst <- drop(do.call(cbind, lapply(bRiskCurveList,"[[","plaRiskCurve")))
  txRiskCurveBootEst <- drop(do.call(cbind, lapply(bRiskCurveList,"[[","txRiskCurve")))

  # the output list
  bList <- list(biomarkerGrid=biomarkerGrid, plaRiskCurveBootEst=plaRiskCurveBootEst, txRiskCurveBootEst=txRiskCurveBootEst)
  if (hinge){
    bList$cpointPbootEst <- sapply(bRiskCurveList,"[[","cpointP")
    bList$cpointTbootEst <- sapply(bRiskCurveList,"[[","cpointT")
  }

  if (!is.null(saveFile)){
    save(bList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n",file.path(saveDir, saveFile),"\n")
  }

  return(invisible(bList))
}

#' Estimation of Conditional Clinical Endpoint Risk under Placebo and Treatment Given Biomarker Response to Treatment in a Three-Phase Sampling Design
#'
#' Estimates \eqn{P\{Y(z)=1|S(1)=s_1\}}, \eqn{z=0,1}, on a grid of \eqn{s_1} values following the estimation method of Juraska, Huang, and Gilbert (2018), where \eqn{Z} is the
#' treatment group indicator (\eqn{Z=1}, treatment; \eqn{Z=0}, placebo), \eqn{S(z)} is a discrete or continuous univariate biomarker under assignment to \eqn{Z=z}
#' measured at fixed time \eqn{t_0} after randomization, and \eqn{Y} is a binary clinical endpoint (\eqn{Y=1}, disease; \eqn{Y=0}, no disease) measured after \eqn{t_0}. The
#' estimator employs the generalized product kernel density estimation method of Hall, Racine, and Li (2004). The risks \eqn{P\{Y(z)=1|S(z)=s_1,X=x\}}, \eqn{z=0,1}, where
#' \eqn{X} is a vector of discrete baseline covariates, are estimated by fitting inverse probability-weighted logistic regression models.
#'
#' @param formula a formula object with the binary clinical endpoint on the left of the \code{~} operator. The first listed variable on the right must be the biomarker response
#' at \eqn{t0} and all variables that follow, if any, are discrete baseline covariates specified in all fitted models that condition on them. Interactions and transformations
#' of the baseline covariates are allowed. All terms in the formula must be evaluable in the data frame \code{data}.
#' @param bsm a character string specifying the variable name in \code{data} representing the baseline surrogate measure
#' @param tx a character string specifying the variable name in \code{data} representing the treatment group indicator
#' @param data a data frame with one row per randomized participant endpoint-free at \eqn{t_0} that contains at least the variables specified in \code{formula}, \code{bsm} and
#' \code{tx}. Biomarker values at baseline and \eqn{t_0} that are unavailable are represented as \code{NA}.
#' @param hinge a logical value (\code{FALSE} by default) indicating whether a hinge model (Fong et al., 2017) shall be used for modeling the effect of \eqn{S(z)} on the
#' clinical endpoint risk. A hinge model specifies that variability in \eqn{S(z)} below the hinge point does not associate with the clinical endpoint risk.
#' @param weights either a numeric vector of weights or a character string specifying the variable name in \code{data} representing weights applied to observations
#' in the phase 2 subset in order to make inference about the target population of all randomized participants endpoint-free at \eqn{t_0}. They reflect that
#' the case:control ratio in the phase 2 subset is different from that in the target population. The weights are passed on to GLMs in the estimation of the hinge point.
#' If \code{NULL} (default), a common weight is calculated separately for all cases and controls, pooling over treatment assignments, in the phase 2 subset.
#' @param biomarkerGrid a numeric vector of \eqn{S(1)} values at which the conditional endpoint risk in each study group is estimated. If \code{NULL} (default), a grid of
#' values spanning the range of observed values of the biomarker will be used.
#' @param saveFile a character string specifying the name of an \code{.RData} file storing the output list. If \code{NULL} (default), the output list will only be returned.
#' @param saveDir a character string specifying a path for the output directory. If \code{NULL} (default), the output list will only be returned; otherwise, if
#' \code{saveFile} is specified, the output list will also be saved as an \code{.RData} file in the specified directory.
#'
#' @return If \code{saveFile} and \code{saveDir} are both specified, the output list (named \code{oList}) is saved as an \code{.RData} file; otherwise it is returned only.
#' The output object is a list with the following components:
#' \itemize{
#' \item \code{biomarkerGrid}: a numeric vector of \eqn{S(1)} values at which the conditional endpoint risk is estimated in the components \code{plaRiskCurve} and
#' \code{txRiskCurve}
#' \item \code{plaRiskCurve}: estimates of \eqn{P\{Y(0)=1|S(1)=s_1\}} for \eqn{s_1} in \code{biomarkerGrid}
#' \item \code{txRiskCurve}: estimates of \eqn{P\{Y(1)=1|S(1)=s_1\}} for \eqn{s_1} in \code{biomarkerGrid}
#' \item \code{fOptBandwidths}: a \code{conbandwidth} object returned by the call of the function \code{npcdensbw} containing the optimal bandwidths, selected by likelihood
#' cross-validation, in the kernel estimation of the conditional density of \eqn{S(1)} given the baseline surrogate measure and any other specified baseline covariates
#' \item \code{gOptBandwidths}: a \code{conbandwidth} object returned by the call of the function \code{npcdensbw} or \code{npudensbw} containing the optimal bandwidths,
#' selected by likelihood cross-validation, in the kernel estimation of the conditional density of \eqn{S(0)} given any specified baseline covariates or the marginal density
#' of \eqn{S(0)} if no baseline covariates are specified in \code{formula}
#' \item \code{cpointP}: if \code{hinge=TRUE}, the estimate of the hinge point in the placebo group
#' \item \code{cpointT}: if \code{hinge=TRUE}, the estimate of the hinge point in the treatment group
#' }
#'
#' @examples
#' n <- 500
#' Z <- rep(0:1, each=n/2)
#' S <- mvrnorm(n, mu=c(2,2,3), Sigma=matrix(c(1,0.9,0.7,0.9,1,0.7,0.7,0.7,1), nrow=3))
#' p <- pnorm(drop(cbind(1,Z,(1-Z)*S[,2],Z*S[,3]) %*% c(-1.2,0.2,-0.02,-0.2)))
#' Y <- sapply(p, function(risk){ rbinom(1,1,risk) })
#' X <- rbinom(n,1,0.5)
#' # delete S(1) in placebo recipients
#' S[Z==0,3] <- NA
#' # delete S(0) in treatment recipients
#' S[Z==1,2] <- NA
#' # generate the indicator of being sampled into the phase 2 subset
#' phase2 <- rbinom(n,1,0.5)
#' # delete Sb, S(0) and S(1) in controls not included in the phase 2 subset
#' S[Y==0 & phase2==0,] <- c(NA,NA,NA)
#' # delete Sb in cases not included in the phase 2 subset
#' S[Y==1 & phase2==0,1] <- NA
#' data <- data.frame(X,Z,S[,1],ifelse(Z==0,S[,2],S[,3]),Y)
#' colnames(data) <- c("X","Z","Sb","S","Y")
#' grid <- with(data, seq(min(S, na.rm=TRUE), max(S, na.rm=TRUE), length.out=5))
#'
#' out <- riskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, biomarkerGrid=grid)
#' # alternatively, to save the .RData output file (no '<-' needed):
#' riskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, saveFile="out.RData", saveDir="./")
#'
#' @seealso \code{\link{bootRiskCurve}}
#' @export
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
  anyBaselineCovar <- length(formulaDecomp[[2]])>1

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
  if (anyBaselineCovar){
    fm.fbw <- as.formula(paste0("S ~ ",paste(c("Sb",formulaDecomp[[2]][-1]),collapse="+")))
    fm.gbw <- as.formula(paste0("S ~ ",paste(formulaDecomp[[2]][-1],collapse="+")))
    gbw <- npcdensbw(fm.gbw, data=dataP2correctRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  } else {
    fm.fbw <- S ~ Sb
    fm.gbw <- ~ S
    # marginal density
    gbw <- npudensbw(fm.gbw, data=dataP2correctRatio, ckertype="epanechnikov")
  }
  fbw <- npcdensbw(fm.fbw, data=dataT2correctRatio, cxkertype="epanechnikov", cykertype="epanechnikov")

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
    if (anyBaselineCovar){
      fm <- as.formula(paste0("Y ~ ",paste(formulaDecomp[[2]][-1],collapse="+")))
    } else {
      fm <- Y ~ 1
    }
    cpointP <- chngptm(formula.1=fm, formula.2=~S, data=subset(data2, Z==0 & !is.na(Y)), family="binomial", type="hinge", prob.weights=weights)$coefficients["chngpt"]
    cpointT <- chngptm(formula.1=Y ~ 1, formula.2=~S, data=subset(data2, Z==1 & !is.na(Y)), family="binomial", type="hinge", prob.weights=weights)$coefficients["chngpt"]
    # use their minimum as the hinge point in the below specified GLMs
    cpoint <- min(cpointP, cpointT)

    # biomarker S left-censored at 'cpoint' for use in the below specified GLMs
    data$Sc <- with(data, ifelse(S>cpoint, S-cpoint, 0))

    # re-extract the subset with phase 2 data
    data2 <- subset(data, !is.na(S))

    # IPW logistic regression model fitted to placebo recipients in the phase 2 subset accounting for two-phase sampling of S
    if (anyBaselineCovar){
      fm <- as.formula(paste0("Y ~ ",paste(c("Sc",formulaDecomp[[2]][-1]),collapse="+")))
    } else {
      fm <- Y ~ Sc
    }
    fit1 <- tps(fm, data=subset(data2, Z==0 & !is.na(Y)), nn0=nPControls, nn1=nPCases, group=rep(1, NROW(subset(data2, Z==0 & !is.na(Y)))), method="PL", cohort=TRUE)
  } else {
    # for passing on to the function 'risk'
    cpoint <- NULL

    # IPW logistic regression model fitted to placebo recipients in the phase 2 subset accounting for two-phase sampling of S
    if (anyBaselineCovar){
      fm <- as.formula(paste0("Y ~ ",paste(c("S",formulaDecomp[[2]][-1]),collapse="+")))
    } else {
      fm <- Y ~ S
    }
    fit1 <- tps(fm, data=subset(data2, Z==0 & !is.na(Y)), nn0=nPControls, nn1=nPCases, group=rep(1, NROW(subset(data2, Z==0 & !is.na(Y)))), method="PL", cohort=TRUE)
  }

  # kernel density estimator for f(s1|Sb=s0, X=x) using the treatment group in the phase 2 subset
  fhat <- npcdens(fbw)

  # kernel density estimator for g(s0|X=x) using the placebo group in the phase 2 subset
  if (anyBaselineCovar){
    ghat <- npcdens(gbw)
  } else {
    ghat <- npudens(gbw)
  }

  # a grid of values of S on which the estimated risk curves are returned
  if (is.null(biomarkerGrid)){
    biomarkerGrid <- seq(min(data2$S), max(data2$S), length.out=200)
  }

  # the first argument of 'risk' is a scalar
  if (anyBaselineCovar){
    fm <- as.formula(paste0("~",paste(formulaDecomp[[2]][-1],collapse="+")))
  } else {
    fm <- ~ 1
  }
  curves <- lapply(biomarkerGrid, function(s){ risk(s, fm, data, data2, fit1, fhat, ghat, cpoint) })
  plaRiskCurve <- sapply(curves, "[[", "plaRisk")
  txRiskCurve <- sapply(curves, "[[", "txRisk")

  # the output list
  oList <- list(biomarkerGrid=biomarkerGrid, plaRiskCurve=plaRiskCurve, txRiskCurve=txRiskCurve, fOptBandwidths=fbw, gOptBandwidths=gbw)
  if (hinge){
    oList$cpointP <- cpointP
    oList$cpointT <- cpointT
  }

  if (!is.null(saveFile) & !is.null(saveDir)){
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
  # patch
  #if (m[2] > 0) { # use data as environment
  #  mf[["formula"]] = eval(mf[[m[1]]], environment(mf[[m[2]]]))
  #} else { # use parent frame
  mf[["formula"]] = eval(mf[[m[1]]], parent.frame())
  #}
  # end of patch
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
