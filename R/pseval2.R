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

#' Bootstrap Estimation of Conditional Clinical Endpoint Risk under Placebo and Treatment Given Biomarker Response to Treatment in a Baseline Surrogate Measure
#' Three-Phase Sampling Design
#'
#' Estimates \eqn{P\{Y(z)=1|S(1)=s_1\}}, \eqn{z=0,1}, on a grid of \eqn{s_1} values in bootstrap resamples (see \code{\link{riskCurve}} for notation introduction). Cases
#' (\eqn{Y=1}) and controls (\eqn{Y=0}) are sampled separately yielding a fixed number of cases and controls in each bootstrap sample. Consequentially, the number of controls
#' with available phase 2 data varies across bootstrap samples.
#'
#' @param formula a formula object with the binary clinical endpoint on the left of the \code{~} operator. The first listed variable on the right must be the biomarker response
#' at \eqn{t0} and all variables that follow, if any, are discrete baseline covariates specified in all fitted models that condition on them. Interactions and transformations
#' of the baseline covariates are allowed. All terms in the formula must be evaluable in the data frame \code{data}.
#' @param bsm a character string specifying the variable name in \code{data} representing the baseline surrogate measure
#' @param tx a character string specifying the variable name in \code{data} representing the treatment group indicator
#' @param data a data frame with one row per randomized participant endpoint-free at \eqn{t_0} that contains at least the variables specified in \code{formula}, \code{bsm} and
#' \code{tx}. Values of \code{bsm} and the biomarker at \eqn{t_0} that are unavailable are represented as \code{NA}.
#' @param hinge a logical value (\code{FALSE} by default) indicating whether a hinge model (Fong et al., 2017) shall be used for modeling the effect of \eqn{S(z)} on the
#' clinical endpoint risk. A hinge model specifies that variability in \eqn{S(z)} below the hinge point does not associate with the clinical endpoint risk. The hinge point
#' is reestimated in each bootstrap sample.
#' @param weights either a numeric vector of weights or a character string specifying the variable name in \code{data} representing weights applied to observations
#' in the phase 2 subset in order to make inference about the target population of all randomized participants endpoint-free at \eqn{t_0}. The weights reflect that
#' the case:control ratio in the phase 2 subset is different from that in the target population and are passed on to GLMs in the estimation of the hinge point.
#' If \code{NULL} (default and recommended), weights for cases and controls are recalculated separately in each study group \emph{within each bootstrap sample}; otherwise the
#' same specified vector of weights is used in each bootstrap sample.
#' @param biomarkerGrid a numeric vector of \eqn{S(1)} values at which the conditional clinical endpoint risk in each study group is estimated. If \code{NULL} (default),
#' a grid of values spanning the range of observed values of the biomarker will be used.
#' @param iter the number of bootstrap iterations
#' @param seed a seed of the random number generator supplied to \code{set.seed} for reproducibility
#' @param saveFile a character string specifying the name of an \code{.RData} file storing the output list. If \code{NULL} (default), the output list will only be returned.
#' @param saveDir a character string specifying a path for the output directory. If \code{NULL} (default), the output list will only be returned; otherwise, if
#' \code{saveFile} is specified, the output list will also be saved as an \code{.RData} file in the specified directory.
#'
#' @return If \code{saveFile} and \code{saveDir} are both specified, the output list (named \code{oList}) is saved as an \code{.RData} file; otherwise it is returned only.
#' The output object is a list with the following components:
#' \itemize{
#' \item \code{biomarkerGrid}: a numeric vector of \eqn{S(1)} values at which the conditional clinical endpoint risk is estimated in the components \code{plaRiskCurveBoot} and
#' \code{txRiskCurveBoot}
#' \item \code{plaRiskCurveBoot}: a \code{length(biomarkerGrid)}-by-\code{iter} matrix of estimates of \eqn{P\{Y(0)=1|S(1)=s_1\}} for \eqn{s_1} in \code{biomarkerGrid},
#' with columns representing bootstrap samples
#' \item \code{txRiskCurveBoot}: a \code{length(biomarkerGrid)}-by-\code{iter} matrix of estimates of \eqn{P\{Y(1)=1|S(1)=s_1\}} for \eqn{s_1} in \code{biomarkerGrid},
#' with columns representing bootstrap samples
#' \item \code{cpointPboot}: if \code{hinge=TRUE}, a numeric vector of estimates of the hinge point in the placebo group in each bootstrap sample
#' \item \code{cpointTboot}: if \code{hinge=TRUE}, a numeric vector of estimates of the hinge point in the treatment group in each bootstrap sample
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
#' phase2 <- rbinom(n,1,0.4)
#' # delete Sb, S(0) and S(1) in controls not included in the phase 2 subset
#' S[Y==0 & phase2==0,] <- c(NA,NA,NA)
#' # delete Sb in cases not included in the phase 2 subset
#' S[Y==1 & phase2==0,1] <- NA
#' data <- data.frame(X,Z,S[,1],ifelse(Z==0,S[,2],S[,3]),Y)
#' colnames(data) <- c("X","Z","Sb","S","Y")
#' qS <- quantile(data$S, probs=c(0.05,0.95), na.rm=TRUE)
#' grid <- seq(qS[1], qS[2], length.out=5)
#'
#' out <- bootRiskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, biomarkerGrid=grid, iter=1, seed=10)
#' # alternatively, to save the .RData output file (no '<-' needed):
#' bootRiskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, iter=1, seed=10, saveFile="out.RData", saveDir="./")
#'
#' @seealso \code{\link{riskCurve}}
#' @export
bootRiskCurve <- function(formula, bsm, tx, data, hinge=FALSE, weights=NULL, biomarkerGrid=NULL, iter, seed=NULL, saveFile=NULL, saveDir=NULL){
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
        bdata2$weights <- NA
        bdata2$weights <- ifelse(bdata2$Z==0 & bdata2$Y==0, nPControls/nPControls2, bdata2$weights)
        bdata2$weights <- ifelse(bdata2$Z==1 & bdata2$Y==0, nTControls/nTControls2, bdata2$weights)
        bdata2$weights <- ifelse(bdata2$Z==0 & bdata2$Y==1, nPCases/nPCases2, bdata2$weights)
        bdata2$weights <- ifelse(bdata2$Z==1 & bdata2$Y==1, nTCases/nTCases2, bdata2$weights)
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
  plaRiskCurveBoot <- sapply(bRiskCurveList,"[[","plaRiskCurve")
  txRiskCurveBoot <- sapply(bRiskCurveList,"[[","txRiskCurve")

  # the output list
  bList <- list(biomarkerGrid=biomarkerGrid, plaRiskCurveBoot=plaRiskCurveBoot, txRiskCurveBoot=txRiskCurveBoot)
  if (hinge){
    bList$cpointPboot <- sapply(bRiskCurveList,"[[","cpointP")
    bList$cpointTboot <- sapply(bRiskCurveList,"[[","cpointT")
  }

  if (!is.null(saveFile)){
    save(bList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n",file.path(saveDir, saveFile),"\n")
  }

  return(invisible(bList))
}

#' Estimation of Conditional Clinical Endpoint Risk under Placebo and Treatment Given Biomarker Response to Treatment in a Baseline Surrogate Measure Three-Phase Sampling Design
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
#' \code{tx}. Values of \code{bsm} and the biomarker at \eqn{t_0} that are unavailable are represented as \code{NA}.
#' @param hinge a logical value (\code{FALSE} by default) indicating whether a hinge model (Fong et al., 2017) shall be used for modeling the effect of \eqn{S(z)} on the
#' clinical endpoint risk. A hinge model specifies that variability in \eqn{S(z)} below the hinge point does not associate with the clinical endpoint risk.
#' @param weights either a numeric vector of weights or a character string specifying the variable name in \code{data} representing weights applied to observations
#' in the phase 2 subset in order to make inference about the target population of all randomized participants endpoint-free at \eqn{t_0}. The weights reflect that
#' the case:control ratio in the phase 2 subset is different from that in the target population and are passed on to GLMs in the estimation of the hinge point.
#' If \code{NULL} (default), weights for cases and controls are calculated separately in each study group.
#' @param biomarkerGrid a numeric vector of \eqn{S(1)} values at which the conditional clinical endpoint risk in each study group is estimated. If \code{NULL} (default),
#' a grid of values spanning the range of observed values of the biomarker will be used.
#' @param saveFile a character string specifying the name of an \code{.RData} file storing the output list. If \code{NULL} (default), the output list will only be returned.
#' @param saveDir a character string specifying a path for the output directory. If \code{NULL} (default), the output list will only be returned; otherwise, if
#' \code{saveFile} is specified, the output list will also be saved as an \code{.RData} file in the specified directory.
#'
#' @return If \code{saveFile} and \code{saveDir} are both specified, the output list (named \code{oList}) is saved as an \code{.RData} file; otherwise it is returned only.
#' The output object (of class \code{"riskCurve"}) is a list with the following components:
#' \itemize{
#' \item \code{biomarkerGrid}: a numeric vector of \eqn{S(1)} values at which the conditional clinical endpoint risk is estimated in the components \code{plaRiskCurve} and
#' \code{txRiskCurve}
#' \item \code{plaRiskCurve}: a numeric vector of estimates of \eqn{P\{Y(0)=1|S(1)=s_1\}} for \eqn{s_1} in \code{biomarkerGrid}
#' \item \code{txRiskCurve}: a numeric vector of estimates of \eqn{P\{Y(1)=1|S(1)=s_1\}} for \eqn{s_1} in \code{biomarkerGrid}
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
#' phase2 <- rbinom(n,1,0.4)
#' # delete Sb, S(0) and S(1) in controls not included in the phase 2 subset
#' S[Y==0 & phase2==0,] <- c(NA,NA,NA)
#' # delete Sb in cases not included in the phase 2 subset
#' S[Y==1 & phase2==0,1] <- NA
#' data <- data.frame(X,Z,S[,1],ifelse(Z==0,S[,2],S[,3]),Y)
#' colnames(data) <- c("X","Z","Sb","S","Y")
#' qS <- quantile(data$S, probs=c(0.05,0.95), na.rm=TRUE)
#' grid <- seq(qS[1], qS[2], length.out=5)
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

  # convert any baseline covariates into factors
  if (anyBaselineCovar){
    # 'vars' is of length >= 3
    vars <- all.vars(formula)[-(1:2)]
    for (i in 1:length(vars)){
      data[,vars[i]] <- as.factor(data[,vars[i]])
    }
  }

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
      data2$weights <- NA
      data2$weights <- ifelse(data2$Z==0 & data2$Y==0, nPControls/nPControls2, data2$weights)
      data2$weights <- ifelse(data2$Z==1 & data2$Y==0, nTControls/nTControls2, data2$weights)
      data2$weights <- ifelse(data2$Z==0 & data2$Y==1, nPCases/nPCases2, data2$weights)
      data2$weights <- ifelse(data2$Z==1 & data2$Y==1, nTCases/nTCases2, data2$weights)
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
  class(oList) <- "riskCurve"

  if (!is.null(saveFile) & !is.null(saveDir)){
    save(oList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n",file.path(saveDir, saveFile),"\n")
  }

  return(invisible(oList))
}

# 'object' is the output list from either 'riskCurve' or 'bootRiskCurve'
contrastRiskCurve <- function(object, contrast){
  if (contrast=="te"){ return(1 - object$txRisk/object$plaRisk) }
  if (contrast=="rr"){ return(object$txRisk/object$plaRisk) }
  if (contrast=="logrr"){ return(log(object$txRisk/object$plaRisk)) }
  if (contrast=="rd"){ return(object$plaRisk - object$txRisk) }
}

# 'object' is the output list from either 'riskCurve' or 'bootRiskCurve'
tContrastRiskCurve <- function(object, contrast){
  if (contrast %in% c("te","rr","logrr")){ return(log(object$txRisk/object$plaRisk)) }
  if (contrast=="rd"){ return(object$plaRisk - object$txRisk) }
}

# 'x' is a numeric vector
invtContrastRiskCurve <- function(x, contrast){
  if (contrast=="te"){ return(1 - exp(x)) }
  if (contrast=="rr"){ return(exp(x)) }
  if (contrast=="logrr"){ return(x) }
  if (contrast=="rd"){ return(x) }
}

#' Summary of Point and Interval Estimation of a Marginal Causal Effect Predictiveness Curve
#'
#' Summarizes point estimates and pointwise and simultaneous Wald-type bootstrap confidence intervals for a specified marginal causal effect predictiveness (mCEP) curve (see,
#' e.g., Juraska, Huang, and Gilbert (2018) for the definition).
#'
#' @param object an object of class \code{"riskCurve"}, typically returned by \code{\link{riskCurve}}
#' @param boot an object returned by \code{\link{bootRiskCurve}}. If \code{NULL} (default), only point estimates are reported.
#' @param contrast a character string specifying the mCEP curve. It must be one of \code{"te"} (treatment efficacy), \code{"rr"} (relative risk), \code{"logrr"} (log relative risk), and \code{"rd"} (risk
#' difference [placebo minus treatment]).
#' @param confLevel the confidence level of pointwise and simultaneous confidence intervals
#' @param \dots for other methods
#'
#' @return The output object is a data frame containing point and possibly interval estimates of the specified mCEP curve.
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
#' phase2 <- rbinom(n,1,0.4)
#' # delete Sb, S(0) and S(1) in controls not included in the phase 2 subset
#' S[Y==0 & phase2==0,] <- c(NA,NA,NA)
#' # delete Sb in cases not included in the phase 2 subset
#' S[Y==1 & phase2==0,1] <- NA
#' data <- data.frame(X,Z,S[,1],ifelse(Z==0,S[,2],S[,3]),Y)
#' colnames(data) <- c("X","Z","Sb","S","Y")
#' qS <- quantile(data$S, probs=c(0.05,0.95), na.rm=TRUE)
#' grid <- seq(qS[1], qS[2], length.out=5)
#'
#' out <- riskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, biomarkerGrid=grid)
#' boot <- bootRiskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, biomarkerGrid=grid, iter=2, seed=10)
#' summary(out, boot, contrast="te")
#'
#' @seealso \code{\link{riskCurve}} and \code{\link{bootRiskCurve}}
#' @export
summary.riskCurve <- function(object, boot=NULL, contrast=c("te", "rr", "logrr", "rd"), confLevel=0.95,...){
  contrast <- match.arg(contrast)

  # point estimates of mCEP(s1)
  MCEP <- contrastRiskCurve(object, contrast)
  out <- data.frame(object$biomarkerGrid, MCEP)
  colnames(out) <- c("biomarkerGrid", contrast)

  # interval estimates of mCEP(s1)
  if(!is.null(boot)){
    # transformed MCEP curve
    tMCEP <- tContrastRiskCurve(object, contrast)
    # transformed bootstrapped MCEP curves
    # assuming the matrices have the same dimensions
    tbMCEP <- tContrastRiskCurve(boot, contrast)

    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)

    # pointwise confidence bounds for MCEP(s1)
    ptLB.MCEP <- invtContrastRiskCurve(tMCEP - qnorm(1-(1-confLevel)/2) * bSE, contrast=contrast)
    ptUB.MCEP <- invtContrastRiskCurve(tMCEP + qnorm(1-(1-confLevel)/2) * bSE, contrast=contrast)

    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=confLevel, na.rm=TRUE)

    smLB.MCEP <- invtContrastRiskCurve(tMCEP - qSupAbsZ * bSE, contrast=contrast)
    smUB.MCEP <- invtContrastRiskCurve(tMCEP + qSupAbsZ * bSE, contrast=contrast)

    if (contrast=="te"){
      tmp <- ptUB.MCEP
      ptUB.MCEP <- ptLB.MCEP
      ptLB.MCEP <- tmp

      tmp <- smUB.MCEP
      smUB.MCEP <- smLB.MCEP
      smLB.MCEP <- tmp
    }

    out$ptLB <- ptLB.MCEP
    out$ptUB <- ptUB.MCEP
    out$smLB <- smLB.MCEP
    out$smUB <- smUB.MCEP
  }

  return(out)
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
testConstancy <- function(object, boot, contrast=c("te", "rr", "logrr", "rd"), null=c("H01", "H02"), overallPlaRisk=NULL, overallTxRisk=NULL, MCEPconstantH02=NULL, limS1=NULL){
  contrast <- match.arg(contrast)
  null <- match.arg(null)

  if (null=="H01"){
    if (is.null(overallPlaRisk) | is.null(overallTxRisk)){ stop("'overallPlaRisk' and 'overallTxRisk' must be specified for the test of H01.") }
  }

  if (null=="H02"){
    if (is.null(MCEPconstantH02)){ stop("'MCEPconstantH02' must be specified for the test of H02.") }

    # trim the risk curves if 'limS1' is specified
    if (!is.null(limS1)){
      object$plaRiskCurve <- object$plaRiskCurve[object$biomarkerGrid>=limS1[1] & object$biomarkerGrid<=limS1[2]]
      object$txRiskCurve <- object$txRiskCurve[object$biomarkerGrid>=limS1[1] & object$biomarkerGrid<=limS1[2]]
      boot$plaRiskCurveBoot <- boot$plaRiskCurveBoot[boot$biomarkerGrid>=limS1[1] & boot$biomarkerGrid<=limS1[2],]
      boot$txRiskCurveBoot <- boot$txRiskCurveBoot[boot$biomarkerGrid>=limS1[1] & boot$biomarkerGrid<=limS1[2],]
    }
  }

  # transformed estimated MCEP curve
  tMCEP <- tContrastRiskCurve(object, contrast)
  # transformed bootstrapped MCEP curves
  tbMCEP <- tContrastRiskCurve(boot, contrast)

  # bootstrap SE of tMCEP estimates
  bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)

  # calculate the supremum statistic for each bootstrap sample
  supAbsZ <- NULL
  for (j in 1:NCOL(tbMCEP)){
    Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
  }

  if (null=="H01"){
    tCEest <- ifelse(contrast=="rd", overallPlaRisk - overallTxRisk, log(overallTxRisk/overallPlaRisk))
    testStat <- max(abs(tMCEP-tCEest)/bSE, na.rm=TRUE)
    return(mean(supAbsZ > testStat))
  }

  if (null=="H02"){
    tMCEPconstantH02 <- switch(contrast, te=log(1-MCEPconstantH02), rr=log(MCEPconstantH02), logrr=MCEPconstantH02, rd=MCEPconstantH02)
    testStat <- max(abs(tMCEP-tMCEPconstantNull)/bSE, na.rm=TRUE)
    return(mean(supAbsZ > testStat))
  }
}

#' Plotting of the Estimated Marginal Causal Effect Predictiveness Curve
#'
#' Plots point estimates and, if available, pointwise and simultaneous Wald-type bootstrap confidence intervals for the specified marginal causal effect predictiveness (mCEP)
#' curve.
#'
#' @param object an object returned by \code{\link{summary.riskCurve}}
#' @param confLevel the confidence level (0.95 by default) of pointwise and simultaneous confidence intervals
#' @param hingePoint the hinge point estimate (\code{NULL} by default)
#' @param xLab a character string specifying the x-axis label (\code{NULL} by default)
#'
#' @return None. The function is called solely for plot generation.
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
#' phase2 <- rbinom(n,1,0.4)
#' # delete Sb, S(0) and S(1) in controls not included in the phase 2 subset
#' S[Y==0 & phase2==0,] <- c(NA,NA,NA)
#' # delete Sb in cases not included in the phase 2 subset
#' S[Y==1 & phase2==0,1] <- NA
#' data <- data.frame(X,Z,S[,1],ifelse(Z==0,S[,2],S[,3]),Y)
#' colnames(data) <- c("X","Z","Sb","S","Y")
#' qS <- quantile(data$S, probs=c(0.05,0.95), na.rm=TRUE)
#' grid <- seq(qS[1], qS[2], length.out=5)
#'
#' out <- riskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, biomarkerGrid=grid)
#' boot <- bootRiskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, biomarkerGrid=grid, iter=2, seed=10)
#' sout <- summary(out, boot, contrast="te")
#' plotMCEPcurve(sout)
#'
#' @seealso \code{\link{riskCurve}}, \code{\link{bootRiskCurve}} and \code{\link{summary.riskCurve}}
#' @export
plotMCEPcurve <- function(object, confLevel=0.95, hingePoint=NULL, xLab=NULL){
  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2

  yLim <- range(object[,-1], na.rm=TRUE)
  if (is.null(xLab)){ xLab <- expression(paste("Biomarker Response at ",t[0])) }
  yLab <- switch(colnames(object)[2], te="Treatment Efficacy", rr="Relative Risk", logrr="Log Relative Risk", rd="Risk Difference (Pla - Tx)")

  par(mar=c(5,5,1,1), cex.lab=cexLab, cex.axis=cexAxis, las=1)

  plot(object[,1], object[,2], type="n", ylim=c(yLim[1]-0.1*(yLim[2]-yLim[1]), yLim[2]), xlab=xLab, ylab=yLab)
  abline(h=ifelse(colnames(object)[2]=="rr", 1, 0), lty="dotted", lwd=2, col="gray50")

  lines(object[,1], object[,2], lwd=3.5)

  # if interval estimates are available in the data frame
  if (NCOL(object) > 2){
    lines(object[,1], object$ptLB, lty="dashed", lwd=3)
    lines(object[,1], object$ptUB, lty="dashed", lwd=3)
    lines(object[,1], object$smLB, lty="dotdash", lwd=3)
    lines(object[,1], object$smUB, lty="dotdash", lwd=3)
  }

  legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c(paste0("Pointwise ",confLevel*100,"% CI"), paste0("Simultaneous ",confLevel*100,"% CI")),
         cex=cexLegend, bty="n")
  if (!is.null(hingePoint)){ legend("bottomright", paste0("Hinge Point = ", round(hingePoint,2),"  "), cex=cexLegend, bty="n") }
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
