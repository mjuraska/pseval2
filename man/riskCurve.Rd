% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssmooth.R
\name{riskCurve}
\alias{riskCurve}
\title{Estimation of Conditional Clinical Endpoint Risk under Placebo and Treatment Given Biomarker Response to Treatment in a Baseline Surrogate Measure Three-Phase Sampling Design}
\usage{
riskCurve(
  formula,
  bsm,
  tx,
  data,
  pstype = c("continuous", "ordered"),
  bsmtype = c("continuous", "ordered"),
  bwtype = c("fixed", "generalized_nn", "adaptive_nn"),
  hinge = FALSE,
  weights = NULL,
  psGrid = NULL,
  saveFile = NULL,
  saveDir = NULL
)
}
\arguments{
\item{formula}{a formula object with the binary clinical endpoint on the left of the \code{~} operator. The first listed variable on the right must be the biomarker response
at \eqn{t0} and all variables that follow, if any, are discrete baseline covariates specified in all fitted models that condition on them. Interactions and transformations
of the baseline covariates are allowed. All terms in the formula must be evaluable in the data frame \code{data}.}

\item{bsm}{a character string specifying the variable name in \code{data} representing the baseline surrogate measure}

\item{tx}{a character string specifying the variable name in \code{data} representing the treatment group indicator}

\item{data}{a data frame with one row per randomized participant endpoint-free at \eqn{t_0} that contains at least the variables specified in \code{formula}, \code{bsm} and
\code{tx}. Values of \code{bsm} and the biomarker at \eqn{t_0} that are unavailable are represented as \code{NA}.}

\item{pstype}{a character string specifying whether the biomarker response shall be treated as a \code{continuous} (default) or \code{ordered} categorical variable in the
kernel density/probability estimation}

\item{bsmtype}{a character string specifying whether the baseline surrogate measure shall be treated as a \code{continuous} (default) or \code{ordered} categorical variable in the
kernel density/probability estimation}

\item{bwtype}{a character string specifying the bandwidth type for continuous variables in the kernel density estimation. The options are \code{fixed} (default) for fixed
bandwidths, \code{generalized_nn} for generalized nearest neighbors, and \code{adaptive_nn} for adaptive nearest neighbors. As noted in the documentation of the function
\code{npcdensbw} in the \code{np} package: "Adaptive nearest-neighbor bandwidths change with each sample realization in the set when estimating the density at the point \eqn{x}.
Generalized nearest-neighbor bandwidths change with the point at which the density is estimated, \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}."}

\item{hinge}{a logical value (\code{FALSE} by default) indicating whether a hinge model (Fong et al., 2017) shall be used for modeling the effect of \eqn{S(z)} on the
clinical endpoint risk. A hinge model specifies that variability in \eqn{S(z)} below the hinge point does not associate with the clinical endpoint risk.}

\item{weights}{either a numeric vector of weights or a character string specifying the variable name in \code{data} representing weights applied to observations
in the phase 2 subset in order to make inference about the target population of all randomized participants endpoint-free at \eqn{t_0}. The weights reflect that
the case:control ratio in the phase 2 subset is different from that in the target population and are passed on to GLMs in the estimation of the hinge point.
If \code{NULL} (default), weights for cases and controls are calculated separately in each study group.}

\item{psGrid}{a numeric vector of \eqn{S(1)} values at which the conditional clinical endpoint risk in each study group is estimated. If \code{NULL} (default),
a grid of values spanning the range of observed values of the biomarker will be used.}

\item{saveFile}{a character string specifying the name of an \code{.RData} file storing the output list. If \code{NULL} (default), the output list will only be returned.}

\item{saveDir}{a character string specifying a path for the output directory. If \code{NULL} (default), the output list will only be returned; otherwise, if
\code{saveFile} is specified, the output list will also be saved as an \code{.RData} file in the specified directory.}
}
\value{
If \code{saveFile} and \code{saveDir} are both specified, the output list (named \code{oList}) is saved as an \code{.RData} file; otherwise it is returned only.
The output object (of class \code{riskCurve}) is a list with the following components:
\itemize{
\item \code{psGrid}: a numeric vector of \eqn{S(1)} values at which the conditional clinical endpoint risk is estimated in the components \code{plaRiskCurve} and
\code{txRiskCurve}
\item \code{plaRiskCurve}: a numeric vector of estimates of \eqn{P\{Y(0)=1|S(1)=s_1\}} for \eqn{s_1} in \code{psGrid}
\item \code{txRiskCurve}: a numeric vector of estimates of \eqn{P\{Y(1)=1|S(1)=s_1\}} for \eqn{s_1} in \code{psGrid}
\item \code{fOptBandwidths}: a \code{conbandwidth} object returned by the call of the function \code{npcdensbw} containing the optimal bandwidths, selected by likelihood
cross-validation, in the kernel estimation of the conditional density of \eqn{S(1)} given the baseline surrogate measure and any other specified baseline covariates
\item \code{gOptBandwidths}: a \code{conbandwidth} object returned by the call of the function \code{npcdensbw} or \code{npudensbw} containing the optimal bandwidths,
selected by likelihood cross-validation, in the kernel estimation of the conditional density of \eqn{S(0)} given any specified baseline covariates or the marginal density
of \eqn{S(0)} if no baseline covariates are specified in \code{formula}
\item \code{cpointP}: if \code{hinge=TRUE}, the estimate of the hinge point in the placebo group
\item \code{cpointT}: if \code{hinge=TRUE}, the estimate of the hinge point in the treatment group
}
}
\description{
Estimates \eqn{P\{Y(z)=1|S(1)=s_1\}}, \eqn{z=0,1}, on a grid of \eqn{s_1} values following the estimation method of Juraska, Huang, and Gilbert (2018), where \eqn{Z} is the
treatment group indicator (\eqn{Z=1}, treatment; \eqn{Z=0}, placebo), \eqn{S(z)} is a continuous or ordered categorical univariate biomarker under assignment to \eqn{Z=z}
measured at fixed time \eqn{t_0} after randomization, and \eqn{Y} is a binary clinical endpoint (\eqn{Y=1}, disease; \eqn{Y=0}, no disease) measured after \eqn{t_0}. The
estimator employs the generalized product kernel density/probability estimation method of Hall, Racine, and Li (2004) implemented in the \code{np} package. The risks
\eqn{P\{Y(z)=1|S(z)=s_1,X=x\}}, \eqn{z=0,1}, where \eqn{X} is a vector of discrete baseline covariates, are estimated by fitting inverse probability-weighted logistic regression
models using the \code{osDesign} package.
}
\examples{
n <- 500
Z <- rep(0:1, each=n/2)
S <- MASS::mvrnorm(n, mu=c(2,2,3), Sigma=matrix(c(1,0.9,0.7,0.9,1,0.7,0.7,0.7,1), nrow=3))
p <- pnorm(drop(cbind(1,Z,(1-Z)*S[,2],Z*S[,3]) \%*\% c(-1.2,0.2,-0.02,-0.2)))
Y <- sapply(p, function(risk){ rbinom(1,1,risk) })
X <- rbinom(n,1,0.5)
# delete S(1) in placebo recipients
S[Z==0,3] <- NA
# delete S(0) in treatment recipients
S[Z==1,2] <- NA
# generate the indicator of being sampled into the phase 2 subset
phase2 <- rbinom(n,1,0.4)
# delete Sb, S(0) and S(1) in controls not included in the phase 2 subset
S[Y==0 & phase2==0,] <- c(NA,NA,NA)
# delete Sb in cases not included in the phase 2 subset
S[Y==1 & phase2==0,1] <- NA
data <- data.frame(X,Z,S[,1],ifelse(Z==0,S[,2],S[,3]),Y)
colnames(data) <- c("X","Z","Sb","S","Y")
qS <- quantile(data$S, probs=c(0.05,0.95), na.rm=TRUE)
grid <- seq(qS[1], qS[2], length.out=3)

out <- riskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, psGrid=grid)
\donttest{
# alternatively, to save the .RData output file (no '<-' needed):
riskCurve(formula=Y ~ S + factor(X), bsm="Sb", tx="Z", data=data, saveFile="out.RData",
          saveDir="./")
}

}
\references{
Fong, Y., Huang, Y., Gilbert, P. B., and Permar, S. R. (2017), chngpt: threshold regression model estimation and inference, \emph{BMC Bioinformatics}, 18.

Hall, P., Racine, J., and Li, Q. (2004), Cross-validation and the estimation of conditional probability densities, \emph{JASA} 99(468), 1015-1026.

Juraska, M., Huang, Y., and Gilbert, P. B. (2020), Inference on treatment effect modification by biomarker response in a three-phase sampling design, Biostatistics, 21(3): 545-560, \url{https://doi.org/10.1093/biostatistics/kxy074}.
}
\seealso{
\code{\link{bootRiskCurve}}, \code{\link{summary.riskCurve}} and \code{\link{plotMCEPcurve}}
}
