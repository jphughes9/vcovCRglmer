\name{vcovCR.glmerMod}
\alias{vcovCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster-robust variance-covariance matrix for an glmerMod object
}
\description{
\code{vcovCR} returns a robust ("sandwich") covariance matrix for glmer objects.
}
\usage{
vcovCR(obj,cluster,type="classic")
}

\arguments{
 \item{obj}{
Fitted glmer model for which to calculate the variance-covariance matrix
}
\item{cluster}{
Optional expression or vector indicating which observations belong to the same cluster. If not specified, will be set to clubSandwich:::get_outer_group(obj)
}
\item{type}{
Character string specifying which small-sample adjustment should be used, with available options "classic", "DF", "KC", "MD", or "FG(r)" where r is a number between 0 and 1 (note: FG(0) is equivalent to classic). See "Details" section for further information. Default is "classic".
}
}

\details{
Write the generalized linear mixed model as 
\deqn{E(Y | \gamma) = g^{-1}(X \beta + Z \gamma) = g^{-1}(\eta)}
\deqn{Var(Y | \gamma) = \Sigma}
where \eqn{Y} is \eqn{n \times 1}, \eqn{X} is \eqn{n \times p}, \eqn{Z} is \eqn{n \times q}, and \eqn{\gamma \sim N(0,R(\theta))}. \eqn{g} is a link function; only identity, log and logit link are accepted at this time.

Define
\deqn{\Delta = \left( \frac{\partial g^{-1}(\eta)}{\partial \eta} \right)_{\beta, \gamma} \mbox{a diagonal matrix}}
\deqn{P = \Delta^{-1} (Y - g^{-1}(\eta)) + X \beta + Z \gamma}
\deqn{V =  Z R(\theta) Z^{\prime} + \Delta^{-1} \Sigma \Delta^{-1}}
\deqn{e_i = P_i - X_i \beta }
\deqn{H_i = X_i  (X^{\prime} V^{-1} X)^{-1} X_i^{\prime} V_i^{-1}}
\deqn{Q_i = X_i^{\prime}  V_i^{-1} X_i (X^{\prime} V^{-1} X)^{-1}}
Note that \eqn{(X^{\prime} V^{-1} X)^{-1}} is the usual model-based variance.

The robust variance of \eqn{\hat{\beta}} is given by the expression\eqn{^1}
\deqn{
c \times (X^{\prime} V^{-1} X)^{-1} \left( \sum_g A_g X_g^{\prime} V_g^{-1} F_g^{\prime} e_g e_g^{\prime} F_g V_g^{-1} X_g A_g \right) (X^{\prime} V^{-1} X)^{-1}}
where the sum is over the independent groups (clusters) and all terms are evaluated at the estimates of \eqn{\hat{\beta}, \hat{\theta}, \mbox{and } \hat{\Sigma}}. See the table for values for c, F and A. Note that for the classical sandwich estimator, c=1 and F = A = I. 

\tabular{lccc}{ 
 \tab c \tab \eqn{A_g} \tab \eqn{F_g} \cr 
 classical \tab 1 \tab I \tab I \cr
DF \tab \eqn{\frac{n}{n-p}} \tab I \tab I \cr
Kauermann and Carroll (2001) \tab 1 \tab I \tab \eqn{(I - H_g^{\prime})^{-1/2}} \cr
Mancl and DeRouen (2001) \tab 1 \tab I \tab \eqn{(I - H_g^{\prime})^{-1}} \cr
Fay and Graubard (2001) \tab 1 \tab \eqn{\mbox{diag} \left\{( 1 - \min(r, [Q_g]_{jj}))^{-1/2} \right\}} \tab I \cr
}
Note that \eqn{0 \leq r < 1}. If type = "FG" with no value for r then a value of R=0.75 is assumed. Also, if type = "DF" is selected but \eqn{n - np \le 0} then a warning is issued and type defaults to "classic".

\eqn{^1} {When type = "KC" is selected then the approximation given by Li et al (2020) is used in place of this expression.}
}

\value{
A matrix of the estimated variance of and covariances between the regression fixed effect coefficient estimates.
}

\references{

Fay MP and Graubard BI (2001) Small-Sample Adjustments for Wald-Type Tests Using Sandwich Estimators. Biometrics 57:1198–1206.

Kauermann G and Carroll RJ (2001) A Note on the Efficiency of Sandwich Covariance Estimation. Journal of the American Statistical Association 96:1387–1396.

Gallis JA, Li F and Turner EL (2020) xtgeebcv: A command for bias-corrected
sandwich variance estimation for GEE analyses of cluster randomized trials. Stata Journal 20(2) 363:381.

Mancl LA and DeRouen TA (2001) A Covariance Estimator for GEE with Improved Small-Sample Properties. Biometrics 57:126–134

SAS GLIMMIX documentation (https://documentation.sas.com/doc/en/statug/15.2/statug_glimmix_syntax01.htm)
}
\author{
James P Hughes
}
%%\note{
%%  ~~further notes~~
%%}

\examples{
rslt = glmer(ct ~ ftime + treat + (1 | fcluster), family=binomial, data=ctdata)
vcovCR(rslt)
vcovCR(rslt,type="FG(0)") # gives identical result
}
\keyword{ robust variance }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line
