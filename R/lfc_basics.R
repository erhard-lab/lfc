
#' Logit transformation to obtain the log fold change representation from the
#' proportion representation.
#'
#' @seealso ltop
#' @family Effect size representations
#' @param p Effect size in proportion representation
#' @return The log2 fold change representation of the effect size
#' @examples
#' ptol(0.5)
#' ptol(2/3)
#' @export
ptol <- function(p) log2(p/(1-p))


#' Inverse logit transformation to obtain proportion representation from the
#' log fold change representation.
#'
#' @seealso ptol
#' @family Effect size representations
#' @param l Effect size in log2 fold change representation
#' @return The proportion representation of the effect size
#' @examples
#' ptol(0)
#' ptol(1)
#' @export
ltop <- function(l) 2^l/(1+2^l)


#' Density, distribution function, quantile function and random
#' generation for the log2 fold change  distribution with parameters
#' ‘a’ and ‘b’ (corresponding to (pseudo-)counts incremented by 1).
#'
#' @title The log2 fold change distribution
#' @describeIn dlfc Density function
#' @family Log2 fold change distribution
#' @param x,q vector of quantiles
#' @param a non-negative parameter
#' @param b non-negative parameter
#' @param log,log.p if TRUE, probabilities p are given as log(p)
#' @examples
#' x <- seq (-5,5,by=0.01)
#' plot(x,dlfc(x,1,1))
#' @export
dlfc <- function(x,a,b,log=FALSE) {
    r <- a*x*log(2)+log(log(2))-lbeta(a,b)-(a+b)*log(1+2^x)
    if (!log) r <- exp(r)
    r
}

#' @describeIn dlfc Distribution function
#' @family Log2 fold change distribution
#' @param lower.tail if TRUE (default), probabilities are P[X <= x],
#'     otherwise, P[X > x].
#' @export
#' @importFrom stats pbeta
plfc <- function(q,a,b,lower.tail = TRUE,log.p=FALSE) pbeta(ltop(q),shape1=a,shape2=b,lower.tail=lower.tail,log.p = log.p)

#' @describeIn dlfc Quantile function
#' @family Log2 fold change distribution
#' @param p vector of probabilities
#' @export
#' @importFrom stats qbeta
qlfc <- function(p,a,b,lower.tail = TRUE,log.p=FALSE) ptol(qbeta(p,shape1=a,shape2=b,lower.tail=lower.tail,log.p = log.p))

#' @describeIn dlfc random generation
#' @family Log2 fold change distribution
#' @param n number of observations
#' @export
#' @importFrom stats rbeta
rlfc <- function(n,a,b) ptol(rbeta(n,shape1=a,shape2=b))


