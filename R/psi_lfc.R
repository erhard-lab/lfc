
#' Computes the prior parameters (i.e. pseudocounts incremented by 1) for
#' the log2 fold change distribution
#'
#' @param A Vector of counts from condition A
#' @param B Vector of counts from condition B
#' @seealso PsiLFC
#' @return A vector of length 2 containing the two parameters
#' @export
#' @importFrom stats median pnorm optim var quantile
EmpiricalBayesPrior=function(A,B) {

    u=A>0 | B>0
    A=A[u]
    B=B[u]

    x=median(log(A)-log(B))
    y=max((quantile(log(A)-log(B),pnorm(1))-x)^2,(-quantile(log(A)-log(B),pnorm(-1))+x)^2)
    if (is.infinite(y)) {
        x=mean(log(A+1)-log(B+1))
        y=var(log(A+1)-log(B+1))
    }
    opt.fun=function(v) (digamma(v[1])-digamma(v[2])-x)^2+(trigamma(v[1])+trigamma(v[2])-y)^2
    optim(c(1,1),opt.fun)$par
}

#' Subtract the median of the given vector (for normalizing log2 fold changes).
#'
#' @param l Vector of effect sizes
#' @seealso PsiLFC
#' @return A vector of length 2 containing the two parameters
#' @export
#' @importFrom stats median
CenterMedian <- function(l) l-median(l)

#' Computes the optimal effect size estimate.
#'
#' @param A Vector of counts from condition A
#' @param B Vector of counts from condition B
#' @param prior Vector of length 2 of the prior parameters
#' @param normalizeFun Function to normalize the obtained effect sizes
#'
#' @return A vector of length 2 containing the two parameters
#' @export
PsiLFC=function(A,B, prior=EmpiricalBayesPrior(A,B), normalizeFun=CenterMedian) {
    lfc=(digamma(A+prior[1])-digamma(B+prior[2]))/log(2)
    CenterMedian(lfc)
}

