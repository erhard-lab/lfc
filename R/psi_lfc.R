
#' Computes the prior parameters (i.e. pseudocounts incremented by 1) for
#' the log2 fold change distribution
#'
#' @param A Vector of counts from condition A
#' @param B Vector of counts from condition B
#' @seealso PsiLFC
#' @return A vector of length 2 containing the two parameters
#' @export
#' @importFrom stats median pnorm optim var quantile
#' @examples
#'   EmpiricalBayesPrior(rnorm(1000,200),rnorm(1000,100))
EmpiricalBayesPrior=function(A,B) {

    u=A>0 | B>0
    A=A[u]
    B=B[u]

    x=median(log(A)-log(B))
    y=max((quantile(log(A)-log(B),pnorm(1))-x)^2,(-quantile(log(A)-log(B),pnorm(-1))+x)^2)
    if (is.infinite(x) || is.infinite(y)) {
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
#' @examples
#'   CenterMedian(rnorm(1000,200))
CenterMedian <- function(l) l-median(l,na.rm=TRUE)

#' Computes the optimal effect size estimate and credible intervals if needed.
#'
#' @title Psi LFC effect size estimator
#' @param A Vector of counts from condition A
#' @param B Vector of counts from condition B
#' @param prior Vector of length 2 of the prior parameters
#' @param normalizeFun Function to normalize the obtained effect sizes
#' @param cre Compute credible intervals as well? (can also be a vector of quantiles)
#'
#' @return Either a vector containing the estimates, or a data frame containing
#'     the credible interval as well
#' @export
#' @examples
#'   PsiLFC(rnorm(1000,200),rnorm(1000,100))
PsiLFC=function(A,B, prior=EmpiricalBayesPrior(A,B), normalizeFun=CenterMedian,cre=FALSE) {
    lfc<-(digamma(A+prior[1])-digamma(B+prior[2]))/log(2)
    r<-normalizeFun(lfc)
    if (all(cre==TRUE)) cre = c(0.05,0.95)
    if (!missing(cre) & any(cre!=FALSE)){
        s<-sapply(cre,function(c) qlfc(c,A+prior[1],B+prior[2]))
        r<-cbind(data.frame(PsiLFC=lfc),s)
        names(r)[-1]<-paste("Credible",cre)
    }
    r
}


#' Computes the optimal effect size estimate and credible intervals if needed
#' for a Bioconductor SummarizedExperiment object
#'
#' @title Psi LFC effect size estimator
#' @param se SummarizedExperiment object
#' @param contrast Vector of length 3 (<name>,<A>,<B>)
#' @param cre Compute credible intervals as well? (can also be a vector of quantiles)
#'
#' @return Either a vector containing the estimates, or a data frame containing
#'     the credible interval as well
#' @export
#' @examples
#'    data(airway, package="airway")
#'    head(PsiLFC.se(airway,contrast=c("dex","untrt","trt")))
PsiLFC.se=function(se,contrast,cre=FALSE) {
    if (!is.character(contrast) | length(contrast)!=3 | !(contrast[1] %in% names(colData(se))) | contrast[2]==contrast[3] ) {
        stop("'contrast' vector should be a character vector of length (colData,A,B)")
    }

    f=do.call("colData",list(se))[,contrast[1]]
    if (!all(contrast[2:3] %in% levels(f))){
        stop("'contrast' vector contains unknown levels!")
    }

    summ<-function(cont) {
        r=do.call("assay",list(se))[,f==cont]
        if (!is.vector(r)) r=apply(r,1,sum)
        r
    }
    A=summ(contrast[2])
    B=summ(contrast[3])
    PsiLFC(A,B,cre=cre)
}


#' Drop-in replacement for DESeq2's results function for simple settings involving
#' a single variable. Appends the PsiLFC estimate.
#'
#' @title Psi LFC effect size estimator for DESeq2
#' @param object the DESeq2DataSet object
#' @param contrast Vector of length 3, specifying the variable and the two levels
#'     to compute effect sizes for (<name>,<A>,<B>)
#' @param cre Compute credible intervals as well? (can also be a vector of quantiles)
#' @param ... Handed over to DESeq2's results function
#'
#' @return Either a vector containing the estimates, or a data frame containing
#'     the credible interval as well
#' @export
results <- function(object, contrast, cre=FALSE, ...) {
    r=PsiLFC.se(object,contrast,cre)
    if ("DESeq2" %in% loadedNamespaces()) {
        if (is.vector(r)) r = data.frame(PsiLFC=r)
        r=cbind(DESeq2::results(object,contrast,...),r)
    }
    r
}


