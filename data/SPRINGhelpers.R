# For simple example for SPRING method

# To implement Kendall correlation estimates on huge.
# original huge function can only take covariance or data matrix.
# huge.mb funciton returns beta values (MB coefficient estimates)
# huge function with method="mb" does not return beta values.
# to do network visualization (for edge color, I need beta)

# data = n by p matrix data. usually through pulsar, data will receive subsamples.
# lambda = a vector of lambda values
# type = a type of variables. "trunc" is default.
# sym = "or" is the symmetrizing rule of the output graphs. If sym = "and", the edge between node i and node j is selected ONLY when both node i and node j are selected as neighbors for each other. If sym = "or", the edge is selected when either node i or node j is selected as the neighbor for each other. The default value is "or". (refer to huge manual)
hugeKmb <- function(data, lambda, type = "trunc", sym = "or", verbose = TRUE) {
  S    <- estimateR(data, type = type)$R
  est  <- huge.mb(S, lambda, sym = sym, verbose = verbose)
  est
}


# SPRING function: Semi-Parametric Rank-based approach for INference in Graphical model. SPRING follows the neighborhood selection methodology outlined in "mb" Meinshausen and Buhlmann and 
# qdata = n by p matrix of quantitative microbiome data 
# type = a type of variables. "trunc" is default. It could be either continuous, binary or trunc.
# nlambda = the number of lambda sequence values
# lambdaseq = a sequence of decreasing positive numbers to control the regularization. Users can specify a sequence to override the default sequence.
# seed = the seed for subsampling
# nc = number of cores to use subsampling
# thresh = threshold for StARS selection criterion. 0.1 is recommended.
# subsample.ratio = 0.8 is default. The recommended values are 10*sqrt(n)/n for n > 144 or 0.8 otherwise.
# rep.num = the repetition number of subsampling

SPRING <- function(qdata, type = "trunc", fun = hugeKmb, lambda.min.ratio = 1e-2, nlambda = 50, lambdaseq = NULL, seed = 10010, ncores = 2, thresh = 0.1, subsample.ratio = 0.8, rep.num = 50){
  p <- ncol(qdata)
  Kcor <- estimateR(qdata, type = type)$R
  
  if(is.null(lambdaseq)){
    # generate lambda sequence
    lambda.max <- max(max(Kcor-diag(p)), -min(Kcor-diag(p)))
    lambda.min <- lambda.min.ratio * lambda.max
    lambdaseq <- exp(seq(log(lambda.max), log(lambda.min), length = nlam))
  }
  
  out1.K_count <- pulsar::pulsar(qdata, fun = fun, fargs = list(lambda = lambdaseq), rep.num = rep.num, criterion = 'stars', seed = seed, ncores = nc, thresh = thresh, subsample.ratio = subsample.ratio)
  
  fit1.K_count <- pulsar::refit(out1.K_count)
  
  return(list(Kcor = Kcor, output = out1.K_count, fit = fit1.K_count, lambdaseq = lambdaseq))
}




# modified central log ratio transformation
# dat = raw count data or compositional data (n by p) does not matter.
# base = exp(1) for natural log
# tol = tolerance for checking zeros
# For eps and atleast, users do not have to specify any values. Default should be enough.
# eps = epsilon in eq (2) of the paper "Yoon, Gaynanova, M\"{u}ller (2019), Frontiers in Genetics". positive shifts to all non-zero compositions. Refer to the paper for more details. eps = absolute value of minimum of log ratio counts plus c.
# atleast = 1: constant c which ensures all nonzero values to be strictly positive. default is 1.
mclr <- function(dat, base = exp(1), tol = 1e-16, MARGIN = 1, eps = NULL, atleast = 1){
  dat <- as.matrix(dat)
  nzero <- (dat >= tol)  # index for nonzero part
  LOG <- ifelse(nzero, log(dat, base), 0.0) # take log for only nonzero values. zeros stay as zeros.
  
  # centralize by the log of "geometric mean of only nonzero part" # it should be calculated by each row.
  if (nrow(dat) > 1){
    clrdat <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0.0)
  } else if (nrow(dat) == 1){
    clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
  }
  
  if (is.null(eps)){
    if(atleast < 0){ 
      warning("atleast should be positive. The functions uses default value 1 instead.")
      atleast = 1
    }
    if( min(clrdat) < 0 ){ # to find the smallest negative value and add 1 to shift all data larger than zero.
      positivecst <- abs(min(clrdat)) + atleast # "atleast" has default 1.
    }else{
      positivecst <- 0
    }
    # positive shift
    ADDpos <- ifelse(nzero, clrdat + positivecst, 0.0) ## make all non-zero values strictly positive.
    return(ADDpos)
  } else if(eps == 0){
    ## no shift. clr transform applied to non-zero proportions only. without pseudo count.
    return(clrdat)
  } else if(eps > 0){
    ## use user-defined eps for additional positive shift.
    ADDpos <- ifelse(nzero, clrdat + eps, 0.0)
    return(ADDpos)
  } else {
    stop("check your eps value for additional positive shift. Otherwise, leave it as NULL.")
  }
} 