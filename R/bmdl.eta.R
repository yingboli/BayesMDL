########################################################################
######     Calculating marginal likelihood for a given eta       #######
########################################################################
###### AR(p); monthly data; no trend; can handle metadata

###### Inputs:
## X: a vector of length N.
##    For monthly data, no need to pre-processed to zero mean.
## month: a vector of length N. Take value in {1, 2, ..., 12}.
##		For annual data, always equal 1.
## eta: changepoint configuration.
##      Length is N, first p elemenet are always 0.
## meta: a vector of 0-1 indicators.
##      Length is N, first p elemenet are always 0.
## p: AR(p)
## fit: 'marlik' marginal likelihood, or 'lik' likelihood.
##      Note, the 'lik' option already included the two-part mdl of mu.
## penalty: 'bmdl' Beta-Binomial prior, 'uniform', 'mdl', 'BIC'
## nu: default value 5
## a: default value 1
## b1: for non-metadata, dafault value 239 (for monthly data)
## b2: for metadata, dafault value 47 (for monthly data)
## period = 12 (for monthly data), or 1 (for annual data)

###### Outputs:
## mdl: Bayesian mdl score for model eta (i.e., lpost)
## mu: conditional posterior mean of mu (Note: the length of mu is m)
## s: EB estimators to seasonalities (Note: the length of s is 12)
## sigmasq: EB estimator
## phi: Yule-Walker estimator
########################################################################

#' Calculate BMDL for a Given Changepoint Model
#'
#' For mean shifts detection in a univariate time series. These functions can
#' handle metadata.
#'
#' \itemize{
#' \item To compute BMDL, use \code{fit = 'marlik'} and \code{penalty = 'bmdl'}.
#' \item To compute automatic MDL, use \code{fit = 'lik'} and
#'   \code{penalty = 'mdl'}.
#' \item To compute BIC, use \code{penalty = 'BIC'} (the \code{fit} argument
#'   does not matter in this case).
#' \item Note, the marginal likelihood \code{fit = 'marlik'} can also be
#'   combined with \code{penalty = 'uniform'}, which assigns uniform prior on
#'   the model space.
#' }
#'
#' @inheritParams bmdl.mean.shift
#' @param eta A 0-1 indicator vector of length \code{N}, changepoint model. The
#'   first \code{p} elemenet are always 0.
#' @param period Period of the seasonal cycle, 12 for monthly data, or 1 for
#'   annual data.
#' @return
#' \item{mdl}{The Bayesian MDL score for the changepoint model \code{eta}.}
#' \item{phi}{A vector of length \code{p}, the Yule-Walker estimator of the
#'   autocorrelation parameter \code{phi}.}
#' \item{sigmasq}{The EB estimator of the white noise variance \code{sigmasq}.}
#' \item{mu}{A vector of length \code{sum(eta)}; the conditional posterior mean
#'   of regime-wise means \code{mu}.}
#' \item{s}{A vector of length \code{period}, the EB estimator of the monthly
#'   means \code{s}. }
#' @export

bmdl.eta = function(X, month, eta, meta, p, fit, penalty, nu, a, b1, b2, period){

  m = sum(eta); N = length(X);

  ## design matrix for seasonality
  A = matrix(as.numeric(matrix(month, ncol = period, nrow = N) == matrix(1:period, ncol = period, nrow = N, byrow = TRUE)), ncol = period);
  ## design matrix for regime means
  D = matrix(as.numeric(matrix(rep(cumsum(eta) + 1, m + 1), ncol = m + 1) == matrix(rep(1:(m + 1), N), ncol = m + 1, byrow = TRUE)), ncol = m + 1);
  D = as.matrix(D[, -1]);
  seg.lengths = apply(D, 2, sum);

  #### compute the Yule-Walker estimate of phi ####
  ## linear model residuals
  Y = lm(X ~ cbind(A, D) - 1)$res;

  ## if AR, not independent, compute phi
  phi = NULL;
  if(p > 0){
    gamma = rep(NA, p + 1); ## gamma_hat of lag h: h = 0, 1, ..., p
    for(h in 0:p){
  	  gamma[h + 1] = sum(c(Y, rep(0, h)) * c(rep(0, h), Y)) / N;
    }
    phi = solve(matrix(gamma[abs(outer(1:p, 1:p, '-')) + 1], nrow = p), gamma[-1]);
  }


  #### compute Ahat, Xhat, and Dhat
  Ahat = as.matrix(A[(p + 1):N, ]);
  Xhat = X[(p + 1):N];
  if(m > 0)
    Dhat = as.matrix(D[(p + 1):N, ]);
  if(p > 0){
  	for(h in 1:p){
  	  Ahat = as.matrix(Ahat - phi[h] * A[(p + 1 - h):(N - h), ]);
  	  Xhat = Xhat - phi[h] * X[(p + 1 - h):(N - h)];
  	  if(m > 0)
  	    Dhat = as.matrix(Dhat - phi[h] * D[(p + 1 - h):(N - h), ]);
  	}
  }

  ## if not the null model
  if(m > 0){

    #### compute EB estimate of sigmasq and conditional posterior mean of mu
    tmp = svd(Dhat);
    UtX = t(tmp$u) %*% Xhat;
    UtA = t(tmp$u) %*% Ahat;
    if(fit == 'marlik')
      d = tmp$d^2 / (tmp$d^2 + 1 / nu);
    if(fit == 'lik')
      d = rep(1, m);

    if(m > 1){
      AtBA = t(Ahat) %*% Ahat - t(UtA) %*% diag(d) %*% UtA;
    }
   if(m == 1){
      AtBA = t(Ahat) %*% Ahat - t(UtA) %*% (d * UtA);
    }
    AtBX = t(Ahat) %*% Xhat - t(UtA) %*% (d * UtX);
    XtBX = sum(Xhat^2) - sum(UtX^2 * d);

    s = solve(AtBA, AtBX);

    if(m > 1){
  	  mu = tmp$v %*% diag(tmp$d / (tmp$d^2 + 1 / nu)) %*% (UtX - UtA %*% s);
    }
    if(m == 1){
  	  mu = tmp$v * (tmp$d / (tmp$d^2 + 1 / nu)) * (UtX - UtA %*% s);
    }

    sigmasq = c((XtBX - t(AtBX) %*% s) / (N - p));

    ## conditional posterior covariance of mu: sigmasq * inv(Dhat' Dhat + I / nu)
    cov.mu = solve(t(Dhat) %*% Dhat + diag(m) / nu) * sigmasq;

    if(fit == 'marlik')
      mdl.data = (N - p) / 2 * log(sigmasq) + m / 2 * log(nu) + sum(log(tmp$d^2 + 1 / nu)) / 2;
    if(fit == 'lik')
      mdl.data = (N - p) / 2 * log(sigmasq) + sum(log(seg.lengths)) / 2;
  }

  ## if the null model
  if(m == 0){
    s = solve(t(Ahat) %*% Ahat, t(Ahat) %*% Xhat);
    mu = NULL; cov.mu = NULL;
    sigmasq = c((t(Xhat) %*% Xhat - t(Xhat) %*% Ahat %*% s) / (N - p));

    mdl.data = (N - p) / 2 * log(sigmasq);
  }

  #### compute bmdl
  if(penalty == 'bmdl'){
    m2 = sum(eta == 1 & meta == 1); N2 = sum(meta);
    m1 = m - m2; N1 = N - p - N2;
    mdl = mdl.data - lgamma(a + m1) - lgamma(b1 + N1 - m1) - lgamma(a + m2) - lgamma(b2 + N2 - m2);
  }

  if(penalty == 'mdl')
  	mdl = mdl.data + log(m + 1) + (m + 1) * log(N - p);

  if(penalty == 'uniform')
  	mdl = mdl.data;

  if(penalty == 'BIC')
    mdl = (N - p) / 2 * log(sigmasq) + m / 2 * log(N - p);

  ## return list "inference"
  return( list(mdl = c(mdl), phi = phi, sigmasq = c(sigmasq), mu = c(mu), s = c(s), cov.mu = cov.mu) );
}

#######################
  # start.time = proc.time();
  # for(i in 1:1e5){
    # #AtBA = AtA - t(UtA) %*% sweep(UtA, 1, d, FUN = "*");
    # AtBA = AtA - t(UtA) %*% diag(d) %*% UtA
  # }
  # proc.time() - start.time;
