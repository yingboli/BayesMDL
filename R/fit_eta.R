###########################################################
#########   Fit an Individual Changepoint Model    ########
###########################################################
#' Fit an individual changepoint model \code{eta}
#'
#' Compute the objective function of the model, i.e., log of posterior proablity
#' for BMDL, or log MDL, and estimate model parameters.
#'
#' @param x The time series data, a numeric vector of length \code{n}.
#' @param A The design matrix for the nuisance coefficients in the linear model.
#'   It is usually the matrix of seasonal indicators, or the design matrix
#'   for harmonic regression with a column of all 1 for intercept.
#' @param eta A multiple changepoint configuration, i.e., model. It is a 0/1
#'   indicator vector of length \code{n}, with the first \code{p} elemenets
#'   always being 0.
#' @param xi Outliers indicator, a 0/1 vector of length \code{n}.
#' @param p The order of the AR process.
#' @param fit For likelihood calculation, \code{'marlik'} for marginal
#'   likelihood, or \code{'lik'} for likelihood. Note that the \code{'lik'}
#'   option already includes the two-part MDL of \code{mu}.
#' @param penalty For penalty function calculation, \code{'bmdl'} for
#'   Beta-Binomial prior, or \code{'mdl'} for MDL.
#' @param nu Prior variance scale of \code{mu}; only used if
#'   \code{fit == 'marlik'}.
#' @param kappa Prior variance scale of outliers.
#' @param a,b The first and second parameters in the Beta-Binomial prior; only
#'   used if \code{penalty == 'bmdl'}.
#' @param scale_trend_design The factor multiplied to the design matrix of trend.
#'   Default is 1/50.
#' @param weights A numeric vector of observation weights, defined the same as
#'   the \code{weights} argument in the function \code{lm}.
#'
#' @return
#' \item{bmdl}{BMDL or MDL.}
#' \item{s}{Estimates of nuisance coefficients in the linear model. It usually
#'   contains the seasonal means or coefficients of harmonic regression
#'   (with intercept).}
#' \item{mu}{Estimates of regime-wise coefficients in the linear model. It
#'   usually contains the regime means and regime trends.}
#' \item{sigmasq}{Estimate of \eqn{\sigma^2} in the AR(p) process.}
#' \item{phi}{Yule-Walker estimates of the AR(p) coefficients.}
#'
#' @export
#' @keywords internal


fit_eta = function(x, A, eta, xi, p = 2, fit = 'marlik', penalty = 'bmdl',
                   nu = 5, kappa = 3, a = 1, b = 1, scale_trend_design = 0.05,
                   weights = NULL){

  m = sum(eta); k = sum(xi); n = length(x);

  ## Design matrix for regime means and regime trend
  tmp = D_eta(x, eta, scale_trend_design);
  D = tmp$D;
  seg_lengths = tmp$seg_lengths;

  ## Weights: incoporate outliers
  if(is.null(weights))
    weights = rep(1, n);
  if(k > 0){
    weights[xi == 1] = weights[xi == 1] / kappa^2;
  }

  ###### Compute the Yule-Walker estimate of phi ######
  ## Linear model residuals
  y = lm(x ~ cbind(A, D) - 1, weights = weights)$res;

  ## If AR, not independent, compute phi
  phi = NULL;
  if(p > 0){
    gamma = rep(NA, p + 1); ## gamma_hat of lag h: h = 0, 1, ..., p
    for(h in 0:p){
      gamma[h + 1] = sum(c(y, rep(0, h)) * c(rep(0, h), y)) / n;
    }
    phi = solve(matrix(gamma[abs(outer(1:p, 1:p, '-')) + 1], nrow = p), gamma[-1]);
  }


  ###### Compute Ahat, Xhat, and Dhat ######
  Ahat = as.matrix(A[(p + 1):n, ]);
  xhat = x[(p + 1):n];
  Dhat = as.matrix(D[(p + 1):n, ]);

  if(p > 0){
    for(h in 1:p){
      Ahat = as.matrix(Ahat - phi[h] * A[(p + 1 - h):(n - h), ]);
      xhat = xhat - phi[h] * x[(p + 1 - h):(n - h)];
      Dhat = as.matrix(Dhat - phi[h] * D[(p + 1 - h):(n - h), ]);
    }
  }

  ## Left multipler each of Ahat, xhat, Dhat by diag(sqrt(weights))
  wsqrt = sqrt(weights)[(p + 1):n];
  Ahat = sweep(Ahat, 1, wsqrt, FUN = '*');
  xhat = xhat * wsqrt;
  Dhat = sweep(Dhat, 1, wsqrt, FUN = '*');

  #### Compute EB estimate of sigmasq and conditional posterior mean of mu
  tmp = svd(Dhat);
  UtX = t(tmp$u) %*% xhat;
  UtA = t(tmp$u) %*% Ahat;
  if(fit == 'marlik')
    d = tmp$d^2 / (tmp$d^2 + 1 / nu);
  if(fit == 'lik')
    d = rep(1, ncol(D));

  if(ncol(D) > 1){
    AtBA = t(Ahat) %*% Ahat - t(UtA) %*% diag(d) %*% UtA;
  }
  if(ncol(D) == 1){
    AtBA = t(Ahat) %*% Ahat - t(UtA) %*% (d * UtA);
  }
  AtBX = t(Ahat) %*% xhat - t(UtA) %*% (d * UtX);
  XtBX = sum(xhat^2) - sum(UtX^2 * d);

  s = solve(AtBA, AtBX);

  if(ncol(D) > 1){
    mu = tmp$v %*% diag(tmp$d / (tmp$d^2 + 1 / nu)) %*% (UtX - UtA %*% s);
  }
  if(ncol(D) == 1){
    mu = tmp$v * (tmp$d / (tmp$d^2 + 1 / nu)) * (UtX - UtA %*% s);
  }

  sigmasq = c((XtBX - t(AtBX) %*% s) / (n - p));

  ## Conditional posterior covariance of mu: sigmasq * inv(Dhat' Dhat + I / nu)
  cov.mu = solve(t(Dhat) %*% Dhat + diag(ncol(D)) / nu) * sigmasq;

  if(fit == 'marlik')
    mdl.data = (n - p) / 2 * log(sigmasq) + ncol(D) / 2 * log(nu) +
               sum(log(tmp$d^2 + 1 / nu)) / 2;
  if(fit == 'lik')
    mdl.data = (n - p) / 2 * log(sigmasq) + sum(log(seg_lengths[-1])) / 2 +
               sum(log(seg_lengths)) / 2;

  ###### Compute bmdl (i.e., log_plik) #####
  if(penalty == 'bmdl')
    bmdl = mdl.data - lgamma(a + m) - lgamma(b + n - p - m) - lgamma(a + k) -
           lgamma(b + n - p - k);

  if(penalty == 'mdl')
    bmdl = mdl.data + log(m + 1) + (m + 1) * log(n - p) + log(k + 1) +
           (k + 1) * log(n - p);

  ## return list "inference"
  return( list(bmdl = c(bmdl), phi = phi, sigmasq = c(sigmasq),
               mu = c(mu), s = c(s), cov.mu = cov.mu) );

}
