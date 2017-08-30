###########################################################
########       Tmax and Tmin joint detection         ######
###########################################################
###### Inputs:
## X: a N by 2 matrix.
## month: a vector of length N. Take value in {1, 2, ..., 12}.
##		For annual data, always equal 1.
## eta: a N by 3 matrix, changepoint configuration.
##      First two columns are 0-1 indicator of changepoints for tmax and tmin, respectively,
##      and the third column is the category 1-4 each time is in.
## meta: a vector of 0-1 indicators.
##      Length is N, first p elemenet are always 0.
## p: VAR(p)
## nu: default value 5
## alpha1, alpha2: default values are c(3/7, 2/7, 2/7, 239) and c(3/7, 2/7, 2/7, 47)
## period = 12 (for monthly data), or 1 (for annual data)

####################################################################
##########     Compute the bmdl for a given model eta     ##########
####################################################################
#' Calculate BMDL for a Given Changepoint Model, in Tmax / Tmin Joint Detection
#'
#' For mean shifts detection in a bivariate time series. These functions can
#' handle metadata.
#'
#' @inheritParams txtn
#' @param eta A \code{N} by 3 matrix, changepoint model. The first two columns
#'   are 0-1 indicators of changepoints for tmax and tmin, respectively, and the
#'   third column is the category 1-4 each time is in. The first \code{p} times
#'   are not allowed to be changepoints.
#' @param period Period of the seasonal cycle, 12 for monthly data and 1 for
#'   annual data.
#' @return
#' \item{mdl}{The Bayesian MDL score for the changepoint model \code{eta}.}
#' \item{Phi}{A 2\code{p} by 2 matrix, the Yule-Walker estimator of the
#'   autocorrelation parameters \code{Phi_1}, ..., \code{Phi_p}, in a block
#'   matrix representation.}
#' \item{Sigma}{A 2 by 2 matrix, the Yule-Walker estimator of the white noise
#'   covariance \code{Sigma}.}
#' \item{mu}{A vector of length \code{sum(eta[ , 1]) + sum(eta[ , 2])}; the
#'   conditional posterior mean of regime-wise means \code{(mu_1, mu_2)}.}
#' \item{s}{A vector of length \code{period}, the EB estimator of the monthly
#' means \eqn{s}. }
#' @export
#' @importFrom MASS ginv


txtn.eta = function(X, month, eta, meta, p, nu, alpha1, alpha2, period){

  N = dim(X)[1]; ## length of time series

  m1 = sum(eta[, 1]);
  m2 = sum(eta[, 2]);

  ## design matrix for seasonality
  A = matrix(as.numeric(matrix(month, ncol = period, nrow = N) == matrix(1:period, ncol = period, nrow = N, byrow = TRUE)), ncol = period);

  ## design matrix for regime means
  D1 = matrix(as.numeric(matrix(rep(cumsum(eta[, 1]) + 1, m1 + 1), ncol = m1 + 1) == matrix(rep(1:(m1 + 1), N), ncol = m1 + 1, byrow = TRUE)), ncol = m1 + 1);
  D1 = as.matrix(D1[, -1]);

  D2 = matrix(as.numeric(matrix(rep(cumsum(eta[, 2]) + 1, m2 + 1), ncol = m2 + 1) == matrix(rep(1:(m2 + 1), N), ncol = m2 + 1, byrow = TRUE)), ncol = m2 + 1);
  D2 = as.matrix(D2[, -1]);

  ## SUR residuals
  H = cbind(A, D1, A, D2);
  invPsi = solve((t(X) %*% X - (t(X) %*% H) %*% ginv(t(H) %*% H) %*% (t(H) %*% X)) / N);

  G1 = cbind(A, D1); G2 = cbind(A, D2);
  invGtG = solve(rbind(cbind(invPsi[1, 1] * t(G1) %*% G1, invPsi[1, 2] * t(G1) %*% G2), cbind(invPsi[2, 1] * t(G2) %*% G1, invPsi[2, 2] * t(G2) %*% G2)));
  GtX = c(t(G1) %*% (invPsi[1, 1] * X[, 1] + invPsi[1, 2] * X[, 2]), t(G2) %*% (invPsi[2, 1] * X[, 1] + invPsi[2, 2] * X[, 2]));
  Y = matrix(c(X) - rbind(cbind(G1, 0 * G2), cbind(0 * G1, G2)) %*% invGtG %*% GtX, ncol = 2);

  #### compute the Yule-Walker estimate of Phi and Sigma ####
  if(p > 0){
  	Gamma = matrix(NA, nrow = 2, ncol = (2 * p + 1) * 2);
  	Gamma[ , p * 2 + 1:2] = t(Y) %*% Y / N;
  	for(h in 1:p){
  	  tmp = t(Y[1:(N - h), ]) %*% Y[(h + 1):N, ] / N;
  	  Gamma[ , (p + h) * 2 + 1:2] = tmp;
  	  Gamma[ , (p - h) * 2 + 1:2] = t(tmp);
  	}
  	Gamma.sym = NULL;
  	for(h in 1:p){
  	  Gamma.sym = rbind(Gamma.sym, Gamma[ , ((p - h) * 2 + 3):((2 * p - h + 1) * 2)]);
  	}
    tmp = matrix(NA, ncol = 2, nrow = p);
    for(h in 1:p){
      tmp[h, ] = c(((p - h) * 2 + 3), ((2 * p - h + 1) * 2));
    }

    Phi = Gamma[ , (p * 2 + 3):((2 * p + 1) * 2)] %*% solve(Gamma.sym);
    Sigma = Gamma[ , p * 2 + 1:2] - Phi %*% t(Gamma[ , (p * 2 + 3):((2 * p + 1) * 2)]);
  }

  if(p == 0){
  	Sigma = t(Y) %*% Y / N;
  	Phi = NULL;
  }

  invSigma = solve(Sigma);

  ## hatDt = Dt - Phi_1 %*% D_(t-1) - ... - Phi_p %*% D_(t-p)
  ## hatAt = At - Phi_1 %&% A_(t-1) - ... - Phi_p %*% A_(t-p)
  ## hatXt = At - Phi_1 %&% X_(t-1) - ... - Phi_p %*% X_(t-p)
  Dv = cbind(D1, 0 * D2, 0 * D1, D2) ## vectorize the two rows of time t into one row
  Av = cbind(A, 0 * A, 0 * A, A) ## vectorize the two rows of time t into one row

  hatDv = Dv[(p + 1):N, ];
  hatAv = Av[(p + 1):N, ];
  hatX = X[(p + 1):N, ];

  if(p > 0){
    for(h in 1:p){
      ## for a given Phi, compute Phi_h %*% dt
      Phi.dt = function(dt){
  	    return(c(t(Phi[ , (h - 1) * 2 + 1:2] %*% matrix(dt, nrow = 2, byrow = TRUE))));
      }
      hatAv = hatAv - t(apply(Av[(p + 1 - h): (N - h), ], 1, Phi.dt)); ## nrow = N - p
      hatX = hatX - t(apply(X[(p + 1 - h): (N - h), ], 1, Phi.dt)); ## length = N - p
      if(m1 + m2 > 0)
        hatDv = hatDv - t(apply(Dv[(p + 1 - h): (N - h), ], 1, Phi.dt)); ## nrow = N - p
    }
  }

  ## for a given invSigma, compute sum(t(hat..) %*% inv(Sigma) %*% hat..)
  tAinvSA = t(matrix(hatAv[1, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% matrix(hatAv[1, ], nrow = 2, byrow = TRUE);
  tXinvSX = t(hatX[1, ]) %*% invSigma %*% hatX[1, ];
  tAinvSX = t(matrix(hatAv[1, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% hatX[1, ];
  if(m1 + m2 > 0){
    tDinvSD = t(matrix(hatDv[1, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% matrix(hatDv[1, ], nrow = 2, byrow = TRUE);
    tDinvSA = t(matrix(hatDv[1, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% matrix(hatAv[1, ], nrow = 2, byrow = TRUE);
    tDinvSX = t(matrix(hatDv[1, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% hatX[1, ];
  }

  for(i in 2:(N - p)){
    tAinvSA = tAinvSA + t(matrix(hatAv[i, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% matrix(hatAv[i, ], nrow = 2, byrow = TRUE);
    tXinvSX = tXinvSX + t(hatX[i, ]) %*% invSigma %*% hatX[i, ];
    tAinvSX = tAinvSX + t(matrix(hatAv[i, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% hatX[i, ];
    if(m1 + m2 > 0){
      tDinvSD = tDinvSD + t(matrix(hatDv[i, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% matrix(hatDv[i, ], nrow = 2, byrow = TRUE);
      tDinvSA = tDinvSA + t(matrix(hatDv[i, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% matrix(hatAv[i, ], nrow = 2, byrow = TRUE);
      tDinvSX = tDinvSX + t(matrix(hatDv[i, ], nrow = 2, byrow = TRUE)) %*% invSigma %*% hatX[i, ];
    }
  }

  if(m1 + m2 > 0){
  	if(m1 + m2 == 1)
      invF = solve(tDinvSD + 1 / nu * c(rep(1 / Sigma[1, 1], m1), rep(1 / Sigma[2, 2], m2)));
  	if(m1 + m2 >= 2)
      invF = solve(tDinvSD + 1 / nu * diag(c(rep(1 / Sigma[1, 1], m1), rep(1 / Sigma[2, 2], m2))));
  	s = solve(tAinvSA - t(tDinvSA) %*% invF %*% tDinvSA) %*% (tAinvSX - t(tDinvSA) %*% invF %*% tDinvSX);
  	mu = invF %*% (tDinvSX - tDinvSA %*% s);
  	mdl.data = 0.5 * ( (N - p) * log(det(Sigma)) + m1 * log(nu * Sigma[1, 1]) + m2 * log(nu * Sigma[2, 2]) - log(det(invF)) +  (tXinvSX - t(tDinvSX) %*% invF %*% tDinvSX) - t(tAinvSX - t(tDinvSA) %*% invF %*% tDinvSX) %*% s);
  }
  if(m1 + m2 == 0){
  	s = solve(tAinvSA) %*% tAinvSX;
  	mu = NULL; invF = NULL;
  	mdl.data = 0.5 * ( (N - p) * log(det(Sigma)) + tXinvSX - t(tAinvSX) %*% s );
  }

  ## compute bmdl
  N2 = sum(meta); N1 = N - p - N2;
  ml1 = ml2 = rep(0, 4);
  tmp1 = table(eta[meta == 0, 3]);
  ml1[as.numeric(names(tmp1))] = tmp1;
  ml1[4] = ml1[4] - p;
  tmp2 = table(eta[meta == 1, 3]);
  ml2[as.numeric(names(tmp2))] = tmp2;

  mdl = mdl.data - sum(lgamma(alpha1 + ml1)) - sum(lgamma(alpha2 + ml2));

  ## return list "inference"
  return(list(mdl = c(mdl), Phi = Phi, Sigma = Sigma, mu = c(mu), s = matrix(s, ncol = 2), cov.mu = invF));

}

############################################################
#### current: list(eta, count.eta, inference, change.eta) ####

# start.time = proc.time();
# for(i in 1:1e2){

  # hatDv = Dv[(p + 1):N, ];
  # hatAv = Av[(p + 1):N, ];
  # hatX = X[(p + 1):N, ];

  # for(i in 1:(N - p)){
  	# for(h in 1:p){
  	  # hatDv[i, ] = hatDv[i, ] - c(t(Phi[ , (h - 1) * 2 + 1:2] %*% matrix(Dv[i + p - h, ], nrow = 2, byrow = TRUE)));
  	  # hatAv[i, ] = hatAv[i, ] - c(t(Phi[ , (h - 1) * 2 + 1:2] %*% matrix(Av[i + p - h, ], nrow = 2, byrow = TRUE)));
  	  # hatX[i, ] = hatX[i, ] - c(t(Phi[ , (h - 1) * 2 + 1:2] %*% X[i + p - h, ]));
  	# }
  # }

# }
# proc.time() - start.time;


