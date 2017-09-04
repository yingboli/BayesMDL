###########################################################
#########             Harmonic Regression          ########
###########################################################
#' Harmonic regression
#'
#' @param x Input time series of length \code{n}.
#' @param k The highest order of the factor inside sin/cos.
#' @param time_ind A index vector of times, or \code{NULL}, in which case
#'   \code{1:n} will be used.
#' @param intercept Logical, whether to include an intercept.
#' @param return_design Logical, whether to return the design matrix.
#'
#' @return A linear regession \code{lm} object.
#'
#' @export
#' @importFrom stats lm
#' @keywords internal
#'
harmonic_lm = function(x, k = 3, time_ind = NULL, intercept = TRUE,
                       return_design = TRUE){

  n = length(x);

  if(is.null(time_ind)){
    time_ind = 1:n;
  }

  ## Design matrix
  if(k > 0){
    design_sin = outer(time_ind, 1:k, function(t, k){sin(2 * pi / 12 * t * k)});
    design_cos = outer(time_ind, 1:k, function(t, k){cos(2 * pi / 12 * t * k)});
    design = data.frame(design_sin, design_cos);
    names(design) = c(paste0('sin', 1:k), paste0('cos', 1:k));

    ## Add intercept
    if(intercept){
      design$intercept = 1;
    }
    ## Also include x
    design$x = x;
  }

  if(k == 0){
    design = data.frame(intercept = rep(1, n), x = x);
  }
  ## Linear regression
  return(lm(x ~ . -1, data = design, x = return_design));
}

###########################################################
######     Create eta from non-zero-one sequence    #######
###########################################################
#' Convert a Location Vector to a 0-1 Indicator Vector
#'
#' Convert a location vector to a 0-1 indicator vector, such that the
#'   corresponding entries in these locations are 1.
#'
#' @param loc An integer vector, length does not exceed \code{N}.
#' @param N The length of the resulting indicator vector.
#' @return A 0-1 indicator vector of length \code{N}.
#' @export
#' @keywords internal

loc2dirac = function(loc, N){
  dirac = rep(0, N);
  dirac[loc] = 1;
  return(dirac);
}

###########################################################
######   Create a design matrix D for a model eta   #######
###########################################################
#' Create a design matrix D for a model eta
#'
#' @inheritParams fit_eta
#' @return
#' \item{D}{A \code{n} by \code{(2*m+1)} matrix. The first \code{m}
#'   columns are for regime means; they are 0/1 indicators for regimes from 2 to
#'   \code{m}. The next \code{(m+1)} columns are for trends.}
#' \item{seg_lengths}{A vector of length \code{m}, regime lengths, including
#'   the first regime. }
#'
#' @export
#' @keywords internal
#'

D_eta = function(x, eta, scale_trend_design){

  m = sum(eta); n = length(x);

  ## Design matrix for regime means: n by m + 1 (eventually: n by m)
  ## Later, the column for the 1st regime will be deleted in D_mu
  D_mu = matrix(as.numeric(matrix(rep(cumsum(eta) + 1, m + 1), ncol = m + 1) ==
                             matrix(rep(1:(m + 1), n), ncol = m + 1, byrow = TRUE)),
                ncol = m + 1);

  seg_lengths = apply(D_mu, 2, sum);

  ## Design matrix for regime trends: n by (m + 1)
  D_alpha = matrix(rep(1:n), nrow = n, ncol = m + 1);
  D_alpha = (D_alpha - matrix(c(0, which(eta == 1) - 1) + (1 + seg_lengths) / 2,
                              nrow = n, ncol = m + 1, byrow = TRUE)) * D_mu;
  D_alpha = D_alpha * scale_trend_design;

  ## Combine the two design matrices of regime parameters
  ## Note that even for the null model, the D matrix still contains 1 column of trend
  D = cbind(D_mu[, -1], D_alpha);

  return(list(D = D, seg_lengths = seg_lengths));
}

###########################################################
#' Location of the last changepoint
#'
#' @inheritParams fit_eta
#' @return An iterger, the location of the last changepoint, or \code{NA} if
#'   there are no changepoints.
#' @export
#' @keywords internal
#'
last_cp = function(eta){
  if(sum(eta) == 0){
    return(NA);
  } else {
    return(max(which(eta == 1)));
  }
}
