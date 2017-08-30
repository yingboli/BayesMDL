## current = list(eta, inference, change.eta)
## inference = list(mdl, phi, sigmasq, mu, s)

####################################################################
##########  Update one random component of eta at a time  ##########
####################################################################
#' One MCMC iteration: propose to flip one random component of eta
#'
#' @inheritParams fit_eta
#' @param current A list object representing a changepoint model. It contains
#'   the following components: \code{eta}, \code{xi}, \code{inference} (output of
#'   the \code{\link{fit_eta}} function), \code{change_eta}, and \code{change_xi}.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @importFrom stats runif
#' @keywords internal
#'

eta_MH_flip = function(x, A, current, p, fit, penalty, nu, kappa, a, b,
                        scale_trend_design, weights){

  eta = current$eta;
  current$change_eta = FALSE;

  ## Randomly select the component to be flipped
  ## Note that if p == 0, the first time cannot be a changepoint!
  n = length(x);
  j = sample((max(p, 1) + 1):n, 1);
  eta2 = eta;
  eta2[j] = 1 - eta[j];

  inference2 = fit_eta(x, A, eta2, current$xi, p, fit, penalty, nu, kappa, a, b,
                       scale_trend_design, weights);

  loga = current$inference$bmdl - inference2$bmdl;
  logu = log(runif(1));
  if(loga > logu){
    current$eta = eta2;
    current$inference = inference2;
    current$change_eta = TRUE;
  }

  return(current);

}

####################################################################
##########       Simple random swapping update of eta     ##########
####################################################################
#' One MCMC iteration: propose to swap a changepoint with a non-changepoint
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'

eta_MH_swap = function(x, A, current, p, fit, penalty, nu, kappa, a, b,
                        scale_trend_design, weights){

  eta = current$eta;
  current$change_eta = FALSE;

  ## Randomly select the component i and j to be swapped
  ## Note that if p == 0, the first time cannot be a changepoint!
  n = length(x);
  if(sum(eta) > 0 && sum(eta) < n - max(p, 1)){
    if(sum(eta) == 1)
      i = which(eta == 1);
    if(sum(eta) > 1)
      i = sample(which(eta == 1), 1);

    ## For p > 0
    if(sum(eta) == n - p - 1 && p > 0)
      j = which(eta == 0)[p + 1];
    if(sum(eta) < n - p - 1 && p > 0)
      j = sample(which(eta == 0)[-(1:p)], 1);
    ## For p == 0
    if(sum(eta) == n - p - 2 && p == 0)
      j = which(eta == 0)[2];
    if(sum(eta) < n - p - 2 && p == 0)
      j = sample(which(eta == 0)[-1], 1);


    ## randomly select the component to be flipped
    eta2 = eta;
    eta2[i] = 1 - eta[i];
    eta2[j] = 1 - eta[j];

    inference2 = fit_eta(x, A, eta2, current$xi, p, fit, penalty, nu, kappa, a,
                         b, scale_trend_design, weights);

    loga = current$inference$bmdl - inference2$bmdl;
    logu = log(runif(1));
    if(loga > logu){
      current$eta = eta2;
      current$inference = inference2;
      current$change_eta = TRUE;
    }
  }

  return(current);

}

####################################################################
##########  Update one random component of xi at a time  ##########
####################################################################
#' One MCMC iteration: propose to flip one random component of xi
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'

xi_MH_flip = function(x, A, current, p, fit, penalty, nu, kappa, a, b,
                      scale_trend_design, weights){

  xi = current$xi;
  current$change_xi = FALSE;

  ## Randomly select the component to be flipped
  ## Note that if p == 0, the first time can be an outlier!
  n = length(x);
  j = sample((p + 1):n, 1);
  xi2 = xi;
  xi2[j] = 1 - xi[j];

  inference2 = fit_eta(x, A, current$eta, xi2, p, fit, penalty, nu, kappa, a, b,
                       scale_trend_design, weights);

  loga = current$inference$bmdl - inference2$bmdl;
  logu = log(runif(1));
  if(loga > logu){
    current$xi = xi2;
    current$inference = inference2;
    current$change_xi = TRUE;
  }

  return(current);

}
