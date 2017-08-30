
## current = list(eta, inference, change.eta)
## inference = list(mdl, phi, sigmasq, mu, s)

####################################################################
##########  Update one random component of eta at a time  ##########
####################################################################

#' Metropolis-Hastings Algorithms for Stoachstic Model Search
#'
#' Two proposals: the function \code{bmdl.MH.flip} proposes to flip a random
#'   dimension, and the function \code{bmdl.MH.swap} propose to swap between a
#'   randomly chosen zero and a randomly chosen one.
#'
#' @inheritParams bmdl.eta
#' @param current A list of current changepoint model. Specifically,
#'   \code{list(eta, inference, change.eta)}, where \code{inference} is the
#'   output of the function \code{\link{bmdl.eta}}.
#' @return
#' A list of current changepoint model. Specifically,
#'   \code{list(eta, inference, change.eta)}, where \code{inference} is the
#'   output of the function \code{\link{bmdl.eta}}.
#' @name bmdl.MH
NULL
####################################################################


#' @rdname bmdl.MH
#' @export

bmdl.MH.flip = function(X, month, current, meta, p, fit, penalty, nu, a, b1, b2, period){

  eta = current$eta;
  current$change.eta = FALSE;

  ## randomly select the component to be flipped
  N = length(X);
  j = sample((p + 1):N, 1);
  eta2 = eta;
  eta2[j] = 1 - eta[j];

  inference2 = bmdl.eta(X, month, eta2, meta, p, fit, penalty, nu, a, b1, b2, period);

  loga = current$inference$mdl - inference2$mdl;
  logu = log(runif(1));
  if(loga > logu){
  	current$eta = eta2;
  	current$inference = inference2;
  	current$change.eta = TRUE;
  }

  return(current);

}

####################################################################
##########       Simple random swapping update of eta     ##########
####################################################################
#' @rdname bmdl.MH
#' @export

bmdl.MH.swap = function(X, month, current, meta, p, fit, penalty, nu, a, b1, b2, period){

  eta = current$eta;
  current$change.eta = FALSE;

  ## randomly select the component i and j to be swapped
  N = length(X);
  if(sum(eta) > 0 && sum(eta) < N - p){
  	if(sum(eta) == 1)
   	  i = which(eta == 1);
   	if(sum(eta) > 1)
   	  i = sample(which(eta == 1), 1);

   	if(sum(eta) == N - p - 1)
   	  j = which(eta == 0)[p + 1];
   	if(sum(eta) < N - p - 1)
      j = sample(which(eta == 0)[-(1:p)], 1);

    ## randomly select the component to be flipped
    eta2 = eta;
    eta2[i] = 1 - eta[i];
    eta2[j] = 1 - eta[j];

    inference2 = bmdl.eta(X, month, eta2, meta, p, fit, penalty, nu, a, b1, b2, period);

    loga = current$inference$mdl - inference2$mdl;
    logu = log(runif(1));
    if(loga > logu){
  	  current$eta = eta2;
  	  current$inference = inference2;
  	  current$change.eta = TRUE;
    }
  }

  return(current);

}
