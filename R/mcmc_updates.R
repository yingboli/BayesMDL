## current = list(eta, inference, change.eta)
## inference = list(mdl, phi, sigmasq, mu, s)

####################################################################
##########  Update one random component of eta at a time  ##########
####################################################################
#' One MCMC iteration: propose to flip one random non-outlier component of eta
#'
#' @inheritParams fit_eta
#' @param current A list object representing a changepoint model. It contains
#'   the following components: \code{eta}, \code{xi}, \code{inference} (output of
#'   the \code{\link{fit_eta}} function), \code{change_eta}, and \code{change_xi}.
#' @param max_changes The maximum number of changepoints, or \code{NULL} if not
#'   specified.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @importFrom stats runif
#' @keywords internal
#'

eta_MH_flip = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                       b_xi, scale_trend_design, weights, max_changes){

  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);

  ## Location of eta that can be flipped
  candidate_eta = 1:n;
  ## An outlier cannot be a changepoint
  if(sum(xi) > 0){ 
    candidate_eta = candidate_eta[-which(xi[candidate_eta] == 1)];
  }
  ## The first max(p, 1) locations cannot be changepoints
  candidate_eta = candidate_eta[candidate_eta > max(p, 1)];
  
  ## The total number of changepoints m cannot exceed max_changes
  if(is.null(max_changes)){
    max_changes = ceiling((length(x) - p - ncol(A) - 1) / 2) - 1;
  }
  if(sum(eta) == max_changes){
    candidate_eta = candidate_eta[eta[candidate_eta] == 1];
  }
  
  ## Randomly select the component to be flipped
  j = sample(candidate_eta, 1);
  eta2 = eta;
  eta2[j] = 1 - eta[j];

  inference2 = fit_eta(x, A, eta2, current$xi, p, fit, penalty, nu, kappa, a, 
                       b_eta, b_xi, scale_trend_design, weights);

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
##########  Update one random component of eta: 0 -> 1    ##########
####################################################################
#' One MCMC iteration: propose to flip one random non-outlier component of eta
#' from 0 to 1
#'
#' @inheritParams fit_eta
#' @param current A list object representing a changepoint model. It contains
#'   the following components: \code{eta}, \code{xi}, \code{inference} (output of
#'   the \code{\link{fit_eta}} function), \code{change_eta}, and \code{change_xi}.
#' @param max_changes The maximum number of changepoints, or \code{NULL} if not
#'   specified.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @importFrom stats runif
#' @keywords internal
#'

eta_MH_birth = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                        b_xi, scale_trend_design, weights, max_changes){
  
  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);
  
  ## The total number of changepoints m cannot exceed max_changes
  if(is.null(max_changes)){
    max_changes = ceiling((length(x) - p - ncol(A) - 1) / 2) - 1;
  }
  if(sum(eta) < max_changes){
    ## Location of eta that can be flipped: among non-changepoints
    candidate_eta = which(eta == 0);
    ## An outlier cannot be a changepoint
    if(sum(xi) > 0){ 
      candidate_eta = candidate_eta[-which(xi[candidate_eta] == 1)];
    }
    ## The first max(p, 1) locations cannot be changepoints
    candidate_eta = candidate_eta[candidate_eta > max(p, 1)];
    
    ## Check candidate_eta
    if(sum(eta[candidate_eta] == 1) > 0){
      stop('Error in eta_MH_birth: wrong candidate_eta.')
    }
  
    ## Randomly select the component to be flipped
    if(length(candidate_eta) == 1){
      j = candidate_eta;
    }
    if(length(candidate_eta) > 1){
      j = sample(candidate_eta, 1);
    }  
    eta2 = eta;
    eta2[j] = 1 - eta[j];
  
    inference2 = fit_eta(x, A, eta2, current$xi, p, fit, penalty, nu, kappa, a, 
                         b_eta, b_xi, scale_trend_design, weights);
  
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
##########  Update one random component of eta: 1 -> 0    ##########
####################################################################
#' One MCMC iteration: propose to flip one random component of eta from 1 to 0
#'
#' @inheritParams fit_eta
#' @param current A list object representing a changepoint model. It contains
#'   the following components: \code{eta}, \code{xi}, \code{inference} (output of
#'   the \code{\link{fit_eta}} function), \code{change_eta}, and \code{change_xi}.
#' @param max_changes The maximum number of changepoints, or \code{NULL} if not
#'   specified.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @importFrom stats runif
#' @keywords internal
#'

eta_MH_death = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                        b_xi, scale_trend_design, weights){
  
  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);
  
  if(sum(eta) > 0){
    ## Location of eta that can be flipped: among non-changepoints
    candidate_eta = which(eta == 1);

    ## Check candidate_eta
    if(sum(eta[candidate_eta] == 0) > 0){
      stop('Error in eta_MH_death: wrong candidate_eta.')
    }
    ## Randomly select the component to be flipped
    if(length(candidate_eta) == 1){
      j = candidate_eta;
    }
    if(length(candidate_eta) > 1){
      j = sample(candidate_eta, 1);
    }
    eta2 = eta;
    eta2[j] = 1 - eta[j];
    
    inference2 = fit_eta(x, A, eta2, current$xi, p, fit, penalty, nu, kappa, a, 
                         b_eta, b_xi, scale_trend_design, weights);
    
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
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'

eta_MH_swap = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                       b_xi, scale_trend_design, weights){

  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  ## Randomly select the component i and j to be swapped
  ## Note that if p == 0, the first time cannot be a changepoint!
  n = length(x);
  if(sum(eta) > 0 && sum(eta) < n - max(p, 1)){
    ## Index i is a current changepoint
    if(sum(eta) == 1)
      i = which(eta == 1);
    if(sum(eta) > 1)
      i = sample(which(eta == 1), 1);
    
    ## Index j is a current non-changepoint and non-outlier
    candidate_eta = which(eta == 0);
    ## An outlier cannot be a changepoint
    if(sum(xi) > 0){ 
      candidate_eta = candidate_eta[-which(xi[candidate_eta] == 1)];
    }
    ## The first max(p, 1) locations cannot be changepoints
    candidate_eta = candidate_eta[candidate_eta > max(p, 1)];
    ## Randomly select the component j 
    if(length(candidate_eta) == 1){
      j = candidate_eta;
    }
    if(length(candidate_eta) > 1){
      j = sample(candidate_eta, 1);
    }  

    ## randomly select the component to be flipped
    eta2 = eta;
    eta2[i] = 1 - eta[i];
    eta2[j] = 1 - eta[j];

    inference2 = fit_eta(x, A, eta2, current$xi, p, fit, penalty, nu, kappa, a,
                         b_eta, b_xi, scale_trend_design, weights);

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
##########   Update one random component of xi at a time  ##########
####################################################################
#' One MCMC iteration: propose to flip one random non-changepoint component of xi
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @param max_outliers The maximum number of outliers, or \code{NULL} if not
#'   specified.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'

xi_MH_flip = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, b_xi,
                      scale_trend_design, weights, max_outliers){

  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);
  candidate_xi = 1:n;
  
  ## If a time is a changepoint, then it cannot be an outlier
  if(sum(eta) > 0){
    candidate_xi = candidate_xi[eta[candidate_xi] == 0];
  }
  ## The total number of outliers l cannot exceed max_outliers
  if(is.null(max_outliers)){
    max_outliers = floor(n / 6);
  }
  if(sum(xi) == max_outliers){
    candidate_xi = candidate_xi[xi[candidate_xi] == 1];
  }
  
  ## Randomly select the component to be flipped
  j = sample(candidate_xi, 1);
  xi2 = xi;
  xi2[j] = 1 - xi[j];

  inference2 = fit_eta(x, A, current$eta, xi2, p, fit, penalty, nu, kappa, a, 
                       b_eta, b_xi, scale_trend_design, weights);

  loga = current$inference$bmdl - inference2$bmdl;
  logu = log(runif(1));
  if(loga > logu){
    current$xi = xi2;
    current$inference = inference2;
    current$change_xi = TRUE;
  }

  return(current);

}

####################################################################
##########    Update one random component of xi: 0 -> 1   ##########
####################################################################
#' One MCMC iteration: propose to flip one random non-changepoint component of 
#' xi from 0 to 1
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @param max_outliers The maximum number of outliers, or \code{NULL} if not
#'   specified.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'

xi_MH_birth = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                       b_xi, scale_trend_design, weights, max_outliers){
  
  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);
  candidate_xi = 1:n;
  
  ## The total number of outliers l cannot exceed max_outliers
  if(is.null(max_outliers)){
    max_outliers = floor(n / 6);
  }
  if(sum(xi) < max_outliers){
    ## Location of xi that can be flipped: among non-outliers
    candidate_xi = which(xi == 0);
    ## An changepoint cannot be an outlier
    if(sum(eta) > 0){ 
      candidate_xi = candidate_xi[-which(eta[candidate_xi] == 1)];
    }
    ## Check candidate_xi
    if(sum(xi[candidate_xi] == 1) > 0){
      stop('Error in xi_MH_birth: wrong candidate_xi.')
    }
    
    ## Randomly select the component to be flipped
    if(length(candidate_xi) == 1){
      j = candidate_xi;
    }
    if(length(candidate_xi) > 1){
      j = sample(candidate_xi, 1);
    }  
    xi2 = xi;
    xi2[j] = 1 - xi[j];
    
    inference2 = fit_eta(x, A, current$eta, xi2, p, fit, penalty, nu, kappa, a, 
                         b_eta, b_xi, scale_trend_design, weights);
    
    loga = current$inference$bmdl - inference2$bmdl;
    logu = log(runif(1));
    if(loga > logu){
      current$xi = xi2;
      current$inference = inference2;
      current$change_xi = TRUE;
    }
  }  
  return(current);
  
}

####################################################################
##########    Update one random component of xi: 1 -> 0   ##########
####################################################################
#' One MCMC iteration: propose to flip one random non-changepoint component of 
#' xi from 1 to 0
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @param max_outliers The maximum number of outliers, or \code{NULL} if not
#'   specified.
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'

xi_MH_death = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                       b_xi, scale_trend_design, weights){
  
  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);
  candidate_xi = 1:n;
  
  if(sum(xi) > 0){
    ## Location of xi that can be flipped: among outliers
    candidate_xi = which(xi == 1);
    ## Check candidate_xi
    if(sum(xi[candidate_xi] == 0) > 0){
      stop('Error in xi_MH_death: wrong candidate_xi.')
    }
    
    ## Randomly select the component to be flipped
    if(length(candidate_xi) == 1){
      j = candidate_xi;
    }
    if(length(candidate_xi) > 1){
      j = sample(candidate_xi, 1);
    }  
    xi2 = xi;
    xi2[j] = 1 - xi[j];
    
    inference2 = fit_eta(x, A, current$eta, xi2, p, fit, penalty, nu, kappa, a, 
                         b_eta, b_xi, scale_trend_design, weights);
    
    loga = current$inference$bmdl - inference2$bmdl;
    logu = log(runif(1));
    if(loga > logu){
      current$xi = xi2;
      current$inference = inference2;
      current$change_xi = TRUE;
    }
  }  
  return(current);
  
}

####################################################################
##########        Flip one changepoint to an outlier      ##########
####################################################################
#' One MCMC iteration: propose to flip one random non-changepoint component of xi
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @inheritParams xi_MH_flip
#' @return A list object representing the (maybe) updated changepoint model.
#'   \item{eta}{The (maybe) updated changepoint model, in the format of
#'     a vector of 0/1 indicators.}
#'   \item{xi}{The outliers, in the format of a vector of 0/1 indicators.}
#'   \item{inference}{Output of the \code{\link{fit_eta}} function on the output
#'     model \code{eta}, \code{xi}.}
#'   \item{change_eta}{Logical, if this \code{eta} is new, i.e.,
#'     the proposed model is accepted.}
#'   \item{change_xi}{Logical, if this \code{xi} is new, i.e.,
#'     the proposed model is accepted.}
#'
#' @export
#' @keywords internal
#'
eta_to_xi_MH_flip = function(x, A, current, p, fit, penalty, nu, kappa, a, b_eta, 
                             b_xi, scale_trend_design, weights, max_outliers){
  
  eta = current$eta;
  xi = current$xi;
  current$change_eta = FALSE;
  current$change_xi = FALSE;
  
  n = length(x);
  
  if(is.null(max_outliers)){
    max_outliers = floor(n / 6);
  }
  
  if(sum(eta) >= 1 && sum(xi) <= max_outliers - 1){
    ## Randomly select the changepoint to be flipped to an outlier
    if(sum(eta) == 1)
      j = which(eta == 1);
    if(sum(eta) > 1){
      candidate_eta = which(eta == 1);
      next_eta = (candidate_eta + 1) %in% candidate_eta;
      prob_eta = rep(1, length(candidate_eta)) + 2 * next_eta;
      prob_eta = prob_eta / sum(prob_eta);
      j = sample(candidate_eta, 1, prob = prob_eta);
    }  
    
    eta2 = eta;
    xi2 = xi;
    eta2[j] = 0;
    xi2[j] = 1;
    
    inference2 = fit_eta(x, A, eta2, xi2, p, fit, penalty, nu, kappa, a,
                         b_eta, b_xi, scale_trend_design, weights);
    
    loga = current$inference$bmdl - inference2$bmdl;
    logu = log(runif(1));
    if(loga > logu){
      current$eta = eta2;
      current$xi = xi2;
      current$inference = inference2;
      current$change_eta = TRUE;
      current$change_xi = TRUE;
    }
  }

  return(current);
  
}