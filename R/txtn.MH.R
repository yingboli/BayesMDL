###########################################################
########       Tmax and Tmin joint detection         ######
###########################################################
###### Inputs:
## X: a N by 2 matrix, each column is re-processed to zero mean.
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


############################################################
#### current: list(eta, count.eta, inference, change.eta) ####

#' Metropolis-Hastings Algorithms for Stoachstic Model Search, in Tmax / Tmin
#' Joint Detection
#'
#' Two proposals: the function \code{txtn.MH.flip} propose to flip a random
#' dimension, and the function \code{txtn.MH.swap} propose to swap between a
#' category 4 (no changes in either series) and a catogory 1-3 (change in at
#' least one series).
#'
#' @inheritParams txtn.eta
#' @inheritParams txtn
#' @param current A list of current changepoint model. Specifically,
#'   \code{list(eta, count.eta, inference, change.eta)}, where \code{inference}
#'   is the output of the function \code{\link{txtn.eta}}.
#' @return
#' A list of current changepoint model. Specifically,
#'   \code{list(eta, count.eta, inference, change.eta)}, where \code{inference}
#'   is the output of the function \code{\link{txtn.eta}}.
#' @export
#' @name txtn.MH

####################################################################
##########  Update one random component of eta at a time  ##########
####################################################################
#' @rdname txtn.MH
#' @export
#'
txtn.MH.flip = function(X, month, current, meta, p, nu, alpha1, alpha2, period, q){

  N = dim(X)[1]; ## length of time series

  eta = current$eta;
  count.eta = current$count.eta;

  ## propose a new eta
  eta2 = eta;
  j = sample((p + 1):N, 1);
  ## New proposal for the j-th time
  indtmp = eta[j, 3];
  eta2[j, 3] = sample((1:4)[-indtmp], 1, prob = q[indtmp, -indtmp]);
  eta2[j, 1:2] = c((4 - eta2[j, 3]) %/% 2, (4 - eta2[j, 3]) %% 2);

  ## change in count.eta
  count.eta2 = count.eta;
  count.eta2[eta[j, 3]] = count.eta2[eta[j, 3]] - 1;
  count.eta2[eta2[j, 3]] = count.eta2[eta2[j, 3]] + 1;

  inference2 = txtn.eta(X, month, eta2, meta, p, nu, alpha1, alpha2, period);
  loga = current$inference$mdl - inference2$mdl + log(q[eta2[j, 3], eta[j, 3]]) - log(q[eta[j, 3], eta2[j, 3]]);
  logu = log(runif(1));
  if(loga > logu){
  	current$eta = eta2;
  	current$count.eta = count.eta2;
  	current$change.eta = TRUE;
  	current$inference = inference2;
  }
  if(loga <= logu)
    current$change.eta = FALSE;

  return(current);

}

####################################################################
##########       Simple random swapping update of eta     ##########
####################################################################
#' @rdname txtn.MH
#' @export
#'
txtn.MH.swap = function(X, month, current, meta, p, nu, alpha1, alpha2, period){

  N = dim(X)[1]; ## length of time series

  eta = current$eta;
  count.eta = current$count.eta;

  ## propose a new eta
  allj = which(eta[, 3] == 4)[-(1:p)];
  alli = which(eta[, 3] != 4);

  if(length(alli) == 0 || length(allj) == 0){
    current$change.eta = FALSE;
  	return( current );
  }
  if(length(alli) == 1)
    i = alli;
  if(length(alli) > 1)
    i = sample(alli, 1);

  if(length(allj) == 1)
    j = allj;
  if(length(allj) > 1)
    j = sample(allj, 1);

  eta2 = eta;
  eta2[i, ] = eta[j, ];
  eta2[j, ] = eta[i, ];

  ## change in count.eta
  count.eta2 = current$count.eta;

  inference2 = txtn.eta(X, month, eta2, meta, p, nu, alpha1, alpha2, period);
  loga = current$inference$mdl - inference2$mdl;
  logu = log(runif(1));
  if(loga > logu){
  	current$eta = eta2;
  	current$count.eta = count.eta2;
  	current$change.eta = TRUE;
  	current$inference = inference2;
  }
  if(loga <= logu)
    current$change.eta = FALSE;

  return(current);
}




