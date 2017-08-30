###################################################################
######    Bayesian MDL for multiple changepoint detection   #######
###################################################################

######## Inputs #######
## X: a vector of length N. NOT pre-processed to zero mean.
## month: a vector of length N. Take value in {1, 2, ..., 12}.
## meta: a vector of 0-1 indicators.
##      Length is N, first p elemenets are always 0.
## nu: default value 5
## a: default value 1
## b1: dafault value 19 (for annual data), 239 (for monthly data)
## b2: dafault value 3 (for annual data), 47 (for monthly data)
## type: 'annual' or 'monthly'

######## Outputs ########
## map200: 200 * (N + 3) matrix
##         200 top models with highest posterior probability that MCMC visited
##         along with mdl, phi, sigmasq
## Eta: (iter / thin + 1) * (N + 3) matrix
##         MCMC samples of change point configuration (model) eta

###################################################################
#' BMDL for Changepoint Detection
#'
#' For univiariate time series. It accomodates metadata (if any), and can handel
#' monthly data and annual data.
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
#' Please see the references for more details.
#'
#' @inheritParams bmdl
#' @param X A vector of length \code{N}, observed time series data.
#' @param month A vector of length \code{N}, month indicators. Takes value in
#'   {1, 2, ..., 12} for monthly data, or {1} for annual data.
#' @param meta A 0-1 indicator vector of length \eqn{N}, metadata indicator.
#'   The first \code{p} elemenet are always 0.
#' @param type Whether data are observed 'annual' or 'monthly'. Default is
#'   'monthly'.
#' @param b1 The \eqn{b} parameter in the Beta-Binomial prior of \eqn{\eta}, for
#'   un-documented times (i.e., times not in metadata). Default value is 239 for
#'   monthly data, and 19 for annual data.
#' @param b2 The \eqn{b} parameter in the Beta-Binomial prior of \eqn{\eta}, for
#'   documented times (i.e., times in metadata). Default value is 47 for monthly
#'   data, and 3 for annual data.
#' @param start.eta A 0-1 indicator vector of length \code{N}, initial
#'   changepoint model for MCMC. The first \code{p} elemenet are always 0.
#'   Default value is \code{NULL}.
#' @param track.time Logical; whether to report processing time (in seconds).
#' @param show.summary Integer; the number of top models to print. Default value
#'   is 10.
#' @param show.month Logical; if show month in the summary. Default is
#'   \code{FALSE}.
#' @param start.year Integer; the year time 1 is in. Default value is 1.
#'
#' @return
#' \item{mcmc}{A (\code{iter / thin} + 1) by (\code{N} + \code{p} + 2) matrix;
#'   each row is a saved MCMC output from an iteration, consisting \code{eta},
#'   bmdl score, \code{phi}, and \code{sigmasq}.}
#' \item{map200}{A 200 by (\code{N} + \code{p} + 2) matrix; the best 200
#'   changepoint models visited. Each row is a model, consisting \code{eta},
#'   bmdl score, \code{phi}, and \code{sigmasq}.}
#' \item{input.parameters}{A list of all the following input parameters:
#'   \code{X}, \code{month}, \code{meta}, \code{iter}, \code{thin}, \code{type},
#'   \code{p}, \code{nu}, \code{a}, \code{b1}, \code{b2}, \code{period}, and
#'   \code{start.year}.}
#' @export
#' @references Yingbo Li, Robert Lund, and Anuradha Hewaarachchi,
#' "Multiple Changepoint Detection with Partial Information on Changepoint Times".
#' Working paper.
#'
#' @examples
#' data(tuscaloosa);
#' X = tuscaloosa[, 3];
#'
#' ## For illustration purpose, here iter is small.
#' ## To get meaningfull inference, use a larger value, e.g., iter = 1e5.
#' results = bmdl.mean.shift(X, month = tuscaloosa[, 2], meta = tuscaloosa[, 7],
#'                           iter = 1e2);
#'

bmdl.mean.shift = function(X, month = NULL, meta = NULL, iter = 1e4,
                           thin = max(1, iter / 1e3), type = 'monthly', p = 3,
                           fit = 'marlik', penalty = 'bmdl', nu = 5, a = 1,
                           b1 = 19 * (type == 'annual') + 239 * (type == 'monthly'),
                           b2 = (b1 - 4) / 5, start.eta = NULL, track.time = TRUE,
                           show.summary = 10, show.month = FALSE, start.year = 1){

  t.start = proc.time();

  change.rate = 0;

  N = length(X);

  ## pre-process meta
  ## for now, non-dirac meta can only be used for annual data!!!
  if(length(meta) < N){
  	meta = loc2dirac(meta, N);
  }
  if(length(meta) > N)
    stop('Length of metadata cannot exceed length of time series.')

  ## output matrices
  map200 = matrix(NA, ncol = N + p + 2, nrow = 200);
  mcmc = matrix(NA, ncol = N + p + 2, nrow = round(iter / thin) + 1);
  colnames(map200) = colnames(mcmc) = c(paste('eta', 1:N, sep = ''), 'mdl', paste('phi', 1:p, sep = ''), 'sigmasq');

  ## initial values of eta
  if(length(start.eta) == 0)
    eta = rbinom(N, 1, sort(runif(2, 0, 0.05))[meta + 1]);	## vector of change point, start value
  if(length(start.eta) == N)
    eta = start.eta;
  if(length(start.eta) < N && length(start.eta) > 0)
    eta = loc2dirac(start.eta, N);
  if(length(start.eta) > N)
    stop('Error in start.eta: length of changepoint configuration cannot exceed length of time series.')
  eta[1:p] = 0;

  ## initial values
  if(type == 'annual'){
  	month = rep(1, N);
  	period = 1;
  }
  if(type == 'monthly'){
  	if(length(month) == 0)
  	  month = rep(1:12, ceiling(N / 12))[1:N];
  	period = 12;
  }
  inference = bmdl.eta(X, month, eta, meta, p, fit, penalty, nu, a, b1, b2, period);
  current = list(eta = eta, inference = inference, change.eta = TRUE);
  mcmc[1, ] = map200[1, ] = unlist(current)[1 : (N + p + 2)];

  ## Start MCMC
  for(it in 1:iter){

  	## Use MH.flip with probablity 80%, use MH.swap with probability 20%
    action = sample(c('flip', 'swap'), 1, prob = c(0.8, 0.2), replace = TRUE);
    ## Metropolis-Hastings update gamma
    if(action == 'flip')
  	  current = bmdl.MH.flip(X, month, current, meta, p, fit, penalty, nu, a, b1, b2, period);
    if(action == 'swap')
  	  current = bmdl.MH.swap(X, month, current, meta, p, fit, penalty, nu, a, b1, b2, period);

    ## try to put the new eta to map200 if good
    if(current$change.eta == TRUE){
      change.rate = change.rate + 1 / iter;

      tmpk = which(c(current$inference$mdl) > map200[, N + 1]);
      if(length(tmpk) == 0)
        k = 0;
      if(length(tmpk) > 0)
  	    k = max(tmpk);
  	  if(k < 200 && is.na(map200[k + 1, N + 1]))
  	    map200[(k + 1), ] = unlist(current)[1 : (N + p + 2)];
  	  if(k < 200 && sum(map200[k + 1, 1:N] != current$eta) > 0){
  	    if(k == 199)
  	      map200[200, ] = unlist(current)[1 : (N + p + 2)];
  	    if(k < 199){
  	      map200[(k + 2):200, ] = map200[(k + 1):199, ];
  	      map200[(k + 1), ] = unlist(current)[1 : (N + p + 2)];
  	    }
  	  }
  	}

  	## save Markov chains
    if(it %% thin == 0){
  	  mcmc[it / thin + 1, ] = unlist(current)[1 : (N + p + 2)];

      ## show: x0% completed
      if( (it * 10) %% iter == 0 && it != iter  && track.time == TRUE)
        cat(paste( (it * 100) / iter), '% completed...\n', sep = '');
      if( it == iter  && track.time == TRUE){
        cat(paste( (it * 100) / iter), '% completed.\n', sep = '');
        cat('\n');
      }
    }
  }

  ## show summary: top 10 models
  if(show.summary > 0){

  	cat('accept eta: ', change.rate, '\n');
    cat('Top changepoint configurations: \n');
    if(type == 'annual'){
      for(r in 1:show.summary)
  	    cat(which(map200[r, 1:N] == 1) + start.year - 1, ' mdl =', map200[r, N + 1], '\n', sep = ' ');
  	}
    if(type == 'monthly' && show.month == TRUE){
      year = (start.year + c(rep(0, 13 - month[1]), c(matrix(rep(1 : ceiling((N - (13 - month[1])) / 12), 12), nrow = 12, byrow = TRUE))))[1 : N];
      yearmonth = paste(year, month.abb[month], sep = '')
      for(r in 1:show.summary)
  	    cat(yearmonth[which(map200[r, 1:N] == 1)], ' mdl =', map200[r, N + 1], '\n', sep = ' ');
  	}
    if(type == 'monthly' && show.month == FALSE){
      for(r in 1:show.summary)
  	    cat(which(map200[r, 1:N] == 1), ' mdl =', map200[r, N + 1], '\n', sep = ' ');
  	}

  }

  input.parameters = list(X = X, month = month, meta = meta, iter = iter, thin = thin, type = type, p = p, nu = nu, a = a, b1 = b1, b2 = b2, period = period, start.year = start.year);

  ## track time
  if(track.time == TRUE){
    t.finish = proc.time();
    cat('Time used (in second): \n')
    print(t.finish - t.start);
    cat('\n');
  }

  return( list(mcmc = mcmc, map200 = map200, input.parameters = input.parameters) );
}





