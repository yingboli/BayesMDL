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

##########################################################################
######  Tmax/Tmin: Bayesian MDL for multiple changepoint detection #######
##########################################################################

#' BMDL for Changepoint Detection, in Tmax / Tmin Joint Detection
#'
#' For biviariate time series. It accomodates metadata (if any), and can handel
#' monthly data and annual data.
#'
#' @param X A \code{N} by 2 matrix, observed time series data. Columns are Tmax
#'   and Tmin, respectively.
#' @param month A vector of length \code{N}, month indicators. Takes value in
#'   {1, 2, ..., 12} for monthly data, or {1} for annual data.
#' @param meta A 0-1 indicator vector of length \code{N}, metadata indicator.
#'   The first \code{p} elemenet are always 0.
#' @param iter Total number of MCMC iterations.
#' @param thin Thinning; save one MCMC iteration for every \code{thin} number of
#'   iterations.
#' @param type Whether data are observed 'annual' or 'monthly'. Default is
#'   'monthly'.
#' @param p The order of the VAR process.
#' @param nu Prior variance scale of \code{mu}; only used if
#'   \code{fit == 'marlik'}.
#' @param alpha1 The parameter in the Multinomial-Dirichlet prior of \code{eta},
#'   for un-documented times (i.e., times not in metadata). Default value is
#'   \code{c(3/7, 2/7, 2/7, 239)} for monthly data, and
#'   \code{c(3/7, 2/7, 2/7, 19)} for annual data.
#' @param alpha2 The parameter in the Multinomial-Dirichlet prior of \code{eta},
#'   for documented times (i.e., times in metadata). Default value is
#'   \code{c(3/7, 2/7, 2/7, 47)} for monthly data, and \code{c(3/7, 2/7, 2/7, 3)}
#'   for annual data.
#' @param q A 4 by 4 transition matrix of the proposal distribution, only used
#'   in the flip proposal.
#' @param start.eta A 0-1 indicator vector of length \code{2*N}, initial
#'   changepoint model for MCMC. The entries correspond to the first \code{p}
#'   times are always 0. Default value is \code{NULL}.
#' @param track.time Logical; whether to report processing time (in seconds).
#' @param show.summary Integer; the number of top models to print. Default value
#'   is 10.
#' @param show.month Logical; if show month in the summary. Default is
#'   \code{FALSE}.
#' @param start.year Integer; the year time 1 is in. Default value is 1.
#'
#' @return
#' \item{mcmc}{A (\code{iter / thin} + 1) by (2 * \code{N} + 4 * \code{p} + 5)
#'   matrix; each row is a saved MCMC output from an iteration, consisting
#'   \code{eta}, bmdl score, \code{Phi}, and \code{Sigma}.}
#' \item{map200}{A 200 by (2 *\code{N} + 4 * \code{p} + 5) matrix; the best 200
#'   changepoint models visited. Each row is a model, consisting \code{eta},
#'   bmdl score, \code{Phi}, and \code{Sigma}.}
#' \item{input.parameters}{A list of all the following input parameters:
#'   \code{X}, \code{month}, \code{meta}, \code{iter}, \code{thin}, \code{type},
#'   \code{p}, \code{nu}, \code{alpha1}, \code{alpha2}, \code{q}, \code{period},
#'   and \code{start.year}.}
#'
#' @export
#' @examples
#' data(tuscaloosa);
#' X = as.matrix(tuscaloosa[, 3:4]);
#'
#' ## For illustration purpose, here iter is small.
#' ## To get meaningfull inference, use a larger value, e.g., iter = 5e4.
#' results = txtn(X, month = tuscaloosa[, 2], meta = tuscaloosa[, 7], iter = 1e1);


txtn = function(X, month = NULL, meta = NULL, iter = 1e4, thin = max(1, iter / 1e3),
				type = 'monthly', p = 3, nu = 5,
				alpha1 = c(3/7, 2/7, 2/7, 19) * (type == 'annual')
				+ c(3/7, 2/7, 2/7, 239) * (type == 'monthly'),
				alpha2 = c(3/7, 2/7, 2/7, 3) * (type == 'annual')
				+ c(3/7, 2/7, 2/7, 47) * (type == 'monthly'),
				q = matrix(c(0, 1/2, 1/2, 1/2,
				1/4, 0, 0, 1/4, 1/4, 0, 0, 1/4, 1/2, 1/2, 1/2, 0), ncol = 4, nrow = 4),
				start.eta = NULL, track.time = TRUE, show.summary = 10,
				show.month = FALSE, start.year = 1){

  t.start = proc.time();

  change.rate = 0;

  N = dim(X)[1]; ## length of time series

  ## pre-process meta
  ## for now, non-dirac meta can only be used for annual data!!!
  if(length(meta) < N){
  	meta = loc2dirac(meta, N);
  }
  if(length(meta) > N)
    stop('Length of metadata cannot exceed length of time series.')

  ## output matrices
  map200 = matrix(NA, ncol = 2 * (N + 2 * p + 2) + 1, nrow = 200);
  mcmc = matrix(NA, ncol = 2 * (N + 2 * p + 2) + 1, nrow = round(iter / thin) + 1);
  colnames(map200) = colnames(mcmc) = c(paste(c(rep('TX', N), rep('TN', N)), rep(1:N, 2), sep = ''), 'mdl', paste('Phi', 1:(4 * p), sep = ''), paste('Sigma', 1:4, sep = ''));

  ## initial values of eta
  if(length(start.eta) == 0){
    eta = matrix(0, ncol = 3, nrow = N);
    eta[meta == 0, 3] = sample(1:4, sum(meta == 0), prob = alpha1 / sum(alpha1), replace = TRUE);
    eta[meta == 1, 3] = sample(1:4, sum(meta == 1), prob = alpha2 / sum(alpha2), replace = TRUE);
    eta[1:p, 3] = 4
    for(t in (p + 1):N){
      if(eta[t, 3] == 2)
        eta[t, 1] = 1;
      if(eta[t, 3] == 3)
        eta[t, 2] = 1;
      if(eta[t, 3] == 1)
        eta[t, 1:2] = 1;
    }
  }
  if(length(start.eta) > 0){
  	eta = matrix(start.eta, ncol = 2);
    eta = cbind(eta, rep(NA, N));
  	for(t in 1:N){
      if(eta[t, 1] == 1 && eta[t, 2] == 1)
        eta[t, 3] = 1;
      if(eta[t, 1] == 1 && eta[t, 2] == 0)
        eta[t, 3] = 2;
      if(eta[t, 1] == 0 && eta[t, 2] == 1)
        eta[t, 3] = 3;
      if(eta[t, 1] == 0 && eta[t, 2] == 0)
        eta[t, 3] = 4;
  	}
  }
  count.eta = rep(0, 4);
  tmp = table(eta[-(1:p), 3]);
  count.eta[as.numeric(names(tmp))] = tmp;

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

  inference = txtn.eta(X, month, eta, meta, p, nu, alpha1, alpha2, period);
  current = list(eta = eta, count.eta = count.eta, inference = inference, change.eta = TRUE);
  mcmc[1, ] = map200[1, ] = c(current$eta[ , 1:2], unlist(current$inference)[1:(2 * (2 * p + 2) + 1)]);

  #### MCMC ####
  for(it in 1:iter){
  	#it = 1;
  	## update eta
  	propose = sample(c('flip', 'swap'), 1, prob = c(0.8, 0.2), replace = TRUE);
  	if(propose == 'flip')
  	  current = txtn.MH.flip(X, month, current, meta, p, nu, alpha1, alpha2, period, q);
  	if(propose == 'swap')
  	  current = txtn.MH.swap(X, month, current, meta, p, nu, alpha1, alpha2, period);

  	if(current$change.eta == TRUE){
      change.rate = change.rate + 1 / iter;

      ## update map200
      tmpk = which(current$inference$mdl > map200[, 2 * N + 1]);
      if(length(tmpk) == 0)
        k = 0;
      if(length(tmpk) > 0)
  	    k = max(tmpk);
  	  if(k < 200 && is.na(map200[k + 1, 2 * N + 1]))
  	    map200[(k + 1), ] = c(current$eta[ , 1:2], unlist(current$inference)[1:(2 * (2 * p + 2) + 1)]);
  	  if(k < 200 && sum(map200[k + 1, 1:(2 * N)] != c(current$eta[, 1:2])) > 0){
  	    if(k == 199)
  	      map200[200, ] = c(current$eta[ , 1:2], unlist(current$inference)[1:(2 * (2 * p + 2) + 1)]);
  	    if(k < 199){
  	      map200[(k + 2):200, ] = map200[(k + 1):199, ];
  	      map200[(k + 1), ] = c(current$eta[ , 1:2], unlist(current$inference)[1:(2 * (2 * p + 2) + 1)]);
  	    }
  	  }
    }

  	## save MCMC outputs
  	if(it %% thin == 0){
  	  mcmc[it / thin + 1, ] = c(current$eta[ , 1:2], unlist(current$inference)[1:(2 * (2 * p + 2) + 1)]);

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
  	    cat('Tmax:', which(map200[r, 1:N] == 1) + start.year - 1, 'Tmin:', which(map200[r, N + 1:N] == 1) + start.year - 1, 'mdl =', map200[r, 2 * N + 1], '\n', sep = ' ');
  	}
    if(type == 'monthly' && show.month == TRUE){
      year = (start.year + c(rep(0, 13 - month[1]), c(matrix(rep(1 : ceiling((N - (13 - month[1])) / 12), 12), nrow = 12, byrow = TRUE))))[1 : N];
      yearmonth = paste(year, month.abb[month], sep = '')
      for(r in 1:show.summary)
  	    cat('Tmax:', yearmonth[which(map200[r, 1:N] == 1)], 'Tmin:', yearmonth[which(map200[r, N + 1:N] == 1)], 'mdl =', map200[r, 2 * N + 1], '\n', sep = ' ');
  	}
  	if(type == 'monthly' && show.month == FALSE){
      for(r in 1:show.summary)
  	    cat('Tmax:', which(map200[r, 1:N] == 1), 'Tmin:', which(map200[r, N + 1:N] == 1), 'mdl =', map200[r, 2 * N + 1], '\n', sep = ' ');
  	}

  }

  input.parameters = list(X = X, month = month, meta = meta, iter = iter, thin = thin, type = type, p = p, nu = nu, alpha1 = alpha1, alpha2 = alpha2, q = q, period = period, start.year = start.year);

  ## track time
  if(track.time == TRUE){
    t.finish = proc.time();
    cat('Time used (in second): \n')
    print(t.finish - t.start);
    cat('\n');
  }

  return( list(mcmc = mcmc, map200 = map200, input.parameters = input.parameters) );

}


