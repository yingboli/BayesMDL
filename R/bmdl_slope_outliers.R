###################################################################
######    Bayesian MDL for multiple changepoint detection   #######
###################################################################
#' Bayesian MDL (or MDL) for multiple changepoint detection
#'
#' Changes in both mean and linear trend, which permitting a global seasonal
#' mean and AR(p) errors.
#'
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @inheritParams xi_MH_flip
#' @param dates The dates records are observed, a \code{POSIXlt} vector. It should
#'   have the \code{mon} component for months, and the \code{year} component
#'   for years. See also \code{\link{strptime}}.
#' @param iter Total number of MCMC iterations.
#' @param thin Thinning; save one MCMC iteration for every \code{thin} number of
#'   iterations.
#' @param time_unit Default is \code{'month'} (\code{'week'} may be added in the
#'   future).
#' @param seasonal_means The seasonal means variables in the linear model.
#'   Either \code{'fixed_effects'} such that each season has a coefficient
#'   (with no intercept), or \code{'harmonic'} for harmonic regression (with
#'   an intercept).
#' @param k The highest degree of harmonic regression. It is only used if
#'   the argument \code{seasonal_means == 'harmonic'}.
#' @param start_eta A vector of 0/1 indicators for the initial model, or
#'   \code{NULL} to randomly sample an initial model.
#' @param start_xi A vector of 0/1 indicators for the initial outliers, or
#'   \code{NULL} to randomly sample an initial model.
#' @param track_time Logical, whether to show runtime on screen.
#' @param detect_outliers Logical, whether to detect outliers. 
#' @return
#' \item{best}{The optimal model which minimizes BMDL or MDL; a list representing
#'   the model, which contains components: \code{eta}, \code{inference} (output
#'   of the \code{\link{fit_eta}} function), and \code{change_eta}.}
#' \item{eta_mcmc}{A matrix to save MCMC iterations of \code{eta}. Each row
#'   is about an iteration, in the format of a vector of length \code{n + p + 2},
#'   containing \code{eta}, \code{bmdl}, \code{phi}, \code{sigmasq}.}
#' \item{A}{The design matrix for the nuisance coefficients in the linear model.
#'   It is usually the matrix of seasonal indicators, if the argument
#'   \code{seasonal_means = 'fixed_effects'}, or the design matrix
#'   for harmonic regression with a column of all 1 for intercept, if the
#'   argument \code{seasonal_means = 'harmonic'}).}
#' \item{runtime}{Runtime, in second.}
#' \item{input_parameters}{A list of input parameters.}
#'
#' @export
#' @importFrom stats rbinom
#' @keywords internal
#'

bmdl = function(x, dates, iter = 1e4, thin = max(1, iter / 1e3), weights = NULL,
                p = 2, time_unit = 'month', seasonal_means = 'harmonic', k = 3,
                scale_trend_design = 0.05, fit = 'marlik', penalty = 'bmdl',
                nu = 5, kappa = 3, a = 1, b_eta = length(x), b_xi = length(x),
                max_changes = NULL, max_outliers = NULL, start_eta = NULL, 
                start_xi = NULL, track_time = TRUE, detect_outliers = TRUE){

  ## Start time. To be used to compute runtime.
  t.start = proc.time();

  change_rate = 0;
  n = length(x);

  ###### Design matrix A for seasonal means ######
  ## If each month gets a different mean
  if(seasonal_means == 'fixed_effects'){
    if(time_unit == 'month'){
      seasons = dates$mon + 1;
      period = 12;

      ## The design matrix A contains 12 columns
      A = matrix(as.numeric(matrix(seasons, ncol = period, nrow = n) ==
                            matrix(1:period, ncol = period, nrow = n,
                                   byrow = TRUE)), ncol = period);
    }
  }
  ## If using harmonic regression of degree k
  if(seasonal_means == 'harmonic'){
    if(time_unit == 'month'){
      time_ind = (dates$year - dates$year[1]) * 12 +
                 (dates$mon - dates$mon[1]) + 1;

      ## Fit a harnomic regression
      hlm1 = harmonic_lm(x, k = k, time_ind = time_ind);
      A = hlm1$x;
    }
  }

  ###### Initial preparation for MCMC ######
  ## Output matrices
  eta_mcmc = xi_mcmc = matrix(NA, ncol = n, nrow = round(iter / thin) + 1);
  para_mcmc = matrix(NA, ncol = p + 2, nrow = round(iter / thin) + 1);

  ## Initial values of eta
  if(length(start_eta) == 0)
    eta = rbinom(n, 1, 0.02);
  if(length(start_eta) == n)
    eta = start_eta;
  if(length(start_eta) < n && length(start_eta) > 0)
    eta = loc2dirac(start_eta, n);
  if(length(start_eta) > n)
    stop('Error in start_eta: length of changepoint configuration cannot exceed
         length of time series.')
  eta[1:max(p, 1)] = 0; ## When p = 0, time 1 cannot be a changpoint

  ## Initial values of xi
  if(length(start_xi) == 0)
    xi = rbinom(n, 1, 0.02);
  if(length(start_xi) == n)
    xi = start_xi;
  if(length(start_xi) < n && length(start_xi) > 0)
    eta = loc2dirac(start_xi, n);
  if(length(start_xi) > n)
    stop('Error in start_xi: length of outlier indicator vector cannot exceed
         length of time series.')
  ## If a time is a changepoint, it cannot be an outlier
  xi[eta == 1] = 0;
  
  ## If not accommodating outliers
  if(detect_outliers == FALSE){
    xi = rep(0, n);
  }

  ## Initial values
  inference = fit_eta(x, A, eta, xi, p, fit, penalty, nu, kappa, a, b_eta, b_xi,
                      scale_trend_design, weights);
  best = current = list(eta = eta, xi = xi, inference = inference,
                        change_eta = FALSE, change_xi = FALSE);
  eta_mcmc[1, ] = current$eta;
  xi_mcmc[1, ] = current$xi;
  para_mcmc[1, ] = c(current$inference$bmdl, current$inference$phi,
                     current$inference$sigmasq);

  ## Start MCMC
  for(it in 1:iter){

    ## Update eta and/or xi
    if(detect_outliers){
      action = sample(c('eta_birth', 'eta_death', 'eta_swap', 'xi_birth', 
                        'xi_death', 'eta_to_xi_flip'), 1, 
                      prob = c(0.25, 0.15, 0.1, 0.25, 0.15, 0.1));
    } else {
      action = sample(c('eta_birth', 'eta_death', 'eta_swap'), 1, 
                      prob = c(0.5, 0.3, 0.2), replace = TRUE);
    }
    
    ## Metropolis-Hastings update eta and/or xi
    if(action == 'eta_birth')
      current = eta_MH_birth(x, A, current, p, fit, penalty, nu, kappa, a, b_eta,
                            b_xi, scale_trend_design, weights, max_changes);
    if(action == 'eta_death')
      current = eta_MH_death(x, A, current, p, fit, penalty, nu, kappa, a, b_eta,
                             b_xi, scale_trend_design, weights);
    if(action == 'eta_swap')
      current = eta_MH_swap(x, A, current, p, fit, penalty, nu, kappa, a, b_eta,
                            b_xi, scale_trend_design, weights);
    if(action == 'xi_birth')
      current = xi_MH_birth(x, A, current, p, fit, penalty, nu, kappa, a, b_eta,
                            b_xi, scale_trend_design, weights, max_outliers);
    if(action == 'xi_death')
      current = xi_MH_death(x, A, current, p, fit, penalty, nu, kappa, a, b_eta,
                            b_xi, scale_trend_design, weights);
    if(action == 'eta_to_xi_flip')
      current = eta_to_xi_MH_flip(x, A, current, p, fit, penalty, nu, kappa, a, 
                                  b_eta, b_xi, scale_trend_design, weights, 
                                  max_outliers);
    
    if(sum(current$eta * current$xi) > 0){
      stop('Error: a time cannot be both a changepoint and an outlier.')
    }
    
    ## try to put the new eta to map200 if good
    if(current$change_eta == TRUE || current$change_xi == TRUE){
      change_rate = change_rate + 1 / iter;

      ## Update the best model is the current model is better
      if(current$inference$bmdl < best$inference$bmdl)
        best = current;
    }

    ## Save MCMC iterations with thinning
    if(it %% thin == 0){
      eta_mcmc[it / thin + 1, ] = current$eta;
      xi_mcmc[it / thin + 1, ] = current$xi;
      para_mcmc[it / thin + 1, ] = c(current$inference$bmdl,
                                     current$inference$phi,
                                     current$inference$sigmasq);

      ## show: x0% completed
      if( (it * 10) %% iter == 0 && it != iter  && track_time == TRUE)
        cat(paste( (it * 100) / iter), '% completed...\n', sep = '');
      if( it == iter  && track_time == TRUE){
        cat(paste( (it * 100) / iter), '% completed.\n', sep = '');
        cat('\n');
      }
    }
  }

  ## Track runtime
  t.finish = proc.time();
  runtime = t.finish - t.start;
  if(track_time == TRUE){
    cat('Time used (in second): \n')
    print(runtime);
    cat('\n');
  }

  ## Save the input parameters as a record
  input_parameters = list(x = x, dates = dates, iter = iter, thin = thin, p = p,
                          time_unit = time_unit, seasonal_means = seasonal_means,
                          k = k, fit = fit, penalty = penalty, nu = nu,
                          kappa = kappa, a = a, b_eta = b_eta, b_xi = b_xi, 
                          start_eta = start_eta, start_xi = start_xi, 
                          track_time = track_time, weights = weights,
                          scale_trend_design = scale_trend_design, 
                          max_changes = max_changes, max_outliers = max_outliers
                          );

  return( list(eta_mcmc = eta_mcmc, xi_mcmc = xi_mcmc, para_mcmc = para_mcmc,
               best = best, A = A, runtime = runtime,
               input_parameters = input_parameters) );
}
