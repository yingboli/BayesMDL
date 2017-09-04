###################################################################
######    Bayesian MDL for multiple changepoint detection   #######
###################################################################
#' Plot the time series data and the changepoints
#'
#' Create a plot, with time series data \code{x} shown in points, linear trends
#' under the selected changepoint model in blue lines, estimated changepoints
#' in verticle dashed blue lines, and fitted value (including seasonal means
#' and linear trends) in a purple line. Optional residual plots can also be
#' plotted.
#'
#' Two sets of inputs can be used:
#' \enumerate{
#'   \item Either \code{current}, \code{A}, \code{scale_trend_design},
#'   \item Or \code{eta}, \code{seasonal_cycles}, \code{mean_and_trend}.
#' }
#' When used, every inputs in that set must but non-NULL, while all inputs in
#' the other set must be NULL.
#'
#' @inheritParams bmdl
#' @inheritParams fit_eta
#' @inheritParams eta_MH_flip
#' @param varname The name of the variable \code{x}, a character.
#' @param seasonal_cycles Estimated seasonal means, a numeric vector of length
#'   \code{n}.
#' @param mean_and_trend Estimated piecewise linear fits, a numeric vector of
#'   length \code{n}.
#' @param plot_res Logical, whether residual plots will be added after the main
#'   plots. If \code{TRUE}, then the following four plots will be added:
#'   scatterplot and acf plot of the linear model residuals, and scatterplot
#'   and acf plot fo the white noise estimates.
#' @param yw_phi Yule walker estimator of AR coefficients, a vector of length
#'   \code{p}; must be non-NULL if \code{plot_res == TRUE} and
#'   \code{is.null(current)}.
#' @param ... Other arguments for the \code{\link{plot}} function.
#' @importFrom graphics abline axis box lines mtext plot
#' @importFrom stats acf
#' @export
#' @keywords internal
#'

plot_eta = function(x, dates, varname = NULL, current = NULL, A = NULL,
                    scale_trend_design = NULL, eta = NULL,
                    seasonal_cycles = NULL, mean_and_trend = NULL,
                    plot_res = FALSE, yw_phi = NULL, ...){

  ## Based on dates, creates some time related objects
  time_ind = (dates$year - dates$year[1]) * 12 + (dates$mon - dates$mon[1]) + 1;
  year_ind = dates$year + 1900;
  jan_ind = which(dates$mon == 0);
  month_number = c(paste0('0', 1:9), '10', '11', '12');

  ## Compute seasonal_cycles, mean_and_trend, and eta if not provided
  if(all(!is.null(current), !is.null(A), !is.null(scale_trend_design))){
    eta = current$eta;
    seasonal_cycles = c(A %*% current$inference$s);
    mu = current$inference$mu;
    D = D_eta(x, eta, scale_trend_design)$D;
    mean_and_trend = c(D %*% mu);
  }

  if(any(is.null(eta), is.null(seasonal_cycles), is.null(mean_and_trend))){
    stop('Missing inputs. See the Details section of the help file.');
    return();
  }

  regime = cumsum(eta) + 1;

  ## Title of the main plot
  if(is.null(varname)){
    title1 = paste0(dates[1]$year + 1900, '/',
                    month_number[dates[1]$mon+1], ' to ',
                    dates[length(dates)]$year + 1900, '/',
                    month_number[dates[length(dates)]$mon+1]);
  } else {
    title1 = paste0(varname, ': from ', dates[1]$year + 1900, '/',
                    month_number[dates[1]$mon+1], ' to ',
                    dates[length(dates)]$year + 1900, '/',
                    month_number[dates[length(dates)]$mon+1]);
  }

  ## Create the main plot
  plot(time_ind, x, xlab = '', ylab = '', axes = FALSE, main = title1, ...);
  axis(1, at = time_ind[jan_ind], labels = year_ind[jan_ind]);
  axis(2); box();
  lines(time_ind, seasonal_cycles + mean_and_trend, col = 6, lty = 1);
  for(r in 1:(sum(eta) + 1)){
    lines(time_ind[regime == r],
          (mean_and_trend + mean(seasonal_cycles))[regime == r],
          col = 4, lwd = 2);
    if(r < (sum(eta) + 1)){
      tau = which(eta == 1)[r];
      mtext(paste(year_ind[tau], month.abb[dates$mon[tau] + 1]), at = tau,
            cex = 1, col = 4)
    }
  }
  abline(v = time_ind[which(eta == 1)], col = 4, lty = 2, lwd = 2);

  ## Create the residual plots
  if(plot_res){
    ## Linear model residuals
    res = x - seasonal_cycles - mean_and_trend;

    ## Scatterplot for residuals
    plot(time_ind, res, xlab = '', ylab = '', main = 'Linear model residuals',
         axes = FALSE, ...)
    axis(1, at = time_ind[jan_ind], labels = year_ind[jan_ind]);
    axis(2); box();
    abline(h = 0, lty = 2);

    ## ACF plot for residuals
    acf(res, main = 'ACF: linear model residuals');

    ## AR white noises
    if(!is.null(current)){
      yw_phi = current$inference$phi;
    }
    p = length(yw_phi);
    n = length(x);
    ar_res = res[(p + 1):n];
    if(p > 0){
      for(j in 1:p){
        ar_res = ar_res - yw_phi[j] * res[((p + 1):n) - j]
      }
    }
    ar_res_na = c(rep(NA, p), ar_res);

    ## Scatterplot for white noises
    plot(time_ind, ar_res_na, xlab = '', ylab = '', main = 'AR white noises',
         axes = FALSE, ...);
    axis(1, at = time_ind[jan_ind], labels = year_ind[jan_ind]);
    axis(2); box();
    abline(h = 0, lty = 2);

    ## ACF plot for white noises
    acf(ar_res, main = 'ACF: white noises');
  }

}

###################################################################
#' Plot eta of the MAP or BMA model of a MCMC fit
#' @inheritParams plot_eta
#' @param results Output object of the function \code{\link{bmdl}}.
#' @param select The selected model to plot, \code{'MAP'} for maximum a
#'   posterior, or \code{'BMA'} for median probability model under Bayesian
#'   model averaging.
#' @export
#' @keywords internal
#'

plot_best = function(results, select = 'MAP', varname = NULL, plot_res = FALSE,
                     ...){

  x = results$input_parameters$x;
  n = length(x);
  dates = results$input_parameters$dates;
  A = results$A;
  scale_trend_design = results$input_parameters$scale_trend_design;

  if(select == 'MAP'){
    current = results$best;
  }
  if(select == 'BMA'){
    results$eta_mcmc;
    iter_save = nrow(results$eta_mcmc);
    keep = ceiling(iter_save / 2):iter_save;
    eta_bma = as.numeric(apply(results$eta_mcmc[keep, 1:n], 2, mean) > 0.5);
    xi_bma = as.numeric(apply(results$xi_mcmc[keep, 1:n], 2, mean) > 0.5);
    
    p = results$input_parameters$p;
    fit = results$input_parameters$fit;
    penalty = results$input_parameters$penalty;
    nu = results$input_parameters$nu;
    a = results$input_parameters$a;
    b_eta = results$input_parameters$b_eta;
    b_xi = results$input_parameters$b_xi;
    
    current = list(eta = eta_bma, inference =
                     fit_eta(x, A, eta_bma, xi_bma, p = p, fit = fit, penalty = penalty,
                             nu = nu, a = a, b_eta = b_eta, b_xi = b_xi, 
                             scale_trend_design = scale_trend_design));
  }


  plot_eta(x, dates, varname = varname, current = current, A = A,
           scale_trend_design = scale_trend_design, plot_res = plot_res, ...);

}

###################################################################
#' Plot the residual
#' @inheritParams plot_eta
#' @inheritParams bmdl
#' @export
#' @keywords internal

plot_res = function(current, x, dates, A, scale_trend_design, p,  ...){

  n = length(x);

  time_ind = (dates$year - dates$year[1]) * 12 +
    (dates$mon - dates$mon[1]) + 1;
  year_ind = dates$year + 1900;
  jan_ind = which(dates$mon == 0)

  seasonal_cycles = A %*% current$inference$s;
  mu = current$inference$mu;
  tmp = D_eta(x, current$eta, scale_trend_design);
  D = tmp$D;
  seg_lengths = tmp$seg_lengths;
  mean_and_trend = D %*% mu;

  regime = cumsum(current$eta) + 1;

  res = x - seasonal_cycles - mean_and_trend;
  ar_res = res[(p + 1):n];
  if(p > 0){
    for(j in 1:p){
      ar_res = ar_res - current$inference$phi[j] * res[((p + 1):n) - j]
    }
    plot(time_ind[-(1:p)], ar_res, xlab = '', axes = FALSE, ...);
  }
  if(p == 0){
    plot(time_ind, ar_res, xlab = '', axes = FALSE, ...);
  }
  axis(1, at = time_ind[jan_ind - p], labels = year_ind[jan_ind]);
  axis(2); box();
  abline(h = 0, col = 2, lty = 2);

  acf(ar_res, main = '');

}

###################################################################
#' Plot the heatmap and barplot of changepoints among multiple variables
#'
#' @inheritParams plot_eta
#' @param best_all A matrix of 0/1 indicators on the selected changepoint models.
#'   Each row is the \code{eta} of a variable.
#' @param type Type of plot: \code{'heatmap'}, \code{'barplot'}, or \code{'both'}.
#'
#' @importFrom grDevices gray
#' @importFrom graphics barplot image
#' @export
#' @keywords internal
#'
plot_bar_heat = function(best_all, dates, type = 'both'){
  n = ncol(best_all);

  month_label = paste((dates$year + 1900), month.abb[dates$mon + 1]);
  ## location of x in the heatmap
  month_at = (1:(n - 1) - 0.5) / (n - 1);
  month_at = c(month_at - month_at[1], 1);

  ## barplot
  if(type == 'barplot' | type == 'both'){
    month_label_tmp = month_label;
    month_label_tmp[-seq(1, n, by = 3)] = '';

    barplot(apply(best_all, 2, sum), names.arg = month_label_tmp, las = 2,
            width = 1/n * 0.8, space = 0.25, xlim = c(0, 1), xaxs = 'i');
  }
  ## heatmap
  if(type == 'heatmap' | type == 'both'){
    image(t(best_all), xlab = '', ylab = 'variable', axes = FALSE,
          col = gray((32:0)/32),
          xlim = c(0 - month_at[2] / 2, 1 + month_at[2] / 2));
    axis(1, at = month_at[seq(1, n, by = 3)],
         labels = month_label[seq(1, n, by = 3)],
         las = 2, cex.lab = 0.5);
    box();
  }

}

