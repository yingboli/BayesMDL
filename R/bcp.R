#' Fit the changepoint detection method to multiple variables
#'
#' Apply Bayesian Minimum Description Length to detect multiple changepoints,
#' individually for a group of time series data (one for each variable to be
#' monitored). Returns the optimal changepoint models, and graphic outputs as a
#' pdf file.
#'
#' @inheritParams bmdl
#' @param X Variables to be monitored, a \code{data.frame} of numerics, with
#'   each row being a variable, in the format of a time series of length
#'   \code{n}. If such information is in a CSV file, \code{X} can also be the
#'   name of this file, with extension name, and a path if needed.
#' @param pdf_file Name of the output pdf file with extension name, a character;
#'   can include path.
#' @param select The selected model to plot, \code{'MAP'} for maximum a
#'   posterior, or \code{'BMA'} for median probability model under Bayesian
#'   model averaging.
#' @param width,height Width and height of the output pdf (per page).
#' @param mar Margin of the output pdf (per page).
#' @param mfrow Number of rows and columns of figures on per pdf page.
#'
#' @return
#'   \item{eta_all}{Estimated changepoints for all variables; a matrix of 0/1
#'     indicators, where each row is for one variable, in the same order as
#'     the input data \code{X}.}
#'   \item{seasonal_cycles_all}{Estimated seasonal means for all variables,
#'     where each row is for one variable.}
#'   \item{mean_and_trend_all}{Estimated linear trends for all variables,
#'     where each row is for one variable.}
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom utils read.csv
#' @export
#' @keywords internal

bcp = function(X, pdf_file = 'bcp.pdf', select = 'MAP', iter = 1e4,
               thin = max(1, iter / 1e3), p = 2, time_unit = 'month',
               seasonal_means = 'harmonic', k = 3, scale_trend_design = 0.05,
               fit = 'marlik', penalty = 'bmdl', nu = 5, a = 1, b = 1,
               width = 8, height = 4, mar = c(5, 5, 4, 1), mfrow = c(1, 1),
               start_eta = NULL, track_time = TRUE){

  ## Read in X if it's in a csv file
  if(is.character(X) == TRUE){
    X = read.csv(file = X);

    rownames(X) = as.character(X[, 1]);
    X = X[, -1];
    colnames(X) = gsub('X', '',  colnames(X));
  }

  ## Compute length of the time series n, and the total number of variables
  n = ncol(X);
  nvars = nrow(X);

  ## Dates and vars
  dates = strptime(paste(colnames(X), '01', sep = ''), format = '%Y%m%d');
  vars = rownames(X);

  ## Matrices to save eta outputs
  eta_all = seasonal_cycles_all = mean_and_trend_all = NULL;

  ## Apply bmdl method
  for(i in 1:nvars){
    x = as.matrix(X)[i, ];

    results = bmdl(x, dates, iter = iter, thin = thin, p = p,
                   time_unit = time_unit, seasonal_means = seasonal_means, k = k,
                   scale_trend_design = scale_trend_design, fit = fit,
                   penalty = penalty, nu = nu, a = a, b = b,
                   start_eta = start_eta, track_time = track_time);
    A = results$A;

    ## Compute eta and current
    if(select == 'MAP'){
      current = results$best;
      eta = current$eta;
    }
    if(select == 'BMA'){
      iter_save = nrow(results$eta_mcmc);
      keep = ceiling(iter_save / 2):iter_save;
      eta = as.numeric(apply(results$eta_mcmc[keep, 1:n], 2, mean) > 0.5);
      current = list(eta = eta, inference =
                       fit_eta(x, A, eta, p = p, fit = fit, penalty = penalty,
                               nu = nu, a = a, b = b,
                               scale_trend_design = scale_trend_design));
    }

    ## Compute seasonal_cycles and mean_and_trend
    seasonal_cycles = c(A %*% current$inference$s);
    mu = current$inference$mu;
    D = D_eta(x, eta, scale_trend_design)$D;
    mean_and_trend = c(D %*% mu);


    ## Save all
    eta_all = rbind(eta_all, eta);
    seasonal_cycles_all = rbind(seasonal_cycles_all, seasonal_cycles);
    mean_and_trend_all = rbind(mean_and_trend_all, mean_and_trend);
  }

  ## Create the pdf file
  pdf(file = pdf_file, width = width, height = height, onefile = TRUE);

  par(mfrow = mfrow, mar = mar);
  plot_bar_heat(eta_all, dates);

  lcp = apply(eta_all, 1, last_cp);
  plot_order = order(lcp, decreasing = TRUE);

  for(i in plot_order){
    x = as.matrix(X)[i, ];

    plot_eta(x = as.matrix(X)[i, ], dates, varname = vars[i], eta = eta_all[i, ],
             seasonal_cycles = seasonal_cycles_all[i, ],
             mean_and_trend = mean_and_trend_all[i, ]);
  }
  dev.off();

  return(list(eta_all = eta_all, seasonal_cycles_all = seasonal_cycles_all,
              mean_and_trend_all = mean_and_trend_all));
}
