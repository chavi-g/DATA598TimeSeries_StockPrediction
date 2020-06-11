generateTs = function(timeSeriesCSV,id){
  
  dataForts = timeSeriesCSV %>% filter(stock_id == id) %>% dplyr::select(Close, dateInDateFormat)
  dataForts = dataForts[dataForts$dateInDateFormat < '2019-11-01', ]
  minmaxdateTime = dataForts %>% summarise(minDate = min(dateInDateFormat), maxDate = max(dateInDateFormat))
  
  timeSeries = ts(dataForts$Close, start = decimal_date(as.Date(minmaxdateTime$minDate))
                  , frequency = 365.25/7)
  
  return(timeSeries)
}

get_data_break = function(timeSeriesCSV, id, test_size) {
  dataForts = timeSeriesCSV %>% filter(stock_id == id) %>% dplyr::select(Close, dateInDateFormat)
  len = nrow(dataForts)
  return(nth(dataForts$dateInDateFormat, len-test_size+1))
}


get_sample_path = function(fc, errors, start_point) {
  sample_path = (fc$mean +  rnorm(length(fc$mean), mean = 0, sd = sd(errors)))
  return(ts(sample_path, start = start_point, freq=365.25/7))
}

get_predictions = function(stock_id, test_size = 14) {
  set.seed(1234)
  dataTimeSeries = read.csv('stock_series_train.csv', header = TRUE, stringsAsFactors = FALSE)
  dataTimeSeries = dataTimeSeries %>% mutate(dateInDateFormat = as.Date(Date, format="%m/%d/%Y"))
  timeSeries = generateTs(dataTimeSeries, stock_id)
  print(end(timeSeries))
  crps_results = read.csv('crps_values.csv', header = TRUE)
  
  h <- 14
  tbats_model <- tbats(timeSeries, use.box.cox = TRUE, 
                       use.trend = TRUE, 
                       use.damped.trend = TRUE)
  arima_model <- auto.arima(timeSeries, lambda=0)
  stl_model <- stlm(timeSeries, lambda=0)
  ets_model <- ets(timeSeries, lambda = 0)
  
  aic = c()
  ks = 19:26
  for (i in ks) {
    fitx = auto.arima(timeSeries, seasonal = FALSE, xreg = forecast::fourier(timeSeries, K = i))
    aic = c(aic, AIC(fitx))
  }
  bestK = ks[which(aic==max(aic))[1]]
  print(paste0("Fourier K ", bestK))
  fourier_model = auto.arima(timeSeries, seasonal = FALSE, xreg = forecast::fourier(timeSeries, K = bestK))
  bestModel = crps_results[stock_id, ]$BestModel
  
  p_mix = c()
  p_arima = c()
  p_tbats = c()
  p_stlm = c()
  p_fourier = c()
  p_ets = c()
  for (i in 1:500) {
    arima_pred = simulate(arima_model, nsim = h)
    tbats_pred = get_sample_path(forecast(tbats_model, h=h), tbats_model$errors, start(arima_pred))
    stl_pred = get_sample_path(forecast(stl_model, h=h), stl_model$residuals, start(arima_pred))
    fourier_pred = simulate(fourier_model, xreg = forecast::fourier(timeSeries, K = bestK, h = h), nsim = h)
    ets_pred = simulate(ets_model, nsim = h)
    
    if (bestModel == 'Mixture') {
      tbats_fitted = fitted(tbats_model)
      arima_fitted = fitted(arima_model)
      fourier_fitted = fitted(fourier_model)
      ets_fitted = fitted(ets_model)
      
      X = cbind(TBATS = tbats_fitted, ARIMA = arima_fitted, FOURIER = fourier_fitted, ETS = ets_fitted)
      X_new = cbind(TBATS = tbats_pred, ARIMA = arima_pred, FOURIER = fourier_pred, ETS = ets_pred)
      
      MLpol0 = mixture(Y = timeSeries, experts = X, model = "Ridge", loss.type = "square")
      p_mix[[i]] = predict(MLpol0, X_new, type = 'response', online = FALSE)
    }
    p_arima[[i]] = as.matrix(arima_pred, c(h,1))
    p_tbats[[i]] = as.matrix(tbats_pred, c(h,1))
    p_stlm[[i]] = as.matrix(stl_pred, c(h,1))
    p_fourier[[i]] = as.matrix(fourier_pred, c(h,1))
    p_ets[[i]] = as.matrix(ets_pred, c(h, 1))
  }
  
  if (bestModel == 'Arima') {
    preds = do.call(cbind, p_arima)
  } else if (bestModel == 'Mixture') {
    preds = do.call(cbind, p_mix)
  } else if (bestModel == 'Stlm') {
    preds = do.call(cbind, p_stlm)
  } else if (bestModel == 'Tbats') {
    preds = do.call(cbind, p_tbats)
  } else if (bestModel == 'Fourier') {
    preds = do.call(cbind, p_fourier)
  } else if (bestModel == 'ETS') {
    preds = do.call(cbind, p_ets)
  }
  
  dates_column = c('11/1/2019', '11/8/2019', '11/15/2019', '11/22/2019', '11/29/2019', '12/6/2019', '12/13/2019', '12/20/2019', '12/27/2019', '1/3/2020', '1/10/2020', '1/17/2020', '1/24/2020', '1/31/2020')
  stock_column = rep(stock_id, h)
  pred_colnames = c()
  for (j in 1:500) {
    pred_colnames[j] = paste0("ysim_", j)
  }

  df = cbind(stock_column, date, preds)
  colnames(df) = c("stock_id", "Date", pred_colnames)
  return(df)
}

final_predictions = data.frame()

for (i in 1:19) {
  print(paste0("Generating results for stock: ", i))
  final_predictions = rbind(final_predictions, get_predictions(i))
  print(dim(final_predictions))
}
#dev.off()

#min.col <- function(m, ...) max.col(-m, ...)
#res2['BestModel'] = colnames(results[c("Mixture", "Arima", "Tbats", "Stlm", "Fourier", "ETS")])[min.col(results[c("Mixture", "Arima", "Tbats", "Stlm")],ties.method="first")]
#print(res2)
#write.csv(res2, 'crps_values3.csv', row.names = FALSE)



