generateTs = function(timeSeriesCSV,id){
  
  dataForts = timeSeriesCSV %>% filter(stock_id == id) %>% dplyr::select(Close, dateInDateFormat)
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
  return(ts(sample_path, start = start_point , freq=365.25/7))
}

pdf("fit_plots2.pdf")
get_crps = function(stock_id, test_size = 14) {
  set.seed(1234)
  dataTimeSeries = read.csv('stock_series_train.csv', header = TRUE, stringsAsFactors = FALSE)
  dataTimeSeries = dataTimeSeries %>% mutate(dateInDateFormat = as.Date(Date, format="%m/%d/%Y"))
  timeSeries = generateTs(dataTimeSeries, stock_id)
  
  data_break = get_data_break(dataTimeSeries, stock_id, test_size)
  train <- window(timeSeries, end=decimal_date(as.Date(data_break)))
  test <- window(timeSeries, start=decimal_date(as.Date(data_break)))
  #print(length(test))
  h <- length(test)
  tbats_model <- tbats(train, use.box.cox = TRUE, 
                       use.trend = TRUE, 
                       use.damped.trend = TRUE)
  arima_model <- auto.arima(train, lambda=0)
  stl_model <- stlm(train, lambda=0)
  ets_model <- ets(train, lambda = 0)
  
  aic = c()
  ks = 19:26
  for (i in ks) {
    fitx = auto.arima(train, seasonal = FALSE, xreg = forecast::fourier(train, K = i))
    aic = c(aic, AIC(fitx))
  }
  bestK = ks[which(aic==max(aic))[1]]
  print(paste0("Fourier K ", bestK))
  fourier_model = auto.arima(train, seasonal = FALSE, xreg = forecast::fourier(train, K = bestK))
  
  p_mix = c()
  p_arima = c()
  p_tbats = c()
  p_stlm = c()
  p_fourier = c()
  p_ets = c()
  for (i in 1:500) {
    tbats_pred = get_sample_path(forecast(tbats_model, h=h), tbats_model$errors, start(test))
    arima_pred = simulate(arima_model, nsim = h)
    stl_pred = get_sample_path(forecast(stl_model, h=h), stl_model$residuals, start(test))
    fourier_pred = simulate(fourier_model, xreg = forecast::fourier(train, K = bestK, h = h), nsim = h)
    ets_pred = simulate(ets_model, nsim = h)
    
    tbats_fitted = fitted(tbats_model)
    arima_fitted = fitted(arima_model)
    fourier_fitted = fitted(fourier_model)
    ets_fitted = fitted(ets_model)
    
    X = cbind(TBATS = tbats_fitted, ARIMA = arima_fitted, FOURIER = fourier_fitted, ETS = ets_fitted)
    X_new = cbind(TBATS = tbats_pred, ARIMA = arima_pred, FOURIER = fourier_pred, ETS = ets_pred)
    
    MLpol0 = mixture(Y = train, experts = X, model = "Ridge", loss.type = "square")
    p_mix[[i]] = predict(MLpol0, X_new, type = 'response', online = FALSE)
    p_arima[[i]] = as.matrix(arima_pred, c(h,1))
    p_tbats[[i]] = as.matrix(tbats_pred, c(h,1))
    p_stlm[[i]] = as.matrix(stl_pred, c(h,1))
    p_fourier[[i]] = as.matrix(fourier_pred, c(h,1))
    p_ets[[i]] = as.matrix(ets_pred, c(h, 1))
  }
  
  preds_mix = do.call(cbind, p_mix)
  preds_arima = do.call(cbind, p_arima)
  preds_tbats = do.call(cbind, p_tbats)
  preds_stl = do.call(cbind, p_stlm)
  preds_fourier = do.call(cbind, p_fourier)
  preds_ets = do.call(cbind, p_ets)
  
  crps_mix = mean(crps_sample(as.numeric(test), preds_mix))
  crps_arima = mean(crps_sample(as.numeric(test), preds_arima))
  crps_tbats = mean(crps_sample(as.numeric(test), preds_tbats))
  crps_stlm = mean(crps_sample(as.numeric(test), preds_stl))
  crps_fourier = mean(crps_sample(as.numeric(test), preds_fourier))
  crps_ets = mean(crps_sample(as.numeric(test), preds_ets))
  
  if (crps_mix < crps_arima) {
    df_sample = cbind(timeSeries, ts(p_mix[[i]],  start=start(test), freq=365.25/7))
    colnames(df_sample) <- c("Data","Mixture")
    p = autoplot(df_sample) + 
      xlab("Year") + ylab("Time Series")  + ggtitle(paste0("Mixture Model fit for Stock ", stock_id))
  } else {
    df_sample = cbind(timeSeries, ts(p_arima[[i]],  start=start(test), freq=365.25/7))
    colnames(df_sample) <- c("Data","Arima")
    p = autoplot(df_sample) +
      xlab("Year") + ylab("Time Series")  + ggtitle(paste0("Arima fit for Stock ", stock_id))
  }
  print(p)
  
  print(data.frame(Stock = stock_id, Mixture = crps_mix, Arima = crps_arima, Tbats = crps_tbats, Stlm = crps_stlm, Fourier = crps_fourier, ETS = crps_ets))
  return(data.frame(Stock = stock_id, Mixture = crps_mix, Arima = crps_arima, Tbats = crps_tbats, Stlm = crps_stlm, Fourier = crps_fourier, ETS = crps_ets))
}

results = data.frame(Stock = numeric(0), Mixture = numeric(0), Arima = numeric(0), Tbats = numeric(0), Stlm = numeric(0))

for (i in 1:19) {
  print(paste0("Generating results for stock: ", i))
  results = rbind(results, get_crps(i))
}
dev.off()

min.col <- function(m, ...) max.col(-m, ...)
results['BestModel'] = colnames(results[c("Mixture", "Arima", "Tbats", "Stlm", "Fourier", "ETS")])[min.col(results[c("Mixture", "Arima", "Tbats", "Stlm")],ties.method="first")]
print(results)
write.csv(results, 'crps_values.csv', row.names = FALSE)

