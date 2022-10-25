
shuffing <- function(data, data_prsim, marginal="empirical") {

  ### vi) reorder the surrogate to the distribution of the original time series
  ### apply daily backtransformation: ensures smoothness of regimes

  #reorder the observation at a day according to the surrogates order
  data_new <- data.frame("seasonal"=data_prsim)
  data_new$rank <- rank(data_new$seasonal)
  data_new$index <- data$index
  data_new$simulated_seasonal <- NA

  d_vec <- unique(data$index)

  # for(di in d_vec[2:(length(d_vec)-1)]){
  #   d <- c(di-1,di,di+1)
  for(d in d_vec){
    #cat(data[which(data$index%in%c(d)),])
    data_day <- data[which(data$index%in%c(d)),]
    #cat(summary(data_day))
    ### use empirical distribution for backtransformation
    if(marginal=="empirical"){
      data_day$rank <- rank(data_day$bcc) # rank of the day in obs

      data_new$rank[which(data$index%in%c(d))] <- rank(data_new[which(data$index%in%c(d)),]$seasonal)

      ### derive corresponding values from the empirical distribution
      ### identify value corresponding to rank in the original time series
      data_ordered <- data_day[order(data_day$rank),] #order of rank from small to high
      data_new$simulated_seasonal[which(data_new$index%in%c(d))] <- data_ordered$bcc[data_new$rank[which(data$index%in%c(d))]]
    }

  }  # end for loop

  # ts.plot(cbind(data_new$simulated_seasonal, data$obs), col=1:4)
  # ts.plot(cbind(data_new$seasonal, data$mod), col=1:4)

  return(data_new$simulated_seasonal)

}
