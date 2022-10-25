#' Title
#'
#' @param data
#' @param subset
#' @param method
#' @param variable
#' @param theta
#' @param PRand
#' @param number_sim
#' @param wavelet
#' @param dt
#' @param dj
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
bc.cwt <- function(data, subset=1:256, method=c("additive","multiplicative"), variable, theta=0.1,
                   QM=TRUE,PRand=FALSE,
                   number_sim=1, wavelet="morlet", dt=1, dj=1/8, upperPeriod, seed=NULL)
  {

  # data <- data.h.df; seed=101; wavelet="morlet"; dt=1; dj=1/8; number_sim=1
  # method <- switch(1, "additive","multiplicative")
  ###===============================###===============================###
  ### data check
  #obs and mod have same length
  if(!is.null(subset)){
    for(l in 1:length(data))data[[l]] <- data[[l]][subset, ]
  }

  ###===============================###===============================###
  ### Generation of white noise for random phases generation
  ### generate random sample of indices for each simulation run
  if(!is.null(seed)) set.seed(seed)
  noise_mat <- list()
  for (r in 1:number_sim){
    ts_wn <- rnorm(n=length(data[[1]]$mod), mean = 0, sd = 1) ### iid time series
    #wt_noise <- wmtsa::wavCWT(x=ts_wn,wavelet=wavelet,n.scale=n_wave)
    wt_noise <- t(WaveletComp::WaveletTransform(x=ts_wn,dt=dt,dj=dj)$Wave)

    noise_mat[[r]] <- as.matrix(wt_noise)
  }

  ### fitting distribution to all stations for which simulations are to be derived
  marginal_list<- list()
  for(l in seq_along(data)) marginal_list[[l]]<-'empirical'

  ###===============================###===============================###
  ###pre-process: bias correction
  ### center_data: subtract mean from values
  for(l in 1:length(data)){
    data[[l]]$norm.o <- data[[l]]$obs-mean(data[[l]]$obs,na.rm=TRUE)
    data[[l]]$norm.m <- data[[l]]$mod-mean(data[[l]]$mod,na.rm=TRUE)
  }

  # spec.pgram(data[[l]]$norm.o)
  # spec.pgram(data[[l]]$norm.o[subset])
  # ts.plot(cbind(data[[l]]$norm.o,data[[l]]$norm.m), col=1:2)

  ### run through all stations
  out_list<-list()
  for(l in 1:length(data)){
    ### list for storing results
    data_sim <- NULL
    data[[l]]$index <- as.numeric(format(data[[l]]$Date,format='%j'))

    ###===============================###===============================###
    ### i) use continuous wavelet transform (CWT) to wavelet transform the data
    # wt_o <- wavCWT(x=data[[l]]$norm.o,wavelet=wavelet,n.scale=n_wave)
    # wt_m <- wavCWT(x=data[[l]]$norm.m,wavelet=wavelet,n.scale=n_wave)

    wt_o <- t(WaveletComp::WaveletTransform(x=data[[l]]$norm.o,dt=dt,dj=dj,upperPeriod = upperPeriod)$Wave)
    wt_m <- t(WaveletComp::WaveletTransform(x=data[[l]]$norm.m,dt=dt,dj=dj,upperPeriod = upperPeriod)$Wave)

    ### return CWT coefficients as a complex matrix with rows and columns representing times and scales, respectively.
    wt_o_mat <- as.matrix(wt_o)
    real.o <- Re(wt_o_mat)
    modulus.o <- Mod(wt_o_mat)       ### derive modulus of complex numbers (radius)
    phases.o <- Arg(wt_o_mat)        ### extract phases (argument)

    wt_m_mat <- as.matrix(wt_m)
    real.m <- Re(wt_m_mat)
    modulus.m <- Mod(wt_m_mat)
    phases.m <- Arg(wt_m_mat)

    ###===============================###===============================###
    ### bias correction
    if(method=="additive"){

      factor.amp <- modulus.o-modulus.m
      factor.arg <- phases.o-phases.m
      phases.c <- phases.m+factor.arg
      modulus.c <- modulus.m+factor.amp

    } else if(method=="multiplicative"){

      factor.amp <- modulus.o/modulus.m
      factor.arg <- phases.o/phases.m
      phases.c <- phases.m*factor.arg
      modulus.c <- modulus.m*factor.amp

    }

    if(QM){
      #modulus.tmp <- MBCr(o.c=modulus.o, m.c=modulus.m, m.p=modulus.p)
      #modulus.tmp <- do.call(paste0("MBC","p"),list(o.c=modulus.o, m.c=modulus.m,m.p=modulus.m,silent=TRUE))
      modulus.tmp <- do.call(paste0("MRS"),list(o.c=modulus.o, m.c=modulus.m,m.p=modulus.m))

      modulus.c <- modulus.tmp$mhat.c
      #modulus.p <- modulus.tmp$mhat.p
      #phases.c <- phases.m
    }

    #plot(modulus.o, modulus.m); abline(a=0, b=1)

    #cat(dim(phases.o))

    # for(i in sample(1:ncol(phases.o),3)) {plot(phases.o[,i], phases.m[,i]); abline(a=0, b=1)}
    # plot(phases.o-phases.m); abline(a=0, b=0)

    if(FALSE){
      if(FALSE){
        ###amplitude correction first
        mat_new <- matrix(complex(modulus=modulus.c,argument=phases.m),ncol=ncol(phases.m))
        rec <- fun_icwt(x=mat_new,dt=dt,dj=dj)

        rec_cwt <- t(WaveletComp::WaveletTransform(x=rec, dt=dt,dj=dj)$Wave)
        factor.arg <- phases.o - Arg(rec_cwt)
        #factor.mod <- modulus.o - Mod(rec_cwt)
        #cat(sum(abs(factor.arg1-factor.arg)))
        #cat(sum(abs(factor.mod1)))

      } else {
        ###phase correction first
        mat_new <- matrix(complex(modulus=modulus.m,argument=phases.c),ncol=ncol(phases.c))
        rec <- fun_icwt(x=mat_new,dt=dt,dj=dj)

        rec_cwt <- t(WaveletComp::WaveletTransform(x=rec, dt=dt,dj=dj)$Wave)
        factor.amp <- modulus.o - Mod(rec_cwt)
      }
    }

    ###amplitude and phase correction
    mat_new <- matrix(complex(modulus=modulus.c,argument=phases.c),ncol=ncol(phases.o))
    #cat(sum(abs(wt_o_mat-mat_new)))

    ###inverse
    rec<- fun_icwt(x=mat_new,dt=dt,dj=dj)

    # wt_m.bc <- analyze.wavelet(data.frame(x=data[[l]]$norm.m), "x",
    #                            loess.span = 0,
    #                            dt = dt, dj = dj, make.pval = TRUE, n.sim = 10)
    # cat(dim(wt_m.bc$Wave))
    # wt_m.bc$Wave <- t(mat_new)
    # rec1 <- reconstruct(wt_m.bc, only.sig = FALSE, plot.rec = FALSE, verbose = FALSE,
    #                     rescale = FALSE)$series$x.r

    if(!missing(variable)) if(variable=="prep") rec[rec<theta] <- 0

    #ts.plot(cbind(data[[l]]$norm.o,rec,rec1), col=1:3, lwd=c(2,2,1))

    rec_random <- (rec - mean(rec))*sd(data[[l]]$mod)/sd(rec) + mean(data[[l]]$mod)

    #cat(quantile(rec_random,0.5))
    #if(!missing(variable)) if(variable=="prep") rec_random[rec_random<theta] <- 0

    #ts.plot(cbind(data[[l]]$obs,rec_random), col=1:2)


    ###===============================###===============================###
    ###post-process: bias correction
    sd.ratio <- sd(data[[l]]$obs)/sd(data[[l]]$mod)
    mean.d <- mean(data[[l]]$obs) - mean(data[[l]]$mod)
    rec_random <- (rec_random - mean(rec_random)) * sd.ratio  + mean(data[[l]]$obs)

    if(!missing(variable)) if(variable=="prep") rec_random[rec_random<theta] <- 0

    data[[l]]$bcc <- rec_random

    #data[[l]]$bcc <- rec

    ###===============================###===============================###
    ###phase randomization
    if(PRand) {
      phases.rand <- lapply(1:number_sim, function(r)  Arg(noise_mat[[r]]))
      data_prsim <- lapply(1:number_sim, function(r) prsim(modulus.c, phases.rand[[r]]))

      data_sim <- lapply(1:number_sim, function(r) shuffing(data, data_prsim[[r]]))
      names(data_sim) <- paste0("r",seq(1:number_sim))
    }


    ### output
    #data_stoch <- data.frame(data[[l]][,c("Date","obs","mod","bcc")])

    out_list[[l]] <- list(data=data[[l]], simulation=data_sim,
                          method=method,
                          factor.amp=factor.amp, factor.arg=factor.arg,
                          factor.sd=sd.ratio, factor.mu=mean.d)
  }

  return(out_list)

}
