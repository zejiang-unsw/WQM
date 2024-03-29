#' CWT based quantile mapping
#'
#' @param data a list of input dataset
#' @param subset a index of number denoting the subset for calibration
#' @param variable a character string denoting the type of variable.
#' @param theta threshold of rainfall.
#' @param QM a character string denoting the qm method used.
#' @param number_sim The total number of realizations.
#' @param wavelet a character string denoting the wavelet filter to use in calculating the CWT.
#' @param dt sampling resolution in the time domain.
#' @param dj sampling resolution in the frequency domain.
#' @param method Shuffling method, M1: non-shuffling and M2: shuffling. M2 by default.
#' @param block Block size.
#' @param seed Seed for shuffling process.
#' @param PR.cal Logical value for phase randomization of calibration.
#' @param do.plot Logical value for ploting.
#' @param ... Additional arguments for QDM.
#'
#' @importFrom graphics abline axis hist legend lines par
#' @importFrom stats fft rnorm
#' @import MBC
#' @import ggplot2
#'
#' @return a list of post-processed data
#' @export
#'
bc_cwt <- function(data, subset, variable, theta=0.1, QM=c("MBC","MRS","QDM"),
                   number_sim=5, wavelet="morlet", dt=1, dj=1,
                   method="M2", block=3, seed=NULL,
                   PR.cal=FALSE, do.plot=FALSE,...)
  {

  flag.wav <- switch(1, "wmtsa", "WaveletComp")
  if(flag.wav=="wmtsa") J <- fun_cwt_J(length(data[[1]]$obs[-subset]), dt, dj) + 1

  ###=================================###====================================###
  ## white noise ----
  # generation of white noise for random phases generation
  #cmat <- cor(sapply(data, function(ls) ls$obs[subset]))
  #cat(cmat)

  noise_mat_cal <- list()
  #ts_wn_mat <- faux::rnorm_multi(length(data[[1]]$obs[subset]), length(data), mu=0, sd=1, r=cmat)
  for (r in 1:number_sim) {
    ts_wn <- rnorm(n=length(data[[1]]$obs[subset]), mean = 0, sd = 1)
    #data.obs <- as.vector(sapply(data, function(ls)ls$obs[subset]))
    #ts_wn <- sample(data.obs, size=length(data[[1]]$obs[subset]), replace=TRUE)

    if(flag.wav=="WaveletComp"){
      wt_noise <- t(WaveletComp::WaveletTransform(x=ts_wn,dt=dt,dj=dj)$Wave)
    } else if(flag.wav=="wmtsa"){
      wt_noise <- wmtsa::wavCWT(x=ts_wn,wavelet=wavelet,n.scale=J)
    }

    noise_mat_cal[[r]] <- as.matrix(wt_noise)
  }

  noise_mat_val <- list()
  #ts_wn_mat <- faux::rnorm_multi(length(data[[1]]$obs[-subset]), length(data), mu=0, sd=1, r=cmat)
  for (r in 1:number_sim) {
    ts_wn <- rnorm(n=length(data[[1]]$obs[-subset]), mean = 0, sd = 1)
    #data.obs <- as.vector(sapply(data, function(ls)ls$obs[subset]))
    #ts_wn <- sample(data.obs, size=length(data[[1]]$obs[-subset]), replace=TRUE)

    if(flag.wav=="WaveletComp"){
      wt_noise <- t(WaveletComp::WaveletTransform(x=ts_wn,dt=dt,dj=dj)$Wave)
    } else if(flag.wav=="wmtsa"){
      wt_noise <- wmtsa::wavCWT(x=ts_wn,wavelet=wavelet,n.scale=J)
    }

    noise_mat_val[[r]] <- as.matrix(wt_noise)
  }

  if(ncol(wt_noise)!=J) message(paste0("No of level: ",ncol(wt_noise)))

  ###=================================###====================================###
  # postprocessing ----
  out_list<-list()
  for(l in 1:length(data)){ # run through all stations
    ## cwt decomposition ----
    # use continuous wavelet transform (CWT) to wavelet transform the data
    if(flag.wav=="wmtsa"){
      wt_o <- wmtsa::wavCWT(x=data[[l]]$obs[subset],wavelet=wavelet,n.scale=J)
      wt_m <- wmtsa::wavCWT(x=data[[l]]$mod[subset],wavelet=wavelet,n.scale=J)
      wt_p <- wmtsa::wavCWT(x=data[[l]]$mod[-subset],wavelet=wavelet,n.scale=J)

      scale <- attr(wt_o,'scale')
    } else if(flag.wav=="WaveletComp"){
      wt_o <- t(WaveletComp::WaveletTransform(x=data[[l]]$obs[subset],dt=dt,dj=dj)$Wave)
      wt_m <- t(WaveletComp::WaveletTransform(x=data[[l]]$mod[subset],dt=dt,dj=dj)$Wave)
      wt_p <- t(WaveletComp::WaveletTransform(x=data[[l]]$mod[-subset],dt=dt,dj=dj)$Wave)

      scale <- NULL
    }

    # return CWT coefficients as a complex matrix with rows and columns
    # representing times and scales, respectively.
    wt_o_mat <- as.matrix(wt_o)
    modulus.o <- Mod(wt_o_mat)       # derive modulus of complex numbers (radius)
    phases.o <- Arg(wt_o_mat)        # extract phases (argument)

    wt_m_mat <- as.matrix(wt_m)
    modulus.m <- Mod(wt_m_mat)
    phases.m <- Arg(wt_m_mat)

    wt_p_mat <- as.matrix(wt_p)
    modulus.p <- Mod(wt_p_mat)
    phases.p <- Arg(wt_p_mat)

    ###================================###===================================###
    ## QM ----
    if((QM=="MBCp") | (QM=="MBCr") | (QM=="MBCn")){
      modulus.tmp <- do.call(QM,list(o.c=modulus.o, m.c=modulus.m, m.p=modulus.p,
                                     ratio.seq=rep(TRUE, ncol(modulus.m)), silent=TRUE))
      modulus.bcc <- modulus.tmp$mhat.c
      modulus.bcf <- modulus.tmp$mhat.p
    } else if(QM=="MRS") {
      modulus.tmp <- do.call(QM,list(o.c=modulus.o, m.c=modulus.m, m.p=modulus.p))
      modulus.bcc <- modulus.tmp$mhat.c
      modulus.bcf <- modulus.tmp$mhat.p
    } else if(QM=="QDM") {
      modulus.tmp <- lapply(1:ncol(modulus.o), function(i)
        MBC::QDM(o.c=modulus.o[,i], m.c=modulus.m[,i], m.p=modulus.p[,i],ratio=TRUE,...))
      modulus.bcc <- sapply(modulus.tmp, function(ls) ls$mhat.c)
      modulus.bcf <- sapply(modulus.tmp, function(ls) ls$mhat.p)
    }

    if(do.plot){
      # calibration
      #summary(modulus.m) %>% print()
      #summary(modulus.bcc) %>% print()
      df.modulus <- rbind(data.frame(mod="obs",no=subset,x=modulus.o),
                          data.frame(mod="cal",no=subset,x=modulus.m),
                          data.frame(mod="bcc",no=subset,x=modulus.bcc)) %>%
        tidyr::gather(lev, amp, 3:11) #%>% mutate(amp=as.numeric(amp), subset=as.numeric(subset))

      df.modulus$lev <- factor(df.modulus$lev, levels = paste0("x.",1:ncol(modulus.o)))
      p.c <- ggplot(df.modulus, aes(x=no, y=amp, color=mod)) +
              geom_line()+
              facet_wrap(.~lev, nrow=3)
      print(p.c)

      # ggplot(df.modulus) +
      #   stat_ecdf(aes(x=amp, color=mod), geom = "step") +
      #   facet_wrap(.~lev, nrow=3)

      # validation
      #summary(modulus.p) %>% print()
      #summary(modulus.bcf) %>% print()

      df.modulus <- rbind(data.frame(mod="val",no=1:nrow(modulus.p),x=modulus.p),
                          data.frame(mod="bcf",no=1:nrow(modulus.p),x=modulus.bcf)) %>%
        tidyr::gather(lev, amp, 3:11) #%>% mutate(amp=as.numeric(amp), subset=as.numeric(subset))

      df.modulus$lev <- factor(df.modulus$lev, levels = paste0("x.",1:ncol(modulus.o)))
      p.f <- ggplot(df.modulus, aes(x=no, y=amp, color=mod)) +
              geom_line()+
              facet_wrap(.~lev, nrow=3)
      print(p.f)

      # ggplot(df.modulus) +
      #   stat_ecdf(aes(x=amp, color=mod), geom = "step") +
      #   facet_wrap(.~lev, nrow=3)
    }

    # post-process
    modulus.bcf[modulus.bcf<0] <- 0
    modulus.bcc[modulus.bcc<0] <- 0

    ###================================###===================================###
    ## bcc----
    mat_new_cal <- matrix(complex(modulus=modulus.bcc,argument=phases.m),ncol=ncol(phases.m))
    rec_cal <- fun_icwt(x.wave=mat_new_cal,dt=dt,dj=dj, flag.wav, scale)
    if(variable=="prep") rec_cal[rec_cal<=theta] <- 0

    if(PR.cal) {
      ### apply wavelet reconstruction to randomized signal----
      #mat_cal_r <- lapply(1:number_sim, function(r) prsim(modulus.bcc, phases.m, noise_mat_cal[[r]]))
      mat_cal_r <- prsim(modulus.bcc, phases.m, noise_mat_cal, method=method, size=block, seed=seed)

      data_sim_cal <- sapply(1:number_sim, function(r) fun_icwt(x.wave=mat_cal_r[[r]], dt=dt, dj=dj, flag.wav, scale))
      if(variable=="prep") data_sim_cal[data_sim_cal<=theta] <- 0
      colnames(data_sim_cal) <- paste0("r",seq(1:number_sim))
    } else{
      data_sim_cal <- data.frame(r=NA)
    }

    data.cal <- data[[l]][subset,]
    data.cal$bcc <-  rec_cal

    ###================================###===================================###
    ## bcf----
    mat_new_val <- matrix(complex(modulus=modulus.bcf,argument=phases.p),ncol=ncol(phases.p))
    rec_val <- fun_icwt(x.wave=mat_new_val,dt=dt,dj=dj, flag.wav, scale)
    if(variable=="prep") rec_val[rec_val<=theta] <- 0

	  ### apply wavelet reconstruction to randomized signal----
    #mat_val_r <- lapply(1:number_sim, function(r) prsim(modulus.bcf, phases.p, noise_mat_val[[r]]))
    mat_val_r <- prsim(modulus.bcf, phases.p, noise_mat_val, method=method, size=block, seed=seed)

	  data_sim_val <- sapply(1:number_sim, function(r) fun_icwt(x.wave=mat_val_r[[r]], dt=dt, dj=dj, flag.wav, scale))
	  if(variable=="prep") data_sim_val[data_sim_val<=theta] <- 0
    colnames(data_sim_val) <- paste0("r",seq(1:number_sim))

    data.val <- data[[l]][-subset,]
    data.val$bcc <-  rec_val

    ### output
    out_list[[l]] <- list(cal=data.frame(data.cal, data_sim_cal),
                          val=data.frame(data.val, data_sim_val))

  }

  return(out_list)

}
