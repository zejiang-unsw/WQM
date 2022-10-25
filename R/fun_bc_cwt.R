#' bc_cwt
#'
#' @param data
#' @param subset
#' @param variable
#' @param theta
#' @param QM
#' @param PRand
#' @param marginal
#' @param number_sim
#' @param wavelet
#' @param dt
#' @param dj
#' @param upperPeriod
#' @param seed
#' @param do.plot
#'
#' @return
#' @export
#'
#' @examples
bc_cwt <- function(data, subset, variable, theta=0.1,
                       QM=c("MBC","MRS","QDM"), PRand=TRUE, marginal='empirical',
                       number_sim=5, wavelet="morlet", dt=1, dj=1/8, upperPeriod, seed=NULL,
                       do.plot=FALSE)
  {
  # data=data.h.df1[112]; subset=subset.df;
  # variable="prep";QM=QM; PRand=TRUE; theta = theta;
  # number_sim=num.sim; wavelet="morlet"; dt=dt; dj=dj; upperPeriod=upperPeriod

  ###===============================###===============================###
  ### Generation of white noise for random phases generation
  ### generate random sample of indices for each simulation run
  if(!is.null(seed)) set.seed(seed)

  if(PRand){
    noise_mat <- list()
    for (r in 1:number_sim) {
      #ts_wn <- rnorm(n=length(data[[1]]$obs[-subset]), mean = 0, sd = 1) ### iid time seris
      data.obs <- as.vector(sapply(data, function(ls)ls$obs[subset]))
      ts_wn <- sample(data.obs, size=length(data[[1]]$obs[-subset]), replace=TRUE)

      #wt_noise <- wavCWT(x=ts_wn,wavelet=wavelet,n.scale=n_wave)
      wt_noise <- t(WaveletComp::WaveletTransform(x=ts_wn,dt=dt,dj=dj,upperPeriod = upperPeriod)$Wave)

      noise_mat[[r]] <- as.matrix(wt_noise)
    }
  }


  ###===============================###===============================###
  ###pre-process: bias correction
  ### center_data: subtract mean from values
  for(l in 1:length(data)){
    data[[l]]$norm.m <-  data[[l]]$mod
    data[[l]]$norm.o <- data[[l]]$obs-mean(data[[l]]$obs[subset],na.rm=TRUE)
    data[[l]]$norm.m[subset] <- data[[l]]$mod[subset]-mean(data[[l]]$mod[subset],na.rm=TRUE)
    data[[l]]$norm.m[-subset] <- data[[l]]$mod[-subset]-mean(data[[l]]$mod[-subset],na.rm=TRUE)

    summary(data[[l]]$norm.o)
    summary(data[[l]]$norm.m)
  }

  # data.h.df1[112][[1]]$obs
  # plot.ts(data[[l]]$obs[subset])
  # plot(data[[l]]$norm.o, ylim=c(0,10), type="l",col=1)
  # lines(data.h.df[[id.test]]$Date,data.h.df[[id.test]]$mod, type="l", col=2)

  ### run through all stations
  out_list<-list()
  for(l in 1:length(data)){
    #cat("Station:", l)

    ### list for storing results
    data_sim <- list()
    data[[l]]$index <- as.numeric(format(data[[l]]$Date,format='%j'))

    ###===============================###===============================###
    ### Use continuous wavelet transform (CWT) to wavelet transform the data
    # wt_o <- wavCWT(x=data[[l]]$norm.o,wavelet=wavelet,n.scale=n_wave)
    # wt_m <- wavCWT(x=data[[l]]$norm.m,wavelet=wavelet,n.scale=n_wave)

    summary(data[[l]]$norm.o[subset])
    wt_o <- t(WaveletComp::WaveletTransform(x=data[[l]]$norm.o[subset],dt=dt,dj=dj,upperPeriod = upperPeriod)$Wave)
    wt_m <- t(WaveletComp::WaveletTransform(x=data[[l]]$norm.m[subset],dt=dt,dj=dj,upperPeriod = upperPeriod)$Wave)
    wt_p <- t(WaveletComp::WaveletTransform(x=data[[l]]$norm.m[-subset],dt=dt,dj=dj,upperPeriod = upperPeriod)$Wave)

    #print(is.na(wt_o)||is.na(wt_m)||is.na(wt_p))
    if(is.na(wt_o)||is.na(wt_m)||is.na(wt_p)) next

    if(l==1) print(paste0("No of decomposition level: ", ncol(wt_o)))

    ### return CWT coefficients as a complex matrix with rows and columns representing times and scales, respectively.
    wt_o_mat <- as.matrix(wt_o)
    real.o <- Re(wt_o_mat)
    modulus.o <- Mod(wt_o_mat)       ### derive modulus of complex numbers (radius)
    phases.o <- Arg(wt_o_mat)        ### extract phases (argument)

    wt_m_mat <- as.matrix(wt_m)
    real.m <- Re(wt_m_mat)
    modulus.m <- Mod(wt_m_mat)
    phases.m <- Arg(wt_m_mat)

    wt_p_mat <- as.matrix(wt_p)
    real.p <- Re(wt_p_mat)
    modulus.p <- Mod(wt_p_mat)
    phases.p <- Arg(wt_p_mat)

    # plot.ts(real.o,xlab=NA)
    # plot.ts(real.m,xlab=NA)
    # colMeans(real.o); colMeans(real.m)
    # colMeans(real.p)
    #
    # plot.ts(modulus.o,xlab=NA)
    # plot.ts(modulus.m,xlab=NA)
    # plot.ts(modulus.p,xlab=NA)
    # cat(dim(modulus.m))

    ###===============================###===============================###
    ### bias correction
    #QM <- switch(4, "MBC","MRS","QDM","QDMF")
    if(QM %like% "MBC"){
      modulus.tmp <- do.call(QM,list(o.c=modulus.o, m.c=modulus.m,
                                                    m.p=modulus.p, ratio.seq=rep(TRUE, ncol(modulus.m)),
                                                    silent=TRUE))
      modulus.bcc <- modulus.tmp$mhat.c
      modulus.bcf <- modulus.tmp$mhat.p
    } else if(QM=="MRS") {
      modulus.tmp <- do.call(QM,list(o.c=modulus.o, m.c=modulus.m, m.p=modulus.p))
      modulus.bcc <- modulus.tmp$mhat.c
      modulus.bcf <- modulus.tmp$mhat.p
    } else if(QM=="QDM") {
      #cat("QDM with ratio=T \n")
      modulus.tmp <- lapply(1:ncol(modulus.o), function(i)
        QDM(o.c=modulus.o[,i], m.c=modulus.m[,i], m.p=modulus.p[,i], ratio=TRUE))#, ratio.max = 1, trace=theta))
      modulus.bcc <- sapply(modulus.tmp, function(ls) ls$mhat.c)
      modulus.bcf <- sapply(modulus.tmp, function(ls) ls$mhat.p)
    } else {
      #cat("QDM with ratio=F \n")
      modulus.tmp <- lapply(1:ncol(modulus.o), function(i)
        QDM(o.c=modulus.o[,i], m.c=modulus.m[,i], m.p=modulus.p[,i], ratio=FALSE))
      modulus.bcc <- sapply(modulus.tmp, function(ls) ls$mhat.c)
      modulus.bcf <- sapply(modulus.tmp, function(ls) ls$mhat.p)

    }
    modulus.bcf[modulus.bcf<0] <-0
    modulus.bcc[modulus.bcc<0] <-0

    phases.tmp <- lapply(1:ncol(phases.o), function(i)
      QDM(o.c=phases.o[,i], m.c=phases.m[,i], m.p=phases.p[,i]))
    phases.bcc <- sapply(phases.tmp, function(ls) ls$mhat.c)
    phases.bcf <- sapply(phases.tmp, function(ls) ls$mhat.p)

    if(do.plot){
    df.modulus <- rbind(data.frame(mod="obs",no=subset,x=modulus.o),
                        data.frame(mod="cur",no=subset,x=modulus.m),
                        data.frame(mod="bcc",no=subset,x=modulus.bcc)) %>%
      gather(lev, amp, 3:11) #%>% mutate(amp=as.numeric(amp), subset=as.numeric(subset))

    df.modulus$lev <- factor(df.modulus$lev, levels = paste0("x.",1:ncol(modulus.o)))
    print(ggplot(df.modulus, aes(x=no, y=amp, color=mod)) +
      geom_line()+
      facet_wrap(.~lev, nrow=3))

    # ggplot(df.modulus) +
    #   stat_ecdf(aes(x=amp, color=mod), geom = "step") +
    #   facet_wrap(.~lev, nrow=3)
    }

    if(do.plot){
      df.modulus <- rbind(data.frame(mod="fut",no=1:nrow(modulus.p),x=modulus.p),
                          data.frame(mod="bcf",no=1:nrow(modulus.p),x=modulus.bcf)) %>%
        gather(lev, amp, 3:11) #%>% mutate(amp=as.numeric(amp), subset=as.numeric(subset))

      df.modulus$lev <- factor(df.modulus$lev, levels = paste0("x.",1:ncol(modulus.o)))
      print(ggplot(df.modulus, aes(x=no, y=amp, color=mod)) +
        geom_line()+
        facet_wrap(.~lev, nrow=3) +
        scale_y_continuous())

      # ggplot(df.modulus) +
      #   stat_ecdf(aes(x=amp, color=mod), geom = "step") +
      #   facet_wrap(.~lev, nrow=3)
    }

    # print(plot.ts(cbind(modulus.p[,9],modulus.bcf[,9]),xlab=NA))
    # print(plot.ts(cbind(modulus.o[,9],modulus.m[,9],modulus.bcc[,9]),xlab=NA))
    ###===============================###===============================###
    ### bcc
    mat_new <- matrix(complex(modulus=modulus.bcc,argument=phases.m),ncol=ncol(phases.m))
    rec<- fun_icwt(x=mat_new,dt=dt,dj=dj)
    #ts.plot(cbind(data[[l]]$norm.o[subset],data[[l]]$norm.m[subset],rec), col=1:3)
    if(!missing(variable)) if(variable=="prep") rec[rec<theta] <- 0

    # ### rescale
    # rec_random <- (rec - mean(rec))*sd(data[[l]]$obs[subset])/sd(rec) + mean(data[[l]]$obs[subset])
    # #rec_random <- rec
    # if(!missing(variable)) if(variable=="prep") rec_random[rec_random<theta] <- 0

    #ts.plot(cbind(data[[l]]$obs[subset],data[[l]]$mod[subset]), col=1:2)
    #ts.plot(cbind(data[[l]]$obs[subset],data[[l]]$mod[subset],rec_random), lwd=c(2,2,1), col=1:3)

    data.cal <- data[[l]][subset,]
    #data.cal$bcc <- rec_random
    data.cal$bcc <- rec

    ###===============================###===============================###
    ###bcf
    # mat_new <- matrix(complex(modulus=modulus.bcf,argument=phases.p),ncol=ncol(phases.p))
    # rec<- fun_icwt(x=mat_new,dt=dt,dj=dj)
    # #ts.plot(cbind(data[[l]]$norm.o[-subset],data[[l]]$norm.m[-subset],rec), col=1:3)
    # if(!missing(variable)) if(variable=="prep") rec[rec<theta] <- 0
    #
    # ### rescale
    # rec_random <- (rec - mean(rec))*sd(data[[l]]$mod[-subset])/sd(rec) + mean(data[[l]]$mod[-subset])
    # if(!missing(variable)) if(variable=="prep") rec_random[rec_random<theta] <- 0
    #
    # #ts.plot(cbind(data[[l]]$obs[-subset],data[[l]]$mod[-subset]), col=1:2)
    # #ts.plot(cbind(data[[l]]$obs[-subset],data[[l]]$mod[-subset],rec_random), lwd=c(2,2,1), col=1:3)
    #
    # ###===============================###===============================###
    # ###post-process: bias correction
    # sd.ratio <- sd(data[[l]]$obs[subset])/sd(data[[l]]$mod[subset])
    # mean.d <- mean(data[[l]]$obs[subset]) - mean(data[[l]]$mod[subset])
    # #rec_random <- (rec - mean(rec)) * sd.ratio + mean(rec) + mean.d
    # rec_random <- (rec_random - mean(rec_random)) * sd.ratio + mean(rec_random) + mean.d
    # #rec_random <- rec
    # if(!missing(variable)) if(variable=="prep") rec_random[rec_random<theta] <- 0

    data.val <- data[[l]][-subset,]
    #data.val$bcc <- rec_random

    ###===============================###===============================###
    ###phase randomization
    if(PRand) {
      #phases.rand <- lapply(1:number_sim, function(r)  Arg(noise_mat[[r]]))

      phases.rand <- lapply(1:number_sim, function(r)  {
        tmp <- Arg(noise_mat[[r]])

        ord.bcf <- apply(phases.p, 2, order)

        tmp.rank <- apply(tmp, 2, sort)
        tmp.n <- sapply(1:ncol(tmp), function(ii) {tmp[ord.bcf[,ii],ii] <- tmp.rank[,ii];
        return(tmp[,ii])})

    		w <- 7
    		groups <- seq(1, nrow(phases.p), by=w)
    		shuff <- do.call(c,lapply(groups, function(x) sample(x:(x+w-1),w)))
    		shuff <- shuff[!shuff>nrow(phases.p)]

    		return(list(mod=modulus.bcf[shuff,], arg=tmp.n[shuff,]))

        #return(tmp.n)
      })

      #phases.rand <- c(list(phases.bcf), phases.rand)

      #data_prsim <- lapply(1:length(phases.rand), function(r) prsim(modulus.bcf, phases.rand[[r]],theta = theta))

	  phases.rand <- c(list(list(mod=modulus.bcf, arg=phases.p)), phases.rand)


	  data_prsim <- lapply(1:length(phases.rand), function(r) prsim(phases.rand[[r]]$mod, phases.rand[[r]]$arg,
                                                              theta=theta))

      # What the distribution parameter should come from?
      #data_sim <- lapply(1:number_sim, function(r) shuffing(data.val, data_prsim[[r]]))

      data_sim <- data_prsim
      names(data_sim) <- c("bcc",paste0("r",seq(1:number_sim)))

      ### output
      out_list[[l]] <- list(cal=data.cal, val=data.frame(data.val, data_sim))

    } else {

      ### output
      out_list[[l]] <- list(cal=data.cal, val=data.val)
    }

  }

  return(out_list)

}
