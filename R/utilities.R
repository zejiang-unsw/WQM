#' Inverse of continuous wavelet transform
#'
#' @param x.wave input complex matrix.
#' @param dt sampling resolution in the time domain.
#' @param dj sampling resolution in the frequency domain.
#'
#' @return reconstruction time series
#' @export
#'
#' @references fun_stoch_sim_wave in PRSim, Brunner and Furrer, 2020.
#'
#' @examples
#' set.seed(100)
#'
#' dt<-1
#' dj<-1/8
#' n <- 100
#' x <- rnorm(n)
#' x.wave <- t(WaveletComp::WaveletTransform(x=x)$Wave)
#' rec <- fun_icwt(x.wave)
#'
#' x.wt <- WaveletComp::analyze.wavelet(data.frame(x=x),"x",dt=dt,dj=dj)
#' rec_orig <- WaveletComp::reconstruct(x.wt,only.sig = FALSE, plot.rec = FALSE)$series$x.r
#'
#' ### compare to original series
#' op <- par(mfrow = c(1, 1), mar=c(3,3,1,1), mgp=c(1, 0.5, 0))
#' plot(1:n, x, type="l", lwd=5, xlab=NA, ylab=NA)
#' lines(1:n, rec, col="red",lwd=3)
#' lines(1:n, rec_orig, col="blue", lwd=1)
#' legend("topright",legend=c("Raw","Inverse","Inverse_orig"),
#'        lwd=c(5,3,1),bg="transparent",bty = "n",
#'        col=c("black","red","blue"),horiz=TRUE)
#' par(op)
fun_icwt<-function(x.wave, dt=1, dj=1/8){

  # dt <- 1
  # dj <- 1/8
  n <- length(x.wave)

  # extract real part of wavelet decomposition
  wt.r <- Re(x.wave)/2 # for wavCWT

  ### define number of scales
  J <- length(wt.r[1,]) - 1

  # Reconstruct as in formula (11), refer to Torrence and Compo, 1998
  dial <- 2*2^((0:J)*dj)
  rec <- rep(NA,(length(wt.r[,1])))
  for(l in 1:(length(wt.r[,1]))){
    rec[l] <- dj*sqrt(dt)/(pi^(-1/4)*0.776)*sum(wt.r[l,]/sqrt(dial)[1:length(wt.r[l,])])
  }

  return(rec)

}

#------------------------------------------------------------------------------#
#' Inverse Fourier transform
#'
#' @param x input time series.
#' @param do.plot Logical value of plot.
#'
#' @return reconstruction time series
#' @export
#'
#' @references fun_stoch_sim in PRSim, Brunner and Furrer, 2020.
#'
#' @examples
#' x <- rnorm(100)
#' x.new <- fun_ifft(x, do.plot=TRUE)
fun_ifft<-function(x, do.plot=FALSE){

  ### compute fast Fourier transform
  x.fft <- fft(x)
  ts_inv <- fft(x.fft, inverse = TRUE)/length(x.fft)

  ### derive modulus of complex numbers (radius)
  modulus <- Mod(x.fft)

  ### extract phases (argument)
  phases <- Arg(x.fft)

  n <- length(x.fft)
  ### determine first half (left part of data)
  first_part <- 2:(floor(n/2)+1)
  ### determine second half (right part of data)
  second_part <- (n+1)-(1:floor(n/2))


  ft_new <- rep(NA, length=n)
  ### add mean value
  ft_new[1] <- x.fft[1]

  ### first half
  ft_new[first_part] <- complex(modulus=modulus[first_part], argument=phases[first_part])
  ### second half with conjugate values (opposite imaginary part)
  ft_new[second_part] <- Conj(ft_new[first_part])

  #ft_new <- complex(modulus=modulus,argument=phases)

  ### this procedure reproduces the original time series
  ### apply the same transformation procedure for the newly generated complex numbers
  ft_inv_new <- fft(ft_new, inverse = TRUE)/length(ft_new)
  ts_inv_new <- Re(ft_inv_new)

  if(do.plot){
    ### compare to original series
    op <- par(mfrow = c(1, 1), mar=c(3,3,1,1), mgp=c(1, 0.5, 0))
    plot(1:n, ts_inv,type="l", lwd=5, xlab=NA, ylab=NA)
    lines(1:n, ts_inv_new,col="red",lwd=3)
    lines(1:n, x,col="blue", lwd=1)
    legend("topright",legend=c("Inverse","Inverse_new", "Raw"),
           lwd=c(5,3,1),bg="transparent",bty = "n",
           col=c("black","red","blue"),horiz=TRUE)
    par(op)
  }

  return(ts_inv_new)
}

#' Function: Total number of decomposition levels
#'
#' @param n
#' @param dt
#' @param dj
#'
#' @return
#' @export
#'
#' @examples
fun_cwt_J <- function(n, dt, dj){
    upperPeriod <- floor(n*dt/3)
    lowerPeriod <- 2*dt

    omega0 = 6
    fourier.factor = (2 * pi)/omega0

    min.scale = lowerPeriod/fourier.factor
    max.scale = upperPeriod/fourier.factor

    J <- as.integer(log2(max.scale/min.scale)/dj)

    return(J)
}
