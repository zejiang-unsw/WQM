#' Function: Inverse of continuous wavelet transform
#'
#' @param x input complex matrix
#'
#' @return  reconstruction time series
#' @export
#'
#' @examples
#'
fun_icwt<-function(x.wave,dt=1,dj=1/8,lowerPeriod=2*dt,scales=NULL){

  #x.wave <- as.matrix(wt_o)
  wt.r<-Re(x.wave)

  ### define number of scales
  J<-length(wt.r[1,])
  # Reconstruct as in formula (11):

  if(is.null(scales)) dial<-2*2^(0:J*dj) else dial<-scales
  rec<-rep(NA,(length(wt.r[,1])))
  for(l in 1:(length(wt.r[,1]))){
    rec[l]<-dj*sqrt(dt)/(pi^(-1/4)*0.776)*sum(wt.r[l,]/sqrt(dial)[1:length(wt.r[l,])]) #refer to Torrence and Compo, 1998
  }

  if(is.null(scales)){
    # Define central angular frequency omega0 and fourier factor:
    omega0 = 6
    #fourier.factor   = (4*pi)/(omega0 + sqrt(2+omega0^2))
    fourier.factor = (2*pi)/omega0

    # Compute scales and periods:
    min.scale = lowerPeriod/fourier.factor             # Convert lowerPeriod to minimum scale
    #max.scale = upperPeriod/fourier.factor             # Convert upperPeriod to maximum scale
    #J = as.integer( log2(max.scale/min.scale) / dj)    # Index of maximum scale -1

    scales = min.scale * 2^((0:J)*dj)        # sequence of scales
    scales.length = length(scales)           # J + 1
    periods = fourier.factor*scales          # sequence of periods

  }

  rec.waves = matrix(0, nrow=length(wt.r[1,]), ncol=length(wt.r[,1]))
  #0.2144548: dj*sqrt(dt)/(pi^(-1/4)*0.776) when dj=0.125, dt=1
  for (s.ind in seq_along(wt.r[1,])) {
    rec.waves[s.ind,] = (Re(wt.r[,s.ind])/sqrt(scales[s.ind]))*dj*sqrt(dt)/(pi^(-1/4)*0.776)
  }

  # reconstructed time series
  x.r  = colSums(rec.waves, na.rm=T)

  #ts.plot(cbind(rec,x.r), col=1:2)
  #cat(sum(abs(x.r-rec)))

  return(x.r)
}

ifft <- function (x) fft(x, inverse = TRUE)/length(x)

fun_ifft<-function(x.fft){

  #x.fft <- mat_new
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
  ft_inv_new <- ifft(ft_new)
  ts_invers_new <- Re(ft_inv_new)

  # plot(data$des,type="l", lwd=2)
  # lines(ts_invers_new,col=2)
  # cat(sum(abs(data$des-ts_invers_new)))

  return(ts_invers_new)
}
