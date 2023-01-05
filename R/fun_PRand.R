#' Phase randomization
#'
#' @param modulus
#' @param phases.p
#' @param noise_mat
#' @param variable
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
prsim <- function(modulus, phases.p, noise_mat, dt, dj, method=c("M1","M2")[2]){

    ###===============================###===============================###
    ### use the noise matrix corresponding to this run
    tmp <- Arg(noise_mat)

    if(method=="M1"){
      phases <- tmp

    } else if(method=="M2"){
      ord.bcf <- apply(phases.p, 2, order)
      tmp.rank <- apply(tmp, 2, sort)
      tmp.n <- sapply(1:ncol(tmp), function(ii) {tmp[ord.bcf[,ii],ii] <- tmp.rank[,ii];
      return(tmp[,ii])})

      if(F){
        #cat("M2-1")
        w <- 6
        groups <- seq(1, nrow(phases.p), by=w)
        shuff <- do.call(c,lapply(groups, function(x) sample(x:(x+w-1),w)))
        shuff <- shuff[!shuff>nrow(phases.p)]
        #cat(shuff,'----')
        #shuff <- 1:nrow(phases.p)
        phases <- tmp.n[shuff,]
        modulus <- modulus[shuff,]
      } else {
        #cat("M2-2")
        shuff <- ShuffleBlocks(1:nrow(phases.p), block=12)
        shuff <- shuff[!shuff>nrow(phases.p)]
        #cat(shuff,'----')
        phases <- tmp.n[shuff,]
        modulus <- modulus[shuff,]
      }
    }

    mat_new <- matrix(complex(modulus=modulus,argument=phases),ncol=ncol(phases))

    ### apply wavelet reconstruction to randomized signal
    rec<- fun_icwt(x=mat_new, dt=dt, dj=dj)

    return(rec)
}

guyrot <- function(v, n)
  {
    l <- length(v)
    n <- n %% l
    if(n == 0)
      return(v)
    tmp <- v[(l - n + 1):l]
    v[(n + 1):l] <- v[1:(l - n)]
    v[1:n] <- tmp
    v
  }

ShuffleBlocks <- function(v, block = 6L) {
  #block <- as.integer(block)
  #stopifnot(length(v) %% block == 0L)
  v <- 1: ((length(v) %/% block+1)*block)

  mat <- matrix(v, nrow = block)
  #out <- as.vector(apply(mat, 2, sample)) # random sample within block

  out <- as.vector(apply(mat, 2, function(x) guyrot(x,sample(1:block,1)))) #rotate within block

  #out <- as.vector(mat[, sample(ncol(mat))]) # block shuffle
  #print(out)
  return(out)
}
