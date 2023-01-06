#' Phase randomization
#'
#' @param modulus
#' @param phases
#' @param noise_mat
#' @param variable
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
prsim <- function(modulus, phases, noise_mat, method=c("M1","M2")[2]){

    ###===============================###===============================###
    ### use the noise matrix corresponding to this run
    tmp <- Arg(noise_mat)

    if(method=="M1"){ # no reshuffle
      phases <- tmp
    } else if(method=="M2"){
      ord.bcf <- apply(phases, 2, order)
      tmp.rank <- apply(tmp, 2, sort)
      tmp.n <- sapply(1:ncol(tmp), function(ii) {tmp[ord.bcf[,ii],ii] <- tmp.rank[,ii];
      return(tmp[,ii])})

      shuff <- ShuffleBlocks(1:nrow(phases), block=12)
      shuff <- shuff[!shuff>nrow(phases)]

      phases <- tmp.n[shuff,]
      modulus <- modulus[shuff,]
    }

    mat_new <- matrix(complex(modulus=modulus,argument=phases),ncol=ncol(phases))

    return(mat_new)
}

guyrot <- function(v, n){
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
