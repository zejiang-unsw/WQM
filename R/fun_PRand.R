#' Phase randomization and shuffling
#'
#' @param modulus  Modulus of complex values.
#' @param phases  Argument of complex values.
#' @param noise_mat Complex matrix from random time series.
#' @param method Shuffling method, M1: non-shuffling and M2: shuffling. M2 by default.
#' @param seed Seed for shuffling process.
#'
#' @return A new complex matrix
#' @export
#'
prsim <- function(modulus, phases, noise_mat, method=c("M1","M2")[2], seed=100){

    if(!is.null(seed)) set.seed(seed)

    ###===============================###===============================###
    ### use the noise matrix corresponding to this run
    mat_new <- vector("list",length(noise_mat))
    for(r in 1:length(noise_mat)){

      tmp <- Arg(noise_mat[[r]])

      if(method=="M1"){ # no reshuffle
        phases <- tmp
      } else if(method=="M2"){

        ord.bcf <- apply(phases, 2, order)
        tmp.rank <- apply(tmp, 2, sort)
        tmp.n <- sapply(1:ncol(tmp), function(ii) {tmp[ord.bcf[,ii],ii] <- tmp.rank[,ii];
        return(tmp[,ii])})

        shuff <- ShuffleBlocks(1:nrow(phases), size=6)
        shuff <- shuff[!shuff>nrow(phases)]

        phases <- tmp.n[shuff,]
        modulus <- modulus[shuff,]
        #if(r==1|r==2) cat('r',r,":",shuff,'-------')
      }

      mat_new[[r]] <- matrix(complex(modulus=modulus,argument=phases),ncol=ncol(phases))

    }

    if(!is.null(seed)) rm(.Random.seed, envir=.GlobalEnv)

    return(mat_new)
}

# shifts (or rotates) the elements of the input vector in a cyclic fashion (end periodicity is used).
# reference: wavethresh
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

# shuffle within a block
# reference https://stackoverflow.com/questions/68576845/shuffle-permutate-vector-block-wise
ShuffleBlocks <- function(v, size = 6L) {
  size <- as.integer(size)
  #stopifnot(length(v) %% size == 0L)
  v <- 1: ((length(v) %/% size+1)*size)

  mat <- matrix(v, nrow = size)
  out <- as.vector(apply(mat, 2, sample)) # random sample within block
  #out <- as.vector(apply(mat, 2, function(x) guyrot(x,sample(1:size,1)))) #rotate within block
  #out <- as.vector(mat[, sample(ncol(mat))]) # block shuffle
  #print(out)
  return(out)
}
