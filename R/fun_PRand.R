#' Phase randomization
#'
#' @param modulus
#' @param phases
#' @param dt
#' @param dj
#' @param variable
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
prsim <- function(modulus, phases, dt=1, dj=1/8, variable="prep", theta){

    ###===============================###===============================###
    ### use the R-package wmtsa, which relates to the book by Percival and Walden
    ### on wavelet methods for time series analysis
    ### allows for a flexible range of different wavelet filters: Morlet, Daubechies, Gaussian,...
    ### A) Produce surrogates using phase randomization as in Chavez and Cazelles 2019
    ###===============================###===============================###
    ### A) Produce surrogates using phase randomization as in Chavez and Cazelles 2019
    ### Requirement: choose a complex values filter: Morlet
    ### later on test alternatives: e.g. Gaussian filter
    ### i) use continuous wavelet transform (CWT) to wavelet transform the data
    ### then, follow the randomization procedure proposed by Chavez and Cazelles 2019
    ### ii) generate a Gaussian white noise time series to match the original data length
    ### iii) derive the wavelet transform of this noise to extract the phase
    ### iv) combine this randomised phase and the WT modulus of the original signal to obtain a surrogate time-frequency distribution
    ### v) inverse wavelet transform
    ### vi) rescale the surrogate to the distribution of the original time series by sorting the data
    ### (after a wavelet filtering in the frequency band of interest) according to the ranking of values of the wavelet-based surrogate

    ###===============================###===============================###
    ### use the noise matrix corresponding to this run
    #noise_mat <- noise_mat_r[[r]]

    ### iv) combine this randomised phase and the WT modulus of the original signal to obtain a surrogate time-frequency distribution
    ### create a new matrix
    ### combine modulus of original series to randomised phase: create new matrix of complex values
    #phases_random <- as.matrix(c(phases.m[1], Arg(noise_mat)[-1,]))
    #phases_random <- Arg(noise_mat)
    # cat(phases_random1-phases_random)


    mat_new <- matrix(complex(modulus=modulus,argument=phases),ncol=ncol(phases))
    #cat(Arg(mat_new)-phases_random)
    #cat(Mod(mat_new)-modulus.c)

    ### plug into the original time-frequency object
    ### wmtsa package does not allow for the inverse transform of a CWT object

    #Extract the real part for the reconstruction: (see Torrence and Campo equation 11)
    ### v) inverse wavelet transform
    ### apply wavelet reconstruction to randomized signal
    rec<- fun_icwt(x=mat_new, dt=dt, dj=dj)
    #rec1 <- ifft(mat_new)
    #ts.plot(cbind(data[[l]]$norm.o,rec,rec1), col=1:3, lwd=c(2,2,1))

    if(variable=="prep") rec[rec<theta] <- 0

    #print(summary(rec))

    ### create new data frame
    #data_new <- data.frame("random"=rec)
    #data_new$index <- data[[l]]$index

    ### use transformed data directly
    #data_new$seasonal <- data_new$random
    #data_new$rank <- rank(data_new$seasonal)
    return(rec)
}
