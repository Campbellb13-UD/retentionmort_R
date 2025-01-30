#' Generate a test dataset that can be input into the `retentionmort` function
#'
#' Using the example code verbatim, this function will produce a series of
#' outputs that can be directly input into the `retentionmort` function.
#'
#' @returns
#'
#' c =     One value for the total number of tagging efforts that is 6 or
#'         greater.
#'
#' err =   A value (between 1 and 26) that represents the weekly mortality rate
#'         and weekly tag loss rate from a preloaded case study listed in the
#'         metadata.
#'
#' n_c1 =  A vector of the number of tagged individuals in one of two
#'         categorical variables.
#'
#' nT =    A vector of the number of tagged individuals for each tagging effort.
#'
#' R =     A vector of the number of recaptured individuals per effort.
#'
#' TaL =   A vector of the cumulative number of tagged individuals following
#'         each effort.
#'
#'
#' @examples
#'
#' #Run this verbatim to produce a single mark-recapture dataset
#' list2env(test_dataset_retentionmort(), envir = .GlobalEnv)
#'
#'@importFrom stats runif
#'
#' @export
test_dataset_retentionmort = function(){

c <- sample(6:100, 1, replace = TRUE) #total number of tagging efforts - minimum 6

nT<- sample(seq(0,500, by = 1), c, replace = TRUE) #number of tagged individuals per tagging effort

TaL<- length(nT) #cumulative number of tagged individuals per effort - will be equal to nT for now
TaL[1] <- nT[1]
for (p in 1:c) {
  TaL[p+1]<- TaL[p]+nT[p+1]
}
TaL<- TaL[-length(TaL)]

n_c1<- nT*runif(1) #proportion of the number of fish in class1 (0-1)


R<- runif(c, min=0*TaL, max=0.2*TaL) #recapture rate for each effort
R<- ceiling(R)
R[1]<-0 #must make the first effort have a recapture rate of 0

err<- sample(1:26, 1, replace = TRUE) #radomly generate a set of error coefficeints from the predetermined list


return(list(c=c, err=err, n_c1=n_c1, nT=nT, R=R, TaL=TaL))

}
