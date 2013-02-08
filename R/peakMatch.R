
#'Function to match peaks across samples
#'
#'Uses the hungarian algorithm to match peaks across samples, one at a time using a template sample.
#'@param peaks is the matrix of peaks in the columns and samples in the rows
#'@target.index is the index of the template sample.
#'@param max.fill is the value to substitute for NAs in the distance matrix. Should be very large, but if too large, will overflow and give an incorrect matching
#'importFrom clue solve_LSAP
#'@export
peakMatch<-function(peaks,target.index,max.fill=1e12){
  if(any(is.na(peaks[target.index,]))){
    stop("template sample must have a full set of peaks detected.")
  }
  result<-matrix(NA,ncol=ncol(peaks),nrow=nrow(peaks))
  sapply(setdiff(1:nrow(peaks),target.index),function(test){
    M<-as.matrix(dist(c(peaks[target.index,],peaks[test,])))
    #browser()
    M[1:ncol(peaks),1:ncol(peaks)]<-1e12
    M[(ncol(peaks)+1):(2*ncol(peaks)),(ncol(peaks)+1):(2*ncol(peaks))]<-1e12
    M[is.na(M)]<-1e12
    sol<-as.vector(solve_LSAP(M))[1:ncol(peaks)]-ncol(peaks)
    result[test,]<<-peaks[test,sol]
  })
  result[target.index,]<-peaks[target.index,]
  result[,order(result[target.index,],decreasing=FALSE)]
}
# 
# #Some test code
# require(clue)
# peaks<-cbind(rnorm(5,12,sd=0.5),rnorm(5,8,sd=0.5),rnorm(5,4,sd=0.5),rnorm(5,1,sd=0.5))
# peaks<-t(apply(peaks,1,sample)) #reorder each column
# peaks[sample(1:length(peaks),5)]<-NA #add some NA
# 
# #The candidate template rows should all give the same output below
# candidates<-which(!apply(peaks,1,function(x)any(is.na(x))))
# sapply(candidates,function(C){
#   list(peakMatch(peaks,C))
# })
