
#' Program Computes neighbors for 4 nearest neighbor structure on a square regular lattice.
#'
#' Locations are assumed to be 1, . . .,k in row 1, k+1, . . .,2k in row two, etc.
#' @param k size of the regular lattice
#' @return Matrix with neighbors information
#' @export

sqrgrid.4nbrs <- function(k){
  us<-1:k
  vs<-1:k
  nbmat<-matrix(0,k^2,4)
  ucnt<-0
  repeat{
    ucnt<-ucnt+1
    tu<-us[ucnt]
    vcnt<-0
    repeat{
      vcnt<-vcnt+1
      tv<-vs[vcnt]
      nbrs<-NULL
      tsite<-(tu-1)*k+tv
      for(ui in (tu-1):(tu+1)){
        if((ui>0) && (ui<=k) && (tv>0) && (tv<=k))
          nbrs<-c(nbrs,(ui-1)*k +tv)
        else nbrs<-c(nbrs,0)}
      for(vi in (tv-1):(tv+1)){
        if((tu>0) && (tu<=k) && (vi>0) && (vi<=k))
          nbrs<-c(nbrs,(tu-1)*k +vi)
        else nbrs<-c(nbrs,0)}
      nbmat[tsite,]<-nbrs[nbrs!=tsite]
      if(vcnt==k) break
    }
    if(ucnt==k) break
  }
  return(nbmat)
}
