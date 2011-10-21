##**********************************************
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##**********************************************

mytimed=function(t){
t <- round(t)
h <- floor(t/3600)
m <- floor((t%%3600)/60)
s <- (t%%60)%%60
return(list(h=h,m=m,s=s))
}
