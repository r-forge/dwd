###############################################################################
## angle.R
##
## Function to compute angle of two directions 
## The angle is between 0 and 180
## 
## Author: Xiaosun Lu
###############################################################################

angle = function(d1,d2,
		degree=TRUE, ## result in degree
		lessThan90 = FALSE ## result less than 90 degree
		){

	l1 = sqrt(sum(d1^2))
	l2 = sqrt(sum(d2^2))
	
	
	if( abs(l1-1)>0.000000000001  ){
		warning(paste('d1 is not a direction. The length is ',l1,'.',sep=""))
	}
	if( abs(l2-1)>0.000000000001  ){
		warning(paste('d2 is not a direction. The length is ',l2,'.',sep=""))
	}
	
	a = acos(sum(d1*d2)/(l1*l2) )
	
	if(lessThan90){ 
		if(a>pi/2){ a= pi-a }
	}
	
	if(degree){ a = a*180/pi }
	
	return(a)
	
}

