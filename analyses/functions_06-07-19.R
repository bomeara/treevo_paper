

extractIntrinsic_from_prcOut<-function(prcOut){
	res <- list(
		intrinsicFn = prcOut$intrinsicFn, 
		intrinsicValues = prcOut$parMeansList$intrinsic, 
		startingValues = prcOut$parMeansList$starting
		)
	return(res)
	}
	
extractExtrinsic_from_prcOut<-function(prcOut){
	res <- list(
		extrinsicFn = prcOut$extrinsicFn, 
		extrinsicValues = prcOut$parMeansList$extrinsic
		)
	return(res)
	}
	
