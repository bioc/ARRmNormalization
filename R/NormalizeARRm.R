#### ARRm Normalization

data(ProbesType)
assign("ProbesType", ProbesType, envir = .GlobalEnv)

################################################################################
################################################################################

getDesignInfo <- 
    function(sampleNames=NULL, chipVector=NULL, positionVector=NULL){
		
		if (!is.null(sampleNames)){
	   		sampleNames= as.character(as.matrix(sampleNames))
	    }
	    
	    
	    ### Extraction of the position numbers and chip numbers
		if (is.null(chipVector) || is.null(positionVector)){
			if (is.null(sampleNames)){
				stop('Both position vector and chip vector must be provided, 
				    or sampleNames must be provided')
			}
			chips<-as.numeric(factor(substr(sampleNames, 1, 10)))        
			positions<-as.numeric(factor(substr(sampleNames, 12, 17)))   
			designInfo = 
				data.frame(chipInfo=chips, positionInfo=positions, 
				    sampleNames=sampleNames, stringsAsFactors=FALSE)
		
		  
		  
		### If chipVector and positionVector are provided:
		} else {
				  
			if (length(positionVector)!=length(chipVector)){
				 stop('positionVector and 
				     chipVector must have the same length')
			}
			if (sum(positionVector>12)>0){
				 stop('Position number cannot be greater than or equal to 13')
			}
			if (sum(is.na(chipVector))>1 || sum(is.na(positionVector))>1){ 
				stop('One or more position or chip numbers missing')
			}
			chips <- as.numeric(chipVector)
			positions <- as.numeric(positionVector)
			designInfo = 
				data.frame(chipInfo=chips, positionInfo=positions, 
				    stringsAsFactors=FALSE)
		}
		
		
		return(designInfo)	
    }
    
################################################################################
################################################################################








################################################################################
################################################################################

getBackground <- 
	    function(greenControlMatrix, redControlMatrix){
			if (ncol(greenControlMatrix) != ncol(redControlMatrix)){
				stop("Both control matrices must
				    have the same number of columns")
			}
			
			### Green background:
			back.green <- 
			    apply(greenControlMatrix, 2, 
			        function(x) 
			            median(x, na.rm=TRUE)) 
			            
			           
			### Red background:             
			back.red <-
				apply(redControlMatrix, 2, 
				    function(x) 
				        median(x,na.rm=TRUE)) 	
				        			        
				        
			backgroundInfo <- data.frame(green=back.green, red=back.red)	
			return (backgroundInfo)
	    }
################################################################################
################################################################################






### To get the quantiles of the beta-value distributions for each probe type
################################################################################
################################################################################

getQuantiles <- 
    function(betaMatrix,goodProbes=NULL){
    	perc<-seq(from=1,to=100,by=1)/100 
    	    ### Percentiles used for the robust regression
		data(ProbesType,package="ARRmNormalization")
		assign("ProbesType", ProbesType, envir = .GlobalEnv)
		
		if (is.null(goodProbes)){
			probes<-as.character(as.matrix(ProbesType$Probe_Name))} else { 
			probes<-goodProbes
		}
	
		probes.design <- 
		    as.character(as.matrix(ProbesType$Design_Type))[
		        match(probes,as.character(as.matrix(ProbesType$Probe_Name)))]
		        
		        
		## The three following subsets of indices refer to the three subsets 
		## of probes to be normalized according to the good probes 
		## given by the user, and the probe design 
		## (Infinium I Green/Red, Infinium II)
		probeTypeIGreen <- 
		    match(probes[which(probes.design=="I Green")],
		        as.character(as.matrix(ProbesType$Probe_Name)))
		        
		probeTypeIRed <- 
		    match(probes[which(probes.design=="I Red")],
		        as.character(as.matrix(ProbesType$Probe_Name)))
		        
		probeTypeII <- 
		    match(probes[which(probes.design=="II")],
		        as.character(as.matrix(ProbesType$Probe_Name)))
		 
		
		### Remark: columns are samples,
		###    rows are percentiles
		quantilesIGreen <- 
		    apply(betaMatrix, 2, 
		        function(x)
		            quantile(x[probeTypeIGreen], probs=perc, na.rm=TRUE))  
		            
		            
		quantilesIRed <- 
		    apply(betaMatrix, 2, 
		        function(x)
		           quantile(x[probeTypeIRed], probs=perc, na.rm=TRUE))  
		           
		           
		quantilesII <- 
		    apply(betaMatrix, 2, 
		        function(x)
		            quantile(x[probeTypeII], probs=perc, na.rm=TRUE))  
		
		
		
		return(list(green=quantilesIGreen,red=quantilesIRed,II=quantilesII))
}

################################################################################
################################################################################







# For each probe type, extract the estimated coefficients associated with biases
################################################################################
################################################################################

getCoefficients <- 
	function(quantiles, designInfo, backgroundInfo, outliers.perc=0.02){
		positions=designInfo$positionInfo
		chips=designInfo$chipInfo
		back.red=backgroundInfo$red
		back.green=backgroundInfo$green
		dye.bias=log(back.red/back.green)
		ns=nrow(designInfo)
	
	
	
	    #### Type I Green Probes
	    green <- 
	        ARRm.regression(quantiles=quantiles$green, positions=positions,
	            chips=chips, background=log(back.green), dye.bias=NULL)
	            
	    res=green$res
	    
	    #### Indices of the outliers to be removed from the regression:
	    outliers<-order(res)[floor(ns*(1-outliers.perc)):ns]   
	        
		green <- 
		    ARRm.regression(quantiles=quantiles$green[,-outliers], 
		    	positions=positions[-outliers], 
		    		chips=chips[-outliers], 
		    		    background=log(back.green)[-outliers], 
		    		        dye.bias=NULL)
	    
	    
	    
	    
	    #### Type I Red Probes
	    red <- 
	        ARRm.regression(quantiles=quantiles$red, positions=positions,
	            chips=chips, background=log(back.red), dye.bias=NULL)
	            
	    res=red$res
	    outliers<-order(res)[floor(ns*(1-outliers.perc)):ns]  
		red <- 
		    ARRm.regression(quantiles=quantiles$red[,-outliers],
		        positions=positions[-outliers],
		            chips=chips[-outliers],
		                background=log(back.red)[-outliers],
		                    dye.bias=NULL)
	    
	    
	    
	    
	    #### Type II Probes
	    II <- 
	        ARRm.regression(quantiles=quantiles$II, positions=positions,
	            chips=chips, background=log(back.red*back.green), 
	                dye.bias=dye.bias)
	            
	    res=II$res
	    outliers<-order(res)[floor(ns*(1-outliers.perc)):ns]  
	  	II <- 
	  	    ARRm.regression(quantiles=quantiles$II[,-outliers],
	  	        positions=positions[-outliers],
	  	            chips=chips[-outliers],
	  	                background=log(back.red*back.green)[-outliers],
	  	                    dye.bias=dye.bias[-outliers])
	  	    
		return(list(green=green, red=red, II=II))
}

################################################################################
################################################################################







### Regression of the quantiles against background intensity, 
###    dye bias, on-chip position and chip index. 
################################################################################
################################################################################

ARRm.regression <- 
    function(quantiles, positions, chips, background, dye.bias){
			c <- C(factor(chips), sum)
			p <- C(factor(positions), sum)
			ns <- length(positions)
			numberOfChips <- length(table(c))
			perc <- seq(from=1, to=100, by=1)/100 
			    ### Percentiles used for the robust regression
			
			
			chip.variations <- matrix(, length(perc), numberOfChips)  	
			    #### Chip deviation from the mean
			position.variations <- matrix(, length(perc), 12)  
			    #### Position deviation from the chip mean 
			background.vector <- matrix(, length(perc), 1) 
			dyebias.vector<-matrix(, length(perc), 1) 
			res<-rep(0, ns)
			
			
			
			if (is.null(dye.bias)){
				for (percindex in 1:length(perc)){
					fit <- lm(unlist(quantiles[percindex,])~c+p+background)
					res<-res+(fit$residuals)^2
					for (j in 1:(numberOfChips-1)){
						chip.variations[percindex,j] <-
					        fit$coefficients[j+1]
					}
					for (j in 1:11){
						position.variations[percindex,j] <- 
						    fit$coefficients[j+numberOfChips]
					}
					
					
					background.vector[percindex] <- 
					    fit$coefficients[numberOfChips+12]
					chip.variations[percindex,numberOfChips] <- 
					   (-sum(chip.variations[percindex,1:(numberOfChips-1)]))
					position.variations[percindex,12] <- 
					    (-sum(position.variations[percindex,1:11]))
					    
					    
				}
			} else {
				for (percindex in 1:length(perc)){
					
					fit <- 
					    lm(unlist(quantiles[percindex,])~c+p+background+dye.bias)
					res<-res+(fit$residuals)^2
					for (j in 1:(numberOfChips-1)){
						chip.variations[percindex,j] <- 
						    fit$coefficients[j+1]
					}
					for (j in 1:11){
						position.variations[percindex,j] <-
						    fit$coefficients[j+numberOfChips]
					}
					
					
					background.vector[percindex] <- 
					    fit$coefficients[numberOfChips+12]    
					dyebias.vector[percindex] <- 
					    fit$coefficients[numberOfChips+13]
					chip.variations[percindex,numberOfChips] <- 
					    (-sum(chip.variations[percindex,1:(numberOfChips-1)]))
					position.variations[percindex,12] <- 
					    (-sum(position.variations[percindex,1:11])) 
					    
				} 
			}
			
			
			
			
			return(list(res=res, background.vector=background.vector, 
			    dyebias.vector=dyebias.vector, chip.variations=chip.variations,
			        position.variations=position.variations))
	}	
	
################################################################################
################################################################################






### Normalization of Type I probes
################################################################################
################################################################################
normalizeI <- 
    function(betaMatrix, background, probeType, designInfo, 
        outliers.perc, chipCorrection){
        	
    	ns=ncol(betaMatrix)
		c<-C(factor(designInfo$chipInfo),sum)
		p<-C(factor(designInfo$positionInfo),sum)
		numberOfChips<-length(table(designInfo$chipInfo))
		perc <- seq(from=1,to=100,by=1)/100 
		    ### Percentiles used for the robust regression
		
		log.background <- log(background)
		log.background <- log.background-median(log.background) 
		    ### Centralization of the background variable
		
		
		quantiles <- 
		    apply(betaMatrix, 2, 
		        function(x) 
		            quantile(x[probeType], probs=perc, na.rm=TRUE))  
		                
		
		
		firstRegression <- 
		    ARRm.regression(quantiles=quantiles,
		        positions=designInfo$positionInfo, 
		        chips=designInfo$chipInfo, 
		            background=log.background, 
		                dye.bias=NULL)
		        
		        
		        
		res=firstRegression$res
		outliers<-order(res)[floor(ns*(1-outliers.perc)):ns]   
		    #### Indices of the outliers to be removed from the regression
		    
		    
		finalEstimation <- 
		    ARRm.regression(quantiles=quantiles[,-outliers],
		        positions=designInfo$positionInfo[-outliers],
		            chips=designInfo$chipInfo[-outliers],
		                background=log.background[-outliers],
		                    dye.bias=NULL)
	
	
	
	
		################  ----- Normalization step ------ ##################
		for (si in 1:ns)
		{
			smoothingParameter=0.1
			
			
			### Individual sample information 
			sample.position <- as.numeric(p[si])   
			sample.chip <- as.numeric(c[si])       
			back.intensity <- log.background[si] 
			
			  
			
			beta <- betaMatrix[probeType,si]            
			beta.reduced <- beta[which(!is.na(beta))]   
			rankedbeta.reduced <-  
			    rank(beta.reduced, ties.method="first", na.last="keep")
			rankedbeta.reduced <- rankedbeta.reduced/length(beta.reduced)
	
	
	        	
	        #####  Interpolation of the coefficients across percentiles
			chip.function <- 
			    smooth.spline(perc, 
			        finalEstimation$chip.variations[,sample.chip], 
			    	    spar=smoothingParameter)
			    	
			position.function <- 
			    smooth.spline(perc, 
			        finalEstimation$position.variations[,sample.position], 
			            spar=smoothingParameter)
			        
			background.function <-
			    smooth.spline(perc, 
			        finalEstimation$background.vector, 
			            spar=smoothingParameter)
	
	
	
			##### Slopes and variations predicted by regression
			background.slopes <- rep(NA,length(beta))
			chipvar <- rep(NA,length(beta))
			posvar <- rep(NA,length(beta))
			background.slopes[which(!is.na(beta))] <- 
			    predict(background.function, rankedbeta.reduced)$y
			chipvar[which(!is.na(beta))] <-
			     predict(chip.function, rankedbeta.reduced)$y
			posvar[which(!is.na(beta))] <- 
			    predict(position.function, rankedbeta.reduced)$y
	
	
			### Normalization:
			if (chipCorrection){
				betaMatrix[probeType,si] <- 
				    beta-background.slopes*(back.intensity)-chipvar-posvar
			} else {
				betaMatrix[probeType,si] <- 
				    beta-background.slopes*(back.intensity)-posvar
			}
		} 
		
		
		return(betaMatrix)
}
################################################################################
################################################################################






### Normalization of Type II probes
################################################################################
################################################################################

normalizeII <- 
	function(betaMatrix, back.green, back.red, probeType, 
	    designInfo, outliers.perc, chipCorrection){
		ns <- ncol(betaMatrix)
		c <- C(factor(designInfo$chipInfo), sum)
		p <- C(factor(designInfo$positionInfo), sum)
		numberOfChips <- length(table(designInfo$chipInfo))
		perc <- seq(from=1, to=100, by=1)/100 
		    ### Percentiles used for the robust regression
		    
		    
		dye.bias <- log(back.green/back.red)
		dye.bias <- dye.bias-median(dye.bias) 
		    ### Centralization of the dye bias variable
		    
		    
		mean.intensity <- log(back.green*back.red)
		mean.intensity <- mean.intensity-median(mean.intensity) 
		    ### Centralization of the background variable
		    
		
	
		quantiles <- 
		    apply(betaMatrix, 2,
		        function(x) 
		            quantile(x[probeType], probs=perc, na.rm=TRUE))
		            
		            
		firstRegression <- 
		     ARRm.regression(quantiles=quantiles, 
		         positions=designInfo$positionInfo,
		             chips=designInfo$chipInfo,
		                 background=mean.intensity,
		                     dye.bias)
		         
		         
		res=firstRegression$res
		outliers<-order(res)[floor(ns*(1-outliers.perc)):ns]   
		finalEstimation <- 
		    ARRm.regression(quantiles=quantiles[,-outliers],
		         positions=designInfo$positionInfo[-outliers],
		              chips=designInfo$chipInfo[-outliers],
		                  background=mean.intensity[-outliers],
		                      dye.bias=dye.bias[-outliers])
		                      
		                      
	
		#################### Normalization step  ##########################
		for (si in 1:ns)
		{
			
			smoothingParameter=0.1
			
			
			### Individual sample information 
			sample.position <- as.numeric(p[si])
			sample.chip <- as.numeric(c[si])	   
			back.intensity <- mean.intensity[si] 
			sample.dyebias <- dye.bias[si]	   
	
		
		
			beta <- betaMatrix[probeType,si]
			beta.reduced <- beta[which(!is.na(beta))] 
			rankedbeta.reduced <- 
			    rank(beta.reduced, ties.method="first", na.last="keep")
			rankedbeta.reduced <- rankedbeta.reduced/length(beta.reduced)
	
	
		        	
	        	
	        #####  Interpolation of the coefficients across percentiles	
			chip.function <- 
			    smooth.spline(perc, 
			        finalEstimation$chip.variations[,sample.chip], 
			    	    spar=smoothingParameter)
			    	
			position.function <-
			    smooth.spline(perc, 
			        finalEstimation$position.variations[,sample.position],
			             spar=smoothingParameter)
			        
			background.function <-
			    smooth.spline(perc, finalEstimation$background.vector, 
			        spar=smoothingParameter)
			        
			dyebias.function <-
			   smooth.spline(perc, finalEstimation$dyebias.vector, 
			       spar=smoothingParameter)
			
			
			
			## Slopes predicted by regression:
			background.slopes <- rep(NA, length(beta))
			dyebias.slopes <- rep(NA, length(beta))
			chipvar <- rep(NA, length(beta))
			posvar <- rep(NA, length(beta))
			background.slopes[which(!is.na(beta))] <- 
			    predict(background.function, rankedbeta.reduced)$y
			dyebias.slopes[which(!is.na(beta))] <- 
			    predict(dyebias.function, rankedbeta.reduced)$y
			chipvar[which(!is.na(beta))] <- 
			    predict(chip.function, rankedbeta.reduced)$y
			posvar[which(!is.na(beta))] <-
			    predict(position.function, rankedbeta.reduced)$y
		
		
		
			## Normalization:
			if (chipCorrection) {
				betaMatrix[probeType,si] <-
			        beta-background.slopes*(back.intensity)-
			            dyebias.slopes*(sample.dyebias)-chipvar-posvar
			} else {
				betaMatrix[probeType,si] <- 
				    beta-background.slopes*(back.intensity)-
				        dyebias.slopes*(sample.dyebias)-posvar
			}
		}
		
		
		return(betaMatrix)
}

################################################################################
################################################################################








################################################################################
################################################################################

normalizeARRm <- 
    function(betaMatrix, designInfo, backgroundInfo, outliers.perc=0.02,
        goodProbes=NULL, chipCorrection=TRUE){
        	
	data(ProbesType,package="ARRmNormalization")
	assign("ProbesType", ProbesType, envir = .GlobalEnv)
	if (ncol(betaMatrix)<13) {
		stop('The number of samples must 
		    be greater than 13 (more than one chip)')
	}
	
	if (ncol(betaMatrix) != length(backgroundInfo$green)){
		stop('Dimension of the control matrices must agree with 
		    the dimension of the Beta value matrix')
	}
	
	if (ncol(betaMatrix) != length(designInfo$chipInfo)){
		stop('The dimension of the design matrix must agree with 
		    the dimension of the Beta value matrix')
	}
	
	if (is.null(goodProbes)){
		probes<-as.character(as.matrix(ProbesType$Probe_Name))} else { 
		probes<-goodProbes}
		
		
	probes.design <- as.character(as.matrix(ProbesType$Design_Type))[
	    match(probes,as.character(as.matrix(ProbesType$Probe_Name)))]
	
	## The three following subsets of indices refer to the three subsets 
	## of probes to be normalized according to the good probes 
	## given by the user, and the probe design 
	## (Infinium I Green/Red, Infinium II)
	probeTypeIGreen <- 
	    match(probes[which(probes.design=="I Green")],
	        as.character(as.matrix(ProbesType$Probe_Name)))
	        
	probeTypeIRed <- 
	    match(probes[which(probes.design=="I Red")],
	        as.character(as.matrix(ProbesType$Probe_Name)))
	        
	probeTypeII <- 
	    match(probes[which(probes.design=="II")],
	        as.character(as.matrix(ProbesType$Probe_Name)))
	
	
	
	normMatrix=betaMatrix
	normMatrix[-match(probes,as.character(as.matrix(ProbesType$Probe_Name))),]=NA 
	    ### (Probes not normalized are set to NA)
	
	normMatrix <- 
	    normalizeI(betaMatrix=normMatrix, background=backgroundInfo$green,
	        probeType=probeTypeIGreen, designInfo=designInfo,
	            outliers.perc=outliers.perc, chipCorrection=chipCorrection)   
	            
	            
	normMatrix <-
	    normalizeI(betaMatrix=normMatrix, background=backgroundInfo$red,
	        probeType=probeTypeIRed, designInfo=designInfo,
	            outliers.perc=outliers.perc, chipCorrection=chipCorrection)
	            
	                   
	normMatrix <- 
	    normalizeII(betaMatrix=normMatrix, back.red=backgroundInfo$red,
	        back.green=backgroundInfo$green, probeType=probeTypeII,
	            designInfo=designInfo, outliers.perc=outliers.perc,
	                chipCorrection=chipCorrection)      
	         
	            
	             
	return(normMatrix)
}
################################################################################
################################################################################








# Visualization tool
################################################################################
################################################################################

quantilePlots <- 
    function(quantiles, backgroundInfo, designInfo,
        percentilesI=NULL, percentilesII=NULL){
        	
	    pdf("QuantilePlots.pdf")
	        	
		myColors=rep(c("black", "red"), 10)
		cex=1
		pch=rep(c(20,18), 10)
		lwd=1.5
		
		
		if (is.null(percentilesI)){
			percentilesI = seq(1, 100, 5)
		} else {
			percentilesI = percentilesI
		}
		
		
		if (is.null(percentilesII)){
			percentilesII = seq(1, 100, 10)
		} else {
			percentilesII = percentilesII
		}
		
		
		y <- seq(from=-0.1, to=1, by=11/100)
		quantilesList = list(quantiles$green, quantiles$red, quantiles$II)
	
		
		
		## Background effects plots
		backgroundList <- 
		    list(log(backgroundInfo$green),
		        log(backgroundInfo$red),
		            log(backgroundInfo$red*backgroundInfo$green))
		            
		            
		titles <- c("Background effects on different percentiles - Inf I Green",
		    "Background effects on different percentiles - Inf I Red",
		        "Background effects on different percentiles - Infinium II")
	
		for (j in 1:3){
			if (j==1 || j==2){
				percentiles=percentilesI} else {
					percentiles=percentilesII}
					
			currentQuantiles = quantilesList[[j]]
	    	background = backgroundList[[j]]
	    	
			x <- seq(from=min(background), to=max(background)+0.5,
			    by=(max(background+0.5)-min(background))/10)
			    
			plot(x, y, type="n", 
			    xlab="Log Background intensity",
			        ylab="Beta value",
			            main=titles[j])
			            
			for (i in 1:length(percentiles)){
				currentP=percentiles[i]
				
				points(unlist(background),
				    unlist(currentQuantiles[currentP,]),
				        col=myColors[i],
				            pch=pch[i],
				                cex=cex)
				                
				abline(lm(unlist(currentQuantiles[currentP,])~background),
				    col="black", lwd=lwd)
			}	
		}
		
		
		
		
		## Dye bias plot
		dyebias = log(backgroundInfo$green/backgroundInfo$red)
		currentQuantile = quantilesList[[3]]
		percentiles = percentilesII
		
		x <- seq(from=min(dyebias), to=max(dyebias)+0.5,
		    by=(max(dyebias+0.5)-min(dyebias))/10)
		
		plot(x, y, type="n", 
		    xlab="Dye bias",
		       ylab="Beta value",
		          main="Dye bias effects on different percentiles -
		               Infinium II")
		            
		            
		for (i in 1:length(percentiles)){
			currentP=percentiles[i]
			
			
			points(unlist(dyebias), 
			    unlist(currentQuantile[currentP,]),
			        col=myColors[i],
			            pch=pch[i],
			                cex=cex)
			                
			                
			abline(lm(unlist(currentQuantile[currentP,])~dyebias),
			    col="black",lwd=lwd)
		}	
	
	    dev.off()
}

################################################################################
################################################################################





################################################################################
################################################################################

positionPlots <- 
    function(quantiles, designInfo, percentiles=c(25,50,75)){
		quantilesList=list(quantiles$green, quantiles$red, quantiles$II)
	
		#### To look at the position deviations from the mean
		grandMeansList <- list()
		centeredQuantiles <- list()
		
		for (i in 1:3){
			grandMeansList[[i]]=apply(quantiles[[i]], 1, mean)
		}
		
		
		
		for (i in 1:3){
			centeredQuantiles[[i]]=quantiles[[i]]-grandMeansList[[i]]
		}
		
		
		chipCenteredQuantiles=centeredQuantiles
		chipInfo=as.numeric(designInfo$chipInfo)
		chipIndices=unique(as.numeric(designInfo$chipInfo))
		for (i in 1:3){
			currentQuantiles=chipCenteredQuantiles[[i]]
			for (j in 1:length(chipInfo)){
				chipQuantiles <- 
				    currentQuantiles[,which(chipInfo==chipIndices[j])]
				mean<- apply(chipQuantiles, 1, mean)
				chipQuantiles <- chipQuantiles-mean
				currentQuantiles[,which(chipInfo==chipIndices[j])] <- 
				    chipQuantiles
			}
			chipCenteredQuantiles[[i]]=currentQuantiles
		}
	
		### Position lines
		order=order(designInfo$positionInfo)
		for (i in 1:3){
			chipCenteredQuantiles[[i]]=chipCenteredQuantiles[[i]][,order]
		}
		newpositions=designInfo$positionInfo[order]
		lines=c()
		for (i in 1:(length(newpositions)-1)){
		if (newpositions[i+1]!=newpositions[i]){lines=c(lines, i+0.5)}
		}
	
		## Tick marks for position axis
		ticks=c()
		tick1=lines[1]/2
		lastTick=(nrow(designInfo)+lines[11])/2
		ticks=tick1
		for (i in 2:11){ticks=c(ticks,(lines[i]+lines[i-1])/2)}
		ticks=c(ticks,lastTick)
		
		pdf("PositionPlots.pdf")
		for (percIndex in 1:length(percentiles)){
			currentPerc = percentiles[percIndex]
			percStr = 
			    paste("Positions Variations - Percentile", currentPerc," - ")
			titles=c(paste(percStr,"Inf I Green"),
			    paste(percStr,"Inf I Red"),
			        paste(percStr,"Inf II"))
			        
			for (k in 1:3){
				plot(chipCenteredQuantiles[[k]][currentPerc,],xaxt='n',
				    ylab="Beta value deviations",
				        xlab="Positions",
				            main=titles[k],
				                cex=1, pch=20)
				                    
				axis (1, at = ticks, labels=1:12)
				abline(h=0)
				for (j in 1:length(lines)){
					abline(v=lines[j], lty=2)
				}
			}
		}
		dev.off()
}

################################################################################
################################################################################





