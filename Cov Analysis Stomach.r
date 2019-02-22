###STEP 1: Make all the L1HS patients line up with the genes
	setwd("/home/quints1/Deseqstuff")
	totaltissue<- read.table("allL1HSmatrix.txt", header=TRUE, row.names=1, sep="\t")
	matchcut_new<- readRDS("patientgene.rds")
	ttissue<- t(totaltissue)
	library("stringr")
	colnames(ttissue)=str_replace_all(colnames(ttissue), "-", ".")
	t_TEcntmatrix = t(ttissue)
	 t_cntmatrix = t(matchcut_new)
	m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
	new_patients = m[,1]
	new_genenames = colnames(m[,-1])
	new_cntmatrix = as.matrix(t(m[,-1]))
	finalpatlist<- data.frame(new_patients)
	rownames(finalpatlist)=finalpatlist[,1]
	library("stringi")
	
###STEP 2: Start on doing the models to find best result
	#read in the L1HS vst and cntvst
	setwd("/home/quints1/Deseqstuff")
	L1HSraw<- read.table("L1HS.VST.cnts.txt", header=TRUE, row.names=1)
	genesraw<- read.table("VSTcnt.txt", header=TRUE, row.names=1)
	tgenesraw<- t(genesraw)
	library("stringi")
	library("stringr")
	
	#stomach
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	stopats<- read.table("Stomachedits.txt", header=TRUE, sep="\t")
	rownames(stopats)=str_replace_all(stopats[,1], "-", ".")
	stogenes<- merge(stopats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	stogenesfinal<- stogenes[,c(4:56207)]
	rownames(stogenesfinal)=stogenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/stomach")
		stocov<- read.table("stomacheditcov.txt", header=TRUE, row.names=1, sep="\t")
	stoall<- merge(stogenesfinal, stocov, by.x="row.names", by.y="row.names")
	finalstocov<- data.frame(stoall[,c(56206:56238)])	
	rownames(finalstocov)=stoall[,1]
	finalstoL1HS<- data.frame(stoall[,c(56205)])
	rownames(finalstoL1HS)= stoall[,1]
	colnames(finalstoL1HS)= c("L1HS")
	finalstogenes<- data.frame(stoall[,c(3:56204)])
	rownames(finalstogenes)=stoall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalstogenes))
		{
		valuepc0 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1])
		valuepc1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,2])
		valuepc2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,3])
		valuein1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,4])
		valuein2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,5])
		valuein3 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,6])
		valuein4 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,7])
		valuein5 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,8])
		valuein6 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,9])
		valuein7 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,10])
		valuein8 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,11])
		valuein9 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,12])
		valuein10 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,13])
		valuein11 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,14])
		valuein12 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,15])
		valuein13 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,16])
		valuein14 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,17])
		valuein15 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,18])
		valuein16 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,19])
		valuein17 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,20])
		valuein18 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,21])
		valuein19 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,22])
		valuein20 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,23])
		valuein21 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,24])
		valuein22 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,25])
		valuein23 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,26])
		valuein24 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,27])
		valuein25 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,28])
		valuein26 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,29])
		valuein27 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,30])
		valuein28 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,31])
		valuein29 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,32])
		valuein30 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,33])
		valueall <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + finalstocov[,2] + 
					finalstocov[,3] + finalstocov[,4] + finalstocov[,5] + finalstocov[,6] + finalstocov[,7] + 
					finalstocov[,8] + finalstocov[,9] + finalstocov[,10] + finalstocov[,11] + finalstocov[,12] +
					finalstocov[,13] + finalstocov[,14] + finalstocov[,15] + finalstocov[,16] + finalstocov[,17] + 
					finalstocov[,18] + finalstocov[,19] + finalstocov[,20] + finalstocov[,21] + finalstocov[,22] + finalstocov[,23] + finalstocov[,24] + 
					finalstocov[,25] + finalstocov[,26] + finalstocov[,27] + finalstocov[,28] + finalstocov[,29] + finalstocov[,30] +
					finalstocov[,31] + finalstocov[,32] + finalstocov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/stomach")
	write.table(finalAICc, "AllAICcresultsstomach.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalstogenes))
		{
		if(! all(is.na(finalstogenes[i]))){
		valuepc0 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1])
		valuepc1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,2])
		valuepc2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,3])
		valuein1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,4])
		valuein2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,5])
		valuein3 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,6])
		valuein4 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,7])
		valuein5 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,8])
		valuein6 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,9])
		valuein7 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,10])
		valuein8 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,11])
		valuein9 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,12])
		valuein10 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,13])
		valuein11 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,14])
		valuein12 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,15])
		valuein13 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,16])
		valuein14 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,17])
		valuein15 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,18])
		valuein16 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,19])
		valuein17 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,20])
		valuein18 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,21])
		valuein19 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,22])
		valuein20 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,23])
		valuein21 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,24])
		valuein22 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,25])
		valuein23 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,26])
		valuein24 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,27])
		valuein25 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,28])
		valuein26 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,29])
		valuein27 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,30])
		valuein28 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,31])
		valuein29 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,32])
		valuein30 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,33])
		valueall <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + finalstocov[,2] + 
					finalstocov[,3] + finalstocov[,4] + finalstocov[,5] + finalstocov[,6] + finalstocov[,7] + 
					finalstocov[,8] + finalstocov[,9] + finalstocov[,10] + finalstocov[,11] + finalstocov[,12] +
					finalstocov[,13] + finalstocov[,14] + finalstocov[,15] + finalstocov[,16] + finalstocov[,17] + 
					finalstocov[,18] + finalstocov[,19] + finalstocov[,20] + finalstocov[,21] + finalstocov[,22] + finalstocov[,23] + finalstocov[,24] + 
					finalstocov[,25] + finalstocov[,26] + finalstocov[,27] + finalstocov[,28] + finalstocov[,29] + finalstocov[,30] +
					finalstocov[,31] + finalstocov[,32] + finalstocov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]}
		else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, pval) 
			}
		pvalues<- data.frame(pvalues)
		tpvalues<- t(pvalues)
		skip=NULL
		
		for (i in colnames(finalstogenes))
		{
			skipped=NULL
			if(!all(is.na(finalstogenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalstogenes))
		{
		if(! all(is.na(finalstogenes[i]))){
		valuepc0 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1])
		valuepc1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,2])
		valuepc2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,3])
		valuein1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,4])
		valuein2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,5])
		valuein3 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,6])
		valuein4 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,7])
		valuein5 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,8])
		valuein6 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,9])
		valuein7 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,10])
		valuein8 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,11])
		valuein9 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,12])
		valuein10 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,13])
		valuein11 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,14])
		valuein12 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,15])
		valuein13 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,16])
		valuein14 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,17])
		valuein15 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,18])
		valuein16 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,19])
		valuein17 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,20])
		valuein18 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,21])
		valuein19 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,22])
		valuein20 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,23])
		valuein21 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,24])
		valuein22 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,25])
		valuein23 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,26])
		valuein24 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,27])
		valuein25 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,28])
		valuein26 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,29])
		valuein27 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,30])
		valuein28 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,31])
		valuein29 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,32])
		valuein30 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,33])
		valueall <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + finalstocov[,2] + 
					finalstocov[,3] + finalstocov[,4] + finalstocov[,5] + finalstocov[,6] + finalstocov[,7] + 
					finalstocov[,8] + finalstocov[,9] + finalstocov[,10] + finalstocov[,11] + finalstocov[,12] +
					finalstocov[,13] + finalstocov[,14] + finalstocov[,15] + finalstocov[,16] + finalstocov[,17] + 
					finalstocov[,18] + finalstocov[,19] + finalstocov[,20] + finalstocov[,21] + finalstocov[,22] + finalstocov[,23] + finalstocov[,24] + 
					finalstocov[,25] + finalstocov[,26] + finalstocov[,27] + finalstocov[,28] + finalstocov[,29] + finalstocov[,30] +
					finalstocov[,31] + finalstocov[,32] + finalstocov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		rval<- finalresult$adj.r.squared}
		else{
				print(paste("Skipping", i))
			}

			rvalues <- cbind(rvalues, rval) 
			}
		rvalues<- data.frame(rvalues)
		trvalues<- t(rvalues)
		
		skip=NULL
		
		
		for (i in colnames(finalstogenes))
		{
			skipped=NULL
			if(!all(is.na(finalstogenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}	
			skip<- t(skip)
			rownames(trvalues)= skip[,1]
			
		#qvalues
		testforq<- tpvalues
		testforq[!is.finite(testforq)] <- NA
		testforq<- data.frame(na.omit(testforq))
		qvalues<- data.frame(qvalue(testforq[,1], pi0=1)$qvalues)
		rownames(qvalues)=rownames(tpvalues)
		results<- merge(tpvalues, trvalues, by.x="row.names", by.y="row.names")	
		rownames(results)=results[,1]
		results<- results[,2:3]
		colnames(results)=c("pvalues", "rsqvalues")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("pvalues", "rsqvalues", "qvalues")
		setwd("/home/quints1/Correlations/stomach")
		write.table(finalresults, "P and Q stomach.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalstogenes))
		{
		if(! all(is.na(finalstogenes[i]))){
			pearson<- cor.test(finalstogenes[,c(i)], finalstoL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalstogenes[,c(i)], finalstoL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalstogenes))
		{
			skipped=NULL
			if(!all(is.na(finalstogenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}	
			skip<- t(skip)
			rownames(trvalueresults)= skip[,1]
			
		results<- merge(tpvalues, trvalueresults, by.x="row.names", by.y="row.names")	
		rownames(results)=results[,1]
		results<- results[,2:4]
		colnames(results)=c("pvalues", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/stomach")
		write.table(finalresults, "PRandQStomach.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		