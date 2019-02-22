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
	
	#breast
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	breastpats<- read.table("Breastedits.txt", header=TRUE, sep="\t")
	rownames(breastpats)=str_replace_all(breastpats[,1], "-", ".")
	breastgenes<- merge(breastpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	breastgenesfinal<- breastgenes[,c(4:56207)]
	rownames(breastgenesfinal)=breastgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/breast")
		breastcov<- read.table("breasteditcov.txt", header=TRUE, row.names=1, sep="\t")
	breastall<- merge(breastgenesfinal, breastcov, by.x="row.names", by.y="row.names")
	finalbreastcov<- data.frame(breastall[,c(56206:56238)])	
	rownames(finalbreastcov)=breastall[,1]
	finalbreastL1HS<- data.frame(breastall[,c(56205)])
	rownames(finalbreastL1HS)= breastall[,1]
	colnames(finalbreastL1HS)= c("L1HS")
	finalbreastgenes<- data.frame(breastall[,c(3:56204)])
	rownames(finalbreastgenes)=breastall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalbreastgenes))
		{
		valuepc0 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1])
		valuepc1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,2])
		valuepc2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,3])
		valuein1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,4])
		valuein2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,5])
		valuein3 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,6])
		valuein4 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,7])
		valuein5 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,8])
		valuein6 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,9])
		valuein7 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,10])
		valuein8 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,11])
		valuein9 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,12])
		valuein10 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,13])
		valuein11 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,14])
		valuein12 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,15])
		valuein13 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,16])
		valuein14 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,17])
		valuein15 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,18])
		valuein16 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,19])
		valuein17 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,20])
		valuein18 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,21])
		valuein19 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,22])
		valuein20 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,23])
		valuein21 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,24])
		valuein22 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,25])
		valuein23 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,26])
		valuein24 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,27])
		valuein25 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,28])
		valuein26 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,29])
		valuein27 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,30])
		valuein28 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,31])
		valuein29 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,32])
		valuein30 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,33])
		valueall <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + finalbreastcov[,2] + 
					finalbreastcov[,3] + finalbreastcov[,4] + finalbreastcov[,5] + finalbreastcov[,6] + finalbreastcov[,7] + 
					finalbreastcov[,8] + finalbreastcov[,9] + finalbreastcov[,10] + finalbreastcov[,11] + finalbreastcov[,12] +
					finalbreastcov[,13] + finalbreastcov[,14] + finalbreastcov[,15] + finalbreastcov[,16] + finalbreastcov[,17] + finalbreastcov[,18]+
					finalbreastcov[,19] + finalbreastcov[,20] + finalbreastcov[,21] + finalbreastcov[,22] + finalbreastcov[,23] + finalbreastcov[,24] +
					finalbreastcov[,25] + finalbreastcov[,26]+ finalbreastcov[,27] + finalbreastcov[,28] + finalbreastcov[,29] + finalbreastcov[,30] + 
					finalbreastcov[,31] + finalbreastcov[,32] + finalbreastcov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/breast")
	write.table(finalAICc, "All AICc results breast.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalbreastgenes))
		{
		if(! all(is.na(finalbreastgenes[i]))){
		valuepc0 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1])
		valuepc1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,2])
		valuepc2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,3])
		valuein1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,4])
		valuein2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,5])
		valuein3 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,6])
		valuein4 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,7])
		valuein5 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,8])
		valuein6 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,9])
		valuein7 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,10])
		valuein8 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,11])
		valuein9 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,12])
		valuein10 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,13])
		valuein11 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,14])
		valuein12 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,15])
		valuein13 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,16])
		valuein14 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,17])
		valuein15 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,18])
		valuein16 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,19])
		valuein17 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,20])
		valuein18 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,21])
		valuein19 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,22])
		valuein20 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,23])
		valuein21 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,24])
		valuein22 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,25])
		valuein23 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,26])
		valuein24 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,27])
		valuein25 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,28])
		valuein26 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,29])
		valuein27 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,30])
		valuein28 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,31])
		valuein29 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,32])
		valuein30 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,33])
		valueall <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + finalbreastcov[,2] + 
					finalbreastcov[,3] + finalbreastcov[,4] + finalbreastcov[,5] + finalbreastcov[,6] + finalbreastcov[,7] + 
					finalbreastcov[,8] + finalbreastcov[,9] + finalbreastcov[,10] + finalbreastcov[,11] + finalbreastcov[,12] +
					finalbreastcov[,13] + finalbreastcov[,14] + finalbreastcov[,15] + finalbreastcov[,16] + finalbreastcov[,17] + finalbreastcov[,18]+
					finalbreastcov[,19] + finalbreastcov[,20] + finalbreastcov[,21] + finalbreastcov[,22] + finalbreastcov[,23] + finalbreastcov[,24] +
					finalbreastcov[,25] + finalbreastcov[,26]+ finalbreastcov[,27] + finalbreastcov[,28] + finalbreastcov[,29] + finalbreastcov[,30] + 
					finalbreastcov[,31] + finalbreastcov[,32] + finalbreastcov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
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
		
		for (i in colnames(finalbreastgenes))
		{
			skipped=NULL
			if(!all(is.na(finalbreastgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalbreastgenes))
		{
		if(! all(is.na(finalbreastgenes[i]))){
		valuepc0 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1])
		valuepc1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,2])
		valuepc2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,3])
		valuein1 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,4])
		valuein2 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,5])
		valuein3 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,6])
		valuein4 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,7])
		valuein5 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,8])
		valuein6 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,9])
		valuein7 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,10])
		valuein8 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,11])
		valuein9 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,12])
		valuein10 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,13])
		valuein11 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,14])
		valuein12 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,15])
		valuein13 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,16])
		valuein14 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,17])
		valuein15 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,18])
		valuein16 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,19])
		valuein17 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,20])
		valuein18 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,21])
		valuein19 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,22])
		valuein20 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,23])
		valuein21 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,24])
		valuein22 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,25])
		valuein23 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,26])
		valuein24 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,27])
		valuein25 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,28])
		valuein26 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,29])
		valuein27 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,30])
		valuein28 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,31])
		valuein29 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,32])
		valuein30 <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,33])
		valueall <-lm(finalbreastL1HS[,1] ~ finalbreastgenes[,c(i)] + finalbreastcov[,1] + finalbreastcov[,2] + 
					finalbreastcov[,3] + finalbreastcov[,4] + finalbreastcov[,5] + finalbreastcov[,6] + finalbreastcov[,7] + 
					finalbreastcov[,8] + finalbreastcov[,9] + finalbreastcov[,10] + finalbreastcov[,11] + finalbreastcov[,12] +
					finalbreastcov[,13] + finalbreastcov[,14] + finalbreastcov[,15] + finalbreastcov[,16] + finalbreastcov[,17] + finalbreastcov[,18]+
					finalbreastcov[,19] + finalbreastcov[,20] + finalbreastcov[,21] + finalbreastcov[,22] + finalbreastcov[,23] + finalbreastcov[,24] +
					finalbreastcov[,25] + finalbreastcov[,26]+ finalbreastcov[,27] + finalbreastcov[,28] + finalbreastcov[,29] + finalbreastcov[,30] + 
					finalbreastcov[,31] + finalbreastcov[,32] + finalbreastcov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
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
		
		
		for (i in colnames(finalbreastgenes))
		{
			skipped=NULL
			if(!all(is.na(finalbreastgenes[i]))){
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
		colnames(finalresult)=c("pvalues", "rsqvalues", "qvalues")
		setwd("/home/quints1/Correlations/breast")
		write.table(finalresults, "P and Q Breast.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalbreastgenes))
		{
		if(! all(is.na(finalbreastgenes[i]))){
			pearson<- cor.test(finalbreastgenes[,c(i)], finalbreastL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalbreastgenes[,c(i)], finalbreastL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalbreastgenes))
		{
			skipped=NULL
			if(!all(is.na(finalbreastgenes[i]))){
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
		setwd("/home/quints1/Correlations/breast")
		write.table(finalresults, "PRandQBreast.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
			