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
	
	#colon - sigmoid
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	colonpats<- read.table("Colonedits.txt", header=TRUE, sep="\t")
	colonsigpats<- subset(colonpats, colonpats[,3] == "Colon - Sigmoid", select=c(1,2,3))
	rownames(colonsigpats)=str_replace_all(colonsigpats[,1], "-", ".")
	colonsiggenes<- merge(colonsigpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	colonsiggenesfinal<- colonsiggenes[,c(4:56207)]
	rownames(colonsiggenesfinal)=colonsiggenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/colonsig")
		colonsigcov<- read.table("colonsigeditcov.txt", header=TRUE, row.names=1, sep="\t")
	colonsigall<- merge(colonsiggenesfinal, colonsigcov, by.x="row.names", by.y="row.names")
	finalcolonsigcov<- data.frame(colonsigall[,c(56206:56238)])	
	rownames(finalcolonsigcov)=colonsigall[,1]
	finalcolonsigL1HS<- data.frame(colonsigall[,c(56205)])
	rownames(finalcolonsigL1HS)= colonsigall[,1]
	colnames(finalcolonsigL1HS)= c("L1HS")
	finalcolonsiggenes<- data.frame(colonsigall[,c(3:56204)])
	rownames(finalcolonsiggenes)=colonsigall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalcolonsiggenes))
		{
		valuepc0 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1])
		valuepc1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,2])
		valuepc2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,3])
		valuein1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,4])
		valuein2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,5])
		valuein3 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,6])
		valuein4 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,7])
		valuein5 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,8])
		valuein6 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,9])
		valuein7 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,10])
		valuein8 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,11])
		valuein9 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,12])
		valuein10 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,13])
		valuein11 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,14])
		valuein12 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,15])
		valuein13 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,16])
		valuein14 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,17])
		valuein15 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,18])
		valuein16 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,19])
		valuein17 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,20])
		valuein18 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,21])
		valuein19 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,22])
		valuein20 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,23])
		valuein21 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,24])
		valuein22 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,25])
		valuein23 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,26])
		valuein24 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,27])
		valuein25 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,28])
		valuein26 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,29])
		valuein27 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,30])
		valuein28 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,31])
		valuein29 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,32])
		valuein30 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,33])
		valueall <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + finalcolonsigcov[,2] + 
					finalcolonsigcov[,3] + finalcolonsigcov[,4] + finalcolonsigcov[,5] + finalcolonsigcov[,6] + finalcolonsigcov[,7] + 
					finalcolonsigcov[,8] + finalcolonsigcov[,9] + finalcolonsigcov[,10] + finalcolonsigcov[,11] + finalcolonsigcov[,12] +
					finalcolonsigcov[,13] + finalcolonsigcov[,14] + finalcolonsigcov[,15] + finalcolonsigcov[,16] + finalcolonsigcov[,17] + finalcolonsigcov[,18]+
					finalcolonsigcov[,19] + finalcolonsigcov[,20] + finalcolonsigcov[,21] + finalcolonsigcov[,22] + finalcolonsigcov[,23] + finalcolonsigcov[,24] +
					finalcolonsigcov[,25] + finalcolonsigcov[,26]+ finalcolonsigcov[,27] + finalcolonsigcov[,28] + finalcolonsigcov[,29] + finalcolonsigcov[,30] + 
					finalcolonsigcov[,31] + finalcolonsigcov[,32] + finalcolonsigcov[,33])
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
	
	setwd("/home/quints1/Correlations/colonsig")
	write.table(finalAICc, "AllAICcresultscolonsig.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalcolonsiggenes))
		{
		if(! all(is.na(finalcolonsiggenes[i]))){
		valuepc0 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1])
		valuepc1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,2])
		valuepc2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,3])
		valuein1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,4])
		valuein2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,5])
		valuein3 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,6])
		valuein4 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,7])
		valuein5 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,8])
		valuein6 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,9])
		valuein7 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,10])
		valuein8 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,11])
		valuein9 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,12])
		valuein10 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,13])
		valuein11 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,14])
		valuein12 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,15])
		valuein13 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,16])
		valuein14 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,17])
		valuein15 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,18])
		valuein16 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,19])
		valuein17 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,20])
		valuein18 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,21])
		valuein19 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,22])
		valuein20 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,23])
		valuein21 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,24])
		valuein22 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,25])
		valuein23 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,26])
		valuein24 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,27])
		valuein25 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,28])
		valuein26 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,29])
		valuein27 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,30])
		valuein28 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,31])
		valuein29 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,32])
		valuein30 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,33])
		valueall <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + finalcolonsigcov[,2] + 
					finalcolonsigcov[,3] + finalcolonsigcov[,4] + finalcolonsigcov[,5] + finalcolonsigcov[,6] + finalcolonsigcov[,7] + 
					finalcolonsigcov[,8] + finalcolonsigcov[,9] + finalcolonsigcov[,10] + finalcolonsigcov[,11] + finalcolonsigcov[,12] +
					finalcolonsigcov[,13] + finalcolonsigcov[,14] + finalcolonsigcov[,15] + finalcolonsigcov[,16] + finalcolonsigcov[,17] + finalcolonsigcov[,18]+
					finalcolonsigcov[,19] + finalcolonsigcov[,20] + finalcolonsigcov[,21] + finalcolonsigcov[,22] + finalcolonsigcov[,23] + finalcolonsigcov[,24] +
					finalcolonsigcov[,25] + finalcolonsigcov[,26]+ finalcolonsigcov[,27] + finalcolonsigcov[,28] + finalcolonsigcov[,29] + finalcolonsigcov[,30] + 
					finalcolonsigcov[,31] + finalcolonsigcov[,32] + finalcolonsigcov[,33])
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
		
		for (i in colnames(finalcolonsiggenes))
		{
			skipped=NULL
			if(!all(is.na(finalcolonsiggenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalcolonsiggenes))
		{
		if(! all(is.na(finalcolonsiggenes[i]))){
		valuepc0 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1])
		valuepc1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,2])
		valuepc2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,3])
		valuein1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,4])
		valuein2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,5])
		valuein3 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,6])
		valuein4 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,7])
		valuein5 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,8])
		valuein6 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,9])
		valuein7 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,10])
		valuein8 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,11])
		valuein9 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,12])
		valuein10 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,13])
		valuein11 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,14])
		valuein12 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,15])
		valuein13 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,16])
		valuein14 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,17])
		valuein15 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,18])
		valuein16 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,19])
		valuein17 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,20])
		valuein18 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,21])
		valuein19 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,22])
		valuein20 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,23])
		valuein21 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,24])
		valuein22 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,25])
		valuein23 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,26])
		valuein24 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,27])
		valuein25 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,28])
		valuein26 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,29])
		valuein27 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,30])
		valuein28 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,31])
		valuein29 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,32])
		valuein30 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,33])
		valueall <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + finalcolonsigcov[,2] + 
					finalcolonsigcov[,3] + finalcolonsigcov[,4] + finalcolonsigcov[,5] + finalcolonsigcov[,6] + finalcolonsigcov[,7] + 
					finalcolonsigcov[,8] + finalcolonsigcov[,9] + finalcolonsigcov[,10] + finalcolonsigcov[,11] + finalcolonsigcov[,12] +
					finalcolonsigcov[,13] + finalcolonsigcov[,14] + finalcolonsigcov[,15] + finalcolonsigcov[,16] + finalcolonsigcov[,17] + finalcolonsigcov[,18]+
					finalcolonsigcov[,19] + finalcolonsigcov[,20] + finalcolonsigcov[,21] + finalcolonsigcov[,22] + finalcolonsigcov[,23] + finalcolonsigcov[,24] +
					finalcolonsigcov[,25] + finalcolonsigcov[,26]+ finalcolonsigcov[,27] + finalcolonsigcov[,28] + finalcolonsigcov[,29] + finalcolonsigcov[,30] + 
					finalcolonsigcov[,31] + finalcolonsigcov[,32] + finalcolonsigcov[,33])
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
		
		
		for (i in colnames(finalcolonsiggenes))
		{
			skipped=NULL
			if(!all(is.na(finalcolonsiggenes[i]))){
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
		setwd("/home/quints1/Correlations/colonsig")
		write.table(finalresults, "P and Q colonsig.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalcolonsiggenes))
		{
		if(! all(is.na(finalcolonsiggenes[i]))){
			pearson<- cor.test(finalcolonsiggenes[,c(i)], finalcolonsigL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalcolonsiggenes[,c(i)], finalcolonsigL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalcolonsiggenes))
		{
			skipped=NULL
			if(!all(is.na(finalcolonsiggenes[i]))){
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
		setwd("/home/quints1/Correlations/colonsig")
		write.table(finalresults, "PRandQColonSig.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		