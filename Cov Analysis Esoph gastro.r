###STEP 1: Make all the L1HS patients line up with the genes
###27373.Colontransj  
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
	
	#esophagus - gastroesophageal junction
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Gastroesophageal Junction", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophgenes<- merge(esophpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	esophgenesfinal<- esophgenes[,c(4:56207)]
	rownames(esophgenesfinal)=esophgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/esophagusgastro")
		esophcov<- read.table("esophagusgastroeditcov.txt", header=TRUE, row.names=1, sep="\t")
	esophall<- merge(esophgenesfinal, esophcov, by.x="row.names", by.y="row.names")
	finalesophcov<- data.frame(esophall[,c(56206:56223)])	
	rownames(finalesophcov)=esophall[,1]
	finalesophL1HS<- data.frame(esophall[,c(56205)])
	rownames(finalesophL1HS)= esophall[,1]
	colnames(finalesophL1HS)= c("L1HS")
	finalesophgenes<- data.frame(esophall[,c(3:56204)])
	rownames(finalesophgenes)=esophall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalesophgenes))
		{
		valuepc0 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1])
		valuepc1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,2])
		valuepc2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,3])
		valuein1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,4])
		valuein2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,5])
		valuein3 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,6])
		valuein4 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,7])
		valuein5 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,8])
		valuein6 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,9])
		valuein7 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,10])
		valuein8 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,11])
		valuein9 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,12])
		valuein10 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,13])
		valuein11 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,14])
		valuein12 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,15])
		valuein13 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,16])
		valuein14 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,17])
		valuein15 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,18])
		valueall <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1] + finalesophcov[,2] + 
					finalesophcov[,3] + finalesophcov[,4] + finalesophcov[,5] + finalesophcov[,6] + finalesophcov[,7] + 
					finalesophcov[,8] + finalesophcov[,9] + finalesophcov[,10] + finalesophcov[,11] + finalesophcov[,12] +
					finalesophcov[,13] + finalesophcov[,14] + finalesophcov[,15] + finalesophcov[,16] + finalesophcov[,17] + finalesophcov[,18])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/esophagusgastro")
	write.table(finalAICc, "AllAICcresultsesophgastro.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalesophgenes))
		{
		if(! all(is.na(finalesophgenes[i]))){
		valuepc0 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1])
		valuepc1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,2])
		valuepc2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,3])
		valuein1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,4])
		valuein2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,5])
		valuein3 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,6])
		valuein4 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,7])
		valuein5 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,8])
		valuein6 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,9])
		valuein7 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,10])
		valuein8 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,11])
		valuein9 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,12])
		valuein10 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,13])
		valuein11 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,14])
		valuein12 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,15])
		valuein13 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,16])
		valuein14 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,17])
		valuein15 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,18])
		valueall <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1] + finalesophcov[,2] + 
					finalesophcov[,3] + finalesophcov[,4] + finalesophcov[,5] + finalesophcov[,6] + finalesophcov[,7] + 
					finalesophcov[,8] + finalesophcov[,9] + finalesophcov[,10] + finalesophcov[,11] + finalesophcov[,12] +
					finalesophcov[,13] + finalesophcov[,14] + finalesophcov[,15] + finalesophcov[,16] + finalesophcov[,17] + finalesophcov[,18])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
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
		
		for (i in colnames(finalesophgenes))
		{
			skipped=NULL
			if(!all(is.na(finalesophgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalesophgenes))
		{
		if(! all(is.na(finalesophgenes[i]))){
		valuepc0 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1])
		valuepc1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,2])
		valuepc2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,3])
		valuein1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,4])
		valuein2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,5])
		valuein3 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,6])
		valuein4 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,7])
		valuein5 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,8])
		valuein6 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,9])
		valuein7 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,10])
		valuein8 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,11])
		valuein9 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,12])
		valuein10 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,13])
		valuein11 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,14])
		valuein12 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,15])
		valuein13 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,16])
		valuein14 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,17])
		valuein15 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,18])
		valueall <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1] + finalesophcov[,2] + 
					finalesophcov[,3] + finalesophcov[,4] + finalesophcov[,5] + finalesophcov[,6] + finalesophcov[,7] + 
					finalesophcov[,8] + finalesophcov[,9] + finalesophcov[,10] + finalesophcov[,11] + finalesophcov[,12] +
					finalesophcov[,13] + finalesophcov[,14] + finalesophcov[,15] + finalesophcov[,16] + finalesophcov[,17] + finalesophcov[,18])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
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
		
		
		for (i in colnames(finalesophgenes))
		{
			skipped=NULL
			if(!all(is.na(finalesophgenes[i]))){
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
		setwd("/home/quints1/Correlations/esophagusgastro")
		write.table(finalresults, "P and Q esophgastro.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalesophgenes))
		{
		if(! all(is.na(finalesophgenes[i]))){
			pearson<- cor.test(finalesophgenes[,c(i)], finalesophL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalesophgenes[,c(i)], finalesophL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalesophgenes))
		{
			skipped=NULL
			if(!all(is.na(finalesophgenes[i]))){
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
		setwd("/home/quints1/Correlations/esophagusgastro")
		write.table(finalresults, "PRandQEsophGastro.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		