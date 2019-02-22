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
	
	#liver
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	livpats<- read.table("Liveredits.txt", header=TRUE, sep="\t")
	rownames(livpats)=str_replace_all(livpats[,1], "-", ".")
	livgenes<- merge(livpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	livgenesfinal<- livgenes[,c(4:56207)]
	rownames(livgenesfinal)=livgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/liver")
		livcov<- read.table("livereditcov.txt", header=TRUE, row.names=1, sep="\t")
	livall<- merge(livgenesfinal, livcov, by.x="row.names", by.y="row.names")
	finallivcov<- data.frame(livall[,c(56206:56223)])	
	rownames(finallivcov)=livall[,1]
	finallivL1HS<- data.frame(livall[,c(56205)])
	rownames(finallivL1HS)= livall[,1]
	colnames(finallivL1HS)= c("L1HS")
	finallivgenes<- data.frame(livall[,c(3:56204)])
	rownames(finallivgenes)=livall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finallivgenes))
		{
		valuepc0 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1])
		valuepc1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,2])
		valuepc2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,3])
		valuein1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,4])
		valuein2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,5])
		valuein3 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,6])
		valuein4 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,7])
		valuein5 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,8])
		valuein6 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,9])
		valuein7 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,10])
		valuein8 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,11])
		valuein9 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,12])
		valuein10 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,13])
		valuein11 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,14])
		valuein12 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,15])
		valuein13 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,16])
		valuein14 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,17])
		valuein15 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,18])
		valueall <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + finallivcov[,2] + 
					finallivcov[,3] + finallivcov[,4] + finallivcov[,5] + finallivcov[,6] + finallivcov[,7] + 
					finallivcov[,8] + finallivcov[,9] + finallivcov[,10] + finallivcov[,11] + finallivcov[,12] +
					finallivcov[,13] + finallivcov[,14] + finallivcov[,15] + finallivcov[,16] + finallivcov[,17] + finallivcov[,18])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/liver")
	write.table(finalAICc, "AllAICcresultsliver.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finallivgenes))
		{
		if(! all(is.na(finallivgenes[i]))){
		valuepc0 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1])
		valuepc1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,2])
		valuepc2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,3])
		valuein1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,4])
		valuein2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,5])
		valuein3 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,6])
		valuein4 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,7])
		valuein5 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,8])
		valuein6 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,9])
		valuein7 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,10])
		valuein8 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,11])
		valuein9 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,12])
		valuein10 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,13])
		valuein11 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,14])
		valuein12 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,15])
		valuein13 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,16])
		valuein14 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,17])
		valuein15 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,18])
		valueall <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + finallivcov[,2] + 
					finallivcov[,3] + finallivcov[,4] + finallivcov[,5] + finallivcov[,6] + finallivcov[,7] + 
					finallivcov[,8] + finallivcov[,9] + finallivcov[,10] + finallivcov[,11] + finallivcov[,12] +
					finallivcov[,13] + finallivcov[,14] + finallivcov[,15] + finallivcov[,16] + finallivcov[,17] + finallivcov[,18])
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
		
		for (i in colnames(finallivgenes))
		{
			skipped=NULL
			if(!all(is.na(finallivgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finallivgenes))
		{
		if(! all(is.na(finallivgenes[i]))){
		valuepc0 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1])
		valuepc1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,2])
		valuepc2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,3])
		valuein1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,4])
		valuein2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,5])
		valuein3 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,6])
		valuein4 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,7])
		valuein5 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,8])
		valuein6 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,9])
		valuein7 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,10])
		valuein8 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,11])
		valuein9 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,12])
		valuein10 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,13])
		valuein11 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,14])
		valuein12 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,15])
		valuein13 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,16])
		valuein14 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,17])
		valuein15 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,18])
		valueall <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + finallivcov[,2] + 
					finallivcov[,3] + finallivcov[,4] + finallivcov[,5] + finallivcov[,6] + finallivcov[,7] + 
					finallivcov[,8] + finallivcov[,9] + finallivcov[,10] + finallivcov[,11] + finallivcov[,12] +
					finallivcov[,13] + finallivcov[,14] + finallivcov[,15] + finallivcov[,16] + finallivcov[,17] + finallivcov[,18])
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
		
		
		for (i in colnames(finallivgenes))
		{
			skipped=NULL
			if(!all(is.na(finallivgenes[i]))){
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
		setwd("/home/quints1/Correlations/liver")
		write.table(finalresults, "P and Q liver.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finallivgenes))
		{
		if(! all(is.na(finallivgenes[i]))){
			pearson<- cor.test(finallivgenes[,c(i)], finallivL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finallivgenes[,c(i)], finallivL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finallivgenes))
		{
			skipped=NULL
			if(!all(is.na(finallivgenes[i]))){
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
		setwd("/home/quints1/Correlations/liver")
		write.table(finalresults, "PRandQLiver.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		