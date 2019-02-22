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
	
	#Prostate
	#GTEX.131XF.5013.SM.7EWFI
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	propats<- read.table("Prostateedits.txt", header=TRUE, sep="\t")
	rownames(propats)=str_replace_all(propats[,1], "-", ".")
	progenes<- merge(propats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	progenesfinal<- progenes[,c(4:56207)]
	rownames(progenesfinal)=progenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/prostate")
		procov<- read.table("prostateeditcov.txt", header=TRUE, row.names=1, sep="\t")
	proall<- merge(progenesfinal, procov, by.x="row.names", by.y="row.names")
	finalprocov<- data.frame(proall[,c(56206:56223)])	
	rownames(finalprocov)=proall[,1]
	finalproL1HS<- data.frame(proall[,c(56205)])
	rownames(finalproL1HS)= proall[,1]
	colnames(finalproL1HS)= c("L1HS")
	finalprogenes<- data.frame(proall[,c(3:56204)])
	rownames(finalprogenes)=proall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalprogenes))
		{
		valuepc0 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,1])
		valuepc1 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,2])
		valuepc2 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,3])
		valuein1 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,4])
		valuein2 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,5])
		valuein3 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,6])
		valuein4 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,7])
		valuein5 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,8])
		valuein6 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,9])
		valuein7 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,10])
		valuein8 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,11])
		valuein9 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,12])
		valuein10 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,13])
		valuein11 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,14])
		valuein12 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,15])
		valuein13 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,16])
		valuein14 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,17])
		valuein15 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,18])
		valueall <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,1] + finalprocov[,2] + 
					finalprocov[,3] + finalprocov[,4] + finalprocov[,5] + finalprocov[,6] + finalprocov[,7] + 
					finalprocov[,8] + finalprocov[,9] + finalprocov[,10] + finalprocov[,11] + finalprocov[,12] +
					finalprocov[,13] + finalprocov[,14] + finalprocov[,15] + finalprocov[,16] + finalprocov[,17] + finalprocov[,18])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/prostate")
	write.table(finalAICc, "AllAICcresultsprostate.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalprogenes))
		{
		if(! all(is.na(finalprogenes[i]))){
		valuepc0 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,1])
		valuepc1 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,2])
		valuepc2 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,3])
		valuein1 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,4])
		valuein2 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,5])
		valuein3 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,6])
		valuein4 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,7])
		valuein5 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,8])
		valuein6 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,9])
		valuein7 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,10])
		valuein8 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,11])
		valuein9 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,12])
		valuein10 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,13])
		valuein11 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,14])
		valuein12 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,15])
		valuein13 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,16])
		valuein14 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,17])
		valuein15 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,18])
		valueall <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,1] + finalprocov[,2] + 
					finalprocov[,3] + finalprocov[,4] + finalprocov[,5] + finalprocov[,6] + finalprocov[,7] + 
					finalprocov[,8] + finalprocov[,9] + finalprocov[,10] + finalprocov[,11] + finalprocov[,12] +
					finalprocov[,13] + finalprocov[,14] + finalprocov[,15] + finalprocov[,16] + finalprocov[,17] + finalprocov[,18])
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
		
		for (i in colnames(finalprogenes))
		{
			skipped=NULL
			if(!all(is.na(finalprogenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalprogenes))
		{
		if(! all(is.na(finalprogenes[i]))){
		valuepc0 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,1])
		valuepc1 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,2])
		valuepc2 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,3])
		valuein1 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,4])
		valuein2 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,5])
		valuein3 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,6])
		valuein4 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,7])
		valuein5 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,8])
		valuein6 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,9])
		valuein7 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,10])
		valuein8 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,11])
		valuein9 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,12])
		valuein10 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,13])
		valuein11 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,14])
		valuein12 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,15])
		valuein13 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,16])
		valuein14 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,17])
		valuein15 <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,18])
		valueall <-lm(finalproL1HS[,1] ~ finalprogenes[,c(i)] + finalprocov[,1] + finalprocov[,2] + 
					finalprocov[,3] + finalprocov[,4] + finalprocov[,5] + finalprocov[,6] + finalprocov[,7] + 
					finalprocov[,8] + finalprocov[,9] + finalprocov[,10] + finalprocov[,11] + finalprocov[,12] +
					finalprocov[,13] + finalprocov[,14] + finalprocov[,15] + finalprocov[,16] + finalprocov[,17] + finalprocov[,18])
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
		
		
		for (i in colnames(finalprogenes))
		{
			skipped=NULL
			if(!all(is.na(finalprogenes[i]))){
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
		setwd("/home/quints1/Correlations/prostate")
		write.table(finalresults, "P and Q prostate.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalprogenes))
		{
		if(! all(is.na(finalprogenes[i]))){
			pearson<- cor.test(finalprogenes[,c(i)], finalproL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalprogenes[,c(i)], finalproL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalprogenes))
		{
			skipped=NULL
			if(!all(is.na(finalprogenes[i]))){
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
		setwd("/home/quints1/Correlations/prostate")
		write.table(finalresults, "PRandQProstate.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		