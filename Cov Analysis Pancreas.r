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
	
	#pancreas
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	panpats<- read.table("Pancreasedits.txt", header=TRUE, sep="\t")
	rownames(panpats)=str_replace_all(panpats[,1], "-", ".")
	pangenes<- merge(panpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	pangenesfinal<- pangenes[,c(4:56207)]
	rownames(pangenesfinal)=pangenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/pancreas")
		pancov<- read.table("pancreaseditcov.txt", header=TRUE, row.names=1, sep="\t")
	panall<- merge(pangenesfinal, pancov, by.x="row.names", by.y="row.names")
	finalpancov<- data.frame(panall[,c(56206:56238)])	
	rownames(finalpancov)=panall[,1]
	finalpanL1HS<- data.frame(panall[,c(56205)])
	rownames(finalpanL1HS)= panall[,1]
	colnames(finalpanL1HS)= c("L1HS")
	finalpangenes<- data.frame(panall[,c(3:56204)])
	rownames(finalpangenes)=panall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalpangenes))
		{
		valuepc0 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1])
		valuepc1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,2])
		valuepc2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,3])
		valuein1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,4])
		valuein2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,5])
		valuein3 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,6])
		valuein4 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,7])
		valuein5 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,8])
		valuein6 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,9])
		valuein7 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,10])
		valuein8 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,11])
		valuein9 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,12])
		valuein10 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,13])
		valuein11 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,14])
		valuein12 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,15])
		valuein13 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,16])
		valuein14 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,17])
		valuein15 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,18])
		valuein16 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,19])
		valuein17 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,20])
		valuein18 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,21])
		valuein19 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,22])
		valuein20 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,23])
		valuein21 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,24])
		valuein22 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,25])
		valuein23 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,26])
		valuein24 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,27])
		valuein25 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,28])
		valuein26 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,29])
		valuein27 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,30])
		valuein28 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,31])
		valuein29 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,32])
		valuein30 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,33])
		valueall <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + finalpancov[,2] + 
					finalpancov[,3] + finalpancov[,4] + finalpancov[,5] + finalpancov[,6] + finalpancov[,7] + 
					finalpancov[,8] + finalpancov[,9] + finalpancov[,10] + finalpancov[,11] + finalpancov[,12] +
					finalpancov[,13] + finalpancov[,14] + finalpancov[,15] + finalpancov[,16] + finalpancov[,17] + 
					finalpancov[,18] + finalpancov[,19] + finalpancov[,20] + finalpancov[,21] + finalpancov[,22] + finalpancov[,23] + finalpancov[,24] + 
					finalpancov[,25] + finalpancov[,26] + finalpancov[,27] + finalpancov[,28] + finalpancov[,29] + finalpancov[,30] +
					finalpancov[,31] + finalpancov[,32] + finalpancov[,33])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/pancreas")
	write.table(finalAICc, "AllAICcresultspancreas.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalpangenes))
		{
		if(! all(is.na(finalpangenes[i]))){
		valuepc0 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1])
		valuepc1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,2])
		valuepc2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,3])
		valuein1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,4])
		valuein2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,5])
		valuein3 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,6])
		valuein4 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,7])
		valuein5 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,8])
		valuein6 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,9])
		valuein7 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,10])
		valuein8 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,11])
		valuein9 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,12])
		valuein10 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,13])
		valuein11 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,14])
		valuein12 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,15])
		valuein13 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,16])
		valuein14 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,17])
		valuein15 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,18])
		valuein16 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,19])
		valuein17 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,20])
		valuein18 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,21])
		valuein19 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,22])
		valuein20 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,23])
		valuein21 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,24])
		valuein22 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,25])
		valuein23 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,26])
		valuein24 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,27])
		valuein25 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,28])
		valuein26 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,29])
		valuein27 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,30])
		valuein28 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,31])
		valuein29 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,32])
		valuein30 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,33])
		valueall <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + finalpancov[,2] + 
					finalpancov[,3] + finalpancov[,4] + finalpancov[,5] + finalpancov[,6] + finalpancov[,7] + 
					finalpancov[,8] + finalpancov[,9] + finalpancov[,10] + finalpancov[,11] + finalpancov[,12] +
					finalpancov[,13] + finalpancov[,14] + finalpancov[,15] + finalpancov[,16] + finalpancov[,17] + 
					finalpancov[,18] + finalpancov[,19] + finalpancov[,20] + finalpancov[,21] + finalpancov[,22] + finalpancov[,23] + finalpancov[,24] + 
					finalpancov[,25] + finalpancov[,26] + finalpancov[,27] + finalpancov[,28] + finalpancov[,29] + finalpancov[,30] +
					finalpancov[,31] + finalpancov[,32] + finalpancov[,33])
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
		
		for (i in colnames(finalpangenes))
		{
			skipped=NULL
			if(!all(is.na(finalpangenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalpangenes))
		{
		if(! all(is.na(finalpangenes[i]))){
		valuepc0 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1])
		valuepc1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,2])
		valuepc2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,3])
		valuein1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,4])
		valuein2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,5])
		valuein3 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,6])
		valuein4 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,7])
		valuein5 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,8])
		valuein6 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,9])
		valuein7 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,10])
		valuein8 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,11])
		valuein9 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,12])
		valuein10 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,13])
		valuein11 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,14])
		valuein12 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,15])
		valuein13 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,16])
		valuein14 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,17])
		valuein15 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,18])
		valuein16 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,19])
		valuein17 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,20])
		valuein18 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,21])
		valuein19 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,22])
		valuein20 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,23])
		valuein21 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,24])
		valuein22 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,25])
		valuein23 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,26])
		valuein24 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,27])
		valuein25 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,28])
		valuein26 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,29])
		valuein27 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,30])
		valuein28 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,31])
		valuein29 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,32])
		valuein30 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,33])
		valueall <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + finalpancov[,2] + 
					finalpancov[,3] + finalpancov[,4] + finalpancov[,5] + finalpancov[,6] + finalpancov[,7] + 
					finalpancov[,8] + finalpancov[,9] + finalpancov[,10] + finalpancov[,11] + finalpancov[,12] +
					finalpancov[,13] + finalpancov[,14] + finalpancov[,15] + finalpancov[,16] + finalpancov[,17] + 
					finalpancov[,18] + finalpancov[,19] + finalpancov[,20] + finalpancov[,21] + finalpancov[,22] + finalpancov[,23] + finalpancov[,24] + 
					finalpancov[,25] + finalpancov[,26] + finalpancov[,27] + finalpancov[,28] + finalpancov[,29] + finalpancov[,30] +
					finalpancov[,31] + finalpancov[,32] + finalpancov[,33])
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
		
		
		for (i in colnames(finalpangenes))
		{
			skipped=NULL
			if(!all(is.na(finalpangenes[i]))){
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
		setwd("/home/quints1/Correlations/pancreas")
		write.table(finalresults, "P and Q pancreas.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalpangenes))
		{
		if(! all(is.na(finalpangenes[i]))){
			pearson<- cor.test(finalpangenes[,c(i)], finalpanL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalpangenes[,c(i)], finalpanL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalpangenes))
		{
			skipped=NULL
			if(!all(is.na(finalpangenes[i]))){
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
		setwd("/home/quints1/Correlations/pancreas")
		write.table(finalresults, "PRandQPancreas.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		