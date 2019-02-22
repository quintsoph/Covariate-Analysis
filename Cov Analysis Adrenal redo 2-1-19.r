setwd("/home/quints1/Deseqstuff")
genesraw<- read.table("normalizedcnt.txt", header=TRUE, row.names=1)
newmatrix=NULL
for (i in colnames(genesraw))
{
	summation<-sum(genesraw[,i])
	newmatrix<- cbind(newmatrix, summation)
	
	}
colnames(newmatrix)=colnames(genesraw)
log2sub<- log2(newmatrix)
tlog2sub<- t(log2sub)
library("stringr")
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
	setwd("/home/quints1/Deseqstuff")
	L1HSraw<- read.table("L1HS.VST.cnts.txt", header=TRUE, row.names=1)
	genesraw<- read.table("VSTcnt.txt", header=TRUE, row.names=1)
	tgenesraw<- t(genesraw)
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	adrenalpats<- read.table("Adrenaledits.txt", header=TRUE, sep="\t")
	rownames(adrenalpats)=str_replace_all(adrenalpats[,1], "-", ".")
	adrenalgenes<- merge(adrenalpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	adrenalgenesfinal<- adrenalgenes[,c(4:56207)]
	rownames(adrenalgenesfinal)=adrenalgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/adrenal")
		adrenalcov<- read.table("adrenaleditcov.txt", header=TRUE, row.names=1, sep="\t")
	adrenalall<- merge(adrenalgenesfinal, adrenalcov, by.x="row.names", by.y="row.names")
	rownames(adrenalall)=adrenalall[,1]
	adrenalall<- adrenalall[,-c(1,2)]
	adrenalallcombo<-  merge(adrenalall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(adrenalallcombo)=adrenalallcombo[,1]
	finaladrenalcov<- data.frame(adrenalallcombo[,c(56205:56222)])	
	rownames(finaladrenalcov)=adrenalallcombo[,1]
	finaladrenalL1HS<- data.frame(adrenalallcombo[,c(56204)])
	rownames(finaladrenalL1HS)= adrenalallcombo[,1]
	colnames(finaladrenalL1HS)= c("L1HS")
	finaladrenalgenes<- data.frame(adrenalallcombo[,c(2:56203)])
	log2sum<- data.frame(adrenalallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(adrenalallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	
	finalAICc=NULL
		for (i in colnames(finaladrenalgenes))
		{
		valuepc0 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + log2sum[,1])
		valuepc1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,2] + log2sum[,1])
		valuepc2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,3] + log2sum[,1])
		valuein1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,4] + log2sum[,1])
		valuein2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,5] + log2sum[,1])
		valuein3 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,6] + log2sum[,1])
		valuein4 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,7] + log2sum[,1])
		valuein5 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,8] + log2sum[,1])
		valuein6 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,9] + log2sum[,1])
		valuein7 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,10] + log2sum[,1])
		valuein8 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,11] + log2sum[,1])
		valuein9 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,12] + log2sum[,1])
		valuein10 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,13] + log2sum[,1])
		valuein11 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,14] + log2sum[,1])
		valuein12 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,15] + log2sum[,1])
		valuein13 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,16] + log2sum[,1])
		valuein14 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,17] + log2sum[,1])
		valuein15 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,18] + log2sum[,1])
		valueall <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15,valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finaladrenalgenes))
		{
		if(! all(is.na(finaladrenalgenes[,c(i)]))){
		valuepc0 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + log2sum[,1])
		valuepc1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,2] + log2sum[,1])
		valuepc2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,3] + log2sum[,1])
		valuein1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,4] + log2sum[,1])
		valuein2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,5] + log2sum[,1])
		valuein3 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,6] + log2sum[,1])
		valuein4 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,7] + log2sum[,1])
		valuein5 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,8] + log2sum[,1])
		valuein6 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,9] + log2sum[,1])
		valuein7 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,10] + log2sum[,1])
		valuein8 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,11] + log2sum[,1])
		valuein9 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,12] + log2sum[,1])
		valuein10 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,13] + log2sum[,1])
		valuein11 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,14] + log2sum[,1])
		valuein12 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,15] + log2sum[,1])
		valuein13 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,16] + log2sum[,1])
		valuein14 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,17] + log2sum[,1])
		valuein15 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,18] + log2sum[,1])
		valueall <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15,valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finaladrenalgenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2 <- effectsize$Effects[1,"pEta-sqr"]
      } else {
        partialeta2 <- NA
		}
		combo<- rbind(pval, partialeta2)}
			pvalues <- cbind(pvalues, combo) 
			}

		pvalues<- data.frame(pvalues)
		tpvalues<- t(pvalues)
		skip=NULL

		
		for (i in colnames(finaladrenalgenes))
		{
			skipped=NULL
			if(!all(is.na(finaladrenalgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	setwd("/home/quints1/Correlations/adrenal")
	write.table(finalAICc, "All AICc results redo.txt", sep="\t", quote=FALSE)
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finaladrenalgenes))
		{
		if(! all(is.na(finaladrenalgenes[i]))){
		valuepc0 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + log2sum[,1])
		valuepc1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,2] + log2sum[,1])
		valuepc2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,3] + log2sum[,1])
		valuein1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,4] + log2sum[,1])
		valuein2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,5] + log2sum[,1])
		valuein3 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,6] + log2sum[,1])
		valuein4 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,7] + log2sum[,1])
		valuein5 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,8] + log2sum[,1])
		valuein6 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,9] + log2sum[,1])
		valuein7 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,10] + log2sum[,1])
		valuein8 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,11] + log2sum[,1])
		valuein9 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,12] + log2sum[,1])
		valuein10 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,13] + log2sum[,1])
		valuein11 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,14] + log2sum[,1])
		valuein12 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,15] + log2sum[,1])
		valuein13 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,16] + log2sum[,1])
		valuein14 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,17] + log2sum[,1])
		valuein15 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,18] + log2sum[,1])
		valueall <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15,valueall)
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
		
		
		for (i in colnames(finaladrenalgenes))
		{
			skipped=NULL
			if(!all(is.na(finaladrenalgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}	
			skip<- t(skip)
			rownames(trvalues)= skip[,1]
			
		#qvalues
		testforq<- tpvalues[,1]
		testforq[!is.finite(testforq)] <- NA
		testforq<- data.frame(na.omit(testforq))
		qvalues<- data.frame(qvalue(testforq[,1], pi0=1)$qvalues)
		rownames(qvalues)=rownames(tpvalues)
		results<- merge(tpvalues, trvalues, by.x="row.names", by.y="row.names")	
		rownames(results)=results[,1]
		results<- results[,2:4]
		colnames(results)=c("pvalues", "Partialeta", "rsqvalues")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "Partialeta", "rsqvalues", "qvalues")
		setwd("/home/quints1/Correlations/adrenal")
		write.table(finalresults, "P and Q Adrenal redo.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finaladrenalgenes))
		{
		if(! all(is.na(finaladrenalgenes[i]))){
			pearson<- cor.test(finaladrenalgenes[,c(i)], finaladrenalL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finaladrenalgenes[,c(i)], finaladrenalL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finaladrenalgenes))
		{
			skipped=NULL
			if(!all(is.na(finaladrenalgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}	
			skip<- t(skip)
			rownames(trvalueresults)= skip[,1]
				
		results<- merge(tpvalues, trvalueresults, by.x="row.names", by.y="row.names")	
		rownames(results)=results[,1]
		results<- results[,2:5]
		colnames(results)=c("pvalues", "partialeta2", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/adrenal")
		write.table(finalresults, "PRandQAdrenal redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
