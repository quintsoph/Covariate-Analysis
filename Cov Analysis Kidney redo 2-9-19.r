library("stringi")

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
	library("stringi")
	library("stringr")

	#kidney
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	kidpats<- read.table("Kidneyedits.txt", header=TRUE, sep="\t")
	rownames(kidpats)=str_replace_all(kidpats[,1], "-", ".")
	kidgenes<- merge(kidpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	kidgenesfinal<- kidgenes[,c(4:56207)]
	rownames(kidgenesfinal)=kidgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/kidney")
		kidneycov<- read.table("kidneyeditcov.txt", header=TRUE, row.names=1, sep="\t")
	kidall<- merge(kidgenesfinal, kidneycov, by.x="row.names", by.y="row.names")
	rownames(kidall)=kidall[,1]
	kidall<- kidall[,-c(1,2)]
	kidallcombo<-  merge(kidall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(kidallcombo)=kidallcombo[,1]
	finalkidcov<- data.frame(kidallcombo[,c(56205:56222)])	
	rownames(finalkidcov)=kidallcombo[,1]
	finalkidL1HS<- data.frame(kidallcombo[,c(56204)])
	rownames(finalkidL1HS)= kidallcombo[,1]
	colnames(finalkidL1HS)= c("L1HS")
	finalkidgenes<- data.frame(kidallcombo[,c(2:56203)])
	log2sum<- data.frame(kidallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(kidallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalkidgenes))
		{
		valuepc0 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,3] + log2sum[,1])
		valuein1 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,4] + log2sum[,1])
		valuein2 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,5] + log2sum[,1])
		valuein3 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,6] + log2sum[,1])
		valuein4 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,7] + log2sum[,1])
		valuein5 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,8] + log2sum[,1])
		valuein6 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,9] + log2sum[,1])
		valuein7 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,10] + log2sum[,1])
		valuein8 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,11] + log2sum[,1])
		valuein9 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,12] + log2sum[,1])
		valuein10 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,13] + log2sum[,1])
		valuein11 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,14] + log2sum[,1])
		valuein12 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,15] + log2sum[,1])
		valuein13 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,16] + log2sum[,1])
		valuein14 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,17] + log2sum[,1])
		valuein15 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,18] + log2sum[,1])
		valueall <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,1] + finalkidcov[,2] + 
					finalkidcov[,3] + finalkidcov[,4] + finalkidcov[,5] + finalkidcov[,6] + finalkidcov[,7] + 
					finalkidcov[,8] + finalkidcov[,9] + finalkidcov[,10] + finalkidcov[,11] + finalkidcov[,12] +
					finalkidcov[,13] + finalkidcov[,14] + finalkidcov[,15] + finalkidcov[,16] + finalkidcov[,17] + finalkidcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/kidney")
	write.table(finalAICc, "AllAICcresultskidney redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalkidgenes))
		{
		if(! all(is.na(finalkidgenes[i]))){
		valuepc0 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,3] + log2sum[,1])
		valuein1 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,4] + log2sum[,1])
		valuein2 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,5] + log2sum[,1])
		valuein3 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,6] + log2sum[,1])
		valuein4 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,7] + log2sum[,1])
		valuein5 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,8] + log2sum[,1])
		valuein6 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,9] + log2sum[,1])
		valuein7 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,10] + log2sum[,1])
		valuein8 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,11] + log2sum[,1])
		valuein9 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,12] + log2sum[,1])
		valuein10 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,13] + log2sum[,1])
		valuein11 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,14] + log2sum[,1])
		valuein12 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,15] + log2sum[,1])
		valuein13 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,16] + log2sum[,1])
		valuein14 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,17] + log2sum[,1])
		valuein15 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,18] + log2sum[,1])
		valueall <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,1] + finalkidcov[,2] + 
					finalkidcov[,3] + finalkidcov[,4] + finalkidcov[,5] + finalkidcov[,6] + finalkidcov[,7] + 
					finalkidcov[,8] + finalkidcov[,9] + finalkidcov[,10] + finalkidcov[,11] + finalkidcov[,12] +
					finalkidcov[,13] + finalkidcov[,14] + finalkidcov[,15] + finalkidcov[,16] + finalkidcov[,17] + finalkidcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finalkidgenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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
		
		for (i in colnames(finalkidgenes))
		{
			skipped=NULL
			if(!all(is.na(finalkidgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalkidgenes))
		{
		if(! all(is.na(finalkidgenes[i]))){
		valuepc0 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,3] + log2sum[,1])
		valuein1 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,4] + log2sum[,1])
		valuein2 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,5] + log2sum[,1])
		valuein3 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,6] + log2sum[,1])
		valuein4 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,7] + log2sum[,1])
		valuein5 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,8] + log2sum[,1])
		valuein6 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,9] + log2sum[,1])
		valuein7 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,10] + log2sum[,1])
		valuein8 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,11] + log2sum[,1])
		valuein9 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,12] + log2sum[,1])
		valuein10 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,13] + log2sum[,1])
		valuein11 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,14] + log2sum[,1])
		valuein12 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,15] + log2sum[,1])
		valuein13 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,16] + log2sum[,1])
		valuein14 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,17] + log2sum[,1])
		valuein15 <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,18] + log2sum[,1])
		valueall <-lm(finalkidL1HS[,1] ~ finalkidgenes[,c(i)] + finalkidcov[,1] + finalkidcov[,2] + 
					finalkidcov[,3] + finalkidcov[,4] + finalkidcov[,5] + finalkidcov[,6] + finalkidcov[,7] + 
					finalkidcov[,8] + finalkidcov[,9] + finalkidcov[,10] + finalkidcov[,11] + finalkidcov[,12] +
					finalkidcov[,13] + finalkidcov[,14] + finalkidcov[,15] + finalkidcov[,16] + finalkidcov[,17] + finalkidcov[,18] + log2sum[,1])
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
		
		
		for (i in colnames(finalkidgenes))
		{
			skipped=NULL
			if(!all(is.na(finalkidgenes[i]))){
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
		setwd("/home/quints1/Correlations/kidney")
		write.table(finalresults, "P and Q kidney redo.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalkidgenes))
		{
		if(! all(is.na(finalkidgenes[i]))){
			pearson<- cor.test(finalkidgenes[,c(i)], finalkidL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalkidgenes[,c(i)], finalkidL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalkidgenes))
		{
			skipped=NULL
			if(!all(is.na(finalkidgenes[i]))){
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
		setwd("/home/quints1/Correlations/kidney")
		write.table(finalresults, "PRandQKidney redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		