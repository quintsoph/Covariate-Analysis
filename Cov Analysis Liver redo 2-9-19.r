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
library("stringi")
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
	rownames(livall)=livall[,1]
	livall<- livall[,-c(1,2)]
	livallcombo<-  merge(livall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(livallcombo)=livallcombo[,1]
	finallivcov<- data.frame(livallcombo[,c(56205:56222)])	
	rownames(finallivcov)=livallcombo[,1]
	finallivL1HS<- data.frame(livallcombo[,c(56204)])
	rownames(finallivL1HS)= livallcombo[,1]
	colnames(finallivL1HS)= c("L1HS")
	finallivgenes<- data.frame(livallcombo[,c(2:56203)])
	log2sum<- data.frame(livallcombo[,56223])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(livallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finallivgenes))
		{
		valuepc0 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + log2sum[,1])
		valuepc1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,2] + log2sum[,1])
		valuepc2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,3] + log2sum[,1])
		valuein1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,4] + log2sum[,1])
		valuein2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,5] + log2sum[,1])
		valuein3 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,6] + log2sum[,1])
		valuein4 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,7] + log2sum[,1])
		valuein5 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,8] + log2sum[,1])
		valuein6 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,9] + log2sum[,1])
		valuein7 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,10] + log2sum[,1])
		valuein8 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,11] + log2sum[,1])
		valuein9 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,12] + log2sum[,1])
		valuein10 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,13] + log2sum[,1])
		valuein11 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,14] + log2sum[,1])
		valuein12 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,15] + log2sum[,1])
		valuein13 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,16] + log2sum[,1])
		valuein14 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,17] + log2sum[,1])
		valuein15 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,18] + log2sum[,1])
		valueall <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + finallivcov[,2] + 
					finallivcov[,3] + finallivcov[,4] + finallivcov[,5] + finallivcov[,6] + finallivcov[,7] + 
					finallivcov[,8] + finallivcov[,9] + finallivcov[,10] + finallivcov[,11] + finallivcov[,12] +
					finallivcov[,13] + finallivcov[,14] + finallivcov[,15] + finallivcov[,16] + finallivcov[,17] + finallivcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/liver")
	write.table(finalAICc, "AllAICcresultsliver redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finallivgenes))
		{
		if(! all(is.na(finallivgenes[i]))){
		valuepc0 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + log2sum[,1])
		valuepc1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,2] + log2sum[,1])
		valuepc2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,3] + log2sum[,1])
		valuein1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,4] + log2sum[,1])
		valuein2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,5] + log2sum[,1])
		valuein3 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,6] + log2sum[,1])
		valuein4 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,7] + log2sum[,1])
		valuein5 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,8] + log2sum[,1])
		valuein6 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,9] + log2sum[,1])
		valuein7 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,10] + log2sum[,1])
		valuein8 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,11] + log2sum[,1])
		valuein9 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,12] + log2sum[,1])
		valuein10 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,13] + log2sum[,1])
		valuein11 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,14] + log2sum[,1])
		valuein12 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,15] + log2sum[,1])
		valuein13 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,16] + log2sum[,1])
		valuein14 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,17] + log2sum[,1])
		valuein15 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,18] + log2sum[,1])
		valueall <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + finallivcov[,2] + 
					finallivcov[,3] + finallivcov[,4] + finallivcov[,5] + finallivcov[,6] + finallivcov[,7] + 
					finallivcov[,8] + finallivcov[,9] + finallivcov[,10] + finallivcov[,11] + finallivcov[,12] +
					finallivcov[,13] + finallivcov[,14] + finallivcov[,15] + finallivcov[,16] + finallivcov[,17] + finallivcov[,18] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finallivgenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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
		valuepc0 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + log2sum[,1])
		valuepc1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,2] + log2sum[,1])
		valuepc2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,3] + log2sum[,1])
		valuein1 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,4] + log2sum[,1])
		valuein2 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,5] + log2sum[,1])
		valuein3 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,6] + log2sum[,1])
		valuein4 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,7] + log2sum[,1])
		valuein5 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,8] + log2sum[,1])
		valuein6 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,9] + log2sum[,1])
		valuein7 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,10] + log2sum[,1])
		valuein8 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,11] + log2sum[,1])
		valuein9 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,12] + log2sum[,1])
		valuein10 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,13] + log2sum[,1])
		valuein11 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,14] + log2sum[,1])
		valuein12 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,15] + log2sum[,1])
		valuein13 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,16] + log2sum[,1])
		valuein14 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,17] + log2sum[,1])
		valuein15 <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,18] + log2sum[,1])
		valueall <-lm(finallivL1HS[,1] ~ finallivgenes[,c(i)] + finallivcov[,1] + finallivcov[,2] + 
					finallivcov[,3] + finallivcov[,4] + finallivcov[,5] + finallivcov[,6] + finallivcov[,7] + 
					finallivcov[,8] + finallivcov[,9] + finallivcov[,10] + finallivcov[,11] + finallivcov[,12] +
					finallivcov[,13] + finallivcov[,14] + finallivcov[,15] + finallivcov[,16] + finallivcov[,17] + finallivcov[,18] + log2sum[,1])
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
		setwd("/home/quints1/Correlations/liver")
		write.table(finalresults, "P and Q liver redo.txt", sep="\t", quote=FALSE)
		
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
		results<- results[,2:5]
		colnames(results)=c("pvalues", "partialeta2", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/liver")
		write.table(finalresults, "PRandQLiver redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		