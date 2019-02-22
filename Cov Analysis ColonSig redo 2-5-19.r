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
	rownames(colonsigall)=colonsigall[,1]
	colonsigall<- colonsigall[,-c(1,2)]
	colonsigallcombo<-  merge(colonsigall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(colonsigallcombo)=colonsigallcombo[,1]
	finalcolonsigcov<- data.frame(colonsigallcombo[,c(56205:56237)])	
	rownames(finalcolonsigcov)=colonsigallcombo[,1]
	finalcolonsigL1HS<- data.frame(colonsigallcombo[,c(56204)])
	rownames(finalcolonsigL1HS)= colonsigallcombo[,1]
	colnames(finalcolonsigL1HS)= c("L1HS")
	finalcolonsiggenes<- data.frame(colonsigallcombo[,c(2:56203)])
	log2sum<- data.frame(colonsigallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(colonsigallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)

	#AICc
	finalAICc=NULL
		for (i in colnames(finalcolonsiggenes))
		{
		valuepc0 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,3] + log2sum[,1])
		valuein1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,4] + log2sum[,1])
		valuein2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,5] + log2sum[,1])
		valuein3 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,6] + log2sum[,1])
		valuein4 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,7] + log2sum[,1])
		valuein5 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,8] + log2sum[,1])
		valuein6 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,9] + log2sum[,1])
		valuein7 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,10] + log2sum[,1])
		valuein8 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,11] + log2sum[,1])
		valuein9 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,12] + log2sum[,1])
		valuein10 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,13] + log2sum[,1])
		valuein11 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,14] + log2sum[,1])
		valuein12 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,15] + log2sum[,1])
		valuein13 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,16] + log2sum[,1])
		valuein14 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,17] + log2sum[,1])
		valuein15 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,18] + log2sum[,1])
		valuein16 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,19] + log2sum[,1])
		valuein17 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,20] + log2sum[,1])
		valuein18 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,21] + log2sum[,1])
		valuein19 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,22] + log2sum[,1])
		valuein20 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,23] + log2sum[,1])
		valuein21 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,24] + log2sum[,1])
		valuein22 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,25] + log2sum[,1])
		valuein23 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,26] + log2sum[,1])
		valuein24 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,27] + log2sum[,1])
		valuein25 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,28] + log2sum[,1])
		valuein26 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,29] + log2sum[,1])
		valuein27 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,30] + log2sum[,1])
		valuein28 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,31] + log2sum[,1])
		valuein29 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,32] + log2sum[,1])
		valuein30 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,33] + log2sum[,1])
		valueall <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + finalcolonsigcov[,2] + 
					finalcolonsigcov[,3] + finalcolonsigcov[,4] + finalcolonsigcov[,5] + finalcolonsigcov[,6] + finalcolonsigcov[,7] + 
					finalcolonsigcov[,8] + finalcolonsigcov[,9] + finalcolonsigcov[,10] + finalcolonsigcov[,11] + finalcolonsigcov[,12] +
					finalcolonsigcov[,13] + finalcolonsigcov[,14] + finalcolonsigcov[,15] + finalcolonsigcov[,16] + finalcolonsigcov[,17] + finalcolonsigcov[,18]+
					finalcolonsigcov[,19] + finalcolonsigcov[,20] + finalcolonsigcov[,21] + finalcolonsigcov[,22] + finalcolonsigcov[,23] + finalcolonsigcov[,24] +
					finalcolonsigcov[,25] + finalcolonsigcov[,26]+ finalcolonsigcov[,27] + finalcolonsigcov[,28] + finalcolonsigcov[,29] + finalcolonsigcov[,30] + 
					finalcolonsigcov[,31] + finalcolonsigcov[,32] + finalcolonsigcov[,33] + log2sum[,1])
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
	write.table(finalAICc, "AllAICcresultscolonsig redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalcolonsiggenes))
		{
		if(! all(is.na(finalcolonsiggenes[i]))){
		valuepc0 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,3] + log2sum[,1])
		valuein1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,4] + log2sum[,1])
		valuein2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,5] + log2sum[,1])
		valuein3 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,6] + log2sum[,1])
		valuein4 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,7] + log2sum[,1])
		valuein5 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,8] + log2sum[,1])
		valuein6 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,9] + log2sum[,1])
		valuein7 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,10] + log2sum[,1])
		valuein8 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,11] + log2sum[,1])
		valuein9 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,12] + log2sum[,1])
		valuein10 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,13] + log2sum[,1])
		valuein11 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,14] + log2sum[,1])
		valuein12 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,15] + log2sum[,1])
		valuein13 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,16] + log2sum[,1])
		valuein14 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,17] + log2sum[,1])
		valuein15 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,18] + log2sum[,1])
		valuein16 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,19] + log2sum[,1])
		valuein17 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,20] + log2sum[,1])
		valuein18 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,21] + log2sum[,1])
		valuein19 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,22] + log2sum[,1])
		valuein20 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,23] + log2sum[,1])
		valuein21 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,24] + log2sum[,1])
		valuein22 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,25] + log2sum[,1])
		valuein23 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,26] + log2sum[,1])
		valuein24 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,27] + log2sum[,1])
		valuein25 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,28] + log2sum[,1])
		valuein26 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,29] + log2sum[,1])
		valuein27 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,30] + log2sum[,1])
		valuein28 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,31] + log2sum[,1])
		valuein29 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,32] + log2sum[,1])
		valuein30 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,33] + log2sum[,1])
		valueall <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + finalcolonsigcov[,2] + 
					finalcolonsigcov[,3] + finalcolonsigcov[,4] + finalcolonsigcov[,5] + finalcolonsigcov[,6] + finalcolonsigcov[,7] + 
					finalcolonsigcov[,8] + finalcolonsigcov[,9] + finalcolonsigcov[,10] + finalcolonsigcov[,11] + finalcolonsigcov[,12] +
					finalcolonsigcov[,13] + finalcolonsigcov[,14] + finalcolonsigcov[,15] + finalcolonsigcov[,16] + finalcolonsigcov[,17] + finalcolonsigcov[,18]+
					finalcolonsigcov[,19] + finalcolonsigcov[,20] + finalcolonsigcov[,21] + finalcolonsigcov[,22] + finalcolonsigcov[,23] + finalcolonsigcov[,24] +
					finalcolonsigcov[,25] + finalcolonsigcov[,26]+ finalcolonsigcov[,27] + finalcolonsigcov[,28] + finalcolonsigcov[,29] + finalcolonsigcov[,30] + 
					finalcolonsigcov[,31] + finalcolonsigcov[,32] + finalcolonsigcov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finalcolonsiggenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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
		valuepc0 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + log2sum[,1])
		valuepc1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,2] + log2sum[,1])
		valuepc2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,3] + log2sum[,1])
		valuein1 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,4] + log2sum[,1])
		valuein2 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,5] + log2sum[,1])
		valuein3 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,6] + log2sum[,1])
		valuein4 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,7] + log2sum[,1])
		valuein5 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,8] + log2sum[,1])
		valuein6 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,9] + log2sum[,1])
		valuein7 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,10] + log2sum[,1])
		valuein8 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,11] + log2sum[,1])
		valuein9 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,12] + log2sum[,1])
		valuein10 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,13] + log2sum[,1])
		valuein11 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,14] + log2sum[,1])
		valuein12 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,15] + log2sum[,1])
		valuein13 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,16] + log2sum[,1])
		valuein14 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,17] + log2sum[,1])
		valuein15 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,18] + log2sum[,1])
		valuein16 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,19] + log2sum[,1])
		valuein17 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,20] + log2sum[,1])
		valuein18 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,21] + log2sum[,1])
		valuein19 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,22] + log2sum[,1])
		valuein20 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,23] + log2sum[,1])
		valuein21 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,24] + log2sum[,1])
		valuein22 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,25] + log2sum[,1])
		valuein23 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,26] + log2sum[,1])
		valuein24 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,27] + log2sum[,1])
		valuein25 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,28] + log2sum[,1])
		valuein26 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,29] + log2sum[,1])
		valuein27 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,30] + log2sum[,1])
		valuein28 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,31] + log2sum[,1])
		valuein29 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,32] + log2sum[,1])
		valuein30 <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,33] + log2sum[,1])
		valueall <-lm(finalcolonsigL1HS[,1] ~ finalcolonsiggenes[,c(i)] + finalcolonsigcov[,1] + finalcolonsigcov[,2] + 
					finalcolonsigcov[,3] + finalcolonsigcov[,4] + finalcolonsigcov[,5] + finalcolonsigcov[,6] + finalcolonsigcov[,7] + 
					finalcolonsigcov[,8] + finalcolonsigcov[,9] + finalcolonsigcov[,10] + finalcolonsigcov[,11] + finalcolonsigcov[,12] +
					finalcolonsigcov[,13] + finalcolonsigcov[,14] + finalcolonsigcov[,15] + finalcolonsigcov[,16] + finalcolonsigcov[,17] + finalcolonsigcov[,18]+
					finalcolonsigcov[,19] + finalcolonsigcov[,20] + finalcolonsigcov[,21] + finalcolonsigcov[,22] + finalcolonsigcov[,23] + finalcolonsigcov[,24] +
					finalcolonsigcov[,25] + finalcolonsigcov[,26]+ finalcolonsigcov[,27] + finalcolonsigcov[,28] + finalcolonsigcov[,29] + finalcolonsigcov[,30] + 
					finalcolonsigcov[,31] + finalcolonsigcov[,32] + finalcolonsigcov[,33] + log2sum[,1])
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
		setwd("/home/quints1/Correlations/colonsig")
		write.table(finalresults, "P and Q colonsig redo.txt", sep="\t", quote=FALSE)
		
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
		results<- results[,2:5]
		colnames(results)=c("pvalues", "partialeta2", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/colonsig")
		write.table(finalresults, "PRandQColonSig redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		