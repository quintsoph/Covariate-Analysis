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
	#stomach
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	stopats<- read.table("Stomachedits.txt", header=TRUE, sep="\t")
	rownames(stopats)=str_replace_all(stopats[,1], "-", ".")
	stogenes<- merge(stopats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	stogenesfinal<- stogenes[,c(4:56207)]
	rownames(stogenesfinal)=stogenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/stomach")
		stocov<- read.table("stomacheditcov.txt", header=TRUE, row.names=1, sep="\t")
	stoall<- merge(stogenesfinal, stocov, by.x="row.names", by.y="row.names")
	rownames(stoall)=stoall[,1]
	stoall<- stoall[,-c(1,2)]
	stoallcombo<-  merge(stoall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(stoallcombo)=stoallcombo[,1]
	finalstocov<- data.frame(stoallcombo[,c(56205:56237)])	
	rownames(finalstocov)=stoallcombo[,1]
	finalstoL1HS<- data.frame(stoallcombo[,c(56204)])
	rownames(finalstoL1HS)= stoallcombo[,1]
	colnames(finalstoL1HS)= c("L1HS")
	finalstogenes<- data.frame(stoallcombo[,c(2:56203)])
	log2sum<- data.frame(stoallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(stoallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalstogenes))
		{
		valuepc0 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + log2sum[,1])
		valuepc1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,2] + log2sum[,1])
		valuepc2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,3] + log2sum[,1])
		valuein1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,4] + log2sum[,1])
		valuein2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,5] + log2sum[,1])
		valuein3 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,6] + log2sum[,1])
		valuein4 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,7] + log2sum[,1])
		valuein5 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,8] + log2sum[,1])
		valuein6 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,9] + log2sum[,1])
		valuein7 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,10] + log2sum[,1])
		valuein8 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,11] + log2sum[,1])
		valuein9 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,12] + log2sum[,1])
		valuein10 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,13] + log2sum[,1])
		valuein11 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,14] + log2sum[,1])
		valuein12 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,15] + log2sum[,1])
		valuein13 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,16] + log2sum[,1])
		valuein14 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,17] + log2sum[,1])
		valuein15 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,18] + log2sum[,1])
		valuein16 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,19] + log2sum[,1])
		valuein17 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,20] + log2sum[,1])
		valuein18 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,21] + log2sum[,1])
		valuein19 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,22] + log2sum[,1])
		valuein20 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,23] + log2sum[,1])
		valuein21 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,24] + log2sum[,1])
		valuein22 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,25] + log2sum[,1])
		valuein23 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,26] + log2sum[,1])
		valuein24 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,27] + log2sum[,1])
		valuein25 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,28] + log2sum[,1])
		valuein26 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,29] + log2sum[,1])
		valuein27 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,30] + log2sum[,1])
		valuein28 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,31] + log2sum[,1])
		valuein29 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,32] + log2sum[,1])
		valuein30 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,33] + log2sum[,1])
		valueall <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + finalstocov[,2] + 
					finalstocov[,3] + finalstocov[,4] + finalstocov[,5] + finalstocov[,6] + finalstocov[,7] + 
					finalstocov[,8] + finalstocov[,9] + finalstocov[,10] + finalstocov[,11] + finalstocov[,12] +
					finalstocov[,13] + finalstocov[,14] + finalstocov[,15] + finalstocov[,16] + finalstocov[,17] + 
					finalstocov[,18] + finalstocov[,19] + finalstocov[,20] + finalstocov[,21] + finalstocov[,22] + finalstocov[,23] + finalstocov[,24] + 
					finalstocov[,25] + finalstocov[,26] + finalstocov[,27] + finalstocov[,28] + finalstocov[,29] + finalstocov[,30] +
					finalstocov[,31] + finalstocov[,32] + finalstocov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/stomach")
	write.table(finalAICc, "AllAICcresultsstomach redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalstogenes))
		{
		if(! all(is.na(finalstogenes[i]))){
		valuepc0 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + log2sum[,1])
		valuepc1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,2] + log2sum[,1])
		valuepc2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,3] + log2sum[,1])
		valuein1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,4] + log2sum[,1])
		valuein2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,5] + log2sum[,1])
		valuein3 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,6] + log2sum[,1])
		valuein4 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,7] + log2sum[,1])
		valuein5 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,8] + log2sum[,1])
		valuein6 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,9] + log2sum[,1])
		valuein7 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,10] + log2sum[,1])
		valuein8 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,11] + log2sum[,1])
		valuein9 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,12] + log2sum[,1])
		valuein10 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,13] + log2sum[,1])
		valuein11 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,14] + log2sum[,1])
		valuein12 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,15] + log2sum[,1])
		valuein13 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,16] + log2sum[,1])
		valuein14 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,17] + log2sum[,1])
		valuein15 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,18] + log2sum[,1])
		valuein16 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,19] + log2sum[,1])
		valuein17 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,20] + log2sum[,1])
		valuein18 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,21] + log2sum[,1])
		valuein19 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,22] + log2sum[,1])
		valuein20 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,23] + log2sum[,1])
		valuein21 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,24] + log2sum[,1])
		valuein22 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,25] + log2sum[,1])
		valuein23 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,26] + log2sum[,1])
		valuein24 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,27] + log2sum[,1])
		valuein25 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,28] + log2sum[,1])
		valuein26 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,29] + log2sum[,1])
		valuein27 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,30] + log2sum[,1])
		valuein28 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,31] + log2sum[,1])
		valuein29 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,32] + log2sum[,1])
		valuein30 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,33] + log2sum[,1])
		valueall <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + finalstocov[,2] + 
					finalstocov[,3] + finalstocov[,4] + finalstocov[,5] + finalstocov[,6] + finalstocov[,7] + 
					finalstocov[,8] + finalstocov[,9] + finalstocov[,10] + finalstocov[,11] + finalstocov[,12] +
					finalstocov[,13] + finalstocov[,14] + finalstocov[,15] + finalstocov[,16] + finalstocov[,17] + 
					finalstocov[,18] + finalstocov[,19] + finalstocov[,20] + finalstocov[,21] + finalstocov[,22] + finalstocov[,23] + finalstocov[,24] + 
					finalstocov[,25] + finalstocov[,26] + finalstocov[,27] + finalstocov[,28] + finalstocov[,29] + finalstocov[,30] +
					finalstocov[,31] + finalstocov[,32] + finalstocov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finalstogenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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
		
		for (i in colnames(finalstogenes))
		{
			skipped=NULL
			if(!all(is.na(finalstogenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalstogenes))
		{
		if(! all(is.na(finalstogenes[i]))){
		valuepc0 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + log2sum[,1])
		valuepc1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,2] + log2sum[,1])
		valuepc2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,3] + log2sum[,1])
		valuein1 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,4] + log2sum[,1])
		valuein2 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,5] + log2sum[,1])
		valuein3 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,6] + log2sum[,1])
		valuein4 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,7] + log2sum[,1])
		valuein5 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,8] + log2sum[,1])
		valuein6 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,9] + log2sum[,1])
		valuein7 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,10] + log2sum[,1])
		valuein8 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,11] + log2sum[,1])
		valuein9 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,12] + log2sum[,1])
		valuein10 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,13] + log2sum[,1])
		valuein11 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,14] + log2sum[,1])
		valuein12 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,15] + log2sum[,1])
		valuein13 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,16] + log2sum[,1])
		valuein14 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,17] + log2sum[,1])
		valuein15 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,18] + log2sum[,1])
		valuein16 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,19] + log2sum[,1])
		valuein17 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,20] + log2sum[,1])
		valuein18 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,21] + log2sum[,1])
		valuein19 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,22] + log2sum[,1])
		valuein20 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,23] + log2sum[,1])
		valuein21 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,24] + log2sum[,1])
		valuein22 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,25] + log2sum[,1])
		valuein23 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,26] + log2sum[,1])
		valuein24 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,27] + log2sum[,1])
		valuein25 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,28] + log2sum[,1])
		valuein26 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,29] + log2sum[,1])
		valuein27 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,30] + log2sum[,1])
		valuein28 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,31] + log2sum[,1])
		valuein29 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,32] + log2sum[,1])
		valuein30 <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,33] + log2sum[,1])
		valueall <-lm(finalstoL1HS[,1] ~ finalstogenes[,c(i)] + finalstocov[,1] + finalstocov[,2] + 
					finalstocov[,3] + finalstocov[,4] + finalstocov[,5] + finalstocov[,6] + finalstocov[,7] + 
					finalstocov[,8] + finalstocov[,9] + finalstocov[,10] + finalstocov[,11] + finalstocov[,12] +
					finalstocov[,13] + finalstocov[,14] + finalstocov[,15] + finalstocov[,16] + finalstocov[,17] + 
					finalstocov[,18] + finalstocov[,19] + finalstocov[,20] + finalstocov[,21] + finalstocov[,22] + finalstocov[,23] + finalstocov[,24] + 
					finalstocov[,25] + finalstocov[,26] + finalstocov[,27] + finalstocov[,28] + finalstocov[,29] + finalstocov[,30] +
					finalstocov[,31] + finalstocov[,32] + finalstocov[,33] + log2sum[,1])
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
		
		
		for (i in colnames(finalstogenes))
		{
			skipped=NULL
			if(!all(is.na(finalstogenes[i]))){
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
		setwd("/home/quints1/Correlations/stomach")
		write.table(finalresults, "P and Q stomach redo.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalstogenes))
		{
		if(! all(is.na(finalstogenes[i]))){
			pearson<- cor.test(finalstogenes[,c(i)], finalstoL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalstogenes[,c(i)], finalstoL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalstogenes))
		{
			skipped=NULL
			if(!all(is.na(finalstogenes[i]))){
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
		setwd("/home/quints1/Correlations/stomach")
		write.table(finalresults, "PRandQStomach redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		
	