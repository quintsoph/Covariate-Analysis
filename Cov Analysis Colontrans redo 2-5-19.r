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
#colon - Transverse
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	colonpats<- read.table("Colonedits.txt", header=TRUE, sep="\t")
	colontranspats<- subset(colonpats, colonpats[,3] == "Colon - Transverse", select=c(1,2,3))
	rownames(colontranspats)=str_replace_all(colontranspats[,1], "-", ".")
	colontransgenes<- merge(colontranspats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	colontransgenesfinal<- colontransgenes[,c(4:56207)]
	rownames(colontransgenesfinal)=colontransgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/colontrans")
		colontranscov<- read.table("colontranseditcov.txt", header=TRUE, row.names=1, sep="\t")
	colontransall<- merge(colontransgenesfinal, colontranscov, by.x="row.names", by.y="row.names")
	rownames(colontransall)=colontransall[,1]
	colontransall<- colontransall[,c(3:56236)]
	colontransallcombo<-  merge(colontransall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(colontransallcombo)=colontransallcombo[,1]
	finalcolontranscov<- data.frame(colontransallcombo[,c(56205:56235)])	
	rownames(finalcolontranscov)=colontransallcombo[,1]
	finalcolontransL1HS<- data.frame(colontransallcombo[,c(56204)])
	rownames(finalcolontransL1HS)= colontransallcombo[,1]
	colnames(finalcolontransL1HS)= c("L1HS")
	finalcolontransgenes<- data.frame(colontransallcombo[,c(2:56203)])
	log2sum<- data.frame(colontransallcombo[,56236])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(colontransallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
#AICc
	finalAICc=NULL
		for (i in colnames(finalcolontransgenes))
		{
		valuepc0 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,1] + log2sum[,1])
		valuepc1 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,2] + log2sum[,1])
		valuepc2 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,3] + log2sum[,1])
		valuein1 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,4] + log2sum[,1])
		valuein2 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,5] + log2sum[,1])
		valuein3 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,6] + log2sum[,1])
		valuein4 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,7] + log2sum[,1])
		valuein5 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,8] + log2sum[,1])
		valuein6 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,9] + log2sum[,1])
		valuein7 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,10] + log2sum[,1])
		valuein8 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,11] + log2sum[,1])
		valuein9 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,12] + log2sum[,1])
		valuein10 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,13] + log2sum[,1])
		valuein11 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,14] + log2sum[,1])
		valuein12 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,15] + log2sum[,1])
		valuein13 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,16] + log2sum[,1])
		valuein14 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,17] + log2sum[,1])
		valuein15 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,18] + log2sum[,1])
		valuein16 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,19] + log2sum[,1])
		valuein17 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,20] + log2sum[,1])
		valuein18 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,21] + log2sum[,1])
		valuein19 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,22] + log2sum[,1])
		valuein20 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,23] + log2sum[,1])
		valuein21 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,24] + log2sum[,1])
		valuein22 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,25] + log2sum[,1])
		valuein23 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,26] + log2sum[,1])
		valuein24 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,27] + log2sum[,1])
		valuein25 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,28] + log2sum[,1])
		valuein26 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,29] + log2sum[,1])
		valuein27 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,30] + log2sum[,1])
		valuein28 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,31] + log2sum[,1])
		valueall <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,1] + finalcolontranscov[,2] + 
					finalcolontranscov[,3] + finalcolontranscov[,4] + finalcolontranscov[,5] + finalcolontranscov[,6] + finalcolontranscov[,7] + 
					finalcolontranscov[,8] + finalcolontranscov[,9] + finalcolontranscov[,10] + finalcolontranscov[,11] + finalcolontranscov[,12] +
					finalcolontranscov[,13] + finalcolontranscov[,14] + finalcolontranscov[,15] + finalcolontranscov[,16] + finalcolontranscov[,17] + finalcolontranscov[,18]+
					finalcolontranscov[,19] + finalcolontranscov[,20] + finalcolontranscov[,21] + finalcolontranscov[,22] + finalcolontranscov[,23] + finalcolontranscov[,24] +
					finalcolontranscov[,25] + finalcolontranscov[,26]+ finalcolontranscov[,27] + finalcolontranscov[,28] + finalcolontranscov[,29] + finalcolontranscov[,30] + 
					finalcolontranscov[,31] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, 
			valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/colontrans")
	write.table(finalAICc, "AllAICcresultscolontrans redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalcolontransgenes))
		{
		if(! all(is.na(finalcolontransgenes[i]))){
		valuepc0 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,1] + log2sum[,1])
		valuepc1 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,2] + log2sum[,1])
		valuepc2 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,3] + log2sum[,1])
		valuein1 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,4] + log2sum[,1])
		valuein2 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,5] + log2sum[,1])
		valuein3 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,6] + log2sum[,1])
		valuein4 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,7] + log2sum[,1])
		valuein5 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,8] + log2sum[,1])
		valuein6 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,9] + log2sum[,1])
		valuein7 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,10] + log2sum[,1])
		valuein8 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,11] + log2sum[,1])
		valuein9 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,12] + log2sum[,1])
		valuein10 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,13] + log2sum[,1])
		valuein11 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,14] + log2sum[,1])
		valuein12 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,15] + log2sum[,1])
		valuein13 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,16] + log2sum[,1])
		valuein14 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,17] + log2sum[,1])
		valuein15 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,18] + log2sum[,1])
		valuein16 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,19] + log2sum[,1])
		valuein17 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,20] + log2sum[,1])
		valuein18 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,21] + log2sum[,1])
		valuein19 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,22] + log2sum[,1])
		valuein20 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,23] + log2sum[,1])
		valuein21 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,24] + log2sum[,1])
		valuein22 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,25] + log2sum[,1])
		valuein23 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,26] + log2sum[,1])
		valuein24 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,27] + log2sum[,1])
		valuein25 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,28] + log2sum[,1])
		valuein26 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,29] + log2sum[,1])
		valuein27 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,30] + log2sum[,1])
		valuein28 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,31] + log2sum[,1])
		valueall <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,1] + finalcolontranscov[,2] + 
					finalcolontranscov[,3] + finalcolontranscov[,4] + finalcolontranscov[,5] + finalcolontranscov[,6] + finalcolontranscov[,7] + 
					finalcolontranscov[,8] + finalcolontranscov[,9] + finalcolontranscov[,10] + finalcolontranscov[,11] + finalcolontranscov[,12] +
					finalcolontranscov[,13] + finalcolontranscov[,14] + finalcolontranscov[,15] + finalcolontranscov[,16] + finalcolontranscov[,17] + finalcolontranscov[,18]+
					finalcolontranscov[,19] + finalcolontranscov[,20] + finalcolontranscov[,21] + finalcolontranscov[,22] + finalcolontranscov[,23] + finalcolontranscov[,24] +
					finalcolontranscov[,25] + finalcolontranscov[,26]+ finalcolontranscov[,27] + finalcolontranscov[,28] + finalcolontranscov[,29] + finalcolontranscov[,30] + 
					finalcolontranscov[,31] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, 
			valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finalcolontransgenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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

		for (i in colnames(finalcolontransgenes))
		{
			skipped=NULL
			if(!all(is.na(finalcolontransgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalcolontransgenes))
		{
		if(! all(is.na(finalcolontransgenes[i]))){
		valuepc0 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,1] + log2sum[,1])
		valuepc1 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,2] + log2sum[,1])
		valuepc2 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,3] + log2sum[,1])
		valuein1 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,4] + log2sum[,1])
		valuein2 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,5] + log2sum[,1])
		valuein3 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,6] + log2sum[,1])
		valuein4 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,7] + log2sum[,1])
		valuein5 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,8] + log2sum[,1])
		valuein6 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,9] + log2sum[,1])
		valuein7 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,10] + log2sum[,1])
		valuein8 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,11] + log2sum[,1])
		valuein9 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,12] + log2sum[,1])
		valuein10 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,13] + log2sum[,1])
		valuein11 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,14] + log2sum[,1])
		valuein12 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,15] + log2sum[,1])
		valuein13 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,16] + log2sum[,1])
		valuein14 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,17] + log2sum[,1])
		valuein15 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,18] + log2sum[,1])
		valuein16 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,19] + log2sum[,1])
		valuein17 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,20] + log2sum[,1])
		valuein18 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,21] + log2sum[,1])
		valuein19 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,22] + log2sum[,1])
		valuein20 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,23] + log2sum[,1])
		valuein21 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,24] + log2sum[,1])
		valuein22 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,25] + log2sum[,1])
		valuein23 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,26] + log2sum[,1])
		valuein24 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,27] + log2sum[,1])
		valuein25 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,28] + log2sum[,1])
		valuein26 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,29] + log2sum[,1])
		valuein27 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,30] + log2sum[,1])
		valuein28 <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,31] + log2sum[,1])
		valueall <-lm(finalcolontransL1HS[,1] ~ finalcolontransgenes[,c(i)] + finalcolontranscov[,1] + finalcolontranscov[,2] + 
					finalcolontranscov[,3] + finalcolontranscov[,4] + finalcolontranscov[,5] + finalcolontranscov[,6] + finalcolontranscov[,7] + 
					finalcolontranscov[,8] + finalcolontranscov[,9] + finalcolontranscov[,10] + finalcolontranscov[,11] + finalcolontranscov[,12] +
					finalcolontranscov[,13] + finalcolontranscov[,14] + finalcolontranscov[,15] + finalcolontranscov[,16] + finalcolontranscov[,17] + finalcolontranscov[,18]+
					finalcolontranscov[,19] + finalcolontranscov[,20] + finalcolontranscov[,21] + finalcolontranscov[,22] + finalcolontranscov[,23] + finalcolontranscov[,24] +
					finalcolontranscov[,25] + finalcolontranscov[,26]+ finalcolontranscov[,27] + finalcolontranscov[,28] + finalcolontranscov[,29] + finalcolontranscov[,30] + 
					finalcolontranscov[,31] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, 
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
		
		
		for (i in colnames(finalcolontransgenes))
		{
			skipped=NULL
			if(!all(is.na(finalcolontransgenes[i]))){
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
		setwd("/home/quints1/Correlations/colontrans")
		write.table(finalresults, "P and Q colontrans redo.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalcolontransgenes))
		{
		if(! all(is.na(finalcolontransgenes[i]))){
			pearson<- cor.test(finalcolontransgenes[,c(i)], finalcolontransL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalcolontransgenes[,c(i)], finalcolontransL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalcolontransgenes))
		{
			skipped=NULL
			if(!all(is.na(finalcolontransgenes[i]))){
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
		setwd("/home/quints1/Correlations/colontrans")
		write.table(finalresults, "PRandQColonTrans redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		
	