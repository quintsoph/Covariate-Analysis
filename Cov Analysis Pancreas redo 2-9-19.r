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
	rownames(panall)=panall[,1]
	panall<- panall[,-c(1,2)]
	panallcombo<-  merge(panall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(panallcombo)=panallcombo[,1]
	finalpancov<- data.frame(panallcombo[,c(56205:56237)])	
	rownames(finalpancov)=panallcombo[,1]
	finalpanL1HS<- data.frame(panallcombo[,c(56204)])
	rownames(finalpanL1HS)= panallcombo[,1]
	colnames(finalpanL1HS)= c("L1HS")
	finalpangenes<- data.frame(panallcombo[,c(2:56203)])
	log2sum<- data.frame(panallcombo[,56238])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(panallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
		#AICc
	finalAICc=NULL
		for (i in colnames(finalpangenes))
		{
		valuepc0 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + log2sum[,1])
		valuepc1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,2] + log2sum[,1])
		valuepc2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,3] + log2sum[,1])
		valuein1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,4] + log2sum[,1])
		valuein2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,5] + log2sum[,1])
		valuein3 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,6] + log2sum[,1])
		valuein4 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,7] + log2sum[,1])
		valuein5 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,8] + log2sum[,1])
		valuein6 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,9] + log2sum[,1])
		valuein7 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,10] + log2sum[,1])
		valuein8 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,11] + log2sum[,1])
		valuein9 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,12] + log2sum[,1])
		valuein10 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,13] + log2sum[,1])
		valuein11 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,14] + log2sum[,1])
		valuein12 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,15] + log2sum[,1])
		valuein13 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,16] + log2sum[,1])
		valuein14 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,17] + log2sum[,1])
		valuein15 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,18] + log2sum[,1])
		valuein16 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,19] + log2sum[,1])
		valuein17 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,20] + log2sum[,1])
		valuein18 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,21] + log2sum[,1])
		valuein19 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,22] + log2sum[,1])
		valuein20 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,23] + log2sum[,1])
		valuein21 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,24] + log2sum[,1])
		valuein22 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,25] + log2sum[,1])
		valuein23 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,26] + log2sum[,1])
		valuein24 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,27] + log2sum[,1])
		valuein25 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,28] + log2sum[,1])
		valuein26 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,29] + log2sum[,1])
		valuein27 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,30] + log2sum[,1])
		valuein28 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,31] + log2sum[,1])
		valuein29 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,32] + log2sum[,1])
		valuein30 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,33] + log2sum[,1])
		valueall <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + finalpancov[,2] + 
					finalpancov[,3] + finalpancov[,4] + finalpancov[,5] + finalpancov[,6] + finalpancov[,7] + 
					finalpancov[,8] + finalpancov[,9] + finalpancov[,10] + finalpancov[,11] + finalpancov[,12] +
					finalpancov[,13] + finalpancov[,14] + finalpancov[,15] + finalpancov[,16] + finalpancov[,17] + 
					finalpancov[,18] + finalpancov[,19] + finalpancov[,20] + finalpancov[,21] + finalpancov[,22] + finalpancov[,23] + finalpancov[,24] + 
					finalpancov[,25] + finalpancov[,26] + finalpancov[,27] + finalpancov[,28] + finalpancov[,29] + finalpancov[,30] +
					finalpancov[,31] + finalpancov[,32] + finalpancov[,33] + log2sum[,1])
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
	write.table(finalAICc, "AllAICcresultspancreas redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalpangenes))
		{
		if(! all(is.na(finalpangenes[i]))){
		valuepc0 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + log2sum[,1])
		valuepc1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,2] + log2sum[,1])
		valuepc2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,3] + log2sum[,1])
		valuein1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,4] + log2sum[,1])
		valuein2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,5] + log2sum[,1])
		valuein3 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,6] + log2sum[,1])
		valuein4 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,7] + log2sum[,1])
		valuein5 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,8] + log2sum[,1])
		valuein6 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,9] + log2sum[,1])
		valuein7 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,10] + log2sum[,1])
		valuein8 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,11] + log2sum[,1])
		valuein9 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,12] + log2sum[,1])
		valuein10 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,13] + log2sum[,1])
		valuein11 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,14] + log2sum[,1])
		valuein12 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,15] + log2sum[,1])
		valuein13 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,16] + log2sum[,1])
		valuein14 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,17] + log2sum[,1])
		valuein15 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,18] + log2sum[,1])
		valuein16 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,19] + log2sum[,1])
		valuein17 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,20] + log2sum[,1])
		valuein18 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,21] + log2sum[,1])
		valuein19 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,22] + log2sum[,1])
		valuein20 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,23] + log2sum[,1])
		valuein21 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,24] + log2sum[,1])
		valuein22 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,25] + log2sum[,1])
		valuein23 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,26] + log2sum[,1])
		valuein24 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,27] + log2sum[,1])
		valuein25 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,28] + log2sum[,1])
		valuein26 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,29] + log2sum[,1])
		valuein27 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,30] + log2sum[,1])
		valuein28 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,31] + log2sum[,1])
		valuein29 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,32] + log2sum[,1])
		valuein30 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,33] + log2sum[,1])
		valueall <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + finalpancov[,2] + 
					finalpancov[,3] + finalpancov[,4] + finalpancov[,5] + finalpancov[,6] + finalpancov[,7] + 
					finalpancov[,8] + finalpancov[,9] + finalpancov[,10] + finalpancov[,11] + finalpancov[,12] +
					finalpancov[,13] + finalpancov[,14] + finalpancov[,15] + finalpancov[,16] + finalpancov[,17] + 
					finalpancov[,18] + finalpancov[,19] + finalpancov[,20] + finalpancov[,21] + finalpancov[,22] + finalpancov[,23] + finalpancov[,24] + 
					finalpancov[,25] + finalpancov[,26] + finalpancov[,27] + finalpancov[,28] + finalpancov[,29] + finalpancov[,30] +
					finalpancov[,31] + finalpancov[,32] + finalpancov[,33] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finalpangenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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
		valuepc0 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + log2sum[,1])
		valuepc1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,2] + log2sum[,1])
		valuepc2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,3] + log2sum[,1])
		valuein1 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,4] + log2sum[,1])
		valuein2 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,5] + log2sum[,1])
		valuein3 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,6] + log2sum[,1])
		valuein4 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,7] + log2sum[,1])
		valuein5 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,8] + log2sum[,1])
		valuein6 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,9] + log2sum[,1])
		valuein7 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,10] + log2sum[,1])
		valuein8 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,11] + log2sum[,1])
		valuein9 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,12] + log2sum[,1])
		valuein10 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,13] + log2sum[,1])
		valuein11 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,14] + log2sum[,1])
		valuein12 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,15] + log2sum[,1])
		valuein13 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,16] + log2sum[,1])
		valuein14 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,17] + log2sum[,1])
		valuein15 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,18] + log2sum[,1])
		valuein16 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,19] + log2sum[,1])
		valuein17 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,20] + log2sum[,1])
		valuein18 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,21] + log2sum[,1])
		valuein19 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,22] + log2sum[,1])
		valuein20 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,23] + log2sum[,1])
		valuein21 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,24] + log2sum[,1])
		valuein22 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,25] + log2sum[,1])
		valuein23 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,26] + log2sum[,1])
		valuein24 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,27] + log2sum[,1])
		valuein25 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,28] + log2sum[,1])
		valuein26 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,29] + log2sum[,1])
		valuein27 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,30] + log2sum[,1])
		valuein28 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,31] + log2sum[,1])
		valuein29 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,32] + log2sum[,1])
		valuein30 <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,33] + log2sum[,1])
		valueall <-lm(finalpanL1HS[,1] ~ finalpangenes[,c(i)] + finalpancov[,1] + finalpancov[,2] + 
					finalpancov[,3] + finalpancov[,4] + finalpancov[,5] + finalpancov[,6] + finalpancov[,7] + 
					finalpancov[,8] + finalpancov[,9] + finalpancov[,10] + finalpancov[,11] + finalpancov[,12] +
					finalpancov[,13] + finalpancov[,14] + finalpancov[,15] + finalpancov[,16] + finalpancov[,17] + 
					finalpancov[,18] + finalpancov[,19] + finalpancov[,20] + finalpancov[,21] + finalpancov[,22] + finalpancov[,23] + finalpancov[,24] + 
					finalpancov[,25] + finalpancov[,26] + finalpancov[,27] + finalpancov[,28] + finalpancov[,29] + finalpancov[,30] +
					finalpancov[,31] + finalpancov[,32] + finalpancov[,33] + log2sum[,1])
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
		setwd("/home/quints1/Correlations/pancreas")
		write.table(finalresults, "P and Q pancreas redo.txt", sep="\t", quote=FALSE)
		
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
		results<- results[,2:5]
		colnames(results)=c("pvalues", "partialeta2", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "partialeta2", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/pancreas")
		write.table(finalresults, "PRandQPancreas redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		
	