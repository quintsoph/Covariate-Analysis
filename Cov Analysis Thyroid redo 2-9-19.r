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
#thyroid
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	thypats<- read.table("Thyroidedits.txt", header=TRUE, sep="\t")
	rownames(thypats)=str_replace_all(thypats[,1], "-", ".")
	thygenes<- merge(thypats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	thygenesfinal<- thygenes[,c(4:56207)]
	rownames(thygenesfinal)=thygenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/thyroid")
		thycov<- read.table("thyroideditcov.txt", header=TRUE, row.names=1, sep="\t")
	thyall<- merge(thygenesfinal, thycov, by.x="row.names", by.y="row.names")
	rownames(thyall)=thyall[,1]
	thyall<- thyall[,-c(1,2)]
	thyallcombo<-  merge(thyall, tlog2sub, by.x="row.names", by.y="row.names")
	rownames(thyallcombo)=thyallcombo[,1]
	finalthycov<- data.frame(thyallcombo[,c(56205:56252)])	
	rownames(finalthycov)=thyallcombo[,1]
	finalthyL1HS<- data.frame(thyallcombo[,c(56204)])
	rownames(finalthyL1HS)= thyallcombo[,1]
	colnames(finalthyL1HS)= c("L1HS")
	finalthygenes<- data.frame(thyallcombo[,c(2:56203)])
	log2sum<- data.frame(thyallcombo[,56253])
	colnames(log2sum)=c("log2sum")
	rownames(log2sum)=rownames(thyallcombo)
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalthygenes))
		{
		valuepc0 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,1] + log2sum[,1])
		valuepc1 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,2] + log2sum[,1])
		valuepc2 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,3] + log2sum[,1])
		valuein1 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,4] + log2sum[,1])
		valuein2 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,5] + log2sum[,1])
		valuein3 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,6] + log2sum[,1])
		valuein4 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,7] + log2sum[,1])
		valuein5 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,8] + log2sum[,1])
		valuein6 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,9] + log2sum[,1])
		valuein7 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,10] + log2sum[,1])
		valuein8 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,11] + log2sum[,1])
		valuein9 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,12] + log2sum[,1])
		valuein10 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,13] + log2sum[,1])
		valuein11 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,14] + log2sum[,1])
		valuein12 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,15] + log2sum[,1])
		valuein13 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,16] + log2sum[,1])
		valuein14 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,17] + log2sum[,1])
		valuein15 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,18] + log2sum[,1])
		valuein16 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,19] + log2sum[,1])
		valuein17 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,20] + log2sum[,1])
		valuein18 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,21] + log2sum[,1])
		valuein19 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,22] + log2sum[,1])
		valuein20 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,23] + log2sum[,1])
		valuein21 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,24] + log2sum[,1])
		valuein22 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,25] + log2sum[,1])
		valuein23 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,26] + log2sum[,1])
		valuein24 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,27] + log2sum[,1])
		valuein25 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,28] + log2sum[,1])
		valuein26 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,29] + log2sum[,1])
		valuein27 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,30] + log2sum[,1])
		valuein28 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,31] + log2sum[,1])
		valuein29 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,32] + log2sum[,1])
		valuein30 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,33] + log2sum[,1])
		valuein31 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,34] + log2sum[,1])
		valuein32 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,35] + log2sum[,1])
		valuein33 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,36] + log2sum[,1])
		valuein34 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,37] + log2sum[,1])
		valuein35 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,38] + log2sum[,1])
		valuein36 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,39] + log2sum[,1])
		valuein37 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,40] + log2sum[,1])
		valuein38 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,41] + log2sum[,1])
		valuein39 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,42] + log2sum[,1])
		valuein40 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,43] + log2sum[,1])
		valuein41 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,44] + log2sum[,1])
		valuein42 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,45] + log2sum[,1])
		valuein43 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,46] + log2sum[,1])
		valuein44 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,47] + log2sum[,1])
		valuein45 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,48] + log2sum[,1])
		valueall <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,1] + finalthycov[,2] + 
					finalthycov[,3] + finalthycov[,4] + finalthycov[,5] + finalthycov[,6] + finalthycov[,7] + 
					finalthycov[,8] + finalthycov[,9] + finalthycov[,10] + finalthycov[,11] + finalthycov[,12] +
					finalthycov[,13] + finalthycov[,14] + finalthycov[,15] + finalthycov[,16] + finalthycov[,17] + finalthycov[,18] +
					finalthycov[,19] + finalthycov[,20] + finalthycov[,21] + finalthycov[,22] + finalthycov[,23] + finalthycov[,24] + 
					finalthycov[,25] + finalthycov[,26] + finalthycov[,27] + finalthycov[,28] + finalthycov[,29] + finalthycov[,30] +
					finalthycov[,31] + finalthycov[,32] + finalthycov[,33] + finalthycov[,34] + finalthycov[,35] + finalthycov[,36] +
					finalthycov[,37] + finalthycov[,38] + finalthycov[,39] + finalthycov[,40] + finalthycov[,41] + finalthycov[,42] +
					finalthycov[,43] + finalthycov[,44] + finalthycov[,45] + finalthycov[,46] + finalthycov[,47] + finalthycov[,48] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valuein31, valuein32, valuein33, valuein34, valuein35, valuein36, valuein37, valuein38, valuein39, valuein40, valuein41,
			valuein42, valuein43, valuein44, valuein45, valueall)
		evaluate<- data.frame(AICc[1,])
		result<- cbind(evaluate, rownames(evaluate))
		name<- paste(i, "topmodel", sep="")
		rownames(result)= name
		finalAICc<- rbind(finalAICc, result)
	}
	
	setwd("/home/quints1/Correlations/thyroid")
	write.table(finalAICc, "AllAICcresultsthyroid redo.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalthygenes))
		{
		if(! all(is.na(finalthygenes[i]))){
		valuepc0 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,1] + log2sum[,1])
		valuepc1 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,2] + log2sum[,1])
		valuepc2 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,3] + log2sum[,1])
		valuein1 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,4] + log2sum[,1])
		valuein2 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,5] + log2sum[,1])
		valuein3 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,6] + log2sum[,1])
		valuein4 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,7] + log2sum[,1])
		valuein5 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,8] + log2sum[,1])
		valuein6 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,9] + log2sum[,1])
		valuein7 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,10] + log2sum[,1])
		valuein8 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,11] + log2sum[,1])
		valuein9 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,12] + log2sum[,1])
		valuein10 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,13] + log2sum[,1])
		valuein11 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,14] + log2sum[,1])
		valuein12 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,15] + log2sum[,1])
		valuein13 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,16] + log2sum[,1])
		valuein14 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,17] + log2sum[,1])
		valuein15 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,18] + log2sum[,1])
		valuein16 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,19] + log2sum[,1])
		valuein17 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,20] + log2sum[,1])
		valuein18 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,21] + log2sum[,1])
		valuein19 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,22] + log2sum[,1])
		valuein20 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,23] + log2sum[,1])
		valuein21 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,24] + log2sum[,1])
		valuein22 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,25] + log2sum[,1])
		valuein23 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,26] + log2sum[,1])
		valuein24 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,27] + log2sum[,1])
		valuein25 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,28] + log2sum[,1])
		valuein26 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,29] + log2sum[,1])
		valuein27 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,30] + log2sum[,1])
		valuein28 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,31] + log2sum[,1])
		valuein29 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,32] + log2sum[,1])
		valuein30 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,33] + log2sum[,1])
		valuein31 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,34] + log2sum[,1])
		valuein32 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,35] + log2sum[,1])
		valuein33 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,36] + log2sum[,1])
		valuein34 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,37] + log2sum[,1])
		valuein35 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,38] + log2sum[,1])
		valuein36 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,39] + log2sum[,1])
		valuein37 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,40] + log2sum[,1])
		valuein38 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,41] + log2sum[,1])
		valuein39 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,42] + log2sum[,1])
		valuein40 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,43] + log2sum[,1])
		valuein41 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,44] + log2sum[,1])
		valuein42 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,45] + log2sum[,1])
		valuein43 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,46] + log2sum[,1])
		valuein44 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,47] + log2sum[,1])
		valuein45 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,48] + log2sum[,1])
		valueall <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,1] + finalthycov[,2] + 
					finalthycov[,3] + finalthycov[,4] + finalthycov[,5] + finalthycov[,6] + finalthycov[,7] + 
					finalthycov[,8] + finalthycov[,9] + finalthycov[,10] + finalthycov[,11] + finalthycov[,12] +
					finalthycov[,13] + finalthycov[,14] + finalthycov[,15] + finalthycov[,16] + finalthycov[,17] + finalthycov[,18] +
					finalthycov[,19] + finalthycov[,20] + finalthycov[,21] + finalthycov[,22] + finalthycov[,23] + finalthycov[,24] + 
					finalthycov[,25] + finalthycov[,26] + finalthycov[,27] + finalthycov[,28] + finalthycov[,29] + finalthycov[,30] +
					finalthycov[,31] + finalthycov[,32] + finalthycov[,33] + finalthycov[,34] + finalthycov[,35] + finalthycov[,36] +
					finalthycov[,37] + finalthycov[,38] + finalthycov[,39] + finalthycov[,40] + finalthycov[,41] + finalthycov[,42] +
					finalthycov[,43] + finalthycov[,44] + finalthycov[,45] + finalthycov[,46] + finalthycov[,47] + finalthycov[,48] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valuein31, valuein32, valuein33, valuein34, valuein35, valuein36, valuein37, valuein38, valuein39, valuein40, valuein41,
			valuein42, valuein43, valuein44, valuein45, valueall)
		bestmodel <- eval(getCall(AICc, 1))
		finalresult<- summary(bestmodel)
		pval<- finalresult$coefficients[1,"Pr(>|t|)"]
		if (finalthygenes[,c(i)] != "L1HS" && !any(finalresult$aliased) && deviance(bestmodel) > 1.0e-08) {
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
		
		for (i in colnames(finalthygenes))
		{
			skipped=NULL
			if(!all(is.na(finalthygenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalthygenes))
		{
		if(! all(is.na(finalthygenes[i]))){
		valuepc0 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,1] + log2sum[,1])
		valuepc1 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,2] + log2sum[,1])
		valuepc2 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,3] + log2sum[,1])
		valuein1 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,4] + log2sum[,1])
		valuein2 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,5] + log2sum[,1])
		valuein3 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,6] + log2sum[,1])
		valuein4 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,7] + log2sum[,1])
		valuein5 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,8] + log2sum[,1])
		valuein6 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,9] + log2sum[,1])
		valuein7 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,10] + log2sum[,1])
		valuein8 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,11] + log2sum[,1])
		valuein9 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,12] + log2sum[,1])
		valuein10 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,13] + log2sum[,1])
		valuein11 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,14] + log2sum[,1])
		valuein12 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,15] + log2sum[,1])
		valuein13 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,16] + log2sum[,1])
		valuein14 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,17] + log2sum[,1])
		valuein15 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,18] + log2sum[,1])
		valuein16 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,19] + log2sum[,1])
		valuein17 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,20] + log2sum[,1])
		valuein18 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,21] + log2sum[,1])
		valuein19 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,22] + log2sum[,1])
		valuein20 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,23] + log2sum[,1])
		valuein21 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,24] + log2sum[,1])
		valuein22 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,25] + log2sum[,1])
		valuein23 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,26] + log2sum[,1])
		valuein24 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,27] + log2sum[,1])
		valuein25 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,28] + log2sum[,1])
		valuein26 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,29] + log2sum[,1])
		valuein27 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,30] + log2sum[,1])
		valuein28 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,31] + log2sum[,1])
		valuein29 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,32] + log2sum[,1])
		valuein30 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,33] + log2sum[,1])
		valuein31 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,34] + log2sum[,1])
		valuein32 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,35] + log2sum[,1])
		valuein33 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,36] + log2sum[,1])
		valuein34 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,37] + log2sum[,1])
		valuein35 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,38] + log2sum[,1])
		valuein36 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,39] + log2sum[,1])
		valuein37 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,40] + log2sum[,1])
		valuein38 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,41] + log2sum[,1])
		valuein39 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,42] + log2sum[,1])
		valuein40 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,43] + log2sum[,1])
		valuein41 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,44] + log2sum[,1])
		valuein42 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,45] + log2sum[,1])
		valuein43 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,46] + log2sum[,1])
		valuein44 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,47] + log2sum[,1])
		valuein45 <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,48] + log2sum[,1])
		valueall <-lm(finalthyL1HS[,1] ~ finalthygenes[,c(i)] + finalthycov[,1] + finalthycov[,2] + 
					finalthycov[,3] + finalthycov[,4] + finalthycov[,5] + finalthycov[,6] + finalthycov[,7] + 
					finalthycov[,8] + finalthycov[,9] + finalthycov[,10] + finalthycov[,11] + finalthycov[,12] +
					finalthycov[,13] + finalthycov[,14] + finalthycov[,15] + finalthycov[,16] + finalthycov[,17] + finalthycov[,18] +
					finalthycov[,19] + finalthycov[,20] + finalthycov[,21] + finalthycov[,22] + finalthycov[,23] + finalthycov[,24] + 
					finalthycov[,25] + finalthycov[,26] + finalthycov[,27] + finalthycov[,28] + finalthycov[,29] + finalthycov[,30] +
					finalthycov[,31] + finalthycov[,32] + finalthycov[,33] + finalthycov[,34] + finalthycov[,35] + finalthycov[,36] +
					finalthycov[,37] + finalthycov[,38] + finalthycov[,39] + finalthycov[,40] + finalthycov[,41] + finalthycov[,42] +
					finalthycov[,43] + finalthycov[,44] + finalthycov[,45] + finalthycov[,46] + finalthycov[,47] + finalthycov[,48] + log2sum[,1])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valuein31, valuein32, valuein33, valuein34, valuein35, valuein36, valuein37, valuein38, valuein39, valuein40, valuein41,
			valuein42, valuein43, valuein44, valuein45, valueall)
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
		
		
		for (i in colnames(finalthygenes))
		{
			skipped=NULL
			if(!all(is.na(finalthygenes[i]))){
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
		setwd("/home/quints1/Correlations/thyroid")
		write.table(finalresults, "P and Q thyroid redo.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalthygenes))
		{
		if(! all(is.na(finalthygenes[i]))){
			pearson<- cor.test(finalthygenes[,c(i)], finalthyL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalthygenes[,c(i)], finalthyL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalthygenes))
		{
			skipped=NULL
			if(!all(is.na(finalthygenes[i]))){
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
		setwd("/home/quints1/Correlations/thyroid")
		write.table(finalresults, "PRandQThyroid redo.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		