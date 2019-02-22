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
	
	#esophagus - muscularis
		#Use the old L1HS patients to get correct patient output
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Muscularis", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophgenes<- merge(esophpats, tgenesraw, by.x="row.names", by.y="row.names", drop=FALSE)
	esophgenesfinal<- esophgenes[,c(4:56207)]
	rownames(esophgenesfinal)=esophgenes[,1]
		#read in the covariate
		setwd("/home/quints1/Correlations/esophagusmuscularis")
		esophcov<- read.table("esophagusmusculariseditcov.txt", header=TRUE, row.names=1, sep="\t")
	esophall<- merge(esophgenesfinal, esophcov, by.x="row.names", by.y="row.names")
	finalesophcov<- data.frame(esophall[,c(56206:56253)])	
	rownames(finalesophcov)=esophall[,1]
	finalesophL1HS<- data.frame(esophall[,c(56205)])
	rownames(finalesophL1HS)= esophall[,1]
	colnames(finalesophL1HS)= c("L1HS")
	finalesophgenes<- data.frame(esophall[,c(3:56204)])
	rownames(finalesophgenes)=esophall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
	
	#AICc
	finalAICc=NULL
		for (i in colnames(finalesophgenes))
		{
		valuepc0 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1])
		valuepc1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,2])
		valuepc2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,3])
		valuein1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,4])
		valuein2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,5])
		valuein3 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,6])
		valuein4 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,7])
		valuein5 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,8])
		valuein6 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,9])
		valuein7 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,10])
		valuein8 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,11])
		valuein9 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,12])
		valuein10 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,13])
		valuein11 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,14])
		valuein12 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,15])
		valuein13 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,16])
		valuein14 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,17])
		valuein15 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,18])
		valuein16 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,19])
		valuein17 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,20])
		valuein18 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,21])
		valuein19 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,22])
		valuein20 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,23])
		valuein21 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,24])
		valuein22 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,25])
		valuein23 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,26])
		valuein24 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,27])
		valuein25 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,28])
		valuein26 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,29])
		valuein27 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,30])
		valuein28 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,31])
		valuein29 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,32])
		valuein30 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,33])
		valuein31 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,34])
		valuein32 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,35])
		valuein33 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,36])
		valuein34 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,37])
		valuein35 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,38])
		valuein36 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,39])
		valuein37 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,40])
		valuein38 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,41])
		valuein39 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,42])
		valuein40 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,43])
		valuein41 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,44])
		valuein42 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,45])
		valuein43 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,46])
		valuein44 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,47])
		valuein45 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,48])
		valueall <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1] + finalesophcov[,2] + 
					finalesophcov[,3] + finalesophcov[,4] + finalesophcov[,5] + finalesophcov[,6] + finalesophcov[,7] + 
					finalesophcov[,8] + finalesophcov[,9] + finalesophcov[,10] + finalesophcov[,11] + finalesophcov[,12] +
					finalesophcov[,13] + finalesophcov[,14] + finalesophcov[,15] + finalesophcov[,16] + finalesophcov[,17] + finalesophcov[,18] +
					finalesophcov[,19] + finalesophcov[,20] + finalesophcov[,21] + finalesophcov[,22] + finalesophcov[,23] + finalesophcov[,24] + 
					finalesophcov[,25] + finalesophcov[,26] + finalesophcov[,27] + finalesophcov[,28] + finalesophcov[,29] + finalesophcov[,30] +
					finalesophcov[,31] + finalesophcov[,32] + finalesophcov[,33] + finalesophcov[,34] + finalesophcov[,35] + finalesophcov[,36] +
					finalesophcov[,37] + finalesophcov[,38] + finalesophcov[,39] + finalesophcov[,40] + finalesophcov[,41] + finalesophcov[,42] +
					finalesophcov[,43] + finalesophcov[,44] + finalesophcov[,45] + finalesophcov[,46] + finalesophcov[,47] + finalesophcov[,48])
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
	
	setwd("/home/quints1/Correlations/esophagusmuscularis")
	write.table(finalAICc, "AllAICcresultsesophmuscularis.txt", sep="\t", quote=FALSE)
	
	#final analysis
	pvalues=NULL
	for (i in colnames(finalesophgenes))
		{
		if(! all(is.na(finalesophgenes[i]))){
		valuepc0 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1])
		valuepc1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,2])
		valuepc2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,3])
		valuein1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,4])
		valuein2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,5])
		valuein3 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,6])
		valuein4 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,7])
		valuein5 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,8])
		valuein6 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,9])
		valuein7 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,10])
		valuein8 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,11])
		valuein9 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,12])
		valuein10 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,13])
		valuein11 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,14])
		valuein12 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,15])
		valuein13 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,16])
		valuein14 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,17])
		valuein15 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,18])
		valuein16 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,19])
		valuein17 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,20])
		valuein18 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,21])
		valuein19 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,22])
		valuein20 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,23])
		valuein21 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,24])
		valuein22 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,25])
		valuein23 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,26])
		valuein24 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,27])
		valuein25 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,28])
		valuein26 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,29])
		valuein27 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,30])
		valuein28 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,31])
		valuein29 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,32])
		valuein30 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,33])
		valuein31 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,34])
		valuein32 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,35])
		valuein33 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,36])
		valuein34 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,37])
		valuein35 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,38])
		valuein36 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,39])
		valuein37 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,40])
		valuein38 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,41])
		valuein39 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,42])
		valuein40 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,43])
		valuein41 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,44])
		valuein42 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,45])
		valuein43 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,46])
		valuein44 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,47])
		valuein45 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,48])
		valueall <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1] + finalesophcov[,2] + 
					finalesophcov[,3] + finalesophcov[,4] + finalesophcov[,5] + finalesophcov[,6] + finalesophcov[,7] + 
					finalesophcov[,8] + finalesophcov[,9] + finalesophcov[,10] + finalesophcov[,11] + finalesophcov[,12] +
					finalesophcov[,13] + finalesophcov[,14] + finalesophcov[,15] + finalesophcov[,16] + finalesophcov[,17] + finalesophcov[,18] +
					finalesophcov[,19] + finalesophcov[,20] + finalesophcov[,21] + finalesophcov[,22] + finalesophcov[,23] + finalesophcov[,24] + 
					finalesophcov[,25] + finalesophcov[,26] + finalesophcov[,27] + finalesophcov[,28] + finalesophcov[,29] + finalesophcov[,30] +
					finalesophcov[,31] + finalesophcov[,32] + finalesophcov[,33] + finalesophcov[,34] + finalesophcov[,35] + finalesophcov[,36] +
					finalesophcov[,37] + finalesophcov[,38] + finalesophcov[,39] + finalesophcov[,40] + finalesophcov[,41] + finalesophcov[,42] +
					finalesophcov[,43] + finalesophcov[,44] + finalesophcov[,45] + finalesophcov[,46] + finalesophcov[,47] + finalesophcov[,48])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15, valuein16, valuein17, valuein18, valuein19,
			valuein20, valuein21, valuein22, valuein23, valuein24, valuein25, valuein26, valuein27, valuein28, valuein29, valuein30,
			valuein31, valuein32, valuein33, valuein34, valuein35, valuein36, valuein37, valuein38, valuein39, valuein40, valuein41,
			valuein42, valuein43, valuein44, valuein45, valueall)
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
		
		for (i in colnames(finalesophgenes))
		{
			skipped=NULL
			if(!all(is.na(finalesophgenes[i]))){
			skipped<- print(paste(i))
			}
			skip<- cbind(skip, skipped)
			}

		skip<- t(skip)
		rownames(tpvalues)= skip[,1]
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finalesophgenes))
		{
		if(! all(is.na(finalesophgenes[i]))){
		valuepc0 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1])
		valuepc1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,2])
		valuepc2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,3])
		valuein1 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,4])
		valuein2 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,5])
		valuein3 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,6])
		valuein4 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,7])
		valuein5 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,8])
		valuein6 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,9])
		valuein7 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,10])
		valuein8 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,11])
		valuein9 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,12])
		valuein10 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,13])
		valuein11 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,14])
		valuein12 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,15])
		valuein13 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,16])
		valuein14 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,17])
		valuein15 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,18])
		valuein16 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,19])
		valuein17 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,20])
		valuein18 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,21])
		valuein19 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,22])
		valuein20 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,23])
		valuein21 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,24])
		valuein22 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,25])
		valuein23 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,26])
		valuein24 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,27])
		valuein25 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,28])
		valuein26 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,29])
		valuein27 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,30])
		valuein28 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,31])
		valuein29 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,32])
		valuein30 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,33])
		valuein31 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,34])
		valuein32 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,35])
		valuein33 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,36])
		valuein34 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,37])
		valuein35 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,38])
		valuein36 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,39])
		valuein37 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,40])
		valuein38 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,41])
		valuein39 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,42])
		valuein40 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,43])
		valuein41 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,44])
		valuein42 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,45])
		valuein43 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,46])
		valuein44 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,47])
		valuein45 <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,48])
		valueall <-lm(finalesophL1HS[,1] ~ finalesophgenes[,c(i)] + finalesophcov[,1] + finalesophcov[,2] + 
					finalesophcov[,3] + finalesophcov[,4] + finalesophcov[,5] + finalesophcov[,6] + finalesophcov[,7] + 
					finalesophcov[,8] + finalesophcov[,9] + finalesophcov[,10] + finalesophcov[,11] + finalesophcov[,12] +
					finalesophcov[,13] + finalesophcov[,14] + finalesophcov[,15] + finalesophcov[,16] + finalesophcov[,17] + finalesophcov[,18] +
					finalesophcov[,19] + finalesophcov[,20] + finalesophcov[,21] + finalesophcov[,22] + finalesophcov[,23] + finalesophcov[,24] + 
					finalesophcov[,25] + finalesophcov[,26] + finalesophcov[,27] + finalesophcov[,28] + finalesophcov[,29] + finalesophcov[,30] +
					finalesophcov[,31] + finalesophcov[,32] + finalesophcov[,33] + finalesophcov[,34] + finalesophcov[,35] + finalesophcov[,36] +
					finalesophcov[,37] + finalesophcov[,38] + finalesophcov[,39] + finalesophcov[,40] + finalesophcov[,41] + finalesophcov[,42] +
					finalesophcov[,43] + finalesophcov[,44] + finalesophcov[,45] + finalesophcov[,46] + finalesophcov[,47] + finalesophcov[,48])
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
		
		
		for (i in colnames(finalesophgenes))
		{
			skipped=NULL
			if(!all(is.na(finalesophgenes[i]))){
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
		setwd("/home/quints1/Correlations/esophagusmuscularis")
		write.table(finalresults, "P and Q esophmuscularis.txt", sep="\t", quote=FALSE)
		
		rvalueresults=NULL	
		for (i in colnames(finalesophgenes))
		{
		if(! all(is.na(finalesophgenes[i]))){
			pearson<- cor.test(finalesophgenes[,c(i)], finalesophL1HS[,1], method="pearson")$estimate
			spearman<- cor.test(finalesophgenes[,c(i)], finalesophL1HS[,1], method="spearman")$estimate
			combinedr<- rbind(pearson, spearman)}
			else{
				print(paste("Skipping", i))
			}
			rvalueresults<- cbind(rvalueresults, combinedr)
			}
		rvalueresults<- data.frame(rvalueresults)
		trvalueresults<- t(rvalueresults)
		
		skip=NULL	
		for (i in colnames(finalesophgenes))
		{
			skipped=NULL
			if(!all(is.na(finalesophgenes[i]))){
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
		setwd("/home/quints1/Correlations/esophagusmuscularis")
		write.table(finalresults, "PRandQEsophMuscularis.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
		