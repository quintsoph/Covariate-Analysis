##STEP 1: organize the covariates and make sure all the patients line up with each other
	#Make all the L1HS patients line up with the genes
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

	#adrenal
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	adrenalpats<- read.table("Adrenaledits.txt", header=TRUE, sep="\t")
	rownames(adrenalpats)=str_replace_all(adrenalpats[,1], "-", ".")
	adrenalcut<- merge(adrenalpats, finalpatlist, by.x="row.names", by.y="row.names")
	finaladrenalpat<- adrenalcut[,5]
	cutnames<- data.frame(stri_sub(finaladrenalpat, 0, -15))
	adrenalnames<- cbind(finaladrenalpat, cutnames)
	rownames(adrenalnames)=adrenalnames[,2]

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	adrenalcovs<- read.table("Adrenal.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tadrenalcovs<- data.frame(t(adrenalcovs))
	covnames<- data.frame(stri_sub(rownames(tadrenalcovs), 0, -15))
	adrenalcovs2<- tadrenalcovs
	adrenalcovs2<- cbind(covnames, adrenalcovs2)
	rownames(adrenalcovs2) = adrenalcovs2[,1]
	adrenalresults <- merge(adrenalcovs2, adrenalnames, by.x="row.names", by.y="row.names", drop=FALSE)
	adrenalcovfinal<- adrenalresults[,c(21, 3:20)]
	setwd("/home/quints1/Correlations/adrenal")
	write.table(adrenalcovfinal, "adrenaleditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Bladder is missing

	#Breast
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	breastpats<- read.table("Breastedits.txt", header=TRUE, sep="\t")
	rownames(breastpats)=str_replace_all(breastpats[,1], "-", ".")
	breastcut<- merge(breastpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalbreastpat<- breastcut[,5]
	cutnames<- data.frame(stri_sub(finalbreastpat, 0, -15))
	breastnames<- data.frame(cbind(finalbreastpat, cutnames))
	colnames(breastnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	breastcovs<- read.table("Breast.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tbreastcovs<- data.frame(t(breastcovs))
	covnames<- data.frame(stri_sub(rownames(tbreastcovs), 0, -15))
	breastcovs2<- tbreastcovs
	breastcovs2<- cbind(covnames, breastcovs2)
	rownames(breastcovs2) = breastcovs2[,1]
	breastresults <- merge(breastcovs2, breastnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	breastcovfinal<- breastresults[,c(36, 3:35)]
	setwd("/home/quints1/Correlations/breast")
	write.table(breastcovfinal, "breasteditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Colon - Sigmoid
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	colonpats<- read.table("Colonedits.txt", header=TRUE, sep="\t")
	colonpats<- subset(colonpats, colonpats[,3] == "Colon - Sigmoid", select=c(1,2,3))
	rownames(colonpats)=str_replace_all(colonpats[,1], "-", ".")
	coloncut<- merge(colonpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalcolonpat<- coloncut[,5]
	cutnames<- data.frame(stri_sub(finalcolonpat, 0, -15))
	colonnames<- data.frame(cbind(finalcolonpat, cutnames))
	colnames(colonnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	coloncovs<- read.table("ColonSigmoid.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tcoloncovs<- data.frame(t(coloncovs))
	covnames<- data.frame(stri_sub(rownames(tcoloncovs), 0, -15))
	coloncovs2<- tcoloncovs
	coloncovs2<- cbind(covnames, coloncovs2)
	rownames(coloncovs2) = coloncovs2[,1]
	colonresults <- merge(coloncovs2, colonnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	coloncovfinal<- colonresults[,c(36, 3:35)]
	setwd("/home/quints1/Correlations/colonsig")
	write.table(coloncovfinal, "colonsigeditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Colon - Transverse
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	colonpats<- read.table("Colonedits.txt", header=TRUE, sep="\t")
	colonpats<- subset(colonpats, colonpats[,3] == "Colon - Transverse", select=c(1,2,3))
	rownames(colonpats)=str_replace_all(colonpats[,1], "-", ".")
	coloncut<- merge(colonpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalcolonpat<- coloncut[,5]
	cutnames<- data.frame(stri_sub(finalcolonpat, 0, -15))
	colonnames<- data.frame(cbind(finalcolonpat, cutnames))
	colnames(colonnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	coloncovs<- read.table("ColonTransverse.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tcoloncovs<- data.frame(t(coloncovs))
	covnames<- data.frame(stri_sub(rownames(tcoloncovs), 0, -15))
	coloncovs2<- tcoloncovs
	coloncovs2<- cbind(covnames, coloncovs2)
	rownames(coloncovs2) = coloncovs2[,1]
	colonresults <- merge(coloncovs2, colonnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	coloncovfinal<- colonresults[,c(34, 3:33)]
	setwd("/home/quints1/Correlations/colontrans")
	write.table(coloncovfinal, "colontranseditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Esophagus - Gastro
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Gastroesophageal Junction", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophcut<- merge(esophpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalesophpat<- esophcut[,5]
	cutnames<- data.frame(stri_sub(finalesophpat, 0, -15))
	esophnames<- data.frame(cbind(finalesophpat, cutnames))
	colnames(esophnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	esophcovs<- read.table("EsophagusGastro.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tesophcovs<- data.frame(t(esophcovs))
	covnames<- data.frame(stri_sub(rownames(tesophcovs), 0, -15))
	esophcovs2<- tesophcovs
	esophcovs2<- cbind(covnames, esophcovs2)
	rownames(esophcovs2) = esophcovs2[,1]
	esophresults <- merge(esophcovs2, esophnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	esophcovfinal<- esophresults[,c(21, 3:20)]
	setwd("/home/quints1/Correlations/esophagusgastro")
	write.table(esophcovfinal, "esophagusgastroeditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Esophagus - Mucosa
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Mucosa", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophcut<- merge(esophpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalesophpat<- esophcut[,5]
	cutnames<- data.frame(stri_sub(finalesophpat, 0, -15))
	esophnames<- data.frame(cbind(finalesophpat, cutnames))
	colnames(esophnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	esophcovs<- read.table("EsophagusMucosa.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tesophcovs<- data.frame(t(esophcovs))
	covnames<- data.frame(stri_sub(rownames(tesophcovs), 0, -15))
	esophcovs2<- tesophcovs
	esophcovs2<- cbind(covnames, esophcovs2)
	rownames(esophcovs2) = esophcovs2[,1]
	esophresults <- merge(esophcovs2, esophnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	esophcovfinal<- esophresults[,c(51, 3:50)]
	setwd("/home/quints1/Correlations/esophagusmucosa")
	write.table(esophcovfinal, "esophagusmucosaeditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Esophagus - Muscularis
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	esophpats<- read.table("Esophagusedits.txt", header=TRUE, sep="\t")
	esophpats<- subset(esophpats, esophpats[,3] == "Esophagus - Muscularis", select=c(1,2,3))
	rownames(esophpats)=str_replace_all(esophpats[,1], "-", ".")
	esophcut<- merge(esophpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalesophpat<- esophcut[,5]
	cutnames<- data.frame(stri_sub(finalesophpat, 0, -15))
	esophnames<- data.frame(cbind(finalesophpat, cutnames))
	colnames(esophnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	esophcovs<- read.table("EsophagusMuscularis.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tesophcovs<- data.frame(t(esophcovs))
	covnames<- data.frame(stri_sub(rownames(tesophcovs), 0, -15))
	esophcovs2<- tesophcovs
	esophcovs2<- cbind(covnames, esophcovs2)
	rownames(esophcovs2) = esophcovs2[,1]
	esophresults <- merge(esophcovs2, esophnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	esophcovfinal<- esophresults[,c(51, 3:50)]
	setwd("/home/quints1/Correlations/esophagusmuscularis")
	write.table(esophcovfinal, "esophagusmusculariseditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Kidney
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	kidneypats<- read.table("Kidneyedits.txt", header=TRUE, sep="\t")
	kidneypats<- subset(kidneypats, kidneypats[,3] == "Kidney - Cortex", select=c(1,2,3))
	rownames(kidneypats)=str_replace_all(kidneypats[,1], "-", ".")
	kidneycut<- merge(kidneypats, finalpatlist, by.x="row.names", by.y="row.names")
	finalkidneypat<- kidneycut[,5]
	cutnames<- data.frame(stri_sub(finalkidneypat, 0, -15))
	kidneynames<- data.frame(cbind(finalkidneypat, cutnames))
	colnames(kidneynames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	kidneycovs<- read.table("Kidney.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tkidneycovs<- data.frame(t(kidneycovs))
	covnames<- data.frame(stri_sub(rownames(tkidneycovs), 0, -15))
	kidneycovs2<- tkidneycovs
	kidneycovs2<- cbind(covnames, kidneycovs2)
	rownames(kidneycovs2) = kidneycovs2[,1]
	kidneyresults <- merge(kidneycovs2, kidneynames, by.x="row.names", by.y="cutnames", drop=FALSE)
	kidneycovfinal<- kidneyresults[,c(21, 3:20)]
	setwd("/home/quints1/Correlations/kidney")
	write.table(kidneycovfinal, "kidneyeditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Liver
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	liverpats<- read.table("Liveredits.txt", header=TRUE, sep="\t")
	liverpats<- subset(liverpats, liverpats[,3] == "Liver", select=c(1,2,3))
	rownames(liverpats)=str_replace_all(liverpats[,1], "-", ".")
	livercut<- merge(liverpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalliverpat<- livercut[,5]
	cutnames<- data.frame(stri_sub(finalliverpat, 0, -15))
	livernames<- data.frame(cbind(finalliverpat, cutnames))
	colnames(livernames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	livercovs<- read.table("Liver.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tlivercovs<- data.frame(t(livercovs))
	covnames<- data.frame(stri_sub(rownames(tlivercovs), 0, -15))
	livercovs2<- tlivercovs
	livercovs2<- cbind(covnames, livercovs2)
	rownames(livercovs2) = livercovs2[,1]
	liverresults <- merge(livercovs2, livernames, by.x="row.names", by.y="cutnames", drop=FALSE)
	livercovfinal<- liverresults[,c(21, 3:20)]
	setwd("/home/quints1/Correlations/liver")
	write.table(livercovfinal, "livereditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Pancreas
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	pancreaspats<- read.table("Pancreasedits.txt", header=TRUE, sep="\t")
	pancreaspats<- subset(pancreaspats, pancreaspats[,3] == "Pancreas", select=c(1,2,3))
	rownames(pancreaspats)=str_replace_all(pancreaspats[,1], "-", ".")
	pancreascut<- merge(pancreaspats, finalpatlist, by.x="row.names", by.y="row.names")
	finalpancreaspat<- pancreascut[,5]
	cutnames<- data.frame(stri_sub(finalpancreaspat, 0, -15))
	pancreasnames<- data.frame(cbind(finalpancreaspat, cutnames))
	colnames(pancreasnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	pancreascovs<- read.table("Pancreas.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tpancreascovs<- data.frame(t(pancreascovs))
	covnames<- data.frame(stri_sub(rownames(tpancreascovs), 0, -15))
	pancreascovs2<- tpancreascovs
	pancreascovs2<- cbind(covnames, pancreascovs2)
	rownames(pancreascovs2) = pancreascovs2[,1]
	pancreasresults <- merge(pancreascovs2, pancreasnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	pancreascovfinal<- pancreasresults[,c(36, 3:35)]
	setwd("/home/quints1/Correlations/pancreas")
	write.table(pancreascovfinal, "pancreaseditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Prostate
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	prostatepats<- read.table("Prostateedits.txt", header=TRUE, sep="\t")
	prostatepats<- subset(prostatepats, prostatepats[,3] == "Prostate", select=c(1,2,3))
	prostatepats2<- cbind(prostatepats, (str_replace_all(prostatepats[,1], "-", ".")))
	colnames(prostatepats2)= c("names", "expression", "type", "newnames")
	prostatecut<- merge(prostatepats2, finalpatlist, by.x="newnames", by.y="row.names")
	finalprostatepat<- prostatecut[,5]
	cutnames<- data.frame(stri_sub(finalprostatepat, 0, -15))
	prostatenames<- data.frame(cbind(finalprostatepat, cutnames))
	colnames(prostatenames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	prostatecovs<- read.table("Prostate.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tprostatecovs<- data.frame(t(prostatecovs))
	covnames<- data.frame(stri_sub(rownames(tprostatecovs), 0, -15))
	prostatecovs2<- tprostatecovs
	prostatecovs2<- cbind(covnames, prostatecovs2)
	rownames(prostatecovs2) = prostatecovs2[,1]
	prostateresults <- merge(prostatecovs2, prostatenames, by.x="row.names", by.y="cutnames", drop=FALSE)
	prostatecovfinal<- prostateresults[,c(21, 3:20)]
	setwd("/home/quints1/Correlations/prostate")
	write.table(prostatecovfinal, "prostateeditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Stomach
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	stomachpats<- read.table("Stomachedits.txt", header=TRUE, sep="\t")
	stomachpats<- subset(stomachpats, stomachpats[,3] == "Stomach", select=c(1,2,3))
	rownames(stomachpats)=str_replace_all(stomachpats[,1], "-", ".")
	stomachcut<- merge(stomachpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalstomachpat<- stomachcut[,5]
	cutnames<- data.frame(stri_sub(finalstomachpat, 0, -15))
	stomachnames<- data.frame(cbind(finalstomachpat, cutnames))
	colnames(stomachnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	stomachcovs<- read.table("Stomach.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tstomachcovs<- data.frame(t(stomachcovs))
	covnames<- data.frame(stri_sub(rownames(tstomachcovs), 0, -15))
	stomachcovs2<- tstomachcovs
	stomachcovs2<- cbind(covnames, stomachcovs2)
	rownames(stomachcovs2) = stomachcovs2[,1]
	stomachresults <- merge(stomachcovs2, stomachnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	stomachcovfinal<- stomachresults[,c(36, 3:35)]
	setwd("/home/quints1/Correlations/stomach")
	write.table(stomachcovfinal, "stomacheditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

	#Thyroid
	setwd("/home/quints1/Correlations/L1HSpatientedits")
	thyroidpats<- read.table("Thyroidedits.txt", header=TRUE, sep="\t")
	thyroidpats<- subset(thyroidpats, thyroidpats[,3] == "Thyroid", select=c(1,2,3))
	rownames(thyroidpats)=str_replace_all(thyroidpats[,1], "-", ".")
	thyroidcut<- merge(thyroidpats, finalpatlist, by.x="row.names", by.y="row.names")
	finalthyroidpat<- thyroidcut[,5]
	cutnames<- data.frame(stri_sub(finalthyroidpat, 0, -15))
	thyroidnames<- data.frame(cbind(finalthyroidpat, cutnames))
	colnames(thyroidnames)=c("realname", "cutnames")

	setwd("/home/quints1/Deseqstuff/covariatefiles")
	thyroidcovs<- read.table("Thyroid.combined_covariates.txt", header=TRUE, sep="\t", row.names=1)
	tthyroidcovs<- data.frame(t(thyroidcovs))
	covnames<- data.frame(stri_sub(rownames(tthyroidcovs), 0, -15))
	thyroidcovs2<- tthyroidcovs
	thyroidcovs2<- cbind(covnames, thyroidcovs2)
	rownames(thyroidcovs2) = thyroidcovs2[,1]
	thyroidresults <- merge(thyroidcovs2, thyroidnames, by.x="row.names", by.y="cutnames", drop=FALSE)
	thyroidcovfinal<- thyroidresults[,c(51, 3:50)]
	setwd("/home/quints1/Correlations/thyroid")
	write.table(thyroidcovfinal, "thyroideditcov.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

###STEP 2: Start on doing the models to find best result
	#read in the L1HS vst and cntvst
	setwd("/home/quints1/Deseqstuff")
	L1HSraw<- read.table("L1HS.VST.cnts.txt", header=TRUE, row.names=1)
	genesraw<- read.table("VSTcnt.txt", header=TRUE, row.names=1)
	tgenesraw<- t(genesraw)
	library("stringi")
	library("stringr")
	
	#adrenal
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
	finaladrenalcov<- data.frame(adrenalall[,c(56206:56223)])	
	rownames(finaladrenalcov)=adrenalall[,1]
	finaladrenalL1HS<- data.frame(adrenalall[,c(56205)])
	rownames(finaladrenalL1HS)= adrenalall[,1]
	colnames(finaladrenalL1HS)= c("L1HS")
	finaladrenalgenes<- data.frame(adrenalall[,c(3:56204)])
	rownames(finaladrenalgenes)=adrenalall[,1]
	library(nlme)
	library(MuMIn)
	library(lmSupport)
	library(qvalue)
	library(matrixStats)
		#linear models - pc0
		testpvaluespc0<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,1])
			anova(value)$Pr[2]
		}
		
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluespc0(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovpc0<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear models - pc1
		testpvaluespc1<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,2])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluespc1(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovpc1<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		# linear models - PC2
		testpvaluespc2<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,3])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluespc2(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovpc2<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov1
		testpvaluesin1<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,4])
			anova(value)$Pr[2]
		}
		
		
		in1value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,4])
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin1(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin1<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov2
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,5])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin2<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov3
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,6])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin3<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov4
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,7])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin4<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov5
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,8])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin5<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov6
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,9])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin6<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov7
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,10])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin7<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov8
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,11])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin8<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov9
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,12])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin9<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov10
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,13])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin10<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov11
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,14])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin11<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov12
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,15])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin12<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov13
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,16])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin13<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov14
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,17])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin14<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#linear model - InferredCov15
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,18])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalcovin15<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		#combined covariates
		testpvaluesin<- function(x)
		{
			value <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(x)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18])
			anova(value)$Pr[2]
		}
		
		pvalues=NULL
		for (i in colnames(finaladrenalgenes))
		{
			p = NULL
			#If not entire row is NAs, run.
			if(! all(is.na(finaladrenalgenes[i]))){
				p <-testpvaluesin(i)
			}
			else{
				print(paste("Skipping", i))
			}

			pvalues <- cbind(pvalues, p) 
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
		pvalallcov<- tpvalues[order(tpvalues[,c(1)]),, drop=FALSE]
		
	
		#library("dplyr")
		#pc0model= NULL
		#for (i in colnames(finaladrenalgenes))
		#{
		#model=NULL
        #model<- lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1])
		#pc0model<- list(pc0model, model)
		#}
		
		finalAICc=NULL
		for (i in colnames(finaladrenalgenes))
		{
		valuepc0 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1])
		valuepc1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,2])
		valuepc2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,3])
		valuein1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,4])
		valuein2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,5])
		valuein3 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,6])
		valuein4 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,7])
		valuein5 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,8])
		valuein6 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,9])
		valuein7 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,10])
		valuein8 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,11])
		valuein9 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,12])
		valuein10 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,13])
		valuein11 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,14])
		valuein12 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,15])
		valuein13 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,16])
		valuein14 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,17])
		valuein15 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,18])
		valueall <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18])
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
		if(! all(is.na(finaladrenalgenes[i]))){
		valuepc0 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1])
		valuepc1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,2])
		valuepc2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,3])
		valuein1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,4])
		valuein2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,5])
		valuein3 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,6])
		valuein4 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,7])
		valuein5 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,8])
		valuein6 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,9])
		valuein7 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,10])
		valuein8 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,11])
		valuein9 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,12])
		valuein10 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,13])
		valuein11 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,14])
		valuein12 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,15])
		valuein13 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,16])
		valuein14 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,17])
		valuein15 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,18])
		valueall <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18])
		AICc<- model.sel(valuepc0, valuepc1, valuepc2, valuein1, valuein2, valuein3, valuein4, valuein5, valuein6, valuein7, valuein8,
			valuein9, valuein10, valuein11, valuein12, valuein13, valuein14, valuein15,valueall)
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
	write.table(finalAICc, "All AICc results.txt", sep="\t", quote=FALSE)
	
	#rvalue
	rvalues=NULL
	for (i in colnames(finaladrenalgenes))
		{
		if(! all(is.na(finaladrenalgenes[i]))){
		valuepc0 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1])
		valuepc1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,2])
		valuepc2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,3])
		valuein1 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,4])
		valuein2 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,5])
		valuein3 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,6])
		valuein4 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,7])
		valuein5 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,8])
		valuein6 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,9])
		valuein7 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,10])
		valuein8 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,11])
		valuein9 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,12])
		valuein10 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,13])
		valuein11 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,14])
		valuein12 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,15])
		valuein13 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,16])
		valuein14 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,17])
		valuein15 <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,18])
		valueall <-lm(finaladrenalL1HS[,1] ~ finaladrenalgenes[,c(i)] + finaladrenalcov[,1] + finaladrenalcov[,2] + 
					finaladrenalcov[,3] + finaladrenalcov[,4] + finaladrenalcov[,5] + finaladrenalcov[,6] + finaladrenalcov[,7] + 
					finaladrenalcov[,8] + finaladrenalcov[,9] + finaladrenalcov[,10] + finaladrenalcov[,11] + finaladrenalcov[,12] +
					finaladrenalcov[,13] + finaladrenalcov[,14] + finaladrenalcov[,15] + finaladrenalcov[,16] + finaladrenalcov[,17] + finaladrenalcov[,18])
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
		colnames(finalresult)=c("pvalues", "rsqvalues", "qvalues")
		setwd("/home/quints1/Correlations/adrenal")
		write.table(finalresults, "P and Q Adrenal.txt", sep="\t", quote=FALSE)
		
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
		results<- results[,2:4]
		colnames(results)=c("pvalues", "pearsonr", "spearmanr")
		finalresults<- merge(results, qvalues, by.x="row.names", by.y="row.names")
		colnames(finalresults)=c("rownames", "pvalues", "pearsonr", "spearmanr", "qvalues")
		setwd("/home/quints1/Correlations/adrenal")
		write.table(finalresults, "PRandQAdrenal.txt", sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE)
			
			
			
			