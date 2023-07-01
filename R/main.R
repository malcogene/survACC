dir.create("data")
# Benchmark dataset for ACC
loadUrl(
 "https://konkukackr-my.sharepoint.com/:u:/g/personal/palelamp_kku_ac_kr/EfFPEElurM5Ptr38fKMBZ4cB2efBy5Zs4z0pgFlew3fDdQ?e=AXBe72"
)
load(file.path(getwd(),
 "../data/bdb.ACC.rda"))
# Pathway processed
loadUrl(
 "https://konkukackr-my.sharepoint.com/:u:/g/personal/palelamp_kku_ac_kr/EQht_pXT7rhOnyjt4_8OXTEB_dFGAWyeynzr3MoZNxccGQ?e=b7Qh83"
)
load(file.path(getwd(),
 "../data/p.KEGG.PID.BioCarta.rda"))
pw.in <- p.KEGG.PID.BioCarta$entrez_gene_ids
names(pw.in) <- p.KEGG.PID.BioCarta$pathway

# ACC batch corrected
loadUrl(
 "https://konkukackr-my.sharepoint.com/:u:/g/personal/palelamp_kku_ac_kr/Eb3wdN_SLfVJt2luFluHcbQBFHpJcGgMdiY--5qRdAmnpQ?e=kgZNDm"
)
load(file.path(getwd(),
 "../data/exp.BC.rda"))

pathwayMapping = F
RFS = F
EarlyStage = F

if (!RFS) {
 pheno.BC = rbind(bdb.ACC[[1]]$y[, 1:5], bdb.ACC[[2]]$y[, 1:5], bdb.ACC[[3]]$y[, 1:5])
 KMplot.ylab = "Overall Survival (OS)"
} else {
 pheno.BC = rbind(bdb.ACC[[1]]$y[, c(8:9, 3:5)], bdb.ACC[[2]]$y[, c(7:8, 3:5)])
 KMplot.ylab = "Recurrence Free Survival (RFS)"
}



# Early stage ACC
if (EarlyStage) {
 stageIandII.inx = grep("Stage I$|Stage II$", pheno.BC$pathologic_stage)
 pheno.BC <- pheno.BC.stageIandII <- pheno.BC[stageIandII.inx,]
 exp.BC <- exp.BC.stageIandII <- exp.BC[, stageIandII.inx]
 KMplot.ylab = "Early Stage Overall Survival (OS)"
}

# Data Cleaning
colnames(pheno.BC)
colnames(pheno.BC)[1:2] <- c("time", "status")
pheno.BC$pathologic_stage <-
 gsub("NA", "", pheno.BC$pathologic_stage)
pheno.BC$pathologic_stage[c(4, 24, 92)] <- NA
pheno.BC.clean <- na.omit(pheno.BC)
exp.BC.clean <-  exp.BC[, rownames(pheno.BC.clean)]

pheno.BC.clean$gender <- as.factor(pheno.BC.clean$gender)
pheno.BC.clean$pathologic_stage <-
 as.factor(pheno.BC.clean$pathologic_stage)

pheno.BC.clean.numeric <- pheno.BC.clean
pheno.BC.clean.numeric$gender <-
 as.numeric(as.factor(pheno.BC.clean$gender))
pheno.BC.clean.numeric$pathologic_stage <-
 as.numeric(as.factor(pheno.BC.clean$pathologic_stage))

# Pathway Mapping
if (pathwayMapping) {
 exp.BC.clean <- exp.BC.clean.pre
 gsesExp = SGSES(exp.BC.clean, predefined.pw.in)
 length(predefined.pw.in)
 exp.BC.clean <- gsesExp
}

# Model building using development cohort.
exp.dev <- exp.BC.clean[, 1:76]
exp.dev <- na.omit(exp.dev)
dim(exp.dev)

vv <-
 add.gs(exp.dev, pheno.BC.clean.numeric[colnames(exp.dev), ], , "lassoCoxRiskScore")
vv$riskScoreFormula$betaS[vv$riskScoreFormula$gm] %*% exp.dev[vv$riskScoreFormula$gm,]
tmp.out.dev <-
 vv$riskScoreFormula$betaS[vv$riskScoreFormula$gm] %*% exp.dev[vv$riskScoreFormula$gm,]
makeVec(vv$riskScoreFormula$betaS[vv$riskScoreFormula$gm])
makeVec(vv$riskScoreFormula$gm)
tmp.pheno <- pheno.BC.clean[colnames(exp.dev),]
tmp.pheno$rs <- factorizeQuantile(tmp.out.dev)
colnames(tmp.pheno)[3:6] <-
 c("Age", "Gender", "Pathologic Stage", "SurvACC")
levels(tmp.pheno$SurvACC) <- list('_Row Risk' = 1, 'High Risk' = 2)

coxph(Surv(time, status) ~ ., tmp.pheno)
ggsurvplot(survfit(Surv(time, status) ~ SurvACC, tmp.pheno))


# External validation conducted utilizing cohort #1.
exp.val1 <- exp.BC.clean[, 77:119]
tmp.out.val1 <-
 vv$riskScoreFormula$betaS[vv$riskScoreFormula$gm] %*% exp.val1[vv$riskScoreFormula$gm,]
tmp.pheno <- pheno.BC.clean[colnames(exp.val1),]
tmp.pheno$rs <- factorizeQuantile(tmp.out.val1)
colnames(tmp.pheno)[3:6] <-
 c("Age", "Gender", "Pathologic Stage", "SurvACC")
levels(tmp.pheno$SurvACC) <- list('_Row Risk' = 1, 'High Risk' = 2)

coxph(Surv(time, status) ~ ., tmp.pheno)
ggsurvplot(survfit(Surv(time, status) ~ SurvACC, data = tmp.pheno))


# External validation conducted utilizing cohort #2.
exp.val2 <- exp.BC.clean[, 120:142]
tmp.out.val2 <-
 vv$riskScoreFormula$betaS[vv$riskScoreFormula$gm] %*% exp.val2[vv$riskScoreFormula$gm,]
tmp.pheno <- pheno.BC.clean[colnames(exp.val2),]
tmp.pheno$rs <- factorizeQuantile(tmp.out.val2)
colnames(tmp.pheno)[3:6] <-
 c("Age", "Gender", "Pathologic Stage", "SurvACC")
levels(tmp.pheno$SurvACC) <- list('_Row Risk' = 1, 'High Risk' = 2)

coxph(Surv(time, status) ~ ., tmp.pheno)
ggsurvplot(survfit(Surv(time, status) ~ SurvACC, data = tmp.pheno))


# Building reference dataset and the ultimate prediction model.
vv.final <-
 add.gs(exp.BC.clean, pheno.BC.clean.numeric, , "lassoCoxRiskScore")
vv.final$riskScoreFormula$Results
tmp.out.final <-
 vv.final$riskScoreFormula$betaS[vv.final$riskScoreFormula$gm] %*% exp.val2[vv.final$riskScoreFormula$gm,]
vv.final$riskScoreFormula
