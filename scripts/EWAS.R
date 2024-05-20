
rm(list=ls())
load("DRN_ML_ReadyforEWAS.Rdata")

# Load necessary libraries
library(sva)

# Perform SVA on U, mC, and hmC datasets
# Ensure the 'design' matrix is properly defined in your environment

# SVA on U dataset
mydf <- na.omit(U.f)
sva.U <- sva(dat = as.matrix(mydf), mod = design, mod0 = design[ , -2])

# SVA on mC dataset
mydf <- na.omit(mC.f)
sva.mC <- sva(dat = as.matrix(mydf), mod = design, mod0 = design[ , -2])

# Filter hmC dataset for relevant CpG sites
hmC.f <- hmC.f[rownames(hmC.f) %in% cpg.hmc, ]

# SVA on hmC dataset
mydf <- na.omit(hmC.f)
sva.hmC <- sva(dat = as.matrix(mydf), mod = design, mod0 = design[ , -2])

# Save the SVA results
save(sva.U, sva.mC, sva.hmC, file = "sva.drn.83.Rdata")
###############################################################
# Load necessary libraries
library(sva)
library(parallel)

# Define the EWAS function
EWAS <- function(x, BraakStage, Ageatdeath, Gender, Chipid, sv1, sv2, sv3, sv4, sv5) {
  fit <- lm(BraakStage ~ as.numeric(x) + Ageatdeath + Gender + Chipid + sv1 + sv2 + sv3 + sv4 + sv5)
  return(cbind(coef(summary(fit))[2, ], shapiro.test(fit$residuals)$p.value))
}

# Set working directory and load data
setwd("/mnt/data1/Ehsan/DRN/")
load("sva.drn.83.Rdata")

# Convert phenotype variables to appropriate types
pheno$BraakStage <- as.numeric(pheno$BraakStage)
pheno$Gender <- as.factor(pheno$Gender)
pheno$Ageatdeath <- as.numeric(pheno$Ageatdeath)
Chipid <- as.factor(substr(pheno$sentrix_ID, 1, 12))
BraakStage <- pheno$BraakStage
Ageatdeath <- pheno$Ageatdeath
Gender <- pheno$Gender

# Perform EWAS for unmodified cytosine (U)
sv1 <- sva.U$sv[, 1]
sv2 <- sva.U$sv[, 2]
sv3 <- sva.U$sv[, 3]
sv4 <- sva.U$sv[, 4]
sv5 <- sva.U$sv[, 5]

cl <- makeCluster(32)
res.U <- t(parApply(cl, U.f, 1, EWAS, BraakStage, Ageatdeath, Gender, Chipid, sv1, sv2, sv3, sv4, sv5))
res.U <- res.U[, c(1:5)]
p.BH.U <- as.data.frame(p.adjust(res.U[, 4], "BH"))
res.U <- cbind(res.U, p.BH.U)
colnames(res.U) <- c("U.Estimate", "U.SE", "U.t", "U.P", "resid.norm.U", "BH.U")
head(res.U)

# Define the EWAS function for mC and hmC
EWAS <- function(x, BraakStage, Ageatdeath, Gender, Chipid, sv1, sv2, sv3) {
  fit <- lm(BraakStage ~ as.numeric(x) + Ageatdeath + Gender + Chipid + sv1 + sv2 + sv3)
  return(cbind(coef(summary(fit))[2, ], shapiro.test(fit$residuals)$p.value))
}

# Perform EWAS for methylated cytosine (mC)
sv1 <- sva.mC$sv[, 1]
sv2 <- sva.mC$sv[, 2]
sv3 <- sva.mC$sv[, 3]
sv4 <- sva.mC$sv[, 4]
sv5 <- sva.mC$sv[, 5]

cl <- makeCluster(64)
res.mC <- t(parApply(cl, mC.f, 1, EWAS, BraakStage, Ageatdeath, Gender, Chipid, sv1, sv2, sv3))
res.mC <- res.mC[, c(1:5)]
p.BH.mC <- as.data.frame(p.adjust(res.mC[, 4], "BH"))
res.mC <- cbind(res.mC, p.BH.mC)
colnames(res.mC) <- c("mC.Estimate", "mC.SE", "mC.t", "mC.P", "resid.norm.mC", "BH.mC")

# Perform EWAS for hydroxymethylated cytosine (hmC)
sv1 <- sva.hmC$sv[, 1]
sv2 <- sva.hmC$sv[, 2]
sv3 <- sva.hmC$sv[, 3]
sv4 <- sva.hmC$sv[, 4]
sv5 <- sva.hmC$sv[, 5]

res.hmC <- t(parApply(cl, hmC.f, 1, EWAS, BraakStage, Ageatdeath, Gender, Chipid, sv1, sv2, sv3))
res.hmC <- res.hmC[, c(1:5)]
colnames(res.hmC) <- c("hmC.Estimate", "hmC.SE", "hmC.t", "hmC.P", "resid.norm.hmC")
stopCluster(cl)

# Combine EWAS results
REULTS.DRN <- cbind(res.hmC, res.mC, res.U)
dim(REULTS.DRN)

# Remove cross-hybridising probes
cross <- read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist <- cross[, 1]
REULTS.DRN.1 <- REULTS.DRN[!rownames(REULTS.DRN) %in% crosslist, ]
dim(REULTS.DRN.1)

# Remove probes with SNPs
snpProbes <- read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
snpProbes <- snpProbes[snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95, ]
SNPlist <- snpProbes[, 1]
REULTS.DRN.2 <- REULTS.DRN.1[!rownames(REULTS.DRN.1) %in% SNPlist, ]

# Remove bad probes
badProbes <- read.csv("/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE, header = TRUE)
badlist <- badProbes[, 1]
REULTS.DRN.3 <- REULTS.DRN.2[!rownames(REULTS.DRN.2) %in% badlist, ]
dim(REULTS.DRN.1)
dim(REULTS.DRN.2)
dim(REULTS.DRN.3)

# Remove probes on sex chromosomes
epicManifest <- read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7)
ANNOT <- data.frame(epicManifest)
sexchrX <- ANNOT[ANNOT$CHR == "X", ]
sexchrY <- ANNOT[ANNOT$CHR == "Y", ]

REULTS.DRN.4 <- REULTS.DRN.3[!rownames(REULTS.DRN.3) %in% sexchrX$IlmnID, ]
REULTS.DRN.5 <- REULTS.DRN.4[!rownames(REULTS.DRN.4) %in% sexchrY$IlmnID, ]

ANNOT.sub <- ANNOT[ANNOT$Name %in% rownames(REULTS.DRN.5), ]

ANNOT.sub <- ANNOT.sub[order(ANNOT.sub$Name), ]
REULTS.DRN.5 <- REULTS.DRN.5[order(rownames(REULTS.DRN.5)), ]

# Check if row names match
identical(as.character(rownames(REULTS.DRN.5)), as.character(ANNOT.sub$Name))

# Combine results with annotations
DRN.RESULTS.ANNOT <- cbind(REULTS.DRN.5, ANNOT.sub)
head(DRN.RESULTS.ANNOT)
table(DRN.RESULTS.ANNOT$CHR)

# Save final annotated results
save(DRN.RESULTS.ANNOT, file = "DRN.RESULTS.ANNOT.Rdata")

