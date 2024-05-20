
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
  cbind(coef(summary(fit))[2, ], shapiro.test(fit$residuals)$p.value)
}

# Set working directory and load data
setwd("/mnt/data1/Ehsan/DRN/")
load("sva.drn.83.Rdata")

# Convert phenotype variables to appropriate types
pheno <- transform(pheno,
                   BraakStage = as.numeric(BraakStage),
                   Gender = as.factor(Gender),
                   Ageatdeath = as.numeric(Ageatdeath),
                   Chipid = as.factor(substr(sentrix_ID, 1, 12)))

# Define SVA variables
sv_vars <- function(sva_obj) list(sva_obj$sv[, 1], sva_obj$sv[, 2], sva_obj$sv[, 3], sva_obj$sv[, 4], sva_obj$sv[, 5])

# Perform EWAS analysis with parallel processing
perform_ewas <- function(data, sva_obj) {
  sv_list <- sv_vars(sva_obj)
  cl <- makeCluster(detectCores())
  res <- t(parApply(cl, data, 1, EWAS, pheno$BraakStage, pheno$Ageatdeath, pheno$Gender, pheno$Chipid, sv_list[[1]], sv_list[[2]], sv_list[[3]], sv_list[[4]], sv_list[[5]]))
  stopCluster(cl)
  res <- res[, c(1:5)]
  p.BH <- p.adjust(res[, 4], "BH")
  res <- cbind(res, p.BH)
  colnames(res) <- c("Estimate", "SE", "t", "P", "resid.norm", "BH")
  res
}

# Run EWAS for U, mC, and hmC
res.U <- perform_ewas(U.f, sva.U)
res.mC <- perform_ewas(mC.f, sva.mC)
res.hmC <- perform_ewas(hmC.f[hmC.f %in% cpg.hmc, ], sva.hmC)

# Combine EWAS results
REULTS.DRN <- cbind(res.hmC, res.mC, res.U)

# Remove cross-hybridising probes, SNPs, and bad probes
remove_probes <- function(results, file_path, column = 1) {
  probes <- read.table(file_path, stringsAsFactors = FALSE)[, column]
  results[!rownames(results) %in% probes, ]
}

REULTS.DRN.1 <- remove_probes(REULTS.DRN, "/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt")
REULTS.DRN.2 <- remove_probes(REULTS.DRN.1, "/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", column = "EUR_AF >= 0.05 & EUR_AF <= 0.95")
REULTS.DRN.3 <- remove_probes(REULTS.DRN.2, "/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv")

# Remove probes on sex chromosomes
epicManifest <- read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7)
ANNOT <- data.frame(epicManifest)
sexchrX <- ANNOT$IlmnID[ANNOT$CHR == "X"]
sexchrY <- ANNOT$IlmnID[ANNOT$CHR == "Y"]

REULTS.DRN.4 <- REULTS.DRN.3[!rownames(REULTS.DRN.3) %in% sexchrX, ]
REULTS.DRN.5 <- REULTS.DRN.4[!rownames(REULTS.DRN.4) %in% sexchrY, ]

# Combine results with annotations
ANNOT.sub <- ANNOT[ANNOT$Name %in% rownames(REULTS.DRN.5), ]
ANNOT.sub <- ANNOT.sub[order(ANNOT.sub$Name), ]
REULTS.DRN.5 <- REULTS.DRN.5[order(rownames(REULTS.DRN.5)), ]
identical(as.character(rownames(REULTS.DRN.5)), as.character(ANNOT.sub$Name))

DRN.RESULTS.ANNOT <- cbind(REULTS.DRN.5, ANNOT.sub)
save(DRN.RESULTS.ANNOT, file = "DRN.RESULTS.ANNOT.Rdata")
