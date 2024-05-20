rm(list=ls())

##Loading libraries

library(methylumi)
library(wateRmelon)
require(gdata)
library(minfi)
library(ggplot2)
require(gridExtra)
library(Biobase)
library(MLML2R)
library(sva)
library(qqman)
```

#Loading data

# Set working directory
setwd("/mnt/data1/Ehsan/DRN")

# Load necessary libraries
library(minfi)
library(methylumi)

# Read methylation array data
rgSet <- read.metharray.exp(base = "/mnt/data1/Ehsan/DRN/IDATS/", extended = TRUE)
save(rgSet, file = "rgSet.LC.Rdata")

# Define path to IDAT files
idatPath <- "/mnt/data1/Ehsan/DRN/IDATS/"

# Read EPIC data
drn.methylumi <- readEPIC(idatPath, barcodes = NULL, oob = FALSE, n = TRUE, parallel = TRUE)
save(drn.methylumi, file = "methylumiset.DRN.Rdata")

# Apply p-filter to methylumi object
drn.pfilters <- pfilter(drn.methylumi, perc = 5)
cpg.pfiltered.drn <- as.character(rownames(drn.pfilters))

# Save filtered CpG sites
save(cpg.pfiltered.drn, file = "cpg.pfiltered.lc.Rdata")

# Load phenotype data
DRN_PHENO <- read.csv("pheno.drn.all.csv")
Basename <- substr(DRN_PHENO$sentrix_ID, 14, 32)
DRN_PHENO <- cbind(DRN_PHENO, Basename)

# Identify samples to be filtered
filtered <- which(DRN_PHENO$Basename %in% c("200772260017_R01C01", "200772260017_R02C01", "200793300066_R06C01"))

# Samples recommended by GenomeScan to be dropped
filtered.GS <- which(DRN_PHENO$Basename %in% c("200772260017 R01C01", "200772260017 R02C01", "200772260017 R04C01", "200793300066 R06C01"))

# Display sample names for the filtered indices
filtered_sample_names <- DRN_PHENO$Sample_Name[filtered]

# Display basename for a specific sample name
specific_basename <- DRN_PHENO$Basename[which(DRN_PHENO$Sample_Name == "7058-001-BS-043")]

# Filter out samples from the phenotype data
DRN_PHENO.f <- DRN_PHENO[!DRN_PHENO$Basename %in% c("200772260017_R01C01", "200772260017_R02C01", "200793300066_R06C01", "200793300066_R05C01"),]

# Filter out samples from the rgSet
rgSet.f <- rgSet[, !colnames(rgSet) %in% c("200772260017_R01C01", "200772260017_R02C01", "200793300066_R06C01", "200793300066_R05C01")]

# Split phenotype data into BS and oxBS pools
pheno.bs <- DRN_PHENO.f[DRN_PHENO.f$Pool_ID == "BS",]
pheno.oxbs <- DRN_PHENO.f[DRN_PHENO.f$Pool_ID == "oxBS",]

# Extract sample IDs from Sample_Name
pheno.bs$id <- substr(pheno.bs$Sample_Name, 13, 15)
pheno.oxbs$id <- substr(pheno.oxbs$Sample_Name, 15, 17)

# Check if the IDs are identical
identical_ids <- identical(pheno.bs$id, pheno.oxbs$id)

# Convert Basename columns to character vectors
BS <- as.character(pheno.bs$Basename)
oxBS <- as.character(pheno.oxbs$Basename)

# Subset rgSet by BS and oxBS samples
RGset.BS <- rgSet[, BS]
RGset.oxBS <- rgSet[, oxBS]

# Verify column names match Basename vectors
identical_BS <- identical(colnames(RGset.BS), BS)
identical_oxBS <- identical(colnames(RGset.oxBS), oxBS)

# Remove unnecessary objects from the environment
rm(rgSet, DRN_PHENO, DRN_PHENO.1)

# Save the final datasets
save(RGset.BS, RGset.oxBS, pheno.bs, pheno.oxbs, file = "DRN_RGsets_phenos.Rdata")

# Set working directory
setwd("/mnt/data1/Ehsan/DRN")

# Load necessary libraries
library(minfi)
library(methylumi)

# Preprocess raw data for BS and oxBS
MSet.BS <- preprocessRaw(RGset.BS)
MSet.oxBS <- preprocessRaw(RGset.oxBS)

# Plot quality control metrics
par(mfrow=c(1,2))
plotQC(getQC(MSet.BS))
plotQC(getQC(MSet.oxBS))

# Plot density of BS and oxBS
par(mfrow=c(1,2))
densityPlot(MSet.BS, main = "BS")
densityPlot(MSet.oxBS, main = "oxBS")

# Extract and plot the intensities of the red/green background and modification signals
## BS
con.red.BS <- getRed(RGset.BS)[getProbeInfo(RGset.BS, type = "Control")$Address, ]
rownames(con.red.BS) <- getProbeInfo(RGset.BS, type = "Control")$ExtendedType
con.green.BS <- getGreen(RGset.BS)[getProbeInfo(RGset.BS, type = "Control")$Address, ]
rownames(con.green.BS) <- getProbeInfo(RGset.BS, type = "Control")$ExtendedType

neg.red.BS <- con.red.BS[getProbeInfo(RGset.BS, type = "Control")$Type == "NEGATIVE", ]
neg.green.BS <- con.green.BS[getProbeInfo(RGset.BS, type = "Control")$Type == "NEGATIVE", ]

neg.red.mean.BS <- colMeans(neg.red.BS)
neg.green.mean.BS <- colMeans(neg.green.BS)
green.mean.BS <- colMeans(getGreen(RGset.BS))
red.mean.BS <- colMeans(getRed(RGset.BS))

## OXBS
con.red.oxBS <- getRed(RGset.oxBS)[getProbeInfo(RGset.oxBS, type = "Control")$Address, ]
rownames(con.red.oxBS) <- getProbeInfo(RGset.oxBS, type = "Control")$ExtendedType
con.green.oxBS <- getGreen(RGset.oxBS)[getProbeInfo(RGset.oxBS, type = "Control")$Address, ]
rownames(con.green.oxBS) <- getProbeInfo(RGset.oxBS, type = "Control")$ExtendedType

neg.red.oxBS <- con.red.oxBS[getProbeInfo(RGset.oxBS, type = "Control")$Type == "NEGATIVE", ]
neg.green.oxBS <- con.green.oxBS[getProbeInfo(RGset.oxBS, type = "Control")$Type == "NEGATIVE", ]

neg.red.mean.oxBS <- colMeans(neg.red.oxBS)
neg.green.mean.oxBS <- colMeans(neg.green.oxBS)
green.mean.oxBS <- colMeans(getGreen(RGset.oxBS))
red.mean.oxBS <- colMeans(getRed(RGset.oxBS))

# Plot the intensities
par(mfrow=c(1,4))

# Plot BS background intensities
plot(neg.red.mean.BS, xlab = "Sample", ylab = "Background Intensity", col = "red", pch = 16, ylim = c(0, 5000))
par(new = TRUE)
plot(neg.green.mean.BS, xlab = "Sample", ylab = "", col = "green", pch = 16, ylim = c(0, 5000))

# Plot BS detected intensities
plot(red.mean.BS, xlab = "Sample", ylab = "Detected Intensity", col = "red", pch = 16, ylim = c(0, 5000))
par(new = TRUE)
plot(green.mean.BS, xlab = "Sample", ylab = "", col = "green", pch = 16, ylim = c(0, 5000))

# Plot oxBS background intensities
plot(neg.red.mean.oxBS, xlab = "Sample", ylab = "Background Intensity", col = "red", pch = 16, ylim = c(0, 5000))
par(new = TRUE)
plot(neg.green.mean.oxBS, xlab = "Sample", ylab = "", col = "green", pch = 16, ylim = c(0, 5000))

# Plot oxBS detected intensities
plot(red.mean.oxBS, xlab = "Sample", ylab = "Detected Intensity", col = "red", pch = 16, ylim = c(0, 5000))
par(new = TRUE)
plot(green.mean.oxBS, xlab = "Sample", ylab = "", col = "green", pch = 16, ylim = c(0, 5000))

# Clean up the environment
rm(con.red.BS, con.green.BS, neg.red.BS, neg.green.BS, neg.red.mean.BS, neg.green.mean.BS, green.mean.BS, red.mean.BS,
   con.red.oxBS, con.green.oxBS, neg.red.oxBS, neg.green.oxBS, neg.red.mean.oxBS, neg.green.mean.oxBS, green.mean.oxBS, red.mean.oxBS)


# Normalisation (Noob) - Applying R/G ratio flip to fix dye bias
MSet.Noob.bs.drn  <- preprocessNoob(RGset.BS)
MSet.Noob.oxbs.drn <- preprocessNoob(RGset.oxBS)

# Plot the normalised data
par(mfrow = c(1, 2))
plotQC(getQC(MSet.Noob.bs.drn))
plotQC(getQC(MSet.Noob.oxbs.drn))

# Check if the column names of BS and oxBS are matched
pheno.bs$id <- substr(pheno.bs$Sample_Name, 13, 15)
pheno.oxbs$id <- substr(pheno.oxbs$Sample_Name, 15, 17)
colnames(MSet.Noob.bs.drn) <- pheno.bs$id
colnames(MSet.Noob.oxbs.drn) <- pheno.oxbs$id
print(identical(colnames(MSet.Noob.bs.drn), colnames(MSet.Noob.oxbs.drn)))

# Extracting methylated and unmethylated values from BS & oxBS
methBS <- getMeth(MSet.Noob.bs.drn)
methOxBS <- getMeth(MSet.Noob.oxbs.drn)
unmethBS <- getUnmeth(MSet.Noob.bs.drn)
unmethOxBS <- getUnmeth(MSet.Noob.oxbs.drn)

# Plot the signal intensities after normalisation
par(mfrow = c(1, 2))
plot(colMeans(methBS), xlab = "Sample", ylab = "Signal Intensity", col = "green", pch = 16, ylim = c(0, 7000))
par(new = TRUE)
plot(colMeans(unmethBS), xlab = "Sample", ylab = "", col = "red", pch = 16, ylim = c(0, 7000))

plot(colMeans(methOxBS), xlab = "Sample", ylab = "Signal Intensity", col = "green", pch = 16, ylim = c(0, 7000))
par(new = TRUE)
plot(colMeans(unmethOxBS), xlab = "Sample", ylab = "", col = "red", pch = 16, ylim = c(0, 7000))

# Getting beta values and subtracting oxBS from BS if the column names of normalised betas match for BS and oxBS data
betas.Noob.bs.drn <- getBeta(MSet.Noob.bs.drn)
betas.Noob.oxbs.drn <- getBeta(MSet.Noob.oxbs.drn)
betas.raw.bs.drn <- getBeta(MSet.BS)
betas.raw.oxbs.drn <- getBeta(MSet.oxBS)

print(identical(colnames(betas.Noob.bs.drn), colnames(betas.Noob.oxbs.drn)))

# Subtract oxBS betas from BS betas
mydf <- betas.Noob.bs.drn - betas.Noob.oxbs.drn

# Plotting the BS-oxBS and beta density plots before and after normalisation
par(mfrow = c(1, 2))

# BS - oxBS density plot
plot(density(mydf), main = "BS-oxBS")
abline(v = 0, col = "red", lty = 1, lwd = 1)

# Density plots before normalisation
par(mfrow = c(1, 2))
plot(density(na.omit(betas.raw.bs.drn)), col = "royalblue3", ylim = c(0, 7), main = "Raw", lwd = 2, xlab = "")
par(new = TRUE)
plot(density(na.omit(betas.raw.oxbs.drn)), col = "red4", ylim = c(0, 7), main = "Raw", lwd = 2, xlab = "")

# Density plots after normalisation
plot(density(na.omit(betas.Noob.bs.drn)), col = "royalblue3", ylim = c(0, 7), main = "Internal Control", lwd = 2, xlab = "")
par(new = TRUE)
plot(density(na.omit(betas.Noob.oxbs.drn)), col = "red4", ylim = c(0, 7), main = "Internal Control", lwd = 2, xlab = "")



 plotting BS-oxBS in percentile BS and oxBS beta value bins

# Define a function to get row names based on beta value ranges
get_rownames_by_range <- function(data, start, end, step) {
  lapply(seq(start, end, by=step), function(i) {
    rownames(data[which(rowMeans(data) > i & rowMeans(data) <= (i + step)), ])
  })
}

# Get row names for both oxbs and bs datasets
oxbs_ranges <- get_rownames_by_range(betas.Noob.oxbs.drn, 0, 1, 0.01)
bs_ranges <- get_rownames_by_range(betas.Noob.bs.drn, 0.1, 1, 0.01)

# Function to generate boxplots
generate_boxplots <- function(data, ranges) {
  boxplot(sapply(ranges, function(r) rowMeans(data[rownames(data) %in% r, ])))
}

# Plotting
par(mfrow=c(1, 2))

# Generate boxplots for oxbs and bs datasets
generate_boxplots(mydf, oxbs_ranges)
abline(h=0, col="red", lty=1, lwd=1)
generate_boxplots(mydf, bs_ranges)
abline(h=0, col="red", lty=1, lwd=1)


# Maximum likelihood estimation of states of cytosine
results_exact <- MLML(T.matrix = methBS, U.matrix = unmethBS, L.matrix = unmethOxBS, M.matrix = methOxBS)

# Extract results
U <- data.frame(results_exact$C)
mC <- data.frame(results_exact$mC)
hmC <- data.frame(results_exact$hmC)

# Filter for rows where mean BS betas are greater than 0.1
hmC.f <- hmC[rownames(betas.Noob.bs.drn[which(rowMeans(betas.Noob.bs.drn) > 0.1), ]), ]

# Identify rows with fewer than half the columns equal to zero
a <- data.frame(apply(hmC.f, 1, function(x) sum(x == 0)))
a <- which(a[, 1] < 0.5 * ncol(hmC.f))

# Subset CpG sites with non-zero hmC values
length(a) #[1] 806673
cpg.hmc <- rownames(hmC.f[a, ])

# Save results
betas <- betas.Noob.bs.drn
save(U, mC, hmC, betas, pheno.bs, cpg.hmc, cpg.pfiltered, file = "DRN_ML_ReadyforEWAS.Rdata")

# Plotting the results
par(mfrow = c(1, 3))

# Density plot for unmethylated cytosine
plot(density(as.matrix(U)), main = "C", ylim = c(0, 10))

# Density plot for methylated cytosine
plot(density(as.matrix(mC)), main = "mC", ylim = c(0, 10))

# Density plot for hydroxymethylated cytosine
plot(density(as.matrix(hmC[cpg.hmc, ])), main = "hmC", ylim = c(0, 10))



































```{r performin sva, include=FALSE}
mydf<-na.omit(U.f)
sva.U<- sva(dat = as.matrix(mydf), mod = design, mod0 = design[ ,-2])
mydf<-na.omit(mC.f)
sva.mC<- sva(dat = as.matrix(mydf), mod = design, mod0 = design[ , -2])
hmC.f.f<-hmC.f[which(rownames(hmC.f)%in% cpg.hmc),]
mydf<-na.omit(hmC.f[b,])
sva.hmC<- sva(dat = as.matrix(mydf), mod = design, mod0 = design[ , -2])
save(sva.U, sva.mC, sva.hmC, file="sva.drn.83.Rdata")
```


```{r EWAS formula, include=TRUE}
 EWAS <- function(x, BraakStage, Ageatdeath, Gender, Chipid,  sv1 , sv2 , sv3, sv4,sv5){
	fit<-lm(BraakStage ~ as.numeric(x) + Ageatdeath + Gender + Chipid + sv1 + sv2 + sv3 + sv4 +sv5)
	return(cbind(coef(summary(fit))[2,], shapiro.test(fit$residuals)$p.value))
} 
```


```{r converting the vars to the correct vectors, include=FALSE}
setwd("/mnt/data1/Ehsan/DRN/")
load("sva.drn.83.Rdata")

pheno$BraakStage<-as.numeric(pheno$BraakStage)
pheno$Gender<-as.factor(pheno$Gender)
pheno$Ageatdeath<-as.numeric(pheno$Ageatdeath)
Chipid<-as.factor(substr((pheno$sentrix_ID),1,12))
BraakStage<-pheno$BraakStage
Ageatdeath<-pheno$Ageatdeath
Gender<-pheno$Gender

```

```{r EWAS.Unmodified , include=FALSE}
sv1<-sva.U$sv[,1]
sv2<-sva.U$sv[,2]
sv3<-sva.U$sv[,3]
sv4<-sva.U$sv[,4]
sv5<-sva.U$sv[,5]

cl<- makeCluster(32)
res.U<-t(parApply(cl, U.f , 1, EWAS, BraakStage, Ageatdeath, Gender , Chipid, sv1 ,sv2 ,sv3, sv4, sv5))

res.U<-res.U[,c(1:5)]
p.BH.U<-as.data.frame(p.adjust(res.U[,4], "BH"))
res.U<-cbind(res.U,p.BH.U)
colnames(res.U)<-c("U.Estimate","U.SE","U.t","U.P","resid.norm.U","BH.U")
head(res.U)
```

```{r EWAS formula, include=TRUE}
 EWAS <- function(x, BraakStage, Ageatdeath, Gender, Chipid,  sv1 , sv2 , sv3){
	fit<-lm(BraakStage ~ as.numeric(x) + Ageatdeath + Gender + Chipid + sv1 + sv2 + sv3 )
	return(cbind(coef(summary(fit))[2,], shapiro.test(fit$residuals)$p.value))
} 
```


```{r EWAS.mC , include=FALSE}
sv1<-sva.mC$sv[,1]
sv2<-sva.mC$sv[,2]
sv3<-sva.mC$sv[,3]
sv4<-sva.mC$sv[,4]
sv5<-sva.mC$sv[,5]
cl<- makeCluster(64)
res.mC<-t(parApply(cl, mC.f , 1, EWAS, BraakStage, Ageatdeath, Gender , Chipid, sv1 ,sv2 ,sv3))

res.mC<-res.mC[,c(1:5)]
p.BH.mC<-as.data.frame(p.adjust(res.mC[,4], "BH"))
res.mC<-cbind(res.mC,p.BH.mC)
colnames(res.mC)<-c("mC.Estimate","mC.SE","mC.t","mC.P","resid.norm.mC","BH.mC")

#```


#```{r EWAS.hmC , include=FALSE}
sv1<-sva.hmC$sv[,1]
sv2<-sva.hmC$sv[,2]
sv3<-sva.hmC$sv[,3]
sv4<-sva.hmC$sv[,4]
sv5<-sva.hmC$sv[,5]

res.hmC<-t(parApply(cl, hmC.f , 1, EWAS, BraakStage, Ageatdeath, Gender , Chipid, sv1 ,sv2 ,sv3))
res.hmC<-res.hmC[,c(1:5)]
colnames(res.hmC)<-c("hmC.Estimate","hmC.SE","hmC.t","hmC.P","resid.norm.hmC")
stopCluster(cl)
#```


#```{r Lambda , include=FALSE}
REULTS.DRN<-cbind(res.hmC,res.mC,res.U)
dim(REULTS.DRN)
#REMOVE CROSS-HYBRIDISING PROBES
cross<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
crosslist<-cross[,1]
length(crosslist)
#[1]44210  
REULTS.DRN.1 <-REULTS.DRN[ !rownames(REULTS.DRN) %in% crosslist,]
dim(REULTS.DRN.1)
##REMOVE PROBES with SNPs
snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
dim(snpProbes)
#[1] 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
dim(snpProbes)
# [1] 10888    18
SNPlist<-snpProbes[,1]
REULTS.DRN.2 <- REULTS.DRN.1[ ! rownames(REULTS.DRN.1) %in% SNPlist,]

##REMOVE bad probes
badProbes<-read.csv("/mnt/data1/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE, header = TRUE)
dim(badProbes)
#[1] 977  48
badlist<-badProbes[,1]
REULTS.DRN.3 <- REULTS.DRN.2[ ! rownames(REULTS.DRN.2) %in% badlist,]
dim(REULTS.DRN.1) # 
dim(REULTS.DRN.2) # 
dim(REULTS.DRN.3) # 

#epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7)
ANNOT<-data.frame(epicManifest)
###### QC 5.  Remove Chr X Y ########
sexchrX<-ANNOT[which(ANNOT$CHR=="X"),]
dim(sexchrX) #[1] 19120    48
sexchrY<- ANNOT[which(ANNOT$CHR== "Y"),]
dim(sexchrY) #[1] 561  48

REULTS.DRN.4 <- REULTS.DRN.3[!(rownames(REULTS.DRN.3) %in% sexchrX$IlmnID),]
REULTS.DRN.5 <- REULTS.DRN.4[!(rownames(REULTS.DRN.4) %in% sexchrY$IlmnID),]

ANNOT.sub<-ANNOT[which(ANNOT$Name%in% rownames(REULTS.DRN.5)),]


ANNOT.sub<-ANNOT.sub[order(ANNOT.sub$Name),]
REULTS.DRN.5<-REULTS.DRN.5[order(rownames(REULTS.DRN.5)),]

identical(as.character(rownames(REULTS.DRN.5)), as.character(ANNOT.sub$Name))

DRN.RESULTS.ANNOT<-cbind(REULTS.DRN.5,ANNOT.sub)
head(DRN.RESULTS.ANNOT)
table(DRN.RESULTS.ANNOT$CHR)
save(DRN.RESULTS.ANNOT, file="DRN.RESULTS.ANNOT.Rdata")
DRN.RESULTS.ANNOT10<-DRN.RESULTS.ANNOT[which(DRN.RESULTS.ANNOT$Name%in%rownames(mydf)),]


col<-c("hmC.Estimate","hmC.SE","hmC.t","hmC.P","mC.Estimate","mC.SE","mC.t","mC.P","U.Estimate","U.SE","U.t","U.P","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","Regulatory_Feature_Group")
hmC.p.ranked<-DRN.RESULTS.ANNOT10[order(DRN.RESULTS.ANNOT10[,4]),]
write.csv(hmC.p.ranked[1:100,col], file="top100.hmC.DRN.csv")


col<-c("hmC.Estimate","hmC.SE","hmC.t","hmC.P","mC.Estimate","mC.SE","mC.t","mC.P","U.Estimate","U.SE","U.t","U.P","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","Regulatory_Feature_Group")
mC.p.ranked<-DRN.RESULTS.ANNOT[order(DRN.RESULTS.ANNOT[,9]),]
write.csv(mC.p.ranked[1:100,col], file="top100.mC.DRN.csv")


col<-c("hmC.Estimate","hmC.SE","hmC.t","hmC.P","mC.Estimate","mC.SE","mC.t","mC.P","U.Estimate","U.SE","U.t","U.P","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","Regulatory_Feature_Group")
C.p.ranked<-DRN.RESULTS.ANNOT[order(DRN.RESULTS.ANNOT[,15]),]
write.csv(C.p.ranked[1:100,col], file="top100.C.DRN.csv")

```


```{r Lambda , include=FALSE}
chisq <- qchisq(1-(DRN.RESULTS.ANNOT10[,4]),1)
lambda = median(chisq)/qchisq(0.5,1)
lambda#1.065208

chisq <- qchisq(1-(DRN.RESULTS.ANNOT[,9]),1)
lambda = median(chisq)/qchisq(0.5,1)
lambda#1.122397

chisq <- qchisq(1-(DRN.RESULTS.ANNOT[,15]),1)
lambda = median(chisq)/qchisq(0.5,1)
lambda#1.043112 

pdf("qqplots.DRN.pdf", width=10, height=4)
par(mfrow=c(1,3))
qq(mC.p.ranked[,9], main="mC", pch=".")
title(sub="lambda=1.12",cex.lab=0.75)
qq(hmC.p.ranked[,4], main="hmC", pch=".")
title(sub="lambda=1.06",cex.lab=0.75)
qq(C.p.ranked[,15], main="C", pch=".")
title(sub="lambda=1.04",cex.lab=0.75)
dev.off()

p.BH.hmC<-as.data.frame(p.adjust(hmC.p.ranked[,4], "BH"))
head(p.BH.hmC,10)
-log10(0.05/nrow(hmC.p.ranked))
-log10(0.05/nrow(mC.p.ranked))

-log10(0.05/nrow(C.p.ranked))

p.BH.mC<-as.data.frame(p.adjust(mC.p.ranked[,9], "BH"))
head(p.BH.mC,10)


p.BH.C<-as.data.frame(p.adjust(C.p.ranked[,15], "BH"))
head(p.BH.C,10)
```

```{r DMR , include=TRUE}
library(dplyr)
DRN.RESULTS.ANNOT$MAPINFO.1<-as.numeric(DRN.RESULTS.ANNOT$MAPINFO)
DRN.RESULTS.ANNOT$MAPINFO.1.1<-(DRN.RESULTS.ANNOT$MAPINFO.1)+1
dmr <- data.frame("chrom" = 	paste0("chr", as.numeric(as.character(DRN.RESULTS.ANNOT$CHR))), 
					        "start" = 	DRN.RESULTS.ANNOT[rownames(DRN.RESULTS.ANNOT), "MAPINFO.1"],
                  "end" 	= 	DRN.RESULTS.ANNOT[rownames(DRN.RESULTS.ANNOT), "MAPINFO.1.1"],
                  "pvalue" = 	DRN.RESULTS.ANNOT[,15])

colnames(dmr) <- c("chrom", "start", "end", "pvalue")
dmr.U.DRN<-dmr[order(dmr[,1],dmr[,2]),]
head(dmr.U.DRN)
setwd("/home/ehsan/")
write.table(dmr.U.DRN, file="DRN.U.DMR.txt",sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)


DRN.RESULTS.ANNOT$MAPINFO.1<-as.numeric(DRN.RESULTS.ANNOT$MAPINFO)
DRN.RESULTS.ANNOT$MAPINFO.1.1<-(DRN.RESULTS.ANNOT$MAPINFO.1)+1
dmr <- data.frame("chrom" = 	paste0("chr", as.numeric(as.character(DRN.RESULTS.ANNOT$CHR))), 
					        "start" = 	DRN.RESULTS.ANNOT[rownames(DRN.RESULTS.ANNOT), "MAPINFO.1"],
                  "end" 	= 	DRN.RESULTS.ANNOT[rownames(DRN.RESULTS.ANNOT), "MAPINFO.1.1"],
                  "pvalue" = 	DRN.RESULTS.ANNOT[,9])

colnames(dmr) <- c("chrom", "start", "end", "pvalue")
dmr.mC.DRN<-dmr[order(dmr[,1],dmr[,2]),]
setwd("/home/ehsan/")
write.table(dmr.mC.DRN, file="DRN.mC.DMR.txt",sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)

DRN.RESULTS.ANNOT10$MAPINFO.1<-as.numeric(DRN.RESULTS.ANNOT10$MAPINFO)
DRN.RESULTS.ANNOT10$MAPINFO.1.1<-(DRN.RESULTS.ANNOT10$MAPINFO.1)+1
dmr <- data.frame("chrom" = 	paste0("chr", as.numeric(as.character(DRN.RESULTS.ANNOT10$CHR))), 
					        "start" = 	DRN.RESULTS.ANNOT10[rownames(DRN.RESULTS.ANNOT10), "MAPINFO.1"],
                  "end" 	= 	DRN.RESULTS.ANNOT10[rownames(DRN.RESULTS.ANNOT10), "MAPINFO.1.1"],
                  "pvalue" = 	DRN.RESULTS.ANNOT10[,4])

colnames(dmr) <- c("chrom", "start", "end", "pvalue")
dmr.mC.DRN<-dmr[order(dmr[,1],dmr[,2]),]
setwd("/home/ehsan/")
write.table(dmr.mC.DRN, file="DRN.hmC.DMR.txt",sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
```

comb-p pipeline -c 4 --dist 1000 --seed 0.01 --anno hg19 -p  out.DRN.U500 DRN.U.DMR.txt
comb-p pipeline -c 4 --dist 1000 --seed 0.01 --anno hg19 -p  out.DRN.mC500 DRN.mC.DMR.txt
comb-p pipeline -c 4 --dist 1000 --seed 0.01 --anno hg19 -p  out.DRN.hmC500 DRN.hmC.DMR.txt

```{r pressure, echo=TRUE}
#select deletion region
32063403	32064838

res.chr.6 <- DRN.RESULTS.ANNOT[which(DRN.RESULTS.ANNOT$CHR =="6"),]
dim(res.chr.6)
TNXB<-res.chr.6[(which(res.chr.6$MAPINFO > 32063403 & res.chr.6$MAPINFO < 32064838)),]
dim(TNXB)
TableResultsTNXB <- data.frame(rownames(TNXB),as.numeric(TNXB[,1]),
                                as.numeric(TNXB[,9]),
                                as.character(TNXB$CHR),
                                as.numeric(as.character(TNXB$MAPINFO)), 
                                rep("TNXB", times = 41))

colnames(TableResultsTNXB)<-c("CpGid","estimate","p-val","CHR","MAPINFO","Gene_Name")
allResult<-makeGRangesFromDataFrame(TableResultsTNXB,    #convert into a Grange object
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("CHR"),
                                    start.field=c("MAPINFO"),
                                    end.field=c("MAPINFO"),
                                    starts.in.df.are.0based=FALSE)
chr <- as.character(unique(seqnames(allResult)))
atrack <- AnnotationTrack(allResult, name = "CpG sites")
axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
knownGenes <- UcscTrack(genome = "hg19", 
                        chromosome = chr,
                        track = "UCSC Genes", 
                        from = 27143234 -1500 , 
                        to = 27143790 +100, 
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", 
                        rends = "exonEnds", 
                        gene = "name",
                        showId=TRUE, 
                        geneSymbol=TRUE, 
                        strand = "strand",
                        fill = "#8282d2", 
                        name = "UCSC Genes")


CPGtrack <- AnnotationTrack(CpGs, name = "CpG Island")
biomTrack <- BiomartGeneRegionTrack(genome = "hg19",chromosome = 7 , start = 27143234-100, end = 27143790 + 100 ,name = "ENSEMBL", 
                                    filter = list(with_ox_refseq_mrna = TRUE),col.line = NULL, col = NULL, stackHeight = 0.5, 
                                    collapseTranscripts ="meta", transcriptAnnotation = "symbol")

HOXA2.betas<-data.frame(betas.all[rownames(HOXA2),])
#all_Diff <- MeanSuicide-MeanCon
extra_info<-data.frame(TableResultsHOXA2$CHR,as.numeric(as.character(TableResultsHOXA2$MAPINFO))) # create extra information
colnames(extra_info)<-c("chr", "start")
extra_info$end<-as.numeric(extra_info$start)+10

colnames(extra_info)<-c("chr", "start", "end")
rownames(extra_info)<-rownames(HOXA2)
extra_info

HOXA2_results<-as.data.frame(cbind(HOXA2.betas,extra_info))

HOXA2_DMR<-makeGRangesFromDataFrame(HOXA2_results, #convert into a Grange object
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo=NULL,
                                     seqnames.field=c("chr"),
                                     start.field=c("start"),
                                     end.field=c("end"),
                                     starts.in.df.are.0based=FALSE)


pheno.all[226,9]="EuGEI"

dTrack <- DataTrack(HOXA2_DMR, 
                    genome = "hg19", 
                    name = "Betas",
                    groups = pheno.all$sampleinfo, 
                    type = c("boxplot","a","g"), 
                    col=c("red", "blue"), 
                    legend = TRUE, lwd=2)

pdf("HOXA2_DMR.pdf", height=5, width=7.5)
plotTracks(list(idxTrack,axTrack,CPGtrack,atrack,dTrack,biomTrack), from = 27143234 -100, to = 27143790 + 100,
           transcriptAnnotation="symbol", 
           stacking="squish",
           background.title="#40464C")
dev.off()
```
