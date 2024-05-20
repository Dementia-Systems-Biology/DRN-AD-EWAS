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



		
