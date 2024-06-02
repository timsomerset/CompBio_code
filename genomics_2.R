library("edgeR")
library("org.Dm.eg.db")
library("splines")
library("DESeq2")
library("gamlss")
library("apeglm")
library("ggvenn")
library("ggpubr")

CountFile <- paste("http://bowtie-bio.sourceforge.net/recount",
                   "countTables",
                   "modencodefly_count_table.txt", sep="/")
Counts <- read.delim(url(CountFile), row.names=1)
SampleFile <- paste("http://bowtie-bio.sourceforge.net/recount",
                    "phenotypeTables",
                    "modencodefly_phenodata.txt", sep="/")
Samples <- read.delim(SampleFile, row.names=1, sep=" ", stringsAsFactors=FALSE)

# Sum counts by stage of sample
PooledCounts <- sumTechReps(Counts, ID = Samples$stage)
# each column is a stage of a sample
# each row is a gene, using FlyBase gene ID's
Hours <- seq(from=2, to=24, by=2)
Time <- paste0(Hours,"hrs")
y <- DGEList(counts=PooledCounts[,1:12], group=Time)
# Gene annotation
multiVals <- function(x) paste(x,collapse=";")
Symbol <- mapIds(org.Dm.eg.db, keys=rownames(y), keytype="FLYBASE",
                 column="SYMBOL", multiVals=multiVals)
y$genes <- data.frame(Symbol=Symbol, stringsAsFactors=FALSE)
HasSymbol <- y$genes$Symbol != "NA"
y <- y[HasSymbol, , keep.lib.sizes=FALSE]
# Filtering and normalisation
keep <- filterByExpr(y) #removes sample with insufficient counts
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y) # normalised counts
norm_counts <- y


#---------------------------------------------
runEdgeR <- function(degree){
  d <- degree
  X <- poly(Hours, degree=d)
  design <- model.matrix(~ X)
  y <- estimateDisp(norm_counts, design)
  fit <- glmQLFit(y,design)
  fit <- glmQLFTest(fit, coef = 2:(d+1)) # perform f-test on all coeffcients (including intercept)
  tab <- as.data.frame(topTags(fit, n = 30))

  AIC.vec <- c()
  chisq.vec <- c()
  for(gene in row.names(norm_counts)){
    AIC.vec[length(AIC.vec)+1] <- 2*(degree+1) +fit[gene,]$deviance
  }


  logCPM.obs <- cpm(y, log=TRUE, prior.count=fit$prior.count)
  logCPM.fit <- cpm(fit, log=TRUE)
  par(mfrow=c(2,2))
  for(i in 1:4) {
    FlybaseID <- row.names(tab)[i]
    Symbol <- tab$Symbol[i]
    logCPM.obs.i <- logCPM.obs[FlybaseID,]
    logCPM.fit.i <- logCPM.fit[FlybaseID,]
    plot(Hours, logCPM.obs.i, ylab="log-CPM", main=Symbol, pch=16)
    lines(Hours, logCPM.fit.i, col="red", lwd=2)
  }
  par(mfrow=c(1,1))

  par(mfrow = c(1,2))
  plotQLDisp(fit)
  plotBCV(y)
  par(mfrow = c(1,1))


  list(fit = fit, AIC = AIC.vec, chisq = chisq.vec)
}
edgeR_linear <- runEdgeR(1)
edgeR_quadratic <- runEdgeR(2)
edgeR_quartic <- runEdgeR(4)

violin_data <- data.frame(log_AIC = log(c(edgeR_linear$AIC,
                                          edgeR_quadratic$AIC, edgeR_quartic$AIC)),
                          Order = factor(rep(c(1,2,4), each = 10973)))
ggplot(violin_data, aes(x = Order, y = log_AIC, color = Order))+
  geom_violin() +
  geom_boxplot(width = 0.1)+
  scale_color_brewer(palette="Dark2")

linear_DE <- edgeR_linear$fit$table[which(edgeR_linear$fit$table$PValue < 0.05),]
quadratic_DE <- edgeR_quadratic$fit$table[which(edgeR_quadratic$fit$table$PValue < 0.05),]
quartic_DE <- edgeR_quartic$fit$table[which(edgeR_quartic$fit$table$PValue < 0.05),]


# DESEQ2----------------------------------------------------------------------------------------------------------------
CountFile <- paste("http://bowtie-bio.sourceforge.net/recount",
                   "countTables",
                   "modencodefly_count_table.txt", sep="/")
Counts <- read.delim(url(CountFile), row.names=1)
SampleFile <- paste("http://bowtie-bio.sourceforge.net/recount",
                    "phenotypeTables",
                    "modencodefly_phenodata.txt", sep="/")
Samples <- read.delim(SampleFile, row.names=1, sep=" ", stringsAsFactors=FALSE)
# only include data sampled from embryos
sample_info <- Samples[grep("Embryos", Samples$stage),]
# create embryo count dataset
embryo_counts <- Counts[, row.names(sample_info)]

PooledCounts <- sumTechReps(Counts, ID = Samples$stage)
embryo_counts2 <- PooledCounts[,1:12]


# create new column in sample dataset called hours
Hours <- seq(from=2, to=24, by=2)
sample_info$Hours <- Hours[factor(sample_info$stage)]
# extract the desired information from sample_info
# write it in a new form that is acceptable to DESeq
colData3 <- data.frame(Hours = Hours,
                       Hours2 = Hours^2,
                       Hours3 = Hours^3,
                       Hours4 = Hours^4
)
rownames(colData3) <- colnames(embryo_counts)

# alternative definition of col data that hard codes the polynomial columns
colData2 <- data.frame(Hours = sample_info[colnames(embryo_counts), 3],
                       Hours2 = (sample_info[colnames(embryo_counts), 3])^2,
                       Hours3 = (sample_info[colnames(embryo_counts), 3])^3,
                       Hours4 = (sample_info[colnames(embryo_counts), 3])^4 )
rownames(colData2) <- colnames(embryo_counts)


runDESeq2 <- function(degree, count_data, colData){
  d <- degree
  embryo_counts <- count_data
  regression <- paste0("~", paste(colnames(colData)[1:d], collapse  = "+"))
  y <- DESeqDataSetFromMatrix(embryo_counts, colData = colData,
                              design = formula(regression))
  y <- y[which(rowSums(counts(y)) >= 10),] # filter out low counts
  y <- estimateSizeFactors(y) # normalise counts
  y <- estimateDispersions(y)
  y <- nbinomWaldTest(y)
  raw_results <- rowData(y)
  result <- results(y, name = "Hours")

  chisq.vec <- c()
  AIC.vec <- c()
  for(gene in row.names(raw_results)){
    fitted <- assays(y)[["mu"]][gene,]
    observed <- as.numeric(embryo_counts[gene,])

    chisq.vec[length(chisq.vec)+1] <- sum((observed - fitted)^2/fitted)/(12-degree-1)

    AIC.vec[length(AIC.vec)+1] <- 2*(d+1) + raw_results[gene,]$deviance
  }


  # plots
  title = paste0("Plynomial Model of order", d)
  par(mfrow = c(2,1))
  plotDispEsts(y, main = title) # DESeq quality plot
  plotMA(result, main = title) # mean expression and LFC needed for significance of linear coef
  par(mfrow = c(1,1))


  list(design = design, y = y, chisq = chisq.vec, AIC = AIC.vec)
}

deseq_linear <- runDESeq2(1, embryo_counts2, colData3)
deseq_quadratic <- runDESeq2(2, embryo_counts2, colData3)
deseq_quartic <- runDESeq2(4, embryo_counts2, colData3)


violin_data <- data.frame(AIC = c(deseq_linear$AIC, deseq_quadratic$AIC, deseq_quartic$AIC),
                          Order = factor(rep(c(1,2,4), each = 12602)))
ggplot(violin_data, aes(x = Order, y = AIC, color = Order))+
  geom_violin() +
  geom_boxplot(width = 0.1)+
  scale_color_brewer(palette="Dark2")

#volcano plots
edgeR_toptags <- topTags(edgeR_linear$fit, n = nrow(edgeR_linear$fit$genes))
edgeR_results <- edgeR_toptags$table
edgeR_results$pass <- as.numeric(edgeR_results$PValue < 0.1
                                   & (edgeR_results$logFC < -9 | edgeR_results$logFC > 9))
deseq_results <- results(deseq_linear$y)
deseq_results$Symbol <- data.frame(Symbol[row.names(deseq_results)], stringsAsFactors = FALSE)
deseq_results$pass <- as.numeric(deseq_results$padj < 0.1
                                   & (deseq_results$log2FoldChange < -0.35 |
  deseq_results$log2FoldChange > 0.35))


plot(edgeR_results$logFC, -10*log(edgeR_results$PValue), col = edgeR_results$pass+1,
     xlab = "LogFC", ylab = "-10log(adj_p_value)", main = "EdgeR: Linear model")
plot(deseq_results$log2FoldChange, -10*log(deseq_results$padj),
     col = deseq_results$pass + 1, xlab = "LogFC",
     ylab = "-10log(adj_p_value)", main = "DESeq2: Linear model")

edgeR_DE <- edgeR_results[which(edgeR_results$pass == 1),]
deseq_DE <- deseq_results[which(deseq_results$pass == 1),]
intersect <- intersect(unlist(deseq_DE$Symbol), unlist(edgeR_DE$Symbol))

edgeR_best <- edgeR_DE[order(edgeR_DE$PValue),]
edgeR_best <- edgeR_best[1:4,]
deseq_best <- deseq_DE[order(deseq_DE$padj), ]
deseq_best <- deseq_best[1:4,]

counts <- c()
par(mfrow = c(2,2))
for(gene in row.names(deseq_best)){
  fitted <- log(assays(deseq_linear$y)[["mu"]][gene,])
  observed <- log(assays(deseq_linear$y)[["counts"]][gene,])

  plot(Hours, observed, main = deseq_best[gene,]$Symbol, ylab = "log(counts)")
  lines(Hours, fitted, col = "red", lwd = 3)
}
par(mfrow = c(1,1))


# comparison of p-value calling
edgeR_toptags <- topTags(edgeR_linear$fit, n = nrow(edgeR_linear$fit$genes))
edgeR_results <- edgeR_toptags$table
edgeR_results$pass <- as.numeric(edgeR_results$PValue < 0.1)

deseq_results <- results(deseq_linear$y)
deseq_results$Symbol <- data.frame(Symbol[row.names(deseq_results)],
                                   stringsAsFactors = FALSE)
deseq_results$pass <- as.numeric(deseq_results$padj < 0.1)

edgeR_results <- edgeR_results[which(edgeR_results$pass == 1),]
deseq_results <- deseq_results[which(deseq_results$pass == 1),]

intersect <- intersect(unlist(deseq_results$Symbol), unlist(edgeR_results$Symbol))

venn_data <- list(DESeq2 = row.names(deseq_results), EdgeR = row.names(edgeR_results))
ggvenn(venn_data)


edgeR_results <- edgeR_results[order(unlist(edgeR_results$Symbol)),]
edgeR_results <- edgeR_results[which(unlist(edgeR_results$Symbol) %in% intersect),]
deseq_results <- deseq_results[order(unlist(deseq_results$Symbol)),]
deseq_results <- deseq_results[which(unlist(deseq_results$Symbol) %in% intersect),]

deseq_counts <- sapply(row.names(deseq_results), function(x)
  mean(assays(deseq_linear$y)[["counts"]][x,]))


deseq_shrink <- lfcShrink(deseq_linear$y, coef = 2)
deseq_shrink$Symbol <- data.frame(Symbol[row.names(deseq_shrink)], stringsAsFactors = FALSE)
deseq_shrink <- deseq_shrink[order(unlist(deseq_shrink$Symbol)),]
deseq_shrink <- deseq_shrink[which(unlist(deseq_shrink$Symbol) %in% intersect),]


plot(edgeR_results$logCPM, log(deseq_counts, base = 2), xlab = "EdgeR log CPM",
     ylab = "DESeq2 log CPM")
x <- seq(-5, 15, 0.1)
lines(x, 5.744 + 1.025*x, col = "red", lwd = 3)

density_data <- data.frame(log_counts = c(edgeR_results$logCPM, log(deseq_counts, base = 2)),
                           tool = rep(c("EdgeR", "DESeq2"), each = 7651))
ggdensity(density_data, x = "log_counts", color = "tool", fill = "tool",
          palette = c("#01FF0D", "#FF0101"), alpha = 0.1)
