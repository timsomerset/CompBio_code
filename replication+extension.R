library("ggplot2")
library("tidyverse")
library("data.table")
library("BiocManager")
library("GEOquery")
library("ComplexHeatmap")
library("magick")
library("circlize")
library("dplyr")

mutation <- fread("data/mutation.tsv")
samples <- fread("data/samples.tsv")

GSE <- getGEO("GSE81540", GSEMatrix = TRUE, destdir = "./data")
data1 <- pData(GSE$`GSE81540-GPL11154_series_matrix.txt.gz`)
data2 <-  pData(GSE$`GSE81540-GPL18573_series_matrix.txt.gz`)

compactify_data <- function(GEO_data){
  compact_data <- data.frame(ID = GEO_data$title,
                             SCAN_B_ID = gsub(".*S0(\\d+)", "S0\\1", GEO_data$characteristics_ch1),
                             PAM50 = GEO_data$`pam50 subtype:ch1`,
                             ER = GEO_data$`er status:ch1`,
                             PgR = GEO_data$`pgr status:ch1`,
                             HER2 = GEO_data$`her2 status:ch1`,
                             Ki67 = GEO_data$`ki67 status:ch1`,
                             NHG = GEO_data$`nhg:ch1`)
  compact_data$SCAN_B_ID <- gsub(".k.a.t", "", compact_data$SCAN_B_ID)

  # compact data includes the 405 samples used for training prediction, so we remove them here
  training <- grep("T", compact_data$ID)
  if(length(training) != 0){
    compact_data <- compact_data[-training,]
  }

  # also have replicates, and not sure what to do with them, so we remove them here
  # a few so far look like they're in agreement with their copies
  replicates <- grep("repl", compact_data$ID)
  if(length(replicates) != 0){
    compact_data <- compact_data[-replicates, ]
  }

  # add histological subtype from website dataset
  histological <- sapply(compact_data$SCAN_B_ID, function(ID){
    row <- which(samples$SAMPLE == ID)
    samples$Histological_Type[row]
  })
  histological <- do.call(c, lapply(histological, function(x){
    if(is.null(x) | length(x) == 0){NA}else{x}
  }))

  compact_data <- cbind(compact_data, data.frame(Hist = histological %>% unlist ))
  in_mutexplorer <- !is.na(compact_data$Hist)
  compact_data <- compact_data[in_mutexplorer,]

  compact_data
}
compact_data1 <- compactify_data(data1)
compact_data2 <- compactify_data(data2)

compact_data <- rbind(compact_data1, compact_data2)

biomark_matrix <- matrix(c(compact_data$Hist, compact_data$PAM50,
                        compact_data$ER, compact_data$PgR, compact_data$HER2,
                        compact_data$Ki67, compact_data$NHG), nrow = 7, byrow = TRUE)
row.names(biomark_matrix) <- c("Hist", "PAM50", "ER", "PgR", "HER2", "Ki67", "NHG")
colnames(biomark_matrix) <- compact_data$SCAN_B_ID


# create another dataset for mutations
# going to do this first, and then re-organise everything for plotting, seems to be simpler to implement this way
mutated_genes <- mutation$gene.symbol %>% unique
n_samples <- compact_data$SCAN_B_ID %>% unique %>% length
mutation_matrix <- sapply(compact_data$SCAN_B_ID, function(sample){
  rows <- which(mutation$SAMPLE == sample)
  mut_genes <- mutation$gene.symbol[rows]
  mut_type <- mutation$SnpEff.Effect.Class[rows]

  mutation_vec <- rep("None", length(mutated_genes))
  names(mutation_vec) <- mutated_genes
  mutation_vec[mut_genes] <- mut_type

  mutation_vec
})


# establish odering of genes by mutation number
mut_type_levels <- c("None", "Frameshift", "Nonsense", "In-frame indel", "Missense",
                     "Splicing", "Synonymous", "UTR", "Other")

mutation_counts <- sapply(mutated_genes, function(gene){
  sample_info <- mutation_matrix[gene,]
  sample_info <- factor(sample_info, levels = mut_type_levels) %>% as.integer
  n_mutated <- which(sample_info != 1) %>% length
  n_mutated
})
gene_order <- mutation_counts %>% sort(decreasing = TRUE) %>% names


# Now, compile ALL the data togther:
#         1. biomark_matrix for biomarkers
#         2. mutation_matrix for mutations

# define a new matrix identifying only if mutations are synonymous/nonsynonymous
# this is used for column ordering only
syn_nsyn_matrix <- sapply(as.vector(mutation_matrix), function(x){
  if(x == "None" | x == "Synonymous"){
    0
  }else{
    1
  }
}) %>% matrix(nrow = length(mutated_genes))
row.names(syn_nsyn_matrix) <- row.names(mutation_matrix)
colnames(syn_nsyn_matrix) <- colnames(mutation_matrix)

# we use the mega_matrix only for column ordering. Row ordering is already done with gene order
mega_matrix <- rbind(syn_nsyn_matrix, biomark_matrix)
# primary ordering by histological subtype, followed by order based on gene order
order_hierarchy <- list(mega_matrix["Hist",])
for(gene in gene_order){
  order_hierarchy[[length(order_hierarchy) + 1]] <- -as.numeric(mega_matrix[gene,])
}
column_order <- do.call("order", order_hierarchy)
column_order <- colnames(mega_matrix)[column_order]
row_order <- gene_order

# double check that ordered by gene and histological subtype is correct
mega_row_order <- c(row_order, row.names(biomark_matrix))
mega_matrix <- mega_matrix[mega_row_order, column_order]

# applying new ordering
mutation_matrix <- mutation_matrix[row_order, column_order]
biomark_matrix <- biomark_matrix[, column_order]

# PLOT-------------------------------
# remove sample names from plot
colnames(mutation_matrix) <- c()
colnames(biomark_matrix) <- c()

# side bar plot dataframe
non_synonymous <- sapply(mutated_genes, function(gene){
  syn_nsyn_matrix[gene,] %>% sum
})
non_synonymous <- non_synonymous[gene_order]
bar_plot_mat <- data.frame(synonymous = (mutation_counts[gene_order] - non_synonymous)/ncol(mutation_matrix),
                          n_synonymous = non_synonymous/ncol(mutation_matrix)) %>% as.matrix
bar_plot_text <- sapply(rowSums(bar_plot_mat), function(mut){
  paste0("(", round(mut*100) %>% as.character, "%)")
})
# colours for alteration type
colors <- structure(c("#ffffff", "#f52516", "#f5ab16", "#43d17c",
                      "#68edc4", "#68c8ed", "#c168ed", "#ed68e9",
                      "#c7c7c7"), names = mut_type_levels)
mutation_annotation <- rowAnnotation(`%` = anno_text(bar_plot_text),
                                     `% Mutant` = anno_barplot(bar_plot_mat,
                                                               gp = gpar(fill = c(synonymous = "red", n_synonymous = "blue")),
                                                               axis_param = list(direction = "reverse")),
                                     empty = anno_empty(border = FALSE)
)
mut_anno_legend <- Legend(labels = c("synonymous", "non-synonymous"),
                          title = "Translational Effect",
                          legend_gp = gpar(fill = c("red", "blue")))

biomark_annotation <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  `Histological Subtype` = biomark_matrix["Hist",],
  `Molecular Subtype` = biomark_matrix["PAM50",],
  `ER/PgR/HER2 Status` = cbind(ER = biomark_matrix["ER",], PgR = biomark_matrix["PgR",],
                               HER2 = biomark_matrix["HER2",]),
  `Ki67 Status` = biomark_matrix["Ki67",],
  `Histological Grade` = biomark_matrix["NHG",],
   col = list(`Histological Subtype` = c(Ductal = "yellow", Lobular = "red", Other = "gray"),
              `Molecular Subtype` = c(LumA = "lightblue", LumB = "cornflowerblue", Her2 = "pink", Basal = "red",
                                      Normal = "green"),
              `ER/PgR/HER2 Status` = c(`0` = "white", `1` = "black", `NA`= "gray"),
              `Ki67 Status` = c(`0` = "white", `1` = "black", `NA` = "gray"),
              `Histological Grade` = c(G1 = "cornflowerblue", G2 = "lightblue", G3 = "blue", `NA` = "gray")),
  annotation_name_side = list(side = "left"),
  annotation_legend_param = list(`Ki67 Status` = list(labels = c("Low", "High", "NA")),
                                 `ER/PgR/HER2 Status` = list(labels = c("Negative", "Positive", "NA")))
)


body <- Heatmap(mutation_matrix, col = colors, left_annotation = mutation_annotation,
                bottom_annotation = biomark_annotation, name = "Alteration type", row_names_side = "left")

draw(body, annotation_legend_list = list(mut_anno_legend))



# SURVIVAL PLOTS -----------------------------------------------------------------------
library(survival)
library(readr)
library(ggplot2)
library(scales)
library(survminer)
library(dplyr)
library(tidyverse)
library(gridExtra)

survival_data <- fread("data/survival_data.tsv")
survival_data <- fread("data/survival_samples.tsv")

# assumes survival_data is already imported
plot_KM <- function(treatment, covar_type){
  if(treatment == "Any"){
    treatment_status <- rep(1, nrow(survival_data)) %>% as.logical
  }else{
    treatment_status <- survival_data[[treatment]] %>% as.logical
  }
  survival_compact <- survival_data[treatment_status,]
  surv_df <- data.frame(SAMPLE = survival_compact$SAMPLE, time = survival_compact$OS_years, surv_status = survival_compact$OS_event,
                        mut_burden = survival_compact$tumor_mutational_burden)
  surv_df$mut_burden <- factor(surv_df$mut_burden, levels = c("Low", "High"))
  surv_samp <- surv_df$SAMPLE
  covariate_df <- cbind(survival_compact[,c(20,21,25)],
                        compact_data[compact_data$SCAN_B_ID %in% surv_samp, c("ER", "PgR", "HER2", "NHG")])
  surv_df <- cbind(surv_df, covariate_df)
  # clean up a litte, turn covariate data into numeric (mostly stored as characters)
  surv_df$ER <- factor(surv_df$ER, exclude = "NA")
  surv_df$HER2 <- factor(surv_df$HER2, exclude = "NA")
  surv_df$PgR <- factor(surv_df$PgR, exclude = "NA")
  surv_df$NHG <- factor(surv_df$NHG, exclude = "NA")
  surv_df$Node_group <- factor(surv_df$Node_group, exclude = "")

  #create survival variable object
  s <- Surv(surv_df$time, surv_df$surv_status)

  # simple single-parameter model
  s_fit <- coxph(s ~ surv_df$mut_burden)
  plot_fit <- surv_fit(s~surv_df$mut_burden, data = surv_df)
  #p_value <- survdiff(s~surv_df$mut_burden)$pvalue
  HR <- exp(s_fit$coefficients[1])
  HR_ci <- summary(s_fit)$conf.int[3:4]


  # cox multi-variable
  covar_set <- list(basic = c("mut_burden", "Age", "TumorSize", "Node_group"),
                    p = c("mut_burden", "Age", "TumorSize", "Node_group", "ER", "PgR", "HER2", "NHG"),
                    sq = c("mut_burden", "Age", "TumorSize", "Node_group", "ER", "PgR", "NHG"),
                    `$` = c("mut_burden", "Age", "TumorSize", "Node_group", "HER2", "NHG"),
                    `#` = c("mut_burden", "Age", "TumorSize", "Node_group", "NHG"))

  covar <- covar_set[[covar_type]]
  formula <- reformulate(response = "s", termlabels = covar)
  mv_fit <- coxph(formula, data = surv_df)
  mv_HR <- exp(mv_fit$coefficients[1])
  mv_HR_ci <- summary(mv_fit)$conf.int[1,3:4]


  my_plot = ggsurvplot(plot_fit, data = surv_df, palette = c("orange", "purple"),
                       title = gsub("Adjuvant_", "", treatment),
                       risk.table = T,
                       risk.table.y.text.col = TRUE, # colour risk table text annotations.
                       risk.table.y.text = FALSE,    # show bars instead of names in text annotations
                       risk.table.height = 0.17,
                       tables.theme = theme_cleantable() +
                         theme(plot.title = element_text(size=13, face="plain",
                                                         color="black", hjust=0,
                                                         margin=margin(b=-0.5, unit="pt"))),
                       tables.y.text = FALSE,
                       fontsize = 4,            # risk table font size

                       pval = TRUE,
                       pval.size = 4,
                       pval.coord = c(0.13, 0.3),

                       xlab = "Time after Diagnosis (years)",
                       ylab = "",
                       ggtheme = theme_classic() +
                         theme(
                           axis.title = element_text(size=13, face="plain", color="black"),
                           axis.text = element_text(size=13, face="plain", colour="black"),
                           axis.ticks.length=unit(.15, "cm"),
                           legend.text=element_text(size=11, face="plain"),
                           plot.title = element_text(size="16", face="bold", colour="black",
                                                     hjust=0.5, margin=margin(b=1, unit="pt"))
                         ),
                       break.time.by = 1,

                       legend.title = ""
  )
  # add text for hazard ratios

  HR_lab <- paste0("HR, ", round(HR, 2), " (", round(HR_ci[1], 2), " to ", round(HR_ci[2], 2), ")" )
  mvHR_lab <- paste0("MV ", covar_type, " HR, ", round(mv_HR, 2),
                     " (", round(mv_HR_ci[1] ,2),
                     " to ", round(mv_HR_ci[2], 2), ")")

  my_plot$plot = my_plot$plot +
    ggplot2::annotate("text", x = 3, y = 0.5, label = HR_lab) +
    ggplot2::annotate("text", x = 3, y = 0.4, label = mvHR_lab)
  my_plot
}

plot_KM("Adjuvant_Endo", "#")

treatment <- c("Any", colnames(survival_data[,2:15])[c(14, 4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9, 13)])
covar_type <- c("p", "p", "$", "$", "basic", "#", "sq", "sq", "$", "basic", "basic", "$", "#", "sq", "#")

for(i in 1:4){
  sublist <- list()
  for(k in (4*(i-1)+1):min((4*i), length(treatment))){
    index <- length(sublist)+1
    sublist[[index]] <- plot_KM(treatment[k], covar_type[k])
  }
  name <- paste0("splot", i, ".jpeg")
  ggsave(name, arrange_ggsurvplots(sublist, print = FALSE, ncol = 2, nrow = 2, width = 20, height = 20))
}


# CO-MUTATION IDENTIFICATION -------------------------------------------------------------------------------------------
library("maftools")
library("dplyr")
library("survival")
library("survminer")
library("ggplot2")


# create a MAF file from the mutation dataset. This doesn't need to be re-run every time, as long as mut_MAF.tsv
# exists in the working directory!
mutations <- fread("data/real_mutation_file.tsv")


send_to_MAF <- function(df, file_name){
  variant_dict <- c(Missense = "Missense_Mutation", UTR = "RNA", Other = "IGR",
                    Synonymous = "Silent", Splicing = "Splice_Site",
                    Nonsense = "Nonsense_Mutation", Frameshift = "Frame_Shift",
                    `In-frame indel` = "In_Frame_")
  variant_type <- sapply(df$TYPE, function(x){
    if(x == "SNV"){
      return("SNP")
    }
    if(x == "Deletion"){
      return("DEL")
    }
    if(x == "Insertion"){
      return("INS")
    }
  })
  variant_class <- sapply(1:nrow(df), function(x){
    original <- df$SnpEff.Effect.Class[x]
    class <- variant_dict[original]
    if(class == "Frame_Shift"){
      if(df$TYPE[x] == "Deletion"){
        class <- "Frame_Shift_Del"
      }else{
        class <- "Frame_Shift_Ins"
      }
    }
    if(class == "In_Frame_"){
      if(df$TYPE[x] == "Deletion"){
        class <- "In_Frame_Del"
      }else{
        class <- "In_Frame_Ins"
      }
    }
    class
  })

  df <- data.frame(`Hugo_Symbol` = df$gene.symbol, `Entrez_Gene_ID` = df$Entrez.ID,
                          Center = rep("SCAN_B", nrow(df)), NCBI_Build = rep("GRch38", nrow(df)),
                          Chromosome = df$CHROM, Start_Position = df$POS, End_Position = as.numeric(df$POS) + (sapply(df$ALT, length)-1),
                          Strand = rep("+", nrow(df)), Variant_Classification = variant_class,
                          Variant_Type = variant_type, Reference_Allele = df$REF, Tumor_Seq_Allele1 = df$ALT,
                          Tumor_Seq_Allele2 = df$ALT, Tumor_Sample_Barcode = df$SAMPLE, Protein_Change = df$SnpEff.Prot.Change)
  file <- paste0(file_name, ".tsv")
  fwrite(df, file = file, sep = "\t")
}


send_to_MAF(mutations, "mut_MAF")

N = samples$SAMPLE[which(samples$HER2 == "NEG" & samples$ER_1perc == "NEG" & samples$PgR_1perc == "NEG")]
HoRp_HERn = samples$SAMPLE[which(samples$HER2 == "NEG" & samples$ER_1perc == "POS" & samples$PgR_1perc == "POS")]
HoRp_HERp = samples$SAMPLE[which(samples$HER2 == "POS" & samples$ER_1perc == "POS" & samples$PgR_1perc == "POS")]
HoRn_HERp = samples$SAMPLE[which(samples$HER2 == "POS" & (samples$ER_1perc == "NEG" | samples$PgR_1perc == "NEG"))]
lob = samples$SAMPLE[which(samples$Histological_Type == "Lobular")]
duc = samples$SAMPLE[which(samples$Histological_Type == "Ductal")]

mutations_nn <- mutations[which(mutations$SAMPLE %in% N),] %>% send_to_MAF(file_name = "mutations_nn")
mutations_pn <- mutations[which(mutations$SAMPLE %in% HoRp_HERn),] %>% send_to_MAF(file_name = "mutations_pn")
mutations_np <- mutations[which(mutations$SAMPLE %in% HoRn_HERp),] %>% send_to_MAF(file_name = "mutations_np")
mutations_pp <- mutations[which(mutations$SAMPLE %in% HoRp_HERp),] %>% send_to_MAF(file_name = "mutations_pp")
mutations_lob <- mutations[which(mutations$SAMPLE %in% lob),] %>% send_to_MAF(file_name = "mutations_lob")
mutations_duc <- mutations[which(mutations$SAMPLE %in% duc),] %>% send_to_MAF(file_name = "mutations_duc")

nn_df <- somaticInteractions(read.maf(maf = "mutations_nn.tsv"), genes = gene_order, pvalue = c(0.05, 0.01))
pn_df <- somaticInteractions(read.maf(maf = "mutations_pn.tsv"), genes = gene_order, pvalue = c(0.05, 0.01))
np_df <- somaticInteractions(read.maf(maf = "mutations_np.tsv"), genes = gene_order, pvalue = c(0.05, 0.01))
pp_df <- somaticInteractions(read.maf(maf = "mutations_pp.tsv"), genes = gene_order, pvalue = c(0.05, 0.01))
lob_df <- somaticInteractions(read.maf(maf = "mutations_lob.tsv"), genes = gene_order, pvalue = c(0.05, 0.01),
                              geneOrder = gene_order)
duc_df <- somaticInteractions(read.maf(maf = "mutations_duc.tsv"), genes = gene_order, pvalue = c(0.05, 0.01),
                              geneOrder = gene_order)


extract_sig_comut <- function(df){
  sig_df <- df[which(df$pAdj < 0.05),]
  outdf <- data.frame(type = sig_df$Event, gene1 = sig_df$gene1, gene2 = sig_df$gene2)
  outdf
}
nn_sig <- extract_sig_comut(nn_df)
pn_sig <- extract_sig_comut(pn_df)
np_sig <- extract_sig_comut(np_df)
pp_sig <- extract_sig_comut(pp_df)
lob_sig <- extract_sig_comut(lob_df)
duc_sig <- extract_sig_comut(duc_df)


# This is where the analysis begins, with mut_MAF.tsv in the work direc already
mut_MAF <- read.maf(maf = "mut_MAF.tsv")
gene_order <- fread("gene_order.tsv", header = FALSE) %>% unlist %>% unname
interaction_df <- somaticInteractions(maf = mut_MAF, genes = gene_order, pvalue = c(0.05, 0.01),
                                      geneOrder = gene_order)

ordered_df <- interaction_df[order(match(interaction_df$gene1, gene_order), match(interaction_df$gene2, gene_order),
                                   interaction_df$pAdj),]
ordered_df <- ordered_df[ordered_df$pAdj < 0.05,]


# read in data from other parts of the analysis
compact_data <- fread("data/compact_data.csv")
survival_data <- fread("data/survival_data.tsv")

# now create mutation dataset recording only genes in gene_order
col_classes <- c("character", rep("integer", 30))
muts_bysample <- lapply(compact_data$SCAN_B_ID, function(sample){
  df_row <- read.table(text = "", colClasses = col_classes, col.names = c("SCAN_B_ID", gene_order))
  row_list <- list()
  row_list[[1]] <- sample
  mut_subset <- mutations[which(mutations$SAMPLE == sample),]
  mut_list <- lapply(gene_order, function(gene){
    if(gene %in% mut_subset$gene.symbol){
      1
    }else{
      0
    }
  })
  row_list <- c(row_list, mut_list)
  names(row_list) <- c("SCAN_B_ID", gene_order)
  df_row <- rbind(df_row, row_list)
  df_row
})
mutation_data <- do.call(rbind, Map(data.frame, muts_bysample))

comut_data <- merge(survival_data, mutation_data, by = "SCAN_B_ID") %>% data.frame
# make sure things are factored properly
comut_data$HER2 <- factor(comut_data$HER2, exclude = "")
comut_data$NHG <- factor(comut_data$NHG, exclude = "")
comut_data$Node_group <- factor(comut_data$Node_group, exclude = "")
comut_data$ER_10perc <- factor(comut_data$ER_10perc, exclude = "")
comut_data$PgR_10perc <- factor(comut_data$PgR_10perc, exclude = "")
comut_data$tumor_mutational_burden <- factor(comut_data$tumor_mutational_burden, levels = c("Low", "High"))
comut_data$`Adjuvant_Any_Treatment` <- (comut_data$Adjuvant_No_Systemic + 1)%%2

# survival analysis + plot function (equivalent to but slightly different from plot_KM)
analyse_survival <- function(gene1, gene2, comut_data, treatment){
  treatment_rows <- comut_data[[treatment]] %>% as.logical
  surv_df <- comut_data[treatment_rows,]

  # append on a column indicating comutation
  # comutation column is a 4-level factor: mutation in none, mutation in gene1, mutation in gene2, comutation
  gene_cols <- c(grep(gene1, colnames(comut_data)), grep(gene2, colnames(comut_data)))
  comutation <- (surv_df[, gene_cols[1]] + 2*surv_df[,gene_cols[2]] ) %>% factor()
  surv_df <- cbind(surv_df, comutation)
  if(surv_df[[treatment]] %>% sum(na.rm = T) < 14){
    print("Insufficient smaple size in treatment group")
    break
  }

  # cox multi-variable covariates
  covar <- c("comutation", "tumor_mutational_burden", "Age", "TumorSize", "Node_group", "ER_10perc", "PgR_10perc",
             "HER2", "NHG")

  # comparison to wild type
  s <- Surv(surv_df$OS_years, surv_df$OS_event)
  s_fit <- coxph(s ~ surv_df$comutation)
  plot_fit <- surv_fit(s~comutation, data = surv_df)
  HR_wt <- exp(s_fit$coefficients[3]) %>% round(digits = 2)
  HR_ci_wt <- summary(s_fit)$conf.int[3,3:4] %>% round(digits = 2)
  p_wt <- survdiff(s ~ surv_df$comutation)$p %>% round(digits = 4)

  formula <- reformulate(response = "s", termlabels = covar)
  mv_fit <- coxph(formula, data = surv_df)
  mv_HR_wt <- exp(mv_fit$coefficients[3]) %>% round(digits = 2)
  mv_HR_ci_wt <- summary(mv_fit)$conf.int[3,3:4] %>% round(digits = 2)
  p_mvwt <- summary(mv_fit)$coefficients[3,5] %>% round(digits = 4)

  # comparison to gene1
  gene1_df <- surv_df[which(surv_df$comutation == 1 | surv_df$comutation == 3),]
  gene1_df$comutation <- factor(gene1_df$comutation)
  s1 <- Surv(gene1_df$OS_years, gene1_df$OS_event)
  s_fit <- coxph(s1 ~ gene1_df$comutation)
  HR_gene1 <- exp(s_fit$coefficients) %>% round(digits = 2)
  HR_ci_gene1 <- summary(s_fit)$conf.int[,3:4] %>% round(digits = 2)
  p_gene1 <- survdiff(s1 ~ gene1_df$comutation)$p %>% round(digits = 4)

  formula <- reformulate(response = "s1", termlabels = covar)
  mv_fit <- coxph(formula, data = gene1_df)
  mv_HR_gene1 <- exp(mv_fit$coefficients[1]) %>% round(digits = 2)
  mv_HR_ci_gene1 <- summary(mv_fit)$conf.int[1,3:4] %>% round(digits = 2)
  p_mvgene1 <- summary(mv_fit)$coefficients[1,5] %>% round(digits = 4)

  # comparison to gene2
  gene2_df <- surv_df[which(surv_df$comutation ==2 | surv_df$comutation == 3),]
  gene2_df$comutation <- factor(gene2_df$comutation)
  s1 <- Surv(gene2_df$OS_years, gene2_df$OS_event)
  s_fit <- coxph(s1 ~ gene2_df$comutation)
  HR_gene2 <- exp(s_fit$coefficients) %>% round(digits = 2)
  HR_ci_gene2 <- summary(s_fit)$conf.int[,3:4] %>% round(digits = 2)
  p_gene2 <- survdiff(s1 ~ gene2_df$comutation)$p %>% round(digits = 4)

  formula <- reformulate(response = "s1", termlabels = covar)
  mv_fit <- coxph(formula, data = gene2_df)
  mv_HR_gene2 <- exp(mv_fit$coefficients[1]) %>% round(digits = 2)
  mv_HR_ci_gene2 <- summary(mv_fit)$conf.int[1,3:4] %>% round(digits = 2)
  p_mvgene2 <- summary(mv_fit)$coefficients[1,5] %>% round(digits = 4)

  title <- gsub("Adjuvant", "", treatment)
  title <- gsub("Cyto", "Chemo", title)
  title <- gsub("_", " ", title)

  my_plot = ggsurvplot(plot_fit, data = surv_df, palette = c("orange", "purple", "lightgreen", "cornflowerblue"),
                       title = title,
                       risk.table = T,
                       risk.table.y.text.col = TRUE, # colour risk table text annotations.
                       risk.table.y.text = FALSE,    # show bars instead of names in text annotations
                       risk.table.height = 0.17,
                       tables.theme = theme_cleantable() +
                         theme(plot.title = element_text(size=13, face="plain",
                                                         color="black", hjust=0,
                                                         margin=margin(b=-0.5, unit="pt"))),
                       tables.y.text = FALSE,
                       fontsize = 4,            # risk table font size


                       xlab = "Time after Diagnosis (years)",
                       ylab = "",
                       ggtheme = theme_classic() +
                         theme(
                           axis.title = element_text(size=13, face="plain", color="black"),
                           axis.text = element_text(size=13, face="plain", colour="black"),
                           axis.ticks.length=unit(.15, "cm"),
                           legend.text=element_text(size=11, face="plain"),
                           plot.title = element_text(size="16", face="bold", colour="black",
                                                     hjust=0.5, margin=margin(b=1, unit="pt"))
                         ),
                       break.time.by = 1,
                       legend.labs = c("wild-type", gene1, gene2, "co-mutation"),
                       legend.title = ""
  )
  # add text for hazard ratios
  gene1_lab <- paste0(gene1, " HR:", HR_gene1, " (", HR_ci_gene1[1],",", HR_ci_gene1[2], "), p = ", p_gene1,"; ",
                      gene1, " mvHR:", mv_HR_gene1, " (", mv_HR_ci_gene1[1],",", mv_HR_ci_gene1[2], "), p = ", p_mvgene1)
  gene2_lab <- paste0(gene2, " HR:", HR_gene2, " (", HR_ci_gene2[1],",", HR_ci_gene2[2], "), p = ", p_gene2, "; ",
                      gene2, " mvHR:", mv_HR_gene2, " (", mv_HR_ci_gene2[1],",", mv_HR_ci_gene2[2], "), p =", p_mvgene2)
  wt_lab <- paste0("wt HR:", HR_wt, " (", HR_ci_wt[1],",", HR_ci_wt[2], "), p = ", p_wt, "; ",
                   "wt mvHR:", mv_HR_wt, " (", mv_HR_ci_wt[1],",", mv_HR_ci_wt[2], "), p = ", p_mvwt)

  my_plot$plot = my_plot$plot +
    ggplot2::annotate("text", x = 3, y = 0.5, label = gene1_lab) +
    ggplot2::annotate("text", x = 3, y = 0.4, label = gene2_lab) +
    ggplot2::annotate("text", x = 3, y = 0.6, label = wt_lab)
  my_plot
}

# select specific gene pairs, and run survival analysis
# preliminary choice of treatments: Endo only, Endo + chemo \pm any, HER2 \pm any
treatments <- c("Adjuvant_Any_Treatment", "Adjuvant_No_Systemic", "Adjuvant_Endo", "Adjuvant_Cyto", "Adjuvant_HER2")
gene1 <- c("AKT1", "TP53", "TP53", "PIK3CA", "TP53")
gene2 <- c("PIK3CA", "MAP3K1", "PIK3CA", "CBFB", "PLEC")
plots <- lapply(1:length(gene1), function(x){
                lapply(treatments, function(treatment){
                  analyse_survival(gene1[x], gene2[x], comut_data, treatment)
                })
})

for(i in 1:length(plots)){
  folder <- paste0("figures/", gene1[i], "-", gene2[i])
  for(j in 1:length(treatments)){
    file <- paste0(folder, "/", treatments[j], ".jpeg")
    ggsave(file, ggsave_workaround(plots[[i]][j][[1]]), width = 10, height = 7)
  }
}