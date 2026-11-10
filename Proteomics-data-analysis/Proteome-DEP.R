library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggforce)
set.seed(1234)

## define a function
DE_protein <- function(condition1_samples_vector=c("DKO1", "DKO2", "DKO3", "DKO4"), 
                          condition2_samples_vector=c("WT1", "WT2", "WT3", "WT4"),
                          condition1_vector=c("DKO", "DKO", "DKO", "DKO"), 
                          condition2_vector=c("WT", "WT", "WT", "WT"),
                          sample_type1="DKO", sample_type2 = "WT",
                          cond1vs2="sample_typeDKO-sample_typeWT",
                          cond_levels=c("sample_typeDKO", "sample_typeWT")
) {
  data <- subset(df1, select = c(
    condition1_samples_vector,
    condition2_samples_vector
  ))
  # Set missing values to 0:
  na.mat <- is.na(data) 
  data[na.mat] <- 0
  df_new <- data[apply(data, 1, function(x) all(x != 0)), ]
  data <- df_new
  # Log transform data:
  dat.trans <- log2(data)
  ## wide to long format for ggplot2
  dat.trans.long <- tidyr::pivot_longer(as.data.frame(dat.trans),
                                        cols      = c(condition1_samples_vector,
                                                      condition2_samples_vector),  
                                        names_to  = "condition",         
                                        values_to = "value"  
  )
  ## boxplot
  dat.trans.long$condition <- factor(dat.trans.long$condition, levels = c(condition1_samples_vector,
                                                                          condition2_samples_vector))
  p1 <- ggplot(dat.trans.long, aes(x = condition, y = value)) +
    geom_boxplot() +
    theme_classic() +
    labs(title = "Log2 transformed data") +
    theme(plot.title = element_text(hjust=0.5))

  dat.filtered <- dat.trans
  coords <- plotMDS(dat.filtered, plot = F)
  plotdat <- data.frame(sample=colnames(dat.filtered), x= coords$x, y=coords$y)
  plotdat$sample_type <- c(
    condition1_vector, condition2_vector
  )
  p3 <- ggplot(plotdat, aes(x = x, y = y, color = sample_type)) +
    geom_point() +  
    geom_mark_hull(aes(fill = sample_type)) +
    # geom_encircle(data = subset(plotdat, sample_type == sample_type1), color = "grey", size = 1, expand=0.02) +  
    # geom_encircle(data = subset(plotdat, sample_type == sample_type2), color = "grey", size = 1, expand=0.02) +
    labs(x = "MDS_1",
         y = "MDS_2") +
    theme_classic() +
    geom_text_repel(data = plotdat, 
                    aes(x = x, y = y, label = sample), hjust= 0.30,
                    max.overlaps=100, min.segment.length = 0
    )
  # p3
  p4 <- ggplot(plotdat, aes(x = x, y = y, color = sample_type)) +
    geom_point() +  
    geom_mark_hull(aes(fill = sample_type)) +
    # geom_encircle(data = subset(plotdat, sample_type == sample_type1), color = "grey", size = 1, expand=0.02) +  
    # geom_encircle(data = subset(plotdat, sample_type == sample_type2), color = "grey", size = 1, expand=0.02) +
    labs(x = paste0("Leading logFC dim 1 (", as.character(as.integer(coords$var.explained[1]*100)), "%)"),
         y = paste0("Leading logFC dim 2 (", as.character(as.integer(coords$var.explained[2]*100)), "%)")) +
    theme_classic() +
    geom_text_repel(data = plotdat, 
                    aes(x = x, y = y, label = sample), hjust= 0.30,
                    max.overlaps=100, min.segment.length = 0
    )
  # p4
  # Set up design matrix:
  mm <- model.matrix(~0 + sample_type, data = plotdat)
  # Fit model to each protein:
  fit <- lmFit(dat.filtered, mm)
  # Specify comparison (“contrast”) being tested:
  # print(cond1vs2)
  contr <- makeContrasts(contrasts= cond1vs2, 
                         levels = cond_levels)
  contr2 <- contrasts.fit(fit, contrasts = contr)
  # Empirical Bayes smoothing adjusts the standard error of each protein to be closer to the average standard error, compensating for unusually high or low variability proteins.
  # Empirical Bayes smoothing:
  contr3 <- eBayes(contr2)
  results <- topTable(contr3, sort.by = "P", n = Inf)
  return(list(p1 = p1,
              p3=p3,
              p4=p4,
              deg = results))
}
