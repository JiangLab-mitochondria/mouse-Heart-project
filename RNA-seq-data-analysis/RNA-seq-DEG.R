library(dplyr)
library(ggplot2)
library(ggalt)
library(ggrepel)
library(edgeR)
set.seed(1234)


## define a function
trne_vs_wt_degs <- function(trne_samples_vector, wt_samples_vector) {
  df <- subset(df1, select = c("Geneid", trne_samples_vector, wt_samples_vector))
  rownames(df) <- df$Geneid
  df <- df[, -1]
  df <- as.matrix(df)
  group <- factor(c(
    "TrnE", "TrnE", "TrnE", "TrnE",
    "WT",  "WT",  "WT",  "WT"
  ))
  group <- relevel(group, ref = "WT")
  y <- DGEList(counts = df, group = group)
  print(y$samples)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  print(y$samples)
  results <- plotMDS(y, gene.selection="common", main = "PCA plot", plot=F)
  data <- data.frame(
    PC1 = results$x,
    PC2 = results$y,
    label = rownames(results$distance.matrix.squared),
    group = c("TrnE", "TrnE", "TrnE", "TrnE",
              "WT",  "WT",  "WT",  "WT")
  )
  print(data)
  pc1_value <- as.character(as.integer(results$var.explained[1] * 100))  
  pc2_value <- as.character(as.integer(results$var.explained[2] * 100))
  p1 <- ggplot(data, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +  
    geom_encircle(data = subset(data, group == "TrnE"), color = "grey", size = 1, expand=0.02) +  
    geom_encircle(data = subset(data, group == "WT"), color = "grey", size = 1, expand=0.02) +
    labs(x = paste0("PC1 (", pc1_value, "%)"),
         y = paste0("PC2 (", pc2_value, "%)")) +
    theme_classic() +
    geom_text_repel(data = data, 
                    aes(x = PC1, y = PC2, label = label), hjust= 0.30,
                    max.overlaps=100, min.segment.length = 0
    )
  
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  DEG <- as.data.frame(topTags(qlf, n=200000))
  DEG$gene_symbol <- rownames(DEG)
  DEG <- DEG %>%
    mutate(gene_type1 = case_when(logFC >= log2(1.2) & PValue <= 0.05 ~ "up",
                                  logFC <= -log2(1.2) & PValue <= 0.05 ~ "down",
                                  TRUE ~ "ns")) 
  print(table(DEG$gene_type1))
  return(list(fig = p1, deg = DEG))
}

nd5_vs_wt_degs <- function(nd5_samples_vector, wt_samples_vector) {
  df <- subset(df1, select = c("Geneid", nd5_samples_vector, wt_samples_vector))
  rownames(df) <- df$Geneid
  df <- df[, -1]
  df <- as.matrix(df)
  group <- factor(c(
    "ND5", "ND5", "ND5", "ND5",
    "WT",  "WT",  "WT",  "WT"
  ))
  group <- relevel(group, ref = "WT")
  y <- DGEList(counts = df, group = group)
  print(y$samples)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  print(y$samples)
  results <- plotMDS(y, gene.selection="common", main = "PCA plot", plot=F)
  data <- data.frame(
    PC1 = results$x,
    PC2 = results$y,
    label = rownames(results$distance.matrix.squared),
    group = c("ND5", "ND5", "ND5", "ND5",
              "WT",  "WT",  "WT",  "WT")
  )
  print(data)
  pc1_value <- as.character(as.integer(results$var.explained[1] * 100))  
  pc2_value <- as.character(as.integer(results$var.explained[2] * 100))
  p1 <- ggplot(data, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +  
    geom_encircle(data = subset(data, group == "ND5"), color = "grey", size = 1, expand=0.02) +  
    geom_encircle(data = subset(data, group == "WT"), color = "grey", size = 1, expand=0.02) +
    labs(x = paste0("PC1 (", pc1_value, "%)"),
         y = paste0("PC2 (", pc2_value, "%)")) +
    theme_classic() +
    geom_text_repel(data = data, 
                    aes(x = PC1, y = PC2, label = label), hjust= 0.30,
                    max.overlaps=100, min.segment.length = 0
    )
  
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  DEG <- as.data.frame(topTags(qlf, n=200000))
  DEG$gene_symbol <- rownames(DEG)
  DEG <- DEG %>%
    mutate(gene_type1 = case_when(logFC >= log2(1.2) & PValue <= 0.05 ~ "up",
                                  logFC <= -log2(1.2) & PValue <= 0.05 ~ "down",
                                  TRUE ~ "ns")) 
  print(table(DEG$gene_type1))
  return(list(fig = p1, deg = DEG))
}

trne_vs_nd5_degs <- function(trne_samples_vector, nd5_samples_vector) {
  df <- subset(df1, select = c("Geneid", trne_samples_vector, nd5_samples_vector))
  rownames(df) <- df$Geneid
  df <- df[, -1]
  df <- as.matrix(df)
  group <- factor(c(
    "TrnE", "TrnE", "TrnE", "TrnE",
    "ND5",  "ND5",  "ND5",  "ND5"
  ))
  group <- relevel(group, ref = "ND5")
  y <- DGEList(counts = df, group = group)
  print(y$samples)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  print(y$samples)
  results <- plotMDS(y, gene.selection="common", main = "PCA plot", plot=F)
  data <- data.frame(
    PC1 = results$x,
    PC2 = results$y,
    label = rownames(results$distance.matrix.squared),
    group = c("TrnE", "TrnE", "TrnE", "TrnE",
              "ND5",  "ND5",  "ND5",  "ND5")
  )
  print(data)
  pc1_value <- as.character(as.integer(results$var.explained[1] * 100))  
  pc2_value <- as.character(as.integer(results$var.explained[2] * 100))
  p1 <- ggplot(data, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +  
    geom_encircle(data = subset(data, group == "TrnE"), color = "grey", size = 1, expand=0.02) +  
    geom_encircle(data = subset(data, group == "ND5"), color = "grey", size = 1, expand=0.02) +
    labs(x = paste0("PC1 (", pc1_value, "%)"),
         y = paste0("PC2 (", pc2_value, "%)")) +
    theme_classic() +
    geom_text_repel(data = data, 
                    aes(x = PC1, y = PC2, label = label), hjust= 0.30,
                    max.overlaps=100, min.segment.length = 0
    )
  
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  DEG <- as.data.frame(topTags(qlf, n=200000))
  DEG$gene_symbol <- rownames(DEG)
  DEG <- DEG %>%
    mutate(gene_type1 = case_when(logFC >= log2(1.2) & PValue <= 0.05 ~ "up",
                                  logFC <= -log2(1.2) & PValue <= 0.05 ~ "down",
                                  TRUE ~ "ns")) 
  print(table(DEG$gene_type1))
  return(list(fig = p1, deg = DEG))
}
