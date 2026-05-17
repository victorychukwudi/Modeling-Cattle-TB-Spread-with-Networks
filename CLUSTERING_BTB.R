

options(bitmapType="cairo")

setwd("C:/Users/Victory/Downloads")
rm(list = ls())

library(tidyverse)
library(lubridate)
library(matrixStats)
library(cluster)
library(factoextra)
library(corrplot)
library(gridExtra)
library(NbClust)
library(clustertend)

bd_df <- read.csv("bd_df_04_Feb_2026_encrypted.csv")

############################################################
# 01a. CORE FILTERING
############################################################

bd_victory <- bd_df %>%
  arrange(herd_no, bd_no) %>%
  filter(all_cases != 0) %>%
  filter(bd_duration_days >= 60) %>%
  filter(feedlot_cfu_ever == 0) %>%
  filter(
    bd_initiated == "Slaughter (non permit animal)" |
      (bd_initiated == "SICTT" &
         bd_first_skin_test_type %in% c("1", "8"))
  )

############################################################
# 01b. BASIC VARIABLES
############################################################

bd_victory <- bd_victory %>%
  mutate(
    bd_start = as.Date(bd_start),

    bd_within_13_months =
      case_when(
        is.na(duration_between_bd) ~ 0,
        duration_between_bd < 390 ~ 1,
        TRUE ~ 0
      ),

    herd_size_mean =
      rowMeans(select(., mean_herd_size_jan,
                      mean_herd_size_may,
                      mean_herd_size_sep))
  )

############################################################
# 01c. HISTORY VARIABLE
############################################################

bd_victory <- bd_victory %>%
  group_by(herd_no) %>%
  arrange(bd_start, .by_group = TRUE) %>%
  mutate(
    no_bds_started_in_prev_5yr =
      sapply(
        bd_start,
        function(x)
          sum(
            bd_start < x &
              bd_start >= x %m-% years(5),
            na.rm = TRUE
          )
      )
  ) %>%
  ungroup()

############################################################
# 01d. CORE EPIDEMIOLOGY
############################################################

bd_victory <- bd_victory %>%
  mutate(

    total_index_reactors =
      cows_positive_index_test +
      bulls_positive_index_test +
      steers_positive_index_test +
      heifers_positive_index_test +
      calves_positive_index_test,

    total_tested_index = no_animals_index_test,

    total_subsequent_reactors =
      (cows_positive + bulls_positive + steers_positive +
         heifers_positive + calves_positive) -
      total_index_reactors,

    prop_positive_index =
      ifelse(total_tested_index > 0,
             total_index_reactors / total_tested_index,
             NA),

    spread_ratio =
      ifelse(total_index_reactors > 0,
             total_subsequent_reactors / total_index_reactors,
             NA),

    spread_any =
      ifelse(total_subsequent_reactors > 0, 1, 0)
  ) %>%
  mutate(
    spread_ratio = ifelse(is.infinite(spread_ratio), NA, spread_ratio)
  )

############################################################
# 02a. SELECT FINAL VARIABLES
############################################################

core_vars <- c(
  "prop_positive_index",
  "spread_ratio",
  "spread_any",
  "bd_duration_days",
  "bd_within_13_months",
  "no_bds_started_in_prev_5yr",
  "herd_size_mean"
)

bd_core <- bd_victory %>%
  select(all_of(core_vars)) %>%
  drop_na()


# 02b. LOG TRANSFORMATION

bd_core_log <- bd_core %>%
  mutate(
    spread_ratio = log1p(spread_ratio),
    bd_duration_days = log1p(bd_duration_days),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )


# 02c. SCALING 

bd_core_scaled <- scale(bd_core_log)

#before exploration just a summary
library(dplyr)

desc_table <- bd_core %>%
  summarise(across(everything(),
                   list(
                     Mean = ~mean(. , na.rm = TRUE),
                     SD   = ~sd(. , na.rm = TRUE),
                     Min  = ~min(. , na.rm = TRUE),
                     Max  = ~max(. , na.rm = TRUE)
                   ))) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = c("Variable", ".value"),
                      names_pattern = "(.*)_(Mean|SD|Min|Max)")

print(desc_table)




############################################################
# ================= EXPLORATION =================
############################################################

############################################################
# 03a. HISTOGRAMS (RAW)
############################################################

png("01a_hist_raw.png",1600,1200)

bd_core %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill="steelblue") +
  facet_wrap(~name, scales="free") +
  theme_bw()

dev.off()


# 03b. HISTOGRAMS (LOG)


png("01b_hist_log.png",1600,1200)

bd_core_log %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill="darkgreen") +
  facet_wrap(~name, scales="free") +
  theme_bw()

dev.off()

############################################################
# 03c. DENSITY
############################################################

png("01c_density.png",1600,1200)

bd_core_log %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value)) +
  geom_density(fill="orange", alpha=0.5) +
  facet_wrap(~name, scales="free") +
  theme_bw()

dev.off()

############################################################
# 03d. SCATTER MATRIX
############################################################

png("01d_scatter_matrix.png",1600,1400)

pairs(as.data.frame(bd_core_log), pch=20, cex=0.4)

dev.off()


# 03e. KEY SCATTER 


png("01e_scatter_spread_vs_severity.png",1200,900)

plot(
  bd_core$prop_positive_index,
  log1p(bd_core$spread_ratio),
  pch=19,
  col=rgb(0,0,1,0.05),
  xlab="Severity",
  ylab="Log Spread"
)

dev.off()

############################################################
# CORRELATION PLOT 
############################################################

png("01f_correlation_full_with_numbers.png",1200,1000)

corr_matrix <- cor(bd_core_log, use = "pairwise.complete.obs")

corrplot(
  corr_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  number.cex = 0.5,
  tl.cex = 0.6
)

dev.off()

############################################################
# 03g. BOXPLOTS
############################################################

png("01g_boxplots.png",1600,1200)

bd_core_log %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x=name, y=value)) +
  geom_boxplot(fill="grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

dev.off()



############################################################
# 07a. SAMPLE Its better after exploration
############################################################

set.seed(123)

sample_size <- min(20000, nrow(bd_core_scaled))

bd_core_sample <- bd_core_scaled[
  sample(nrow(bd_core_scaled), sample_size),
]


# ================= PCA ANALYSIS =================


library(factoextra)
library(ggplot2)
library(gridExtra)
library(corrplot)

bd_pca_df <- as.data.frame(bd_core_sample)


# 07b. DISTRIBUTION CHECK


png("04a_pca_hist.png",1600,1200)

bd_pca_df %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ name, scales = "free") +
  theme_bw()

dev.off()


# 07c. CORRELATION MATRIX


png("04b_pca_correlation.png",1200,1000)

corrplot(
  cor(bd_pca_df, use = "pairwise.complete.obs"),
  method = "color",
  type = "upper",
  order = "hclust"
)

dev.off()


# 07d. RUN PCA 

bd_pca <- prcomp(
  bd_pca_df,
  center = FALSE,
  scale. = FALSE
)


res.pca <- prcomp(bd_core_sample, scale = TRUE)

print(res.pca)
summary(res.pca)

# 07e. VARIANCE EXPLAINED

png("04c_variance_explained.png",1200,900)

fviz_eig(
  bd_pca,
  addlabels = TRUE,
  ylim = c(0, 60)
)

dev.off()

# 07f. PCA INDIVIDUALS 
png("04d_pca_individuals.png",1200,900)

fviz_pca_ind(
  bd_pca,
  geom = "point",
  pointsize = 1,
  alpha.ind = 0.4,
  ggtheme = theme_bw()
)

dev.off()


# ================= CLUSTERING ON PCA =================



# 08a. EXTRACT PCA SCORES

bd_pca_scores <- as.data.frame(bd_pca$x)

# OPTIONAL: 
bd_pca_reduced <- bd_pca_scores[, 1:3]

############################################################
# PCA CLUSTER VISUALS (K = 2–5)
############################################################

png("pca_kmeans_clusters.png", width = 1600, height = 1200)

p2 <- fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(bd_K2$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
) + ggtitle("K = 2")

p3 <- fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(bd_K3$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
) + ggtitle("K = 3")

p4 <- fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(bd_K4$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
) + ggtitle("K = 4")

p5 <- fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(kmeans(bd_pca_reduced, centers = 5, nstart = 25)$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
) + ggtitle("K = 5")

gridExtra::grid.arrange(p2, p3, p4, p5, ncol = 2)

dev.off()

############################################################
# 08c. PCA CLUSTER VISUALS
############################################################

png("05a_pca_clusters_k2.png",1200,900)

fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(bd_K2$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
)

dev.off()

png("05b_pca_clusters_k3.png",1200,900)

fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(bd_K3$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
)

dev.off()

png("05c_pca_clusters_k4.png",1200,900)

fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = factor(bd_K4$cluster),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
)

dev.off()

############################################################
# 08d. VARIABLE CONTRIBUTIONS
############################################################

png("05d_variable_contribution_PC1.png",1200,900)
fviz_contrib(bd_pca, choice = "var", axes = 1)
dev.off()

png("05e_variable_contribution_PC2.png",1200,900)
fviz_contrib(bd_pca, choice = "var", axes = 2)
dev.off()

############################################################
# 08e. BIPLOT
############################################################

library(ggfortify)

png("05f_pca_biplot.png",1200,900)

autoplot(
  bd_pca,
  loadings = TRUE,
  loadings.colour = "darkblue",
  loadings.label = TRUE,
  loadings.label.size = 3
)

dev.off()




############################################################
# 07g. OPTIMAL NUMBER OF PCA COMPONENTS
############################################################

# Run PCA (ensure consistency)
res.pca <- prcomp(bd_core_sample, center = TRUE, scale. = TRUE)

# Extract variance
pca_var <- summary(res.pca)$importance

pca_var_df <- data.frame(
  PC = paste0("PC", 1:ncol(pca_var)),
  Variance_Explained = pca_var[2, ],
  Cumulative_Variance = pca_var[3, ]
)

print(pca_var_df)

# Save
write.csv(pca_var_df, "04c_pca_variance_table.csv", row.names = FALSE)

# Select optimal PCs (≥80% variance)
threshold <- 0.80
optimal_pcs <- which(pca_var_df$Cumulative_Variance >= threshold)[1]

cat("Optimal number of PCs:", optimal_pcs, "\n")


# Scree plot with selection
png("04c_variance_explained.png",1200,900)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 60)) +
  geom_vline(xintercept = optimal_pcs, linetype = "dashed", color = "red") +
  ggtitle("Scree Plot with Optimal Components")

dev.off()


############################################################
# 07h. PCA COMPONENT INTERPRETATION
############################################################

# Extract loadings
loadings <- as.data.frame(res.pca$rotation)

print(round(loadings, 3))

write.csv(loadings, "04d_pca_loadings.csv")

# Optional: contributions plots
png("04d_PC1_contributions.png",1200,900)
fviz_contrib(res.pca, choice = "var", axes = 1)
dev.off()

png("04e_PC2_contributions.png",1200,900)
fviz_contrib(res.pca, choice = "var", axes = 2)
dev.off()

png("04f_PC3_contributions.png",1200,900)
fviz_contrib(res.pca, choice = "var", axes = 3)
dev.off()




############################################################
# 07h. PCA COMPONENT INTERPRETATION
############################################################

# Extract loadings
loadings <- as.data.frame(res.pca$rotation)

print(round(loadings, 3))

write.csv(loadings, "04d_pca_loadings.csv")

# Optional: contributions plots
png("04d_PC1_contributions.png",1200,900)
fviz_contrib(res.pca, choice = "var", axes = 1)
dev.off()

png("04e_PC2_contributions.png",1200,900)
fviz_contrib(res.pca, choice = "var", axes = 2)
dev.off()

png("04f_PC3_contributions.png",1200,900)
fviz_contrib(res.pca, choice = "var", axes = 3)
dev.off()





# ================= CLUSTER TENDENCY =================


#set.seed(123)

#sample_size <- min(20000, nrow(bd_core_scaled))

#bd_core_sample <- bd_core_scaled[
 # sample(nrow(bd_core_scaled), sample_size),
#]

#ToCluster <- as.matrix(bd_core_sample)

############################################################
# 04a. HOPKINS
############################################################

hopkins_stat <- hopkins(
  ToCluster,
  n = floor(0.1 * nrow(ToCluster))
)

hopkins_stat$H

############################################################
# 04b. VAT
############################################################

set.seed(123)

vat_index <- sample(1:nrow(ToCluster), 1000)
vat_data <- ToCluster[vat_index, ]

png("02a_vat_plot.png",1200,1000)

fviz_dist(dist(vat_data))

dev.off()


# ================= CLUSTER NUMBER =================

# 05a. WSS

png("02b_wss.png",1200,900)

fviz_nbclust(bd_core_sample, kmeans, method="wss", k.max=10)

dev.off()

##ELBOW METHOD using variance 
library(ClusterR)
set.seed(123)

opt <- Optimal_Clusters_KMeans(
  bd_core_sample,
  max_clusters = 10,
  plot_clusters = TRUE
)
dev.off()


# 05a ELBOW (ClusterR)

png("elbow07v_plot.png",height=800,width=600)
library(ClusterR)

set.seed(123)

opt <- Optimal_Clusters_KMeans(
  bd_core_sample,
  max_clusters = 10,
  plot_clusters = TRUE
)
dev.off()

#for curves not ellipses
png("02bc_wss_curve.png", 1200, 900)

fviz_nbclust(bd_core_sample, kmeans, method = "wss", k.max = 10) +
  ggtitle("Elbow Method (WSS vs Number of Clusters)")

dev.off()

# 05b. SILHOUETTE METHOD


png("02c_silhouette_method.png",1200,900)

fviz_nbclust(bd_core_sample, kmeans, method="silhouette", k.max=10)

dev.off()


# 05c. GAP STAT 


png("02d_gap_stat.png",1200,900)

fviz_nbclust(
  bd_core_sample,
  kmeans,
  method="gap_stat",
  k.max=10,
  nboot=10   
)

dev.off()


# 05d. NBCLUST 
set.seed(123)

nb_all <- NbClust(
  data = bd_core_sample,
  distance = "euclidean",
  min.nc = 2,
  max.nc = 10,
  method = "ward.D2"
)

png("02e_nbclust_barplot.png",1200,900)

barplot(table(nb_all$Best.nc[1, ]),
        col="grey",
        main="Optimal K")

dev.off()




#Hierachical clustering

# Distance matrix
dist_mat <- dist(bd_core_sample, method = "euclidean")

# Hierarchical clustering (choose linkage)
hc <- hclust(dist_mat, method = "single")  # or "complete", "average"

# Plot dendrogram
plot(hc, main = "Hierarchical Clustering Dendrogram (Euclidean Distance)")

library(factoextra)

fviz_dend(hc, 
          k = 3,              # show 3 clusters
          rect = TRUE,        # draw boxes
          main = "Dendrogram (Euclidean Distance)")

# ================= K COMPARISON =================

# 06a. K-MEANS MODELS 

set.seed(123)

bd_clusters2 <- kmeans(bd_core_sample, 2, nstart=25)
bd_clusters3 <- kmeans(bd_core_sample, 3, nstart=25)
bd_clusters4 <- kmeans(bd_core_sample, 4, nstart=25)
bd_clusters5 <- kmeans(bd_core_sample, 5, nstart=25)



############################################################
# 06d. K-MEANS CLUSTER MEANS (7 VARIABLES)
############################################################

bd_sample_df <- as.data.frame(bd_core_sample)

# K = 2
k2_means <- bd_sample_df %>%
  mutate(cluster = bd_clusters2$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(k2_means, "07v_kmeans_k2.csv", row.names = FALSE)

# K = 3
k3_means <- bd_sample_df %>%
  mutate(cluster = bd_clusters3$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(k3_means, "07v_kmeans_k3.csv", row.names = FALSE)

# K = 4
k4_means <- bd_sample_df %>%
  mutate(cluster = bd_clusters4$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(k4_means, "07v_kmeans_k4.csv", row.names = FALSE)

# K = 5
k5_means <- bd_sample_df %>%
  mutate(cluster = bd_clusters5$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(k5_means, "07v_kmeans_k5.csv", row.names = FALSE)


# 06b. CLUSTER VISUALS 
fviz_cluster(
  list(data = bd_core_sample, cluster = bd_clusters3$cluster),
  geom = "point",
  ellipse.type = "norm",
  ellipse.level = 0.95,
  repel = TRUE
)


png("03a_cluster_plots.png", width = 1600, height = 1200)

p2 <- fviz_cluster(
  list(data = bd_core_sample, cluster = bd_clusters2$cluster),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("K = 2")

p3 <- fviz_cluster(
  list(data = bd_core_sample, cluster = bd_clusters3$cluster),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("K = 3")

p4 <- fviz_cluster(
  list(data = bd_core_sample, cluster = bd_clusters4$cluster),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("K = 4")

p5 <- fviz_cluster(
  list(data = bd_core_sample, cluster = bd_clusters5$cluster),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("K = 5")

gridExtra::grid.arrange(p2, p3, p4, p5, ncol = 2)

dev.off()


############################################################
# K = 3 CLUSTER MEANS ONLY
############################################################

bd_sample_df <- as.data.frame(bd_core_sample)

k3_means <- bd_sample_df %>%
  mutate(cluster = bd_clusters3$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(k3_means, "kmeans_k3_means.csv", row.names = FALSE)

############################################################
# K = 3 VISUAL ONLY (NO GRID)
############################################################

png("kmeans_k3_plot.png", width = 1200, height = 900)

fviz_cluster(
  list(data = bd_core_sample, cluster = bd_clusters3$cluster),
  geom = "point",
  ellipse.type = "norm",
  ellipse.level = 0.95,
  repel = TRUE
) + ggtitle("K = 3")

dev.off()








# 06c. SILHOUETTE COMPARISON FOR KMEANS
library(cluster)

dist_sample <- dist(bd_core_sample)

sil2 <- silhouette(bd_clusters2$cluster, dist_sample)
sil3 <- silhouette(bd_clusters3$cluster, dist_sample)
sil4 <- silhouette(bd_clusters4$cluster, dist_sample)
sil5 <- silhouette(bd_clusters5$cluster, dist_sample)

s2 <- fviz_silhouette(sil2)
s3 <- fviz_silhouette(sil3)
s4 <- fviz_silhouette(sil4)
s5 <- fviz_silhouette(sil5)

png("03b_silhouette_plots.png",1600,1200)
grid.arrange(s2,s3,s4,s5,ncol=2)
dev.off()

# 06c. SILHOUETTE (K = 3 ONLY) FOR KMEANS

library(cluster)
library(factoextra)

dist_sample <- dist(bd_core_sample)

sil3 <- silhouette(bd_clusters3$cluster, dist_sample)

p_sil3 <- fviz_silhouette(sil3) + ggtitle("Silhouette Plot (K = 3)")

png("03b_silhouette_k3.png", 800, 600)
print(p_sil3)
dev.off()












##PAM

############################################################
# 06e. PAM CLUSTERING (7 VARIABLES)
############################################################

library(cluster)

# K = 2
pam2 <- pam(bd_core_sample, k = 2)
pam2_medoids <- as.data.frame(pam2$medoids)
write.csv(pam2_medoids, "07v_pam_k2_medoids.csv", row.names = FALSE)

# K = 3
pam3 <- pam(bd_core_sample, k = 3)
pam3_medoids <- as.data.frame(pam3$medoids)
write.csv(pam3_medoids, "07v_pam_k3_medoids.csv", row.names = FALSE)

# K = 4
pam4 <- pam(bd_core_sample, k = 4)
pam4_medoids <- as.data.frame(pam4$medoids)
write.csv(pam4_medoids, "07v_pam_k4_medoids.csv", row.names = FALSE)

# K = 5
pam5 <- pam(bd_core_sample, k = 5)
pam5_medoids <- as.data.frame(pam5$medoids)
write.csv(pam5_medoids, "07v_pam_k5_medoids.csv", row.names = FALSE)



# PAM cluster means (for comparison)

pam3_means <- bd_core_sample %>%
  mutate(cluster = pam3$clustering) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(pam3_means, "07v_pam_k3_means.csv", row.names = FALSE)


############################################################
# 06f. PAM CLUSTER VISUALS
############################################################

p_pam2 <- fviz_cluster(list(data = bd_core_sample, cluster = pam2$clustering)) + ggtitle("PAM K=2")
p_pam3 <- fviz_cluster(list(data = bd_core_sample, cluster = pam3$clustering)) + ggtitle("PAM K=3")
p_pam4 <- fviz_cluster(list(data = bd_core_sample, cluster = pam4$clustering)) + ggtitle("PAM K=4")
p_pam5 <- fviz_cluster(list(data = bd_core_sample, cluster = pam5$clustering)) + ggtitle("PAM K=5")



p_pam2 <- fviz_cluster(
  list(data = bd_core_sample, cluster = pam2$clustering),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("PAM K = 2")

p_pam3 <- fviz_cluster(
  list(data = bd_core_sample, cluster = pam3$clustering),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("PAM K = 3")

p_pam4 <- fviz_cluster(
  list(data = bd_core_sample, cluster = pam4$clustering),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("PAM K = 4")

p_pam5 <- fviz_cluster(
  list(data = bd_core_sample, cluster = pam5$clustering),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("PAM K = 5")

png("07v_pam_cluster_plots.png", width = 1600, height = 1200)

gridExtra::grid.arrange(p_pam2, p_pam3, p_pam4, p_pam5, ncol = 2)

dev.off()




# PAM CLUSTERING (K = 3)


library(cluster)
library(dplyr)
library(factoextra)

# Run PAM with K = 3
pam3 <- pam(bd_core_sample, k = 3)

# Save medoids
pam3_medoids <- as.data.frame(pam3$medoids)
write.csv(pam3_medoids, "07v_pam_k3_medoids.csv", row.names = FALSE)

# Cluster means (standardized variables)
pam3_means <- bd_core_sample %>%
  mutate(cluster = pam3$clustering) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(pam3_means, "07v_pam_k3_means.csv", row.names = FALSE)


############################################################
# PAM VISUAL (K = 3 ONLY)
############################################################

p_pam3 <- fviz_cluster(
  list(data = bd_core_sample, cluster = pam3$clustering),
  geom = "point",
  ellipse.type = "norm"
) + ggtitle("PAM (K = 3)")

# Save plot
png("07v_pam_k3_plot.png", width = 800, height = 600)
print(p_pam3)
dev.off()












############################################################
# 06g. PAM SILHOUETTE ANALYSIS
############################################################

library(cluster)

dist_sample <- dist(bd_core_sample)

pam_sil2 <- silhouette(pam2$clustering, dist_sample)
pam_sil3 <- silhouette(pam3$clustering, dist_sample)
pam_sil4 <- silhouette(pam4$clustering, dist_sample)
pam_sil5 <- silhouette(pam5$clustering, dist_sample)

sp2 <- fviz_silhouette(pam_sil2) + ggtitle("PAM K=2")
sp3 <- fviz_silhouette(pam_sil3) + ggtitle("PAM K=3")
sp4 <- fviz_silhouette(pam_sil4) + ggtitle("PAM K=4")
sp5 <- fviz_silhouette(pam_sil5) + ggtitle("PAM K=5")

png("07v_pam_silhouette.png", 1600, 1200)
grid.arrange(sp2, sp3, sp4, sp5, ncol = 2)
dev.off()


############################################################
# 06h. K-MEANS vs PAM SILHOUETTE COMPARISON (07 VARIABLES)
############################################################

# Function to extract average silhouette width
avg_sil <- function(sil_obj) {
  mean(sil_obj[, 3])
}

# K-means average silhouette
kmeans_sil_scores <- data.frame(
  Method = "KMEANS",
  K = c(2, 3, 4, 5),
  Avg_Silhouette = c(
    avg_sil(sil2),
    avg_sil(sil3),
    avg_sil(sil4),
    avg_sil(sil5)
  )
)

# PAM average silhouette
pam_sil_scores <- data.frame(
  Method = "PAM",
  K = c(2, 3, 4, 5),
  Avg_Silhouette = c(
    avg_sil(pam_sil2),
    avg_sil(pam_sil3),
    avg_sil(pam_sil4),
    avg_sil(pam_sil5)
  )
)

# Combine results
sil_comparison <- rbind(kmeans_sil_scores, pam_sil_scores)

# Save CSV
write.csv(sil_comparison, "07v_kmeans_vs_pam_silhouette.csv", row.names = FALSE)




# 06g. PAM SILHOUETTE (K = 3 ONLY)

library(cluster)
library(factoextra)

dist_sample <- dist(bd_core_sample)

pam_sil3 <- silhouette(pam3$clustering, dist_sample)

p_sil3 <- fviz_silhouette(pam_sil3) + ggtitle("PAM K = 3")

png("07v_pam_silhouette_k3.png", 800, 600)
print(p_sil3)
dev.off()


############################################################
# 06i. SILHOUETTE COMPARISON PLOT
############################################################

png("07v_kmeans_vs_pam_silhouette_plot.png", 1200, 900)

ggplot(sil_comparison, aes(x = factor(K), y = Avg_Silhouette, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(
    title = "K-means vs PAM Silhouette Comparison (7 Variables)",
    x = "Number of Clusters (K)",
    y = "Average Silhouette Width"
  )

dev.off()

#PAM HERD_CLUSTER CSV



library(dplyr)
library(cluster)

set.seed(123)

sample_size <- min(20000, nrow(bd_core_scaled))

bd_core_sample <- bd_core_scaled[
  sample(nrow(bd_core_scaled), sample_size),
]

bd_core_sample <- as.data.frame(bd_core_sample)

# Run PAM
pam3 <- pam(bd_core_sample, k = 3)

# Attach clusters
bd_core_sample_clustered <- bd_core_sample %>%
  mutate(cluster = pam3$clustering)

# Save
write.csv(
  bd_core_sample_clustered,
  "C:/Users/Victory/Downloads/07v_pam_k3_cluster_assignments.csv",
  row.names = FALSE
)



library(dplyr)
multi_data <- read.csv("Data_for_MN_regression_6-5-26_encrypted.csv")
pam_full <- read.csv("pam_k3_herd_clusters_sampled.csv")
multi_data$herd_no <- as.character(multi_data$herd_no)
pam_full$herd_no <- as.character(pam_full$herd_no)

multi_data_joined <- multi_data %>%
  left_join(pam_full, by = "herd_no")

analysis_data <- multi_data_joined %>%
  filter(!is.na(cluster))

print(dim(analysis_data))
print(table(analysis_data$cluster))

# Save final dataset
write.csv(
  analysis_data,
  "C:/Users/Victory/Downloads/multi_data_with_clusters_FINAL.csv",
  row.names = FALSE
)




library(dplyr)
library(cluster)

set.seed(123)

bd_clustered <- bd_victory %>%
  select(herd_no, bd_no, all_of(core_vars)) %>%
  drop_na() %>%
  sample_n(min(20000, n())) %>% 
  mutate(
    spread_ratio = log1p(spread_ratio),
    bd_duration_days = log1p(bd_duration_days),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )

scaled_data <- bd_clustered %>%
  select(-herd_no, -bd_no) %>%
  scale() %>%
  as.data.frame()

# PAM clustering
pam_model <- pam(scaled_data, k = 3)

# Attach clusters 
bd_clustered$cluster <- pam_model$clustering

final_clusters <- bd_clustered %>%
  select(herd_no, bd_no, cluster)

# Save
write.csv(
  final_clusters,
  "C:/Users/Victory/Downloads/pam_k3_herd_bd_clusters_sampled.csv",
  row.names = FALSE
)








library(mclust)
set.seed(123)

bd_mclust <- Mclust(bd_core_sample)
summary(bd_mclust)
png("04a_mclust_BIC.png",1200,900)

plot(bd_mclust, what = "BIC")

dev.off()

png("04b_mclust_classification.png",1200,900)

plot(bd_mclust, what = "classification")

dev.off()


bd_mclust_clusters <- bd_mclust$classification


head(bd_mclust_clusters)
bd_mclust_prob <- bd_mclust$z
head(bd_mclust_prob)

bd_uncertainty <- bd_mclust$uncertainty

summary(bd_uncertainty)

bd_df <- as.data.frame(bd_core_sample)

gmm_means <- bd_df %>%
  mutate(cluster = bd_mclust_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

write.csv(gmm_means, "07v_gmm_means.csv", row.names = FALSE)
############################################################
# 09b. GMM UNCERTAINTY EXPORT
############################################################

gmm_uncertainty_df <- data.frame(
  cluster = bd_mclust_clusters,
  uncertainty = bd_uncertainty
)

write.csv(gmm_uncertainty_df, "07v_gmm_uncertainty.csv", row.names = FALSE)


############################################################
# 09c. GMM CLASSIFICATION + UNCERTAINTY PLOT (SAVE)
############################################################

png("07v_gmm_uncertainty_plot.png", 1200, 900)

plot(bd_mclust, what = "uncertainty")

dev.off()



############################################################
# 10a. POST-CLUSTER PCA (K = 3 EXAMPLE)
############################################################

bd_post <- bd_sample_df %>%
  mutate(cluster = factor(bd_clusters3$cluster))

png("post_cluster_PCA_k3.png", 1200, 900)

fviz_pca_ind(
  bd_pca,
  geom = "point",
  habillage = bd_post$cluster,
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 1,
  alpha.ind = 0.5
)

dev.off()


############################################################
# 10b. POST-CLUSTER DISTRIBUTIONS
############################################################

png("post_cluster_distributions.png", 1600, 1200)

bd_post %>%
  pivot_longer(-cluster) %>%
  ggplot(aes(value, fill = cluster)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~name, scales = "free") +
  theme_bw()

dev.off()


############################################################
# 10c. POST-CLUSTER SCATTER MATRIX
############################################################

png("post_cluster_scatter_matrix.png", 1600, 1400)

pairs(
  bd_post[, -ncol(bd_post)],
  col = bd_post$cluster,
  pch = 19,
  cex = 0.5
)

dev.off()



# Add cluster assignments to your dataset
pam3_clusters <- bd_core_sample %>%
  mutate(cluster = pam3$clustering)

# View herd number with cluster
pam3_clusters %>%
  select(herd_no, cluster)

# Optional: save to CSV
write.csv(pam3_clusters %>% select(herd_no, cluster),
          "07v_pam_k3_herd_clusters.csv",
          row.names = FALSE)


##proportion of each cluster per year 
############################################################
# 11a. Full data for pam no sampling lets try if it works
############################################################

bd_core_full <- bd_victory %>%
  select(
    bd_start,
    bd_start_yr,
    all_of(core_vars)
  ) %>%
  drop_na()

bd_core_log_full <- bd_core_full %>%
  mutate(
    spread_ratio = log1p(spread_ratio),
    bd_duration_days = log1p(bd_duration_days),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )

bd_scaled_full <- bd_core_log_full %>%
  select(all_of(core_vars)) %>%
  scale() %>%
  as.data.frame()

############################################################
# 11b. PAM K = 3 (FULL DATA)
############################################################

pam3_full <- pam(bd_scaled_full, k = 3)



############################################################
# 11c. ATTACH CLUSTERS
############################################################

bd_clustered <- bd_core_full %>%
  mutate(cluster = pam3_full$clustering)

############################################################
# 11d. PROPORTION PER YEAR
############################################################

cluster_year_prop <- bd_clustered %>%
  group_by(bd_start_yr, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bd_start_yr) %>%
  mutate(
    proportion = n / sum(n)
  )

print(cluster_year_prop)

write.csv(cluster_year_prop,
          "08_cluster_proportion_per_year.csv",
          row.names = FALSE)


############################################################
# 11e. PLOT
############################################################

png("08_cluster_proportion_per_year.png", 1200, 900)

ggplot(cluster_year_prop,
       aes(x = bd_start_yr,
           y = proportion,
           fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "Proportion of PAM (K = 3) Clusters by Year",
    x = "Year",
    y = "Proportion",
    fill = "Cluster"
  )

dev.off()


##what are mine testing?
##do herds in different cluster have different
# bd duration proportion overtime?

library(dplyr)
library(survival)
library(survminer)

# Recreating the SAME dataset used for clustering 
bd_analysis <- bd_victory %>%
  select(
    bd_duration_days,
    prop_positive_index,
    spread_ratio,
    spread_any,
    bd_within_13_months,
    no_bds_started_in_prev_5yr,
    herd_size_mean
  ) %>%
  drop_na()

# Apply SAME transformations as clustering
bd_analysis_log <- bd_analysis %>%
  mutate(
    spread_ratio = log1p(spread_ratio),
    bd_duration_days = log1p(bd_duration_days),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )

bd_analysis_scaled <- scale(bd_analysis_log)

# Run PAM on FULL dataset
pam3_full <- pam(bd_analysis_scaled, k = 3)

# Attach clusters
bd_analysis$cluster <- factor(pam3_full$clustering)

# Define survival object
# All breakdowns ended → event = 1
surv_obj <- Surv(
  time = bd_analysis$bd_duration_days,
  event = rep(1, nrow(bd_analysis))
)

#  Kaplan-Meier
km_fit <- survfit(surv_obj ~ cluster, data = bd_analysis)

#  Plot
ggsurvplot(
  km_fit,
  data = bd_analysis,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Cluster",
  xlab = "Breakdown duration (days)",
  ylab = "Probability of remaining in breakdown"
)



###NEXT BREAKDOWN
library(survival)
library(survminer)
library(dplyr)
library(cluster)

# 1. Build dataset 
bd_surv <- bd_victory %>%
  select(
    duration_between_bd,
    prop_positive_index,
    spread_ratio,
    spread_any,
    bd_within_13_months,
    no_bds_started_in_prev_5yr,
    herd_size_mean
  ) %>%
  drop_na(prop_positive_index, spread_ratio, herd_size_mean)

# 2. Transform variables (same as clustering)
bd_surv_log <- bd_surv %>%
  mutate(
    spread_ratio = log1p(spread_ratio),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )

# 3. Scale (exclude duration_between_bd)
bd_surv_scaled <- scale(bd_surv_log %>% select(-duration_between_bd))

# 4. PAM clustering
pam3_full <- pam(bd_surv_scaled, k = 3)

# 5. Attach clusters
bd_surv$cluster <- factor(pam3_full$clustering)

############################################################
# 6. DEFINE SURVIVAL VARIABLES 
############################################################

# Choose a follow-up limit (same as your plot)
max_followup <- 1000

time <- ifelse(
  is.na(bd_surv$duration_between_bd),
  max_followup,                      # censored at study limit
  bd_surv$duration_between_bd
)

event <- ifelse(
  is.na(bd_surv$duration_between_bd),
  0,                                 # no next breakdown → censored
  1                                  # next breakdown occurred
)

# Survival object
surv_obj <- Surv(time = time, event = event)

# 7. Kaplan-Meier fit
km_fit <- survfit(surv_obj ~ cluster, data = bd_surv)

# 8. Plot
ggsurvplot(
  km_fit,
  data = bd_surv,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Cluster",
  xlab = "Time to next breakdown (days)",
  ylab = "Probability of remaining breakdown-free",
  xlim = c(0, 1000),
  break.time.by = 100
)



############################################################
# Hierarchical clustering 
############################################################

library(factoextra)

# Distance matrix (use scaled data ideally)
dist_mat <- dist(bd_core_sample, method = "euclidean")

# Hierarchical clustering (Ward = recommended)
hc <- hclust(dist_mat, method = "ward.D2")

# Plot dendrogram (clean)
fviz_dend(hc, 
          k = 3,                 # number of clusters
          rect = TRUE,           # draw cluster boxes
          show_labels = FALSE,   # remove messy labels
          main = "Hierarchical Clustering Dendrogram (Ward's Method)")



##MULTINOMIAL REGRESSION

#What characteristics/attributes drives this herds to be 
#in each cluster label?
setwd("C:/Users/Victory/Downloads")
rm(list = ls())

library(tidyverse)
library(lubridate)
library(matrixStats)
library(cluster)
library(factoextra)
library(corrplot)
library(gridExtra)
library(NbClust)
library(clustertend)

bd_df <- read.csv("bd_df_04_Feb_2026_encrypted.csv")

############################################################
# 01a. CORE FILTERING
############################################################

bd_victory <- bd_df %>%
  arrange(herd_no, bd_no) %>%
  filter(all_cases != 0) %>%
  filter(bd_duration_days >= 60) %>%
  filter(feedlot_cfu_ever == 0) %>%
  filter(
    bd_initiated == "Slaughter (non permit animal)" |
      (bd_initiated == "SICTT" &
         bd_first_skin_test_type %in% c("1", "8"))
  )

############################################################
# 01b. BASIC VARIABLES
############################################################

bd_victory <- bd_victory %>%
  mutate(
    bd_start = as.Date(bd_start),
    
    bd_within_13_months =
      case_when(
        is.na(duration_between_bd) ~ 0,
        duration_between_bd < 390 ~ 1,
        TRUE ~ 0
      ),
    
    herd_size_mean =
      rowMeans(select(., mean_herd_size_jan,
                      mean_herd_size_may,
                      mean_herd_size_sep))
  )

############################################################
# 01c. HISTORY VARIABLE
############################################################

bd_victory <- bd_victory %>%
  group_by(herd_no) %>%
  arrange(bd_start, .by_group = TRUE) %>%
  mutate(
    no_bds_started_in_prev_5yr =
      sapply(
        bd_start,
        function(x)
          sum(
            bd_start < x &
              bd_start >= x %m-% years(5),
            na.rm = TRUE
          )
      )
  ) %>%
  ungroup()

############################################################
# 01d. CORE EPIDEMIOLOGY
############################################################

bd_victory <- bd_victory %>%
  mutate(
    
    total_index_reactors =
      cows_positive_index_test +
      bulls_positive_index_test +
      steers_positive_index_test +
      heifers_positive_index_test +
      calves_positive_index_test,
    
    total_tested_index = no_animals_index_test,
    
    total_subsequent_reactors =
      (cows_positive + bulls_positive + steers_positive +
         heifers_positive + calves_positive) -
      total_index_reactors,
    
    prop_positive_index =
      ifelse(total_tested_index > 0,
             total_index_reactors / total_tested_index,
             NA),
    
    spread_ratio =
      ifelse(total_index_reactors > 0,
             total_subsequent_reactors / total_index_reactors,
             NA),
    
    spread_any =
      ifelse(total_subsequent_reactors > 0, 1, 0)
  ) %>%
  mutate(
    spread_ratio = ifelse(is.infinite(spread_ratio), NA, spread_ratio)
  )

############################################################
# 02a. SELECT FINAL VARIABLES
############################################################

core_vars <- c(
  "prop_positive_index",
  "spread_ratio",
  "spread_any",
  "bd_duration_days",
  "bd_within_13_months",
  "no_bds_started_in_prev_5yr",
  "herd_size_mean"
)

bd_core <- bd_victory %>%
  select(all_of(core_vars)) %>%
  drop_na()
# Log transforming them first
bd_core_log <- bd_core %>%
  mutate(
    spread_ratio = log1p(spread_ratio),
    bd_duration_days = log1p(bd_duration_days),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )

# THEN scale
bd_core_scaled <- scale(bd_core_log)

# THEN PAM
pam3_full <- pam(bd_core_scaled, k = 3)


bd_core$cluster <- factor(pam3_full$clustering)
bd_core_log$cluster <- bd_core$cluster
##Using cluster 1 as the baseline

bd_core$cluster <- relevel(bd_core$cluster, ref = "1")


library(nnet)
model <- multinom(
  cluster ~ prop_positive_index + spread_ratio + spread_any +
    bd_duration_days + bd_within_13_months +
    no_bds_started_in_prev_5yr + herd_size_mean,
  data = bd_core_log
)

summary(model)

exp(coef(model))
z <- summary(model)$coefficients / summary(model)$standard.errors
p_values <- 2 * (1 - pnorm(abs(z)))

p_values


pred <- fitted(model)
head(pred)

##EFFECT OF HERDSIZE-What if herdsize changes and everything remains
#the same.
new_data <- data.frame(
  herd_size_mean = quantile(bd_core_log$herd_size_mean,
                            probs = seq(0.1, 0.9, length.out = 5)),

  spread_ratio = mean(bd_core_log$spread_ratio),
  spread_any = 0,
  bd_duration_days = mean(bd_core_log$bd_duration_days),
  bd_within_13_months = 0,
  no_bds_started_in_prev_5yr = mean(bd_core_log$no_bds_started_in_prev_5yr),
  prop_positive_index = mean(bd_core_log$prop_positive_index)
)

predict(model, newdata = new_data, type = "probs")

# What if prop_positive_index changes while others stay constant

# Step 1: Create baseline (average herd)
base <- data.frame(
  prop_positive_index = mean(bd_core_log$prop_positive_index),
  spread_ratio = mean(bd_core_log$spread_ratio),
  spread_any = 0,
  bd_duration_days = mean(bd_core_log$bd_duration_days),
  bd_within_13_months = 0,
  no_bds_started_in_prev_5yr = mean(bd_core_log$no_bds_started_in_prev_5yr),
  herd_size_mean = mean(bd_core_log$herd_size_mean)
)

# Step 2: Replicate baseline 5 times
test_prop <- base[rep(1, 5), ]

# Step 3: Replace with realistic values (quantiles, not arbitrary range)
test_prop$prop_positive_index <- quantile(
  bd_core_log$prop_positive_index,
  probs = seq(0.1, 0.9, length.out = 5)
)

# Step 4: Predict probabilities
probs <- predict(model, newdata = test_prop, type = "probs")

# Step 5: View results
probs

# Optional: Barplot (visualise effect)
barplot(t(probs),
        beside = TRUE,
        legend = TRUE,
        ylab = "Probability",
        ylim = c(0, 1),
        names.arg = round(test_prop$prop_positive_index, 2),
        xlab = "prop_positive_index levels")


##Spread_ratio effect

test_spread <- base[rep(1,5), ]

test_spread$spread_ratio <- quantile(
  bd_core_log$spread_ratio,
  probs = seq(0.1, 0.9, length.out = 5)
)

predict(model, newdata = test_spread, type = "probs")

##Spread_any
test_any <- base[rep(1,2), ]
test_any$spread_any <- c(0,1)

predict(model, newdata = test_any, type = "probs")

#RECURRENCE+DURATION

test_rec <- expand.grid(
  bd_within_13_months = c(0,1),
  bd_duration_days = quantile(
    bd_core_log$bd_duration_days,
    probs = seq(0.1, 0.9, length.out = 5)
  )
)

# Fill others with mean
test_rec$prop_positive_index <- mean(bd_core_log$prop_positive_index)
test_rec$spread_ratio <- mean(bd_core_log$spread_ratio)
test_rec$spread_any <- 1
test_rec$no_bds_started_in_prev_5yr <- mean(bd_core_log$no_bds_started_in_prev_5yr)
test_rec$herd_size_mean <- mean(bd_core_log$herd_size_mean)

predict(model, newdata = test_rec, type = "probs")

##BARPLOT

new_data <- data.frame(
  prop_positive_index = quantile(
    bd_core_log$prop_positive_index,
    probs = c(0.1, 0.9)
  ),
  spread_ratio = mean(bd_core_log$spread_ratio),
  spread_any = 0,
  bd_duration_days = mean(bd_core_log$bd_duration_days),
  bd_within_13_months = 0,
  no_bds_started_in_prev_5yr = mean(bd_core_log$no_bds_started_in_prev_5yr),
  herd_size_mean = mean(bd_core_log$herd_size_mean)
)

probs <- predict(model, newdata = new_data, type = "probs")

barplot(t(probs),
        beside = TRUE,
        legend = TRUE,
        names.arg = c("Low index", "High index"),
        ylab = "Probability",
        xlab = "Scenario")

##BARPLOT Spread_any
new_data <- data.frame(
  prop_positive_index = mean(bd_core_log$prop_positive_index),
  spread_ratio = mean(bd_core_log$spread_ratio),
  spread_any = c(0,1),
  bd_duration_days = mean(bd_core_log$bd_duration_days),
  bd_within_13_months = 0,
  no_bds_started_in_prev_5yr = mean(bd_core_log$no_bds_started_in_prev_5yr),
  herd_size_mean = mean(bd_core_log$herd_size_mean)
)

probs <- predict(model, newdata = new_data, type = "probs")

barplot(t(probs),
        beside = TRUE,
        legend = TRUE,
        names.arg = c("No spread", "Spread"),
        ylab = "Probability")


#RECURRENCE +DURATION
new_data <- data.frame(
  prop_positive_index = mean(bd_core_log$prop_positive_index),
  spread_ratio = mean(bd_core_log$spread_ratio),
  spread_any = 1,
  bd_duration_days = quantile(
    bd_core_log$bd_duration_days,
    probs = c(0.25, 0.75)
  ),
  bd_within_13_months = c(0,1),
  no_bds_started_in_prev_5yr = mean(bd_core_log$no_bds_started_in_prev_5yr),
  herd_size_mean = mean(bd_core_log$herd_size_mean)
)

probs <- predict(model, newdata = new_data, type = "probs")

barplot(t(probs),
        beside = TRUE,
        legend = TRUE,
        names.arg = c("No recurrence", "Recurrence"),
        ylab = "Probability")
###SO lets get started#
library(dplyr)
library(cluster)
library(tidyr)

set.seed(123)


# LOAD MULTINOMIAL DATA


multi_data <- read.csv("Data_for_MN_regression_6-5-26_encrypted.csv")

multi_data$herd_no <- as.character(multi_data$herd_no)


# BUILD BREAKDOWN CLUSTER DATA


bd_cluster_data <- bd_victory %>%
  select(herd_no, bd_no, all_of(core_vars)) %>%
  drop_na() %>%
  mutate(
    spread_ratio = log1p(spread_ratio),
    bd_duration_days = log1p(bd_duration_days),
    herd_size_mean = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr = log1p(no_bds_started_in_prev_5yr)
  )


# SCALE VARIABLES


scaled_data <- bd_cluster_data %>%
  select(-herd_no, -bd_no) %>%
  scale() %>%
  as.data.frame()


# PAM CLUSTERING


pam_model <- pam(scaled_data, k = 3)


# DATASET 1 = HERD + BREAKDOWN + CLUSTER


dataset1 <- bd_cluster_data %>%
  mutate(cluster = pam_model$clustering) %>%
  select(herd_no, bd_no, cluster)

# JOIN WITH RISK FACTOR DATA


data_bd <- dataset1 %>%
  left_join(multi_data, by = c("herd_no", "bd_no"))


# IMPORTANT # KEEP ONLY NON-CONTROL ROWS


data_bd <- data_bd %>%
  filter(status != "Control")


#  TRUE CONTROLS ONLY


control_herds <- multi_data %>%
  group_by(herd_no) %>%
  summarise(
    all_control = all(status == "Control"),
    .groups = "drop"
  ) %>%
  filter(all_control == TRUE)

controls <- multi_data %>%
  filter(herd_no %in% control_herds$herd_no) %>%
  sample_n(30000)


#ASSIGNING CLUSTER 0 TO CONTROLS


controls$cluster <- 0


# FINAL DATASET

final_data <- bind_rows(data_bd, controls)

final_data$cluster <- as.factor(final_data$cluster)

#Saving the data

write.csv(
  final_data,
  "C:/Users/Victory/Downloads/finale_dataset_with_clusters.csv",
  row.names = FALSE
)


table(final_data$cluster)

table(final_data$status, final_data$cluster)
##MULTINOMIAL REGRESSION FOR THE RISK FACTORS
rm(list = ls())

library(nnet)
library(dplyr)
main_multi_data=read.csv("finale_dataset_with_clusters.csv")
names(main_multi_data)
#converting cluster to factor

main_multi_data$cluster <- as.factor(main_multi_data$cluster)

# Set controls (cluster 0) as baseline
main_multi_data$cluster <- relevel(
  main_multi_data$cluster,
  ref = "0"
)
library(nnet)
#Fit a multimonial model

multi_model <- multinom(
  cluster ~
    Control_moves +
    NewBD_moves +
    Burnout_moves +
    Unknown_moves +
    Burnout_neighbours +
    Control_neighbours +
    NewBD_neighbours,
  
  data = main_multi_data
)

#Model summary

summary(multi_model)

#ODD RATIOS
odds_ratios <- exp(coef(multi_model))

odds_ratios

#P-VALUES

z <- summary(multi_model)$coefficients /
  summary(multi_model)$standard.errors

p_values <- 2 * (1 - pnorm(abs(z)))

p_values

#Combine results in table

results_table <- data.frame(
  Variable = colnames(summary(multi_model)$coefficients),
  
  Cluster1_OR = odds_ratios[1, ],
  Cluster1_p = p_values[1, ],
  
  Cluster2_OR = odds_ratios[2, ],
  Cluster2_p = p_values[2, ],
  
  Cluster3_OR = odds_ratios[3, ],
  Cluster3_p = p_values[3, ]
)

results_table

# Save results
write.csv(
  results_table,
  "C:/Users/Victory/Downloads/multinomial_results.csv",
  row.names = FALSE
)


# PREDICTED PROBABILITIES


predicted_probs <- fitted(multi_model)

head(predicted_probs)


#AVERAGE PREDICTED PROBABILITIES


colMeans(predicted_probs)


#CLUSTER COUNTS


table(main_multi_data$cluster)


# STATUS vs CLUSTER CHECK


table(
  main_multi_data$status,
  main_multi_data$cluster
)

############################################################
# EFFECT OF NewBD_moves
############################################################

# Baseline herd
base <- data.frame(
  Control_moves = mean(main_multi_data$Control_moves, na.rm = TRUE),
  NewBD_moves = mean(main_multi_data$NewBD_moves, na.rm = TRUE),
  Burnout_moves = mean(main_multi_data$Burnout_moves, na.rm = TRUE),
  Unknown_moves = mean(main_multi_data$Unknown_moves, na.rm = TRUE),
  Burnout_neighbours = mean(main_multi_data$Burnout_neighbours, na.rm = TRUE),
  Control_neighbours = mean(main_multi_data$Control_neighbours, na.rm = TRUE),
  NewBD_neighbours = mean(main_multi_data$NewBD_neighbours, na.rm = TRUE)
)

# Repeat 5 rows
test_moves <- base[rep(1, 5), , drop = FALSE]

# Change ONLY NewBD_moves
test_moves$NewBD_moves <- as.numeric(
  quantile(
    main_multi_data$NewBD_moves,
    probs = seq(0.1, 0.9, length.out = 5),
    na.rm = TRUE
  )
)

# Predict probabilities
probs <- predict(
  multi_model,
  newdata = test_moves,
  type = "probs"
)

# View
probs

# CREATE CLEAN DATASET USED BY BOTH MODELS

model_data <- main_multi_data %>%
  select(
    cluster,
    Control_moves,
    NewBD_moves,
    Burnout_moves,
    Unknown_moves,
    Burnout_neighbours,
    Control_neighbours,
    NewBD_neighbours
  ) %>%
  na.omit()

# FULL MODEL
multi_model <- multinom(
  cluster ~
    Control_moves +
    NewBD_moves +
    Burnout_moves +
    Unknown_moves +
    Burnout_neighbours +
    Control_neighbours +
    NewBD_neighbours,
  
  data = model_data
)

# NULL MODEL
null_model <- multinom(
  cluster ~ 1,
  data = model_data
)

# LIKELIHOOD RATIO TEST
library(lmtest)

lrtest(null_model, multi_model)

# McFadden pseudo-R2
pseudo_r2 <- 1 - (
  as.numeric(logLik(multi_model)) /
    as.numeric(logLik(null_model))
)

pseudo_r2


# FUNCTION TO TEST EFFECT OF ONE PREDICTOR


test_predictor_effect <- function(var_name) {
  
  # Baseline herd
  base <- data.frame(
    Control_moves = mean(main_multi_data$Control_moves, na.rm = TRUE),
    NewBD_moves = mean(main_multi_data$NewBD_moves, na.rm = TRUE),
    Burnout_moves = mean(main_multi_data$Burnout_moves, na.rm = TRUE),
    Unknown_moves = mean(main_multi_data$Unknown_moves, na.rm = TRUE),
    Burnout_neighbours = mean(main_multi_data$Burnout_neighbours, na.rm = TRUE),
    Control_neighbours = mean(main_multi_data$Control_neighbours, na.rm = TRUE),
    NewBD_neighbours = mean(main_multi_data$NewBD_neighbours, na.rm = TRUE)
  )
  
  # Repeat rows
  test_data <- base[rep(1, 5), , drop = FALSE]
  
  # Increase ONLY one variable
  test_data[[var_name]] <- as.numeric(
    quantile(
      main_multi_data[[var_name]],
      probs = seq(0.1, 0.9, length.out = 5),
      na.rm = TRUE
    )
  )
  
  # Predicted probabilities
  probs <- predict(
    multi_model,
    newdata = test_data,
    type = "probs"
  )
  
  # Final table
  results <- data.frame(
    Predictor = var_name,
    Changed_Value = test_data[[var_name]],
    Prob_Cluster0 = probs[,1],
    Prob_Cluster1 = probs[,2],
    Prob_Cluster2 = probs[,3],
    Prob_Cluster3 = probs[,4]
  )
  
  return(results)
}



all_results <- rbind(
  test_predictor_effect("Control_moves"),
  test_predictor_effect("NewBD_moves"),
  test_predictor_effect("Burnout_moves"),
  test_predictor_effect("Unknown_moves"),
  test_predictor_effect("Burnout_neighbours"),
  test_predictor_effect("Control_neighbours"),
  test_predictor_effect("NewBD_neighbours")
)



print(all_results)



write.csv(
  all_results,
  "C:/Users/Victory/Downloads/predictor_effects.csv",
  row.names = FALSE
)