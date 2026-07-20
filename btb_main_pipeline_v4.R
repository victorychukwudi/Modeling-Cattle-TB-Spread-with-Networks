############################################################
# bTB BREAKDOWN CLUSTERING ANALYSIS


options(bitmapType = "cairo")
setwd("C:/Users/Victory/Downloads")


# LIBRARIES


library(tidyverse)
library(lubridate)
library(matrixStats)
library(cluster)
library(factoextra)
library(corrplot)
library(gridExtra)
library(NbClust)
library(clustertend)
library(ggfortify)
library(mclust)
library(fpc)
library(lmtest)
library(ClusterR)
library(nnet)
library(survival)
library(survminer)
library(pROC)
library(ppclust)

############################################################
# 01a. LOAD DATA
############################################################

bd_df <- read.csv("bd_df_04_Feb_2026_encrypted.csv")

############################################################
# 01b. CORE FILTERING
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
# 01c. BASIC VARIABLES
############################################################

bd_victory <- bd_victory %>%
  mutate(
    bd_start = as.Date(bd_start),

    bd_within_13_months = case_when(
      is.na(duration_between_bd) ~ 0,
      duration_between_bd < 390  ~ 1,
      TRUE                       ~ 0
    ),

    herd_size_mean = rowMeans(
      select(., mean_herd_size_jan,
             mean_herd_size_may,
             mean_herd_size_sep)
    )
  )

############################################################
# 01d. HISTORY VARIABLE
############################################################

bd_victory <- bd_victory %>%
  group_by(herd_no) %>%
  arrange(bd_start, .by_group = TRUE) %>%
  mutate(
    no_bds_started_in_prev_5yr = sapply(
      bd_start,
      function(x)
        sum(bd_start < x & bd_start >= x %m-% years(5),
            na.rm = TRUE)
    )
  ) %>%
  ungroup()

############################################################
# 01e. CORE EPIDEMIOLOGY
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

    prop_positive_index = ifelse(
      total_tested_index > 0,
      total_index_reactors / total_tested_index,
      NA
    ),

    spread_ratio = ifelse(
      total_index_reactors > 0,
      total_subsequent_reactors / total_index_reactors,
      NA
    ),

    spread_any = ifelse(total_subsequent_reactors > 0, 1, 0)
  ) %>%
  mutate(
    spread_ratio = ifelse(is.infinite(spread_ratio), NA, spread_ratio)
  )

############################################################
# 02a. DEFINE CLUSTERING VARIABLES
############################################################

cluster_vars <- c(
  "prop_positive_index",
  "spread_ratio",
  "bd_duration_days",
  "no_bds_started_in_prev_5yr",
  "herd_size_mean"
)

bd_core <- bd_victory %>%
  select(
    herd_no, bd_no, bd_start, bd_start_yr,
    duration_between_bd, bd_within_13_months, spread_any,
    all_of(cluster_vars)
  ) %>%
  drop_na(all_of(cluster_vars))

############################################################
# 02b. LOG TRANSFORMATION
############################################################

bd_core <- bd_core %>%
  mutate(
    spread_ratio_log               = log1p(spread_ratio),
    bd_duration_days_log           = log1p(bd_duration_days),
    herd_size_mean_log             = log1p(herd_size_mean),
    no_bds_started_in_prev_5yr_log = log1p(no_bds_started_in_prev_5yr)
  )

cluster_vars_log <- c(
  "prop_positive_index",
  "spread_ratio_log",
  "bd_duration_days_log",
  "no_bds_started_in_prev_5yr_log",
  "herd_size_mean_log"
)

############################################################
# 02c. SCALE → CANONICAL CLUSTERING MATRIX
############################################################

bd_scaled_full <- bd_core %>%
  select(all_of(cluster_vars_log)) %>%
  scale() %>%
  as.data.frame()

############################################################
# 02d. DESCRIPTIVE TABLE
############################################################

desc_table <- bd_core %>%
  select(all_of(cluster_vars), bd_within_13_months, spread_any) %>%
  summarise(across(
    everything(),
    list(
      Mean = ~mean(., na.rm = TRUE),
      SD   = ~sd(.,   na.rm = TRUE),
      Min  = ~min(.,  na.rm = TRUE),
      Max  = ~max(.,  na.rm = TRUE)
    )
  )) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to  = c("Variable", ".value"),
    names_pattern = "(.*)_(Mean|SD|Min|Max)"
  )

print(desc_table)
write.csv(desc_table, "00_descriptive_table.csv", row.names = FALSE)

############################################################
# ============= EXPLORATION =============
############################################################

png("01a_hist_raw.png", 1600, 1200)
bd_core %>%
  select(all_of(cluster_vars)) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill = "steelblue") +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  ggtitle("Raw distributions (clustering variables)")
dev.off()

png("01b_hist_log.png", 1600, 1200)
bd_core %>%
  select(all_of(cluster_vars_log)) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill = "darkgreen") +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  ggtitle("Log-transformed distributions (clustering variables)")
dev.off()

png("01c_density.png", 1600, 1200)
bd_core %>%
  select(all_of(cluster_vars_log)) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value)) +
  geom_density(fill = "orange", alpha = 0.5) +
  facet_wrap(~name, scales = "free") +
  theme_bw()
dev.off()

png("01d_scatter_matrix.png", 1600, 1400)
pairs(bd_scaled_full, pch = 20, cex = 0.4)
dev.off()

png("01f_correlation.png", 1200, 1000)
corrplot(
  cor(bd_scaled_full, use = "pairwise.complete.obs"),
  method = "color", type = "upper", order = "hclust",
  addCoef.col = "black", number.cex = 0.5, tl.cex = 0.6
)
dev.off()

png("01g_boxplots.png", 1600, 1200)
bd_scaled_full %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot(fill = "grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Scaled clustering variables")
dev.off()

# Severity vs spread scatter
png("01e_scatter_severity_vs_spread.png", 1200, 900)
plot(bd_core$prop_positive_index, log1p(bd_core$spread_ratio),
     pch=19, col=rgb(0,0,1,0.05),
     xlab="Proportion positive at index test (severity)",
     ylab="Log spread ratio",
     main="Index severity vs within-herd spread")
dev.off()



############################################################
# ============= PCA =============
############################################################

bd_pca <- prcomp(bd_scaled_full, center = FALSE, scale. = FALSE)
print(summary(bd_pca))

png("04a_variance_explained.png", 1200, 900)
fviz_eig(bd_pca, addlabels = TRUE, ylim = c(0, 70))
dev.off()

set.seed(123)
pca_plot_idx <- sample(nrow(bd_scaled_full), min(5000, nrow(bd_scaled_full)))

png("04b_pca_individuals.png", 1200, 900)
fviz_pca_ind(
  bd_pca, geom = "point", pointsize = 0.8, alpha.ind = 0.3,
  select.ind = list(idx = pca_plot_idx), ggtheme = theme_bw()
)
dev.off()

png("04c_var_contrib_PC1.png", 1200, 900)
fviz_contrib(bd_pca, choice = "var", axes = 1)
dev.off()

png("04d_var_contrib_PC2.png", 1200, 900)
fviz_contrib(bd_pca, choice = "var", axes = 2)
dev.off()

png("04e_pca_biplot.png", 1200, 900)
autoplot(
  bd_pca, loadings = TRUE, loadings.colour = "darkblue",
  loadings.label = TRUE, loadings.label.size = 3
)
dev.off()



############################################################
# K-MEANS ON PCA SPACE (first 3 PCs)
# Exploratory — shows cluster separation in reduced space
############################################################



bd_pca_scores  <- as.data.frame(bd_pca$x)
bd_pca_reduced <- bd_pca_scores[, 1:3]

# Subsample for plotting
set.seed(123)
pca_plot_idx   <- sample(nrow(bd_pca_reduced), min(5000, nrow(bd_pca_reduced)))
bd_pca_sub     <- bd_pca_reduced[pca_plot_idx, ]

# Run k-means on the subsampled PCA scores
set.seed(123)
pca_km2 <- kmeans(bd_pca_sub, centers = 2, nstart = 25)
pca_km3 <- kmeans(bd_pca_sub, centers = 3, nstart = 25)
pca_km4 <- kmeans(bd_pca_sub, centers = 4, nstart = 25)
pca_km5 <- kmeans(bd_pca_sub, centers = 5, nstart = 25)

# Use fviz_cluster directly on the subsampled PCA scores
pp2 <- fviz_cluster(list(data = bd_pca_sub, cluster = pca_km2$cluster),
                    geom = "point", ellipse.type = "norm") + ggtitle("K-means on PCA - K=2")
pp3 <- fviz_cluster(list(data = bd_pca_sub, cluster = pca_km3$cluster),
                    geom = "point", ellipse.type = "norm") + ggtitle("K-means on PCA - K=3")
pp4 <- fviz_cluster(list(data = bd_pca_sub, cluster = pca_km4$cluster),
                    geom = "point", ellipse.type = "norm") + ggtitle("K-means on PCA - K=4")
pp5 <- fviz_cluster(list(data = bd_pca_sub, cluster = pca_km5$cluster),
                    geom = "point", ellipse.type = "norm") + ggtitle("K-means on PCA - K=5")

png("pca_kmeans_clusters.png", 1600, 1200)
gridExtra::grid.arrange(pp2, pp3, pp4, pp5, ncol = 2)
dev.off()

















############################################################
# ============= CLUSTER TENDENCY =============
############################################################

set.seed(123)
tend_idx  <- sample(nrow(bd_scaled_full), min(2000, nrow(bd_scaled_full)))
tend_data <- as.matrix(bd_scaled_full[tend_idx, ])

hopkins_stat <- hopkins(tend_data, n = floor(0.1 * nrow(tend_data)))
cat("Hopkins statistic:", hopkins_stat$H, "\n")

png("02a_vat_plot.png", 1200, 1000)
fviz_dist(dist(tend_data))
dev.off()


############################################################
# ============= OPTIMAL K SELECTION =============
############################################################

set.seed(123)
nb_idx    <- sample(nrow(bd_scaled_full), min(5000, nrow(bd_scaled_full)))
nb_sample <- bd_scaled_full[nb_idx, ]

png("02b_wss.png", 1200, 900)
fviz_nbclust(nb_sample, kmeans, method = "wss", k.max = 10) +
  ggtitle("Elbow Method (WSS) – subsample n = 5,000")
dev.off()

png("02c_silhouette_method.png", 1200, 900)
fviz_nbclust(nb_sample, kmeans, method = "silhouette", k.max = 10) +
  ggtitle("Silhouette Method – subsample n = 5,000")
dev.off()

png("02d_gap_stat.png", 1200, 900)
fviz_nbclust(nb_sample, kmeans, method = "gap_stat", k.max = 10, nboot = 50) +
  ggtitle("Gap Statistic – subsample n = 5,000")
dev.off()

set.seed(123)
nb_all <- NbClust(
  data = nb_sample, distance = "euclidean",
  min.nc = 2, max.nc = 10, method = "ward.D2"
)

png("02e_nbclust_barplot.png", 1200, 900)
barplot(table(nb_all$Best.nc[1, ]), col = "grey",
        main = "NbClust: Majority vote for optimal K")
dev.off()



# CLUSTERR ELBOW — variance explained curve

library(ClusterR)

png("02b2_clusterR_elbow.png", 1200, 900)
set.seed(123)
opt <- Optimal_Clusters_KMeans(
  nb_sample,
  max_clusters  = 10,
  plot_clusters = TRUE,
  verbose       = FALSE
)
dev.off()




############################################################
# ============= HIERARCHICAL CLUSTERING (exploratory) =============
############################################################

set.seed(123)
hc_idx    <- sample(nrow(bd_scaled_full), min(3000, nrow(bd_scaled_full)))
hc_sample <- bd_scaled_full[hc_idx, ]
dist_hc   <- dist(hc_sample, method = "euclidean")
hc        <- hclust(dist_hc, method = "ward.D2")

png("02f_dendrogram.png", 1600, 1000)
fviz_dend(hc, k = 3, rect = TRUE, show_labels = FALSE,
          main = "Hierarchical Clustering – Ward's Method (subsample n = 3,000)")
dev.off()




############################################################
# ============= CANONICAL PAM CLUSTERING (K = 3) =============
############################################################

set.seed(123)
pam3_canonical <- pam(bd_scaled_full, k = 3)

bd_core$cluster <- factor(pam3_canonical$clustering)

cat("PAM K=3 cluster sizes:\n")
print(table(bd_core$cluster))

write.csv(as.data.frame(pam3_canonical$medoids), "pam_k3_medoids.csv", row.names = FALSE)

pam3_means_raw <- bd_core %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    across(c(all_of(cluster_vars), bd_within_13_months, spread_any), mean, na.rm = TRUE),
    .groups = "drop"
  )

print(pam3_means_raw)
write.csv(pam3_means_raw, "pam_k3_cluster_means.csv", row.names = FALSE)


# ========= PAM K COMPARISON (K = 2 to 5) ======


set.seed(123)
pam2 <- pam(bd_scaled_full, k = 2)
pam4 <- pam(bd_scaled_full, k = 4)
pam5 <- pam(bd_scaled_full, k = 5)

set.seed(123)
vis_idx  <- sample(nrow(bd_scaled_full), min(5000, nrow(bd_scaled_full)))
vis_data <- bd_scaled_full[vis_idx, ]

p_pam2 <- fviz_cluster(list(data=vis_data, cluster=pam2$clustering[vis_idx]),
                       geom="point", ellipse.type="norm") + ggtitle("PAM K = 2")
p_pam3 <- fviz_cluster(list(data=vis_data, cluster=pam3_canonical$clustering[vis_idx]),
                       geom="point", ellipse.type="norm") + ggtitle("PAM K = 3")
p_pam4 <- fviz_cluster(list(data=vis_data, cluster=pam4$clustering[vis_idx]),
                       geom="point", ellipse.type="norm") + ggtitle("PAM K = 4")
p_pam5 <- fviz_cluster(list(data=vis_data, cluster=pam5$clustering[vis_idx]),
                       geom="point", ellipse.type="norm") + ggtitle("PAM K = 5")

png("03a_pam_cluster_plots.png", 1600, 1200)
gridExtra::grid.arrange(p_pam2, p_pam3, p_pam4, p_pam5, ncol = 2)
dev.off()

png("03b_pam_k3_plot.png", 1200, 900)
print(p_pam3)
dev.off()

p_pca2 <- fviz_pca_ind(bd_pca, geom="point", habillage=factor(pam2$clustering)[vis_idx],
                       addEllipses=TRUE, ellipse.level=0.95, pointsize=0.8, alpha.ind=0.4,
                       select.ind=list(idx=vis_idx)) + ggtitle("PAM K = 2 on PCA")
p_pca3 <- fviz_pca_ind(bd_pca, geom="point", habillage=factor(pam3_canonical$clustering)[vis_idx],
                       addEllipses=TRUE, ellipse.level=0.95, pointsize=0.8, alpha.ind=0.4,
                       select.ind=list(idx=vis_idx)) + ggtitle("PAM K = 3 on PCA")
p_pca4 <- fviz_pca_ind(bd_pca, geom="point", habillage=factor(pam4$clustering)[vis_idx],
                       addEllipses=TRUE, ellipse.level=0.95, pointsize=0.8, alpha.ind=0.4,
                       select.ind=list(idx=vis_idx)) + ggtitle("PAM K = 4 on PCA")
p_pca5 <- fviz_pca_ind(bd_pca, geom="point", habillage=factor(pam5$clustering)[vis_idx],
                       addEllipses=TRUE, ellipse.level=0.95, pointsize=0.8, alpha.ind=0.4,
                       select.ind=list(idx=vis_idx)) + ggtitle("PAM K = 5 on PCA")

png("03c_pca_pam_clusters.png", 1600, 1200)
gridExtra::grid.arrange(p_pca2, p_pca3, p_pca4, p_pca5, ncol = 2)
dev.off()

#PAM SILHOUETTE
set.seed(123)
sil_idx  <- sample(nrow(bd_scaled_full), min(5000, nrow(bd_scaled_full)))
sil_data <- bd_scaled_full[sil_idx, ]
dist_sil <- dist(sil_data)

pam_sil2 <- silhouette(pam2$clustering[sil_idx],           dist_sil)
pam_sil3 <- silhouette(pam3_canonical$clustering[sil_idx], dist_sil)
pam_sil4 <- silhouette(pam4$clustering[sil_idx],           dist_sil)
pam_sil5 <- silhouette(pam5$clustering[sil_idx],           dist_sil)

sp2 <- fviz_silhouette(pam_sil2) + ggtitle("PAM K = 2")
sp3 <- fviz_silhouette(pam_sil3) + ggtitle("PAM K = 3")
sp4 <- fviz_silhouette(pam_sil4) + ggtitle("PAM K = 4")
sp5 <- fviz_silhouette(pam_sil5) + ggtitle("PAM K = 5")

png("03d_pam_silhouette_all.png", 1600, 1200)
grid.arrange(sp2, sp3, sp4, sp5, ncol = 2)
dev.off()

png("03e_pam_silhouette_k3.png", 900, 700)
print(sp3)
dev.off()

avg_sil <- function(s) mean(s[, 3])

sil_summary <- data.frame(
  K = c(2, 3, 4, 5),
  Avg_Silhouette = c(avg_sil(pam_sil2), avg_sil(pam_sil3),
                     avg_sil(pam_sil4), avg_sil(pam_sil5))
)

print(sil_summary)
write.csv(sil_summary, "pam_silhouette_summary.csv", row.names = FALSE)


############################################################
# K-MEANS MODELS (K = 2 to 5) — comparator method
############################################################

set.seed(123)
bd_clusters2 <- kmeans(bd_scaled_full, centers = 2, nstart = 25)
bd_clusters3 <- kmeans(bd_scaled_full, centers = 3, nstart = 25)
bd_clusters4 <- kmeans(bd_scaled_full, centers = 4, nstart = 25)
bd_clusters5 <- kmeans(bd_scaled_full, centers = 5, nstart = 25)

############################################################
# K-MEANS CLUSTER MEANS (saved on subsample for speed)
############################################################

km_df <- as.data.frame(bd_scaled_full)

k2_means <- km_df %>%
  mutate(cluster = bd_clusters2$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

k3_means <- km_df %>%
  mutate(cluster = bd_clusters3$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

k4_means <- km_df %>%
  mutate(cluster = bd_clusters4$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

k5_means <- km_df %>%
  mutate(cluster = bd_clusters5$cluster) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

write.csv(k2_means, "kmeans_k2_means.csv", row.names = FALSE)
write.csv(k3_means, "kmeans_k3_means.csv", row.names = FALSE)
write.csv(k4_means, "kmeans_k4_means.csv", row.names = FALSE)
write.csv(k5_means, "kmeans_k5_means.csv", row.names = FALSE)

############################################################
# K-MEANS CLUSTER VISUALS (subsample for plot speed)
############################################################

km_p2 <- fviz_cluster(
  list(data = vis_data, cluster = bd_clusters2$cluster[vis_idx]),
  geom = "point", ellipse.type = "norm"
) + ggtitle("K-means K = 2")

km_p3 <- fviz_cluster(
  list(data = vis_data, cluster = bd_clusters3$cluster[vis_idx]),
  geom = "point", ellipse.type = "norm"
) + ggtitle("K-means K = 3")

km_p4 <- fviz_cluster(
  list(data = vis_data, cluster = bd_clusters4$cluster[vis_idx]),
  geom = "point", ellipse.type = "norm"
) + ggtitle("K-means K = 4")

km_p5 <- fviz_cluster(
  list(data = vis_data, cluster = bd_clusters5$cluster[vis_idx]),
  geom = "point", ellipse.type = "norm"
) + ggtitle("K-means K = 5")

png("kmeans_cluster_plots.png", 1600, 1200)
gridExtra::grid.arrange(km_p2, km_p3, km_p4, km_p5, ncol = 2)
dev.off()


# K-MEANS SILHOUETTE 

km_sil2 <- silhouette(bd_clusters2$cluster[sil_idx], dist_sil)
km_sil3 <- silhouette(bd_clusters3$cluster[sil_idx], dist_sil)
km_sil4 <- silhouette(bd_clusters4$cluster[sil_idx], dist_sil)
km_sil5 <- silhouette(bd_clusters5$cluster[sil_idx], dist_sil)

ks2 <- fviz_silhouette(km_sil2) + ggtitle("K-means K = 2")
ks3 <- fviz_silhouette(km_sil3) + ggtitle("K-means K = 3")
ks4 <- fviz_silhouette(km_sil4) + ggtitle("K-means K = 4")
ks5 <- fviz_silhouette(km_sil5) + ggtitle("K-means K = 5")

png("kmeans_silhouette_plots.png", 1600, 1200)
grid.arrange(ks2, ks3, ks4, ks5, ncol = 2)
dev.off()


# K-MEANS vs PAM SILHOUETTE COMPARISON


avg_sil <- function(s) mean(s[, 3])

sil_comparison <- data.frame(
  Method = c(rep("K-means", 4), rep("PAM", 4)),
  K      = rep(c(2, 3, 4, 5), 2),
  Avg_Silhouette = c(
    avg_sil(km_sil2), avg_sil(km_sil3),
    avg_sil(km_sil4), avg_sil(km_sil5),
    avg_sil(pam_sil2), avg_sil(pam_sil3),
    avg_sil(pam_sil4), avg_sil(pam_sil5)
  )
)

print(sil_comparison)
write.csv(sil_comparison, "kmeans_vs_pam_silhouette.csv", row.names = FALSE)

png("kmeans_vs_pam_silhouette_plot.png", 1200, 900)
ggplot(sil_comparison, aes(x = factor(K), y = Avg_Silhouette, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(
    title = "K-means vs PAM Silhouette Comparison",
    x     = "Number of Clusters (K)",
    y     = "Average Silhouette Width"
  )
dev.off()


# CLUSTER STABILITY (clusterboot) -K = 3

set.seed(123)
cb_pam3 <- clusterboot(
  bd_scaled_full, B = 200, clustermethod = pamkCBI,
  krange = 3, seed = 123, multipleboot = FALSE, showplots = FALSE
)

cat("Cluster stability (Jaccard, B = 200):\n")
print(cb_pam3$bootmean)


# ============= CLUSTER PROPORTION BY YEAR =============

cluster_year_prop <- bd_core %>%
  group_by(bd_start_yr, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bd_start_yr) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

write.csv(cluster_year_prop, "cluster_proportion_by_year.csv", row.names = FALSE)

png("04f_cluster_proportion_by_year.png", 1200, 900)
ggplot(cluster_year_prop, aes(x = bd_start_yr, y = proportion, fill = factor(cluster))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Proportion of PAM (K = 3) Clusters by Year",
       x = "Year", y = "Proportion", fill = "Cluster")
dev.off()


############################################################
# ============= SURVIVAL ANALYSIS 1: BD DURATION =============
############################################################

surv_dur <- bd_core %>% select(bd_duration_days, cluster) %>% drop_na()

surv_obj_dur <- Surv(time = surv_dur$bd_duration_days, event = rep(1, nrow(surv_dur)))
km_fit_dur   <- survfit(surv_obj_dur ~ cluster, data = surv_dur)

png("05a_km_duration.png", 1400, 1000)
print(ggsurvplot(km_fit_dur, data = surv_dur, pval = TRUE, conf.int = TRUE,
                 risk.table = TRUE, legend.title = "Cluster",
                 xlab = "Breakdown duration (days)",
                 ylab = "Probability of remaining in breakdown"))
dev.off()


# ============= SURVIVAL ANALYSIS 2: TIME TO NEXT BD =============


surv_next <- bd_core %>%
  select(duration_between_bd, cluster) %>%
  mutate(
    time  = ifelse(is.na(duration_between_bd), 1000, duration_between_bd),
    event = ifelse(is.na(duration_between_bd), 0, 1)
  )

surv_obj_next <- Surv(time = surv_next$time, event = surv_next$event)
km_fit_next   <- survfit(surv_obj_next ~ cluster, data = surv_next)

png("05b_km_next_breakdown.png", 1400, 1000)
print(ggsurvplot(km_fit_next, data = surv_next, pval = TRUE, conf.int = TRUE,
                 risk.table = TRUE, legend.title = "Cluster",
                 xlab = "Time to next breakdown (days)",
                 ylab = "Probability of remaining breakdown-free",
                 xlim = c(0, 1000), break.time.by = 100))
dev.off()



############################################################
# ============= INTERNAL MULTINOMIAL REGRESSION =============
############################################################

mn_internal_data <- bd_core %>%
  select(cluster, all_of(cluster_vars_log), bd_within_13_months, spread_any) %>%
  drop_na()

mn_internal_data$cluster <- relevel(mn_internal_data$cluster, ref = "1")

set.seed(123)
model_internal <- multinom(
  cluster ~ prop_positive_index + spread_ratio_log + spread_any +
    bd_duration_days_log + bd_within_13_months +
    no_bds_started_in_prev_5yr_log + herd_size_mean_log,
  data = mn_internal_data, trace = FALSE
)

cat("\n--- INTERNAL MULTINOMIAL (descriptive only) ---\n")
print(summary(model_internal))

or_internal <- exp(coef(model_internal))
z_internal  <- summary(model_internal)$coefficients / summary(model_internal)$standard.errors
p_internal  <- 2 * (1 - pnorm(abs(z_internal)))

cat("Odds ratios:\n"); print(or_internal)
cat("P-values:\n");    print(p_internal)

base_internal <- data.frame(
  prop_positive_index            = mean(mn_internal_data$prop_positive_index,            na.rm=TRUE),
  spread_ratio_log               = mean(mn_internal_data$spread_ratio_log,               na.rm=TRUE),
  spread_any                     = 0,
  bd_duration_days_log           = mean(mn_internal_data$bd_duration_days_log,           na.rm=TRUE),
  bd_within_13_months            = 0,
  no_bds_started_in_prev_5yr_log = mean(mn_internal_data$no_bds_started_in_prev_5yr_log, na.rm=TRUE),
  herd_size_mean_log             = mean(mn_internal_data$herd_size_mean_log,             na.rm=TRUE)
)

vary_one <- function(var_name, data_log, base_df, model) {
  vals  <- quantile(data_log[[var_name]], probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
  td    <- base_df[rep(1, 3), ]
  td[[var_name]] <- vals
  probs <- predict(model, newdata = td, type = "probs")
  data.frame(Variable = var_name, Level = c("10th","50th","90th"), as.data.frame(probs))
}

int_vars <- c("prop_positive_index","spread_ratio_log","bd_duration_days_log",
              "no_bds_started_in_prev_5yr_log","herd_size_mean_log")

int_effects <- bind_rows(lapply(int_vars, vary_one,
                                data_log = mn_internal_data, base_df = base_internal, model = model_internal))

write.csv(int_effects, "internal_mn_marginal_effects.csv", row.names = FALSE)




############################################################
# INTERNAL MULTINOMIAL — MARGINAL EFFECT BARPLOTS
############################################################

# Severity (prop_positive_index)
nd_sev <- base_internal[rep(1, 2), ]
nd_sev$prop_positive_index <- quantile(
  mn_internal_data$prop_positive_index, probs = c(0.10, 0.90), na.rm = TRUE
)
probs_sev <- predict(model_internal, newdata = nd_sev, type = "probs")

png("internal_mn_severity_barplot.png", 900, 700)
barplot(t(probs_sev), beside = TRUE, legend = TRUE,
        names.arg = c("Low severity", "High severity"),
        ylab = "Predicted probability", col = 2:4,
        main = "Effect of index severity on cluster membership")
dev.off()

# Spread any (0 vs 1)
nd_spread <- base_internal[rep(1, 2), ]
nd_spread$spread_any <- c(0, 1)
probs_spread <- predict(model_internal, newdata = nd_spread, type = "probs")

png("internal_mn_spreadany_barplot.png", 900, 700)
barplot(t(probs_spread), beside = TRUE, legend = TRUE,
        names.arg = c("No spread", "Spread occurred"),
        ylab = "Predicted probability", col = 2:4,
        main = "Effect of within-herd spread on cluster membership")
dev.off()

# Recurrence (bd_within_13_months) x Duration
nd_rec <- expand.grid(
  bd_within_13_months = c(0, 1),
  bd_duration_days_log = quantile(
    mn_internal_data$bd_duration_days_log,
    probs = c(0.25, 0.75), na.rm = TRUE
  )
)
nd_rec$prop_positive_index            <- mean(mn_internal_data$prop_positive_index,            na.rm = TRUE)
nd_rec$spread_ratio_log               <- mean(mn_internal_data$spread_ratio_log,               na.rm = TRUE)
nd_rec$spread_any                     <- 1
nd_rec$no_bds_started_in_prev_5yr_log <- mean(mn_internal_data$no_bds_started_in_prev_5yr_log, na.rm = TRUE)
nd_rec$herd_size_mean_log             <- mean(mn_internal_data$herd_size_mean_log,             na.rm = TRUE)

probs_rec <- predict(model_internal, newdata = nd_rec, type = "probs")

png("internal_mn_recurrence_barplot.png", 900, 700)
barplot(t(probs_rec), beside = TRUE, legend = TRUE,
        names.arg = c("No recur/Short", "No recur/Long",
                      "Recur/Short",    "Recur/Long"),
        ylab = "Predicted probability", col = 2:4,
        main = "Effect of recurrence and duration on cluster membership",
        las = 2)
dev.off()

############################################################
# ============= MAIN ANALYSIS: RISK FACTOR MULTINOMIAL =============
############################################################

multi_data <- read.csv("Data_for_MN_regression_6-5-26_encrypted.csv")
multi_data$herd_no <- as.character(multi_data$herd_no)
bd_core$herd_no    <- as.character(bd_core$herd_no)

bd_with_clusters <- bd_core %>% select(herd_no, bd_no, cluster, bd_start_yr)

data_bd <- bd_with_clusters %>%
  left_join(multi_data, by = c("herd_no", "bd_no")) %>%
  filter(!is.na(cluster)) %>%
  filter(status != "Control")

control_herds <- multi_data %>%
  group_by(herd_no) %>%
  summarise(all_control = all(status == "Control"), .groups = "drop") %>%
  filter(all_control == TRUE)

set.seed(123)
controls <- multi_data %>%
  filter(herd_no %in% control_herds$herd_no) %>%
  sample_n(min(30000, n()))

############################################################
# FINAL DATASET
############################################################

data_bd$cluster  <- as.character(data_bd$cluster)
controls$cluster <- "0"

final_data         <- bind_rows(data_bd, controls)
final_data$cluster <- as.factor(final_data$cluster)
final_data$cluster <- relevel(final_data$cluster, ref = "0")

cat("\nFinal dataset dimensions:", dim(final_data), "\n")
cat("Cluster distribution:\n");  print(table(final_data$cluster))
cat("Status x Cluster:\n");      print(table(final_data$status, final_data$cluster))

write.csv(final_data, "final_dataset_with_clusters.csv", row.names = FALSE)

############################################################
# MAIN MODEL DATASET
# No bd_start_yr — controls have no breakdown year so including
# it in na.omit() would silently drop all 30,000 controls.
# Neighbour NAs filled with 0 for controls (no neighbours
# in breakdown since they never had one).
############################################################

model_data <- final_data %>%
  select(
    cluster,
    Control_moves, NewBD_moves, Burnout_moves, Unknown_moves,
    Burnout_neighbours, Control_neighbours, NewBD_neighbours
  ) %>%
  mutate(across(
    c(Burnout_neighbours, Control_neighbours, NewBD_neighbours),
    ~replace_na(., 0)
  )) %>%
  na.omit()

cat("\nmodel_data cluster table:\n")
print(table(model_data$cluster))

############################################################
# MAIN MULTINOMIAL MODEL
############################################################

set.seed(123)
multi_model <- multinom(
  cluster ~ Control_moves + NewBD_moves + Burnout_moves + Unknown_moves +
    Burnout_neighbours + Control_neighbours + NewBD_neighbours,
  data = model_data, trace = FALSE
)

cat("\n--- RISK FACTOR MULTINOMIAL ---\n")
print(summary(multi_model))

odds_ratios <- exp(coef(multi_model))
z_main      <- summary(multi_model)$coefficients / summary(multi_model)$standard.errors
p_values    <- 2 * (1 - pnorm(abs(z_main)))

cat("Odds ratios:\n"); print(odds_ratios)
cat("P-values:\n");    print(p_values)

results_table <- data.frame(
  Variable    = colnames(summary(multi_model)$coefficients),
  Cluster1_OR = odds_ratios["1", ], Cluster1_p = p_values["1", ],
  Cluster2_OR = odds_ratios["2", ], Cluster2_p = p_values["2", ],
  Cluster3_OR = odds_ratios["3", ], Cluster3_p = p_values["3", ]
)

print(results_table)
write.csv(results_table, "multinomial_results.csv", row.names = FALSE)

############################################################
# MODEL FIT — McFadden pseudo-R2
############################################################

null_model <- multinom(cluster ~ 1, data = model_data, trace = FALSE)

lrt <- lrtest(null_model, multi_model)
cat("\nLikelihood Ratio Test (null vs main model):\n"); print(lrt)

pseudo_r2 <- 1 - (as.numeric(logLik(multi_model)) / as.numeric(logLik(null_model)))
cat("McFadden pseudo-R2:", round(pseudo_r2, 4), "\n")

############################################################
# YEAR-INTERACTION DATASET
# Controls excluded — they have no breakdown year.
# Neighbour NAs filled with 0 for consistency.
############################################################

model_data_yr <- final_data %>%
  filter(cluster != "0") %>%
  select(
    cluster, bd_start_yr,
    Control_moves, NewBD_moves, Burnout_moves, Unknown_moves,
    Burnout_neighbours, Control_neighbours, NewBD_neighbours
  ) %>%
  mutate(across(
    c(Burnout_neighbours, Control_neighbours, NewBD_neighbours),
    ~replace_na(., 0)
  )) %>%
  na.omit()

model_data_yr$cluster <- droplevels(model_data_yr$cluster)
model_data_yr$cluster <- relevel(model_data_yr$cluster, ref = "1")

cat("\nmodel_data_yr cluster table (clusters 1-3 only):\n")
print(table(model_data_yr$cluster))

############################################################
# YEAR INTERACTION MODEL
############################################################

multi_model_yr <- multinom(
  cluster ~
    NewBD_neighbours * bd_start_yr +
    Burnout_neighbours * bd_start_yr +
    Control_moves + NewBD_moves + Burnout_moves +
    Unknown_moves + Control_neighbours,
  data = model_data_yr, trace = FALSE
)

or_yr <- exp(coef(multi_model_yr))
z_yr  <- summary(multi_model_yr)$coefficients / summary(multi_model_yr)$standard.errors
p_yr  <- 2 * (1 - pnorm(abs(z_yr)))

cat("\n--- YEAR INTERACTION MODEL (clusters 1-3 only) ---\n")
print(summary(multi_model_yr))
cat("Odds ratios:\n"); print(or_yr)
cat("P-values:\n");    print(p_yr)

write.csv(as.data.frame(or_yr), "multinomial_year_interaction_OR.csv", row.names = TRUE)

cat("\nNote: multi_model and multi_model_yr answer different questions\n")
cat("on different populations so are reported separately, not compared\n")
cat("via lrtest.\n")

############################################################
# PREDICTED PROBABILITIES
############################################################

predicted_probs <- fitted(multi_model)
cat("Average predicted probabilities:\n")
print(colMeans(predicted_probs))

############################################################
# MARGINAL EFFECTS
############################################################

base_rf <- data.frame(
  Control_moves      = mean(model_data$Control_moves,      na.rm=TRUE),
  NewBD_moves        = mean(model_data$NewBD_moves,        na.rm=TRUE),
  Burnout_moves      = mean(model_data$Burnout_moves,      na.rm=TRUE),
  Unknown_moves      = mean(model_data$Unknown_moves,      na.rm=TRUE),
  Burnout_neighbours = mean(model_data$Burnout_neighbours, na.rm=TRUE),
  Control_neighbours = mean(model_data$Control_neighbours, na.rm=TRUE),
  NewBD_neighbours   = mean(model_data$NewBD_neighbours,   na.rm=TRUE)
)

rf_predictors <- c("Control_moves","NewBD_moves","Burnout_moves","Unknown_moves",
                   "Burnout_neighbours","Control_neighbours","NewBD_neighbours")

rf_effects <- bind_rows(
  lapply(rf_predictors, function(pred) {
    vals  <- quantile(model_data[[pred]], probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
    td    <- base_rf[rep(1, 3), ]
    td[[pred]] <- vals
    probs <- predict(multi_model, newdata = td, type = "probs")
    data.frame(Predictor = pred, Level = c("10th","50th","90th"), as.data.frame(probs))
  })
)

print(rf_effects)
write.csv(rf_effects, "predictor_effects.csv", row.names = FALSE)

# R auto-prefixes numeric column names with "X" (e.g. "0" becomes "X0")
cat("rf_effects column names before rename:\n")
print(names(rf_effects))

names(rf_effects)[names(rf_effects) == "X0"] <- "Control"
names(rf_effects)[names(rf_effects) == "X1"] <- "Cluster 1"
names(rf_effects)[names(rf_effects) == "X2"] <- "Cluster 2"
names(rf_effects)[names(rf_effects) == "X3"] <- "Cluster 3"

cat("rf_effects column names after rename:\n")
print(names(rf_effects))

plot_long <- rf_effects %>%
  pivot_longer(
    cols      = c("Control", "Cluster 1", "Cluster 2", "Cluster 3"),
    names_to  = "Cluster",
    values_to = "Probability"
  )

p_effects <- ggplot(plot_long, aes(x = Level, y = Probability, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Predictor) +
  labs(title = "Effect of Risk Factor Predictors on Cluster Membership",
       x = "Predictor Percentile", y = "Predicted Probability") +
  theme_minimal()

print(p_effects)
ggsave("predictor_effects_barplot.png", plot = p_effects, width = 14, height = 10, dpi = 300)

cat("\n=== MAIN PIPELINE COMPLETE ===\n")
cat("Objects in memory for Andrew feedback script:\n")
cat("bd_scaled_full, bd_core, cluster_vars_log, pam3_canonical,\n")
cat("model_data, multi_model, multi_model_yr, null_model,\n")
cat("data_bd, multi_data, control_herds\n")

# Sections
#   A. Cluster stability – clusterboot (B=100)
#   B. Cluster stability – split-sample validation
#   C. Cluster stability – sequential variable removal
#   D. Multiple-seed control sampling robustness
#   E. Multinomial discrimination – multiclass AUC
#   F. Cox proportional hazards model
#   G. Fuzzy C-means (exploratory)
#   H. Methods text snippet for K-means description
#   I. Consensus Clustering (Monti et al. 2003) — distinct from
#      clusterboot; produces a pairwise co-occurrence matrix and
#      heatmap, the specific method Andrew named explicitly


options(bitmapType = "cairo")
setwd("C:/Users/Victory/Downloads")

library(tidyverse)
library(cluster)
library(fpc)
library(factoextra)
library(nnet)
library(survival)
library(survminer)
library(gridExtra)
library(lmtest)
library(pROC)
library(ppclust)
library(pheatmap)



required_objects <- c(
  "bd_scaled_full", "bd_core", "cluster_vars_log", "pam3_canonical",
  "model_data", "multi_model", "null_model", "data_bd", "multi_data",
  "control_herds"
)

missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects from main script: ",
    paste(missing_objects, collapse = ", "),
    ". Run btb_main_pipeline_v3.R first in this session."
  )
}


# A. BOOTSTRAP CLUSTER STABILITY (clusterboot, B = 200)
# Jaccard > 0.85 = stable; 0.60-0.85 = moderately stable


cat("\n============================================================\n")
cat("A. BOOTSTRAP CLUSTER STABILITY (clusterboot, B = 200)\n")
cat("============================================================\n")

set.seed(123)
cb_pam3 <- clusterboot(
  bd_scaled_full,
  B             = 200,
  clustermethod = pamkCBI,
  krange        = 3,
  seed          = 123,
  multipleboot  = FALSE,
  showplots     = FALSE
)

cat("Jaccard stability per cluster (mean over 200 bootstraps):\n")
print(cb_pam3$bootmean)

stability_results <- data.frame(
  Cluster         = 1:3,
  Jaccard_Mean    = cb_pam3$bootmean,
  Times_Dissolved = cb_pam3$bootbrd,
  Times_Recovered = 200 - cb_pam3$bootbrd,
  Stability       = ifelse(
    cb_pam3$bootmean > 0.85, "Stable",
    ifelse(cb_pam3$bootmean > 0.60, "Moderately stable", "Unstable")
  )
)

print(stability_results)
write.csv(stability_results, "A_cluster_stability_bootstrap.csv", row.names = FALSE)

png("A_jaccard_stability.png", 1000, 700)
barplot(
  cb_pam3$bootmean,
  names.arg = paste("Cluster", 1:3),
  ylim = c(0, 1),
  col  = ifelse(cb_pam3$bootmean > 0.85, "steelblue",
                ifelse(cb_pam3$bootmean > 0.60, "orange", "red")),
  ylab = "Mean Jaccard Coefficient",
  main = "PAM K=3 Bootstrap Stability (B = 200)",
  las  = 1
)
abline(h = 0.85, lty = 2, col = "darkgreen", lwd = 2)
abline(h = 0.60, lty = 2, col = "red", lwd = 2)
legend("topright",
       legend = c("Stable threshold (0.85)", "Minimum threshold (0.60)"),
       lty = 2, col = c("darkgreen", "red"), lwd = 2)
dev.off()


# B. SPLIT-SAMPLE VALIDATION (5 × 70/30 splits)

cat("\n============================================================\n")
cat("B. SPLIT-SAMPLE VALIDATION (5 × 70/30 splits)\n")
cat("============================================================\n")

ari <- function(true, pred) mclust::adjustedRandIndex(true, pred)

set.seed(123)
n_obs       <- nrow(bd_scaled_full)
n_splits    <- 5
ari_results <- numeric(n_splits)

for (i in seq_len(n_splits)) {
  train_idx <- sample(n_obs, size = floor(0.70 * n_obs))
  test_idx  <- setdiff(seq_len(n_obs), train_idx)
  
  train_data <- bd_scaled_full[train_idx, ]
  test_data  <- bd_scaled_full[test_idx, ]
  
  pam_train <- pam(train_data, k = 3)
  
  medoids     <- train_data[pam_train$id.med, ]
  dists       <- as.matrix(dist(rbind(medoids, test_data)))
  n_med       <- nrow(medoids)
  n_test      <- nrow(test_data)
  dist_to_med <- dists[(n_med + 1):(n_med + n_test), 1:n_med]
  pred_cluster <- apply(dist_to_med, 1, which.min)
  
  true_cluster <- pam3_canonical$clustering[test_idx]
  
  ari_results[i] <- ari(true_cluster, pred_cluster)
  cat(sprintf("  Split %d  ARI = %.4f\n", i, ari_results[i]))
}

split_summary <- data.frame(Split = 1:n_splits, ARI = round(ari_results, 4))
split_summary$Interpretation <- ifelse(
  split_summary$ARI > 0.70, "Strong agreement",
  ifelse(split_summary$ARI > 0.40, "Moderate agreement", "Weak agreement")
)

cat("\nSplit-sample ARI summary:\n")
print(split_summary)
cat(sprintf("Mean ARI: %.4f  (SD: %.4f)\n", mean(ari_results), sd(ari_results)))

write.csv(split_summary, "B_split_sample_validation.csv", row.names = FALSE)


# C. SEQUENTIAL VARIABLE REMOVAL

cat("\n============================================================\n")
cat("C. SEQUENTIAL VARIABLE REMOVAL\n")
cat("============================================================\n")

var_removal_results <- data.frame(
  Variable_Removed = character(),
  ARI_vs_canonical = numeric(),
  stringsAsFactors = FALSE
)

for (drop_var in cluster_vars_log) {
  
  remaining_vars <- setdiff(cluster_vars_log, drop_var)
  
  reduced_data <- bd_core %>%
    select(all_of(remaining_vars)) %>%
    scale() %>%
    as.data.frame()
  
  set.seed(123)
  pam_reduced <- pam(reduced_data, k = 3)
  
  ari_val <- ari(pam3_canonical$clustering, pam_reduced$clustering)
  cat(sprintf("  Removed: %-35s  ARI = %.4f\n", drop_var, ari_val))
  
  var_removal_results <- rbind(
    var_removal_results,
    data.frame(Variable_Removed = drop_var, ARI_vs_canonical = round(ari_val, 4),
               stringsAsFactors = FALSE)
  )
}

var_removal_results$Interpretation <- ifelse(
  var_removal_results$ARI_vs_canonical > 0.80, "Robust (variable not dominant)",
  ifelse(var_removal_results$ARI_vs_canonical > 0.50, "Moderate change",
         "Large change (variable influential)")
)

cat("\nVariable removal summary:\n")
print(var_removal_results)
write.csv(var_removal_results, "C_variable_removal_stability.csv", row.names = FALSE)

png("C_variable_removal_ari.png", 1200, 800)
par(mar = c(10, 5, 4, 2))
barplot(
  var_removal_results$ARI_vs_canonical,
  names.arg = var_removal_results$Variable_Removed,
  las = 2, ylim = c(0, 1), col = "steelblue",
  ylab = "ARI vs Canonical PAM K=3",
  main = "Sequential Variable Removal – Cluster Stability"
)
abline(h = 0.80, lty = 2, col = "darkgreen", lwd = 2)
abline(h = 0.50, lty = 2, col = "red", lwd = 2)
dev.off()


# D. MULTIPLE-SEED CONTROL SAMPLING ROBUSTNESS


cat("\n============================================================\n")
cat("D. MULTIPLE-SEED CONTROL SAMPLING ROBUSTNESS\n")
cat("============================================================\n")

all_controls_pool <- multi_data %>%
  filter(herd_no %in% control_herds$herd_no)

seeds        <- c(1, 2, 3)
seed_or_list <- list()

# data_bd cluster must be character for safe binding (do this once)
data_bd_char <- data_bd
data_bd_char$cluster <- as.character(data_bd_char$cluster)

for (s in seeds) {
  
  set.seed(s)
  controls_s <- sample_n(all_controls_pool, min(30000, nrow(all_controls_pool)))
  controls_s$cluster <- "0"   # character, matches data_bd_char$cluster type
  
  fd_s <- bind_rows(data_bd_char, controls_s)
  fd_s$cluster <- as.factor(fd_s$cluster)
  fd_s$cluster <- relevel(fd_s$cluster, ref = "0")
  
  md_s <- fd_s %>%
    select(cluster, Control_moves, NewBD_moves, Burnout_moves,
           Unknown_moves, Burnout_neighbours, Control_neighbours,
           NewBD_neighbours) %>%
    na.omit()
  
  m_s <- multinom(
    cluster ~ Control_moves + NewBD_moves + Burnout_moves +
      Unknown_moves + Burnout_neighbours + Control_neighbours +
      NewBD_neighbours,
    data  = md_s,
    trace = FALSE
  )
  
  or_s <- exp(coef(m_s))
  z_s  <- summary(m_s)$coefficients / summary(m_s)$standard.errors
  p_s  <- 2 * (1 - pnorm(abs(z_s)))
  
  seed_or_list[[paste0("seed_", s)]] <- data.frame(
    Seed     = s,
    Variable = colnames(summary(m_s)$coefficients),
    C1_OR    = as.numeric(or_s["1", ]), C1_p = as.numeric(p_s["1", ]),
    C2_OR    = as.numeric(or_s["2", ]), C2_p = as.numeric(p_s["2", ]),
    C3_OR    = as.numeric(or_s["3", ]), C3_p = as.numeric(p_s["3", ])
  )
  
  cat(sprintf("\n  Seed %d  – Cluster sizes: %s\n",
              s, paste(table(md_s$cluster), collapse = " / ")))
}

seed_or_df <- bind_rows(seed_or_list)
write.csv(seed_or_df, "D_multiseed_control_OR_comparison.csv", row.names = FALSE)

or_range <- seed_or_df %>%
  group_by(Variable) %>%
  summarise(
    C1_OR_min = min(C1_OR), C1_OR_max = max(C1_OR),
    C2_OR_min = min(C2_OR), C2_OR_max = max(C2_OR),
    C3_OR_min = min(C3_OR), C3_OR_max = max(C3_OR),
    .groups = "drop"
  )

cat("\nOR range across seeds:\n")
print(or_range)
write.csv(or_range, "D_multiseed_OR_range.csv", row.names = FALSE)

newbd_rows <- seed_or_df %>%
  filter(Variable == "NewBD_neighbours") %>%
  select(Seed, C1_OR, C2_OR, C3_OR) %>%
  pivot_longer(cols = -Seed, names_to = "Cluster", values_to = "OR")

png("D_multiseed_NewBD_neighbours_OR.png", 1000, 700)
ggplot(newbd_rows, aes(x = factor(Seed), y = OR, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  labs(title = "NewBD_neighbours OR across control sampling seeds",
       x = "Random seed", y = "Odds Ratio")
dev.off()


# E. MULTINOMIAL DISCRIMINATION – MULTICLASS AUC


cat("\n============================================================\n")
cat("E. MULTINOMIAL DISCRIMINATION – MULTICLASS AUC\n")
cat("============================================================\n")

pred_probs <- as.data.frame(fitted(multi_model))
true_class <- model_data$cluster

cat("\nOne-vs-rest AUC per class:\n")
ovr_aucs <- sapply(levels(true_class), function(cls) {
  binary_truth <- as.numeric(true_class == cls)
  auc_val <- pROC::auc(pROC::roc(binary_truth, pred_probs[[cls]], quiet = TRUE))
  cat(sprintf("  Class %s:  AUC = %.4f\n", cls, auc_val))
  auc_val
})

ovr_df <- data.frame(Class = levels(true_class), AUC = round(ovr_aucs, 4))
write.csv(ovr_df, "E_one_vs_rest_AUC.csv", row.names = FALSE)

cat("\nHand-Till multiclass AUC:\n")
mc_roc <- pROC::multiclass.roc(response = true_class, predictor = pred_probs)
cat(sprintf("  Multiclass AUC = %.4f\n", as.numeric(mc_roc$auc)))

multiclass_auc_df <- data.frame(
  Metric = "Hand-Till multiclass AUC",
  AUC    = round(as.numeric(mc_roc$auc), 4),
  McFadden_pseudo_R2 = round(
    1 - (as.numeric(logLik(multi_model)) / as.numeric(logLik(null_model))), 4
  )
)
write.csv(multiclass_auc_df, "E_multiclass_AUC.csv", row.names = FALSE)

png("E_one_vs_rest_ROC.png", 1400, 1000)
par(mfrow = c(2, 2))
for (cls in levels(true_class)) {
  binary_truth <- as.numeric(true_class == cls)
  roc_obj <- pROC::roc(binary_truth, pred_probs[[cls]], quiet = TRUE)
  plot(roc_obj,
       main = sprintf("Class %s vs Rest  (AUC = %.3f)", cls, pROC::auc(roc_obj)),
       col = "steelblue", lwd = 2)
}
dev.off()


# F. COX PROPORTIONAL HAZARDS MODELS

cat("\n============================================================\n")
cat("F. COX PROPORTIONAL HAZARDS MODELS\n")
cat("============================================================\n")

# F1. Breakdown duration
cox_dur_data <- bd_core %>%
  select(bd_duration_days, cluster, herd_size_mean,
         no_bds_started_in_prev_5yr, bd_within_13_months) %>%
  drop_na()

cox_dur_data$cluster <- relevel(cox_dur_data$cluster, ref = "1")

cox_fit_dur <- coxph(
  Surv(bd_duration_days, rep(1, nrow(cox_dur_data))) ~
    cluster + log1p(herd_size_mean) + log1p(no_bds_started_in_prev_5yr) +
    bd_within_13_months,
  data = cox_dur_data
)

cat("\nCox model – breakdown duration:\n")
print(summary(cox_fit_dur))

cox_dur_summary <- as.data.frame(summary(cox_fit_dur)$conf.int)
cox_dur_summary$Variable <- rownames(cox_dur_summary)
cox_dur_summary$p_value  <- summary(cox_fit_dur)$coefficients[, "Pr(>|z|)"]
write.csv(cox_dur_summary, "F1_cox_duration_HR.csv", row.names = FALSE)

png("F1_cox_duration_schoenfeld.png", 1200, 900)
cox.zph_dur <- cox.zph(cox_fit_dur)
plot(cox.zph_dur, main = "Schoenfeld Residuals – Breakdown Duration")
dev.off()

cat("\nProportional hazards test:\n"); print(cox.zph_dur)

# F2. Time to next breakdown
cox_next_data <- bd_core %>%
  select(duration_between_bd, cluster, herd_size_mean,
         no_bds_started_in_prev_5yr, bd_within_13_months) %>%
  mutate(
    time  = ifelse(is.na(duration_between_bd), 1000, duration_between_bd),
    event = ifelse(is.na(duration_between_bd), 0, 1)
  ) %>%
  drop_na(cluster, herd_size_mean)

cox_next_data$cluster <- relevel(cox_next_data$cluster, ref = "1")

cox_fit_next <- coxph(
  Surv(time, event) ~
    cluster + log1p(herd_size_mean) + log1p(no_bds_started_in_prev_5yr) +
    bd_within_13_months,
  data = cox_next_data
)

cat("\nCox model – time to next breakdown:\n")
print(summary(cox_fit_next))

cox_next_summary <- as.data.frame(summary(cox_fit_next)$conf.int)
cox_next_summary$Variable <- rownames(cox_next_summary)
cox_next_summary$p_value  <- summary(cox_fit_next)$coefficients[, "Pr(>|z|)"]
write.csv(cox_next_summary, "F2_cox_next_breakdown_HR.csv", row.names = FALSE)

png("F2_cox_next_schoenfeld.png", 1200, 900)
cox.zph_next <- cox.zph(cox_fit_next)
plot(cox.zph_next, main = "Schoenfeld Residuals – Time to Next Breakdown")
dev.off()

cat("\nProportional hazards test:\n"); print(cox.zph_next)

hr_df <- cox_next_summary %>%
  rename(HR = `exp(coef)`, Lower = `lower .95`, Upper = `upper .95`) %>%
  mutate(Variable = factor(Variable, levels = rev(Variable)))

png("F2_cox_next_forest.png", 1200, 800)
ggplot(hr_df, aes(x = HR, y = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
  scale_x_log10() +
  theme_bw() +
  labs(title = "Hazard Ratios – Time to Next Breakdown (Cox model)",
       x = "Hazard Ratio (log scale)", y = "")
dev.off()


# G. FUZZY C-MEANS (exploratory)


cat("\n============================================================\n")
cat("G. FUZZY C-MEANS (exploratory)\n")
cat("============================================================\n")

set.seed(123)
fcm_result <- ppclust::fcm(as.matrix(bd_scaled_full), centers = 3, nstart = 5)

cat("Fuzzy membership summary (first 6 rows):\n")
print(head(round(fcm_result$u, 3)))

fcm_hard <- apply(fcm_result$u, 1, which.max)
fcm_ari  <- mclust::adjustedRandIndex(pam3_canonical$clustering, fcm_hard)
cat(sprintf("\nARI between PAM K=3 and Fuzzy C-means (hard): %.4f\n", fcm_ari))

max_membership <- apply(fcm_result$u, 1, max)
cat("\nDistribution of maximum membership probability:\n")
print(quantile(max_membership, probs = c(0.10, 0.25, 0.50, 0.75, 0.90)))

png("G_fuzzy_membership_histogram.png", 1000, 700)
hist(max_membership, breaks = 40, col = "steelblue",
     xlab = "Maximum cluster membership probability",
     main = "Fuzzy C-means: Confidence of cluster assignment", xlim = c(0, 1))
abline(v = 0.80, lty = 2, col = "red", lwd = 2)
legend("topleft", legend = "0.80 threshold", lty = 2, col = "red", lwd = 2)
dev.off()

fcm_u <- as.data.frame(fcm_result$u)
names(fcm_u) <- paste0("Cluster_", 1:3)

png("G_fuzzy_membership_scatter.png", 1200, 900)
ggplot(fcm_u, aes(x = Cluster_1, y = Cluster_2, colour = Cluster_3)) +
  geom_point(alpha = 0.2, size = 0.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme_bw() +
  labs(title = "Fuzzy C-means membership: Cluster 1 vs Cluster 2",
       x = "Membership probability – Cluster 1",
       y = "Membership probability – Cluster 2",
       colour = "Cluster 3\nmembership")
dev.off()

fcm_export <- bd_core %>%
  select(herd_no, bd_no) %>%
  bind_cols(fcm_u) %>%
  mutate(
    PAM_cluster    = pam3_canonical$clustering,
    FCM_hard       = fcm_hard,
    Max_membership = max_membership
  )

write.csv(fcm_export, "G_fuzzy_membership_assignments.csv", row.names = FALSE)

############################################################
# H. METHODS TEXT SNIPPET – K-MEANS DESCRIPTION
############################################################

cat("\n============================================================\n")
cat("H. METHODS TEXT SNIPPET \n")
cat("============================================================\n")

methods_text <- "
K-means clustering was applied as a comparator method to evaluate the
consistency of the cluster solution. In k-means clustering, k cluster
centroids are initialised and iteratively updated by assigning each
observation to the nearest centroid and recalculating centroids as the
mean of assigned observations, continuing until within-cluster
sum-of-squares (WSS) converges. To avoid local minima, 25 random
restarts (nstart = 25) were used and the solution with lowest total
WSS was retained. K-means was applied for k = 2, 3, 4, and 5 on the
same scaled matrix used for PAM. Cluster number was informed by the
elbow method (WSS vs k), the average silhouette width, the gap
statistic, and majority vote across 26 indices using NbClust.
Agreement between k-means and PAM cluster assignments was quantified
using the Adjusted Rand Index (ARI), where ARI = 1 indicates perfect
agreement. The final cluster solution was based on PAM, which is more
robust to outliers than k-means, as PAM uses actual observations
(medoids) rather than means as cluster centres.
"

cat(methods_text)
writeLines(methods_text, "H_methods_snippet_kmeans.txt")

############################################################
# I. CONSENSUS CLUSTERING (Monti et al. 2003)
# ------------------------------------------------------------
# It was explicitly stated that Consensus Clustering as a method to
# "kick the tyres" on the K=3 PAM solution. This is DIFFERENT
# from clusterboot (Section A): clusterboot checks whether the
# SAME CLUSTERS reform after resampling; consensus clustering
# checks, for every PAIR of herds, how often they land in the
# same cluster together across many resamples — producing a
# full pairwise co-occurrence matrix and a visual heatmap.
#
# Algorithm (manual implementation, same logic as
# ConsensusClusterPlus / Monti et al. 2003):
#   1. Resample 80% of herds (no replacement), B = 200 times
#   2. Run PAM K=3 on each resample
#   3. Build an N x N consensus matrix:
#        M[i,j] = (# resamples where i & j were both drawn AND
#                   placed in the same cluster) /
#                  (# resamples where i & j were both drawn)
#   4. Reorder the matrix by hierarchical clustering and plot
#      as a heatmap. Sharp, well-separated blocks = stable,
#      genuine clusters. Smeared/fuzzy blocks = weak structure.
#   5. Compute the consensus CDF and area-under-CDF (a single
#      number summarising overall stability: closer to 1 = more
#      block-like = more stable).
#
# NOTE ON SCALE: consensus clustering is O(n^2) in memory and
# compute for the matrix step, so this runs on a 2,000-row
# subsample, consistent with how other expensive stability
# checks (VAT, Hopkins) are handled elsewhere in this pipeline.
# This is documented, not hidden.
############################################################

cat("\n============================================================\n")
cat("I. CONSENSUS CLUSTERING (Monti et al. 2003)\n")
cat("============================================================\n")

set.seed(123)
cc_idx  <- sample(nrow(bd_scaled_full), min(2000, nrow(bd_scaled_full)))
cc_data <- bd_scaled_full[cc_idx, ]
n_cc    <- nrow(cc_data)

B_cc          <- 200   # number of resampling iterations
resample_frac <- 0.80  # proportion of herds drawn each iteration

# connectivity_count[i,j]: number of times i & j were BOTH sampled
#                          AND placed in the same cluster
# total_count[i,j]:        number of times i & j were BOTH sampled
connectivity_count <- matrix(0, n_cc, n_cc)
total_count         <- matrix(0, n_cc, n_cc)

cat(sprintf("Running %d resampling iterations on n = %d herds...\n", B_cc, n_cc))

set.seed(123)
for (b in seq_len(B_cc)) {
  
  sampled_idx <- sample(seq_len(n_cc), size = floor(resample_frac * n_cc))
  
  pam_b <- pam(cc_data[sampled_idx, ], k = 3)
  cluster_b <- pam_b$clustering
  names(cluster_b) <- sampled_idx
  
  # Update total_count: all pairs within this resample were "both sampled"
  total_count[sampled_idx, sampled_idx] <- total_count[sampled_idx, sampled_idx] + 1
  
  # Update connectivity_count: pairs placed in the SAME cluster this round
  for (k in 1:3) {
    members_k <- sampled_idx[cluster_b == k]
    if (length(members_k) > 1) {
      connectivity_count[members_k, members_k] <-
        connectivity_count[members_k, members_k] + 1
    }
  }
  
  if (b %% 50 == 0) cat(sprintf("  ...completed %d / %d iterations\n", b, B_cc))
}

# Consensus matrix: proportion of co-sampled iterations where pair
# landed in the same cluster. Diagonal forced to 1 (self-agreement).
consensus_matrix <- connectivity_count / total_count
consensus_matrix[total_count == 0] <- NA
diag(consensus_matrix) <- 1

cat("\nConsensus matrix built. Summary of off-diagonal values:\n")
print(summary(consensus_matrix[upper.tri(consensus_matrix)]))

############################################################
# CONSENSUS HEATMAP — reordered by hierarchical clustering
# of (1 - consensus) as a distance measure
############################################################

consensus_dist <- as.dist(1 - consensus_matrix)
consensus_hc   <- hclust(consensus_dist, method = "average")

# Also get a K=3 grouping FROM the consensus matrix itself, for
# comparison against the canonical PAM labels
consensus_groups <- cutree(consensus_hc, k = 3)

png("I_consensus_heatmap.png", 1200, 1200)
pheatmap(
  consensus_matrix,
  cluster_rows   = consensus_hc,
  cluster_cols   = consensus_hc,
  show_rownames  = FALSE,
  show_colnames  = FALSE,
  color          = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
  main           = sprintf("Consensus Clustering Matrix (K=3, B=%d resamples, n=%d)",
                           B_cc, n_cc)
)
dev.off()

############################################################
# CONSENSUS CDF + AREA UNDER CDF
# A single summary number: higher = more block-like = more
# stable clustering structure. This is the standard Monti et al.
# stability metric used alongside the heatmap.
############################################################

off_diag_vals <- consensus_matrix[upper.tri(consensus_matrix)]
off_diag_vals <- off_diag_vals[!is.na(off_diag_vals)]

consensus_cdf <- ecdf(off_diag_vals)
cdf_grid      <- seq(0, 1, by = 0.01)
cdf_values    <- consensus_cdf(cdf_grid)

# Area under the CDF curve (trapezoidal integration)
auc_cdf <- sum(diff(cdf_grid) * (head(cdf_values, -1) + tail(cdf_values, -1)) / 2)

cat(sprintf("\nConsensus CDF area-under-curve: %.4f\n", auc_cdf))
cat("(Closer to 1 = sharper bimodal consensus = more stable clusters;\n")
cat(" closer to 0.5 = consensus values spread uniformly = weak/ambiguous structure)\n")

png("I_consensus_cdf.png", 1000, 700)
plot(
  cdf_grid, cdf_values, type = "l", lwd = 2, col = "steelblue",
  xlab = "Consensus index value",
  ylab = "CDF",
  main = sprintf("Consensus CDF (Area under curve = %.3f)", auc_cdf)
)
abline(v = 0.5, lty = 2, col = "grey50")
dev.off()

############################################################
# AGREEMENT: Consensus-derived groups vs canonical PAM
############################################################

# Map canonical PAM labels for the same subsample
canonical_labels_subsample <- pam3_canonical$clustering[cc_idx]

consensus_vs_canonical_ari <- mclust::adjustedRandIndex(
  canonical_labels_subsample,
  consensus_groups
)

cat(sprintf(
  "\nARI between canonical PAM K=3 and consensus-derived groups: %.4f\n",
  consensus_vs_canonical_ari
))

############################################################
# SAVE RESULTS
############################################################

consensus_summary <- data.frame(
  Metric = c(
    "Subsample size (n)",
    "Resampling iterations (B)",
    "Resample fraction",
    "Consensus CDF area-under-curve",
    "ARI: consensus groups vs canonical PAM K=3"
  ),
  Value = c(
    n_cc,
    B_cc,
    resample_frac,
    round(auc_cdf, 4),
    round(consensus_vs_canonical_ari, 4)
  )
)

print(consensus_summary)
write.csv(consensus_summary, "I_consensus_clustering_summary.csv", row.names = FALSE)

cat("\n=== SECTION I (CONSENSUS CLUSTERING) COMPLETE ===\n")

#VALIDATION SUMMARY

cat("\n============================================================\n")
cat("VALIDATION SUMMARY\n")
cat("============================================================\n")

summary_bootstrap <- data.frame(
  Method = paste0("Bootstrap Jaccard – Cluster ", 1:3),
  Result = round(cb_pam3$bootmean, 4)
)

summary_split <- data.frame(
  Method = paste0("Split-sample ARI – Split ", 1:n_splits),
  Result = round(ari_results, 4)
)
summary_split <- rbind(
  summary_split,
  data.frame(Method = "Split-sample ARI – Mean", Result = round(mean(ari_results), 4))
)

summary_varremoval <- data.frame(
  Method = paste0("Variable removal ARI – remove ", var_removal_results$Variable_Removed),
  Result = var_removal_results$ARI_vs_canonical
)

summary_other <- data.frame(
  Method = c(
    "Fuzzy C-means ARI vs PAM (hard assignment)",
    "Multiclass AUC (Hand-Till)",
    "McFadden pseudo-R2 (multinomial)",
    "Consensus clustering CDF area-under-curve",
    "Consensus clustering ARI vs canonical PAM K=3"
  ),
  Result = c(
    round(fcm_ari, 4),
    round(as.numeric(mc_roc$auc), 4),
    round(1 - (as.numeric(logLik(multi_model)) / as.numeric(logLik(null_model))), 4),
    round(auc_cdf, 4),
    round(consensus_vs_canonical_ari, 4)
  )
)

validation_summary <- bind_rows(
  summary_bootstrap, summary_split, summary_varremoval, summary_other
)

print(validation_summary)
write.csv(validation_summary, "VALIDATION_SUMMARY.csv", row.names = FALSE)

cat("\n===COMPLETED ===\n")
cat("Output files written to:", getwd(), "\n")
