# Modelling Bovine Tuberculosis Spread Through Cattle Trade and Neighbourhood Networks

> **PhD Research — University College Dublin | VistaMilk Research Ireland**
> *Stochastic metapopulation modelling · Unsupervised machine learning · Network epidemiology · One Health*

---

## Overview

Bovine tuberculosis (bTB), caused by *Mycobacterium bovis*, remains one of the most persistent and economically devastating infectious diseases affecting cattle in Ireland  costing approximately €84 million annually despite decades of statutory testing, movement restriction, and compulsory slaughter programmes.

This repository contains the full analytical and simulation pipeline underpinning my PhD research into **why bTB breakdowns are not all the same**  and what that heterogeneity means for disease control.

The core contribution of this work is the discovery, through unsupervised machine learning applied to over **2.8 million surveillance observations**, that Irish cattle herds experience bTB through three fundamentally distinct outbreak phenotypes. These phenotypes  self-limiting, chronically recurrent, and acutely explosive  reflect qualitatively different transmission mechanisms, and demand qualitatively different intervention responses.

This repository provides:
- The complete **data processing and feature engineering pipeline**
- **Unsupervised clustering** (PAM, k-means, GMM, Fuzzy C-means, Consensus Clustering)
- **Network-based transmission modelling** incorporating trade movement and neighbourhood structure
- **Stochastic metapopulation simulation** via continuous-time Markov chains (CTMC)
- **Multinomial regression** and **survival analysis** for cluster risk factor characterisation
- The mathematical foundation for **cluster-recovery validation** asking whether stochastic simulations can reproduce machine-learning-discovered structure from first principles

---

## The Central Research Question

*Can stochastic simulation models, parameterised from epidemiological surveillance data and grounded in network theory, reproduce the three outbreak phenotypes identified through unsupervised machine learning  and are those simulated clusters consistent with the empirically observed clusters?*

A positive answer provides strong evidence that the discovered phenotypes reflect genuine biological structure rather than statistical artefact. A negative answer identifies precisely which transmission mechanism the model has misspecified. Either way, the result advances understanding.

---

## Repository Structure

```
├── data/
│   ├── preprocessing/          # Data cleaning, linkage, exclusion criteria
│   ├── feature_engineering/    # Clustering variable construction and transformation
│   └── README.md               # Data dictionary and variable definitions
│
├── clustering/
│   ├── tendency/               # Hopkins statistic, VAT plot
│   ├── pam/                    # Primary PAM clustering (K=3)
│   ├── kmeans/                 # K-means comparator
│   ├── gmm/                    # Gaussian Mixture Model
│   ├── fuzzy/                  # Fuzzy C-means
│   ├── consensus/              # Consensus clustering (B=200 resamples)
│   └── validation/             # Bootstrap stability, split-sample ARI
│
├── pca/
│   └── pca_analysis.R          # Dimensionality reduction and visualisation
│
├── networks/
│   ├── trade/                  # Cattle movement network construction (AIMS data)
│   ├── neighbourhood/          # Spatial neighbourhood network
│   └── metapopulation/         # Combined network transmission model
│
├── simulation/
│   ├── ode/                    # Deterministic ODE system (deSolve)
│   ├── ctmc/                   # Stochastic CTMC tau-leap (SimInf-compatible)
│   └── validation/             # ODE vs CTMC ensemble mean comparison
│
├── regression/
│   └── multinomial/            # Multinomial logistic regression (nnet)
│
├── survival/
│   ├── kaplan_meier/           # KM curves for breakdown duration and recurrence
│   └── cox/                    # Cox proportional hazards models
│
├── figures/                    # All publication-quality figures
└── manuscript/                 # LaTeX source for submitted manuscript
```

---

## The Three Outbreak Phenotypes

| | **Cluster 1** | **Cluster 2** | **Cluster 3** |
|---|---|---|---|
| **Label** | Low-risk, self-limiting | Chronic recurrence | Acute high-intensity |
| **Share of breakdowns** | 65.8% | 18.3% | 16.5% |
| **Mean spread ratio** | 0.004 | 0.073 | 1.975 |
| **Mean duration (days)** | 149 | 177 | 333 |
| **Prior breakdowns (5yr)** | 0.000 | 1.109 | 0.126 |
| **Mechanism** | Early detection, rapid resolution | Residual infection or persistent wildlife/environmental exposure | Naïve herd encountering intense infection pressure |
| **Control implication** | Standard test-and-removal | Intensified environmental and wildlife investigation | Rapid response to limit within-herd amplification |

---

## Model Architecture

### Force of Infection

$$\lambda_j(t) = \beta \cdot T_j(t) \cdot V_j(t) \cdot \left[\frac{1}{N} \sum_{i \neq j} I_i(t)\right] + \omega W(t)$$

where $T_j(t)$ is the trade-movement intensity for herd $j$, $V_j(t)$ is within-herd prevalence, the summation term captures metapopulation infection pressure from neighbouring herds, and $W(t)$ is the wildlife reservoir index.

### Compartmental Structure

Each herd $j$ transitions through seven epidemiological states:

$$S_j \xrightarrow{p_K \lambda_j} I_{Kj} \xrightarrow{\gamma_K} D_{Kj} \xrightarrow{\mu_K} S_j \quad \text{or} \quad D_{Kj} \xrightarrow{\phi_K} I_{Kj}$$

for $K \in \{A, B, C\}$, where $p_A + p_B + p_C = 1$ and the relapse pathway $D_{Kj} \rightarrow I_{Kj}$ captures the biological reality that detection does not guarantee complete clearance.

### ODE System

$$\frac{dS_j}{dt} = -\lambda_j S_j + \mu_A D_{Aj} + \mu_B D_{Bj} + \mu_C D_{Cj}$$

$$\frac{dI_{Kj}}{dt} = p_K \lambda_j S_j - \gamma_K I_{Kj} + \phi_K D_{Kj}, \quad K \in \{A,B,C\}$$

$$\frac{dD_{Kj}}{dt} = \gamma_K I_{Kj} - (\mu_K + \phi_K) D_{Kj}, \quad K \in \{A,B,C\}$$

### Discrete-Time Stochastic Extension

$$p_{\text{inf}}(j,t) = 1 - \exp\!\left(-\lambda_j(t)\,\Delta t\right)$$

$$\Delta I_{Kj}(t) \sim \mathrm{Binomial}\!\left(S_j(t),\; p_K \cdot p_{\text{inf}}(j,t)\right)$$

---

## Key Findings

**Neighbourhood infection pressure dominates.** Each additional new-breakdown neighbour was associated with a **37–47% increase in the odds** of belonging to a breakdown phenotype relative to control herds substantially larger than the effect of direct trade-mediated exposure from infected source herds.

**Three phenotypes are stable and reproducible.** Bootstrap resampling across 200 iterations produced mean Jaccard coefficients of **0.996, 0.998, and 0.982** for Clusters 1, 2, and 3 respectively, with zero dissolution events. Split-sample validation across five 70/30 partitions produced a mean ARI of **0.999**.

**Cluster 2 recurs at twice the rate of Cluster 1.** Cox proportional hazards modelling confirmed HR = 2.370 (95% CI: 2.177–2.580) for subsequent breakdown after full adjustment and this excess risk persists even after controlling for prior breakdown history, pointing to unmeasured persistent environmental or wildlife exposure.

**Cluster 3 takes twice as long to clear.** Mean breakdown duration of 333 days versus 149 days for Cluster 1, with a Cox HR of 0.114 for clearance rate — reflecting the sustained within-herd amplification that characterises the acute high-intensity phenotype.

---

## Technical Stack

```r
# Core simulation and modelling
library(SimInf)       # Stochastic metapopulation simulation (CTMC)
library(deSolve)      # Deterministic ODE integration (lsoda)

# Clustering
library(cluster)      # PAM clustering
library(ClusterR)     # K-means optimisation
library(mclust)       # Gaussian Mixture Models
library(NbClust)      # 26 cluster validity indices
library(fpc)          # Bootstrap stability (clusterboot)

# Network analysis
library(igraph)       # Cattle trade and neighbourhood network construction

# Statistical modelling
library(nnet)         # Multinomial logistic regression
library(survival)     # Kaplan-Meier and Cox proportional hazards

# Visualisation and data
library(ggplot2)
library(dplyr)
library(factoextra)
```

---

## Reproducibility

All analyses use a fixed random seed (`set.seed(123)`). The full pipeline is designed to be reproducible from raw surveillance data through to publication figures.

```r
# Entry point — run full pipeline
source("btb_main_pipeline_v4.R")
```

---

## Publications and Outputs

**Chukwudi V.**, Murphy T.B., Gormley C.I., Tratalos J., Madden J., McGrath G., Barrett D., Byrne A., Quinn D., McAloon C.G., Harvey N., Farrell D.
*Quantifying Heterogeneity in Bovine Tuberculosis Breakdown Dynamics Using Clustering and Finite Mixture Modelling.*
University College Dublin, School of Veterinary Medicine — **Manuscript under review**

---

## About

**Victory Chukwudi**
PhD Researcher — Stochastic Disease Modelling and AI Safety
University College Dublin · VistaMilk Research Ireland · Mastercard Foundation Scholar

[![LinkedIn](https://img.shields.io/badge/LinkedIn-Victory%20Chukwudi-blue?logo=linkedin)](https://www.linkedin.com/in/victory-chukwudi-531579241)
[![Portfolio](https://img.shields.io/badge/Portfolio-victory--chukwudi-lightgrey)](https://victory-chukwudi-portfolio.lovable.app/)

---
*"The goal is not to describe what we observe — it is to understand why it happens, and to build models rigorous enough to prove it."*
