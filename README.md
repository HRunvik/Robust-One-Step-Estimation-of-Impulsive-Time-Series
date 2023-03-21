# Robust-One-Step-Estimation-of-Impulsive-Time-Series
Code for parameter and input estimation according to the article with the same name.

## Authors
* Håkan Runvik (hakan.runvik 'at' it.uu.se),
* Alexander Medvedev (alexander.medvedev 'at' it.uu.se)

# Code
Examples demonstrating the estimation algorithms.

## startup.m
Setup paths

## firstorderdynamics.m
Monte carlo runs with estimation of b of first order plant and input impulses from synthetic data. Corresponds to the experiment in Section 4.1.1 in the article.

## secondorderdynamics.m
Monte carlo runs with estimation of b1, b2 of second order plant and input impulses from synthetic data. Corresponds to the experiment in Section 4.1.2 in the article.

## secondorderdynamics_basal.m
Monte carlo runs with estimation of b1, b2 of second order plant, basal level and input impulses from synthetic data. Corresponds to the experiment in Section 4.1.3 in the article.

## secondorderdynamics_robust.m
Monte carlo runs with estimation of b1, b2 of second order plant and input impulses from synthetic data corrupted with outliers. Corresponds to the experiment in Section 4.1.4 in the article.