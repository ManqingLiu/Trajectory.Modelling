# Trajectory Modelling
## Objectives
1. Identify hidden clusters of trajectories of longitudinal mortality risks from gradient boosting machine (GBM) algorithm. This algorithm was designed to predict 180-day mortality among outpatients with cancer (Parikh et al. 2019). 
2. Evaluate the association between identified trajectories and End-of-Life care patients received. 

## FPCA analysis 
R codes/FPCA.R 
```
## Perform FPCA analysis 
library(fdapace)
library(EMCluster)

## 1. Data with longitudinal mortality risks (pred). Interval is the days from each appointment to date of death. 
pred <- pred[with(pred, order(PAT_ID, interval, decreasing = c(FALSE, FALSE))),]

## 2. Turn vector inputs to the list so they can be used in FPCA
pred_fpca <- MakeFPCAInputs(pred$PAT_ID, pred$interval, pred$prediction_score)

## 3. Run FPCA. Ly is the observed moratlity risks for each patient, Lt is days from each appointment to death corresponding to each Ly
fpcaObjpred <- FPCA(pred_fpca$Ly, pred_fpca$Lt, list(dataType='Sparse', plot = T, methodMuCovEst = 'smooth'))

save(fpcaOBJpred, file="fpcaOBJ.RData")

## 4. Selects number of functional principal components for given FPCA output using a FVE threshold of 95%
SelectK(fpcaObjpred, criterion = 'FVE', FVEthreshold = 0.95) 

## 5. Cluster functional data using the functional principal component (FPC) scores obtained
## from the data using EMCluster (Chen and Maitra, 2015)
newClust <- FClust(pred_fpca$Ly, pred_fpca$Lt, k = 2, optnsFPCA =
                     list(methodMuCovEst = 'smooth', FVEthreshold = 0.95))

```

Package used: 
1. FPCA: _fdapace_ R package (Carroll et al. 2020).  
Github source of this package: https://github.com/functionaldata/tPACE
2. EM Algorithm for Model-Based Clustering: _EMCluster_ R package (Chen & Maitra. 2020). 

## Figures
R codes/Figures.R

Package used:
Multiple Imputation: _mice_ R package (Buuren & Groothuis-Oudshoorn 2011).

## References
Parikh, R. B. et al. Machine Learning Approaches to Predict 6-Month Mortality Among Patients With Cancer. JAMA Netw Open 2, e1915997 (2019). doi: 10.1001/jamanetworkopen.2019.15997  
Carroll C. et al. fdapace: Functional Data Analysis and Empirical Dynamics. R package. (2020). version 0.5.3. URL: https://CRAN.R-project.org/package=fdapace  
Chen, W.C. & Maitra, R. EMCluster: EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian Distribution. R Package. (2020). URL: https://CRAN.R-project.org/package=EMCluster  
Buuren, S. & Groothuis-Oudshoorn, K.  mice: Multivariate Imputation by Chained Equations in R. Journal of
  Statistical Software, 45(3), 1-67. (2011). URL https://www.jstatsoft.org/v45/i03/.
