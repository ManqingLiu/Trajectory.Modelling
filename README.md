# Trajectory Modelling
## Objectives
1. Identify hidden clusters of trajectories of longitudinal mortality risks from gradient boosting machine (GBM) algorithm. This algorithm was designed to predict 180-day mortality among outpatients with cancer (Parikh et al. 2019). 
2. Evaluate the association between identified trajectories and End-of-Life care patients received. 

## FPCA analysis 
R codes/FPCA.R 

Package used: 
1. FPCA: _fdapace_ R package (Carroll et al. 2020).  
Github source of this package: https://github.com/functionaldata/tPACE
2. EM Algorithm for Model-Based Clustering: _EMCluster_ R package (Chen & Maitra. 2020). 

## Figures


## References
Parikh, R. B. et al. Machine Learning Approaches to Predict 6-Month Mortality Among Patients With Cancer. JAMA Netw Open 2, e1915997 (2019). doi: 10.1001/jamanetworkopen.2019.15997  
Carroll C. et al. fdapace: Functional Data Analysis and Empirical Dynamics. R package. (2020). version 0.5.3. URL: https://CRAN.R-project.org/package=fdapace  
Chen, W.C. & Maitra, R. EMCluster: EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian Distribution. R Package. (2020). URL: https://CRAN.R-project.org/package=EMCluster
