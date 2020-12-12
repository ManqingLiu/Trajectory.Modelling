# Trajectory Modelling
## Objectives
1. Identify hidden clusters of trajectories of longitudinal mortality risks from gradient boosting machine (GBM) algorithm. This algorithm was designed to predict 180-day mortality among outpatients with cancer (Parikh et al. 2019). 
2. Evaluate the association between identified trajectories and End-of-Life care patients received. 

## FPCA analysis 
R codes/FPCA.R 

Key packages used: 
1. FPCA: _fdapace_ R package (Carroll et al. 2020).  
Github source of this package: https://github.com/functionaldata/tPACE
2. EM Algorithm for Model-Based Clustering: _EMCluster_ R package (Chen & Maitra. 2020). 

## Figures
R codes/Figures.R

![Figure1AB](https://github.com/ManqingLiu/Trajectory.Modelling/blob/main/Figures/Figure1AB.png)
![Figure1CD](https://github.com/ManqingLiu/Trajectory.Modelling/blob/main/Figures/Figure1CD.png)
![Figure1EF](https://github.com/ManqingLiu/Trajectory.Modelling/blob/main/Figures/Figure1EF.png)
![Figure1G](https://github.com/ManqingLiu/Trajectory.Modelling/blob/main/Figures/Figure1G.jpg)
![Figure2](https://github.com/ManqingLiu/Trajectory.Modelling/blob/main/Figures/Figure2.jpg)

Key package used:  
1. FPCA: _fdapace_ R package (Carroll et al. 2020).  
Github source of this package: https://github.com/functionaldata/tPACE
2. EM Algorithm for Model-Based Clustering: _EMCluster_ R package (Chen & Maitra. 2020). 
3. Multiple Imputation: _mice_ R package (Buuren & Groothuis-Oudshoorn 2011).

## References
Parikh, R. B. et al. Machine Learning Approaches to Predict 6-Month Mortality Among Patients With Cancer. JAMA Netw Open 2, e1915997 (2019). doi: 10.1001/jamanetworkopen.2019.15997  
Carroll C. et al. fdapace: Functional Data Analysis and Empirical Dynamics. R package. (2020). version 0.5.3. URL: https://CRAN.R-project.org/package=fdapace  
Chen, W.C. & Maitra, R. EMCluster: EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian Distribution. R Package. (2020). URL: https://CRAN.R-project.org/package=EMCluster  
Buuren, S. & Groothuis-Oudshoorn, K.  mice: Multivariate Imputation by Chained Equations in R. Journal of
  Statistical Software, 45(3), 1-67. (2011). URL: https://www.jstatsoft.org/v45/i03/.
