
## Perform FPCA analysis 
library(fdapace)
library(EMCluster)

## 1. Data with longitudinal mortality risks (pred). Interval is the days from each appointment to date of death. 
pred <- pred[with(pred, order(PAT_ID, interval, decreasing = c(FALSE, FALSE))),]

## 2. Turn vector inputs to the list so they can be used in FPCA
pred_fpca <- MakeFPCAInputs(pred$PAT_ID, pred$interval, pred$prediction_score)

## 3. Run FPCA. Ly is the observed moratlity risks for each patient, Lt is days from each appointment to death corresponding to each Ly
fpcaObjpred <- FPCA(pred_fpca$Ly, pred_fpca$Lt, list(dataType='Sparse', plot = T, methodMuCovEst = 'smooth'))

## 4. Selects number of functional principal components for given FPCA output using a FVE threshold of 95%
SelectK(fpcaObjpred, criterion = 'FVE', FVEthreshold = 0.95) 

## 5. Cluster functional data using the functional principal component (FPC) scores obtained
## from the data FPC analysis using EMCluster (Chen and Maitra, 2015)
newClust <- FClust(pred_fpca$Ly, pred_fpca$Lt, k = 2, optnsFPCA =
                     list(methodMuCovEst = 'smooth', FVEthreshold = 0.95))


## 6. Assign cluster indicator to patient-encounter data 
library(dplyr)

cluster <- newClust$cluster

pred_id <- pred[c('PAT_ID')]

pred_id <- pred_id %>% distinct(PAT_ID, .keep_all = TRUE)

pred_id$cluster <- cluster

predf <- merge(pred, pred_id, by = "PAT_ID")

predf<- predf %>%
  group_by(cluster) %>%
  mutate(mean_pred = mean(prediction_score))

predf$cluster_c<- with(predf, ifelse(mean_pred > mean(prediction_score), "2nd Cluster",
                                     "1st Cluster"))

save(predf, file = 'predf.Rdata')




