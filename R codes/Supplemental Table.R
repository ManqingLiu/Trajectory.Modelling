
library(readxl)
library(dplyr)

load('predf.Rdata')
load('bl_f.Rdata')

tab <- bl_f

tab$clus <- with(tab, ifelse(cluster_c == 'Cluster identified by 2nd FPC score', 1, 0))

vars <- c("PAT_AGE", "insurance", "count", "elix_count", "ECOG_strict", "treat_loc")
impute <- tab[vars]

library(mice)
tempData <- mice(impute,m=50,maxit=5,meth='pmm',seed=500)

completeD <- complete(tempData, "all")

cols <- c('Cancer', 'race','clus','treat_loc','female','marital_status','elixhauser_cat','stage2',
          'year_death', 'Insurance','age_cat')

#rdata = factor(data,labels=c("I","II","III"))

tab[cols] <- lapply(tab[cols], factor)

vars <- c("clus",'age_cat','female','race','marital_status','Insurance','elixhauser_cat',
          'stage2','year_death','treat_loc','Cancer')

keep <- tab[vars]

completeD2 = list()
for (i in 1:50) {
  completeD2[[i]] <- cbind(completeD[[i]], keep)
}

for (i in 1:50) {
  completeD2[[i]]$ECOG_strict_c <- with(completeD2[[i]], ifelse(ECOG_strict == 0 | ECOG_strict == 1, "0-1", "2+"))
  completeD2[[i]]$ECOG_strict_c <- as.factor(completeD2[[i]]$ECOG_strict_c)
}


modelFit = list()
for (i in 1:50) {
  modelFit[[i]] <- glm(clus ~ age_cat + count + female + race + marital_status + Insurance +
                         elixhauser_cat + ECOG_strict_c + stage2 + year_death+treat_loc+Cancer, data = completeD2[[i]], family = "binomial")
}


sum <- summary(pool(modelFit))

est1 <- sum[,2]
low_ci1 <- sum[,2] - qt(0.975, df = sum[,5]) * sum[, 3]
high_ci1 <- sum[,2] + qt(0.975, df = sum[,5]) * sum[, 3]


ci1 <- paste0( "[", sprintf('%.2f',exp(low_ci1)), "-", sprintf('%.2f',exp(high_ci1)), "]")
Names <- data.frame(sum[,1])

row1 <- cbind(Names = Names, OR = round(exp(est1),2), '95% CI' = ci1, 'PValue' = sprintf('%.3f',sum$p.value))

row1 <- as.data.frame(row1)