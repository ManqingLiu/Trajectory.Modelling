
## 1. load packages and design theme of plot 
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(fdapace)
library(EMCluster)

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

## 2. load data
load('predf.Rdata')
load('fpcaOBJ.Rdata')


## 3. Draw figures
### 3.1 Figure 1. Trajectories of mortality risk
# Figure 1A. Spaghetti plot  with smoothed line of FPC1
# Figure 1B. Spaghetti plot  with smoothed line of FPC2

fig <- predf

fig1A <- subset(fig, cluster_c == '1st Cluster')
fig1B <- subset(fig, cluster_c == '2nd Cluster')


fig1A$cut <- with(fig1A, ifelse(interval <= 30, 1,
                                ifelse(interval > 30 & interval <= 60, 2,
                                       ifelse(interval > 60 & interval <= 90, 3,
                                              ifelse(interval > 90 & interval <= 120, 4,
                                                     ifelse(interval > 120 & interval <= 150, 5,6))))))


countA <- fig1A %>% group_by(cut) %>% summarize(N = n_distinct(PAT_ID))

number <- countA[,2][1]

names(unlist(number[1:6]))

fig1a <- ggplot(fig1A, aes(x=interval, y=prediction_score, group = PAT_ID)) + 
  geom_line(alpha = 0.5, colour = 'slategray1') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 180, by = 30))+
  geom_smooth(aes(group=1),se=FALSE, colour="royalblue", size=1, method = 'loess')+
  theme_Publication()+
  labs(tag = "A")+annotate("text", x = c(15, 45,75,105,135,165), y =0.95, angle = 90,label =paste0("N=", as.character(countA$N)))

fig1B$cut <- with(fig1B, ifelse(interval <= 30, 1,
                                ifelse(interval > 30 & interval <= 60, 2,
                                       ifelse(interval > 60 & interval <= 90, 3,
                                              ifelse(interval > 90 & interval <= 120, 4,
                                                     ifelse(interval > 120 & interval <= 150, 5,6))))))


countB <- fig1B %>% group_by(cut) %>% summarize(N = n_distinct(PAT_ID))

paste0("N=", as.character(countB$N))

fig1b <- ggplot(fig1B, aes(x=interval, y=prediction_score, group = PAT_ID)) + 
  geom_line(alpha = 0.5, colour = 'lightpink') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(0, 1, by = 0.2))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 180, by = 30))+
  geom_smooth(aes(group=1),se=FALSE, colour="firebrick1", size=1, method = 'loess')+
  theme_Publication()+
  annotate("text", x = c(15, 45,75,105,135,165), y =0.93, angle = 90,label = paste0("N=", as.character(countB$N)))+
  labs(tag = "B")

library(gridExtra)
grid.arrange(fig1a, fig1b,nrow=1)

# Figure 1C. Eigenfunctions of FPC1
# Figure 1D. Eigenfunction of FPC2

fittedY <- fitted(fpcaOBJ, ciOptns = list(alpha=0.05), cvgMethod = 'band')
workGrid <- fpcaOBJ$workGrid

fpc <- as.data.frame(fpcaOBJ$xiEst)

first2eigen <- as.data.frame(fpcaOBJ$phi[, 1:2])
first2eigen$x.axis <- seq(0,180, by = 3.6)

eigen1 <- first2eigen[c('V1', 'x.axis')]
eigen2 <- first2eigen[c('V2', 'x.axis')]

fig1c <- ggplot(eigen1, aes(x=x.axis, y=V1)) + 
  geom_line(colour = 'royalblue') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Eigenfunction', breaks = seq(-0.1, 0.12, by = 0.04), limits = c(-0.1,0.12))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 180, by = 30))+
  theme_Publication()+
  labs(tag = "C")

fig1d <- ggplot(eigen2, aes(x=x.axis, y=V2)) + 
  geom_line(colour = 'firebrick1') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Eigenfunction',  breaks = seq(-0.1, 0.12, by = 0.04), limits = c(-0.1,0.12))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 180, by = 30))+
  theme_Publication()+
  labs(tag = "D")

grid.arrange(fig1c, fig1d,nrow=1)

# Figure 1E. Representative individual patient with maximal FPC score for FPC1 among 1st cluster with >= 10 encounters
# Figure 1F. Representative individual patient with maximal FPC score for FPC2 among 2nd cluster with >= 10 encounters

cvgUpper <- fittedY$cvgUpper
cvgLower <- fittedY$cvgLower

fig1e <- as.data.frame(cbind(workGrid, fittedY$fitted[1611,], cvgUpper[1611,], cvgLower[1611,],
                             fpcaOBJ$inputData$Lt[[1611]],fpcaOBJ$inputData$Ly[[1611]]))

fig1f <- as.data.frame(cbind(workGrid, fittedY$fitted[988,], cvgUpper[988,], cvgLower[988,],
                             fpcaOBJ$inputData$Lt[[988]],fpcaOBJ$inputData$Ly[[988]]))


fig1E <- ggplot() + 
  geom_line(data = fig1e, aes(x = workGrid, y = V2), colour = 'royalblue') + 
  geom_line(data = fig1e, aes(x = workGrid, y = V3), colour = 'dodgerblue', linetype = "dashed") +
  geom_line(data = fig1e, aes(x = workGrid, y = V4), colour = 'dodgerblue', linetype = "dashed") +
  geom_point(data = fig1e, aes(x = V5, y = V6), colour = 'blue', shape = 1)+
  guides(colour=FALSE) +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(-0.1,0.9, by = 0.1), limits = c(-0.1,0.9))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 180, by = 30))+
  theme_Publication()+
  labs(tag = "E")

fig1F <- ggplot() + 
  geom_line(data = fig1f, aes(x = workGrid, y = V2), colour = 'firebrick1') + 
  geom_line(data = fig1f, aes(x = workGrid, y = V3), colour = 'lightpink', linetype = "dashed") +
  geom_line(data = fig1f, aes(x = workGrid, y = V4), colour = 'lightpink', linetype = "dashed") +
  geom_point(data = fig1f, aes(x = V5, y = V6), colour = 'red', shape = 1)+
  guides(colour=FALSE) +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(-0.1,1, by = 0.1), limits = c(-0.1,0.9))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 180, by = 30))+
  theme_Publication()+
  labs(tag = "F")

grid.arrange(fig1E, fig1F,nrow=1)


# Figure 1G. Validation plot of FPC1 vs. FPC2
set.seed(152)
train <- pred %>% 
  group_by(PAT_ID) %>% 
  summarise() %>% 
  sample_frac(0.75) %>%
  left_join(pred)

test <- anti_join(pred, train)
train <- train[with(train, order(PAT_ID, interval, decreasing = c(FALSE, FALSE))),]
train_fpca <- MakeFPCAInputs(train$PAT_ID, train$interval, train$prediction_score)
tr <- FPCA(train_fpca$Ly, train_fpca$Lt, list(dataType='Sparse', plot = FALSE, methodMuCovEst = 'smooth'))
tr_fpc <- tr$xiEst[,1:2]
test <- test[with(test, order(PAT_ID, interval, decreasing = c(FALSE, FALSE))),]
test_fpca <- MakeFPCAInputs(test$PAT_ID, test$interval, test$prediction_score)
test_fpc <- predict(tr, test_fpca$Ly, test_fpca$Lt, K = 7, methodMuCovEst = 'smooth')
tr_fpc <- as.data.frame(tr_fpc)
test_fpc <- as.data.frame(test_fpc)
tr_fpc$type = "Training set (N = 2460)"
tr_fpc$lab = 1
test_fpc$type = "Validation set (N = 820)"
test_fpc$lab = 16

ggplot() + 
  geom_point(data = tr_fpc, aes(x = V1, y = V2, shape = type))+
  geom_point(data = test_fpc, aes(x = V1, y = V2, shape = type))+
  scale_shape_manual(values = c(1,16))+
  guides(shape=guide_legend(" ", nrow = 2, byrow = T)) +
  scale_y_continuous(name = paste0("FPC 2", " (", round(tr$cumFVE[2]-tr$cumFVE[1]), "%)"), breaks = seq(-4.3, 4, by = 1), limits = c(-4.3,4))+
  scale_x_continuous(name = paste0("FPC 1"," (", round(tr$cumFVE[1]), "%)"), breaks = seq(-6.4,7, by = 2),limits = c(-6.4, 7))+
  theme_Publication()+
  theme(legend.position = c(0.15,0.95))+
  labs(tag = "G")


# 3.2 Figure 2. Association between trajectories and EOL outcomes
# forest plot of adjusted odds ratios for each EOL outcome.

## load baseline patient characteristics data 
load(file = "bl_eol.Rdata")

## 3.2.1 Run multiple imputation via MICE for each of the EOL outcome 
cols <- c('Cancer', 'race','cluster','treat_loc','female','marital_status','elixhauser_cat','stage',
          'year_death','Insurance','age_cat','InpatientDeath','enroll','IcuLast30Days','ChemoLast14Days')
bl_eol[cols] <- lapply(bl_eol[cols], factor)
vars <- c("cluster",'age_cat','female','race','marital_status','Insurance','elixhauser_cat',
          'stage','year_death','treat_loc','Cancer','InpatientDeath','enroll','IcuLast30Days','ChemoLast14Days')
keep <- bl_eol[vars]
completeD2 = list()
for (i in 1:50) {
  completeD2[[i]] <- cbind(completeD[[i]], keep)
}
for (i in 1:50) {
  completeD2[[i]]$ECOG_strict_c <- with(completeD2[[i]], ifelse(ECOG_strict == 0 | ECOG_strict == 1, "0-1", "2+"))
  completeD2[[i]]$ECOG_strict_c <- as.factor(completeD2[[i]]$ECOG_strict_c)
}

### Inpatient death 
modelFit = list()
for (i in 1:50) {
  modelFit[[i]] <- glm(InpatientDeath ~ clus + age_cat + count + female + race + marital_status + Insurance +
                         elixhauser_cat + ECOG_strict_c + stage + year_death+treat_loc+Cancer, data = completeD2[[i]], family = "binomial")
}

sum <- summary(pool(modelFit))

sum2 <- sum[2, ]
est1 <- sum2[,1]
low_ci1 <- sum2[,1] - qt(0.975, df = sum2[,4]) * sum2[, 2]
high_ci1 <- sum2[,1] + qt(0.975, df = sum2[,4]) * sum2[, 2]
ci1 <- paste0( "[", round(exp(low_ci1),2), "-", round(exp(high_ci1),2), "]")
row1 <- cbind('OR' = exp(est1), 'll' = low_ci1, 'hl' =high_ci1, 'P Value' = round(sum2$p.value,3))

### ICU in last 30 days
modelFit = list()
for (i in 1:50) {
  modelFit[[i]] <- glm(IcuLast30Days ~ clus + age_cat + count + female + race + marital_status + Insurance +
                         elixhauser_cat + ECOG_strict_c + stage + year_death+treat_loc+Cancer, data = completeD2[[i]], family = "binomial")
}


sum <- summary(pool(modelFit))
sum2 <- sum[2, ]
est1 <- sum2[,1]
low_ci1 <- sum2[,1] - qt(0.975, df = sum2[,4]) * sum2[, 2]
high_ci1 <- sum2[,1] + qt(0.975, df = sum2[,4]) * sum2[, 2]
ci1 <- paste0( "[", round(exp(low_ci1),2), "-", round(exp(high_ci1),2), "]")
row2 <- cbind('OR' = exp(est1), 'll' = low_ci1, 'hl' =high_ci1, 'P Value' = round(sum2$p.value,3))

### Enrolled to hospice
modelFit = list()
for (i in 1:50) {
  modelFit[[i]] <- glm(enroll ~ clus + age_cat + count + female + race + marital_status + Insurance +
                         elixhauser_cat + ECOG_strict_c + stage + year_death+treat_loc+Cancer, data = completeD2[[i]], family = "binomial")
}
sum <- summary(pool(modelFit))
sum2 <- sum[2, ]
est1 <- sum2[,1]
low_ci1 <- sum2[,1] - qt(0.975, df = sum2[,4]) * sum2[, 2]
high_ci1 <- sum2[,1] + qt(0.975, df = sum2[,4]) * sum2[, 2]
ci1 <- paste0( "[", round(exp(low_ci1),2), "-", round(exp(high_ci1),2), "]")
row3 <- cbind('OR' = exp(est1), 'll' = low_ci1, 'hl' =high_ci1, 'P Value' = round(sum2$p.value,3))

### Chemotherapy in last 14 days
modelFit = list()
for (i in 1:50) {
  modelFit[[i]] <- glm(ChemoLast14Days ~ clus + age_cat + count + female + race + marital_status + Insurance +
                         elixhauser_cat + ECOG_strict_c + stage + year_death+treat_loc+Cancer, data = completeD2[[i]], family = "binomial")
}
sum <- summary(pool(modelFit))
sum2 <- sum[2, ]
est1 <- sum2[,1]
low_ci1 <- sum2[,1] - qt(0.975, df = sum2[,4]) * sum2[, 2]
high_ci1 <- sum2[,1] + qt(0.975, df = sum2[,4]) * sum2[, 2]
ci1 <- paste0( "[", round(exp(low_ci1),2), "-", round(exp(high_ci1),2), "]")
row4 <- cbind('OR' = exp(est1), 'll' = low_ci1, 'hl' =high_ci1, 'P Value' = round(sum2$p.value,3))

### 3.2.2 Draw forest plot using the odds ratios from MI 
library(knitr)
library(kableExtra)

tab <- rbind(row1, row2, row3, row4)

rownames(tab) <- c("Inpatient Death", "Admitted to ICU in last 30 days", 
                   "Enrolled to hospice or not","Chemotherapy in last 14 days")
tab[,2:3] <- exp(tab[,2:3])
tab <- as.data.frame(tab)
tab <- tab %>% mutate_at(vars(OR, ll, hl), funs(sprintf('%.2f',.)))

library(forestplot)
rmeta <- 
  structure(list(
    mean  = c(NA, tab[1,1], tab[2,1], tab[3,1], tab[4,1]), 
    lower = c(NA, tab[1,2], tab[2,2], tab[3,2], tab[4,2]),
    upper = c(NA, tab[1,3], tab[2,3], tab[3,3], tab[4,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L), 
    class = "data.frame")


tabletext<-cbind(
  c("EOL Metrics", "Inpatient Death", "Admitted to ICU in last 30 days", 
    "Enrolled to hospice","Chemotherapy in last 14 days"),
  # c("N", "251", "681", "215", "196", "140", "176", "173", "487"),
  c("OR", toString(tab[1,1]), toString( tab[2,1]), toString(tab[3,1]), toString(tab[4,1])),
  c("95% CI",paste0( tab[1,2], "-", tab[1,3]),  paste0( tab[2,2], "-",  tab[2,3]), 
    paste0( tab[3,2], "-", tab[3,3]),  paste0( tab[4,2], "-",tab[4,3]))
)

rmeta$mean <- as.numeric(rmeta$mean)
rmeta$lower <- as.numeric(rmeta$lower)
rmeta$upper <- as.numeric(rmeta$upper)

forestplot(tabletext, 
           rmeta,
           xlog=T, 
           #graphwidth = unit(6,"cm"),
           lineheight = unit(1,"cm"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10, fontfamily = "serif"), ticks=gpar(cex=0.8)),
           is.summary=c(TRUE,rep(FALSE,4)),
           col=fpColors(box="black",line="black"),
           xlim = c(0.5, 2.5),
           xticks = seq(from = .5, to = 2.5, by = 0.5))
