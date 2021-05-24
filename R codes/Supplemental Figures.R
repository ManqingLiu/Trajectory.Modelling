library(dplyr)

library(ggplot2)
library(RColorBrewer)

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



### Supplemental Figure 1
load("pred1Y.Rdata")

library(fdapace)
library(EMCluster)

pred1Y <- pred1Y[with(pred1Y, order(PAT_ID, interval, decreasing = c(FALSE, FALSE))),]

#pred1Y_long <- pred1Y[c('PAT_ID', 'pred1Yiction_score', 'interval')]

pred1Y_fpca <- MakeFPCAInputs(pred1Y$PAT_ID, pred1Y$interval, pred1Y$prediction_score)

fpcaObjpred1Y <- FPCA(pred1Y_fpca$Ly, pred1Y_fpca$Lt, list(dataType='Sparse', plot = T, methodMuCovEst = 'smooth'))

# SelectK(fpcaObjpred1Y, criterion = 'FVE', FVEthreshold = 0.95) 

newClust <- FClust(pred1Y_fpca$Ly, pred1Y_fpca$Lt, k = 2, optnsFPCA =
                     list(methodMuCovEst = 'smooth', FVEthreshold = 0.95))

# cluster <- as.vector(newClust$cluster)
cluster <- newClust$cluster

pred1Y_id <- pred1Y[c('PAT_ID')]

pred1Y_id <- pred1Y_id %>% distinct(PAT_ID, .keep_all = TRUE)

pred1Y_id$cluster <- cluster

pred1Yf <- merge(pred1Y, pred1Y_id, by = "PAT_ID")

pred1Yf<- pred1Yf %>%
  group_by(cluster) %>%
  mutate(mean_pred1Y = mean(prediction_score))

pred1Yf$cluster_c<- with(pred1Yf, ifelse(mean_pred1Y > mean(prediction_score), "Cluster identified by 2nd FPC score",
                                     "Cluster identified by 1st FPC score"))


fig1A <- subset(pred1Yf, cluster_c == 'Cluster identified by 1st FPC score')
fig1B <- subset(pred1Yf, cluster_c == 'Cluster identified by 2nd FPC score')


fig1A$cut <- with(fig1A, ifelse(interval <= 60, 1,
                                ifelse(interval > 60 & interval <= 120, 2,
                                       ifelse(interval > 120 & interval <= 180, 3,
                                              ifelse(interval > 180 & interval <= 240, 4,
                                                     ifelse(interval > 240 & interval <= 300, 5,6))))))


countA <- fig1A %>% group_by(cut) %>% summarize(N = n_distinct(PAT_ID))

number <- countA[,2][1]

names(unlist(number[1:6]))

fig1a <- ggplot(fig1A, aes(x=interval, y=prediction_score, group = PAT_ID)) + 
  geom_line(alpha = 0.5, colour = 'slategray1') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 360, by = 60))+
  geom_smooth(aes(group=1),se=FALSE, colour="royalblue", size=1, method = 'loess')+
  theme_Publication()+
  labs(tag = "A")+annotate("text", x = c(30, 90,150,210,270,330), y =0.95, angle = 90,label =paste0("N=", as.character(countA$N)))

fig1B$cut <- with(fig1B, ifelse(interval <= 60, 1,
                                ifelse(interval > 60 & interval <= 120, 2,
                                       ifelse(interval > 120 & interval <= 180, 3,
                                              ifelse(interval > 180 & interval <= 240, 4,
                                                     ifelse(interval > 240 & interval <= 300, 5,6))))))


countB <- fig1B %>% group_by(cut) %>% summarize(N = n_distinct(PAT_ID))

paste0("N=", as.character(countB$N))

fig1b <- ggplot(fig1B, aes(x=interval, y=prediction_score, group = PAT_ID)) + 
  geom_line(alpha = 0.5, colour = 'lightpink') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(0, 1, by = 0.2))+
  scale_x_continuous(name = 'Time to death (days)', breaks = seq(0, 360, by = 60))+
  geom_smooth(aes(group=1),se=FALSE, colour="firebrick1", size=1, method = 'loess')+
  theme_Publication()+
  annotate("text", x =  c(30, 90,150,210,270,330), y =0.93, angle = 90,label = paste0("N=", as.character(countB$N)))+
  labs(tag = "B")

library(gridExtra)
grid.arrange(fig1a, fig1b,nrow=1)

#### Supplemental Figure 2
load(file = "bl_eol.Rdata")

TAB <- bl_eol

TAB <- TAB[c("InpatientDeath", "IcuLast30Days", "enroll", "PennHospiceEligible",
             "ChemoLast14Days", "cluster_c")]

vars <- c("InpatientDeath", "IcuLast30Days", "enroll", "PennHospiceEligible", "ChemoLast14Days")

TAB[vars] <- lapply(TAB[vars], as.numeric)

# substract <- function(x) {
#   x-1
# }
# 
# TAB[vars] <- apply(TAB[vars],2, substract)


sum_tab <- TAB %>% 
  group_by(cluster_c) %>%
  summarize(SumInpdeath = sum(InpatientDeath),
            SumICU = sum(IcuLast30Days),
            SumEligible = sum(PennHospiceEligible),
            SumEnroll = sum(enroll, na.rm = T),
            SumChemo = sum(ChemoLast14Days))

N = c(873, 2027)

sum_tab <- cbind(sum_tab,N)

library(reshape)
sum_tab2 <- melt(sum_tab, id=c("cluster_c", "N"))

z = qnorm(0.975)

tab_p1 <- subset(sum_tab2, !(variable %in% c('SumEnroll', 'SumEligible')) )

sum_tab3 <- tab_p1 %>%
  mutate(Rate = (value/N)*100,
         Ymin = ((value+(z^2)/2)/(N+z^2)
                 -(z/(N+z^2))*sqrt(value*(N-value)/N+z^2/4))*100,
         Ymax = ((value+(z^2)/2)/(N+z^2)
                 +(z/(N+z^2))*sqrt(value*(N-value)/N+z^2/4))*100)

tab_p2 <- subset(sum_tab2, (variable %in% c('SumEnroll')) )
tab_p2$N <- c(526,1393)

sum_tab4 <- tab_p2 %>%
  mutate(Rate = (value/N)*100,
         Ymin = ((value+(z^2)/2)/(N+z^2)
                 -(z/(N+z^2))*sqrt(value*(N-value)/N+z^2/4))*100,
         Ymax = ((value+(z^2)/2)/(N+z^2)
                 +(z/(N+z^2))*sqrt(value*(N-value)/N+z^2/4))*100)

sum_tabf <- rbind(sum_tab3, sum_tab4)

sum_tabf$cluster_label <- with(sum_tabf, ifelse(cluster_c == 'Cluster identified by 1st FPC score', 'Unpredictable trajectories',
                                                'Predictable trajectories'))

sum_tabf$cluster_label <- factor(sum_tabf$cluster_label,
                levels = c("Unpredictable trajectories", "Predictable trajectories"))

ggplot(sum_tabf, aes(x=variable, y=Rate, fill=cluster_label)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Ymin, ymax=Ymax), width=.2,
                position=position_dodge(.9))+
  labs(tag="", x="", y = "% of Decedents")+
  theme_Publication()+
  theme(legend.title=element_blank(),axis.text.x = element_text(angle=0, size = 12))+
  guides(guide_legend(" ")) +
  scale_y_continuous(breaks = seq(0,63, by = 10), limits = c(0,63))+
  scale_fill_manual(values=c('royalblue','firebrick1'))+
  scale_x_discrete(labels=c("SumInpdeath" = "Inpatient Death", "SumICU" = "Admitted to ICU \n in last 30 days",
                            "SumEnroll" = "Enrolled to hospice", "SumChemo" = "Chemotherapy \n in last 14 days"))

#### Suplemental Table 3

dat <- read.csv("cohort_traj_2018Jan_2020May.csv", header  = T)

vars <- c("PAT_ID", "prediction_score", "PAT_AGE", "female", "race", "marital_status", "insurance", "ecog_strict_c", "elixhauser_cat",
          "ECOG_strict", "elix_count", "SICP_SITE_TEAM", "treat_loc", "label", "DEATH_DATE", 
          "APPT_TIME", "STAGE_GROUP",'n_is_METS',"n_is_Lymp",'NPI','EDIT_DATETIME','CODE_TIME_DATE')

dat <- dat[vars]

dat <- dat %>% 
  mutate(
    appt_date = as.Date(APPT_TIME, format = "%Y-%m-%d"), 
    death_date = as.Date(DEATH_DATE, format = "%Y-%m-%d"),
    stage_date = as.Date(EDIT_DATETIME, format = "%Y-%m-%d"),
    ecog_date = as.Date(CODE_TIME_DATE, format = "%Y-%m-%d")
  )


control <- anti_join(dat, predf, by = 'PAT_ID')

control <- control[!duplicated(control[c('PAT_ID','appt_date')]),]


control <- control[with(control, order(PAT_ID, appt_date, decreasing = c(FALSE, TRUE))),]

control <- control %>% group_by(PAT_ID) %>%
  mutate(days_to_death = first(appt_date)-appt_date, days_to_death)

control$days_to_death <- as.numeric(control$days_to_death)

control <- control[!(control$days_to_death >= 180), ]


control$status <- 'control'

library(fdapace)
library(EMCluster)

control_long <- control[c('PAT_ID', 'prediction_score', 'days_to_death')]

control_fpca <- MakeFPCAInputs(control_long$PAT_ID, control_long$days_to_death, control_long$prediction_score)

fpcaObjcontrol <- FPCA(control_fpca$Ly, control_fpca$Lt, list(dataType='Sparse', plot = TRUE, methodMuCovEst = 'smooth'))

SelectK(fpcaObjcontrol, criterion = 'FVE', FVEthreshold = 0.95) ## K = 2

newClust <- FClust(control_fpca$Ly, control_fpca$Lt, k = 2, optnsFPCA =
                     list(methodMuCovEst = 'smooth', FVEthreshold = 0.95))

# cluster <- as.vector(newClust$cluster)
cluster <- newClust$cluster

pred_id <- control[c('PAT_ID')]

pred_id <- pred_id %>% distinct(PAT_ID, .keep_all = TRUE)

pred_id$cluster <- cluster

controlf <- merge(control, pred_id, by = "PAT_ID")

controlf<- controlf %>%
  group_by(cluster) %>%
  mutate(mean_control = mean(prediction_score))

controlf$cluster_c<- with(controlf, ifelse(mean_control > mean(prediction_score), "Cluster identified by 2nd FPC score",
                                     "Cluster identified by 1st FPC score"))

save(controlf, file = 'controlf.Rdata')


fig <- controlf

#fig$cluster <- with(fig, ifelse(cluster_pred = '1st FPC score', 'Cluster identified by 1st FPC score', 'Cluster identified by 2nd FPC score'))

fig1A <- subset(fig, cluster_c == 'Cluster identified by 1st FPC score')
fig1B <- subset(fig, cluster_c == 'Cluster identified by 2nd FPC score')


fig1A$cut <- with(fig1A, ifelse(days_to_death <= 30, 1,
                                ifelse(days_to_death > 30 & days_to_death <= 60, 2,
                                       ifelse(days_to_death > 60 & days_to_death <= 90, 3,
                                              ifelse(days_to_death > 90 & days_to_death <= 120, 4,
                                                     ifelse(days_to_death > 120 & days_to_death <= 150, 5,6))))))


countA <- fig1A %>% group_by(cut) %>% summarize(N = n_distinct(PAT_ID))

number <- countA[,2][1]


fig1a <- ggplot(fig1A, aes(x=days_to_death, y=prediction_score, group = PAT_ID)) + 
  geom_line(alpha = 0.5, colour = 'slategray1') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  scale_x_continuous(name = 'Time to last appointment (days)', breaks = seq(0, 180, by = 30))+
  geom_smooth(aes(group=1),se=FALSE, colour="royalblue", size=1, method = 'loess')+
  theme_Publication()+
  labs(tag = "A")+annotate("text", x = c(15, 45,75,105,135,165), y =0.95, angle = 90,label =paste0("N=", as.character(countA$N)))

fig1B$cut <- with(fig1B, ifelse(days_to_death <= 30, 1,
                                ifelse(days_to_death > 30 & days_to_death <= 60, 2,
                                       ifelse(days_to_death > 60 & days_to_death <= 90, 3,
                                              ifelse(days_to_death > 90 & days_to_death <= 120, 4,
                                                     ifelse(days_to_death > 120 & days_to_death <= 150, 5,6))))))


countB <- fig1B %>% group_by(cut) %>% summarize(N = n_distinct(PAT_ID))

paste0("N=", as.character(countB$N))

fig1b <- ggplot(fig1B, aes(x=days_to_death, y=prediction_score, group = PAT_ID)) + 
  geom_line(alpha = 0.5, colour = 'lightpink') + guides(colour=FALSE) + xlab("Days to death") +
  scale_y_continuous(name = 'Mortality Risk', breaks = seq(0, 1, by = 0.2))+
  scale_x_continuous(name = 'Time to last appointment (days)', breaks = seq(0, 180, by = 30))+
  geom_smooth(aes(group=1),se=FALSE, colour="firebrick1", size=1, method = 'loess')+
  theme_Publication()+
  annotate("text", x = c(15, 45,75,105,135,165), y =0.93, angle = 90,label = paste0("N=", as.character(countB$N)))+
  labs(tag = "B")

library(gridExtra)
grid.arrange(fig1a, fig1b,nrow=1)
