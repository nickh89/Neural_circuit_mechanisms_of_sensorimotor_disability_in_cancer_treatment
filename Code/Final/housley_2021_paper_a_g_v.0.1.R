# Version of Jan 6, 2021.
# Stephen N. Housley 
# housley.nick@gmail.com
# 

# This work is licensed under the licenses 
# Paper: Creative Commons Attribution 3.0 Unported License 
# Code: GPL-3 License 
# Depends: R (>= 3.5.0)
# Version: 0.1
# Description: code to run analytics and graphic functions associated with:
#         ZZZZZZZZZZZZZ
#
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 
#

## version Hx
# v0.1- original
# Commit comments



## TEMPLATE ##

###################################################################
########################### Figure ZZZZ ###########################
###################################################################

########################### load dependencies ###########################
########################### Custom Functions ###########################
########################### Load Data ###########################
########################### Data Wrangling ###########################
########################### quick visualization ###########################
########################### analyses/modeling ###########################
########################### saving data ###########################
########################### saving figures ###########################
########################### Clean up ###########################


## TEMPLATE ##



########################### prelims ########################### 

invisible(rm(list = ls()))
invisible(gc())

########################### timer start ########################
ptm <- proc.time()

########################### set dirs ########################### 
data_Name<- 'ZZZZZZZZ'   ### fill with name of file
mainDir <- "~/Dropbox/papers_dropbox/circuits_2021" ### set this filepath to the main directory 
figDir <- "Figures" ## this is the subdirectory where data is held
dataDir <- "Data"
dataFinalFold <- "Final"
figFolder <-"figFolder" ## this will be created if not already in existence
saveFolder <- "saveFolder" ## this will be created if not already in existence
invisible(  ifelse(!dir.exists(file.path(mainDir, figDir, figFolder)), dir.create(file.path(mainDir, figDir,figFolder)), FALSE))
invisible(  ifelse(!dir.exists(file.path(mainDir, dataDir, saveFolder)), dir.create(file.path(mainDir, dataDir,saveFolder)), FALSE))


########################### load general dependencies ########################### 
packagesS<-c("devtools",
             "dplyr",
             "parallel",
             "ggplot2",
             "readxl",
             "crayon",
             "PerformanceAnalytics",
             "Hmisc",
             "tidyr",
             "GGally",
             "tibble",
             "ggplot2",
             "leaps",
             "caret",
             "rstan",
             "rstanarm",
             "coda",
             "broom",
             "tidybayes",
             "emmeans",
             "loo",
             "bayesplot",
             "magrittr",
             "ggpubr",
             "modelr",
             "broom.mixed",
             "EnvStats",
             "ggrepel",
             "plotrix",
             "STAR",
             "R.matlab",
             "data.table",
             "FactoMineR",
             "factoextra"
             
)

invisible(suppressWarnings(suppressMessages(package.check <- lapply(
  packagesS,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
))))

## do not exclude these options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################################################################
########################### Figure 1 ########################### 
################################################################

########################### load dependencies ###########################
source("~/Dropbox/papers_dropbox/circuits_2021/Code/Final/cust_ggplot.R")
source("~/Dropbox/papers_dropbox/circuits_2021/Code/Final/SummarySE.R")
########################### Custom Functions ###########################
########################### Load Data ###########################
data_Name<- 'All_Neurons.csv'   ### fill with name of file
Original_Data<-read.csv(file.path(mainDir,dataDir,dataFinalFold,data_Name), 
                        header = TRUE, sep = ',',check.names = FALSE)

cols_filter<-c("treatment","treatment_text","autoClass","treatment_class","Cancer","Chemotherapy","Dyn.pfr","Dyn.IB","Stat.mSfr","Stat.afr", "Stat.HzSTD","Stat.HzEnd", "Thr.T","Thr.F", "Thr.L","Stat.LSpkT","Stat.LSpkF","Dyn.spkNum","Stat.spkNum","Dyn.DI","Dyn.slp","Stat.slp","Dyn.F",  "Dyn.pfr1","Dyn.pfr3","Thr.T1","Thr.F1","Thr.L1", "Thr.T3","Thr.F3","Thr.L3","Dyn.slp1","Dyn.slp3" ,"Dyn.spkNum1", "Dyn.spkNum3","Dyn.RDR", "Dyn.IFRdrop")
cols<-c("treatment_class","Dyn.pfr","Dyn.IB","Stat.mSfr","Stat.afr", "Stat.HzSTD","Stat.HzEnd", "Thr.T","Thr.F", "Thr.L","Stat.LSpkT","Stat.LSpkF","Dyn.spkNum","Stat.spkNum","Dyn.DI","Dyn.slp","Stat.slp","Dyn.F",  "Dyn.pfr1","Dyn.pfr3","Thr.T1","Thr.F1","Thr.L1", "Thr.T3","Thr.F3","Thr.L3","Dyn.slp1","Dyn.slp3" ,"Dyn.spkNum1", "Dyn.spkNum3","Dyn.RDR", "Dyn.IFRdrop")

df<-Original_Data[cols_filter]

########################### Data Wrangling ###########################
df<-df %>% filter(treatment_text %in% c("POX", "WT")) %>%
  select(treatment_text, autoClass, Thr.L,Dyn.pfr,Stat.spkNum,Thr.L1 ,Dyn.spkNum,Stat.mSfr) %>%
  mutate(autoClass = factor(autoClass, levels=c("Ia","II","1.5", "Ib")),
         treatment_text = factor(treatment_text, levels = c("WT", "POX")))%>%
  as.data.table()
encoding_params <-c("Thr.L", "Dyn.pfr", "Stat.spkNum", "Thr.L1", "Dyn.spkNum", "Stat.mSfr")

########################### quick visualization ###########################
########################### analyses/modeling ###########################
########################### saving figures ###########################
fig1_plots <- list()
for (i in seq_along(encoding_params)){
  fig1_plots[[i]]<- cust_ggplot(df, "autoClass", encoding_params[i], "treatment_text")
}
names(fig1_plots) <- encoding_params

setwd(file.path(mainDir,figDir,figFolder))
lapply(names(fig1_plots), 
       function(x) ggsave(filename=paste("fig1_",x,".eps",sep=""), plot=fig1_plots[[x]]))

########################### saving data ###########################
df<-Original_Data[cols_filter]
sum_behaviours <- summarySE((gather(df[cols], variable, value, -treatment_class)), measurevar =  "value",
                            groupvar = c("treatment_class", "variable"), na.rm = TRUE)
setwd(file.path(mainDir,dataDir,saveFolder))
write.csv(sum_behaviours, "Summary_Statistics_Figure_1.csv")
########################### Clean up ###########################

rm(fig1_plots, cols, cols_filter, sum_behaviours, i, encoding_params)


################################################################
########################### Figure 2 ########################### 
################################################################


## filter out COX and Pirc
# df_WT_POX <- df%>%filter(treatment_text %in% c('WT', 'POX'))
# 
# df_WT_POX_t<-t(df_WT_POX[cols_filter][,-(1:6)]) ## just WT and POX neurons
# colnames(df_WT_POX_t) <- df_WT_POX$treatment_text
# 
# res.pca <- PCA(t(df_WT_POX_t), graph = FALSE)
# var <- get_pca_var(res.pca)
# eig.val <- get_eigenvalue(res.pca)
# 
# fviz_pca_ind(res.pca,
#              geom.ind = "point", # show points only (nbut not "text")
#              col.ind = rownames(t(df_WT_POX_t)), # color by groups
#              paletteNature2018 = c("#808080","#3e5ce0","#db0415","#933bd3"),
#              addEllipses = TRUE, # Concentration ellipses
#              legend.title = "Groups"
# )
# res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# res.desc$Dim.1
# 
# setwd(file.path(mainDir,figDir,figFolder))
# cairo_ps("Fig2_b.eps")
# fviz_pca_var(res.pca, col.var = "black", repel=T, select.var= list(name = c("Dyn.pfr",  "Stat.spkNum",  "Thr.L", "Thr.L1", "Dyn.spkNum" ,"Stat.mSfr","Dyn.slp3","Dyn.pfr1")))
# invisible(suppressMessages(suppressWarnings(dev.off())))
# 

####### PCA updated 5/31/21 originally retraced steps on 9/28/21 ###### 
df<-Original_Data[cols_filter]

df_WT_POX <- df%>%filter(treatment_text %in% c('WT', 'POX'))
test<-t(df_WT_POX[cols_filter][,-(1:6)]) ## just WT and POX neurons
colnames(test) <- df_WT_POX$treatment_text
require(FactoMineR)
require(factoextra)
res.pca <- PCA(t(test), graph = FALSE)
var <- get_pca_var(res.pca)
eig.val <- get_eigenvalue(res.pca)
### write data
setwd(file.path(mainDir,dataDir,saveFolder))
write.csv(cbind(res.pca$ind$coord, df_WT_POX[,(1:6)]), "pox_wt_res.pca.csv")


## biplot --> must look at the specifics Fig 2c
fviz_pca_var(res.pca, col.var = "black", repel=T, select.var= list(name = c("Dyn.pfr",  "Stat.spkNum",  "Thr.L", "Thr.L1", "Dyn.spkNum" ,"Stat.mSfr","Dyn.slp3","Dyn.pfr1")))
dimdesc(res.pca, axes = c(1,2), proba = 0.05)
PCA(t(test), graph = FALSE)
fviz_pca_biplot(res.pca,
                col.ind = rownames(t(test)), palette = "jco",
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species")

### read data for plot 2a
data_Name<- 'pox_wt_res.pca.csv'   ### fill with name of file
wt_pox_pca_decoding_1<-read.csv(file.path(mainDir,dataDir,saveFolder,data_Name), 
                        header = TRUE, sep = ',',check.names = FALSE)

df_1<-wt_pox_pca_decoding %>% rename(
  V1 = Dim.1,
  V2 = Dim.2,
  V3 = Dim.3,
  merged = treatment_class) %>%
  select(merged, treatment_text, V1, V2) 
# %>% filter(merged !='WT_II') #### used to check the 7th group because shape discrimination only allows 6

### plot 2a 
ggscatter(df_1, x = "V1", y = "V2",
          color = "merged",
          # palette = "npg",
          palette = c("#703B96", "#703B96", "#703B96", "#703B96", "#010101", "#010101", "#010101", "#010101"),
          # shape = "merged",
          # shape = c(15,16,17,18,15,16,17,18),
          ellipse = TRUE,
          ellipse.level = 0.95,
          ellipse.type = "confidence",
          mean.point = TRUE,
          star.plot = F,
          ggtheme = theme_minimal())

### plot Fig 2b
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))


################################################################
########################### Figure 3 ########################### 
################################################################

##### define helper functions #####
flattenlist<- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}

######
WT_data<-readMat("/Users/nickhousley/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/new data?/spike_times_WT.mat") ## control

spike_times_WT_data<-WT_data$spike.times
spike_times_WT_data_flat<-flattenlist(spike_times_WT_data)
spike_times_WT_data_flat<-list(spike_times_WT_data_flat)

res1 <- lapply(spike_times_WT_data_flat, function(x){
  as.repeatedTrain(x)
})

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig3_c.eps")
psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)
invisible(suppressMessages(suppressWarnings(dev.off())))


POX_data<-readMat("/Users/nickhousley/Dropbox/papers_dropbox_/decoding_cancer_chemo/GT_2019_vEnv_1_Cancer_Chemo_Decoding_Project_Folder/Data/spike_times_POX_fixed.mat") ## pox
spike_times_POX_data<-POX_data$spike.times
spike_times_POX_data_flat<-flattenlist(spike_times_POX_data)
spike_times_POX_data_flat<-list(spike_times_POX_data_flat)

res1 <- lapply(spike_times_POX_data_flat, function(x){
  as.repeatedTrain(x)
})

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig3_d.eps")
psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)
invisible(suppressMessages(suppressWarnings(dev.off())))
########################### Clean up ###########################

rm(POX_data,spike_times_POX_data,spike_times_POX_data_flat, WT_data,spike_times_WT_data, spike_times_WT_data_flat, res1)


################################################################
########################### Figure 4 ########################### 
################################################################


################################################################
########################### Figure 5 ########################### 
################################################################


########################### load dependencies ###########################
########################### Custom Functions ###########################
########################### Load Data ###########################
data_Name<- 'ssp.xlsx'   ### fill with name of file
ssp_data<-read_excel(file.path(mainDir,dataDir,dataFinalFold,data_Name), 
                     na = "NA")
########################### Data Wrangling ###########################
########################### quick visualization ###########################
########################### analyses/modeling ###########################
sum_behaviours<-ssp_data %>%
  group_by(treatment) %>%
  filter(ssp_type == "vib_ssp" 
         & delay_type == "first"
         & flag_dont_use != "1"
  ) %>%
  dplyr::select(treatment, mV,delay_type) %>%
  dplyr::mutate(
    mV = mV * 100) %>%
  dplyr::summarise(n = n(),
            n_gt6 = sum(mV < 0.5))

########################### saving data ###########################
setwd(file.path(mainDir,dataDir,saveFolder))
write.csv(sum_behaviours, "summary_ssp_detection.csv")
########################### saving figures ###########################

### b
fig5_b<-ssp_data %>% filter(!is.na(delay) & !delay_type == "delay_ramp")%>%
  ggboxplot( x = "treatment", y = "delay",
             width=0.3,
             bxp.errorbar.width = 1,
             outlier.shape = NA,
             palette = c("#7d4a9c", "#808080"),
             boxwex = 10,  fill = "treatment")+
  facet_wrap(.~delay_type)+
  stat_compare_means(label.y = 550, method = "anova")+
  stat_summary(fun.data = function(x) data.frame(y=-60, label = paste("Mean=",round(mean(x,na.rm=T),digit = 1),"ms")), geom="text") +
  stat_n_text()+ 
  labs(title="Detection Threshold SSPs (4mm/s movement)", y = "delay (ms)", x = "treatment")
setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_b.eps")
fig5_b
invisible(suppressMessages(suppressWarnings(dev.off())))

### c
fig5_c<-ssp_data%>%
  filter(ssp_type == "dynamic_ssp" 
         & delay_type == "delay_tri_first"
  ) %>%
  select(treatment, mV,delay_type) %>%
  mutate(
    mV = mV * 100,
    treatment = factor(treatment, levels=c("control", "pox"))) %>%
  
  ### build the plot
  ggboxplot( x = "treatment", 
             y = "mV",
             width=0.3,
             bxp.errorbar.width = 1,
             palette = c("#808080", "#7d4a9c"),
             outlier.shape = NA,
             boxwex = 10,  
             # ylim = c(-.5,2.5),
             fill = "treatment")+
  stat_compare_means(label.y = 2.5, method = "anova")+
  stat_summary(fun.data = function(x) data.frame(y=0.5, label = paste("Mean=",round(mean(x,na.rm=T),digit = 2),"mV")), geom="text") +
  stat_n_text()+ 
  labs(title="In vivo Decoding", y = "sp (mV)", x = "treatment")

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_c.eps")
fig5_c
invisible(suppressMessages(suppressWarnings(dev.off())))

### e
fig5_e<-ssp_data %>% filter(!is.na(delay) & delay_type == "delay_ramp")%>%
  ggboxplot( x = "treatment", y = "delay",
             width=0.3,
             bxp.errorbar.width = 1,
             outlier.shape = NA,
             palette = c("#808080", "#7d4a9c"),
             boxwex = 10,  fill = "treatment")+
  facet_wrap(.~delay_type)+
  stat_compare_means(label.y = 60, method = "anova")+
  # stat_summary(fun.data = fun_mean, geom="text")+
  stat_summary(fun.data = function(x) data.frame(y=-10, label = paste("Mean=",round(mean(x,na.rm=T),digit = 1),"ms")), geom="text") +
  stat_n_text()+ 
  labs(title="Detection Threshold SSPs (20mm/s movement)", y = "delay (ms)", x = "treatment")

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_e.eps")
fig5_e
invisible(suppressMessages(suppressWarnings(dev.off())))

### f
fig5_f<-ssp_data[!(ssp_data$treatment == "control" & ssp_data$ssp_type == "static_ssp"& ssp_data$mV > 3
             | ssp_data$ssp_type == "dynamic_ssp"& ssp_data$mV > 9 ) ,] %>%
filter(ssp_type == "static_ssp"| ssp_type == "dynamic_ssp" & delay_type != "delay_tri_first") %>%
  select(treatment, mV, ssp_type) %>%
  ### build the plot
  ggboxplot( x = "treatment", 
             y = "mV",
             width=0.3,
             bxp.errorbar.width = 1,
             palette = c("#808080", "#7d4a9c"),
             #add = "jitter",
             outlier.shape = NA,
             boxwex = 10,  
             ylim = c(-.5,2.5),
             fill = "treatment")+
  facet_wrap(.~ssp_type)+
  stat_compare_means(label.y = 2.5, method = "anova")+
  stat_summary(fun.data = function(x) data.frame(y=-0.5, label = paste("Mean=",round(mean(x,na.rm=T),digit = 2),"mV")), geom="text") +
  stat_n_text()+ 
  labs(title="In vivo Decoding", y = "ssp (mV)", x = "treatment")

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_f.eps")
fig5_f
invisible(suppressMessages(suppressWarnings(dev.off())))

### g
fig5_g<-ssp_data[!(ssp_data$treatment == "control" & ssp_data$ssp_type == "static_ssp"& ssp_data$mV > 3
           | ssp_data$ssp_type == "dynamic_ssp"& ssp_data$mV > 9 ) ,] %>%
  filter(ssp_type == "static_ssp"| ssp_type == "dynamic_ssp" & delay_type != "delay_tri_first") %>%
  select(treatment, mV, ssp_type,cell) %>%
  spread(key = ssp_type, value = mV) %>%
  
  
  ggscatter(x = "dynamic_ssp", y = "static_ssp", color = "treatment",
            palette = c("#808080", "#7d4a9c"),
            add = "reg.line",
            conf.int = TRUE)+
  stat_cor(aes(color = treatment)) 
setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_g.eps")
fig5_g
invisible(suppressMessages(suppressWarnings(dev.off())))


### i
fig5_i<-ssp_data%>%
  filter(ssp_type == "vib_ssp" 
         # & delay_type == "first"
         & flag_dont_use != "1"
  ) %>%
  select(treatment, mV,delay_type, animal_simplified, signal_minus_noise) %>%
  mutate(
    signal_minus_noise = signal_minus_noise * 100) %>%
  ### build the plot
  ggstripchart("treatment", "signal_minus_noise",  size = 2, 
               color = "treatment", palette = c("#7d4a9c","#808080" ),
               add = "mean_sd")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")

  facet(fig5_i,facet.by = "delay_type") ## need to remove the "first" dplyr filter

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_i.eps")
# fig5_i
facet(fig5_i,facet.by = "delay_type") ## need to remove the "first" dplyr filter

invisible(suppressMessages(suppressWarnings(dev.off())))

#### get avgs for manuscript narrative ###

ssp_data%>%
  filter(ssp_type == "vib_ssp" 
         # & delay_type == "first"
         & flag_dont_use != "1"
  ) %>%
  select(treatment, mV,delay_type, animal_simplified, signal_minus_noise) %>%
  mutate(
    signal_minus_noise = signal_minus_noise * 100) %>%
  group_by(treatment) %>%
  summarise(mean = mean(signal_minus_noise))



### j #### velocity dependence mv peak amp
data_Name<- 'velocity_dependence.xlsx'   ### fill with name of file
velocity_dependence_data<-read_excel(file.path(mainDir,dataDir,dataFinalFold,data_Name), 
                                     na = "NA")

velocity_dependence_data$stim <- as.factor(velocity_dependence_data$stim)

fig5_j<-velocity_dependence_data%>%
  mutate(stim = factor(stim, levels=c("tri","vib","ramp", "tenTap"))) %>%
  ggstripchart("stim", "mV_ratio",  size = 2, 
               color = "stim", 
               add = "mean_sd")+
  scale_colour_brewer(palette = "Purples")

setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig5_j.eps")
fig5_j
invisible(suppressMessages(suppressWarnings(dev.off())))


### k   #### velocity dependence delay
select_crit <- c("tri", "ramp", "tenTap")
fig5_k<-velocity_dependence_data%>%
  filter(stim %in% select_crit) %>%
  mutate(stim = factor(stim, levels=c("tri", "ramp", "tenTap"))) %>%
  ggstripchart("stim", "delay_ratio",  size = 2, 
               color = "stim", 
               add = "mean_sd")+
  scale_colour_brewer(palette = "Purples")

cairo_ps("fig5_k.eps")
fig5_k
invisible(suppressMessages(suppressWarnings(dev.off())))

########################### Clean up ###########################



################################################################
########################### Figure 6 ########################### 
################################################################

## c
setwd(file.path(mainDir,figDir,figFolder))
controlData <- data.frame(group=c("Miss 1.8 ±  2%", "Hit 98.2 ±  2%"), FR=c(1.8, 98.2))
cairo_ps("fig6_c.eps")
pie3D(controlData$FR, labels = controlData$group, main = "An exploded 3D pie chart", explode=0.1, radius=.8, labelcex = 1.2,  start=.5)
invisible(suppressMessages(suppressWarnings(dev.off())))

## d
PoxData <- data.frame(group=c("Miss 19.5 ±  5%", "Hit 80.5 ±  2%"), FR=c(19.5, 80.5))
cairo_ps("fig6_d.eps")
pie3D(PoxData$FR, labels = PoxData$group, main = "An exploded 3D pie chart", explode=0.1, radius=.8, labelcex = 1.2,  start=.2)
invisible(suppressMessages(suppressWarnings(dev.off())))
rm(controlData, PoxData)

## e
source("~/Dropbox/papers_dropbox/circuits_2021/Code/Final/kinematic_replacement_errors.R")
cairo_ps("fig6_e.eps")
fig6_e
invisible(suppressMessages(suppressWarnings(dev.off())))

## f
cairo_ps("fig6_f.eps")
fig6_f
invisible(suppressMessages(suppressWarnings(dev.off())))

################################################################
########################### Figure 7 ########################### 
################################################################

## c1
source("~/Dropbox/papers_dropbox/circuits_2021/Code/Final/biomechanical_cor_with_pop.R")
setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig7_c1.eps")
fig7_c1
invisible(suppressMessages(suppressWarnings(dev.off())))


## c2
source("~/Dropbox/papers_dropbox/circuits_2021/Code/Final/SPs_Correlation_with_population_code.R")
setwd(file.path(mainDir,figDir,figFolder))
cairo_ps("fig7_c2.eps")
fig7_c2
invisible(suppressMessages(suppressWarnings(dev.off())))

################################################################
########################### Supplemental #######################
################################################################



################################################################
########################### References #########################
################################################################


################################################################
########################### testing/un-used ####################
################################################################


# ### signal to noise ratio from Figure 5
# data_Name<- 'ssp.xlsx'   ### fill with name of file
# ssp_data<-read_excel(file.path(mainDir,dataDir,dataFinalFold,data_Name), 
#                      na = "NA")
# fig5_unused_1<-ssp_data%>%
#   filter(ssp_type == "vib_ssp" 
#          & delay_type == "first"
#          & flag_dont_use != "1"
#   ) %>%
#   select(treatment, sig_to_noise) %>%
#   ggboxplot( x = "treatment", 
#              y = "sig_to_noise",
#              width=0.3,
#              bxp.errorbar.width = 1,
#              palette = c("#7d4a9c", "#808080"),
#              #add = "jitter",
#              outlier.shape = NA,
#              boxwex = 10,  
#              # ylim = c(-.5,2.5),
#              fill = "treatment")+
#   stat_compare_means(label.y = 2.5, method = "anova")
# setwd(file.path(mainDir,figDir,figFolder))
# cairo_ps("fig5_unused_1.eps")
# fig5_unused_1
# invisible(suppressMessages(suppressWarnings(dev.off())))
# 
# ### tendon Taps from Figure 5
# fig5_unused_2<-ssp_data%>%
#   filter(ssp_type == "tendon_tap_ssp" 
#   ) %>%
#   select(treatment, mV, animal, delay_type, delay) %>%
#   
#   ## build the plot
#   ggstripchart("treatment", "mV",  size = 2,
#                # shape = "treatment",
#                color = "treatment",
#                palette = c("#7d4a9c","#808080" ),
#                add = "mean_sd")+
#   geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")
# cairo_ps("fig5_unused_2.eps")
# fig5_unused_2
# invisible(suppressMessages(suppressWarnings(dev.off())))




########################### timer stop ###########################
proc.time() - ptm


