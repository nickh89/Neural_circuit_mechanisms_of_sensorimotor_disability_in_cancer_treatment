#####    this is to extract the population level encoding (afferents) from the 40 neurons 
#####   to run a correlation analysis with the SPs on both control and pox data 
# rm(list = ls())
# gc()

invisible(suppressWarnings(suppressMessages(lapply(c("readxl","tidyr", "ggplot2", "dplyr","magrittr", "readr","PerformanceAnalytics","Hmisc","data.table","ggpubr","EnvStats", "R.matlab", "STAR"), require, character.only = TRUE))))

######## ######## Read in spike data from matlab for POX
c<-readMat("/Users/nickhousley/Dropbox/papers_dropbox_/decoding_cancer_chemo/GT_2019_vEnv_1_Cancer_Chemo_Decoding_Project_Folder/Data/spike_times_POX_fixed.mat") ## pox

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

######## ######## reshape the df for 'STAR' to run PSTH on. 
d<-c$spike.times
ddd<-flattenlist(d)
neuron<-list(ddd)

res1 <- lapply(neuron, function(x){
  as.repeatedTrain(x)
})

# psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)   ### preferred plot lower bin size
# psth(res1[[1]],breaks=c(bw=0.5,step=0.05),colCI=2)

psth1<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01), plot = F)
meanFreq <- psth1$freq
rm(c,d,ddd,neuron, res1, psth1)

######################## read in POX ssp voltage ramp data #########################

setwd("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/POX/")
file.ls <- list.files(path='~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/POX/untitled folder/', pattern=c('_r','.txt'))
times <- read.table(file.ls[1], header=TRUE, sep="\t")[,1]     # gene names
df_ramp    <- do.call(cbind,lapply(file.ls,function(fn)read.table(fn,header=TRUE, sep="\t")[,2]))
colnames(df_ramp) <- sub(".txt", "", file.ls)
df_ramp    <- as.data.frame(cbind(times,df_ramp))
rm(file.ls, times)


######################## read in Stim data #########################
ramp_stimuli <- read_excel("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/length_force_stiffness.xlsx",col_names = c("time","motor_l","motor_f","motor_s"),skip=1)

suppressMessages(suppressWarnings(avg_length_tri <- read_delim("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/avg_length_tri.txt", 
                                                               "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("time","motor_l"),skip=1)))
suppressMessages(suppressWarnings(avg_force_tri <- read_delim("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/avg_force_tri.txt", 
                                                              "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("time","motor_f"),skip=1)))
suppressMessages(suppressWarnings(tri_stimuli <- cbind(avg_length_tri, avg_force_tri[2])))
suppressMessages(suppressWarnings(tri_stimuli <- cbind(tri_stimuli, motor_s=tri_stimuli$motor_f/tri_stimuli$motor_l)))

rm(avg_length_tri,avg_force_tri)



######## ######## ramps POX
ramp_popHz<-meanFreq[c(522:670)]
plot(ramp_popHz,type="l",col="red")
par(new=TRUE)
plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

# ramp_popHz<-meanFreq[c(715:863)]
# plot(ramp_popHz,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

######## ######## need to upsample the stim data because it was low freq calculation: use spline interpolation
n <- length(ramp_popHz)
n1<-nrow(ramp_stimuli)

splineData_ramp <- data.frame(
  with(as.data.frame(ramp_popHz), 
       spline(ramp_popHz, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

# ######## ######## plot the mechanicas and population code to make sure they are aligned
# plot(splineData_ramp$y,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

######## ######## cor of stim to pop code should be ~40%
res_ramp<-rcorr(as.matrix(bind_cols(splineData_ramp[-1,c(2)],df_ramp[, c(2:length(df_ramp))])))
res_ramp[["r"]][2:12]
mean(res_ramp[["r"]][2:12])

######## ######## clean up the df and rename the variables for easy binding later
pox_ssp_pop_ramp_r2<-as.data.frame(res_ramp$r[c(2:ncol(res_ramp$r)),c(1:1)])
names(pox_ssp_pop_ramp_r2)[1] <- "r2"
pox_ssp_pop_ramp_r2['treatment'] = "pox"
pox_ssp_pop_ramp_r2['stim'] = "ramp"


######################## read in POX ssp voltage tri data #########################

setwd("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/POX/")
file.ls <- list.files(path='~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/POX/untitled folder/', pattern=c('_t','.txt'))
times <- read.table(file.ls[1], header=TRUE, sep="\t")[,1]     # gene names
df_tri    <- do.call(cbind,lapply(file.ls,function(fn)read.table(fn,header=TRUE, sep="\t")[,2]))
colnames(df_tri) <- sub(".txt", "", file.ls)
df_tri    <- as.data.frame(cbind(times,df_tri))


######## ######## triangles 
tri_popHz<-meanFreq[c(0:470)]
plot(tri_popHz,type="l",col="red")
par(new=TRUE)
plot(tri_stimuli$motor_l,col="green", axes=T, type="l")


######## ######## need to upsample the stim data because it was low freq calculation: use spline interpolation
n <- length(tri_popHz)
n1<-nrow(tri_stimuli)
splineData_tri <- data.frame(
  with(as.data.frame(tri_popHz), 
       spline(tri_popHz, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

# ######## ######## plot the mechanicas and population code to make sure they are aligned
# plot(splineData_tri$y,type="l",col="red")
# par(new=TRUE)
# plot(tri_stimuli$motor_f,col="green", axes=T, type="l")

######## ######## cor of stim to pop code should be  ~52%
res_tri<-rcorr(as.matrix(bind_cols(splineData_tri[-1,c(2)],df_tri[, c(2:length(df_tri))])))
res_tri[["r"]][2:4]
mean(res_tri[["r"]][2:4])

######## ######## clean up the df and rename the variables for easy binding later
pox_ssp_pop_tri_r2<-as.data.frame(res_tri$r[c(2:ncol(res_tri$r)),c(1:1)])
names(pox_ssp_pop_tri_r2)[1] <- "r2"
pox_ssp_pop_tri_r2['treatment'] = "pox"
pox_ssp_pop_tri_r2['stim'] = "tri"

######## ######## bind both ramp and triangles together
pox_r2<-rbind(pox_ssp_pop_ramp_r2,pox_ssp_pop_tri_r2)
# rm(list=setdiff(ls(), "pox_r2"))










########################## Read in spike data from matlab for Control ##################
c<-readMat("/Users/nickhousley/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/new data?/spike_times_WT.mat") ## control

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

######## ######## reshape the df for 'STAR' to run PSTH on. 
res1 <- lapply(list(flattenlist(c$spike.times)), function(x){
  as.repeatedTrain(x)
})

######## ######## graphic for population code 
# psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)   ### preferred plot lower bin size
# psth(res1[[1]],breaks=c(bw=0.5,step=0.05),colCI=2)

psth1<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01), plot = F)
meanFreq <- psth1$freq
rm(c,d,ddd,neuron, res1, psth1)

######################## read in ramp control data #########################

setwd("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/Control/")
file.ls <- list.files(path='~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/Control/untitled folder/', pattern=c('_r','.txt'))
times <- read.table(file.ls[1], header=TRUE, sep="\t")[,1]     # gene names
suppressMessages(suppressWarnings(df    <- do.call(cbind,lapply(file.ls,function(fn)read.table(fn,header=TRUE, sep="\t")[,2]))))
colnames(df) <- sub(".txt", "", file.ls)
df_ramp_control    <- as.data.frame(cbind(times,df))

######################## read in ramp control data #########################
setwd("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/Control/")
file.ls <- list.files(path='~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/Control/untitled folder/', pattern=c('_t','.txt'))
suppressMessages(suppressWarnings(times <- read.table(file.ls[1], header=TRUE, sep="\t")[,1]))   
suppressMessages(suppressWarnings(df    <- do.call(cbind,lapply(file.ls,function(fn)read.table(fn,header=TRUE, sep="\t")[,2]))))
colnames(df) <- sub(".txt", "", file.ls)
df_tri_control    <- as.data.frame(cbind(times,df))



######################## read in Stim data #########################
ramp_stimuli <- read_excel("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/length_force_stiffness.xlsx",col_names = c("time","motor_l","motor_f","motor_s"),skip=1)

suppressMessages(suppressWarnings(avg_length_tri <- read_delim("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/avg_length_tri.txt", 
                                                               "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("time","motor_l"),skip=1)))
suppressMessages(suppressWarnings(avg_force_tri <- read_delim("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/avg_force_tri.txt", 
                                                              "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("time","motor_f"),skip=1)))
suppressMessages(suppressWarnings(tri_stimuli <- cbind(avg_length_tri, avg_force_tri[2])))
suppressMessages(suppressWarnings(tri_stimuli <- cbind(tri_stimuli, motor_s=tri_stimuli$motor_f/tri_stimuli$motor_l)))

rm(avg_length_tri,avg_force_tri)





######## ######## ramps for population code CONTROL
ramp_popHz<-meanFreq[c(522:670)]
plot(ramp_popHz,type="l",col="red")
par(new=TRUE)
plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")


######## ######## need to upsample the stim data because it was low freq calculation: use spline interpolation
n <- length(ramp_popHz)
n1<-nrow(ramp_stimuli)

splineData_ramp <- data.frame(
  with(as.data.frame(ramp_popHz), 
       spline(ramp_popHz, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)


######## ######## cor of stim to pop code should be ~88%
res_ramp<-rcorr(as.matrix(bind_cols(splineData_ramp[-1,c(2)],df_ramp_control[, c(2:length(df_ramp_control))])))
res_ramp[["r"]][2:11]
mean(res_ramp[["r"]][2:11])

######## ######## clean up the df and rename the variables for easy binding later
control_ssp_pop_ramp_r2<-as.data.frame(res_ramp[["r"]][2:11])
names(control_ssp_pop_ramp_r2)[1] <- "r2"
control_ssp_pop_ramp_r2['treatment'] = "control"
control_ssp_pop_ramp_r2['stim'] = "ramp"


######## ######## triangles 
tri_popHz<-meanFreq[c(0:480)]
plot(tri_popHz,type="l",col="red")
par(new=TRUE)
plot(tri_stimuli$motor_l,col="green", axes=T, type="l")

######## ######## need to upsample the stim data because it was low freq calculation: use spline interpolation
n <- length(tri_popHz)
n1<-nrow(tri_stimuli)
splineData_tri <- data.frame(
  with(as.data.frame(tri_popHz), 
       spline(tri_popHz, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

######## ######## cor of stim to pop code should be  ~81%
res_tri<-rcorr(as.matrix(bind_cols(splineData_tri[-1,c(2)],df_tri_control[, c(2:length(df_tri_control))])))
res_tri[["r"]][2:6]
mean(res_tri[["r"]][2:6])

######## ######## clean up the df and rename the variables for easy binding later
control_ssp_pop_tri_r2<-as.data.frame(res_tri[["r"]][2:6])
names(control_ssp_pop_tri_r2)[1] <- "r2"
control_ssp_pop_tri_r2['treatment'] = "control"
control_ssp_pop_tri_r2['stim'] = "tri"

######## ######## bind both ramp and triangles together
control_r2<-rbind(control_ssp_pop_ramp_r2,control_ssp_pop_tri_r2)

# rm(list=ls()[! ls() %in% c("pox_r2","control_r2")])


######## ######## build graphic ######## ########

fig7_c2<-rbind(pox_r2,control_r2)%>%
  ### build the plot
  ggboxplot( x = "treatment", 
             y = "r2",
             width=0.3,
             bxp.errorbar.width = 1,
             palette = c("#7d4a9c","#808080"),
             #add = "jitter",
             outlier.shape = NA,
             boxwex = 10,  
             # ylim = c(0,1.2),
             fill = "treatment")+
  # facet_wrap(.~mechan_param)+
  stat_compare_means(label.y = 1, method = "anova")+
  # stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) + #### put a point where the mean is
  #stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7)+          #### use the above fun_mean function
  stat_summary(fun.data = function(x) data.frame(y=0, label = paste("Mean=",round(mean(x,na.rm=T),digit = 2),"%")), geom="text") +
  # stat_summary(fun.data = function(x) data.frame(y=(mean(x)+1),label=mean(x,na.rm=T)),digit=2),geom="text") +   #### alternative to immediately above with dynamically controlled y axis placement "(mean(x)+1)" 
  stat_n_text()+ 
  labs(title="Population Code Correlation with ssp", y = "R2 (%)", x = "treatment")



