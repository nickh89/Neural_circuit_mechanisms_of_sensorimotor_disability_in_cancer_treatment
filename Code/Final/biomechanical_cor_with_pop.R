### biomechanical cor with pop code

# rm(list = ls())
# gc()

invisible(suppressWarnings(suppressMessages(lapply(c("readxl","tidyr", "ggplot2", "dplyr","magrittr", "readr","PerformanceAnalytics","Hmisc","data.table","ggpubr","EnvStats"), require, character.only = TRUE))))

#### Read in spike data
require(R.matlab)
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

d<-c$spike.times
ddd<-flattenlist(d)
neuron<-list(ddd)

require(STAR)
res1 <- lapply(neuron, function(x){
  as.repeatedTrain(x)
})

# psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)   ### preferred plot lower bin size
# psth(res1[[1]],breaks=c(bw=0.5,step=0.05),colCI=2)

psth1<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01), plot = F)
meanFreq <- psth1$freq
rm(c,d,ddd,neuron, res1, psth1)

######################## read in Stim data #########################
ramp_stimuli <- read_excel("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/length_force_stiffness.xlsx",col_names = c("time","motor_l","motor_f","motor_s"),skip=1)

suppressMessages(suppressWarnings(avg_length_tri <- read_delim("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/avg_length_tri.txt", 
                                                               "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("time","motor_l"),skip=1)))
suppressMessages(suppressWarnings(avg_force_tri <- read_delim("~/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/Data/ssp_textFiles/stimulus_ground_truth/avg_force_tri.txt", 
                                                              "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("time","motor_f"),skip=1)))
suppressMessages(suppressWarnings(tri_stimuli <- cbind(avg_length_tri, avg_force_tri[2])))
suppressMessages(suppressWarnings(tri_stimuli <- cbind(tri_stimuli, motor_s=tri_stimuli$motor_f/tri_stimuli$motor_l)))

rm(avg_length_tri,avg_force_tri)

##148 

ramp_popHz_1<-meanFreq[c(522:670)]
# plot(ramp_popHz_1,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")
# 

n <- length(ramp_popHz_1)
n1<-nrow(ramp_stimuli)

ramp_popHz_1 <- data.frame(
  with(as.data.frame(ramp_popHz_1), 
       spline(ramp_popHz_1, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)



ramp_popHz_2<-meanFreq[c(715:863)]
# plot(ramp_popHz_2,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_2 <- data.frame(
  with(as.data.frame(ramp_popHz_2), 
       spline(ramp_popHz_2, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_3<-meanFreq[c(908:1056)]
# plot(ramp_popHz_3,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_3 <- data.frame(
  with(as.data.frame(ramp_popHz_3), 
       spline(ramp_popHz_3, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_4<-meanFreq[c(1090:1238)]
# plot(ramp_popHz_4,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_4 <- data.frame(
  with(as.data.frame(ramp_popHz_4), 
       spline(ramp_popHz_4, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_5<-meanFreq[c(1290:1438)]
# plot(ramp_popHz_5,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_5 <- data.frame(
  with(as.data.frame(ramp_popHz_5), 
       spline(ramp_popHz_5, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_6<-meanFreq[c(1480:1628)]
# plot(ramp_popHz_6,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_6 <- data.frame(
  with(as.data.frame(ramp_popHz_6), 
       spline(ramp_popHz_6, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)


ramp_popHz_7<-meanFreq[c(1673:1821)]
# plot(ramp_popHz_7,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_7 <- data.frame(
  with(as.data.frame(ramp_popHz_7), 
       spline(ramp_popHz_7, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_multipl<-cbind(ramp_popHz_1$y,ramp_popHz_2$y, ramp_popHz_3$y, ramp_popHz_4$y, ramp_popHz_5$y, ramp_popHz_6$y, ramp_popHz_7$y)
ramp_multipl<-as.data.frame(ramp_multipl)
rm(ramp_popHz_1,ramp_popHz_2,ramp_popHz_3,ramp_popHz_4, ramp_popHz_5, ramp_popHz_6, ramp_popHz_7)





## cor of stim to pop code  55%
res_biomech_pop<-rcorr(as.matrix(bind_cols(ramp_multipl[-1,],ramp_stimuli[, c(2:length(ramp_stimuli))])))
res_biomech_pop[["r"]][1:7, c(8:10)]
mean(res_biomech_pop[["r"]][1:7, c(8:10)])

res_biomech_pop_r2<-as.data.frame(res_biomech_pop$r[c(1:7),c(8:10)])
res_biomech_pop_r2<-gather(res_biomech_pop_r2, value = "r2" ,key = "mechan_param" )
res_biomech_pop_r2['treatment'] = "POX"














### |||||||||||| ########
c<-readMat("/Users/nickhousley/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/new data?/spike_times_WT.mat") ## control


d<-c$spike.times
ddd<-flattenlist(d)
neuron<-list(ddd)

require(STAR)
res1 <- lapply(neuron, function(x){
  as.repeatedTrain(x)
})

# psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)   ### preferred plot lower bin size
# psth(res1[[1]],breaks=c(bw=0.5,step=0.05),colCI=2)

psth1<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01), plot = F)
meanFreq <- psth1$freq
rm(c,d,ddd,neuron, res1, psth1)


ramp_popHz_1<-meanFreq[c(520:670)]
# plot(ramp_popHz_1,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")


n <- length(ramp_popHz_1)
n1<-nrow(ramp_stimuli)

ramp_popHz_1 <- data.frame(
  with(as.data.frame(ramp_popHz_1), 
       spline(ramp_popHz_1, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)


ramp_popHz_2<-meanFreq[c(713:863)]
# plot(ramp_popHz_2,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_2 <- data.frame(
  with(as.data.frame(ramp_popHz_2), 
       spline(ramp_popHz_2, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_3<-meanFreq[c(900:1050)]
# plot(ramp_popHz_3,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_3 <- data.frame(
  with(as.data.frame(ramp_popHz_3), 
       spline(ramp_popHz_3, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_4<-meanFreq[c(1090:1240)]
# plot(ramp_popHz_4,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_4 <- data.frame(
  with(as.data.frame(ramp_popHz_4), 
       spline(ramp_popHz_4, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_5<-meanFreq[c(1280:1430)]
# plot(ramp_popHz_5,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_5 <- data.frame(
  with(as.data.frame(ramp_popHz_5), 
       spline(ramp_popHz_5, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_popHz_6<-meanFreq[c(1472:1622)]
# plot(ramp_popHz_6,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_6 <- data.frame(
  with(as.data.frame(ramp_popHz_6), 
       spline(ramp_popHz_6, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)


ramp_popHz_7<-meanFreq[c(1663:1813)]
# plot(ramp_popHz_7,type="l",col="red")
# par(new=TRUE)
# plot(ramp_stimuli$motor_f,col="green", axes=T, type="l")

ramp_popHz_7 <- data.frame(
  with(as.data.frame(ramp_popHz_7), 
       spline(ramp_popHz_7, xout = seq(1, n, by = ((n-1)/n1)))
  ),
  method = "spline()"
)

ramp_multipl<-cbind(ramp_popHz_1$y,ramp_popHz_2$y, ramp_popHz_3$y, ramp_popHz_4$y, ramp_popHz_5$y, ramp_popHz_6$y, ramp_popHz_7$y)
ramp_multipl<-as.data.frame(ramp_multipl)
rm(ramp_popHz_1,ramp_popHz_2,ramp_popHz_3,ramp_popHz_4, ramp_popHz_5, ramp_popHz_6, ramp_popHz_7)





## cor of stim to pop code  91%
res_biomech_pop<-rcorr(as.matrix(bind_cols(ramp_multipl[-1,],ramp_stimuli[, c(2:length(ramp_stimuli))])))
res_biomech_pop[["r"]][1:7, c(8:10)]
mean(res_biomech_pop[["r"]][1:7, c(8:10)])

res_biomech_pop_r2_control<-as.data.frame(res_biomech_pop$r[c(1:7),c(8:10)])
res_biomech_pop_r2_control<-gather(res_biomech_pop_r2_control, value = "r2" ,key = "mechan_param" )
res_biomech_pop_r2_control['treatment'] = "control"





fig7_c1<-rbind(res_biomech_pop_r2,res_biomech_pop_r2_control)%>%
  ### build the plot
  ggboxplot( x = "treatment", 
             y = "r2",
             width=0.3,
             bxp.errorbar.width = 1,
             palette = c("#7d4a9c","#808080"),
             #add = "jitter",
             outlier.shape = NA,
             boxwex = 10,  
             # ylim = c(0,1.1),
             fill = "treatment")+
  # facet_wrap(.~mechan_param)+
  stat_compare_means(label.y = 1, method = "anova")+
  # stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) + #### put a point where the mean is
  #stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7)+          #### use the above fun_mean function
  stat_summary(fun.data = function(x) data.frame(y=0, label = paste("Mean=",round(mean(x,na.rm=T),digit = 2),"%")), geom="text") +
  # stat_summary(fun.data = function(x) data.frame(y=(mean(x)+1),label=mean(x,na.rm=T)),digit=2),geom="text") +   #### alternative to immediately above with dynamically controlled y axis placement "(mean(x)+1)" 
  stat_n_text()+ 
  labs(title="Population Code Correlation with Mechanics", y = "R2 (%)", x = "treatment")

rm(ramp_multipl, ramp_stimuli,res_biomech_pop, res_biomech_pop_r2, res_biomech_pop_r2_control, tri_stimuli, meanFreq, n, n1,flattenlist)

