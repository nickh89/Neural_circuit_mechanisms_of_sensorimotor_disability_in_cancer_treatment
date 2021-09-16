#### Read in data
require(R.matlab)
require(STAR)

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
c<-readMat("/Users/nickhousley/Dropbox (GaTech)/papers_dropbox_/decoding_cancer_chemo/new data?/spike_times_WT.mat") ## control

d<-c$spike.times
ddd<-flattenlist(d)
neuron<-list(ddd)

res1 <- lapply(neuron, function(x){
  as.repeatedTrain(x)
})

fig3_c<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)


c<-readMat("/Users/nickhousley/Dropbox/papers_dropbox_/decoding_cancer_chemo/GT_2019_vEnv_1_Cancer_Chemo_Decoding_Project_Folder/Data/spike_times_POX_fixed.mat") ## pox
d<-c$spike.times
ddd<-flattenlist(d)
neuron<-list(ddd)

res1 <- lapply(neuron, function(x){
  as.repeatedTrain(x)
})

fig3_d<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01),colCI=2)



# 
# psth(res1[[1]],breaks=c(bw=0.5,step=0.05),colCI=2)
# 
# psth1<-psth(res1[[1]],breaks=c(bw=0.1,step=0.01), plot = F)
# 
# maxFreq <- psth1$freq
# 
# 
# ## First get lists containing PSTHs from each neuron
# psth1 <- psth(CAL1V[["neuron 1"]],breaks=c(bw=0.5,step=0.05),plot=FALSE)
# psth2 <- psth(CAL1V[["neuron 2"]],breaks=c(bw=1,step=0.1),plot=FALSE)
# psth3 <- psth(CAL1V[["neuron 3"]],breaks=c(bw=0.5,step=0.05),plot=FALSE)
# psth4 <- psth(CAL1V[["neuron 4"]],breaks=c(bw=2,step=0.2),plot=FALSE)
# 
# ## Get the maximal frequency to display
# maxFreq <- max(max(psth1$ciUp),max(psth2$ciUp),max(psth3$ciUp),max(psth4$ciUp))
# ## Build plot
# plot(c(0,10),c(0,75),type="n",
#      xaxs="i",yaxs="i",xlab="Time (s)",
#      ylab="Freq. (Hz)",
#      main="PSTHs from 4 simultaneously recorded neurons",
#      sub="20 stimulations with vanillin were used.")
# ## Add rectangle corresponding to stimulation command
# rect(4.49,0,4.99,75,col="grey80",lty=0)
# ## Add the neurons PSTHs as confidence bands
# polygon(c(psth1$mids,rev(psth1$mids)),c(psth1$ciLow,rev(psth1$ciUp)),col=1,border=NA)
# polygon(c(psth2$mids,rev(psth2$mids)),c(psth2$ciLow,rev(psth2$ciUp)),col=2,border=NA)
# polygon(c(psth3$mids,rev(psth3$mids)),c(psth3$ciLow,rev(psth3$ciUp)),col=3,border=NA)
# polygon(c(psth4$mids,rev(psth4$mids)),c(psth4$ciLow,rev(psth4$ciUp)),col=4,border=NA)
# legend(0.1,maxFreq,legend=paste("neuron",1:4),lty=1,col=1:4,bty="n")
# 
# c$spike.times
# 
# ddd[1] %>% filter [ddd[1] > 3]
# 
# test <-
# 
# ddd[[1]][!rowSums(ddd[[1]] >20),]
# x[[i]][!rowSums(x[[i]] >20),]
# 
# 
# res2 <- lapply(ddd, function(x){
#   x[[]][!rowSums(x[[i]] >20),]
# })
# 
# 
# 
# 
# ff1 = function(x) 
# {
#   ans = numeric(length(x))
#   
#   for(i in seq_along(x)) ans[i] = x[[i]][!rowSums(x[[i]] >20),]
# 
#   return(ans)
# }
# 
# for(i in ddd) 
# {
#   ddd[[i]][!rowSums(ddd[[i]] >20),]
#   
# }
# 
# 
