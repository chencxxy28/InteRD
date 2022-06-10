#setwd("/Users/chixiangchen/Desktop/postdoc/deconvolution/multi-robust/from_cluster/multirobust/toast/lessmg10_seger/ww_50/")
#setwd("/storage/home/bzs438/chixiang/deconvolution/interd2_new/baronscrna/ref_free/")
#setwd("/storage/home/bzs438/chixiang/deconvolution/interd2/baronscrna/n40marker20/true_pbar")

#setwd("/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar")
setwd("/Users/chixiang.chen/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/postdoc/postdoc/deconvolution/ref_based_rd/simulations/simulation_09272021/interd2_new/xinscrna/true_pbar/")

meth.names<-c("interd1","interd2","scdc_baron","scdc_xin","scdc_ensemble_mad_all","toast")
#comb.names<-c("n25m20","n40m20","n25m40","n40m40","n80m20","n80m40")
#comb.names.true<-c("n25mr20","n40mr20","n25mr40","n40mr40","n80mr20","n80mr40")
comb.names<-c("n80m40")
comb.names.true<-c("n80mr40")
all_mad_list<-list()
all_cor_list<-list()
for (i in 1:length(comb.names))
{
  comb.names.used<-comb.names[i]
  comb.names.used.true<-comb.names.true[i]
  mad_eachcell_mean_all<-rep()
  cor_eachcell_mean_all<-rep()
  mad_each_all<-rep()
  cor_each_all<-rep()
  for(j in meth.names )
  {
    names.used<-paste0(".*",j,comb.names.used,".*\\.csv$")
    temp = list.files(pattern=names.used,full.names=TRUE)
    myfiles = lapply(temp, read.csv)
    data_est<-do.call(cbind,myfiles)
    data_est<-data_est[,-which(colnames(data_est)=="X")]
    names.used.true<-paste0(".*truep",comb.names.used.true,".*\\.csv$")
    temp = list.files(pattern=names.used.true,full.names=TRUE)
    myfiles = lapply(temp, read.csv)
    data_truep<-do.call(cbind,myfiles)
    data_truep<-data_truep[,-which(colnames(data_truep)=="X")]
    
    #mad
    mad_eachrun<-apply(abs(data_est-data_truep),2,mean)
    names(mad_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=100)
    # mad_eachcell_mean<-tapply(mad_eachrun,names(mad_eachrun),mean)
    # mad_eachcell_summary<-tapply(mad_eachrun,names(mad_eachrun),summary)
    # mad_eachcell_mean<-c(mad_eachcell_mean,mean(mad_eachcell_mean))
    # mad_eachcell_mean_all<-rbind(mad_eachcell_mean_all,i=mad_eachcell_mean)
    mad_each_all<-rbind(mad_each_all,mad_eachrun)
    
    cor_eachrun<-unlist(sapply(1:ncol(data_est),function (x)
    {
      cor(data_truep[,x],data_est[,x],method = c("kendall"))
    }
    ))
    names(cor_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=100)
    # cor_eachcell_mean<-tapply(cor_eachrun,names(cor_eachrun),mean)
    # cor_eachcell_summary<-tapply(cor_eachrun,names(cor_eachrun),summary)
    # cor_eachcell_mean<-c(cor_eachcell_mean,mean(cor_eachcell_mean))
    # cor_eachcell_mean_all<-rbind(cor_eachcell_mean_all,i=cor_eachcell_mean)
    cor_each_all<-rbind(cor_each_all,cor_eachrun)
  }
  #colnames(mad_eachcell_mean_all)<-colnames(cor_eachcell_mean_all)<-c("alpha","beta","delta","gamma","averaged")
  rownames(mad_each_all)<-rownames(cor_each_all)<-c("InteRD1",
                                                                      "InteRD2",
                                                                      "SCDC-Baron",
                                                                      "SCDC-Xin",
                                                                      "SCDC-ENSEMBLE",
                                                                      "TOAST")
}


setwd("/Users/chixiang.chen/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/postdoc/postdoc/deconvolution/ref_based_rd/simulations/simulation_09272021/interd2_new/xinscrna/wrong_pbar/")

meth.names<-c("interd2","toast")
#comb.names<-c("n25m20","n40m20","n25m40","n40m40","n80m20","n80m40")
#comb.names.true<-c("n25mr20","n40mr20","n25mr40","n40mr40","n80mr20","n80mr40")
# comb.names<-c("n25m20")
# comb.names.true<-c("n25mr20")
all_mad_list<-list()
all_cor_list<-list()
for (i in 1:length(comb.names))
{
  comb.names.used<-comb.names[i]
  comb.names.used.true<-comb.names.true[i]
  mad_eachcell_mean_all<-rep()
  cor_eachcell_mean_all<-rep()
  mad_each_all2<-rep()
  cor_each_all2<-rep()
  for(j in meth.names )
  {
    names.used<-paste0(".*",j,comb.names.used,".*\\.csv$")
    temp = list.files(pattern=names.used,full.names=TRUE)
    myfiles = lapply(temp, read.csv)
    data_est<-do.call(cbind,myfiles)
    data_est<-data_est[,-which(colnames(data_est)=="X")]
    names.used.true<-paste0(".*truep",comb.names.used.true,".*\\.csv$")
    temp = list.files(pattern=names.used.true,full.names=TRUE)
    myfiles = lapply(temp, read.csv)
    data_truep<-do.call(cbind,myfiles)
    data_truep<-data_truep[,-which(colnames(data_truep)=="X")]
    
    #mad
    mad_eachrun<-apply(abs(data_est-data_truep),2,mean)
    names(mad_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=100)
    # mad_eachcell_mean<-tapply(mad_eachrun,names(mad_eachrun),mean)
    # mad_eachcell_summary<-tapply(mad_eachrun,names(mad_eachrun),summary)
    # mad_eachcell_mean<-c(mad_eachcell_mean,mean(mad_eachcell_mean))
    # mad_eachcell_mean_all<-rbind(mad_eachcell_mean_all,i=mad_eachcell_mean)
    mad_each_all2<-rbind(mad_each_all2,mad_eachrun)
    
    cor_eachrun<-unlist(sapply(1:ncol(data_est),function (x)
    {
      cor(data_truep[,x],data_est[,x],method = c("kendall"))
    }
    ))
    names(cor_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=100)
    # cor_eachcell_mean<-tapply(cor_eachrun,names(cor_eachrun),mean)
    # cor_eachcell_summary<-tapply(cor_eachrun,names(cor_eachrun),summary)
    # cor_eachcell_mean<-c(cor_eachcell_mean,mean(cor_eachcell_mean))
    # cor_eachcell_mean_all<-rbind(cor_eachcell_mean_all,i=cor_eachcell_mean)
    cor_each_all2<-rbind(cor_each_all2,cor_eachrun)
  }
  #colnames(mad_eachcell_mean_all)<-colnames(cor_eachcell_mean_all)<-c("alpha","beta","delta","gamma","averaged")
  rownames(mad_each_all2)<-rownames(cor_each_all2)<-c("Wrong-InteRD2",
                                                    "Wrong-TOAST")
}


setwd("/Users/chixiang.chen/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/postdoc/postdoc/deconvolution/ref_based_rd/simulations/simulation_09272021/interd2_new/xinscrna/ref_free/")

meth.names<-c("ref_free")
#comb.names<-c("n25m20","n40m20","n25m40","n40m40","n80m20","n80m40")
#comb.names.true<-c("n25mr20","n40mr20","n25mr40","n40mr40","n80mr20","n80mr40")
# comb.names<-c("n25m20")
# comb.names.true<-c("n25mr20")
all_mad_list<-list()
all_cor_list<-list()
for (i in 1:length(comb.names))
{
  comb.names.used<-comb.names[i]
  comb.names.used.true<-comb.names.true[i]
  mad_eachcell_mean_all<-rep()
  cor_eachcell_mean_all<-rep()
  mad_each_all3<-rep()
  cor_each_all3<-rep()
  for(j in meth.names)
  {
    names.used<-paste0(".*",j,comb.names.used,".*\\.csv$")
    temp = list.files(pattern=names.used,full.names=TRUE)
    myfiles = lapply(temp, read.csv)
    data_est<-do.call(cbind,myfiles)
    data_est<-data_est[,-which(colnames(data_est)=="X")]
    names.used.true<-paste0(".*truep",comb.names.used.true,".*\\.csv$")
    temp = list.files(pattern=names.used.true,full.names=TRUE)
    myfiles = lapply(temp, read.csv)
    #data_truep<-do.call(cbind,myfiles)
    #data_truep<-data_truep[,-which(colnames(data_truep)=="X")]
    
    #mad
    mad_eachrun<-apply(abs(data_est-data_truep),2,mean)
    names(mad_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=100)
    # mad_eachcell_mean<-tapply(mad_eachrun,names(mad_eachrun),mean)
    # mad_eachcell_summary<-tapply(mad_eachrun,names(mad_eachrun),summary)
    # mad_eachcell_mean<-c(mad_eachcell_mean,mean(mad_eachcell_mean))
    # mad_eachcell_mean_all<-rbind(mad_eachcell_mean_all,i=mad_eachcell_mean)
    mad_each_all3<-rbind(mad_each_all3,mad_eachrun)
    
    cor_eachrun<-unlist(sapply(1:ncol(data_est),function (x)
    {
      cor(data_truep[,x],data_est[,x],method = c("kendall"))
    }
    ))
    names(cor_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=100)
    # cor_eachcell_mean<-tapply(cor_eachrun,names(cor_eachrun),mean)
    # cor_eachcell_summary<-tapply(cor_eachrun,names(cor_eachrun),summary)
    # cor_eachcell_mean<-c(cor_eachcell_mean,mean(cor_eachcell_mean))
    # cor_eachcell_mean_all<-rbind(cor_eachcell_mean_all,i=cor_eachcell_mean)
    cor_each_all3<-rbind(cor_each_all3,cor_eachrun)
  }
  #colnames(mad_eachcell_mean_all)<-colnames(cor_eachcell_mean_all)<-c("alpha","beta","delta","gamma","averaged")
  rownames(mad_each_all3)<-rownames(cor_each_all3)<-c("Ref-free")
}
#names(all_mad_list)<-names(all_cor_list)<-comb.names
mad_each_all_pool<-rbind(mad_each_all,mad_each_all2,mad_each_all3)

data.mad.plot<-data.frame(MAD=c(mad_each_all_pool))
data.mad.plot$Celltype<-rep(rep(c("Alpha","Beta","Delta","Gamma"),times=100),each=9)
data.mad.plot$method<-rep(rownames(mad_each_all_pool),times=ncol(mad_each_all_pool))

ggplot(data.mad.plot, aes(x=Celltype, y=MAD,fill= method)) + 
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian() +
  xlab("Cell Type")

cor_each_all_pool<-rbind(cor_each_all,cor_each_all2,cor_each_all3)
data.cor.plot<-data.frame(Kendall=c(cor_each_all_pool))
data.cor.plot$Celltype<-rep(rep(c("Alpha","Beta","Delta","Gamma"),times=100),each=9)
data.cor.plot$method<-rep(rownames(cor_each_all_pool),times=ncol(cor_each_all_pool))

ggplot(data.cor.plot, aes(x=Celltype, y=Kendall,fill= method)) + 
  geom_boxplot(outlier.size=0.5) +
  coord_cartesian()+
  xlab("Cell Type")

#calculate p values for mad
wilcox.test(MAD ~ method, data = data.mad.plot[(data.mad.plot$Celltype=="Alpha" & 
                                                  (data.mad.plot$method %in% c("InteRD1","SCDC-ENSEMBLE"))),],exact = FALSE)

#calculate p values for mad
wilcox.test(Kendall ~ method, data = data.cor.plot[(data.cor.plot$Celltype=="Gamma" & 
                                                  (data.cor.plot$method %in% c("InteRD2","Wrong-TOAST"))),],exact = FALSE)


  