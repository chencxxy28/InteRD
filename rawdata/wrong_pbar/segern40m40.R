library(SCDC)
library(xbioc)
#library(pheatmap)
library(limSolve)
library(dplyr)
#library(Seurat)
library(MultiRD)
#library(magic)
#library(psych) 
#library(ggplot2)
#library(reshape2)

#read single cell files
qc.seger <- readRDS("/storage/home/bzs438/chixiang/deconvolution/singlecell/qc_segerstolpe.rds") #after qc
qc.baron <- readRDS("/storage/home/bzs438/chixiang/deconvolution/singlecell/qc_baron.rds") #after qc
#qc.xin.all<-readRDS("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/postdoc/postdoc/deconvolution/real data scdc/qc.xin.all.rds") #after qc
#qc.xin <- readRDS("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/postdoc/postdoc/deconvolution/real data scdc/qc_xin.rds") #after qc
xin<-readRDS("/storage/home/bzs438/chixiang/deconvolution/singlecell/Xin_nonD.rds")

# qc.seger <- readRDS("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/postdoc/postdoc/deconvolution/real data scdc/qc_segerstolpe.rds") #after qc
# qc.baron <- readRDS("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/postdoc/postdoc/deconvolution/real data scdc/qc_baron.rds") #after qc
# xin<-readRDS("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/postdoc/postdoc/deconvolution/real data scdc/Xin_nonD.rds")

#change here
list_marker<-readRDS("/storage/home/bzs438/chixiang/deconvolution/markergenes/list_markerxin40.rds")
#list_marker<-readRDS("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/postdoc/postdoc/deconvolution/ref_based_rd/simulations/simulation_08252021/interb1/list_markerbaron20.rds")

set.seed(3155421)  #change here
#num.markers<-20
samplesize<-40 #change here

####criterion fct
criteria_onegroup<-function (bulk_data, prop_used)
{
  message("calculate criteria")
  bulk_matrix_raw<-(exprs(bulk_data))
  find_zero_index<-which(rowSums(bulk_matrix_raw,na.rm=T)*1e5==0)
  if(length(find_zero_index)>1)
  {
    bulk_nozero<-getCPM0(bulk_matrix_raw[-find_zero_index,])
  }else{
    bulk_nozero<-getCPM0(bulk_matrix_raw)
  }
  all_nozero_index<-which(apply(bulk_nozero,1,function (x) sum(x==0))==0)
  bulk_all_nozero<-bulk_nozero[all_nozero_index,]
  check_rowsum<-rowSums(bulk_all_nozero)
  bulk_all_nozero<-bulk_all_nozero[check_rowsum<quantile(check_rowsum,0.95) & check_rowsum>quantile(check_rowsum,0.15),]
  sigma_subject<-apply(bulk_all_nozero,1,sd)
  bulk_all_nozero<-bulk_all_nozero[sigma_subject!=0,]
  gene_validate<-rownames(bulk_all_nozero)
  
  
  prop_new<-as.matrix(prop_used)
  ##do training
  
  bulk_nozero_train<-bulk_all_nozero
  prop_new_train<-prop_new
  X_all<-t(sapply(1:length(gene_validate),function (xx)
  {
    gene_xx<-gene_validate[xx]
    bulk_xx<-bulk_nozero_train[gene_xx,]*1e5
    prop_xx<-prop_new_train
    fit<-nnls(A=prop_xx,B=bulk_xx)
    fit$X
  }
  ))
  rownames(X_all)<-gene_validate
  
  bulk_nozero_test<-bulk_all_nozero
  prop_new_test<-prop_new
  prop_new_withscaler<-t(sapply(1:ncol(bulk_nozero_test), function (x)
  {
    #each subject estimation
    Y_initial<-bulk_nozero_test[,x]
    zero_index<-which(Y_initial==0) #find genes with zero expression in current subject
    if(length(zero_index)>0)
    {
      Y_initial<-Y_initial[-zero_index] #delate zero expressed genes
    }
    genes_left<-names(Y_initial)
    
    #construct Y 
    gene_used_now<-intersect(genes_left,gene_validate)
    Y<-Y_initial[gene_used_now]*1e5
    #Y<-Y[Y<quantile(Y,0.85) & Y>quantile(Y,0.15)] #remove the outliers
    gene_used_final<-names(Y)
    
    #construct X
    X<-X_all[gene_used_final,]
    
    values_final<-sort(abs(Y-X%*%(prop_new_test[x,])),decreasing = T)[1:50]
    
    c(median(values_final),mean(values_final),sqrt(mean(values_final^2)))
    #print(x)
  })
  )
  apply(prop_new_withscaler,2,mean)
}


#######################################################
####select markder genes by using 
#######################################################

#create seurat object and do marker selection from baron
cell_type_unique<-c("alpha","beta","delta","gamma")
# test.ident<-qc.baron[["sc.eset.qc"]]@phenoData@data[["cluster"]]
# sc.counts<-exprs(qc.baron[["sc.eset.qc"]])[,test.ident %in% cell_type_unique]
# test.cell.names <- colnames(sc.counts)
# test.ident.sub<-test.ident[test.ident %in% cell_type_unique]
# names(test.ident.sub) <- test.cell.names
# seurat_baron_obj <- CreateSeuratObject(
#     counts=sc.counts,
#     project = "CreateSeuratObject",
#     assay = "RNA",
#     names.field = 1,
#     names.delim = "_",
#     meta.data = NULL,
#     min.cells = 50, 
#     min.features = 200
# )
# Idents(seurat_baron_obj)<-test.ident.sub
# 
# #do normalization
# seurat_baron_obj <- NormalizeData(seurat_baron_obj, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# #cluster5.markers <- FindAllMarkers(seurat_baron_obj, ident.1 = beta, ident.2=gamma, min.pct = 0.25)
# marker.selection <- FindAllMarkers(seurat_baron_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# marker.pool<-marker.selection %>%
#     group_by(cluster) %>%
#     top_n(n = num.markers, wt = avg_log2FC)
# 
# list_marker<-lapply(seq_len(length(cell_type_unique)),function (x)
# {
#     marker.pool$gene[marker.pool$cluster== cell_type_unique[x]]
# }
# )

# singleref_matrix<-exprs(qc.seger[["sc.eset.qc"]])
# cell_type_vector<-qc.seger[["sc.eset.qc"]]@phenoData@data[["cluster"]]
# #cell_type_unique<-unique(cell_type_vector)
# cell_type_unique<-c("alpha","beta","delta","gamma")
# mean_profile_ct<-NULL
# for(ct in cell_type_unique)
# {
#     ct_index<-which(cell_type_vector==ct)
#     matrix_ct<-singleref_matrix[,ct_index]
#     vector_ct<-apply(matrix_ct,1,mean)
#     mean_profile_ct<-cbind(mean_profile_ct,vector_ct)
# }
# 
# mean_profile_ct_screen<-mean_profile_ct[apply(mean_profile_ct,1,sum)>=10,]
# ra_mean_profile_ct_screen<-as.matrix(t(apply(mean_profile_ct_screen,1,function (x) x/sum(x))))
# colnames(ra_mean_profile_ct_screen)<-as.vector(cell_type_unique)
# head(ra_mean_profile_ct_screen)
# 
# 
# list_marker<-lapply(1:ncol(ra_mean_profile_ct_screen),function (x)
# {
#     gene_select<-rownames(ra_mean_profile_ct_screen)[which(ra_mean_profile_ct_screen[,x]<=1 & ra_mean_profile_ct_screen[,x]>0.8)]
# }
# )
# 
# max_expr<-apply(ra_mean_profile_ct_screen,1,max)
# markergene_matrix<-ra_mean_profile_ct_screen[max_expr<=1 & max_expr>0.8,]
#heatmap(markergene_matrix)
#markergenes<-rownames(ra_mean_profile_ct_screen)[max_expr>=0.8]

#######################################################
####simulate data
#######################################################


overall<-100
toast<-NULL
interd1<-NULL
interd2<-NULL
ref_free<-NULL
scdc_baron<-NULL
scdc_xin<-NULL
scdc_ensemble_mad_all<-NULL
truep_all<-NULL

for(iter_rep in 1:overall)
{
  pseudo.seger<-generateBulk(qc.seger[["sc.eset.qc"]], ct.varname = "cluster", sample = "sample", ct.sub = c("alpha","beta","delta","gamma"), nbulk = samplesize, low_s = 0.3, upp_s = 0.7)
  colnames(pseudo.seger[["pseudo_eset"]]@phenoData@data)[1]<-"sample"
  names(pseudo.seger)[1]<-"truep"
  truep<-pseudo.seger$truep[complete.cases(pseudo.seger$truep),]
  
  pancreas.sc <- list(baron = qc.baron$sc.eset.qc,
                      xin   = xin
  )
  fit_music_semi<-SCDC_ENSEMBLE(bulk.eset = pseudo.seger$pseudo_eset, sc.eset.list = pancreas.sc, ct.varname = "cluster",
                                sample = "sample", weight.basis = T,truep = truep, ct.sub =  c("alpha","beta","delta","gamma"), search.length = 0.02,grid.search=T)
  comb<-fit_music_semi$prop.only
  weight_matrix<-fit_music_semi$w_table["mAD_Y",1:2]
  SCDC_ENSEMBLE_MAD<-SCDC:::wt_prop(weight_matrix,comb)
  #SCDC_ENSEMBLE_MAD<-read.table("/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/wmarkers/SCDC_ENSEMBLE_MAD.txt")
  
  #calculate avaerage of truep
  
  P1<-comb[[1]][,cell_type_unique]
  P2<-comb[[2]][,cell_type_unique]
  # P1<-read.table("/Users/chixiangchen/Documents/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/cmarkers/baronref.txt") 
  # P1<-P1[,cell_type_unique]
  # P1<-as.matrix(P1[,cell_type_unique])
  # P2<-read.table("/Users/chixiangchen/Documents/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/cmarkers/xinref.txt")
  # P2<-P2[,cell_type_unique]
  # P2<-as.matrix(P2[,cell_type_unique])
  #P3<-read.table("/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_methods/seger_bulk/used_method_ensemble/scdc_xinref.txt")
  
  comb_used<-list(as.matrix(P1),as.matrix(P2))
  
  ref_based_ensemble.onegroup<-function (pseudo.seger,list_marker,cell_type_unique,comb_used,lambda_option){
    bulk.eset<-pseudo.seger
    bulk_data<-(exprs(bulk.eset))
    marker_gene_used<-unlist(list_marker)
    gene_test<-setdiff(rownames(bulk_data),marker_gene_used)
    bulk_data<-bulk_data[gene_test,]
    exprs_data<-as.matrix(bulk_data)
    pdata<-data.frame(sample=colnames(bulk_data))
    fdata<-data.frame(genes=rownames(bulk_data))
    rownames(pdata)<-colnames(bulk_data)
    rownames(fdata)<-rownames(bulk_data)
    pseudo.seger_test<-ExpressionSet(exprs_data,
                                     AnnotatedDataFrame(pdata),
                                     AnnotatedDataFrame(fdata))
    
    #run reference free approach
    ct.sub<-cell_type_unique
    bulk.eset<-pseudo.seger
    bulk_matrix_raw<-(exprs(bulk.eset))
    find_zero_index<-which(rowSums(bulk_matrix_raw,na.rm=T)*1e5==0)
    if(length(find_zero_index)>0)
    {
      bulk_nozero<-getCPM0(bulk_matrix_raw[-find_zero_index,])
    }else{
      bulk_nozero<-getCPM0(bulk_matrix_raw)
    }
    marker_all<-unlist(list_marker)
    #lambda_option<-c(seq(from=0,to=0.075,length=15),1,10,30,50,1000)
    #lambda_option<-c(1,10,30,50,1000)
    
    
    
    criterion<-NULL
    est_all<-list(NULL)
    weights_all<-list(NULL)
    
    marker_list_sub<-lapply(1:length(list_marker),function (xx)
    {
      marker_xx<-intersect(list_marker[[xx]],rownames(bulk_nozero)) #create new marker list
    }
    )
    names(marker_list_sub)<-ct.sub
    gene_used<-unlist(marker_list_sub) #vectorize the values
    cluster_identifier<-unlist(lapply(1:length(marker_list_sub),function (xx)
    {
      rep(xx,length(marker_list_sub[[xx]])) #create the identifier to discriminate which cell type are markers coming from
    }
    ))
    
    
    ll<-1
    for (lambda in lambda_option)
    {
      num_est<-length(comb_used)
      prop_old<-0
      for (i in 1:num_est)
      {
        prop_old<-prop_old+comb_used[[i]]/num_est
      }
      iter<-1
      weights_old<-rep(1/num_est,num_est)
      
      repeat{  
        X<-t(sapply(1:length(gene_used),function (xx)
        {
          x_xx<-rep(0,length(ct.sub))
          gene_xx<-gene_used[xx]
          bulk_xx<-bulk_nozero[gene_xx,]
          Y_aug<-c(bulk_xx,rep(0,length(ct.sub)))
          
          cell_type_index<-cluster_identifier[xx]
          #prop_xx<-prop_old[,cell_type_index]
          penalty_matrix<-diag(1,length(ct.sub))
          penalty_matrix[cell_type_index,cell_type_index]<-0
          #colnames(penalty_matrix)<-ct.sub
          X_aug<-rbind(prop_old,lambda*penalty_matrix)
          fit<-nnls(A=(X_aug),B=Y_aug)
          fit$X*1e5
          #x_xx[cell_type_index]<-sum((bulk_xx*prop_xx))/sum(prop_xx^2)*1e5
          #x_xx[cell_type_index]<-sum(bulk_xx)/sum(prop_xx)*1e5 #this formula is weird, I think the previous one is better
        }
        ))
        
        Y_all<-c(bulk_nozero[gene_used,])*1e5
        
        X_all<-rep()
        for(i in 1:ncol(bulk_nozero))
        {
          Matrix_i<-matrix(0,nrow=length(gene_used),ncol=num_est)
          for (j in 1:num_est)
          {
            Matrix_i[,j]<-X%*%(comb_used[[j]][i,])
          }
          X_all<-rbind(X_all,Matrix_i)
        }
        fit<-nnls(A=X_all,B=Y_all)
        weights_new<-fit$X/sum(fit$X)
        
        prop_new<-0
        for (i in 1:num_est)
        {
          prop_new_i<-weights_new[i]*comb_used[[i]]
          prop_new<-prop_new+prop_new_i
        }
        
        if(sum(abs(weights_new-weights_old))<1e-06)
        {break}else{
          prop_old<-prop_new
          weights_old<-weights_new
          iter<-iter+1
        }
      }
      criterion_sub<-criteria_onegroup(bulk_data=pseudo.seger_test, prop_used=prop_new)[2]
      criterion<-c(criterion,criterion_sub)
      est_all[[ll]]<-prop_new
      weights_all[[ll]]<-weights_new
      #print(ll)
      ll<-ll+1
    }
    
    prop_est<-est_all[[which.min(criterion)]]
    list(est=est_all,musec=criterion,weights_list=weights_all)
  }
  
  lambda_option<-c(0,0.01,0.05,0.1,1,5,100)
  results_ref_based_ensemble<-ref_based_ensemble.onegroup(pseudo.seger=pseudo.seger$pseudo_eset,list_marker,cell_type_unique,comb_used=comb_used,lambda_option)
  criteria_values<-results_ref_based_ensemble$musec
  weight_final<-results_ref_based_ensemble$weights_list[[which.min(criteria_values)]]
  
  comb_sampled<-results_ref_based_ensemble$est[[which.min(criteria_values)]]
  
  sampleid<-qc.seger[["sc.eset.qc"]]@phenoData@data[["sample"]]
  clusterid<-qc.seger[["sc.eset.qc"]]@phenoData@data[["cluster"]]
  sc_proportions<-sapply(1:length(unique(qc.seger[["sc.eset.qc"]]@phenoData@data[["sample"]])),function (x)
  {
    ct_x<-qc.seger[["sc.eset.qc"]]@phenoData@data[["cluster"]][sampleid %in% unique(sampleid)[x] & (clusterid %in% cell_type_unique)]
    table(ct_x)/sum(table(ct_x))
  }
  )
  
  
  # sampleid<-xin@phenoData@data[["sample"]]
  # clusterid<-xin@phenoData@data[["cluster"]]
  # sc_proportions<-sapply(1:length(unique(xin@phenoData@data[["sample"]])),function (x)
  # {
  #     ct_x<-xin@phenoData@data[["cluster"]][sampleid %in% unique(sampleid)[x] & (clusterid %in% cell_type_unique)]
  #     table(ct_x)/sum(table(ct_x))
  # }
  # )
  
  ave_est<-apply(sc_proportions[cell_type_unique,],1,mean)
  ave_sd<-apply(sc_proportions[cell_type_unique,],1,sd)

  ave_est<-apply(P2,2,mean)
  ave_sd<-apply(P2,2,sd)
  
  #ave_est<-apply(truep[,cell_type_unique],2,mean)
  #ave_sd<-1
  
  #lambda_option<-c(seq(from=0,to=0.05,length=10),1,10,50,100,1000)
  #lambda_option<-c(0.000000e+00, 5.357143e-03, 0.03,10,50,100,500,1000)
  #lambda_option<-c(0,0.005,0.01,0.03,1,5,10,100) #choose 600 is the best maybe
  #lambda_option<-c(10000)
  #lambda_option<-c(seq(from=0,to=0.5,length=10))
  #lambda_option<-c(0.000000e+00, 5.357143e-03, 0.03,10,50,100,500,1000)
  #lambda_option<-c(0,0.005,0.01,0.02,0.03) #choose 600 is the best maybe
  lambda_option<-c(0,seq(from=1,to=20,length=4),seq(from=30,to=100,length=4),200,500,1000000^2)
  #lambda_option<-c(0,1,50,100,200,500,100,1000000^2)
  MURD.onegroup<-function (pseudo.seger,list_marker,cell_type_unique,comb_sampled,ave_est,ave_sd,lambda_option1){
    bulk.eset<-pseudo.seger
    bulk_data<-(exprs(bulk.eset))
    marker_gene_used<-unlist(list_marker)
    gene_test<-setdiff(rownames(bulk_data),marker_gene_used)
    bulk_data<-bulk_data[gene_test,]
    exprs_data<-as.matrix(bulk_data)
    pdata<-data.frame(sample=colnames(bulk_data))
    fdata<-data.frame(genes=rownames(bulk_data))
    rownames(pdata)<-colnames(bulk_data)
    rownames(fdata)<-rownames(bulk_data)
    pseudo.seger_test<-ExpressionSet(exprs_data,
                                     AnnotatedDataFrame(pdata),
                                     AnnotatedDataFrame(fdata))
    
    #run reference free approach
    comb_sampled<-comb_sampled[,cell_type_unique]
    ct.sub<-cell_type_unique
    bulk.eset<-pseudo.seger
    bulk_matrix_raw<-(exprs(bulk.eset))
    find_zero_index<-which(rowSums(bulk_matrix_raw,na.rm=T)*1e5==0)
    if(length(find_zero_index)>0)
    {
      bulk_nozero<-getCPM0(bulk_matrix_raw[-find_zero_index,])
    }else{
      bulk_nozero<-getCPM0(bulk_matrix_raw)
    }
    marker_genes<-unlist(list_marker) #match the marker gene 
    #check marker gene availability for each cell types
    marker_list_sub_i<-lapply(1:length(list_marker),function (xx)
    {
      marker_xx<-intersect(list_marker[[xx]],rownames(bulk_nozero)) #create new marker list
    }
    )
    names(marker_list_sub_i)<-ct.sub
    gene_used<-unlist(marker_list_sub_i) #vectorize the values
    cluster_identifier<-unlist(lapply(1:length(marker_list_sub_i),function (xx)
    {
      rep(xx,length(marker_list_sub_i[[xx]])) #create the identifier to discriminate which cell type are markers coming from
    }
    ))
    #setdiff(unlist(list_marker),intersect(rownames(bulk_nozero),unlist(list_marker)))
    marker_all<-unlist(list_marker)
    #lambda_option<-c(seq(from=0,to=0.075,length=15),1,10,30,50,1000)
    #lambda_option<-c(1)
    group<-pseudo.seger@phenoData@data[["groups"]]
    criterion<-NULL
    est_all<-rep()
    for(lambda1 in lambda_option1)
    {
        iter<-0
        prop_old<-matrix(1/length(ct.sub),nrow=ncol(bulk_matrix_raw),ncol=length(ct.sub))
        X<-t(sapply(1:length(gene_used),function (xx)
        {
          x_xx<-rep(0,length(ct.sub))
          gene_xx<-gene_used[xx]
          bulk_xx<-bulk_nozero[gene_xx,]
          cell_type_index<-cluster_identifier[xx]
          prop_xx<-prop_old[,cell_type_index]
          #x_xx[cell_type_index]<-sum((bulk_xx*prop_xx))/sum(prop_xx^2)*1e5
          x_xx[cell_type_index]<-sum(bulk_xx)/sum(prop_xx)*1e5 #this formula is weird, I think the previous one is better
          x_xx
        }
        ))
        rownames(X)<-gene_used
        lambda_adjust<-mean((bulk_nozero[gene_used,]*1e5-X%*%t(prop_old))^2)
        repeat{
          X<-t(sapply(1:length(gene_used),function (xx)
          {
            x_xx<-rep(0,length(ct.sub))
            gene_xx<-gene_used[xx]
            bulk_xx<-bulk_nozero[gene_xx,]
            cell_type_index<-cluster_identifier[xx]
            prop_xx<-prop_old[,cell_type_index]
            #x_xx[cell_type_index]<-sum((bulk_xx*prop_xx))/sum(prop_xx^2)*1e5
            x_xx[cell_type_index]<-sum(bulk_xx)/sum(prop_xx)*1e5 #this formula is weird, I think the previous one is better
            x_xx
          }
          ))
          rownames(X)<-gene_used
          
          prop_new<-t(sapply(1:ncol(bulk_nozero), function (x)
          {
            #each subject estimation
            Y_initial<-bulk_nozero[,x]
            zero_index<-which(Y_initial==0) #find genes with zero expression in current subject
            if (length(zero_index)>0)
            {
              Y_initial<-Y_initial[-zero_index] #delate zero expressed genes !!!!!!!!!!!!check this!!!!!!!!!!!
            }
            #Y_initial<-Y_initial[Y_initial<quantile(Y_initial,0.99)] #remove the outliers
            gene_names<-names(Y_initial)
            marker_genes<-intersect(gene_names,marker_all) #match the marker gene
            #check marker gene availability for each cell types
            marker_list_sub_i<-lapply(1:length(list_marker),function (xx)
            {
              marker_xx<-intersect(list_marker[[xx]],marker_genes) #create new marker list
            }
            )
            names(marker_list_sub_i)<-ct.sub
            gene_used<-unlist(marker_list_sub_i) #vectorize the values
            cluster_identifier<-unlist(lapply(1:length(marker_list_sub_i),function (xx)
            {
              rep(xx,length(marker_list_sub_i[[xx]])) #create the identifier to discriminate which cell type are markers coming from
            }
            ))
            find_zero_counts<-which(table(cluster_identifier)==0)
            
            #construct Y and corresponding augmented iterms (major change)
            Y<-Y_initial[gene_used]*1e5
            lambda_adjust1<-lambda_adjust*lambda1
            
            ave_truep<-as.matrix(comb_sampled)[x,]*sqrt(lambda_adjust1)
            ave_truep2<-ave_est*sqrt(lambda_adjust)/ave_sd
            
            X<-X[gene_used,]
            X_aug<-sqrt(lambda_adjust1)*diag(1,length(ct.sub))
            X_aug2<-sqrt(lambda_adjust)*diag(1,length(ct.sub))/ave_sd
            
            Y_aug<-ave_truep
            y_aug2<-ave_truep2
            #Y_aug<-c(sapply(1:length(combo), function (xx) {sqrt(lambda_adjust*weight[xx])*(combo[[xx]][x,])}))
            Y_comb<-c(Y,Y_aug,y_aug2)
            
            #X_aug<-do.call(rbind,lapply(1:length(combo), function (xx) {sqrt(lambda_adjust*weight[xx])*diag(1,length(ct.sub))}))
            X_comb<-rbind(X,X_aug,X_aug2)
            #heatmap(X)
            ##construct the constrain matrix
            E_used<-rep(1,length(ct.sub))
            F_used<-1
            G_used<-diag(1,length(ct.sub))
            H_used<-rep(0,length(ct.sub))
            #fit<-lsei(A=X,B=Y,E=E_used,F=F_used,G=G_used,H=H_used)
            #fit$X
            fit<-nnls(A=X_comb,B=Y_comb)
            
            #calcualte GCV values
            #x_x<-t(X)%*%X
            #inverse_term<-solve(x_x+sum(lambda_adjust*weight))
            #P_lambda<-X%*%inverse_term%*%t(X)
            #penalty_1<-(1-(tr(P_lambda)-1)/length(Y))^2
            #penalty_2<-
            #gcv_values<-(Y%*%((diag(1,nrow(P_lambda))-P_lambda)^2)%*%Y)/penalty_1+
            fit$X/sum(fit$X)
            #print(x)
          })
          )
          if(mean(abs(prop_new[,1:length(ct.sub)]-prop_old))<0.0005 | iter>1000)
          {
            break
          }else{
            prop_old<-prop_new[,1:length(ct.sub)]
            lambda_adjust<-mean((bulk_nozero[gene_used,]*1e5-X%*%t(prop_old))^2)
            iter<-iter+1
          }
        }
        est_meta_all<-list(prop_new)
        criterion_sub<-criteria_onegroup(bulk_data=pseudo.seger_test, prop_used=prop_new)[2]
        criterion<-c(criterion,criterion_sub)
        est_all<-append(est_all,est_meta_all)
        #print(ll)
    }
    
    prop_est<-est_all[[which.min(criterion)]]
    list(est=est_all,musec=criterion)
  }
  
  
  results_1<-MURD.onegroup(pseudo.seger=pseudo.seger$pseudo_eset,list_marker,cell_type_unique,comb_sampled=comb_sampled,ave_est,ave_sd,lambda_option1=lambda_option)
  musec_1<-results_1$musec
  prop_est<-results_1$est[[which.min(musec_1)]]
  prop_est_ref<-results_1$est[[1]]
  
  toast<-cbind(toast,prop_est_ref)
  interd2<-cbind(interd2,prop_est)
  interd1<-cbind(interd1,comb_sampled)
  scdc_baron<-cbind(scdc_baron,comb_used[[1]])
  scdc_xin<-cbind(scdc_xin,comb_used[[2]])
  scdc_ensemble_mad_all<-cbind(scdc_ensemble_mad_all,SCDC_ENSEMBLE_MAD)
  truep_all<-cbind(truep_all,truep)
  print(iter_rep)
}

#change here
write.csv(toast,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/toastn40m40.csv")
write.csv(interd2,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/interd2n40m40.csv")
write.csv(interd1,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/interd1n40m40.csv")
write.csv(scdc_baron,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/scdc_baronn40m40.csv")
write.csv(scdc_xin,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/scdc_xinn40m40.csv")
write.csv(scdc_ensemble_mad_all,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/scdc_ensemble_mad_alln40m40.csv")
write.csv(truep_all,file="/storage/home/bzs438/chixiang/deconvolution/interd2_new/xinscrna/wrong_pbar/truepn40mr40.csv")

# prop_est<-prop_est_ref
# prop_est<-comb_sampled
# mean(apply(abs(truep-prop_est),2,mean))
# # mean(apply(abs(truep-prop_est),2,mean)/apply(truep,2,mean))
# mean(sapply(1:length(cell_type_unique),function (x)
# {
#   cor(truep[,x],prop_est[,x],method = c("kendall"))
# }
# ))
# write.table(prop_est,file="/Users/chixiangchen/Documents/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/cmarkers/InteRB_zscore.txt",quote = FALSE,sep="\t")
# 
# 
# write.table(P1,file="/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/wmarkers/baronref.txt",quote = FALSE,sep="\t")
# write.table(P2,file="/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/wmarkers/xinref.txt",quote = FALSE,sep="\t")
# write.table(prop_est,file="/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/wmarkers/ours_ciborsortx.txt",quote = FALSE,sep="\t")
# write.table(SCDC_ENSEMBLE_MAD,file="/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/wmarkers/SCDC_ENSEMBLE_MAD.txt",quote = FALSE,sep="\t")
# write.table(truep,file="/Users/chixiangchen/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/wmarkers/truep.txt",quote = FALSE,sep="\t")
# write.table(prop_est,file="/Users/chixiangchen/Documents/Dropbox/My Mac (Chixiang’s MacBook Pro (2))/Desktop/postdoc/deconvolution/ref_based_rd/simulations/ensemble_ref/seger_bulk/two_ref/cmarkers/InteRB_zscore.txt",quote = FALSE,sep="\t")

# data_est<-interd2
# data_truep<-truep_all
# mad_eachrun<-apply(abs(data_est-data_truep),2,mean)
# names(mad_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=5)
# mad_eachcell_mean<-tapply(mad_eachrun,names(mad_eachrun),mean)
# mad_eachcell_summary<-tapply(mad_eachrun,names(mad_eachrun),summary)
# mad_eachcell_mean
# mean(mad_eachcell_mean)
# mad_eachcell_summary
# mad_eachcell_sd

#pearson correlation
# cor_eachrun<-unlist(sapply(1:ncol(data_est),function (x)
# {
#   CCC_all<-CCC(data_est[,x],data_truep[,x])
#   CCC_all$rho[1]
# }
#   ))

# cor_eachrun<-unlist(sapply(1:ncol(data_est),function (x)
# {
#   cor(data_truep[,x],data_est[,x],method = c("kendall"))
# }
# ))
# names(cor_eachrun)<-rep(c("alpha","beta","delta","gamma"),times=5)
# cor_eachcell_mean<-tapply(cor_eachrun,names(cor_eachrun),mean)
# cor_eachcell_summary<-tapply(cor_eachrun,names(cor_eachrun),summary)
# cor_eachcell_mean
# mean(cor_eachcell_mean)









