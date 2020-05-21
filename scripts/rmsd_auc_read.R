rmsd_read<-function(pc,method,n){
  rmsd<-c()
  for (i in 1:(pc+1)){
    
    name<-paste0("rmsd_",method,"_n_",n,"_pcs_",i-1,".txt")
    rmsd_new<-read.table(name)
    auc_new<-unlist(rmsd_new)
    rmsd<-c(rmsd,rmsd_new)
  }
  write.table(unlist(rmsd), file = paste0("rmsd_pca_n_", n, "_pcs_",pc,"_collection.txt") )
}

auc_read<-function(pc,method,n){
  auc<-c()
  for (i in 1:(pc+1)){
    
    name<-paste0("auc_",method,"_n_",n,"_pcs_",i-1,".txt")
    auc_new<-read.table(name)
    auc_new<-unlist(auc_new)
    auc<-c(auc,auc_new)
  }
  write.table(unlist(auc), file = paste0("auc_pca_n_", n, "_pcs_",pc,"_collection.txt") )
}

#read rmsd of pca for pcs from 1 to 90 when sample szie is 100
#example: rmsd_read(pc=90, method="pca", n=100)