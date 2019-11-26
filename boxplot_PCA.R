#####
#Boxplot of RMSD and AUC when N=100 &k=10 & no family strcuture 
####
setwd("C:/Users/yiqiy/Desktop/PCA_data")
par(mfrow=c(2,1))

rmsd_100<-read.table("rmsd_k_fixed_k_10_pcs_1_90_n_100_10_30_repeat_10.txt")
M_rmsd_LMM_n_100<-read.table("M_rmsd_LMM_n_100.txt")
M_rmsd_LM_n_100<-read.table("M_rmsd_LM_n_100.txt")
rmsd_100<-rmsd_100[,1:90]
rmsd_100<-cbind(M_rmsd_LM_n_100,rmsd_100)

colnames(rmsd_100)[1]<-"LM"
colnames(rmsd_100)[2:91]<-1:90
rmsd_100<-cbind(rmsd_100,M_rmsd_LMM_n_100)
colnames(rmsd_100)[92]<-"LMM"
boxplot(rmsd_100,main="Boxplot of RMSD when n=100", xlab="The number of PCs", ylab="RMSD")



auc_100<-read.table("auc_k_fixed_k_10_pcs_1_90_n_100_10_30_repeat_10.txt")
M_auc_LMM_n_100<-read.table("M_auc_LMM_n_100.txt")
M_auc_LM_n_100<-read.table("M_auc_LM_n_100.txt")
auc_100<-auc_100[,1:90]
auc_100<-cbind(M_auc_LM_n_100,auc_100)

colnames(auc_100)[1]<-"LM"
colnames(auc_100)[2:91]<-1:90
auc_100<-cbind(auc_100,M_auc_LMM_n_100)
colnames(auc_100)[92]<-"LMM"
boxplot(auc_100,main="Boxplot of AUC when n=100", xlab="The number of PCs", ylab="AUC")

#####
#Boxplot of RMSD and AUC when N=1000 &k=10 & no family strcuture 
####
par(mfrow=c(2,1))
rmsd_1000<-read.table("rmsd_k_fixed_k_10_pcs_1_90_n_1000_10_30_repeat_10.txt")
M_rmsd_LMM_n_1000<-read.table("M_rmsd_LMM_n_1000.txt")
M_rmsd_LM_n_1000<-read.table("M_rmsd_LM_n_1000.txt")
rmsd_1000<-rmsd_1000[,1:90]
rmsd_1000<-cbind(M_rmsd_LM_n_100,rmsd_1000)

colnames(rmsd_1000)[1]<-"LM"
colnames(rmsd_1000)[2:91]<-1:90
rmsd_1000<-cbind(rmsd_1000,M_rmsd_LMM_n_1000)
colnames(rmsd_1000)[92]<-"LMM"
boxplot(rmsd_1000,main="Boxplot of RMSD when n=1000", xlab="The number of PCs", ylab="RMSD")



auc_1000<-read.table("auc_k_fixed_k_10_pcs_1_90_n_1000_10_30_repeat_10.txt")
M_auc_LMM_n_1000<-read.table("M_auc_LMM_n_1000.txt")
M_auc_LM_n_1000<-read.table("M_auc_LM_n_1000.txt")
auc_1000<-auc_1000[,1:90]
auc_1000<-cbind(M_auc_LM_n_100,auc_1000)

colnames(auc_1000)[1]<-"LM"
colnames(auc_1000)[2:91]<-1:90
auc_1000<-cbind(auc_1000,M_auc_LMM_n_1000)
colnames(auc_1000)[92]<-"LMM"
boxplot(auc_1000,main="Boxplot of AUC when n=1000", xlab="The number of PCs", ylab="AUC")




#####
#Boxplot of RMSD and AUC when N=1000 &k=10 & family exists
####
par(mfrow=c(2,1))
M_rmsd_LM_n_1000_family<-read.table("M_rmsd_LM_n_1000_family.txt")
M_auc_LM_n_1000_family<-read.table("M_auc_LM_n_1000_family.txt")
M_rmsd_LMM_n_1000_family<-read.table("M_rmsd_LMM_n_1000_family.txt")
M_auc_LMM_n_1000_family<-read.table("M_auc_LMM_n_1000_family.txt")
rmsd_1000_family<-read.table("rmsd_k_fixed_k_10_pcs_2_90_n_1000_10_31_family.txt")
rmsd_miss<-read.table("M_rmsd_LM_n_1000_family_miss.txt")
auc_miss<-read.table("M_auc_LM_n_1000_family_miss.txt")
rmsd_1000_family<-rmsd_1000_family[,1:90]

rmsd_1000_family<-cbind(M_rmsd_LM_n_1000_family,rmsd_1000_family)

colnames(rmsd_1000_family)[1]<-"LM"
colnames(rmsd_1000_family)[2:91]<-1:90
rmsd_1000_family<-cbind(rmsd_1000_family,M_rmsd_LMM_n_1000_family)
colnames(rmsd_1000_family)[92]<-"LMM"

boxplot(rmsd_1000_family,main="Boxplot of RMSD when n=1000 & family structure exists", xlab="The number of PCs", ylab="RMSD")



auc_1000_family<-read.table("auc_k_fixed_k_10_pcs_2_90_n_1000_10_31_family.txt")
auc_1000_family<-auc_1000_family[,1:90]
auc_1000_family<-cbind(M_auc_LM_n_1000_family,auc_1000_family)

colnames(auc_1000_family)[1]<-"LM"
colnames(auc_1000_family)[2:91]<-1:90
auc_1000_family<-cbind(auc_1000_family,M_auc_LMM_n_1000_family)
colnames(auc_1000_family)[92]<-"LMM"
boxplot(auc_1000_family,main="Boxplot of AUC when n=1000 & family structure exists", xlab="The number of PCs", ylab="AUC")
