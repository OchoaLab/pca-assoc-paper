
setwd("C:/Users/yiqiy/Desktop/gas-rgls/scripts")
source("myFig.R")

setwd("C:/Users/yiqiy/Desktop/gas-rgls/PCA_gwas/PCA_data/PCA_GCTApcs")
#read data from different repliction 
read_data<-function(name,rep){
  data<-c()
  for (i in 1:rep){
    name_read<-paste0(name,"_",i,".","txt")
    data_temp<-read.table(name_read)
    data<-rbind(data,data_temp)
  }
  return(data)
}

read_data_family<-function(name,rep){
  data<-c()
  for (i in 1:rep){
    name_read<-paste0(name,"_",i,"_family",".","txt")
    data_temp<-read.table(name_read)
    data<-rbind(data,data_temp)
  }
  return(data)
}

rmsd_pca_100<-read_data('rmsd_k_fixed_k_10_pcs_1_90_n_100_1_28_repeat_10',5)
auc_pca_100<-read_data('auc_k_fixed_k_10_pcs_1_90_n_100_1_28_repeat_10',5)

rmsd_pca_1000<-read_data('rmsd_k_fixed_k_10_pcs_1_90_n_1000_1_28_repeat_10',5)
auc_pca_1000<-read_data('auc_k_fixed_k_10_pcs_1_90_n_1000_1_28_repeat_10',5)

rmsd_gcta_1000<-read_data("rmsd_k_fixed_k_10_pcs_1_90_n_1000_1_28_repeat_10_GCTAPC",5)
auc_gcta_1000<-read_data("auc_k_fixed_k_10_pcs_1_90_n_1000_1_28_repeat_10_GCTAPC",5)

rmsd_gcta_100<-read_data("rmsd_k_fixed_k_10_pcs_1_90_n_100_1_28_repeat_10_GCTAPC",5)
auc_gcta_100<-read_data("auc_k_fixed_k_10_pcs_1_90_n_100_1_28_repeat_10_GCTAPC",5)

rmsd_lm_1000<-read_data("auc_k_fixed_k_10_pcs_1_90_n_1000_1_28_repeat_10_LM",5)
auc_lm_1000<-read_data("auc_k_fixed_k_10_pcs_1_90_n_1000_1_28_repeat_10_lm",5)

rmsd_lm_100<-read_data("auc_k_fixed_k_10_pcs_1_90_n_100_1_28_repeat_10_LM",5)
auc_lm_100<-read_data("auc_k_fixed_k_10_pcs_1_90_n_100_1_28_repeat_10_lm",5)


rmsd_pca_1000_family<-read_data_family('rmsd_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10',5)
auc_pca_1000_family<-read_data_family('auc_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10',5)

rmsd_gcta_1000_family<-read_data_family("rmsd_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10_GCTAPC",5)
auc_gcta_1000_family<-read_data_family("auc_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10_GCTAPC",5)

rmsd_lm_1000_family<-read_data_family("rmsd_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10_LM",5)
auc_lm_1000_family<-read_data_family("auc_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10_LM",5)



lab_rmsd <- expression( bold( RMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )

auc_k_fixed_k_10_pcs_1_90_n_1000_2_22_repeat_10_3_family

# main plotting code!
boxplots_rmsd_auc <- function(name_out, rmsd, rmsd_lm, rmsd_lmm, auc, auc_lm, auc_lmm, r_max = 90) {
  # inspection shows that 91-100 are blank, remove those here
  rs <- 1 : r_max
  rmsd <- rmsd[, rs]
  auc <- auc[, rs]
  rmsd_lmm <- rmsd_lmm[, rs]
  auc_lmm <- auc[, rs]
  # blank all non-multiples of 10, for plot
  rs[ rs %% 10 != 0 ] <- NA
  # add nice names prior to combining
  colnames( rmsd_lmm ) <- rs
  colnames( auc_lmm ) <- rs
  colnames( rmsd_lm ) <- 0
  colnames( auc_lm ) <- 0
  colnames( rmsd ) <- rs
  colnames( auc ) <- rs
  # combine
  rmsd <- cbind( rmsd_lmm, rmsd_lm, rmsd )
  auc <- cbind( auc_lmm, auc_lm, auc )
  
  # start PDF
  name_out <- paste0('../', name_out)
  myFig(name_out, width = 7, height = 5, botMar = 1.5)
  # add lower margin, so inner margins can be even smaller
  par( oma = c(1.5, 0, 0, 0) )
  # two panels
  par( mfrow = c(2, 1) )
  ## # width hack
  ## widths <- rep.int( 1, ncol(auc) ) # all are the same
  ## #widths[1] <- 8 # but LMM (first) is this much larger!
  at <- 1:ncol(auc)
  at[1] <- -5
  # boxplots!
  boxplot(rmsd, main = "", xlab = "", ylab = lab_rmsd, at = at) # , width = widths
  boxplot(auc, main = "", xlab = '', ylab = lab_auc, at = at) # , width = widths
  # add outer margin
  mtext(
    "Number of PCs (r)",
    side = 1,
    line = 0.5,
    adj = 0.55,
    outer = TRUE
  )
  invisible( dev.off() )
}


boxplots_rmsd_auc(name_out="boxplot_n_100", rmsd=rmsd_pca_100, rmsd_lm=rmsd_lm_100, rmsd_lmm=rmsd_gcta_100, auc=auc_pca_100, auc_lm
                  =auc_lm_100, auc_lmm=auc_gcta_100, r_max = 90)
boxplots_rmsd_auc(name_out="boxplot_n_1000", rmsd=rmsd_pca_1000, rmsd_lm=rmsd_lm_1000, rmsd_lmm=rmsd_gcta_1000, auc=auc_pca_1000, auc_lm
                  =auc_lm_1000, auc_lmm=auc_gcta_1000, r_max = 90)
boxplots_rmsd_auc(name_out="boxplot_n_1000_family", rmsd=rmsd_pca_1000_family, rmsd_lm=rmsd_lm_1000_family, rmsd_lmm=rmsd_gcta_1000_family, auc=auc_pca_1000_family, auc_lm
                  =auc_lm_1000_family, auc_lmm=auc_gcta_1000_family, r_max = 90)