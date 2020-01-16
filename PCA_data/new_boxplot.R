# code to make PDFs
source('../scripts/myFig.R')

# shared stuff
# fancy y-axis lab!
lab_rmsd <- expression( bold( RMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )

# lazy shortcut
read_mat <- function(name) {
  # add extension
  name <- paste0(name, '.txt')
  # read dataframe, convert to numerical matrix right away so column names can be NA later
  as.matrix( read.table( name ) )
}

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
  # close PDF
  invisible( dev.off() )
}

# move to data location
setwd("PCA_data")

#####
#Boxplot of RMSD and AUC when N=100 &k=10 & no family strcuture 
####
#n = 100 and m = 10
# read data
rmsd <- read_mat( "rmsd_k_fixed_k_10_pcs_1_90_n_1000_1_2_repeat_10" )
rmsd_lm <- read_mat( "rmsd_k_fixed_k_10_pcs_1_90_n_1000_1_2_repeat_10_LM" )

rmsd_lm <-as.data.frame(rmsd_lm[,1])
rmsd_lmm <- read_mat( "rmsd_k_fixed_k_10_pcs_1_90_n_1000_1_2_repeat_10_GCTAPC" )
auc <- read_mat( "auc_k_fixed_k_10_pcs_1_90_n_1000_1_2_repeat_10" )
auc_lm <- read_mat( "auc_k_fixed_k_10_pcs_1_90_n_1000_1_2_repeat_10_LM" )
auc_lm <-as.data.frame(auc_lm[,1])
auc_lmm <- read_mat( "auc_k_fixed_k_10_pcs_1_90_n_1000_1_2_repeat_10_GCTAPC" )
# make plot!



