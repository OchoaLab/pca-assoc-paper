library(scales) # for transparency

source('../../../scripts/myFig.R')

#read data from different repliction 
read_data<-function(name, rep, family = FALSE, verbose = TRUE){
  if (verbose)
    message('loading repeats for: ', name)
  # figure out file suffixes
  suffix <- '.txt'
  if (family) {
    if (verbose)
      message('(FAMILY)')
    suffix <- paste0( '_family', suffix )
  }
  # start loading files
  data <- c()
  for (i in 1:rep){
    file_in <- paste0(name, "_", i, suffix)
    # don't just die when things are missing...
    if ( file.exists( file_in ) ) {
      data_temp <- read.table( file_in )
    } else {
      warning('expected file missing: ', file_in)
      # HACK add missing rows for one particular missing file (this won't work in other situations...)
      data_temp <- matrix(NA, nrow = 2, ncol = 90)
    }
    data <- rbind(data, data_temp)
  }
  if (verbose)
    message("\tReplicates: ", nrow(data))
  return(data)
}



# large sample size
rmsd_pca_1000 <- read_data('rmsd_pca_n_1000_3_19', 5)
auc_pca_1000 <- read_data('auc_pca_n_1000_3_19', 5)

rmsd_gcta_1000 <- read_data('rmsd_gcta_n_1000_3_19', 5)
auc_gcta_1000 <-  read_data('auc_gcta_n_1000_3_19', 5)


# small sample size
rmsd_pca_100 <- read_data('rmsd_pca_n_100_3_19', 4)
auc_pca_100 <- read_data('auc_pca_n_100_3_19', 4)

rmsd_gcta_100 <- read_data('rmsd_gcta_n_100_3_19', 4)
auc_gcta_100 <-  read_data('auc_gcta_n_100_3_19', 4)

# large sample size & family structure

rmsd_pca_1000_family <- read_data('rmsd_pca_n_1000_3_19_family', 5)
auc_pca_1000_family <- read_data('auc_pca_n_1000_3_19_family', 5)

rmsd_gcta_1000_family <- read_data('rmsd_gcta_n_1000_3_19_family', 5)
auc_gcta_1000_family <-  read_data('auc_gcta_n_1000_3_19_family', 5)


# plotting labels
lab_rmsd <- expression( bold( RMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )

# remove data beyond r_max, but first checking that it is truly zero
trim_matrix <- function( data, r_max ) {
  nr <- ncol(data)
  if (nr > r_max) {
    # trimming is needed in this case
    # indexes that should be zero
    rs <- ( r_max + 1 ) : nr
    stopifnot( all( data[, rs] == 0 ) )
    # now we trim as needed
    rs <- 1 : r_max
    data <- data[, rs]
  } else if (nr < r_max)
    stop('data did not have enough columns: wanted ', r_max, ', had ', nr)
  # else there's no changes and no errors
  # in all cases return data
  return( data )
}

# main plotting code!
boxplots_rmsd_auc <- function(
  name_out,
  rmsd_pca,
  rmsd_lmm,
  auc_pca,
  auc_lmm,
  r_max = 90,
  legend_pos = 'topright'
) {
  # inspection shows that 91-100 are blank, remove those here
  rmsd_pca <- trim_matrix( rmsd_pca, r_max=91 )
  auc_pca <- trim_matrix( auc_pca, r_max=91 )
  rmsd_lmm <- trim_matrix( rmsd_lmm, r_max=91 )
  auc_lmm <- trim_matrix( auc_lmm, r_max=91 )
  # concatenate the fixed effects series
  rmsd_fixed <-  rmsd_pca
  auc_fixed <-  auc_pca 
  # same for mixed effects, though the LMM with zero PCs is currently missing
  rmsd_mixed <-  rmsd_lmm 
  auc_mixed <-  auc_lmm 
  
  # get common range of data plotted
  range_rmsd <- range(rmsd_fixed, rmsd_mixed, na.rm = TRUE)
  range_auc <- range(auc_fixed, auc_mixed, na.rm = TRUE)
  
  # make labels
  # blank all non-multiples of 10, for plot
  rs <- 0 : r_max
  rs[ rs %% 10 != 0 ] <- NA
  # add nice names to all
  colnames( rmsd_fixed ) <- rs
  colnames( auc_fixed ) <- rs
  colnames( rmsd_mixed ) <- rs
  colnames( auc_mixed ) <- rs
  
  # colors
  col_fixed <- 'red'
  col_mixed <- 'blue'
  # add transarency!!!
  alpha <- 0.5
  col_fixed_alpha <- alpha(col_fixed, alpha)
  col_mixed_alpha <- alpha(col_mixed, alpha)
  
  # other shared params (same as standard boxplot)
  outline <- FALSE # no outliers plotted separately (so busy as it is)
#  range <- 0 # make whiskers extend to full range
  range <- 1.5 # default whiskers range
  whisklty <- 1 # whisker line type (default 2?)
  
  # start PDF
  myFig(name_out, width = 7, height = 5, botMar = 1.5)
  # add lower margin, so inner margins can be even smaller
  par( oma = c(1.5, 0, 0, 0) )
  # two panels
  par( mfrow = c(2, 1) )
  # boxplots!
  # top panel
  boxplot(
    rmsd_fixed,
    xlab = "",
    ylab = lab_rmsd,
    ylim = range_rmsd,
    border = col_fixed_alpha,
    col = NA,
    outline = outline,
    range = range,
    whisklty = whisklty
  )
  boxplot(
    rmsd_mixed,
    ylim = range_rmsd,
    add = TRUE,
    border = col_mixed_alpha,
    col = NA,
    outline = outline,
    range = range,
    whisklty = whisklty
  )
  # add legend to top panel only
  legend(
    legend_pos,
    c('Fixed effects (PCA)', 'Mixed effects (LMM+PCA)'),
    text.col = c( col_fixed, col_mixed ),
    bty = 'n'
  )
  # bottom panel
  boxplot(
    auc_fixed,
    xlab = '',
    ylab = lab_auc,
    ylim = range_auc,
    border = col_fixed_alpha,
    col = NA,
    outline = outline,
    range = range,
    whisklty = whisklty
  )
  boxplot(
    auc_mixed,
    ylim = range_auc,
    add = TRUE,
    border = col_mixed_alpha,
    col = NA,
    outline = outline,
    range = range,
    whisklty = whisklty
  )
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



boxplots_rmsd_auc(
  name_out = "boxplot_n_1000",
  rmsd_pca = rmsd_pca_1000,
  rmsd_lmm = rmsd_gcta_1000,
  auc_pca =  auc_pca_1000,
  auc_lmm = auc_gcta_1000,
  r_max = 90
)


boxplots_rmsd_auc(
  name_out = "boxplot_n_1000_family",
  rmsd_pca = rmsd_pca_1000_family,
  rmsd_lmm = rmsd_gcta_1000_family,
  auc_pca =  auc_pca_1000_family,
  auc_lmm = auc_gcta_1000_family,
  r_max = 90
)


boxplots_rmsd_auc(
  name_out = "boxplot_n_100",
  rmsd_pca = rmsd_pca_100,
  rmsd_lmm = rmsd_gcta_100,
  auc_pca =  auc_pca_100,
  auc_lmm = auc_gcta_100,
  r_max = 90,
  legend_pos = 'topleft'
)

