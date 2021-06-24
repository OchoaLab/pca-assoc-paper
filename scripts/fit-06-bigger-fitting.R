# a slower, empirical approach for deciding on F for undifferentiating MAFs to result in more reasonable structures in terms of popkin and also MAF distribution.
# computes table of sum of squared errors for both kinship and MAF, so we're explicitly trying to fit both/either
# MAF ascertainment is built into test, which is actually what makes this so hard (without it, popkin fits are generally great, so it would just be MAF distributions to worry about)

library(optparse)
library(bnpsd)
library(genio)
library(readr)
library(popkin)
library(simtrait) # for allele_freqs
library(tibble)

# hardcoded params
verbose <- TRUE # to print messages
# hardcoded maf_real = TRUE

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 10000, 
                help = "number of loci", metavar = "int"),
    make_option("--maf_min", type = "double", default = 0.01, 
                help = "Minimum MAF, to simulate MAF-based ascertainment bias", metavar = "double"),
    make_option("--beta", action = "store_true", default = FALSE, 
                help = "Fit symmetric Beta distributions instead of `bnpsd::undiff_af`"),
    make_option("--distr", action = "store_true", default = FALSE, 
                help = "Pass AFs from `bnpsd::undiff_af` to `draw_all_admix` as distribution instead of fixed vector.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
m_loci <- opt$m_loci
maf_min <- opt$maf_min
fit_beta <- opt$beta
as_distr <- opt$distr

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# don't try to do two things at once!
if ( fit_beta && as_distr )
    stop( 'Cannot use `--beta` and `--distr` together!' )

##############################
### LOAD PRECALCUATED DATA ###
##############################

# directory above must already exist
setwd( '../data/' )

# deduce original file path
name_real <- sub( '_sim$', '', name )
setwd( name_real )

# load MAF
load( 'maf.RData' )
# loads: maf
# rename for consistency
maf_real <- maf
# subsample and sort for a later comparison of distributions
maf_real_m <- sort( sample( maf_real, m_loci ) )

# also load kinship estimate from real data
kinship_real <- read_grm( 'popkin' )$kinship

# go back down so we can find the simulated path
setwd( '..' )

setwd( name )

# load precalculated bnpsd data from RData file
load( 'bnpsd.RData' )
# loads: admix_proportions, tree_subpops, fam

sim_data_measure_fit <- function(dat, i, p_anc = NULL, p_anc_distr = NULL, beta = NA) {
    # draw allele freqs and genotypes
    if (verbose)
        message('draw_all_admix')
    # simulate data from this
    # only difference is here number of loci will be smaller than real data, on purpose (and their MAF distribution will be practically uniform too)
    X <- draw_all_admix(
        admix_proportions = admix_proportions,
        tree_subpops = tree_subpops,
        m_loci = m_loci,
        p_anc = p_anc,
        p_anc_distr = p_anc_distr,
        maf_min = maf_min,
        verbose = verbose,
        beta = beta
    )$X

    # now estimate kinship
    message( 'popkin...' )
    kinship <- popkin( X, fam$fam )
    
    # calculate desired allele frequencies
    # sorted cause we only want to compare distributions
    message( 'allele_freqs...' )
    maf <- sort( allele_freqs( X, fold = TRUE ) )
    
    # calculate and store mean errors
    dat$kin[i] <- mean( ( kinship_real - kinship )^2 )
    dat$maf[i] <- mean( ( maf_real_m - maf )^2 )

    return( dat )
}


if ( fit_beta ) {
    # specify in log space
    alphas <- 10^seq( -2, 2, by = 0.1 )
    
    # desired tibble, with parameter of interest and error measures
    dat <- tibble( alpha = alphas, kin = NA, maf = NA )

    # slowest part: simulate data and store error measurements
    for ( i in seq_len( length( alphas ) ) ) {
        message( 'alpha = ', alphas[i] )
        
        # use built-in ability to draw from Beta!
        dat <- sim_data_measure_fit( dat, i, beta = alphas[i] )
    }
} else {
    # need FST of this simulation, which has been calculated already
    # this is sort of like a minimum, in that we might want to un-differentiate even more
    fst_tree <- as.numeric( read_lines( 'fst.txt' ) )
    # round for simplicity
    fst_min <- round( fst_tree, 2 )
    # this is the max allowed by `undiff_af`:
    fst_max <- 4 * mean( ( maf_real - 0.5 )^2 )

    fsts <- seq( fst_min, fst_max, by = 0.01 )
    message( 'Fst seq: ', toString( fsts ) )

    # desired tibble, with parameter of interest and error measures
    dat <- tibble( fst = fsts, kin = NA, maf = NA )

    # slowest part: simulate data and store error measurements
    for ( i in seq_len( length( fsts ) ) ) {
        message( 'F = ', fsts[i] )
        
        # undifferentiate whole original distribution
        p_anc <- undiff_af( maf_real, fsts[i] )$p

        if ( as_distr ) {
            # don't pass subsampled, instead draw from p_anc inside function as distribution (with replacement)
            dat <- sim_data_measure_fit( dat, i, p_anc_distr = p_anc )
        } else {
            # select random allele frequencies for a smaller collection of loci
            p_anc <- sample( p_anc, m_loci )
            # pass as fixed p_anc vector (re-draws don't change it!)
            dat <- sim_data_measure_fit( dat, i, p_anc = p_anc )
        }
    }
}

# hopefully results are similar for smaller m to run faster tests!
file_out <- paste0( if ( fit_beta ) 'B' else if ( as_distr ) 'Fd' else 'F', '-fit_m-', m_loci, '.txt' )
write_tsv( dat, file_out )
