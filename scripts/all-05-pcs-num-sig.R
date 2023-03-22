# test PCs against the trait under a "null" model with no loci

library(optparse)
library(readr)
library(genio)
library(qvalue)
library(ochoalabtools)

# the name is for dir only, actual file is just "data"
name_in <- 'data'
# navigate all reps
rep_max <- 50
# load file with all PCs
n_pcs_max <- 90
# use q-values for threshold
q_cut <- 0.05 # a standard value

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2

# do this consistency check early
if ( !is.na( env1 ) && is.na( env2 ) )
    stop( 'If --env1 is specified, must also specify --env2!' )

# this is the same for all datasets, reps!
# adjust paths if using fes model or non-default m_causal_fac or herit
dir_phen <- '' # current dir
# use subdir instead in this case
if ( fes )
    dir_phen <- 'fes/'
if ( m_causal_fac != 10 )
    dir_phen <- paste0( dir_phen, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_phen <- paste0( dir_phen, 'h', herit, '/' )
if ( !is.na( env1 ) )
    dir_phen <- paste0( dir_phen, 'env', env1, '-', env2, '/' )

# move to where the data is
setwd( '../data/' )

# read datasets info (names for inputs and output, colors, line types)
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )

# output list
qvals_num_sig <- list()

for ( i in 1 : nrow( datasets ) ) {
    # gets reused often
    name_paper <- datasets$name_paper[ i ]
    # go where the data is
    message( datasets$name_dir[ i ] )
    setwd( datasets$name_dir[ i ] )

    # matrix of p-values, every rep and PC
    pvals <- matrix( NA, nrow = rep_max, ncol = n_pcs_max )
    # estimate q-values too!
    qvals <- matrix( NA, nrow = rep_max, ncol = n_pcs_max )

    # navigate reps
    for ( rep in 1L : rep_max ) {

        # move higher to the "reps" location
        # this is so GCTA's temporary files don't overwrite files from other parallel runs
        dir_out <- paste0( 'rep-', rep )
        setwd( dir_out )

        # genotypes, PCs:
        # - in real data, are all in lower level (shared across reps)
        # - in simulated data, are all in current level (not shared across reps)
        name_in_lower <- if ( datasets$type[ i ] != 'Real' ) name_in else paste0( '../', name_in )

        # load all PCs
        file_pcs <- paste0( name_in_lower, '-plink-n_pcs_', n_pcs_max )
        pcs <- read_eigenvec( file_pcs )$eigenvec

        # load phenotype
        name_phen <- paste0( dir_phen, name_in )
        phen <- read_phen( name_phen )$pheno

        # start statistical tests
        # tests current PC against all previous PCs as null
        # first p-value is different
        pvals[ rep, 1L ] <- anova( lm( phen ~ pcs[,1]) )[1L,'Pr(>F)']
        # all rest
        for ( pc in 2L : n_pcs_max ) {
            pvals[ rep, pc ] <- anova( lm( phen ~ pcs[ , 1L:(pc-1L) ] ), lm( phen ~ pcs[ , 1L : pc ] ) )[2L,'Pr(>F)']
        }

        # use q-values to be more permissive

        # set pi0=1 since number of PCs is relatively small, don't want wild qvalue estimates due to that
        # also solves this error
        ## Error in smooth.spline(lambda, pi0, df = smooth.df) : 
        ##   missing or infinite values in inputs are not allowed

        # apply to each rep separately
        # easier to explain conceptually, and for boxplots/etc
        qvals[ rep, ] <- qvalue( pvals[ rep, ], pi0 = 1 )$qvalues

        # go back down
        setwd( '..' )
    }

    # set threshold, count per row/rep
    qvals_num_sig[[ name_paper ]] <- rowSums( qvals < q_cut )
    
    # go back down
    setwd( '..' )
}

# save data!

# place it in a good location that matches the parameters we looked at for the trait
# these directories should exist already if we're doing this towards the end of our analysis
if ( dir_phen != '' )
    setwd( dir_phen )

# being lazy here, instead of transforming into a nice table, just save list as Rdata
save( qvals_num_sig, file = 'pcs-num-sig.RData' )

# make figure!

# first dataset should have fill be white, not black (default)!
datasets$col[1] <- 0

# get max width
width <- fig_width() / 2
fig_start(
    'pcs-num-sig',
    width = width,
    height = width,
    mar_b = 9
)
boxplot(
    qvals_num_sig,
    #xlab = 'Dataset',
    ylab = 'Num. Sig. PCs',
    las = 2,
    col = datasets$col
)
mtext( 'Dataset', side = 1, line = 8 )
fig_end()
