# makes sure all p-value vectors have the right length (number of loci in data)

library(optparse)
library(genio)

# constants
methods <- c('pca-plink-pure', 'gcta') # 'pca-plink', 
# the name is for dir only, actual file is just "data"
name_in <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--n_pcs", type = "integer", default = 90,
                help = "Max number of PCs", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 50,
                help = "Max replicates", metavar = "int"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters location only)"),
    make_option("--final", action = "store_true", default = FALSE, 
                help = "Dies if files are missing (otherwise they are skipped silently)"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option(c('-l', "--labs"), action = "store_true", default = FALSE, 
                help = "Include LMM with labels data")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2
labs <- opt$labs

# do this consistency check early
if ( !is.na( env1 ) && is.na( env2 ) )
    stop( 'If --env1 is specified, must also specify --env2!' )

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# add a third method in this case
if ( labs )
    methods <- c( methods, 'gcta-labs' )

# move to where the data is
setwd( '../data/' )
setwd( name )

# remember this location, using global path, for easily returning across multiple levels when done
dir_base <- getwd()

# in simulations, each rep has its own bim file, but they all have the same number of loci, so let's just read it from rep-1
if (opt$sim)
    name_in <- paste0( 'rep-1/', name_in )
# get number of loci
m_loci <- count_lines( name_in, 'bim' )

for ( rep in 1 : rep_max ) {
    # specify location of files to process, as many levels as needed
    dir_out <- paste0( 'rep-', rep, '/' )
    if ( fes )
        dir_out <- paste0( dir_out, 'fes/' )
    if ( m_causal_fac != 10 )
        dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
    if ( herit != 0.8 )
        dir_out <- paste0( dir_out, 'h', herit, '/' )
    if ( !is.na( env1 ) )
        dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )
    
    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) ) {
        if ( opt$final ) {
            # fatal if file is missing and we set `--final` option
            stop( 'whole rep ', rep, ' MISSING!' )
        } else 
            next
    }
    setwd( dir_out )
    
    # start a big loop
    for ( method in methods ) {
        # only one method is not expected to have PCs (important for final validation)
        n_pcs_max_method <- if ( method == 'gcta-labs' ) 0 else n_pcs_max
        for ( n_pcs in 0 : n_pcs_max_method ) {
            # file to read
            file_pvals <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )

            if ( !file.exists( file_pvals ) ) {
                if ( opt$final ) {
                    # fatal if file is missing and we set `--final` option
                    stop(
                        'rep: ', rep,
                        ', method: ', method,
                        ', pcs: ', n_pcs,
                        ', MISSING!'
                    )
                } else 
                    next
            }

            # count lines using linux system terminal (fastest)
            # run wc on terminal
            n_lines <- system(
                paste0(
                    'zcat ',
                    file_pvals,
                    ' | wc -l'
                ),
                intern = TRUE
            )
            # returns string, turn into numeric
            n_lines <- as.numeric( n_lines )

            # data is either complete (m_loci) or result was NULL (1 line), nothing else is acceptable
            if ( n_lines != m_loci && n_lines != 1 ) {
                # this is what we wanted to find!
                stop(
                    'rep: ', rep,
                    ', method: ', method,
                    ', pcs: ', n_pcs,
                    ', lines: ', n_lines,
                    ', expected: ', m_loci
                )
            }
        }
    }
    
    # return to base when done with this rep
    setwd( dir_base )
}
