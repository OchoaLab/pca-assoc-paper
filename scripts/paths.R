# this will make things work on DCC as well as all my local machines!
dir_bin <- paste0( path.expand("~"), '/bin/' )

# another trick to use local plink installation
# (will be a slightly different version, but the one I normally install doesn't take advantage of vectorizations because that one doesn't work on my home desktop, so this one on DCC should be better)

# tests if a path leads to a functioning binary!
# works for paths in PATH, and for full paths regardless
which_bin <- function( bin ) {
    # don't want verbosity, just a yes/no answer
    ret <- system(
        paste0( 'which ', bin ),
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
    )
    # TRUE if all good, FALSE otherwise
    return ( ret == 0 )
}

# define some binary paths
# 1) different on DCC than local machines
# actually just asks if the bare name is in the PATH, otherwise tries fuller path on my directory
plink2_bin    <- 'plink2'
if ( !which_bin( plink2_bin ) )
    plink2_bin <- paste0( dir_bin, plink2_bin )
# don't check that the second binary exists, should only check if we need it

# 2) same on all machines
# (don't check that binaries exist always, only check when needed)
plink1_bin    <- paste0( dir_bin, 'plink1/plink' )
emmax_bin     <- paste0( dir_bin, 'emmax-beta-07Mar2010/emmax' )
emmax_kin_bin <- paste0( dir_bin, 'emmax-beta-07Mar2010/emmax-kin' )
gemma_bin     <- paste0( dir_bin, 'gemma-0.98.1-linux-static' )
gcta_bin      <- paste0( dir_bin, 'gcta_1.93.2beta/gcta64' )
bolt_bin      <- paste0( dir_bin, 'BOLT-LMM_v2.3.4/bolt' )
smartpca_bin  <- paste0( dir_bin, 'EIG-7.2.1/bin/smartpca' )
twstats_bin   <- paste0( dir_bin, 'EIG-7.2.1/bin/twstats' )
twstats_tab   <- paste0( dir_bin, 'EIG-7.2.1/POPGEN/twtable' )
