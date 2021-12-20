library(popkin)
library(genio)
library(readr)
library(ochoalabtools)

# go where the data is
setwd( '../data/' )

# load MAF data
load( 'mafs.RData' )
# loads data `datasets`, `mafs` and function `plot_mafs_panel`

# for the first simulated datasets we just need these
sim1 <- read_grm( 'sim-n1000-k10-f0.1-s0.5-g1/rep-1/popkin', verbose = FALSE )
sim3 <- read_grm( 'sim-n1000-k10-f0.1-s0.5-g20/rep-1/popkin', verbose = FALSE )

# list of data to plot so far
data <- list(
    inbr_diag( sim1$kinship ),
    inbr_diag( sim3$kinship ),
    # maf plot requires different CEX in multipanel context commpared to single panel version
    function() plot_mafs_panel( datasets, mafs, leg_cex = 0.45, leg_seg_len = 5 )
)
# only used for real datasets
labs <- list(
    NULL,
    NULL,
    NULL
)
# hacked for now, could be more elegant and/or impose more consistency with `datasets`
titles <- c(
    'Admix. Large/Small sim.',
    'Admix. Family sim.',
    '',
    'Human Origins',
    'Human Origins sim.',
    'Human Origins sim.',
    'HGDP',
    'HGDP sim.',
    'HGDP sim.',
    '1000 Genomes',
    '1000 Genomes sim.',
    '1000 Genomes sim.'
)

# for a special loop for each real dataset
# rows get added in this order!
names_real <- c(
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
)

for ( name_real in names_real ) {

    #####################
    ### READ ALL DATA ###
    #####################
    
    setwd( name_real )

    # read full fam data
    fam <- read_fam( 'data', verbose = FALSE )
    
    # read existing popkin data
    kinship_real <- read_grm( 'popkin', verbose = FALSE )$kinship

    # read annotations
    subpop_info <- read_tsv( 'pops-annot.txt', comment = '#', show_col_types = FALSE )

    # also read tree from here (there's a copy in the next _sim dataset, but meh)
    load( 'tree.RData' )

    # now switch dirs to simulated data
    setwd( '..' )
    setwd( paste0( name_real, '_sim' ) )
    # data is always in first rep
    setwd( 'rep-1' )

    # read existing popkin data from simulation
    kinship_sim <- read_grm( 'popkin', verbose = FALSE )$kinship

    # go back down for next dataset (two levels because of "rep-1/")
    setwd( '../..' )

    ###############
    ### REORDER ###
    ###############

    # NOTE: this takes the order of the human papers, which disagree a bit with the tree!
    
    # map subpopulations using sub-subpopulations
    fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]

    # confirm matching orders here too
    stopifnot( fam$id == rownames( kinship_real ) )
    stopifnot( fam$id == rownames( kinship_sim ) )
    # reorder individuals so subpopulations come out in desired order:
    indexes <- order( match( fam$fam, subpop_info$pop ) )
    fam <- fam[ indexes, ]
    kinship_real <- kinship_real[ indexes, indexes ]
    kinship_sim <- kinship_sim[ indexes, indexes ]

    # transform diagonal just before plotting
    kinship_real <- inbr_diag( kinship_real )
    kinship_sim <- inbr_diag( kinship_sim )

    ## # DEBUG report ranges cause I don't know what these are
    ## message( name_real, ': ', max( kinship_real ) )
    ## message( paste0( name_real, '_sim' ), ': ', max( kinship_sim ) )

    # shorten labels for nicer plots
    # this is an HGDP case
    fam$superpop[ fam$superpop == 'CENTRAL_SOUTH_ASIA' ] <- 'SOUTH_ASIA'
    
    # dump all current data into list to plot
    data <- c( data, list( kinship_real, kinship_sim, tree ) )
    # extend labels too
    labs <- c( labs, list( fam$superpop, fam$superpop, NULL ) )
}

# good kinship estimate has zero min (no capping necessary), only look at upper range and trust diagonal
# NOTE: based on HO only (had highest range)
alpha <- 0.01
capHi <- quantile( diag( data[[4]] ), probs = 1 - alpha )
# apply cap to all kinship matrices
# this hand-picks the kinship matrices
for ( i in c( 1, 2, 4, 5, 7, 8, 10, 11 ) ) {
    data[[ i ]][data[[ i ]] > capHi] <- capHi
}

# scaling factors for labels in the different datasets
cex_ho <- 0.4
cex_hgdp <- 0.3
cex_tgp <- 0.5
labs_cex <- list(
    1, 1, 1,
    cex_ho, cex_ho, cex_ho,
    cex_hgdp, cex_hgdp, cex_hgdp,
    cex_tgp, cex_tgp, cex_tgp
)

# same but this one applies to tree tip labels only
names_cex <- rep.int( 1, 12 )
names_cex[  6 ] <- 0.1 # HO
names_cex[  9 ] <- 0.3 # HGDP
names_cex[ 12 ] <- 0.5 # TGP

# x-axis labels (per panel)
ylab <- c( rep.int('', 9), "Individuals", "Individuals", "Kinship" )

# get dimensions to stretch to full page fig retaining desired ratio
wh <- fig_scale( 3/4, 'plos' )

fig_start(
    'kinship',
    width = wh[ 1 ],
    height = wh[ 2 ],
    mar_t = 2
)

# shrink everything a little bit, titles and other things get too messed up
par( cex = 0.8 )

plot_popkin(
    data,
    titles = titles,
    names_cex = names_cex, # in this case affects tree tips only
    labs = labs,
    labs_line = 0.2,
    labs_las = 2,
    labs_cex = labs_cex,
    labs_lwd = 0.1,
#    labs_even = TRUE,
    ylab = ylab,
    ylab_line = 2.5, # place below `labs_line` (above)
    ylab_side = 1, # x-axis instead of y-axis
    layout_rows = 4,
    leg_width = 0.3,
    leg_column = 3,
    panel_letters_adj = -0.25, # have bigger margins here, can make more negative!
    oma = c(0.5, 0, 0, 0) # move outer margin to bottom only (default is 1.5 for left only)
)

fig_end()


### subsets for slides ###


# this one is bigger than original, since titles are big
fig_start(
    'kinship-admix-sims',
    width = wh[ 1 ],
    height = wh[ 2 ] * 3 / 8,
    mar_t = 2,
    mar_b = 0,
    mar_l = 0
)
plot_popkin(
    data[ 1:2 ],
    titles = titles[ 1:2 ],
    leg_width = 0.25,
    panel_letters_adj = 0 # reduced margins here
)
fig_end()

# covers real cases

plot_kinship_real_sim_tree <- function( name, indexes ) {
    fig_start(
        paste0( 'kinship-real-sim-tree_', name ),
        width = wh[ 1 ],
        height = wh[ 2 ] * 1.1 / 4, # same as main figure but only one row
        mar_t = 2
    )
    # shrink everything a little bit, titles and other things get too messed up
    par( cex = 0.8 )
    plot_popkin(
        data[ indexes ],
        titles = titles[ indexes ],
        names_cex = names_cex[ indexes ], # in this case affects tree tips only
        labs = labs[ indexes ],
        labs_line = 0.2,
        labs_las = 2,
        labs_cex = labs_cex[ indexes ],
        labs_lwd = 0.1,
        ylab = ylab[ 10:12 ], # always use last row here
        ylab_line = 3.1, # place below `labs_line` (above)
        ylab_side = 1, # x-axis instead of y-axis
        leg_width = 0.3,
        leg_column = 3,
        panel_letters_adj = -0.25, # have bigger margins here, can make more negative!
        oma = c(1, 0, 0, 0) # move outer margin to bottom only (default is 1.5 for left only)
    )
    fig_end()
}

plot_kinship_real_sim_tree( 'ho', 4:6 )
plot_kinship_real_sim_tree( 'hgdp', 7:9 )
plot_kinship_real_sim_tree( 'tgp', 10:12 )
