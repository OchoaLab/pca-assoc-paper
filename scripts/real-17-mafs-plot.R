# plots MAF distributions gathered in previous script

library(tibble)
library(ochoalabtools)

# names of datasets (paired list)
# data matches rmsd-vs-lambdas fig, including short names and colors
datasets <- tibble(
    name_short = c(
        'Large sample size sim.',
        'Small sample size sim.',
        'Family structure sim.',
        'Human Origins',
        'HGDP',
        'HGDP sim.',
        '1000 Genomes'
    ),
    name_long = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        'sim-n100-k10-f0.1-s0.5-g1',
        'sim-n1000-k10-f0.1-s0.5-g20',
        'HoPacAll_ld_prune_1000kb_0.3',
        #'hgdp_wgs_autosomes_ld_prune_1000kb_0.3',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_sim',
        #'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01'
    ),
    col = 1:7
)

# move to where the data is
setwd( '../data/' )

# loads individual `maf` vectors into named list for all datasets
mafs <- list()
for ( name in datasets$name_long ) {
    setwd( name )
    load( 'maf.RData' )
    mafs[[ name ]] <- maf
    setwd( '..' )
}

# start plot
fig_start( 'mafs' )

# navigate datasets
for ( i in 1 : nrow( datasets ) ) {
    # get name used in key
    name_long <- datasets$name_long[ i ]
    # extract MAF from list, convert immediately to CDF
    x <- ecdf( mafs[[ name_long ]] )
    # plot, with add as needed
    plot(
        x,
        col = datasets$col[ i ],
        add = i > 1,
        # axis labels for first time
        xlab = 'Minor Allele Frequency',
        ylab = 'Cumulative Probability',
        main = '',
        # get rid of ugly limit lines
        col.01line = NA,
        # get rid of points at knots
        do.points = FALSE,
        # connect lines in stepwise function (otherwise they're all horizontal)
        verticals = TRUE,
        # reduce large end gaps in x direction
        xaxs = 'i',
        # don't like default axis marks/locations, will do manually at end
        xaxt = 'n'
    )
}

# add x axis marks, where they make the most sense for this data
axis( 1, at = c(0, 0.25, 0.5) )

# legend, same as RMSD-vs-lambda fig
# reverse order because it matches curve appearance better
legend(
    'bottomright',
    # these two are reversed together
    rev( datasets$name_short ),
    text.col = rev( datasets$col ),
    title = 'Dataset',
    title.col = 'black',
    pch = NA,
    bty = 'n',
    cex = 0.5
)

fig_end()
