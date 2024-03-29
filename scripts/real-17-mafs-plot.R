# plots MAF distributions gathered in previous script

library(readr)
library(ochoalabtools)

# only thing compared to other plots is we want this alternate order
name_dir_order <- c(
    'sim-n1000-k10-f0.1-s0.5-g1',
    'sim-n100-k10-f0.1-s0.5-g1',
    'sim-n1000-k10-f0.1-s0.5-g20',
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01_sim',
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_sim',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_sim',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
)

# move to where the data is
setwd( '../data/' )

# read datasets info (names for inputs and output, colors, line types)
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# reorder as desired
indexes <- match( name_dir_order, datasets$name_dir )
datasets <- datasets[ indexes, ]

# loads individual `maf` vectors into named list for all datasets
mafs <- list()
for ( name in datasets$name_dir ) {
    setwd( name )
    load( 'maf.RData' )
    mafs[[ name ]] <- maf
    setwd( '..' )
}

plot_mafs_panel <- function( datasets, mafs, leg_cex = 0.5, leg_seg_len = 3 ) {
    # set up plot area
    plot(
        NA,
        # axis labels for first time
        xlab = 'Minor Allele Frequency',
        ylab = 'Cumulative Probability',
        main = '',
        xlim = c(0, 0.5),
        ylim = c(0, 1)
    )

    # navigate datasets
    for ( i in 1 : nrow( datasets ) ) {
        # get name used in key
        name_dir <- datasets$name_dir[ i ]
        
        # a nicer alternative to `ecdf`, more flexibility for me
        # https://stat.ethz.ch/pipermail/r-help/2016-June/439416.html
        x <- sort( mafs[[ name_dir ]] )
        y <- ppoints( x )
        
        # plot, with add as needed
        lines(
            x,
            y,
            col = datasets$col[ i ],
            lty = datasets$lty[ i ]
        )
    }

    # legend, (almost the) same as RMSD-vs-lambda fig
    # reverse order because it matches curve appearance better
    legend(
        'bottomright',
        # these two are reversed together
        rev( datasets$name_paper ),
        text.col = rev( datasets$col ),
        title = 'Dataset',
        title.col = 'black',
        pch = NA,
        col = rev( datasets$col ),
        lty = rev( datasets$lty ),
        bty = 'n',
        cex = leg_cex,
        seg.len = leg_seg_len
    )
}

# make independent plot
fig_start( 'mafs' )
plot_mafs_panel( datasets, mafs )
fig_end()

# save data and panel function for replotting as part of bigger, multipanel figure
save( plot_mafs_panel, datasets, mafs, file = 'mafs.RData' )
