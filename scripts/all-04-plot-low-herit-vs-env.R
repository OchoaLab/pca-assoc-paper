# aggregates replicates into just mean, plots low herit vs env versions to show increase in AUC for both LMM and PCA (all r)
# - include all sims and real (but real-sim don't have these data)
# - FES is greater interest, but RC could also be included

library(readr)
library(dplyr)
library(tidyr)
library(ochoalabtools)

#################
### CONSTANTS ###
#################

# output file name (big table)
file_table <- 'sum.txt'
# directory name, needed in one mode (weird var name just stuck)
dir_phen <- 'fes'
dir_low <- 'm_causal_fac-27/h0.3'
dir_env <- paste0( dir_low, '/env0.3-0.2' )
# methods to keep in analysis
method_to_label <- list(
    'pca-plink-pure' = 'PCA',
    gcta = 'LMM'
)
methods <- names( method_to_label )
# labels of metrics
lab_rmsd <- expression( SRMSD[p] )
lab_auc <- expression( AUC[PR] )
# and simulation types (env or not)
lab_low <- ' in low herit. sims.'
lab_env <- ' in env. sims.'


############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read datasets info (names for inputs and output, colors, line types)
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# do not process "Tree" sims, which do not have the data of interest
datasets <- datasets[ datasets$type != 'Tree', ]

# big table of interest
# initialize this way, it'll grow correctly
data <- NULL

# load each dataset
for ( i in 1 : nrow( datasets ) ) {
    # enter dir
    setwd( datasets$name_dir[ i ] )

    for ( fes in c(FALSE, TRUE) ) {
        # move in one more level in this case
        if ( fes )
            setwd( dir_phen )

        for ( trait in c('lo', 'env') ) {
            # treat low as default
            dir_trait <- dir_low
            if ( trait == 'env' )
                dir_trait <- dir_env
            
            # to get back easily
            dir_pre_trait <- getwd()
            # actually go there
            setwd( dir_trait )
            
            # read the big table!
            tib <- read_tsv(
                file_table,
                col_types = 'ciiddd'
            )

            # subset to use only the two methods we talk about in the paper
            tib <- tib[ tib$method %in% methods, ]

            # summarize reps by just mean
            tib <- group_by( tib, method, pc )
            tib <- summarize( tib, rmsd = mean( rmsd ), auc = mean( auc ) )
            tib <- ungroup( tib )

            # recall the dataset of origin
            tib$dataset <- datasets$name_paper[ i ]
            # and the trait simulation type
            tib$fes <- fes
            tib$env <- trait == 'env' # true or false
            
            # concatenate into bigger table
            data <- bind_rows( data, tib )

            # get back
            setwd( dir_pre_trait )
        }
        
        # move back one more level in this case
        if ( fes )
            setwd( '..' )
    }
    # go back down
    setwd( '..' )
}

# to simplify plotting further, pivot metrics
metrics <- c('rmsd', 'auc')
data <- pivot_longer( data, all_of( metrics ), names_to = 'metric', values_to = 'value' )

###########
### FIG ###
###########

# use full page width
width <- fig_width()
fig_start(
    'low-herit-vs-env',
    width = width,
    height = width,
    mar_t = 1.5 # because there's titles and panel letters
)

# create 4 panels, to just include everything here
par( mfrow = c(2, 2) )

for ( metric in metrics ) {
    for ( fes in c(TRUE, FALSE) ) {
        # make panel for one evaluation
        # create area
        data2 <- data[ data$metric == metric & data$fes == fes, ]
        lim <- range( data2$value, na.rm = TRUE ) # range of all data to plot
        # wow, using expression variables nested in expressions is complicated!  but this works
        # https://stackoverflow.com/questions/4973898/combining-paste-and-expression-functions-in-plot-labels
        # followed solution from daniel.heydebreck (wasn't top one but includes crucial [[1]] for variables which are themselves expressions)
        lab_metric <- if (metric == 'auc') lab_auc else lab_rmsd
        # label column, first row only
        title <- ''
        if ( metric == metrics[1] )
            title <- if ( fes ) 'FES traits' else 'RC traits'
        plot(
            NA,
            xlab = eval( bquote( expression( bold( .(lab_metric[[1]]) ~ .(lab_low) ) ) ) ),
            ylab = eval( bquote( expression( bold( .(lab_metric[[1]]) ~ .(lab_env) ) ) ) ),
            xlim = lim,
            ylim = lim,
            main = title
        )
        abline( 0, 1, lty = 2, col = 'gray' )
        # navigate datasets
        # each dataset/method combo is one line, so must plot separately
        for ( dataset in datasets$name_paper ) {
            for ( method in methods ) {
                # subset as needed
                data_i <- data2[ data2$dataset == dataset & data2$method == method, ]
                # compare env to non-env!
                data_env <- data_i[ data_i$env, ]
                data_low <- data_i[ !data_i$env, ]
                # confirm alignment (pcs only)
                stopifnot( all( data_env$pc == data_low$pc ) )
                # and that order is increasing
                stopifnot( all( data_env$pc == 0L : ( nrow( data_env ) - 1L ) ) )
                # this is what we want to plot!
                lines(
                    data_low$value,
                    data_env$value,
                    col = datasets$col[ match( dataset, datasets$name_paper ) ],
                    lty = if ( method == 'gcta' ) 1 else 2 # datasets$lty[ match( dataset, datasets$name_paper ) ]
                )
                # in all these cases, mark the r=0 point with a special character
                # (we've already checked it's the first one, so subset is easy)
                points(
                    data_low$value[1],
                    data_env$value[1],
                    col = datasets$col[ match( dataset, datasets$name_paper ) ],
                    pch = 'x'
                )
            }
        }
        # add legends first panel only
        if ( metric == metrics[1] && fes ) {
            legend(
                'bottomright',
                datasets$name_paper,
                title = 'Dataset',
                col = datasets$col,
                lty = 1,
                cex = 0.8
            )
            legend(
                'right', # slightly above first legend
                unlist( method_to_label ),
                title = 'Model',
                lty = 2:1,
                cex = 0.8
            )
        }
        # add panel letter first row only
        if ( metric == metrics[1] )
            panel_letter( if ( fes ) 'A' else 'B' )
    }
}

fig_end()
