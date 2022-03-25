# install CRAN packages (include only packages not already installed)
install.packages( c('devtools', 'tibble', 'readr', 'dplyr', 'parallel', 'optparse', 'scales', 'ape', 'BEDMatrix', 'genio', 'popkin', 'bnpsd', 'simfam', 'simtrait') )

# github-only packages
library(devtools)
install_github( c('OchoaLab/ochoalabtools', 'OchoaLab/genbin') )
