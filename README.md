# WMCapacity (WoMMBAT)
A GUI R implementation of hierarchical Bayesian models of working memory, used for analyzing change detection data.

See Morey (2011) for modelling details (paper here: http://www.sciencedirect.com/science/article/pii/S0022249610001094), and Morey & Morey (2011) for the software manual (paper here, open access: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218285/). 

### Installing

To install the latest stable version from CRAN, use `install.packages`:

```R
install.packages('WMCapacity', dependencies = TRUE)
```
or your R graphical user interface's install packages menu.

To install the latest development version, you can use `install_github` from the `devtools` package:

```R
## install devtools if necessary
install.packages('devtools')
## Load devtools package for install_github()
library(devtools)
## install from github
install_github('richarddmorey/WMCapacity', subdir='WMCapacity', dependencies = TRUE)
```
