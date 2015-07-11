# WMCapacity
A GUI R implementation of hierarchical Bayesian models of working memory, used for analyzing change detection data.


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
