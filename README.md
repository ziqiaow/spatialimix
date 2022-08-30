# Spatial IMIX
Spatial IMIX: A Mixture Model Approach to Spatially Correlated Multi-Omics Data Integration


### Installation
Install the spatialimix package via github
```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("ziqiaow/spatialimix", build_vignettes = TRUE) 
```
To save time, download the package without building the vignettes by setting 'build_vignettes = FALSE'.

### To get started with the package
Load the package and open the package vignette:
```
library("spatialimix")
vignette("spatialimix")
```
There will be examples for all the functions in spatialimix package in [vignette](vignettes/spatialimix.pdf). If you have any question, use help().

**Reference**

* Ziqiao Wang and Peng Wei. Spatial IMIX: A Mixture Model Approach to Spatially Correlated Multi-Omics Data Integration.
