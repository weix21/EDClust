# EDClust: An EM-MM hybrid method for cell clustering

`EDClust` is an R package designed for cell clustering in multiple-subject single-cell RNA sequencing. 
It utilizes [EDClust.jl](https://github.com/weix21/EDClust.jl) for its core routines 
to achieve a reasonable computational performance of clustering directly in R.
  
`EDClust` adopts a Dirichlet-multinomial mixture model 
and explicitly accounts for cell type heterogeneity, subject heterogeneity, and clustering uncertainty. 
EDClust pipeline includes three steps: (A). input data. (B). initialize parameters. (C). clustering.
Based on an EM and MM hybrid framework, 
`EDClust` offers functions for predicting cell type labels,
estimating parameters of effects from different sources,
and posterior probabilities for cells being in each cluster.  

![image](https://github.com/weix21/EDClust/blob/main/vignettes/flowchart.jpg)
  
In this readme file, we briefly present how to install `EDClust` package through GitHub and the basic functionalities. 
One might need pre-install dependent R packages such as JuliaCall, SAHRP, and mclust.

All the copyrights are explained by Xin Wei <xwei44@emory.edu> from [Dr. Wu's lab](http://www.haowulab.org/). 
Any EDClust-related questions should be posted to the GitHub Issue section of `EDClust`
homepage at https://github.com/weix21/EDClust/issues.

## Installation

### Install EDClust

```{r install, message=FALSE, warning=FALSE}
library(devtools)
install_github("weix21/EDClust")
library(EDClust)
```
### Install Julia

`EDClust` requires a working installation of Julia. This can be easily done via:

```{r install, message=FALSE, warning=FALSE}
julia <- setup_julia() 
```

which will automatically install both Julia and the required Julia packages if they are missing. 

If you want to setup Julia manually, you can download a generic binary from https://julialang.org/downloads/. 
Before using an existing Julia binary, make sure that Julia is found in the path. Then you can do 

```{r install, message=FALSE, warning=FALSE}
julia <- setup_julia(path = "the folder that contains Julia binary") 
```

For more information about Julia setup, 
please see the julia_setup() function from [JuliaCall](https://github.com/Non-Contradiction/JuliaCall) package, 
which provides an R interface to Julia.

## Quick start

Here we show the key steps for parameter initialization and clustering.
This code chunk assumes you have a pre-processed expression count matrix called count_all_notna, an array of subject ID information called subject_all_notna.

Here's we show an example in Mlung_sub dataset step by step.

### (1) Setup the package 

```{r quick_start, eval = FALSE}
library(EDClust)
julia <- setup_julia()
## You can use `Julia` at a specific location by
## julia <- setup_julia(path = "the folder that contains Julia binary") 
```

### (2) Initialize parameters

```{r quick_start, eval = FALSE}
data("Mlung_sub")
alpha_0 <- InitVal_S(count_all_notna, subject_all_notna, Ncluster = 6, ID = 2, seed = 2345) 
```

### (3) Clustering

```{r quick_start, eval = FALSE}
result <- FitPolya(count_all_notna, subject_all_notna, alpha_0)

library(mclust)
adjustedRandIndex(result$mem, annot_all_notna)
```

Based on the clustering result, The t-SNE plot can be shown as: 

![image](https://github.com/weix21/EDClust/blob/main/vignettes/TSNE.jpg)

### For detailed usage of EDClust, please refer to the vignette file through

```{r vignettes, eval = FALSE}
vignette("EDClust")
## You might need to use
## install_github("weix21/EDClust", build_vignettes = TRUE)
## to build vignettes, which could be time consuming and may require additional packages.
```

The content in this README file is essentially the same as the package vignette.



























