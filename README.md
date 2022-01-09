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

If you plan to update `EDClust` to newest version, after updating, please also specify `Update = TRUE` to update the related `EDClust` Julia package version.

```{r updatejulia, eval = FALSE}
julia <- setup_julia(Update = TRUE) 
```

## Feature selection

A newly developed feature selection tool, [FEAST](https://github.com/suke18/FEAST), shows superior performance in substantial real data analyses. FEAST is embedded in EDClust for feature selection. 

```{r feature_selection, message=FALSE, warning=FALSE}
data <- FEAST_select(count, subject, Ncluster = 6, Nfeature = 500)
```
FEAST bascially provides two functions for feature selection: FEAST() and FEAST_fast(), and FEAST_select() inherits most of their arguments. For extreme large dataset (sample size >5000), FEAST_fast() will be automatically applied. In our real data analyses, we could obtain satisfactory results with only 500 features, and thus by Nfeature = 500 by default.

## Example in Mlung_sub Dataset step-by-step

Here we show the key steps for baseline selection, parameter initialization and clustering.
This code chunk assumes you have a pre-processed expression count matrix called count_all_notna, an array of subject ID information called subject_all_notna.

Here's we show an example in Mlung_sub dataset step by step.

### (1) Setup the package 

```{r quick_start, eval = FALSE}
library(EDClust)
julia <- setup_julia()
## You can use `Julia` at a specific location by
## julia <- setup_julia(path = "the folder that contains Julia binary") 
```

### (2) Baseline selection

```{r quick_start, eval = FALSE}
data("Mlung_sub")
Baseline_select <- function(count_all_notna, subject_all_notna, Ncluster = 6)
## We found subject 2 has the greatest score 
## And thus we use subject 2 as the baseline subject for parameter initialization.
```

### (3) Initialize parameters

```{r quick_start, eval = FALSE}
alpha_0 <- InitVal_S(count_all_notna, subject_all_notna, Ncluster = 6, ID = 2, seed = 2345) 
```

### (4) Clustering

```{r quick_start, eval = FALSE}
result <- FitPolya(count_all_notna, subject_all_notna, alpha_0, BaseID=2L)

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



























