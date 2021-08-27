# EDClust: An EM-MM hybrid method for cell clustering

`EDClust` is an R package designed for cell clustering in multiple-subject single-cell RNA sequencing. 
It utilizes [EDClust.jl](https://github.com/weix21/EDClust.jl) for its core routines 
to achieve a reasonable computational performance of clustering directly in R.
  
`EDClust` adopts a Dirichlet-multinomial mixture model 
and explicitly accounts for cell type heterogeneity, subject heterogeneity, and clustering uncertainty. 
Based on a rigorous staitstical framework, 
`EDClust` offers functions for predicting cell type labels,
estimating parameters of effects from different sources,
and posterior probabilities for cells being in each cluster.  
  
One might need pre-install dependent R packages such as JuliaCall, SAHRP, and mclust.

All the copyrights are explained by Xin Wei [xwei44@emory.edu](xwei44@emory.edu) from [Dr. Wu's lab](http://www.haowulab.org/). 

## Installation

### Install EDClust

```{r install, message=FALSE, warning=FALSE}
library(devtools)
install_github("weix21/EDClust", dependencies=T)
library(EDClust)
```
### Install Julia

`EDClust` requires a working installation of Julia. This can be easily done via:

```{r install, message=FALSE, warning=FALSE}
julia <- setup_julia() 
```

which will automatically install both Julia and the required Julia packages if they are missing. 

If you want to use an existing Julia binary, make sure that julia is found in the path. Then you can do 

```{r install, message=FALSE, warning=FALSE}
julia <- setup_julia(path = "the folder that contains Julia binary") 
```

To setup Julia manually, you can download a generic binary from https://julialang.org/downloads/. 
For more information about Julia setup, 
please see the julia_setup() function from [JuliaCall](https://github.com/Non-Contradiction/JuliaCall) package, 
which provides an R interface to Julia.


























