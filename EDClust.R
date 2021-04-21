#Using JuliaCall for calling Julia from R

if (!requireNamespace("JuliaCall", quietly = TRUE))
  install.packages("JuliaCall")

#JuliaCall Guidence
library(JuliaCall)
julia = julia_setup(JULIA_HOME = "Your Julia Path")

#or using julia_setup(installJulia = TRUE), which will invoke install_julia automatically if Julia is not found and also do initialization of JuliaCall.

#Install Julia Packages 
julia$install_package_if_needed("Distributions")
julia$install_package_if_needed("SpecialFunctions")
julia$install_package_if_needed("StatsBase")

julia$source("EDClust.jl")
#source .jl document, which includes several useful packages and EDClust function

source("SHARP_init.R")

load("Mlung_sub.RData")

#Collect the information of first individual
count = count_all_notna[,which(sample_all_notna==levels(factor(sample_all_notna))[1])]
celllabel = annot_all_notna[which(sample_all_notna==levels(factor(sample_all_notna))[1])]
sampleid = sample_all_notna[which(sample_all_notna==levels(factor(sample_all_notna))[1])]

##Normalization for SHARP
data = normalizeSC.total(count)
data = log2(data$normdata+1)
Ncluster = length(unique(celllabel))

#Using cluster results given by SHARP to compute initial alpha0
#Number of cell types is predetermined in this steps
SLabel = SHARP(data, N.cluster=Ncluster,rN.seed=1234)
alpha_0 = computeInitVal(count,SLabel$pred_clusters,sampleid)$alpha[[1]]


EDClust_result = julia$call("fitPolya",count_all_notna,sample_all_notna,alpha_0,100L,3L,1e-4)
#function fitPolya(count,sampleid,initial_alpha0,...)
#100L: upper limit of EM iterations 
#3L: upper limit of MM iterations 
#1e-4: stopping criterion
#100L, 3L, 1e-4 are the default settings

#EMMM_result is a list contains the following:
#1. cell labels
#2. Likelihood
#3. alpha_0: cell type effect
#4. delta_l: individual l's effect
#5. alpha_l: alpha_0 + delta_i
#6. log(p): p_kl indicates the probability of assigning a cell from individual l to group k.

#ARI
adjustedRandIndex(EDClust_result[[1]], annot_all_notna)






