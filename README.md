5 main files corresponding to Tables 1-5 in the main body respectively:

Vstat_gamma.m
VstatSin.m
OptBartlett.m
MVopt.m
MVoptChi.m

Explanation of parameters:

ns: the list of sample size considered. For example, for Vstat_gamma.m, it is set as [15,30] which are the sample sizes considered in Table 1.
qs: the quantiles of chi-square for the interested nominal level. It is set as [1.642,2.7055,3.8415] which are the 80%, 90%, and 95% quantiles of the chi-square distribution with 1 degree of freedom.
N: the number of replications to estimate the coverage probability and CI half widths.


Outputs:
cover: array of size (length(ns),length(qs), N, 3). cover(i,j,k,l) refers to whether the CI covers the truth in the k-th replication when using the i-th sample size in ns and the j-th quantile in qs. l=1 stands for EL, l=2 stands for EB, l=3 stands for TB. 

length: array of the same size as cover. cover(i,j,k,l) refers to the width of the CI in the k-th replication when using the i-th sample size in ns and the j-th quantile in qs. l=1 stands for EL, l=2 stands for EB, l=3 stands for TB. 
