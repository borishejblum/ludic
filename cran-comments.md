# This is an update of the ludic package updating citations 

## Test environments  
 * local R installation, R 4.4.3
 * Linux (Ubuntu 24.04), macOS (14.7) and Windows (Server 2022 10.0), R devel and release (through GitHub Actions)

## R CMD check results  
0 ERRORs | 0 WARNINGs | 0 NOTEs


## Reverse dependencies  
There are no reverse dependencies.

## Installed size  
On some architecture, the CHECK returns one NOTE because 
installed size is too big. My understanding is that this 
inflation of the libs subdirectory is due to the use of 
Rcpp. Indeed, some functions of the ludic package have 
been written in C++ using Rcpp. They are needed to perform 
efficient matching. Without the speed up gained from 
those C++ functions, this package would become impractical.


Thanks, Boris Hejblum

---
