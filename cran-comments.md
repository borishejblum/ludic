# This is an update of the ludic package fixing errors in r-devel (4.0.0) Debian checks 

## Test environments  
* local macOS 10.15.1 install, R 3.6.1
* ubuntu 16.04.6 (on travis-ci), R devel and release
* Windows Visual Studio 2015 (on appveyor), R devel and release

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
