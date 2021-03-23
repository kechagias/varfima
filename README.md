Purpose:
To fit bivariate vector autoregressive fractionally integrated time series with general phase. The functions of this repo
can be used to replicate the results of https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504

Contents:
The repository contains three files:
  VARFIMAModules.sas:    Definitions of necessary functions
  DemoLogLik.sas         Demo on how to compute a bivariate VARFIMA(p,D,q) likelihood with p,q<2
  DemoFitLogLik1D0.sas   Demo on how to fit a bivariate VARFIMA(1,D,0) likelihood.
  
Usage: Pull the repository to your local directory and take the following steps:
Change the path of the localDirectory macro variable in the beginning of the DemoLogLik.sas and DemoFitLogLik1D0.sas 
files to the path you used to save the pulled repo. Then simply run the demo programs.

Disclaimer: This repo is creeated only for demonstration purposes and to allow researchers to replicate the results of the 
paper https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504. As such the macros and programs contained in the repo
are not fully tested. Please feel free to reach out to me at stefanoskeh@gmail.com with any questions or suggestions.
