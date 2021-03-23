**Purpose**:<br>
&nbsp; To fit bivariate vector autoregressive fractionally integrated time series with general phase. The functions of this repo 
&nbsp; can be used to replicate the results of https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504

**Contents**:<br>
&nbsp;   The repository contains three files: <br>
&nbsp; &nbsp; &nbsp;    1. <i>VARFIMAModules.sas</i>:    Definitions of necessary functions <br>
&nbsp; &nbsp; &nbsp;    2. DemoLogLik.sas         Demo on how to compute a bivariate VARFIMA(p,D,q) likelihood with p,q<2. <br>
&nbsp; &nbsp; &nbsp;    3. DemoFitLogLik1D0.sas   Demo on how to fit a bivariate VARFIMA(1,D,0) likelihood.
  
**Usage**: <br>
&nbsp;  Pull the repository to your local directory and take the following steps:
&nbsp;  Change the path of the localDirectory macro variable in the beginning of the DemoLogLik.sas and DemoFitLogLik1D0.sas 
&nbsp;  files to the path you used to save the pulled repo. Then simply run the demo programs.

**Disclaimer**: <br>
&nbsp;  This repo is creeated only for demonstration purposes and to allow researchers to replicate the results of the 
&nbsp;  paper https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504. As such the macros and programs contained in the repo
&nbsp;  are not fully tested. Please feel free to reach out to me at stefanoskeh@gmail.com with any questions or suggestions.

