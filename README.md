**Repository Version 1.0** <br>

**Purpose**:<br>
&nbsp; Provide SAS code to fit bivariate vector autoregressive fractionally integrated time series with general phase.  <br>
&nbsp; The SAS functions in the VARFIMAModules.sas file implement the methods of <a href="https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504">Kechagias and Pipiras (2019)</a> <br> 
&nbsp; and of <a href="https://www.sciencedirect.com/science/article/abs/pii/S0165168410004019">Helgason etal (2011)</a>.


**Contents**:<br>
&nbsp;   The repository contains three files: <br>
&nbsp; &nbsp; &nbsp;    1. <i>VARFIMAModules.sas</i>  &nbsp; &nbsp; &nbsp; Definitions of necessary functions <br>
&nbsp; &nbsp; &nbsp;    2. <i>DemoLogLik.sas</i>      &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Demo on how to compute a bivariate VARFIMA(p,D,q) likelihood with p,q<2. <br>
&nbsp; &nbsp; &nbsp;    3. <i>DemoFitLogLik1D0.sas</i> &nbsp; &nbsp; Demo on how to fit a bivariate VARFIMA(1,D,0) likelihood.
  
**Usage**: <br>
&nbsp;  1. Pull the repository to your local directory.  <br>
&nbsp;  2. Open the DemoLogLik.sas and DemoFitLogLik1D0.sas files.  <br>
&nbsp;  3. Change the path of the <i>localDirectory</i> macro variable in the beginning of these files to your local directory.  <br>
&nbsp;  4. Run the demo files.

**Disclaimer**: <br>
&nbsp;  This repo is created only for demonstration purposes. As a result the macros and programs contained in the repo <br>
&nbsp;  are not fully tested and may contain bugs or defects. If you ecnounter any problems or if you have any questions <br> 
&nbsp;  or comments please do not hesitate to reach out to me at stefanoskeh@gmail.com.

**SAS Version**: <br>
&nbsp;  All functions were developed under SAS 9.4.

**References**: <br>
&nbsp;  Kechagias S. and Pipiras V. (2020), Modeling bivariate long‐range dependence with general phase. <i>Journal of Time <br> 
&nbsp;  Series Analysis</i> 41: 268-292. <br>
&nbsp;  Helgason H., Pipiras V. and Abry P. (2011), Fast and exact synthesis of stationary multivariate Gaussian time series <br>
&nbsp;  using circulant embedding. <i>Signal Processing</i> 91(5): 1123–1133.
