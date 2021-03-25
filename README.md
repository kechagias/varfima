**Repository Version 1.0** <br>

**Purpose**:<br>
&nbsp; Provide SAS code for synthesizing stationary multivariate time series from a given covariance structure and for fitting  <br>
&nbsp; bivariate vector autoregressive fractionally integrated time series with general phase. The SAS functions in the  <br> 
&nbsp; VARFIMAModules.sas file implement the methods of <a href="https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504">Kechagias and Pipiras (2019)</a> and of <a href="https://www.sciencedirect.com/science/article/abs/pii/S0165168410004019">Helgason etal (2011)</a>.


**Contents**:<br>
&nbsp;   The repository contains three files: <br>
&nbsp; &nbsp; &nbsp;    1. <i>VARFIMAModules.sas</i>  &nbsp; &nbsp; &nbsp; Contains definitions of necessary functions. <br>
&nbsp; &nbsp; &nbsp;    2. <i>DemoLogLik.sas</i>      &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Demo on how to compute a bivariate VARFIMA(p,D,q) likelihood with p,q<2. <br>
&nbsp; &nbsp; &nbsp;    3. <i>DemoFitLogLik1D0.sas</i> &nbsp; &nbsp; Demo on how to fit a bivariate VARFIMA(1,D,0) likelihood.
  
**Usage**: <br>
&nbsp;  1. Pull the repository to your local directory. If you are working with SAS On Demand for Academics you will need to <br>
&nbsp;  &nbsp; &nbsp;  upload the files to the SAS server via SAS Studio. Check out this helpful <a href="https://support.sas.com/ondemand/manuals/UploadingDataUsers.pdf">guide</a> for more details. <br>
&nbsp;  2. Open the DemoLogLik.sas and DemoFitLogLik1D0.sas files.  <br>
&nbsp;  3. Change the path of the macro variable <i>localDirectory</i> located in the beginning of these files to the path of your <br> 
&nbsp; &nbsp; &nbsp; local directory (where you saved the repo after you downloaded it) or to the path of the SAS Server directory <br> 
&nbsp; &nbsp; &nbsp;  (the SAS Studio folder that you uploaded your files into). To find the directory of the SAS Studio folder, <br> 
&nbsp; &nbsp; &nbsp;  right click on it and press properties. <br>
&nbsp;  4. Run the demo files.

**Disclaimer**: <br>
&nbsp;  This repo is created for demonstration purposes. As a result the functions and programs contained in the repo <br>
&nbsp;  are not fully tested and may contain bugs or defects. If you ecnounter any problems or if you have any questions <br> 
&nbsp;  or comments please do not hesitate to reach out to me at stefanoskeh@gmail.com.

**SAS Version and Browser**: <br>
&nbsp;  All functions were developed under SAS 9.4 using proc iml (a SAS procedure that allows matrix programming). We used <br> 
&nbsp;  Google Chrome to verify that the Demo files run clean.

**References**: <br>
&nbsp;  Kechagias S. and Pipiras V. (2020), Modeling bivariate long‐range dependence with general phase. <i>Journal of Time <br> 
&nbsp;  Series Analysis</i> 41: 268-292. <br>
&nbsp;  Helgason H., Pipiras V. and Abry P. (2011), Fast and exact synthesis of stationary multivariate Gaussian time series <br>
&nbsp;  using circulant embedding. <i>Signal Processing</i> 91(5): 1123–1133.
