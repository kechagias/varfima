/*=======================================================================
PURPOSE:    Demo to evaluate bivariate VARFIMA log likelihoods.

Notes:      The following code compiles the VARFIMA functions and computes 
			likelihoods of four varfima models for a bivariate Gaussian
			white noise series.

Author:     Stefanos Kechagias
Date:       August 2020 
========================================================================*/

/*local directory where you saved the git repo*/
%let localDirectory = C:/Users/statha/Desktop/varfima;

/*create a SAS library in the specified local directory*/
libname VARFIMA %tslit(&localDirectory);

/*compile the VARFIMA functions*/
%include %tslit(&localDirectory/VARFIMAModules.sas);


proc iml;

/*catalog where the modules are located*/
reset storage = VARFIMA.VARFIMAModules;  

/*load modules that are needed to compute the log-likelihood*/
load module = MyToeplitz;
load module = ComputeOmegaBlocks;
load module = ComputeGammas;
load module = LargeCovar0DqVec;
load module = EvalInvQuadForm;
load module = ApplyAR;
load module = LogLik0D0;
load module = LogLik0D1;
load module = LogLik1D0;
load module = LogLik1D1;


/*generate a bivariate normal time series (no depedence)*/
T = 100;
Y = j(T, 2, .);         
call randgen(Y, "Normal");
YY = Y[,1] // Y[,2];

/*compute log-ikelihood value for a VARFIMA 0D0 model at a point x0*/
x0 = {0.2 0.4 0.5 1 0 1};
loglik1 = LogLik0D0(x0);
print loglik1;

/*compute log-ikelihood value for a VARFIMA 0D1 model at a point x0*/
x0 = {0.2 0.4 0.5 1 0 1 0.7 -0.2 0 0.5};
loglik2 = LogLik0D1(x0);
print loglik2;

/*compute log-ikelihood value for a VARFIMA 1D0 model at a point x0*/
x0 = {0.2 0.4 0.5 1 0.3 1 0.7 -0.2 0 0.5};
loglik3 = LogLik1D0(x0);
print loglik3;


/*compute log-ikelihood value at x0 - 1D1 model*/
x0 = {0.2 0.4 0.5 1 0.3 1 0.7 -0.2 0 0.5 0.7 -0.2 0 0.5};
loglik4 = LogLik1D1(x0);
print loglik4;
