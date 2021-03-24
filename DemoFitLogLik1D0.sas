/*=======================================================================
Purpose: Demo file to fenerate and fit a VARFIMA model

Notes:   The following code compiles the VARFIMA functions and fits a
		 VARFIMA(1,D,0) model to a synthetic bivariate VARFIMA(1,D,0)
		 series. See this paper for more info:
		 https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504

Author:  Stefanos Kechagias
Date     August 2020
=======================================================================*/

/*local directory where you saved the git repo*/
%let localDirectory = C:/Users/statha/Desktop/varfima;

/*create a SAS library in the specified local directory*/
libname VARFIMA %tslit(&localDirectory);

/*compile the VARFIMA functions*/
%include %tslit(&localDirectory/VARFIMAModules.sas);


proc iml;

/*----------------    STEP  1: Compile functions   ---------------------*/

/*load modules that are needed to compute likelihood*/
reset storage=VARFIMA.VARFIMAModules;  

/*load modules to synthesize data*/
load module = embedAutoCov;
load module = embedCrossCov;
load module = circembed;
load module = InitMultivarGauss;
load module = SynthStepMultivarGauss;

/*load modules necessary for applying a bivariate AR filter*/
load module = PowerOfComplex;
load module = cplxMult;
load module = ComplexInv;
load module = ApplyInvArConv;
load module = conv;


/*load modules for log likelihood computation*/
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


/*----------------    STEP  2: Generate a bivariate VARFIMA(0,D,1) series  -----------------*/
/*Specify VARFIMA(1,D,0) parameters*/
M      = 200;
Sigmae = {3 2, 2 3};
U      = root(Sigmae);
d      = {0.2 0.4};
c 	   = 0.6;
Thetaq = I(2) ;
Phi1   = {0.5 0 , 0 -0.8};

/*True Parameter Values*/
TrueParam = d[1] || d[2] || c || U[1,1] || U[1,2] || U[2,2] || Phi1[1,1] || Phi1[1,2] || Phi1[2,1] || Phi1[2,2];

/*generate covariances and synthesize FARIMA(0,d,q) process*/
run CovarVARFIMA0Dq(M,d,c,Sigmae,Thetaq,  Rneg, R);

/*generate a bivariate VARFIMA(0,D,q)series*/
run SynthStepMultivarGauss(X, Yiid,   R);

/*apply the inverse AR filter to obtain VARFIMA(1,D,0) series*/
Y = ApplyInvArConv(Phi1,X);
YY = Y[,1] // Y[,2];


/*----------------    STEP  3: Fit a bivariate VARFIMA(0,D,1) model    -----------------*/
/*select an initial point: d1,d2,c,U11,U12,U22,Phi11,Phi12,Phi21,Phi22*/
x0 = {0.25 0.25 0 1 0 1 1 0 0 1};

/*optimization options*/
optn = j(1,11,.);
optn[2] = 3;
optn[1] = 1;
optn[6] = 3;
optn[10] = 0; /*number of nonlinear constraints*/
  
/*boundary constraints*/
con = {-0.49 -0.49 . . . . . . . .,
        0.49 0.49 . . . . . . . .};

/*run the optimization*/
call nlpqn(rc,ParamEst,"LogLik1D0",x0,optn) blc=con; 


/*print parameter estimates*/
ParamNames = {'d1', 'd2', 'c', 'U11', 'U12', 'U22', 'Phi11', 'Phi12', 'Phi21', 'Phi22'};
RowNames   = {'True Parameters','Parameter Estimates'};
EstimationResults = TrueParam // ParamEst;
print EstimationResults[rowname=RowNames colname=ParamNames];






