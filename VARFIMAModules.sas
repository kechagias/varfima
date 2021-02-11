/*========================================================================================================
Purpose:      Define all the functions for the VARFIMA package.

Notes:        1. From a high level all the functions here serve one of two purposes. The first is to 
			     generate data with a given covariance structure, and the second is to compute a
			     Gaussian likelihood. The end of the file has a categorization of all defined 
                 functions.
              2. The methodology to synthesize the data is taken from Reference 1, and the methodo-
                 logy to compute the likelihood is taken from Reference 2.
			  3. The functions for synthesis are translated from the matlab package hermir (see 
				 http://hermir.org/) with some adhoc adaptations to deal with complex numbers.
       

Reference:    1. Data Synthesis: https://www.sciencedirect.com/science/article/abs/pii/S0165168410004019
			  2. VARFIMA model:  https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.12504

Author:       Stefanos Kechagias
Contributors  Vladas Pipiras, Rick Wicklin, Xilong Chen
Date:         August 2020 (development started in 2015)
========================================================================================================*/

/*local directory*/
libname VARFIMA 'C:\Users\statha\Desktop\VARFIMA';


proc iml;

/* The modules: embedAutoCov, embedCrossCov, circembed, InitMultivarGauss, SynthStepMultivarGauss,
are used for simulating the stationary Gaussian series from a given covariance function. They are 
translated from hermir.org. The function that actually synmthesizes data is SynthStepMultivarGauss.*/


/*----------------------------------------------------------------------------------------------*/
start embedAutoCov(r);                                    
	N = nrow(r);
	tmp = r[2:N-1];
	z = r[1] // tmp // r[N] // tmp[N-2:1]; 
return(z);
finish;
/*----------------------------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------------------------*/
start embedCrossCov(rxy,ryx);
	N = nrow(ryx);
	tmp = ryx[2:N];
	rcirc = rxy[1:N-1] // tmp[N-1:1];
return(rcirc);
finish;
/*----------------------------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------------------------*/
start circembed(R);
/*circular embedding of covariance*/
/*works only for bivariate*/

n = nrow(R);
Rcirc = repeat(0,2*n-2,3);

Rcirc[,1] = embedAutoCov(R[,1]); 
Rcirc[,2] = embedAutoCov(R[,4]); 

/*assuming assymetry*/
Rcirc[,3]=2*embedCrossCov(R[,2],R[,3]);

return(Rcirc);
finish;
/*----------------------------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------------------------*/

start InitMultivarGauss(sqrtL, crossSpec, L,    R);

Rcirc =  circembed(R);
/*works for bivariate*/
P  = 2;
M = 2*nrow(Rcirc)-2;

N = nrow(Rcirc)/2+1;
L = repeat(0,N,6);

L[,1:2]=fft(Rcirc[,1]);
L[,3:4]=fft(Rcirc[,2]);
L[,5:6]=fft(Rcirc[,3]);

L =   L //  (L[N-1:2,1:5] || -L[N-1:2,6]);

sqrtL = repeat(0,2*N-2,2);
/* check if circularized-imaginary part has to be 0*/
do k = 1 to P;
		if norm(L[,2*k])>0.001 then do;
			print 'Imaginary part of FFT seems non-zero.';
		end;

		if sum(L[,2*k-1]<0)>0 then do;
			print 'FFT has to be positive.';
		end;
    sqrtL[,k]=sqrt(L[,2*k-1]);
end;


LxyR = L[,5];
LxyI = L[,6];
tmp = sqrtL[,1]#sqrtL[,2];
nonzeroInd = loc(tmp);
crossSpec = LxyR[nonzeroInd]/tmp[nonzeroInd] || LxyI[nonzeroInd]/tmp[nonzeroInd];
finish;
/*----------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------------*/
start SynthStepMultivarGauss(X, Xiid,    R);

run  InitMultivarGauss(sqrtL, crossSpec, L,    R);

CmR = repeat(0,2,2);
CmI = repeat(0,2,2);
Z1 = repeat(0,nrow(sqrtL),2);
Z2 = repeat(0,nrow(sqrtL),2);

/* The code here is a bit different than matlab bc SAS dosent support 
complex matrices. Therefore, to find the Cholesky decomposition of 
a complex matrix A, I need to use the Cholesky decompositions
of real(A) and imag(A)*/

do m = 1 to nrow(sqrtL);
	CupperR = crossSpec[m,1];
	CupperI = crossSpec[m,2];

	/* real and imaginary parts of Cm */ 
    CmR = (2 || CupperR) // (CupperR || 2);
    CmI = (2 || CupperI) // (CupperI || 2);
	/*find Cholesky decomposition of Cm */

	C1 = root(CmR);
    C2 = root(CmI);

	l1R = C1[1,1] || C1[1,2];
	l2R = {0} || sqrt(C1[1,1]**2-C1[1,2]**2-C2[1,2]**2);


	AmR = (l1R // l2R)`;
	AmI = {0 0} // (C2[1,2] || 0);


/*    call randseed(1);*/
	A  = randnormal(2, {0 0}, {1 0,0 1})/sqrt(2);
    WR = A[1,]`;
    WI = A[2,]`;
/*	WR = {0.2,-0.9};*/       /*this is super useful to compare with matlab code!*/
/*	WI = {-0.3,0.1};*/
	tmpWR = AmR[1,1]*WR[1,1] // AmR[2,1]*WR[1,1]-AmI[2,1]*WI[1,1]+AmR[2,2]*WR[2,1] ;
	tmpWI = AmR[1,1]*WI[1,1] // AmI[2,1]*WR[1,1]+AmR[2,1]*WI[1,1]+AmR[2,2]*WI[2,1];


Z1[m,] = tmpWR[1,1] || tmpWI[1,1] ;
Z2[m,] = tmpWR[2,1] || tmpWI[2,1] ;

end;

tmp1 = sqrtL[,1]#Z1;
tmp2 = sqrtL[,2]#Z2;

/*Note: Since tmp1 and tmp2 are complex valued, we have to express their fft through the ffts
of their real and imaginary part. More specifically we have the following:

Real(fft(tmp1)) = Real(fft(Real(tmp1))+Imag(fft(Real(tmp1)))   (1)

Also note that in SAS IML fft is defined
differently than in matlab-the complex exponential has no minus*/


R1 = fft(tmp1[,1])[,1]+fft(tmp1[,2])[,2]; /*corresponds to (1) above*/
I1 = -fft(tmp1[,1])[,2]+fft(tmp1[,2])[,1];
R2 = fft(tmp2[,1])[,1]+fft(tmp2[,2])[,2];
I2 = -fft(tmp2[,1])[,2]+fft(tmp2[,2])[,1];
F1 = R1 || I1;
F2 = R2 || I2;


Xtilde = 1/sqrt(nrow(sqrtL))*F1 || 1/sqrt(nrow(sqrtL))*F2;
X = Xtilde[,1] || Xtilde[,3];

Xiid = Xtilde[,2] || Xtilde[,4];

finish;
/*-------------------------------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------*/
start conv(u,v);
/* w = conv(u,v) convolves vectors u and v.
 * Algebraically, convolution is the same operation as
 * multiplying the polynomials whose coefficients are the
 * elements of u and v. Straight convolution is too slow,
 * so use the FFT.
 *
 * Both of u and v are column vectors.
 */
   m = nrow(u);
   n = nrow(v);

   wn = m + n - 1;
   /* find p so that 2##(p-1) < wn <= 2##p */
   p = ceil( log(wn)/ log(2) );
   nice = 2##p;

   a = fft( u // j(nice-m,1,0) );
   b = fft( v // j(nice-n,1,0) );
   /* complex multiplication of a and b */
   wReal = a[,1]#b[,1] - a[,2]#b[,2];
   wImag = a[,1]#b[,2] + a[,2]#b[,1];
   w = wReal || wImag;
   z=ifft(w);
   z = z[1:wn,1] / nice;  /* take real part and first wn elements */
   return (z);
finish;
/*--------------------------------------------------------------------*/


/*--------------------------------------*/
start PowerOfComplex(x,n);
/*
Compute nth power of a complex number
 Input
	x    1x2 vector, real and imaginary part 
	n    nth power
 Output
    y    1x2 vector, real and imag part of 
         x  raised to n*/

/*polar coordinates*/
r = sqrt(x[1]##2+ x[2]##2);
phi = atan2(x[2],x[1]);

y = (r##n#cos(n*phi) || r##n#sin(n*phi));
return(y);
finish ;
/*--------------------------------------*/

/*--------------------------------------*/
start cplxMult(A, B);
/*Author: Rick Wicklin from SAS communities*/
/* Complex multiplication of A*B. 
   A vector of complex numbers is a two-column 
   matrix where Re(A)=A[,1] and Im(A)=A[,2]. 
   If A = x + iy and B = u + iv then 
   C = A*B = (x#u - y#v) + i(x#v + y#u)*/
   C = j(nrow(B),2);
   C[,1] = A[,1]#B[,1] - A[,2]#B[,2];
   C[,2] = A[,1]#B[,2] + A[,2]#B[,1];
   return( C );
finish;
/*--------------------------------------*/


/*----------------------------------*/
start ComplexInv(W);
pi = constant('PI');
/*Consider a 2x2 matrix of complex-valued 
eigenvectors (one is conjugate of the other)
The matrix is W = [w1 w1';w2 w2'] where prime 
is conjugate operator. This function computes inv(W)
using the formula inv(A) = transpose(adjoint(A))/detA)

Input
	first column of W   The second is conjugate
Output
	first row of inv(W) The second rrow will be the conjugate 
*/

/*work with polar coordinates */
r1 = sqrt(W[1,1]##2+ W[1,2]##2);
f1 = atan2(W[1,2],W[1,1]);

r2 = sqrt(W[2,1]##2+ W[2,2]##2);
f2 = atan2(W[2,2],W[2,1]);

s2 = 0.5/(r2*sin(f1-f2));
s1 = 0.5/(r1*sin(f1-f2));
invW = j(2,2,.);
invW[1,] = cos(-pi/2-f2)*s1 || sin(-pi/2-f2)*s1;
invW[2,] = -cos(-pi/2-f1)*s2 || -sin(-pi/2-f1)*s2;

return(invW);
finish;
/*----------------------------------*/


/*--------------------------------------------------------------------*/
start ApplyInvArConv(Phi1,Y);
/* Phi is a 2x2 matrix and Y is a bivariate series are given. 
This module calculates X such that Phi(B)X = Y. */

M = nrow(Y);
/*eigendecomposition of Phi to compute its powers faster*/
call eigen(L,W, Phi1);

k = (0:M-1)`; /*can probably replace with 0:min(M,200)-1*/

/*two cases: real or complex eigenvalues eigenvalues for Phi*/
if (ncol(L)>1) then do;
	if (L[1,2]=0)  then do; /*real eigenvalues*/
		Ytilde = inv(W)*Y`;
		C = conv((L[1,1]**k),Ytilde[1,]`)` // conv((L[2,1]**k),Ytilde[2,]`)`;
		X = (W*C)`;
		X = X[1:M,];
	end;
	else do;
		/*write data in real and complex parts*/
		Y1 = Y[,1] || j(M,1,0);
		Y2 = Y[,2] || j(M,1,0);

		/*compute inverse of W*/
		invW =  ComplexInv(W);
	
		/*compute inv(W)*Y*/
		Ytilde = (cplxMult(invW[1,], Y1)+ cplxMult(invW[2,], Y2));
	
		/*compute the powers of the eigenvalues*/
		Lk = PowerOfComplex(L,k);

		/*compute product L^k*inv(W)*Y */
		tmp = conv(Lk[,1],Ytilde[,1])- conv(Lk[,2],Ytilde[,2]) || conv(Lk[,1],Ytilde[,2]) + conv(Lk[,2],Ytilde[,1]);
		tmpConj = tmp[,1] || -tmp[,2];
		WConj = W[,1] || -W[,2];

		/*compute final product*/
		a1 = cplxMult(W[1,],tmp[1:M,]) + cplxMult(WConj[1,],tmpConj[1:M,]);
		a2 = cplxMult(W[2,],tmp[1:M,]) + cplxMult(WConj[2,],tmpConj[1:M,]);
		X = a1[,1] || a2[,1];
	end;
end;
else do;
	/*when eigenvalue are real, sometimes the eigenvalues magtrix L has a column of zero 
	imaginary parts, other times there is only the real part */
	Ytilde = inv(W)*Y`;
	C = conv((L[1]**k),Ytilde[1,]`)` // conv((L[2]**k),Ytilde[2,]`)`;
	X = (W*C)`;
	X = X[1:M,];
end;

return(X);
finish;
/*--------------------------------------------------------------------*/



/*-------------------------------------------------------------------------------------------*/
start CovarVARFIMA0Dq(M,d,c,Sigmae,Theta,  GG, GGNeg);

/* Covariance of identifiable general phase bivariate VARFIMA(0,D,q)
Inputs:
	M		 length of the covariance sequence to generate
	d   	 vector of 2 LRD parameters
	c        1-dim parameter related to phase--c-0 for 1-sided VARFIMA
	Theta    coefficient 2x2 matrices in MA part

Output 
	GG GGNeg Mx4 matrices containing the entries of the covariance matrices
			 at positive and negative lags
          

For example, The first row of GG will have the (1,1), (1,2), (2,1) and
(2,2) entries of the 2x2 covariance matrix Gamma(0). The second row
will have the entries of G(1) and so on. 
*/


P = ncol(d); /* dimension of the series*/
q = ncol(Theta)/P-1; /*MA order*/

/* ensure elements in d are inside principal range */
ind = (d<=-1/2 | d>=1/2);
if sum(ind[1,])>0 then do;
    print 'd does not satisfy -1/2<d(k)<1/2 for all k';
end;

pi = constant('PI');

/*allocate memory*/
GG = repeat(0,M,P*P);
GGNeg = repeat(0,M,P*P);
R = repeat(0,P,P);


/*For every lag, that is for every n, I need to calculate a 2x2 covariance matrix*/
/*j and k in the loop below go over the 4 elements of the matrix*/

/**The following is implementing relations (3.3)-(3.5) in  
http://stefanos.web.unc.edu/files/2018/12/Modeling-bivariate-long-range-dependence-with-general-phase.pdf*/

do n = 0 to M-1;
	do j = 1 to P;
		do k = 1 to P;
		/*see relation (3.4)*/
			a1 = c##2*(-1)##((j+k));
			a4 = c*(-1)##(k+1);	
			a2 = c*(-1)##(j+1);
     
			/*sums in (3.3)*/
			Gsum1 = 0;
			Gsum2 = 0;
			Gsum3 = 0;
			Gsum4 = 0;

			/*numerators in (3.5) that depend on on j and k*/
			num1 = 2*exp(lgamma(1-d[j]-d[k]))*sin(d[k]*pi);
			num3 = 2*exp(lgamma(1-d[k]-d[j]))*sin(d[j]*pi);
			if (abs(d[j]+d[k])>10##(-10)) then do;
				num4 = 2*pi/gamma(d[j]+d[k]);
				num2 = 2*pi/gamma(d[k]+d[j]); /*gamma_jk(n) = gamma_kj(-n)*/
			end;
			else do;
				num4 = 0;
				num2 = 0;
			end;

			/*need to break the loop into cases--If I have large n I ll have to use lgamma function
			so fractions dont explode. When n+t-s = 0 , however, the case d = 0 cannot be handled by lgamma*/
			do s = 0 to q;		
				do t = 0 to s-n-1;  	
					if ((abs(d[k])<10##-8 && 1+n+t-s<=0) || (abs(d[j])<10##-8 && n+t-s<=0)) then do;
						G1 = 0;
						G3 = 0;
					end;
					else do;
						G1 = num1*gamma(n+t-s+d[k])/gamma(n+t-s+1-d[j]);
						G3 = num3*gamma(n+t-s+d[j])/gamma(n+t-s+1-d[k]);
					end;
					do u = 1 to P;	/*P is the dimension so it would be 2*/												   		   			/**/   	    /**/
						do v = 1 to P;	
						Gsum1 = Gsum1 + Theta[j,2*s+u]*Theta[k,2*t+v]*a1*Sigmae[u,v]*G1;
						Gsum3 = Gsum3 + Theta[j,2*s+u]*Theta[k,2*t+v]*Sigmae[u,v]*G3;
						end;
					end;
				end;

				do t = max(s-n,0) to q;
					if (n+t-s=0) then do;
						G1 = 2*gamma(1-d[j]-d[k])*pi/(gamma(1-d[j])*gamma(1-d[k]));
						G3 = G1;
						G4 = 2*pi;
					end;
					else do;
						tmp1 = lgamma(n+t-s+1-d[j]) - lgamma(n+t-s+d[k]);
						denom1 = exp(tmp1);
						G1 = num1/denom1;

						tmp3 = lgamma(n+t-s+1-d[k]) - lgamma(n+t-s+d[j]);
						denom3 = exp(tmp3);
						G3 = num3/denom3;

						tmp4 = lgamma(n+1+t-s) - lgamma(n+t-s+d[j]+d[k]);
						denom4 = exp(tmp4);
						G4 = num4/denom4;
					end;
					do u = 1 to P;													   		   			
						do v = 1 to P;
							Gsum1 = Gsum1 + Theta[j,2*s+u]*Theta[k,2*t+v]*a1*Sigmae[u,v]*G1;
							Gsum3 = Gsum3 + Theta[j,2*s+u]*Theta[k,2*t+v]*Sigmae[u,v]*G3;
							Gsum4 = Gsum4 + Theta[j,2*s+u]*Theta[k,2*t+v]*a4*Sigmae[u,v]*G4;
						end;
					end;
				end;
				do t = 0 to s-n;
					if (n+t-s=0) then do; 
						G2 = 2*pi;
					end;
					else do;
						G2 = 2*pi*exp(lgamma(d[j]+d[k]-n+s-t))/(exp(lgamma(d[j]+d[k]))*exp(lgamma(1-n+s-t)));
					end;
					do u = 1 to P;													   		   			
						do v = 1 to P;
							Gsum2 = Gsum2 + Theta[j,2*s+u]*Theta[k,2*t+v]*a2*Sigmae[u,v]*G2;
						end;
					end;
				end;
			end;

			/*final 2x2 covariance matrix in relation (3.3) for current value of n (lag)*/
			R[j,k] = (1/(2*pi))*(Gsum1 + Gsum2 + Gsum3 + Gsum4);

			/*save the current covariance matrix*/
			GG[n+1, P*(j-1)+k] = R[j,k];
			GGNeg[n+1, P*(k-1)+j] = R[j,k];
		end;		      		  				              												
	end;				  		        		              												
end;

finish;
/*---------------------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------------------*/
start ComputeGammas(M,q,d, GMinus,LGMinus, GPlus, LGPlus);
/*
PURPOSE:
This function computes gamma functions of  the form gamma(k+d),
where k is an integer ranging from -q to M+q and d is a real in (-0.5 0.5).

Notes            M typically will range from 1 to 1000 and q = 0,1,2,3.

INPUT 
	M 			 sample size
	q			 MA order
	d  			 2-dimensional parameter
   j,k	         indices of d
OUTPUT
	GMinus   	 gamma(k-d) when k-d<0
	GPlus    	 gamma(k+d) when k-d<0
	LGMinus      log(gamma(k-d)) when k-d>0
	LGPlus       log(gamma(k+d)) when k-d>0*/
	
/*range of k*/
T = M+q+1;
indexRange = -q:T;

/*I will split indices according to the sign of k+d and k-d for both d[1] and d[2]--
So overall I have 8 categories*/

T = (M+2*q+2);
IndexMinus = j(1,T,.);
AllGammasMinus = j(1,T,0);
AllGammasPlus = j(1,T,0);

IndexPlus= indexRange + d;
IndexMinus = indexRange - d;


/*CASE 1: For k-d negative arguments*/
NegIndexMinus = loc(IndexMinus<1); /*locations of entries with negative arguments*/

/*provided there is at least one such location compute gammas*/
if mod(IndexMinus,1) then 
GMinus =  gamma(IndexMinus[NegIndexMinus]);

/*CASE 2:For k-d positive arguments*/
PosIndexMinus = loc(IndexMinus>0);
z = sum(IndexMinus<=0);
if (z) then do;
	LGMinus= j(z,1,.) // lgamma(IndexMinus[PosIndexMinus]);
end;
else do;
	LGMinus= lgamma(IndexMinus[PosIndexMinus]);
end;

/*CASE 3:For k+d negative arguments*/
NegIndexPlus = loc(IndexPlus<1);

if mod(IndexPlus,1) then 
GPlus =    gamma(IndexPlus[NegIndexPlus]);

/*CASE 4:For k+d positive arguments*/
PosIndexPlus = loc(IndexPlus>0);
z = sum(IndexPlus<=0);
if (z) then do;
	LGPlus =    j(z,1,.) // lgamma(IndexPlus[PosIndexPlus]);
end;
else do;
	LGPlus =     lgamma(IndexPlus[PosIndexPlus]);
end;

finish;
/*--------------------------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------------*/
start LargeCovar0DqVec(M,d,c,Sigmae,Theta,  trows, tcols);

/* Purpose: Compute the KM*KM variance-covariance matrix of a*/
/* VARFIMA0Dq series. If X is a bivariate series of size 100, 
then K=2, M =100.

Let G(h) = E X_0 X_h be the 2x2 autocovariance function at lag h. 
(zero mean)

INPUTS
	 M         Sample size

	 Parameters
	 d         2-dim LRD (long-range dependence) parameters
	 c 	       1-dim parameter for two sided model
	 Sigmae    2x2 Innovation variance matrix
	 Theta     list of MA matrices

OUTPUT
	trows, 
    tcols      first rows and first columns of the 4 Toeplitz matrices that build the 
			   large coavariance matrix of a bivariate time series

In more detail:
There are two ways to put the elements of the M many 2x2 matrixes G(h) together and form
one large covariance matrix (M=T):

CASE 1:
 Omega   G(0) G(-1)...G(-T+1)
         G(1) G(0)... G(-T+2)
         ....................
         ....................
         G(T-1)..........G(0)

CASE2
 OmBlock  [G11 G12
           G21 G22]

/* where G11 is toeplitz with first row the (1,1) elements*/
/* of the matrices G(0) G(-1), G12 is toeplitz with first row*/
/* the (1,2) elements of the matrices G(0) G(-1), etc.*/

/* OUTPUT*/
/* trows, tcols 1st rows and 1st columns of matrices Gij. Since these are toeplitz matrices
all I need to know to have OmBlock is the first rows and first columns of G11, G12, G21, G22*/

/* Stefanos Kechagias--August 20 2017.*/

P = ncol(d); /* dimensionality */
q = ncol(Theta)/P-1;

/* ensure elements in d are inside principal range */
ind = (d<=-1/2 | d>=1/2);
if sum(ind[1,])>0 then do;
    print 'd does not satisfy -1/2<d(k)<1/2 for all k';
end;

pi = constant('PI');

/* calculate R */
trows = j(P,P*M,.);
tcols = j(P*M,P,.);

R = repeat(0,P,P);

/*will use these when computing the double sums from u to P and from v to P
see covariance expressions in  relations (3.3)-(3.5) in  
http://stefanos.web.unc.edu/files/2018/12/Modeling-bivariate-long-range-dependence-with-general-phase.pdf*/
u = 1:2;
v = 1:2;
/*For every lag, that is for every n, I need to calculate a 2x2 covariance matrix*/
/*j and k in the loop below go over the 4 elements of the matrix*/

/**The following is implementing relations (3.3)-(3.5) in  
http://stefanos.web.unc.edu/files/2018/12/Modeling-bivariate-long-range-dependence-with-general-phase.pdf*/

/*First I will compute and store all gamma and log gamma functions that enter these formulas--and contain n*/

call ComputeGammas(M,q,0, GMinus,LGMinus, GPlus, LGPlus); 
call ComputeGammas(M,q,d[1], GMinus_1,LGMinus_1, GPlus_1, LGPlus_1);
call ComputeGammas(M,q,d[2], GMinus_2,LGMinus_2, GPlus_2, LGPlus_2);

call ComputeGammas(M,q,d[1]+d[1], GMinus_11,LGMinus_11, GPlus_11, LGPlus_11);
call ComputeGammas(M,q,d[1]+d[2], GMinus_12,LGMinus_12, GPlus_12, LGPlus_12);
call ComputeGammas(M,q,d[2]+d[2], GMinus_22,LGMinus_22, GPlus_22, LGPlus_22);

GMinusRows = j(1,5,.);
if (GMinus_1) then GMinusRows[1] = nrow(GMinus_1);
if (GMinus_2) then GMinusRows[2] = nrow(GMinus_2);
if (GMinus_11) then GMinusRows[3] = nrow(GMinus_11);
if (GMinus_12) then GMinusRows[4] = nrow(GMinus_12);
if (GMinus_22) then GMinusRows[5] = nrow(GMinus_22);

GPlusRows = j(1,5,.);
if (GPlus_1) then GPlusRows[1] = nrow(GPlus_1);
if (GPlus_2) then GPlusRows[2] = nrow(GPlus_2);
if (GPlus_11) then GPlusRows[3] = nrow(GPlus_11);
if (GPlus_12) then GPlusRows[4] = nrow(GPlus_12);
if (GMPlus_22) then GPlusRows[5] = nrow(GPlus_22);


GMRows = max(GMinusRows);
GPRows = max(GPlusRows);

GMinus_All = j(GMRows,5,.);

if (GMinus_1) then GMinus_All[1:nrow(GMinus_1),1] = GMinus_1;
if (GMinus_2) then GMinus_All[1:nrow(GMinus_2),2] = GMinus_2;
if (GMinus_11) then GMinus_All[1:nrow(GMinus_11),3] = GMinus_11;
if (GMinus_12) then GMinus_All[1:nrow(GMinus_12),4] = GMinus_12;
if (GMinus_22) then GMinus_All[1:nrow(GMinus_22),5] = GMinus_22;

GPlus_All = j(GPRows,5,.);
if (GPlus_1) then GPlus_All[1:nrow(GPlus_1),1] = GPlus_1;
if (GPlus_2) then GPlus_All[1:nrow(GPlus_2),2] = GPlus_2;
if (GPlus_11) then GPlus_All[1:nrow(GPlus_11),3] = GPlus_11;
if (GPlus_12) then GPlus_All[1:nrow(GPlus_12),4] = GPlus_12;
if (GPlus_22) then GPlus_All[1:nrow(GPlus_22),5] = GPlus_22;

LGMinus_All = LGMinus_1 || LGMinus_2 || LGMinus_11 || LGMinus_12 || LGMinus_22;
LGPlus_All = LGPlus_1 || LGPlus_2 || LGPlus_11|| LGPlus_12|| LGPlus_22;

TGp = nrow(GPlus_All);
TGm = nrow(GMinus_All);
TG = min(TGp,TGm);

TLGp = nrow(LGPlus_All);
TLGm = nrow(LGMinus_All);
TLG = min(TLGp,TLGm);

u = 1:2;
v = 1:2;

ThetaSigmaTheta = j(4,4,.);


do k = 1 to P; /*P=2*/
	Theta_k = Theta[k,];
	do j = 1 to P;
		Theta_j = Theta[j,];  
		/*Compute and save Theta_j[, 2*s+u] * Sigmae*Theta_k[,2*t+v]` */
		do s = 0 to q;
			do t = 0 to q;
				ThetaSigmaTheta[2*s+t+1,2*k+j-2] = Theta_j[, 2*s+u] * Sigmae*Theta_k[,2*t+v]`; 
			end;
		end;

		/*Compute the four gamma ratios*/
		num1 = 2*exp(lgamma(1-d[j]-d[k]))*sin(d[k]*pi);
		if (TG>1) then G1Ratio = num1*GPlus_All[1:TG-1,k]/GMinus_All[2:TG,j];
		LG1Ratio = num1/exp(LGMinus_All[2:TLG,j] - LGPlus_All[1:TLG-1,k]); 

		num3 = 2*exp(lgamma(1-d[k]-d[j]))*sin(d[j]*pi);
		if (TG>1) then G3Ratio = num3*GPlus_All[1:TG-1,j]/GMinus_All[2:TG,k];
		LG3Ratio = num3/exp(LGMinus_All[2:TLG,k] - LGPlus_All[1:TLG-1,j]);
		if (abs(d[j]+d[k])>10##(-10)) then do;
			num4 = 2*pi/gamma(d[j]+d[k]);
			num2 = 2*pi/gamma(d[k]+d[j]); /*gamma_jk(n) = gamma_kj(-n)*/
		end;
		else do;
			num4 = 0;
			num2 = 0;
		end;
		
		LG4Ratio = num4/exp(LGMinus[2:M+2*q+2,] - LGPlus_All[1:TLG-1,k+j+1]);
		LG2Ratio = num2*exp(LGPlus_All[1:TLG-1,j+k+1] - LGMinus[2:M+2*q+2,]);


			a1 = c##2*(-1)##((j+k));
			a4 = c*(-1)##(k+1);	
			a2 = c*(-1)##(j+1);

			G1Small = 2*gamma(1-d[j]-d[k])*pi/(gamma(1-d[j])*gamma(1-d[k]));

		do n = 0 to M-1; /*sample size--ranging from 50-1000*/

			Gsum1 = 0;
			Gsum2 = 0;
			Gsum3 = 0;
			Gsum4 = 0;

		   /*----------------------------------------------------------------------------------------------------------*/
			do s = 0 to q;	/*q is the MA order-so it would be 0 or 1*/	
				do t = 0 to s-n-1;  	
					if ((abs(d[k])<10##-8 & 1+n+t-s<=0) | (abs(d[j])<10##-8 & n+t-s<=0)) then do;
						G1 = 0;
						G3 = 0;
					end;
					else do;
						index = n+t-s+q+1;
						G1 = G1Ratio[index];
						G3 = G3Ratio[index];
					end;
					tmp = ThetaSigmaTheta[2*s+t+1,2*k+j-2];
					Gsum1 = GSum1 + a1*G1*tmp;
					Gsum3 = GSum3 +    G3*tmp;
				end;
				
				do t = max(s-n,0) to q;
					if (n+t-s=0) then do;
						G1 = G1Small;
						G3 = G1;
						G4 = 2*pi;
					end;
					else do;
						index = n+t-s+q+1;
						G1 = LG1Ratio[index];
						G3 = LG3Ratio[index];
						G4 = LG4Ratio[index];
					end;
					tmp = ThetaSigmaTheta[2*s+t+1,2*k+j-2];
					Gsum1 = GSum1 + a1*G1*tmp;
					Gsum3 = GSum3 +    G3*tmp;
					Gsum4 = GSum4 + a4*G4*tmp;
				end;
				do t = 0 to s-n;
					if (n+t-s=0) then do; 
						G2 = 2*pi;
					end;
					else do;
						if (abs(d[j]+d[k])<10##-8) then do;
							G2 = 0;
						end;
						else do;
							index = -n+s-t+q+1;
							G2 = LG2Ratio[index];
						end;
					end;
					tmp = ThetaSigmaTheta[2*s+t+1,2*k+j-2]; 
					Gsum2 = GSum2 + a2*G2*tmp;
				end;
			end;
			R[j,k] = (1/(2*pi))*(Gsum1 + Gsum2 + Gsum3 + Gsum4);
			trows[k,(j-1)*M+n+1] = R[j,k];
			tcols[(j-1)*M+n+1,k] = R[j,k];
		end;		      		  				              												
	end;				  		        		              												
end;
finish;
/*---------------------------------------------------------------------------------------------*/



/*-------------------------------------------------------------------------------------------*/
start LargeCovar0Dq(M,d,c,Sigmae,Theta,  trows, tcols);
/*NOTE: This fucntion is an older and slower version of LargeCovar0DqVec. I leave it here 
as it is a bit easier to read it and udnerstand the code. (I dont precopmpute the Gamma functions)*/
/* Purpose: Compute the KM*KM variance-covariance matrix of a*/
/* VARFIMA0Dq series. If X is a bivariate series of size 100, 
then K=2, M =100.

Let G(h) = E X_0 X_h be the 2x2 autocovariance function at lag h. 
(zero mean)

INPUTS
	 M         Sample size

	 Parameters
	 d         2-dim LRD (long-range dependence) parameters
	 c 	       1-dim parameter for two sided model
	 Sigmae    2x2 Innovation variance matrix
	 Theta     list of MA matrices

OUTPUT
	trows, 
    tcols      first rows and first columns of the 4 Toeplitz matrices that build the 
			   large coavariance matrix of a bivariate time series

In more detail:
There are two ways to put the elements of the M many 2x2 matrixes G(h) together and form
one large covariance matrix (M=T):

CASE 1:
 Omega   G(0) G(-1)...G(-T+1)
         G(1) G(0)... G(-T+2)
         ....................
         ....................
         G(T-1)..........G(0)

CASE2
 OmBlock  [G11 G12
           G21 G22]

/* where G11 is toeplitz with first row the (1,1) elements*/
/* of the matrices G(0) G(-1), G12 is toeplitz with first row*/
/* the (1,2) elements of the matrices G(0) G(-1), etc.*/

/* OUTPUT*/
/* trows, tcols 1st rows and 1st columns of matrices Gij. Since these are toeplitz matrices
all I need to know to have OmBlock is the first rows and first columns of G11, G12, G21, G22*/

/* Stefanos Kechagias--August 20 2017.*/



P = ncol(d); /* dimensionality */
q = ncol(Theta)/P-1;

/* ensure elements in d are inside principal range */
ind = (d<=-1/2 | d>=1/2);
if sum(ind[1,])>0 then do;
    print 'd does not satisfy -1/2<d(k)<1/2 for all k';
end;

pi = constant('PI');


/* calculate R */
trows = j(P,P*M,.);
tcols = j(P*M,P,.);

R = repeat(0,P,P);

/*For every lag, that is for every n, I need to calculate a 2x2 covariance matrix*/
/*j and k in the loop below go over the 4 elements of the matrix*/

/**The following is implementing relations (3.3)-(3.5) in  
http://stefanos.web.unc.edu/files/2018/12/Modeling-bivariate-long-range-dependence-with-general-phase.pdf*/
do k = 1 to P; /*P=2*/
	do j = 1 to P;
		do n = 0 to M-1; /*sample size--ranging from 50-1000*/
			a1 = c##2*(-1)##((j+k));
			a4 = c*(-1)##(k+1);	
			a2 = c*(-1)##(j+1);

			Gsum1 = 0;
			Gsum2 = 0;
			Gsum3 = 0;
			Gsum4 = 0;

		   /*----------------------------------------------------------------------------------------------------------*/
			num1 = 2*exp(lgamma(1-d[j]-d[k]))*sin(d[k]*pi);
			num3 = 2*exp(lgamma(1-d[k]-d[j]))*sin(d[j]*pi);
			if (abs(d[j]+d[k])>10##(-10)) then do;
				num4 = 2*pi/gamma(d[j]+d[k]);
				num2 = 2*pi/gamma(d[k]+d[j]); /*gamma_jk(n) = gamma_kj(-n)*/
			end;
			else do;
				num4 = 0;
				num2 = 0;
			end;
			
			do s = 0 to q;	/*q is the MA order-so it would be 0 or 1*/	
				do t = 0 to s-n-1;  	
					if ((abs(d[k])<10##-8 && 1+n+t-s<=0) || (abs(d[j])<10##-8 && n+t-s<=0)) then do;
						G1 = 0;
						G3 = 0;
					end;
					else do;
						G1 = num1*gamma(n+t-s+d[k])/gamma(n+t-s+1-d[j]);
						G3 = num3*gamma(n+t-s+d[j])/gamma(n+t-s+1-d[k]);
					end;
					do u = 1 to P;	/*P is the dimension so it would be 2*/												   		   			/**/   	    /**/
						do v = 1 to P;	
						Gsum1 = Gsum1 + Theta[j,2*s+u]*Theta[k,2*t+v]*a1*Sigmae[u,v]*G1;
						Gsum3 = Gsum3 + Theta[j,2*s+u]*Theta[k,2*t+v]*Sigmae[u,v]*G3;
						end;
					end;
				end;
				
				do t = max(s-n,0) to q;
					if (n+t-s=0) then do;
						G1 = 2*gamma(1-d[j]-d[k])*pi/(gamma(1-d[j])*gamma(1-d[k]));
						G3 = G1;
						G4 = 2*pi;
					end;
					else do;
						tmp1 = lgamma(n+t-s+1-d[j]) - lgamma(n+t-s+d[k]);
						denom1 = exp(tmp1);
						G1 = num1/denom1;

						tmp3 = lgamma(n+t-s+1-d[k]) - lgamma(n+t-s+d[j]);
						denom3 = exp(tmp3);
						G3 = num3/denom3;

						tmp4 = lgamma(n+1+t-s) - lgamma(n+t-s+d[j]+d[k]);
						denom4 = exp(tmp4);
						G4 = num4/denom4;
					end;
					do u = 1 to P;													   		   			
						do v = 1 to P;
							Gsum1 = Gsum1 + Theta[j,2*s+u]*Theta[k,2*t+v]*a1*Sigmae[u,v]*G1;
							Gsum3 = Gsum3 + Theta[j,2*s+u]*Theta[k,2*t+v]*Sigmae[u,v]*G3;
							Gsum4 = Gsum4 + Theta[j,2*s+u]*Theta[k,2*t+v]*a4*Sigmae[u,v]*G4;
						end;
					end;
				end;
				do t = 0 to s-n;
					if (n+t-s=0) then do; 
						G2 = 2*pi;
					end;
					else do;
						G2 = 2*pi*exp(lgamma(d[j]+d[k]-n+s-t))/(exp(lgamma(d[j]+d[k]))*exp(lgamma(1-n+s-t)));
					end;
					do u = 1 to P;													   		   			
						do v = 1 to P;
							Gsum2 = Gsum2 + Theta[j,2*s+u]*Theta[k,2*t+v]*a2*Sigmae[u,v]*G2;
						end;
					end;
				end;
			end;
			R[j,k] = (1/(2*pi))*(Gsum1 + Gsum2 + Gsum3 + Gsum4);
			trows[k,(j-1)*M+n+1] = R[j,k];
			tcols[(j-1)*M+n+1,k] = R[j,k];
		end;		      		  				              												
	end;				  		        		              												
end;


finish;
/*---------------------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------------------*/
start EvalInvQuadForm(A, v, logdetA, QuadForm);
/* Evaluate quadratic form v`*inv(A)*v where
   A is symmetric positive definite
   Want   QuadForm = v` * inv(A) * v = v` * w  where A*w = v
   ==> (U`U)*w = v             Cholesky decomp of A
   ==>  First solve U` z = v   for z,
        then solve   U w = z   for w */
/*Input */
   U = root(A);           /* Cholesky root */
   z = trisolv(2, U, v);  /* solve linear system */
   w = trisolv(1, U, z);  
   QuadForm = v` * w;            /* dot product of vectors */
   logdetA = 2*sum(log(vecdiag(U)));
 
   /*initial version Rick Wicklin, April 2019--
   modified to also compute logdetA*/
finish;
/*---------------------------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------------------------*/
start MyToeplitz(c,r);
/*create a toeplitz matrix with first row r 
and first column c. Similar to the matlab function toeplitz.*/

if c[1]^=r[1] then do;
	print 'First element of input row and first element 
	of input col are not equal';
end;

p = nrow(r);
m = nrow(c);


ind1 = do(p,2,-1);

x = r[ind1] // c;

ind2 = do(p,1,-1);

ij = repeat((0:m-1)`,1,p) + ind2;
t = shape(x[ij],p,m);

return(t);
finish;
/*---------------------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
start ComputeOmegaBlocks(trows,tcols,   T1,T2,T3,T4);
/*
Purpose:    Compute the blocks of the large covariance matrix

Inputs  
	trows
	tcols	First row and column for block matrices T1-T4
OUTPUT:
	T1-T4 	

Author:     Stefanos Kechagias
Date:       August 2020
*/

T = ncol(trows)/2;
T1 = MyToeplitz(tcols[1:T,1],trows[1,1:T]`);
tcols[1,2]=trows[1,T+1];

T2 = MyToeplitz(tcols[1:T,2],trows[1,T+1:2*T]`);
T3 = MyToeplitz(tcols[T+1:2*T,2],trows[2,T+1:2*T]`);
trows[2,1] = tcols[T+1,1];

T4 = MyToeplitz(tcols[T+1:2*T,1],trows[2,1:T]`);
finish;
/*------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------*/
start LogLik0D1(x) global(YY);
/* PURPOSE:
	compute Gaussian log-likelihood of a VARFIMA(0,D,1) series using brute force

   INPUT:
	 x = d1,d2,c,U11,U12,U22,Th11,Th12,Th21,Th22    parameter vector, where d1,d2
     are the long memory parameters, c is a phase related parameter, U11, U12, U22
     are the cholesky factos of Sigmae (variance of errors) and Th11, Th12, Th21, Th22
     are the MA coeffcients.

     YY is a vector of the individual series stacked one under the other and enters
	 the function as a global variable.
	
   OUTPUT:
     log-likelihood value at x (without the -(T/2)*log(2pi)) 

Author:    Stefanos Kechagias
Date:      August 2020
*/

/*retrieve the sample size*/
T = nrow(YY)/2;

/*retrieve the parameters*/
d = x[1:2]`;
c = x[3];
U = (x[4] || x[5])//(0 || x[6]);
Sigmae = U`*U;
tmp = (x[7] || x[8]) // (x[9] || x[10]);
Thetaq = {1 0, 0 1} || tmp;


/*compute the covariance for the given parameters*/
run LargeCovar0DqVec(T,d,c,Sigmae,Thetaq,  trows, tcols);

/*put the elements of the covariance in a block covariance matrix*/
run ComputeOmegaBlocks(trows,tcols,   Tp1,Tp2,Tp3,Tp4);
OmBlock = (Tp1 ||Tp2 ) // (Tp4 || Tp3);

/*compute the two parts of the Gaussian log-likelihood*/
call EvalInvQuadForm(OmBlock, YY, LogLik1, LogLik2);

/*obtain the value of the objective function*/
LogLik = -0.5*LogLik1 - 0.5*LogLik2;
return(LogLik);
finish;
/*---------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------*/
start LogLik0D0(x) global(YY);
/* PURPOSE:
	compute Gaussian log-likelihood of a VARFIMA(0,D,0) series using brute force

   INPUT:
	 x = d1,d2,c,U11,U12,U22 parameter vector, where d1,d2
     are the long memory parameters, c is a phase related parameter, U11, U12, U22
     are the cholesky factos of Sigmae (variance of errors).

     YY is a vector of the individual series stacked one under the other and enters
	 the function as a global variable.
	
   OUTPUT:
     log-likelihood value at x (without the -(T/2)*log(2pi)) 

Author:    Stefanos Kechagias
Date:      August 2020
*

/*retrieve the sample size*/
T = nrow(YY)/2;

/*retrieve the parameters*/
d = x[1:2]`;
c = x[3];
U = (x[4] || x[5])//(0 || x[6]);
Sigmae = U`*U;
Thetaq = I(2);

/*compute the covariance for the given parameters*/
run LargeCovar0DqVec(T,d,c,Sigmae,Thetaq,  trows, tcols);

/*put the elements of the covariance in a block covariance matrix*/
run ComputeOmegaBlocks(trows,tcols,   Tp1,Tp2,Tp3,Tp4);
OmBlock = (Tp1 ||Tp2 ) // (Tp4 || Tp3);

/*compute the two parts of the Gaussian log-likelihood*/
call EvalInvQuadForm(OmBlock, YY, LogLik1, LogLik2);

/*obtain the value of the objective function*/
LogLik = -0.5*LogLik1 - 0.5*LogLik2;
return(LogLik);
finish;
/*---------------------------------------------------------------------------------*/


start ApplyAR(Phi,Z,Y);
/*This function will be needed when we need to have AR filter apply to the process*/

/*INPUT :
  -Phi          matrix coefficientsin AR polynomial
  -Z            given process
 
  OUTPOUT :     Y = Phi(B)Z, where B is backshift operator*/
/*works for AR(1) only at this time*/

M = nrow(Z);
k = ncol(Z);					/*dimensionality of the series*/
p = ncol(Phi)/k; 			/*AR order*/

Y = j(M,k,.);
Y[1:p,] = Z[1:p,];
Y[p+1:M,] = (Z[p+1:M,]-Z[p:M-1,]*Phi);
finish;
/*-------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
start ApplyInvArConv(Phi1,Y);
/* Phi is a 2x2 matrix and Y is a bivariate series are given. 
This module calculates X such that Phi(B)X = Y. */

M = nrow(Y);
/*eigendecomposition of Phi to compute its powers faster*/
call eigen(L,W, Phi1);

k = (0:M-1)`; /*can probably replace with 0:min(M,200)-1*/

/*two cases: real or complex eigenvalues eigenvalues for Phi*/
if (ncol(L)>1) then do;
	if (L[1,2]=0)  then do; /*real eigenvalues*/
		Ytilde = inv(W)*Y`;
		C = conv((L[1,1]**k),Ytilde[1,]`)` // conv((L[2,1]**k),Ytilde[2,]`)`;
		X = (W*C)`;
		X = X[1:M,];
	end;
	else do;
		/*write data in real and complex parts*/
		Y1 = Y[,1] || j(M,1,0);
		Y2 = Y[,2] || j(M,1,0);

		/*compute inverse of W*/
		invW =  ComplexInv(W);
	
		/*compute inv(W)*Y*/
		Ytilde = (cplxMult(invW[1,], Y1)+ cplxMult(invW[2,], Y2));
	
		/*compute the powers of the eigenvalues*/
		Lk = PowerOfComplex(L,k);

		/*compute product L^k*inv(W)*Y */
		tmp = conv(Lk[,1],Ytilde[,1])- conv(Lk[,2],Ytilde[,2]) || conv(Lk[,1],Ytilde[,2]) + conv(Lk[,2],Ytilde[,1]);
		tmpConj = tmp[,1] || -tmp[,2];
		WConj = W[,1] || -W[,2];

		/*compute final product*/
		a1 = cplxMult(W[1,],tmp[1:M,]) + cplxMult(WConj[1,],tmpConj[1:M,]);
		a2 = cplxMult(W[2,],tmp[1:M,]) + cplxMult(WConj[2,],tmpConj[1:M,]);
		X = a1[,1] || a2[,1];
	end;
end;
else do;
	/*when eigenvalue are real, sometimes the eigenvalues magtrix L has a column of zero 
	imaginary parts, other times there is only the real part */
	Ytilde = inv(W)*Y`;
	C = conv((L[1]**k),Ytilde[1,]`)` // conv((L[2]**k),Ytilde[2,]`)`;
	X = (W*C)`;
	X = X[1:M,];
end;

return(X);
finish;
/*--------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------*/
start LogLik1D0(x0) global(YY);
/* PURPOSE:
	compute Gaussian log-likelihood of a VARFIMA(1,D,0) series using brute force

   INPUT:
	 x = d1,d2,c,U11,U12,U22 parameter vector, where d1,d2
     are the long memory parameters, c is a phase related parameter, U11, U12, U22
     are the cholesky factos of Sigmae (variance of errors).

     YY is a vector of the individual series stacked one under the other and enters
	 the function as a global variable.
	
   OUTPUT:
     log-likelihood value at x (without the -(T/2)*log(2pi)) 

Author:    Stefanos Kechagias
Date:      August 2020
*

/*retrieve the sample size*/
T = nrow(YY)/2;

/*retrieve the parameters*/
d = x0[1:2]`;
c = x0[3];
U = (x0[4] || x0[5])//(0 || x0[6]);
Sigmae = U`*U;
Thetaq = I(2);
Phip = (x0[7] || x0[8]) // (x0[9] || x0[10]);
p = ncol(Phip)/2;


/*compute the covariance for the given parameters*/
run LargeCovar0DqVec(T,d,c,Sigmae,Thetaq,  trows, tcols);


/*put the elements of the covariance in a block covariance matrix*/
run ComputeOmegaBlocks(trows,tcols,   Tp1,Tp2,Tp3,Tp4);
OmBlock = (Tp1 ||Tp2 ) // (Tp4 || Tp3);
ArT1 = Tp1[2:T,2:T];
ArT2 = Tp2[2:T,2:T];
ArT3 = Tp3[2:T,2:T];
ArT4 = Tp4[2:T,2:T];

OmBlockP = (ArT1 ||ArT2 ) // (ArT4 || ArT3);

/*apply Ar filter*/
Y = YY[1:T] || YY[T+1:2*T];
run ApplyAR(Phip, Y, Z);
ZZ = Z[(p+1):T,1] // Z[(p+1):T,2];


/*compute the two parts of the Gaussian log-likelihood*/
call EvalInvQuadForm(OmBlockP, ZZ, LogLik1, LogLik2);

/*obtain the value of the objective function*/
LogLik = -0.5*LogLik1 - 0.5*LogLik2;
return(LogLik);
finish;
/*---------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------*/
start LogLik1D1(x0) global(YY);
/* PURPOSE:
	compute Gaussian log-likelihood of a VARFIMA(1,D,1) series using brute force

   INPUT:
	 x = d1,d2,c,U11,U12,U22 parameter vector, where d1,d2
     are the long memory parameters, c is a phase related parameter, U11, U12, U22
     are the cholesky factos of Sigmae (variance of errors).

     YY is a vector of the individual series stacked one under the other and enters
	 the function as a global variable.
	
   OUTPUT:
     log-likelihood value at x (without the -(T/2)*log(2pi)) 

Author:    Stefanos Kechagias
Date:      August 2020
*

/*retrieve the sample size*/
T = nrow(YY)/2;

/*retrieve the parameters*/
d = x0[1:2]`;
c = x0[3];
U = (x0[4] || x0[5])//(0 || x0[6]);
Sigmae = U`*U;

Phip = (x0[7] || x0[8]) // (x0[9] || x0[10]);
p = ncol(Phip)/2;

tmp = (x0[11] || x0[12]) // (x0[13] || x0[14]);
Thetaq = I(2)  || tmp;

/*compute the covariance for the given parameters*/
run LargeCovar0DqVec(T,d,c,Sigmae,Thetaq,  trows, tcols);

/*put the elements of the covariance in a block covariance matrix*/
run ComputeOmegaBlocks(trows,tcols,   Tp1,Tp2,Tp3,Tp4);

ArT1 = Tp1[2:T,2:T];
ArT2 = Tp2[2:T,2:T];
ArT3 = Tp3[2:T,2:T];
ArT4 = Tp4[2:T,2:T];

OmBlockP = (ArT1 ||ArT2 ) // (ArT4 || ArT3);

/*apply Ar filter*/
Y = YY[1:T] || YY[T+1:2*T];
run ApplyAR(Phip, Y, Z);
ZZ = Z[(p+1):T,1] // Z[(p+1):T,2];


/*compute the two parts of the Gaussian log-likelihood*/
call EvalInvQuadForm(OmBlockP, ZZ, LogLik1, LogLik2);

/*obtain the value of the objective function*/
LogLik = -0.5*LogLik1 - 0.5*LogLik2;
return(LogLik);
finish;
/*---------------------------------------------------------------------------------*/


/* save the modules to a SAS catalog */
reset storage = VARFIMA.VARFIMAModules;  

/*load modules that synthesize data from a given covariance structure*/
store module = embedAutoCov;
store module = embedCrossCov;
store module = circembed;
store module = InitMultivarGauss;
store module = SynthStepMultivarGauss;
store module = CovarVARFIMA0Dq;


/*load modules necessary for applying an inverse bivariate AR filter*/
store module = PowerOfComplex;
store module = cplxMult;
store module = ComplexInv;
store module = ApplyInvArConv;
store module = conv;


/*load modules related to likeliuhood computation*/
store module = ComputeGammas;
store module = LargeCovar0DqVec;
store module = MyToeplitz;
store module = ComputeOmegaBlocks;
store module = EvalInvQuadForm;
store module = ApplyAR;
store module = LogLik0D0;
store module = LogLik0D1;
store module = LogLik1D0;
store module = LogLik1D1;

quit;
