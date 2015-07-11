// Analysis of Visual Array Data
// Richard D. Morey (richarddmorey@gmail.com)

#include <R.h>
#include <Rmath.h>  
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include <Rversion.h>
#include <Rconfig.h>
#include <R_ext/Constants.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>



SEXP WM2_GibbsSampler(SEXP iters, SEXP starteffects, SEXP nCovMatR, SEXP obsCovMatR, SEXP sizeCovMatR, SEXP parStartR, SEXP covEffSlopeR, SEXP nHit, SEXP nMiss, SEXP nFA,
		      SEXP nCR,  SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, 
		      SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, 
		      SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, 
		      SEXP useA, SEXP Ktype, SEXP epsLoR, SEXP epsRngR, SEXP LFstepsR, SEXP compWeightsR, SEXP dfWish, SEXP progress, SEXP pBar, SEXP rho,
              SEXP storePred, SEXP metrop, SEXP metropSD, SEXP metropThin);

SEXP WM2_GibbsSamplerNoCov(SEXP iters, SEXP starteffects, SEXP nHit, SEXP nMiss, SEXP nFA,
		      SEXP nCR,  SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, 
		      SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, 
		      SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, 
			   SEXP useA, SEXP Ktype, SEXP epsLoR, SEXP epsRngR, SEXP LFstepsR, SEXP compWeightsR, SEXP progress, SEXP pBar, SEXP rho, SEXP storePred, SEXP metrop, SEXP metropSD, SEXP metropThin);

double WM2_quadform(double *x, double *A, int N, int incx);
double *std_rWishart_factor(double df, int p, int upper, double ans[]);
SEXP WM2_rWishartR(SEXP dfp, SEXP scal, SEXP Inv);
void WM2_rWishartC(double df, double *scal, int p, int Inv);
SEXP WM2_rmvGaussianR(SEXP mu, SEXP Sigma);
void WM2_rmvGaussianC(double *mu, double *Sigma, int p);

void WM2_hybridMC(double *ystart, int N, double epsLo, double epsRng, int LFsteps, double *compWeights, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMats, int *obsCovMats, int *sizeCovMats, int *parStart, double *means, double **pointerCovMat, double *logLike, int maxSizeCovMat, double *predVals, double *logPostCurrEff);


void WM2_metrop(double *ystart, int N, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMats, int *obsCovMats, int *sizeCovMats, int *parStart, double *means, double **pointerCovMat, double *logLike, int maxSizeCovMat, double *predVals, double *metropSD, double *logPostCurrEff);

SEXP getListElement(SEXP list, const char *str);
SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface);
int InvMatrixUpper(double *A, int p);


double LogPosterior(double *params, int p, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMat, int *obsCovMat, int *sizeCovMat, int *parStart, double *means, double **pointerCovMat, double *logLike, int maxSizeCovMat, double *predVals);
SEXP RLogPosteriorNoCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype);

SEXP RLogPosteriorWithCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype, SEXP covList, SEXP means, SEXP obsCovMat, SEXP sizeCovMat, SEXP parStart, SEXP covEffSlope);


double LogLikelihood(double *params, int p, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *Aeffects, int *Geffects, int useA, int Ktype, double *predVals);

SEXP RLogLikelihood(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP Aeffects, SEXP Geffects, SEXP useA, SEXP Ktype);

void gradLogPosterior(double *params, double *grad, int p, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMat, int *obsCovMat, int *sizeCovMat, int *parStart, double *means, double **pointerCovMat, int maxSizeCovMat);

SEXP RgradLogPosteriorNoCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype);

SEXP RgradLogPosteriorWithCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype, SEXP covList, SEXP means, SEXP obsCovMat, SEXP sizeCovMat, SEXP parStart, SEXP covEffSlope);

SEXP RPredictedProbabilities(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP Aeffects, SEXP Geffects, SEXP useA, SEXP Ktype);


#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static R_INLINE double*
internal_symmetrize(double *a, int nc)
{
    int i,j;
    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}


double WM2_quadform(double *x, double *A, int N, int incx)
{
  
  int Nsqr = N*N,info,i=0;
  double *B = Calloc(Nsqr,double);
  //double one=1;
  //double zero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;

  for(i=0;i<N;i++){
    y[i] = x[i*incx];
    //printf("%d %f\n",i,y[i]);
  }
  Memcpy(B,A,Nsqr);
  
  F77_NAME(dpotrf)("U", &N, B, &N, &info);
  F77_NAME(dtrmv)("U","N","N", &N, B, &N, y, &iOne);
  
  for(i=0;i<N;i++){
    sumSq += y[i]*y[i];
  }
  
  Free(B);
  
  return(sumSq);
}

// Shamelessly ripped from the source of the lme4 package

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and df degrees of freedom.
 *
 * @param df degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
double *std_rWishart_factor(double df, int p, int upper, double ans[])
{
    int i, j, pp1 = p + 1;

    if (df < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");
    for (j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(df - (double) j));
	for (i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = norm_rand();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

/**
 * Simulate a sample of random matrices from a Wishart distribution
 *
 * @param ns Number of samples to generate
 * @param dfp Degrees of freedom
 * @param scal Positive-definite scale matrix
 *
 * @return
 */
SEXP WM2_rWishartR(SEXP dfp, SEXP scal, SEXP Inv)
{

    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), p, psqr, iInv = INTEGER_VALUE(Inv);
    double df = asReal(dfp); // one = 1, zero = 0;
    SEXP scCp;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("rWishart: scal must be a square, real matrix");
    p = dims[0];
    psqr = p * p;
    
    PROTECT(scCp = allocMatrix(REALSXP, p, p));
    Memcpy(REAL(scCp), REAL(scal), psqr);
    
	GetRNGstate();
    WM2_rWishartC(df, REAL(scCp), p, iInv);
	PutRNGstate();
	
    UNPROTECT(1); 
    return scCp;
}

void WM2_rWishartC(double df, double *scal, int p, int Inv)
{

    int j, psqr;
    double *scCp, *tmp, one = 1, zero = 0;
    
    if(df<p)
      error("df lower than dimension in Wishart rng.");

    psqr = p * p;
    tmp = Calloc(psqr, double);
    
    AZERO(tmp, psqr);
    scCp = Memcpy(Calloc(psqr, double), scal, psqr);
    F77_NAME(dpotrf)("U", &p, scCp, &p, &j);
    if (j)
	error("scal matrix is not positive-definite.");

    //GetRNGstate();	
    
    std_rWishart_factor(df, p, 1, tmp);
    
    F77_NAME(dtrmm)("R", "U", "N", "N", &p, &p,
			&one, scCp, &p, tmp, &p);
    F77_NAME(dsyrk)("U", "T", &p, &p,
			&one, tmp, &p,
			&zero, scCp, &p);
    
    if(Inv){
      InvMatrixUpper(scCp,p);
    }
    internal_symmetrize(scCp, p);
    
    //PutRNGstate();
    
    Memcpy(scal, scCp, psqr);
    Free(scCp); Free(tmp);
}

int InvMatrixUpper(double *A, int p)
{
      int info1, info2;
      F77_NAME(dpotrf)("U", &p, A, &p, &info1);
      F77_NAME(dpotri)("U", &p, A, &p, &info2);      
      //make sure you make it symmetric later...
      return(info1);
}


SEXP WM2_rmvGaussianR(SEXP mu, SEXP Sigma)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(Sigma, R_DimSymbol)), psqr, p;
    double *scCp, *ansp;//, one = 1, zero = 0;

    if (!isMatrix(Sigma) || !isReal(Sigma) || dims[0] != dims[1])
	error("Sigma must be a square, real matrix");
    
    p = dims[0];
    psqr=p*p;

    PROTECT(ans = allocVector(REALSXP, p));
    ansp = REAL(ans);
    Memcpy(ansp, REAL(mu), p);

    scCp = Memcpy(Calloc(psqr,double), REAL(Sigma), psqr);

    GetRNGstate();
    WM2_rmvGaussianC(ansp, scCp, p);
    PutRNGstate();
	
    Free(scCp);

    UNPROTECT(1);
    return ans;
}

void WM2_rmvGaussianC(double *mu, double *Sigma, int p)
{
  double ans[p];
  int info, psqr,j=0, intOne=1;
  double *scCp, one = 1; //zero = 0;
  
  psqr = p * p;
  scCp = Memcpy(Calloc(psqr,double), Sigma, psqr);

  F77_NAME(dpotrf)("L", &p, scCp, &p, &info);
  if (info)
    error("Sigma matrix is not positive-definite");
  
  //GetRNGstate();
  for(j=0;j<p;j++)
    {
      ans[j] = rnorm(0,1);
    }
  F77_NAME(dtrmv)("L","N","N", &p, scCp, &p, ans, &intOne);
  F77_NAME(daxpy)(&p, &one, ans, &intOne, mu, &intOne);
  //PutRNGstate();
  Free(scCp);
}


/* function HybridMC arguments */

/* ystart      : initial value  (after call will contain new values)                    */
/* n           : length of ystart                                                       */
/* epsLo       : the lower bound of the leapfrog step size                              */
/* epsRng      : the width of the leapfrog step size range                              */
/* LFsteps     : the number of leapfrog steps to take                                   */
/* compWeights : the "weights" or "masses" of the individual variables                  */


void WM2_hybridMC(double *ystart, int N, double epsLo, double epsRng, int LFsteps, double *compWeights, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMats, int *obsCovMats, int *sizeCovMats, int *parStart, double *means, double **pointerCovMat, double *logLike, int maxSizeCovMat, double *predVals, double *logPostCurrEff)
{
	int j=0,i=0;
	double e=0,epsilon=0,b=0;
	double H=0,cande=0,candH=0;
	double like1, like2;
	//Rprintf("Starting hybridMC\n");
	
	double predVals1[2*designDims[0]];
	double predVals2[2*designDims[0]];

	GetRNGstate();

	double m[N],sumM2=0;//mStart[N];
	double x[N], g[N];


	// make a copy of the starting value
	Memcpy(x,ystart,N);
	

	  // Sample new momenta, and compute the sum of squares
	  sumM2=0;
	  for(i=0;i<N;i++){
	    m[i] = rnorm(0,pow(compWeights[i],.5));
	    sumM2+=m[i]*m[i]/(2*compWeights[i]);
	  }
	  
	  // Copy the starting momenta for later
	  //memcpy(mStart,m,N*sizeof(double));

	  //printf("First Eval Hybrid\n");
	  // Evaluate the log-density and the derivative of the log density
	  
	  if(nCovMats==0){
	    e = *logPostCurrEff;
	  }else{
	    //printf("First Log Post\n");
	    e = LogPosterior(x, N, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, KeffectsSlope, KeffectsCov, Aeffects, AeffectsSlope, AeffectsCov, Geffects, GeffectsSlope, GeffectsCov, invGammaPriorA, invGammaPriorB, Kmu0, Ksig20, Amu0, Asig20, Gmu0, Gsig20, useA, Ktype, nCovMats, obsCovMats, sizeCovMats, parStart, means, pointerCovMat, &like1, maxSizeCovMat, predVals1);
	  }
	
		
	  gradLogPosterior(x, g, N, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, KeffectsSlope, KeffectsCov, Aeffects, AeffectsSlope, AeffectsCov, Geffects, GeffectsSlope, GeffectsCov, invGammaPriorA, invGammaPriorB, Kmu0, Ksig20, Amu0, Asig20, Gmu0, Gsig20, useA, Ktype, nCovMats, obsCovMats, sizeCovMats, parStart, means, pointerCovMat, maxSizeCovMat);

	  // Compute the hamiltonian
	  H = sumM2 + e;
	  sumM2 = 0;
	  
	  
	  // Determine the value of epsilon, the time discretization parameter
	  if(epsRng==0){
	    epsilon = epsLo;
	  }else{
	    epsilon = runif(epsLo,epsLo+epsRng);
	  }
	  
	  
	  //printf("Starting Leapfrog steps\n");

	  // Leapfrog steps
	  for(i=0;i<LFsteps;i++){
	    for(j=0;j<N;j++){
	      m[j] = m[j] - epsilon * g[j]/2;
	      x[j] = x[j] + epsilon * m[j]/compWeights[j];
	      //printf("i: %d,j: %d, epsilon: %f, m: %f, compWeights: %f, x: %f, g: %f, ystart: %f\n",i,j,epsilon,m[j],compWeights[j],x[j], g[j], ystart[j]);
	    }
	    gradLogPosterior(x, g, N, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, KeffectsSlope, KeffectsCov, Aeffects, AeffectsSlope, AeffectsCov, Geffects, GeffectsSlope, GeffectsCov, invGammaPriorA, invGammaPriorB, Kmu0, Ksig20, Amu0, Asig20, Gmu0, Gsig20, useA, Ktype, nCovMats, obsCovMats, sizeCovMats, parStart, means, pointerCovMat, maxSizeCovMat);

	    for(j=0;j<N;j++){
	      m[j] = m[j] - epsilon * g[j]/2;
	      if(i == (LFsteps-1)) sumM2+=m[j]*m[j]/(2*compWeights[j]);
	    }
	  }
	  

	  // Determine hamiltonian for candidate
	  cande = LogPosterior(x, N, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, KeffectsSlope, KeffectsCov, Aeffects, AeffectsSlope, AeffectsCov, Geffects, GeffectsSlope, GeffectsCov, invGammaPriorA, invGammaPriorB, Kmu0, Ksig20, Amu0, Asig20, Gmu0, Gsig20, useA, Ktype, nCovMats, obsCovMats, sizeCovMats, parStart, means, pointerCovMat, &like2, maxSizeCovMat, predVals2);
	  candH = cande + sumM2;

	  //printf("H: %f, candH: %f, candx[0]: %f\n",H,candH,x[0]);

	  // Metropolis-Hastings decision
	  b=-log(runif(0,1));
	  if(b>(candH-H)){
	    Memcpy(ystart,x,N);
	    *logLike = like2;
	    Memcpy(predVals,predVals2,2*designDims[0]);
	    *logPostCurrEff=cande;
	  }else{
	    if(nCovMats>0) *logLike = like1;
	    Memcpy(predVals,predVals1,2*designDims[0]);
	  }
	 
	//Rprintf("Cur logPost: %f\n",logPostCurrEff[0]);
	PutRNGstate();	
	
	//printf("Ending HybridMC\n");
}

void WM2_metrop(double *ystart, int N, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, 
double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, 
int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, 
double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, 
int Ktype, int nCovMats, int *obsCovMats, int *sizeCovMats, int *parStart, double *means, double **pointerCovMat, 
double *logLike, int maxSizeCovMat, double *predVals, double *metropSD,double *logPostCurrEff)
{
	int i=0;
	double e=0,b=0;
	double cande=0;
	double like1, like2;
	//printf("Starting metrop\n");
	
	double predVals1[2*designDims[0]];
	double predVals2[2*designDims[0]];
	double x[N];

	// make a copy of the starting value
	Memcpy(x,ystart,N);

	GetRNGstate();
	
	if(nCovMats==0){
	  e = *logPostCurrEff;
	}else{
	  e = LogPosterior(x, N, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, KeffectsSlope, KeffectsCov, Aeffects, AeffectsSlope, AeffectsCov, Geffects, GeffectsSlope, GeffectsCov, invGammaPriorA, invGammaPriorB, Kmu0, Ksig20, Amu0, Asig20, Gmu0, Gsig20, useA, Ktype, nCovMats, obsCovMats, sizeCovMats, parStart, means, pointerCovMat, &like1, maxSizeCovMat, predVals1);
	}
	for(i=0;i<N;i++)
	  x[i] = x[i] + rnorm(0,metropSD[i]);

	  cande = LogPosterior(x, N, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, KeffectsSlope, KeffectsCov, Aeffects, AeffectsSlope, AeffectsCov, Geffects, GeffectsSlope, GeffectsCov, invGammaPriorA, invGammaPriorB, Kmu0, Ksig20, Amu0, Asig20, Gmu0, Gsig20, useA, Ktype, nCovMats, obsCovMats, sizeCovMats, parStart, means, pointerCovMat, &like2, maxSizeCovMat, predVals2);




	  // Metropolis-Hastings decision
	  b=-log(runif(0,1));
	  if(b>(cande-e)){
	    Memcpy(ystart,x,N);
	    *logLike = like2;
	    Memcpy(predVals,predVals2,2*designDims[0]);
	    *logPostCurrEff = cande;
	  }else{
	    if(nCovMats>0) *logLike = like1;
	    Memcpy(predVals,predVals1,2*designDims[0]);
	    
	  }
	 
	PutRNGstate();	
	
	//printf("Ending metrop\n");
}


/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
SEXP 
alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error("negative extents to 3D array");
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error("alloc3Darray: too many elements specified");
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}



/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str)
{
SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
int i;
for (i = 0; i < length(list); i++)
if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
elmt = VECTOR_ELT(list, i);
break;
}
return elmt;
}

SEXP RLogLikelihood(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP Aeffects, SEXP Geffects, SEXP useA, SEXP Ktype)
{
  int p = length(params);
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
  double logLike;
  SEXP returnLogLike;
  double predVals[2*designDims[0]];

  PROTECT(returnLogLike = allocVector(REALSXP,1));

  logLike = LogLikelihood(REAL(params), p, INTEGER_POINTER(nHit), INTEGER_POINTER(nMiss), INTEGER_POINTER(nFA), INTEGER_POINTER(nCR),
			 INTEGER_POINTER(setSize), INTEGER_POINTER(design), REAL(designCont), designDims, INTEGER_POINTER(Keffects),
			  INTEGER_POINTER(Aeffects), INTEGER_POINTER(Geffects), INTEGER_VALUE(useA), INTEGER_VALUE(Ktype),predVals);
  

  *REAL(returnLogLike) = logLike;
  UNPROTECT(1);
  return(returnLogLike);
}

SEXP RPredictedProbabilities(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP Aeffects, SEXP Geffects, SEXP useA, SEXP Ktype)
{
  int p = length(params);
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
  
  SEXP predVals;
  PROTECT(predVals = allocMatrix(REALSXP,designDims[0],2));
  double *pPredVals = REAL(predVals);

  double logLike;


  logLike = LogLikelihood(REAL(params), p, INTEGER_POINTER(nHit), INTEGER_POINTER(nMiss), INTEGER_POINTER(nFA), INTEGER_POINTER(nCR),
			 INTEGER_POINTER(setSize), INTEGER_POINTER(design), REAL(designCont), designDims, INTEGER_POINTER(Keffects),
			  INTEGER_POINTER(Aeffects), INTEGER_POINTER(Geffects), INTEGER_VALUE(useA), INTEGER_VALUE(Ktype), pPredVals);
  

  UNPROTECT(1);
  return(predVals);
}


SEXP RLogPosteriorNoCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype)
{
  int p = length(params);
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
  double logPost=0;
  double logLike;
  SEXP returnLogPost;
  double predVals[2*designDims[0]];

  PROTECT(returnLogPost = allocVector(REALSXP,1));

  logPost = LogPosterior(REAL(params), p, INTEGER_POINTER(nHit), INTEGER_POINTER(nMiss), INTEGER_POINTER(nFA), INTEGER_POINTER(nCR),
			 INTEGER_POINTER(setSize), INTEGER_POINTER(design), REAL(designCont), designDims, INTEGER_POINTER(Keffects),
			 INTEGER_POINTER(KeffectsSlope), INTEGER_POINTER(KeffectsCov), INTEGER_POINTER(Aeffects),
			 INTEGER_POINTER(AeffectsSlope), INTEGER_POINTER(AeffectsCov), INTEGER_POINTER(Geffects),
			 INTEGER_POINTER(GeffectsSlope), INTEGER_POINTER(GeffectsCov), REAL(invGammaPriorA)[0], REAL(invGammaPriorB)[0], 
			 REAL(Kmu0)[0], REAL(Ksig20)[0], REAL(Amu0)[0], REAL(Asig20)[0], REAL(Gmu0)[0], REAL(Gsig20)[0], INTEGER_VALUE(useA), 
			 INTEGER_VALUE(Ktype), 0, (int *)(0), (int *)(0), (int *)(0), (double *)(0), (double **)(0), &logLike, 0, predVals);
  

  *REAL(returnLogPost) = logPost;
  UNPROTECT(1);
  return(returnLogPost);
}


SEXP RgradLogPosteriorNoCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype)
{
  int p = length(params);
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
 
  double *grad;
  SEXP returnLogPost;
  PROTECT(returnLogPost = allocVector(REALSXP,p));
  grad = REAL(returnLogPost);

  gradLogPosterior(REAL(params), grad, p, INTEGER_POINTER(nHit), INTEGER_POINTER(nMiss), INTEGER_POINTER(nFA), INTEGER_POINTER(nCR),
			 INTEGER_POINTER(setSize), INTEGER_POINTER(design), REAL(designCont), designDims, INTEGER_POINTER(Keffects),
			 INTEGER_POINTER(KeffectsSlope), INTEGER_POINTER(KeffectsCov), INTEGER_POINTER(Aeffects),
			 INTEGER_POINTER(AeffectsSlope), INTEGER_POINTER(AeffectsCov), INTEGER_POINTER(Geffects),
			 INTEGER_POINTER(GeffectsSlope), INTEGER_POINTER(GeffectsCov), REAL(invGammaPriorA)[0], REAL(invGammaPriorB)[0], 
			 REAL(Kmu0)[0], REAL(Ksig20)[0], REAL(Amu0)[0], REAL(Asig20)[0], REAL(Gmu0)[0], REAL(Gsig20)[0], INTEGER_VALUE(useA), 
			 INTEGER_VALUE(Ktype), 0, (int *)(0), (int *)(0), (int *)(0), (double *)(0), (double **)(0), 0);
  

  UNPROTECT(1);
  return(returnLogPost);
}

SEXP RLogPosteriorWithCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype, SEXP covList, SEXP means, SEXP obsCovMat, SEXP sizeCovMat, SEXP parStart, SEXP covEffSlope)
{
  //printf("Calling log likelihood...\n");
  int i=0,q;
  int p = length(params);
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
  double logPost=0;
  double logLike;
  double predVals[2*designDims[0]];

  SEXP returnLogPost;
  PROTECT(returnLogPost = allocVector(REALSXP,1));
  
  //printf("Define some variables...\n");
  int nCovMat = length(covList);
  int maxSizeCov=0;
  double *pointerCovMat[nCovMat];

  //printf("Getting pointers...\n");
  for(i=0;i<nCovMat;i++)
    {
      pointerCovMat[i]=REAL(VECTOR_ELT(covList, i));
      //printf("in loop %d of %d\n",i,nCovMat);
      //printf("%f %f %f %f\n",(pointerCovMat[i])[0],(pointerCovMat[i])[1],(pointerCovMat[i])[2],(pointerCovMat[i])[3]);
      q = INTEGER_POINTER(sizeCovMat)[i];
      //printf("q: %d\n",q);
      if(q>maxSizeCov) maxSizeCov = q;
    }

  //printf("Compute log posterior...\n");
  logPost = LogPosterior(REAL(params), p, INTEGER_POINTER(nHit), INTEGER_POINTER(nMiss), INTEGER_POINTER(nFA), INTEGER_POINTER(nCR),
			 INTEGER_POINTER(setSize), INTEGER_POINTER(design), REAL(designCont), designDims, INTEGER_POINTER(Keffects),
			 INTEGER_POINTER(KeffectsSlope), INTEGER_POINTER(KeffectsCov), INTEGER_POINTER(Aeffects),
			 INTEGER_POINTER(AeffectsSlope), INTEGER_POINTER(AeffectsCov), INTEGER_POINTER(Geffects),
			 INTEGER_POINTER(GeffectsSlope), INTEGER_POINTER(GeffectsCov), REAL(invGammaPriorA)[0], REAL(invGammaPriorB)[0], 
			 REAL(Kmu0)[0], REAL(Ksig20)[0], REAL(Amu0)[0], REAL(Asig20)[0], REAL(Gmu0)[0], REAL(Gsig20)[0], INTEGER_VALUE(useA), 
			 INTEGER_VALUE(Ktype), nCovMat, INTEGER_POINTER(obsCovMat), INTEGER_POINTER(sizeCovMat), INTEGER_POINTER(parStart), REAL(means), pointerCovMat, &logLike, maxSizeCov, predVals);
  

  *REAL(returnLogPost) = logPost;
  UNPROTECT(1);
  return(returnLogPost);
}


SEXP RgradLogPosteriorWithCov(SEXP params, SEXP nHit, SEXP nMiss, SEXP nFA, SEXP nCR, SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, SEXP useA, SEXP Ktype, SEXP covList, SEXP means, SEXP obsCovMat, SEXP sizeCovMat, SEXP parStart, SEXP covEffSlope)
{
  int p = length(params),q,i=0;
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
 
  double *grad;
  SEXP returnLogPost;
  PROTECT(returnLogPost = allocVector(REALSXP,p));
  grad = REAL(returnLogPost);


  int nCovMat = length(covList);
  int maxSizeCov=0;
  double *pointerCovMat[nCovMat];

  for(i=0;i<nCovMat;i++)
    {
      pointerCovMat[i]=REAL(VECTOR_ELT(covList, i));
      q = INTEGER(sizeCovMat)[i];
      if(q>maxSizeCov) maxSizeCov = q;
    }


  gradLogPosterior(REAL(params), grad, p, INTEGER_POINTER(nHit), INTEGER_POINTER(nMiss), INTEGER_POINTER(nFA), INTEGER_POINTER(nCR),
			 INTEGER_POINTER(setSize), INTEGER_POINTER(design), REAL(designCont), designDims, INTEGER_POINTER(Keffects),
			 INTEGER_POINTER(KeffectsSlope), INTEGER_POINTER(KeffectsCov), INTEGER_POINTER(Aeffects),
			 INTEGER_POINTER(AeffectsSlope), INTEGER_POINTER(AeffectsCov), INTEGER_POINTER(Geffects),
			 INTEGER_POINTER(GeffectsSlope), INTEGER_POINTER(GeffectsCov), REAL(invGammaPriorA)[0], REAL(invGammaPriorB)[0], 
			 REAL(Kmu0)[0], REAL(Ksig20)[0], REAL(Amu0)[0], REAL(Asig20)[0], REAL(Gmu0)[0], REAL(Gsig20)[0], INTEGER_VALUE(useA), 
			 INTEGER_VALUE(Ktype), nCovMat, INTEGER_POINTER(obsCovMat), INTEGER_POINTER(sizeCovMat), INTEGER_POINTER(parStart), 
			 REAL(means), pointerCovMat, maxSizeCov);
  
  

  UNPROTECT(1);
  return(returnLogPost);
}



double LogPosterior(double *params, int p, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMat, int *obsCovMat, int *sizeCovMat, int *parStart, double *means, double **pointerCovMat, double *logLike, int maxSizeCovMat, double *predVals)
{
  int i=0,j=0,paramCounter=0,nDesign=designDims[1];
  double logPost=0, sumSq=0,mean=0;
  double a0=invGammaPriorA;
  double b0=invGammaPriorB;

  //printf("Starting Log Posterior\n");
 
  
  logPost += LogLikelihood(params, p, nHit, nMiss, nFA, nCR, setSize, design, designCont, designDims, Keffects, Aeffects, Geffects, useA, Ktype, predVals);
  logLike[0] = logPost;
  //Rprintf("loglike: %f\n",logPost);
  
  logPost += -.5*pow(params[0]-Kmu0,2)/Ksig20;
  if(useA){ 
	logPost += -.5*pow(params[1]-Amu0,2)/Asig20; 
	}
  logPost += -.5*pow(params[2]-Gmu0,2)/Gsig20;
  
  if(useA){
	paramCounter=3;
  }else{
	paramCounter=2;
  }
  for(i=0;i<nDesign;i++)
    {
      if(Keffects[i]>0 && KeffectsCov[i]==0)
		{
			sumSq=0;
			mean=0;
			for(j=0;j<Keffects[i];j++)
			{
				sumSq += params[paramCounter+j]*params[paramCounter+j]; 
				if(KeffectsSlope[i])
				  {
				    mean+=params[paramCounter+j]/Keffects[i];
				  }
			}
			
			if(KeffectsSlope[i])
			  {
			    logPost += -(a0+.5*Keffects[i]-.5)*log(b0 + .5*(sumSq - Keffects[i]*mean*mean ));
			  }else{
			  logPost += -(a0+.5*Keffects[i])*log(b0 + .5*sumSq);
			}
		}
		paramCounter += Keffects[i];
      
      
	  if(Aeffects[i]>0 && useA && AeffectsCov[i]==0)
		{
			sumSq=0;
			mean=0;
			for(j=0;j<Aeffects[i];j++)
			{
				sumSq += params[paramCounter+j]*params[paramCounter+j]; 
				if(AeffectsSlope[i])
				  {
				    mean+=params[paramCounter+j]/Aeffects[i];
				  }
			}
			if(AeffectsSlope[i])
			  {
			    logPost += -(a0+.5*Aeffects[i]-.5)*log(b0 + .5*(sumSq - Aeffects[i]*mean*mean ));
			  }else
			  {
			    logPost += -(a0+.5*Aeffects[i])*log(b0 + .5*sumSq);
			}
		}
		paramCounter += Aeffects[i];
      
      
	  if(Geffects[i]>0 && GeffectsCov[i]==0)
		{
			sumSq=0;
			mean=0;
			for(j=0;j<Geffects[i];j++)
			{
				sumSq += params[paramCounter+j] * params[paramCounter+j]; 
				if(GeffectsSlope[i])
				  {
				    mean+=params[paramCounter+j]/Geffects[i];
				  }
			}
			if(GeffectsSlope[i])
			  {
			    logPost += -(a0+.5*Geffects[i]-.5)*log(b0 + .5*(sumSq - Geffects[i]*mean*mean ));
			  }else
			  {
			    logPost += -(a0+.5*Geffects[i])*log(b0 + .5*sumSq);
			  }
		}
		paramCounter += Geffects[i];

    }
  
  double *vecWorkspace, val;
  int m=0;

  //Do covariance matrices, if any
  if(nCovMat>0){
    for(i=0;i<nCovMat;i++)
      {
	sumSq=0;
	vecWorkspace = Calloc(obsCovMat[i]*sizeCovMat[i], double);
	
	
	//printf("LP: means: %f %f \n",means[0*nCovMat + i],means[1*nCovMat + i]);
	//Rprintf("LP: Prec %d: %f %f %f %f\n",i,pointerCovMat[i][0],pointerCovMat[i][1],pointerCovMat[i][2],pointerCovMat[i][3]);
	for(j=0;j<obsCovMat[i];j++)
	  {
	    for(m=0;m<sizeCovMat[i];m++)
	      {
		val = params[parStart[m*nCovMat + i] + j];
		val = val - means[m*nCovMat + i];
		*(vecWorkspace+m*obsCovMat[i]+j)=val;
	      }
	  }
	 
	for(j=0;j<obsCovMat[i];j++)
	  {
	    sumSq += WM2_quadform(vecWorkspace+j,pointerCovMat[i],sizeCovMat[i],obsCovMat[i]);
	  }
	logPost += -.5*sumSq;
	Free(vecWorkspace);
      }
  }
  
  //for(i=0;i<p;i++) Rprintf("%f ",params[i]);
  //Rprintf("\n%f %f %f %f\n",pointerCovMat[0][0],pointerCovMat[0][1],pointerCovMat[0][2],pointerCovMat[0][3]);
  //Rprintf("\nEnding Log Posterior: %f\n",logPost);

  return(-logPost);

}


double LogLikelihood(double *params, int p, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *Aeffects, int *Geffects, int useA, int Ktype, double *predVals)
{
  
  int h=0,i=0,N=designDims[0],paramCounter=0,nDesign=designDims[1];
  double K=0,G=0,A=0,PH=0,PF=0,logLike=0,d=0,d2=0;
  double effectiveG=0;
  
  //Rprintf("Starting Log Likelihood\n");
	void R_CheckUserInterrupt(void);
  for(h=0;h<N;h++)
    {
    K = params[0];
    if(useA){
		A = params[1];
		G = params[2];
		paramCounter = 3;
    }else{
		G = params[1];
		paramCounter = 2;

    }
    //printf("RESET parametercount: %i\n",paramCounter);
    
    for(i=0;i<nDesign;i++)
      {
	//printf("OUT: h=%i N=%i i=%i ELT=%i\n",h,N,i,i*N+h);
	//printf("EFF: K=%i A=%i G=%i\n",pKeffects[i],pAeffects[i],pGeffects[i]);
	if(Keffects[i]>0)
	  {
 	    //printf("K: pc=%i ELT=%i\n",paramCounter,paramCounter+pdesign[i*N + h]-1);
	    K = K + designCont[i*N + h]*params[paramCounter+design[i*N + h]-1];
		paramCounter += Keffects[i];	
	  }
	
	
	if(Aeffects[i]>0 && useA)
	  {
	    //printf("A: pc=%i ELT=%i\n",paramCounter,paramCounter+pdesign[i*N + h]-1);
	    A = A + designCont[i*N + h]*params[paramCounter+design[i*N + h]-1];
		paramCounter += Aeffects[i];	
	  }
	
	
	if(Geffects[i]>0)
		{
	    //printf("G: pc=%i ELT=%i\n",paramCounter,paramCounter+pdesign[i*N + h]-1);
	    G = G + designCont[i*N + h]*params[paramCounter+design[i*N + h]-1];
		paramCounter += Geffects[i];	
		}
	
    }		
    
    //printf("DONE: pc=%i h=%i\n",paramCounter,h);

    G = plogis(G,0,1,TRUE,FALSE);
    effectiveG=G;
    if(K<0)
      {
		K=0; 
      }else if(K>setSize[h])
      {
		K=setSize[h];
		effectiveG=0;
      }
    d = K/setSize[h];
    if(useA)
      {
	A = plogis(A,0,1,TRUE,FALSE);
      }else{
	A = 1;
      }
    
    if(Ktype==0)
      {
		PH = (1-A)*G + A*d + A*(1-d)*G;
		PF = (1-A)*G + A*(1-d)*G;
      }else if(Ktype==1){
		PH = (1-A)*G + A*d + A*(1-d)*G;
		PF = (1-A)*G + A*effectiveG;
	  }else if(Ktype==2){
		if(K>(setSize[h]-1)){
			d2 = 1;
		}else{
			d2 = 1 - (1-K/(1.0*setSize[h]))*(1 - K/(1.0*setSize[h]-1));
		}
		PH = (1-A)*G + A*d2 + A*(1-d2)*G;
		PF = (1-A)*G + A*(1-d)*G;
	  }
	logLike += nHit[h]*log(PH) + nMiss[h]*log(1-PH) + nFA[h]*log(PF) + nCR[h]*log(1-PF);
    predVals[h] = PH;
    predVals[h + N] = PF;
    //printf("Computed loglike\n");	
    }

  //printf("Ending Log Likelihood\n");
  return(logLike);
  
}
  
void gradLogPosterior(double *params, double *grad, int p, int *nHit, int *nMiss, int *nFA, int *nCR, int *setSize, int *design, double *designCont, int *designDims, int *Keffects, int *KeffectsSlope, int *KeffectsCov, int *Aeffects, int *AeffectsSlope, int *AeffectsCov, int *Geffects, int *GeffectsSlope, int *GeffectsCov, double invGammaPriorA, double invGammaPriorB, double Kmu0, double Ksig20, double Amu0, double Asig20, double Gmu0, double Gsig20, int useA, int Ktype, int nCovMat, int *obsCovMat, int *sizeCovMat, int *parStart, double *means, double **pointerCovMat, int maxSizeCovMat)
{
  int h=0,i=0,j=0,paramCounter=0,N=designDims[0],nDesign=designDims[1];
  double sumSq=0,mean=0;
  double tf=0,th=0,dK=0,dA=0,dG=0;
  double k=0,a=0,g=0,d=0,d2=0,PH=0,PF=0;
  double A=0,G=0;
  double a0=invGammaPriorA;
  double b0=invGammaPriorB;
  double effectiveG=0;
  double gradLogPost[p];
 
  //printf("Starting grad\n");
void R_CheckUserInterrupt(void);
  for(i=0;i<p;i++)
    {
      gradLogPost[i]=0;
    }

  //printf("First**************************g[14]: %f\n",-gradLogPost[14]);


  for(h=0;h<N;h++)
    {
    k = params[0];
    if(useA){
		a = params[1];
		g = params[2];
		paramCounter = 3;
    }else{
		g = params[1];
		paramCounter = 2;
	}
    for(i=0;i<nDesign;i++)
      {
	
	if(Keffects[i]>0)
	  {
	    k = k + designCont[i*N + h]*params[paramCounter+design[N*i+h]-1];
	    //printf("h: %d, paramCounter: %d, N: %d, i: %d, k0: %f, designCont: %f, design:%d, param: %f\n",h,paramCounter, N, i,params[0],designCont[i*N + h],design[N*i+h],params[paramCounter+design[N*i+h]-1]);
	    paramCounter += Keffects[i];
	  }
	
	
	if(Aeffects[i]>0 && useA)
	  {
	    a = a + designCont[i*N + h]*params[paramCounter+design[N*i+h]-1];
		paramCounter += Aeffects[i];
	  }
	
	
	if(Geffects[i]>0)
	  {
	    g = g + designCont[i*N + h]*params[paramCounter+design[N*i+h]-1];
	    //printf("h: %d, paramCounter: %d, N: %d, i: %d, g0: %f, designCont: %f, design:%d, param: %f\n",h,paramCounter, N, i,params[3],designCont[i*N + h],design[N*i+h],params[paramCounter+design[N*i+h]-1]);

		paramCounter += Geffects[i];
	  }
	
	
      }
   
      G = plogis(g,0,1,TRUE,FALSE);
      if(k<0)
      {
		d=0; 
		effectiveG=G;
      }else if(k>setSize[h])
      {
		d=1;
		effectiveG=0;
      }else
      {
		d=k/setSize[h];
		effectiveG=G;
	  }
    if(useA){
		A = plogis(a,0,1,TRUE,FALSE);
    }else{
		A = 1;
     }
	
    if(Ktype==0)
      {
		PH = (1-A)*G + A*d + A*(1-d)*G;
		PF = (1-A)*G + A*(1-d)*G;
      }else if(Ktype==1){
		PH = (1-A)*G + A*d + A*(1-d)*G;
		PF = (1-A)*G + A*effectiveG;
	  }else if(Ktype==2){
		if(k>(setSize[h]-1)){
			d2 = 1;
		}else{
			d2 = 1 - (1-k/(1.0*setSize[h]))*(1 - k/(1.0*setSize[h]-1));
		}
		PH = (1-A)*G + A*d2 + A*(1-d2)*G;
		PF = (1-A)*G + A*(1-d)*G;
	  }
	//PF = (1-A)/2 + A*(1-d)*G;
    
    th = (nHit[h] - PH*(nHit[h]+nMiss[h])) / (PH*(1-PH));
    tf = (nFA[h] - PF*(nFA[h]+nCR[h])) / (PF*(1-PF));
    
	/* Changed model 
	if(Ktype==0){
		dK = (A*(1-G)/setSize[h]*th + -A*G/setSize[h]*tf)*(int)(k>0 && k<setSize[h]);
		dG = A*(1-d)*(th + tf)*dlogis(g,0,1,FALSE)*(int)(k<setSize[h]);
		if(useA){
			dA = ((-.5 + d + (1-d)*G)*th + (-.5 + (1-d)*G)*tf)*dlogis(a,0,1,FALSE);
		}
	}else if(Ktype==1){
		dK = (A*(1-G)/setSize[h]*th)*(int)(k>0 && k<setSize[h]);
		dG = (A*(1-d)*th + A*tf)*dlogis(g,0,1,FALSE)*(int)(k<setSize[h]);
		if(useA){
			dA = ((-.5 + d + (1-d)*G)*th + (-.5 + G)*tf)*dlogis(a,0,1,FALSE);
		}
	}
	*/
	if(Ktype==0){
		dK = ((1-G)*th - G*tf)*A/setSize[h]*(int)(k>0 && k<setSize[h]); 
		dG = (1-A*d)*(th + tf)*dlogis(g,0,1,FALSE);
		if(useA){
			dA = (d*(1-G)*th - d*G*tf)*dlogis(a,0,1,FALSE);
		}
	}else if(Ktype==1){
		dK = (A*(1-effectiveG)/setSize[h]*th)*(int)(k>0 && k<setSize[h]);
		dG = ((1-A + A*(1-d)*(int)(k<setSize[h]))*th + (1-A + A*(int)(k<setSize[h]))*tf)*dlogis(g,0,1,FALSE);
		if(useA){
			dA = ((-G + d + (1-d)*effectiveG)*th + (-G + effectiveG)*tf)*dlogis(a,0,1,FALSE);
		}
	}else if(Ktype==2){
		if(k>(setSize[h]-1)){
			d2 = 1;
		}else{
			d2 = 1 - (1-k/(1.0*setSize[h]))*(1 - k/(1.0*setSize[h]-1));
		}
		dK = (1-G)*A*th*(2.0*setSize[h]-1-2*k)/(1.0*setSize[h]*(setSize[h]-1)) * (int)(k>0 && k<(setSize[h]-1)) +
			  -G*A*tf/(1.0*setSize[h]) * (int)(k>0 && k<setSize[h]);		
		dG = ((1-A*d2)*th + (1-A*d)*tf)*dlogis(g,0,1,FALSE);
		if(useA){
			dA = (d2*(1-G)*th - d*G*tf)*dlogis(a,0,1,FALSE);
		}
	}
	
	//printf("h=%i H=%i M=%i F=%i C=%i pH=%f pF=%f k=%f A=%f G=%f N=%i dK=%f\n",h,nHit[h],nMiss[h],nFA[h],nCR[h],PH,PF,k,A,G,setSize[h],dK);
    
	
	//printf("2**************************g[14]: %f\n",-gradLogPost[14]);

    gradLogPost[0] += dK;
    if(useA){
		gradLogPost[1] += dA;
		gradLogPost[2] += dG;
		paramCounter=3;
	}else{
		gradLogPost[1] += dG;
		paramCounter=2;
	}
	

   
    for(i=0;i<nDesign;i++)
      {
	if(Keffects[i]>0)
	  {
	    gradLogPost[paramCounter+design[N*i+h]-1] += dK * designCont[i*N + h];
	    //printf("h=%i i=%i pc=%i dK=%f ELT=%i pd=%i\n",h,i,paramCounter,dK,paramCounter+design[N*i+h]-1,design[N*i+h]);
		paramCounter += Keffects[i];
	  }
	
	
	if(Aeffects[i]>0 && useA)
	  {
	    gradLogPost[paramCounter+design[N*i+h]-1] += dA * designCont[i*N + h];
	    paramCounter += Aeffects[i];
	  }
	
	
	if(Geffects[i]>0)
	  {
	    gradLogPost[paramCounter+design[N*i+h]-1] += dG * designCont[i*N + h];
	    paramCounter += Geffects[i];
	  }
	
      }

    }
    
  //printf("3**************************g[14]: %f\n",-gradLogPost[14]);

  gradLogPost[0] += -(params[0] - Kmu0)/Ksig20;
  if(useA){
	gradLogPost[1] += -(params[1] - Amu0)/Asig20;
	gradLogPost[2] += -(params[2] - Gmu0)/Gsig20;
	paramCounter=3;
  }else{
	gradLogPost[1] += -(params[1] - Gmu0)/Gsig20;
	paramCounter=2;
  }
  
  for(i=0;i<nDesign;i++)
    {
      if(Keffects[i]>0 && KeffectsCov[i]==0)
	{
	  sumSq=0;
	  mean=0;
	  for(j=0;j<Keffects[i];j++)
	    {
	      sumSq += params[paramCounter+j] * params[paramCounter+j]; 
	      if(KeffectsSlope[i])
		{
		  mean += params[paramCounter+j]/Keffects[i];
		}
	    }
	  for(j=0;j<Keffects[i];j++)
	    {
	      if(KeffectsSlope[i])
		{
		  gradLogPost[paramCounter+j] += -(a0+.5*Keffects[i]-.5) / (b0 + .5*(sumSq - Keffects[i]*mean*mean)) * (params[paramCounter+j] - mean);
		}else
		{
		  gradLogPost[paramCounter+j] += -(a0+.5*Keffects[i]) / (b0 + .5*sumSq) * params[paramCounter+j];
		}
	    }
	}
	paramCounter += Keffects[i];

      
      
      if(Aeffects[i]>0 && useA && AeffectsCov[i]==0)
	{
	  sumSq=0;
	  mean=0;
	  for(j=0;j<Aeffects[i];j++)
	    {
	      sumSq += params[paramCounter+j] * params[paramCounter+j]; 
	      if(AeffectsSlope[i])
		{
		  mean += params[paramCounter+j]/Aeffects[i];
		}
	    }
	  for(j=0;j<Aeffects[i];j++)
	    {
	      if(AeffectsSlope[i])
		{
		  gradLogPost[paramCounter+j] += -(a0+.5*Aeffects[i]-.5) / (b0 + .5*(sumSq - Aeffects[i]*mean*mean)) * (params[paramCounter+j]-mean);
		}else
		{
		  gradLogPost[paramCounter+j] += -(a0+.5*Aeffects[i]) / (b0 + .5*sumSq) * params[paramCounter+j];
		}
	    }
	}
	paramCounter += Aeffects[i];
            
      if(Geffects[i]>0 && GeffectsCov[i]==0)
	{
	  sumSq=0;
	  mean=0;
	  for(j=0;j<Geffects[i];j++)
	    {
	      sumSq += params[paramCounter+j] * params[paramCounter+j]; 
	      if(GeffectsSlope[i])
		{
		  mean += params[paramCounter+j]/Geffects[i];
		}
	    }
	  for(j=0;j<Geffects[i];j++)
	    {
	      if(GeffectsSlope[i])
		{
		  gradLogPost[paramCounter+j] += -(a0+.5*Geffects[i]-.5) / (b0 + .5*(sumSq - Geffects[i]*mean*mean)) * (params[paramCounter+j] - mean);
		}else
		{
		  gradLogPost[paramCounter+j] += -(a0+.5*Geffects[i]) / (b0 + .5*sumSq) * params[paramCounter+j];
		}
	    }
	    
	}
        paramCounter += Geffects[i];
      
    }
  
  //printf("4**************************g[14]: %f\n",-gradLogPost[14]);


  double *tmparray, *vecWorkspace, val=0;
  int m=0, iOne=1;
  double minusone = -1, zero=0;

  //Do covariance matrices, if any
  if(nCovMat>0){
    for(i=0;i<nCovMat;i++)
      {
	//allocate some workspace
	tmparray=Calloc(sizeCovMat[i],double);
	vecWorkspace = Calloc(obsCovMat[i]*sizeCovMat[i],double);

	
	for(j=0;j<obsCovMat[i];j++)
	  {
	    for(m=0;m<sizeCovMat[i];m++)
	      {
		val = params[parStart[m*nCovMat + i] + j];
		val = val - means[m*nCovMat + i];
		*(vecWorkspace+m*obsCovMat[i]+j)=val;
	      }
	  }
	for(j=0;j<obsCovMat[i];j++)
	  {
	    F77_CALL(dsymv)("U", &(sizeCovMat[i]), &minusone, pointerCovMat[i], &(sizeCovMat[i]), vecWorkspace+j, &(obsCovMat[i]), &zero, tmparray, &iOne);
	    
	    for(m=0;m<sizeCovMat[i];m++)
	      {
		gradLogPost[parStart[m*nCovMat + i] + j] += tmparray[m];
	      }
	  }
	Free(tmparray);
	Free(vecWorkspace);
      }

  }

  //Memcpy(grad, gradLogPost, p);
  
  for(i=0;i<p;i++)
    grad[i] = -gradLogPost[i];
  
  //printf("Last**************************g[14]: %f\n",grad[14]);
  //printf("Last2**************************g[14]: %f\n",gradLogPost[14]);

  //printf("Ending grad\n");

}

SEXP WM2_GibbsSampler(SEXP iters, SEXP starteffects, SEXP nCovMatR, SEXP obsCovMatR, SEXP sizeCovMatR, SEXP parStartR, SEXP covEffSlopeR, SEXP nHit, SEXP nMiss, SEXP nFA,
		      SEXP nCR,  SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, 
		      SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, 
		      SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, 
		      SEXP useA, SEXP Ktype, SEXP epsLoR, SEXP epsRngR, SEXP LFstepsR, SEXP compWeightsR, SEXP dfWish, SEXP progress, SEXP pBar, SEXP rho,              SEXP storePred, SEXP metrop, SEXP metropSD, SEXP metropThin)
{

  //printf("Starting C call...\n");

  int iIters = INTEGER_VALUE(iters),m, i,j,q,k,n, p,pSub,counter=0, counter2=0;
  int nEffectPars = length(starteffects), nCovMat=INTEGER_VALUE(nCovMatR), iMetropThin = INTEGER_VALUE(metropThin);
  SEXP effectPars, meanPars, covMatsList;
  SEXP returnList, logLikes;
  double *pLogLikes, tmp, *tmpmeans, tmpmean=0;
  int iStorePred = INTEGER_VALUE(storePred);
  int iMetrop = INTEGER_VALUE(metrop); 
  double *dMetropSD = REAL(metropSD);

  if(!iMetrop) iMetropThin = 1;
  if(iMetropThin>iIters) error("Thinning interval larger than number of iterations.");
  int saveIters = iIters/iMetropThin;


  double dfWishPrior = REAL(dfWish)[0];
  
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2],nMeanPars=0;
  designDims[0]=dims[0];
  designDims[1]=dims[1];


  int *pHit = INTEGER_POINTER(nHit);
  int *pMiss = INTEGER_POINTER(nMiss);
  int *pFA = INTEGER_POINTER(nFA);
  int *pCR = INTEGER_POINTER(nCR);
  int *pSetSize = INTEGER_POINTER(setSize);
  int *pDesign = INTEGER_POINTER(design);
  double *pDesignCont = REAL(designCont);
  int *pKeffects = INTEGER_POINTER(Keffects);
  int *pKeffectsCov = INTEGER_POINTER(KeffectsCov);
  int *pKeffectsSlope = INTEGER_POINTER(KeffectsSlope);
  int *pAeffects = INTEGER_POINTER(Aeffects);
  int *pAeffectsCov = INTEGER_POINTER(AeffectsCov);
  int *pAeffectsSlope = INTEGER_POINTER(AeffectsSlope);
  int *pGeffects = INTEGER_POINTER(Geffects);
  int *pGeffectsCov = INTEGER_POINTER(GeffectsCov);
  int *pGeffectsSlope = INTEGER_POINTER(GeffectsSlope);

  double params[nEffectPars];
  Memcpy(params,REAL(starteffects),nEffectPars);
  double logPostCurrEff=0;

  double epsLo = REAL(epsLoR)[0];
  double epsRng = REAL(epsRngR)[0];
  double *compWeights = REAL(compWeightsR);
  int LFsteps = INTEGER_VALUE(LFstepsR);

  double a0 = REAL(invGammaPriorA)[0];
  double b0 = REAL(invGammaPriorB)[0];
  double dKmu0 = REAL(Kmu0)[0];
  double dKsig20 = REAL(Ksig20)[0];
  double dAmu0 = REAL(Amu0)[0];
  double dAsig20 = REAL(Asig20)[0];
  double dGmu0 = REAL(Gmu0)[0];
  double dGsig20 = REAL(Gsig20)[0];
  int iUseA = INTEGER_VALUE(useA);
  int iKtype = INTEGER_VALUE(Ktype);

  double *S, *S0, *S1;

  GetRNGstate();

  // This is the high precision value for random effects - the reciprocal is the noninformative variance
  //double precHigh=1000000;
  SEXP predVals;
  double *pPredVals;
  if(iStorePred){
	PROTECT(predVals = alloc3Darray(REALSXP,designDims[0],2,saveIters));
	pPredVals = REAL(predVals);
  }else{
	PROTECT(predVals = allocVector(REALSXP,1));
	pPredVals  = Calloc(2*designDims[0],double); 
  }
  
  PROTECT(effectPars = allocMatrix(REALSXP, nEffectPars, saveIters));
  PROTECT(logLikes   = allocVector(REALSXP, iIters));
  pLogLikes = REAL(logLikes);

    //printf("Initialize Matrix stuff\n");

    int *sizeCovMat = INTEGER_POINTER(sizeCovMatR);
    int *obsCovMat  = INTEGER_POINTER(obsCovMatR);
    int maxSizeCovMat=0;
    int *parStart = INTEGER_POINTER(parStartR);
    int *covEffSlope = INTEGER_POINTER(covEffSlopeR);
    

    int nSlopes[nCovMat];
    for(j=0;j<nCovMat;j++){
		if(sizeCovMat[j]>maxSizeCovMat) maxSizeCovMat = sizeCovMat[j];		
		nSlopes[j]=0;
		for(k=0;k<sizeCovMat[j];k++)
			if(covEffSlope[k*nCovMat + j]) nSlopes[j]++; 
	}
	
    nMeanPars = nCovMat*maxSizeCovMat;  
  
    double *means = Calloc(nMeanPars, double);
    

    //initialize all means to zero
    //AZERO(means,nMeanPars);

    PROTECT(returnList = allocVector(VECSXP, 5));
    PROTECT(covMatsList    = allocVector(VECSXP,nCovMat));
    PROTECT(meanPars   = allocMatrix(REALSXP, nMeanPars, saveIters));
    

    //initialize starting values
 
    double *pMeanPars   = REAL(meanPars);
    double *pointerCovMat[nCovMat];
    double *pointerChainCovMat[nCovMat];


    SEXP pCovMats[nCovMat];
    SEXP chainCovMats[nCovMat];

    SEXP R_fcall,sampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
    PROTECT(sampCounter = NEW_INTEGER(1));
    int iProgress = INTEGER_VALUE(progress);
    int *pSampCounter = INTEGER_POINTER(sampCounter);

    //printf("Protecting Cov matrices\n");

    for(j=0;j<nCovMat;j++)
      {
	PROTECT(pCovMats[j]=allocMatrix(REALSXP, sizeCovMat[j],sizeCovMat[j]));
	PROTECT(chainCovMats[j] = alloc3Darray(REALSXP,sizeCovMat[j],sizeCovMat[j],saveIters));
	pointerCovMat[j]=REAL(pCovMats[j]);
	pointerChainCovMat[j]=REAL(chainCovMats[j]);
	// Initialize all covariance matrices to the identity matrix
	AZERO(pointerCovMat[j],sizeCovMat[j]*sizeCovMat[j]);	
	for(m=0;m<sizeCovMat[j];m++){
	  *(pointerCovMat[j] + m + m*sizeCovMat[j]) = 1;
	}
      }  
    
    
    //End Covariance matrix stuff

   double *pEffectPars = REAL(effectPars);

  //printf("Starting Gibbs Sampler\n");

  // start Gibbs sampler  
  for(i=0;i<iIters;i++){
  R_CheckUserInterrupt();

    //Check to see if we need to update the progress bar
    if(iProgress && !((i+1)%iProgress)){
      pSampCounter[0]=i+1;
      SETCADR(R_fcall, sampCounter);
      eval(R_fcall, rho); //Update the progress bar
    }
    //printf("Starting iteration %d...\n",i);

    //Sample effect parameters
    if(iStorePred){
      if(iMetrop){
      	WM2_metrop(params, nEffectPars, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, nCovMat, obsCovMat, sizeCovMat, 
		   parStart, means, pointerCovMat, &(pLogLikes[i]), maxSizeCovMat, &(pPredVals[i*2*designDims[0]]), dMetropSD, &logPostCurrEff);
      }else{
	WM2_hybridMC(params, nEffectPars, epsLo, epsRng, LFsteps, compWeights, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, nCovMat, obsCovMat, sizeCovMat, 
		     parStart, means, pointerCovMat, &(pLogLikes[i]), maxSizeCovMat, &(pPredVals[i*2*designDims[0]]),&logPostCurrEff);
      }    
    }else{
      if(iMetrop){
	WM2_metrop(params, nEffectPars, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, nCovMat, obsCovMat, sizeCovMat, 
		   parStart, means, pointerCovMat, &(pLogLikes[i]), maxSizeCovMat, pPredVals,dMetropSD,&logPostCurrEff);
      }else{
	WM2_hybridMC(params, nEffectPars, epsLo, epsRng, LFsteps, compWeights, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, nCovMat, obsCovMat, sizeCovMat, 
		     parStart, means, pointerCovMat, &(pLogLikes[i]), maxSizeCovMat, pPredVals,&logPostCurrEff);
      }
    }

    //Copy to chain
    //printf("Copy effect pars\n");
    
    if( !( (i+1)%iMetropThin ))
	Memcpy(pEffectPars + (((i+1)/iMetropThin)-1)*nEffectPars, params, nEffectPars);
    
    
    //printf("Sampling Cov Matrix\n");
    //Sample covariance matrices
      for(j=0;j<nCovMat;j++)
	{
	  //printf("Init Cov matrix %d\n",j);
	  S = pointerCovMat[j];
	  
	  p=sizeCovMat[j];
	  n=obsCovMat[j];
	  AZERO(S,p*p);
	  
	  for(k=0;k<p;k++)
	  {
	    S[k + p*k]=1;
	  }
	  
	  
	  for(m=0;m<n;m++) //iterate through observations
	  {
		for(k=0;k<p;k++) //iterate through parameters
		{
			for(q=0;q<=k;q++){
				//printf("k:%d q:%d   -   %f %f\n",k,q,params[parStart[k*nCovMat + j] + m],params[parStart[q*nCovMat + j] + m]);
				tmp = (params[parStart[k*nCovMat + j] + m]-means[k*nCovMat + j])*(params[parStart[q*nCovMat + j] + m]-means[q*nCovMat + j]);
				S[q + p*k]+=tmp;
				if(k!=q) S[k + p*q]+=tmp;
			}
		}
	  }

	//printf("%f %f %f %f\n",S[0],S[1],S[2],S[3]);

	//printf("Sample Wishart\n");
	//Sample Wishart
	InvMatrixUpper(S,p);
	internal_symmetrize(S,p);

	//printf("dfWishPrior: %f, n: %d\n",dfWishPrior,n);
	WM2_rWishartC(dfWishPrior + n, S, p, 0);
	//Copy to chain?
	//printf("copy to chain\n");
	
	//printf("GS: COV:  %f %f %f %f\n",S[0],S[1],S[2],S[3]);

	
	if( !( (i+1)%iMetropThin ))    
	  Memcpy(pointerChainCovMat[j] + (((i+1)/iMetropThin)-1)*p*p, S, p*p);
	}
	
	
    //Sample means
    for(j=0;j<nCovMat;j++)
	{
	  //printf("Init  mean sample %d\n",j);
	   if(nSlopes[j]>1){	
			S = pointerCovMat[j];
			p=sizeCovMat[j];
			n=obsCovMat[j];
			pSub = nSlopes[j];
		
			S1 = Calloc(p*p,double);
			S0 = Calloc(pSub*pSub,double);
			tmpmeans = Calloc(pSub,double);
			
			Memcpy(S1,S,p*p);
			InvMatrixUpper(S1,p);
			internal_symmetrize(S1,p);
			
			//print S matrix
			/*for(k=0;k<p;k++)
			  {
			    for(m=0;m<p;m++)
			      {
				printf("%f ",S[m+k*p]);
			      }
			    printf("\n");
			  }
			*/

			
			counter=0;
			for(k=0;k<p;k++)
			{  
			  //printf("Is k%d of %d a slope? %d\n",k,p,covEffSlope[k*nCovMat + j]);
			  if(covEffSlope[k*nCovMat + j])
			    {
			      S0[counter+pSub*counter] = S1[k + p*k]/n;
			      counter2=counter;
			      if(k!=p){
				for(m=(k+1);m<p;m++){
				  if(covEffSlope[m*nCovMat + j]){
				    //printf("c1:%d, c2:%d\n",counter,counter2);
				    S0[++counter2+pSub*counter]=S1[k+p*m]/n;
				    S0[counter + pSub*counter2]=S1[k+p*m]/n;
				  }
				}
				}
			      counter++;
			      //printf("counter: %d\n",counter);
			    }
			}
			//printf("i:%d %f %f %f %f\n",i,S0[0],S0[1],S0[2],S0[3]);
			//error("Stopping.\n");
			//InvMatrixUpper(S0,p);
			//internal_symmetrize(S0,p);
			
			for(m=0;m<n;m++) //iterate through observations
			{
			  counter=0;
			  for(k=0;k<p;k++) //iterate through parameters
			    {
			      if(covEffSlope[k*nCovMat + j]) tmpmeans[counter++]+=params[parStart[k*nCovMat + j]+m]/n;
			    }
			}
			//Sample Normal
			//printf("Sample Normal\n");

			WM2_rmvGaussianC(tmpmeans,S0,pSub);
			counter=0;
			//printf("Copy\n");
			for(k=0;k<p;k++)
			{
				if(covEffSlope[k*nCovMat + j])
				{
					means[k*nCovMat + j]=tmpmeans[counter++];
				}else
				{
					means[k*nCovMat + j]=0;
				}
			}
			//printf("Free memory...\n");
			Free(S0);
			Free(tmpmeans);
			Free(S1);
	   }else if(nSlopes[j]==1){	
	     S = pointerCovMat[j];
	     
		 p=sizeCovMat[j];
	     n=obsCovMat[j];
		 S1 = Calloc(p*p,double);
		 Memcpy(S1,S,p*p);
		 InvMatrixUpper(S1,p);
		 internal_symmetrize(S1,p);
		 
	     tmpmean=0;
	     
		 
		 
		 for(k=0;k<p;k++)
	       {
		 if(covEffSlope[k*nCovMat + j])
		   {
		     for(m=0;m<n;m++) //iterate through observations
		       tmpmean+=params[parStart[k*nCovMat + j]+m]/n;
		     means[k*nCovMat + j]=rnorm(tmpmean,sqrt(S1[k+k*p]/n));
		   }else
		   {
		     means[k*nCovMat + j]=0;
		   }
	       }
	    Free(S1);
	   }else{
	     p=sizeCovMat[j];
	     for(k=0;k<p;k++) means[k*nCovMat + j]=0;
	   }
	}
    //printf("Copy to chain\n");
	 //copy to chain
	 
	 if( !( (i+1)%iMetropThin ))
		Memcpy(pMeanPars + (((i+1)/iMetropThin)-1)*nMeanPars,means,nMeanPars);
	 
	 //printf("Done copying.\n");
  
  }  
  

  //printf("Prepare to return...\n");
    //Prepare things to return

    SET_VECTOR_ELT(returnList, 0, effectPars);
    SET_VECTOR_ELT(returnList, 1, logLikes);
  
    for(j=0;j<nCovMat;j++)
      {
	SET_VECTOR_ELT(covMatsList, j, chainCovMats[j]);         
      }  
  
     SET_VECTOR_ELT(returnList, 2, meanPars);       
     SET_VECTOR_ELT(returnList, 3, covMatsList);       
     SET_VECTOR_ELT(returnList, 4, predVals);

     UNPROTECT(8+2*nCovMat);
     
     Free(means);
	 if(!iStorePred) Free(pPredVals); 
  
	 PutRNGstate();
      return(returnList);
}
 

SEXP WM2_GibbsSamplerNoCov(SEXP iters, SEXP starteffects, SEXP nHit, SEXP nMiss, SEXP nFA,
		      SEXP nCR,  SEXP setSize, SEXP design, SEXP designCont, SEXP Keffects, SEXP KeffectsSlope, SEXP KeffectsCov, 
		      SEXP Aeffects, SEXP AeffectsSlope, SEXP AeffectsCov, SEXP Geffects, SEXP GeffectsSlope, SEXP GeffectsCov, 
		      SEXP invGammaPriorA, SEXP invGammaPriorB, SEXP Kmu0, SEXP Ksig20, SEXP Amu0, SEXP Asig20, SEXP Gmu0, SEXP Gsig20, 
			   SEXP useA, SEXP Ktype, SEXP epsLoR, SEXP epsRngR, SEXP LFstepsR, SEXP compWeightsR, SEXP progress, SEXP pBar, SEXP rho, SEXP storePred, SEXP metrop, SEXP metropSD, SEXP metropThin)
{
  //printf("Start C call.\n");
  int iIters = INTEGER_VALUE(iters), i=0, nEffectPars = length(starteffects);
  int iMetrop = INTEGER_VALUE(metrop), iMetropThin=INTEGER_VALUE(metropThin); 
  
  if(!iMetrop) iMetropThin = 1;
  if(iMetropThin>iIters) error("Thinning interval larger than number of iterations.");
  int saveIters = iIters/iMetropThin;

  double *dMetropSD = REAL(metropSD);
  SEXP effectPars;
  SEXP returnList, logLikes;
  SEXP sampCounter, R_fcall;
  int *pSampCounter, iProgress;
  double *pLogLikes;
  int iStorePred = INTEGER_VALUE(storePred);
  
  int *dims = INTEGER(getAttrib(design, R_DimSymbol));
  int designDims[2];
  designDims[0]=dims[0];
  designDims[1]=dims[1];
  //printf("dims1: %d, dims2: %d\n",dims[0],dims[1]);


  int *pHit = INTEGER_POINTER(nHit);
  int *pMiss = INTEGER_POINTER(nMiss);
  int *pFA = INTEGER_POINTER(nFA);
  int *pCR = INTEGER_POINTER(nCR);


  int *pSetSize = INTEGER_POINTER(setSize);
  int *pDesign = INTEGER_POINTER(design);
  double *pDesignCont = REAL(designCont);
  double logPostCurrEff = 0;

  int *pKeffects = INTEGER_POINTER(Keffects);
  int *pKeffectsCov = INTEGER_POINTER(KeffectsCov);
  int *pKeffectsSlope = INTEGER_POINTER(KeffectsSlope);
  int *pAeffects = INTEGER_POINTER(Aeffects);
  int *pAeffectsCov = INTEGER_POINTER(AeffectsCov);
  int *pAeffectsSlope = INTEGER_POINTER(AeffectsSlope);
  int *pGeffects = INTEGER_POINTER(Geffects);
  int *pGeffectsCov = INTEGER_POINTER(GeffectsCov);
  int *pGeffectsSlope = INTEGER_POINTER(GeffectsSlope);


  double params[nEffectPars];
  Memcpy(params,REAL(starteffects),nEffectPars);
 
  double epsLo = REAL(epsLoR)[0];
  double epsRng = REAL(epsRngR)[0];
  double *compWeights = REAL(compWeightsR);
  int LFsteps = INTEGER_VALUE(LFstepsR);

  double a0 = REAL(invGammaPriorA)[0];
  double b0 = REAL(invGammaPriorB)[0];
  double dKmu0 = REAL(Kmu0)[0];
  double dKsig20 = REAL(Ksig20)[0];
  double dAmu0 = REAL(Amu0)[0];
  double dAsig20 = REAL(Asig20)[0];
  double dGmu0 = REAL(Gmu0)[0];
  double dGsig20 = REAL(Gsig20)[0];
  int iUseA = INTEGER_VALUE(useA);
  int iKtype = INTEGER_VALUE(Ktype);

  
  
  SEXP predVals;
  double *pPredVals;
  if(iStorePred){
	PROTECT(predVals = alloc3Darray(REALSXP,designDims[0],2,iIters));
	pPredVals = REAL(predVals);
  }else{
	PROTECT(predVals = allocVector(REALSXP,1));
	pPredVals  = Calloc(2*designDims[0],double); 
  }
  
  PROTECT(R_fcall = lang2(pBar, R_NilValue));
  PROTECT(sampCounter = NEW_INTEGER(1));
  PROTECT(effectPars = allocMatrix(REALSXP, nEffectPars, saveIters));
  PROTECT(logLikes   = allocVector(REALSXP, iIters));
  PROTECT(returnList = allocVector(VECSXP, 3));
  pLogLikes = REAL(logLikes);
  iProgress = INTEGER_VALUE(progress);
  pSampCounter = INTEGER_POINTER(sampCounter);

   double *pEffectPars = REAL(effectPars);

   //printf("Starting Gibbs Sampler\n");


   // Get Current Log Posterior
   logPostCurrEff = LogPosterior(params, nEffectPars, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, 0, (int *)(0), (int *)(0), (int *)(0), (double *)(0), (double **)(0), &(pLogLikes[i]), 0, pPredVals);
 

  // start Gibbs sampler  
  for(i=0;i<iIters;i++){
  R_CheckUserInterrupt();
  
  //printf("Iteration %d\n",i);
  //Check to see if we need to update the progress bar
  if(iProgress && !((i+1)%iProgress)){
    pSampCounter[0]=i+1;
    SETCADR(R_fcall, sampCounter);
    eval(R_fcall, rho); //Update the progress bar
  }

    //Sample effect parameters
    if(iStorePred){
      if(iMetrop){
	WM2_metrop(params, nEffectPars, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, 0, (int *)(0), (int *)(0), (int *)(0), (double *)(0), (double **)(0), &(pLogLikes[i]), 0, &(pPredVals[i*2*designDims[0]]), dMetropSD, &logPostCurrEff);
      }else{
	WM2_hybridMC(params, nEffectPars, epsLo, epsRng, LFsteps, compWeights, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, 0, (int *)(0), (int *)(0), 
		     (int *)(0), (double *)(0), (double **)(0), &(pLogLikes[i]), 0, &(pPredVals[i*2*designDims[0]]),&logPostCurrEff);
      }
    }else{
      if(iMetrop){
      	WM2_metrop(params, nEffectPars, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, 0, (int *)(0), (int *)(0), 
		   (int *)(0), (double *)(0), (double **)(0), &(pLogLikes[i]), 0, pPredVals, dMetropSD,&logPostCurrEff);
      }else{
	WM2_hybridMC(params, nEffectPars, epsLo, epsRng, LFsteps, compWeights, pHit, pMiss, pFA, pCR, pSetSize, pDesign, pDesignCont, 
		 designDims, pKeffects, pKeffectsSlope, pKeffectsCov, pAeffects, pAeffectsSlope, pAeffectsCov, pGeffects, pGeffectsSlope, 
		 pGeffectsCov, a0, b0, dKmu0, dKsig20, dAmu0, dAsig20, dGmu0, dGsig20, iUseA, iKtype, 0, (int *)(0), (int *)(0), 
		     (int *)(0), (double *)(0), (double **)(0), &(pLogLikes[i]), 0, pPredVals,&logPostCurrEff);
      }
    }
    //printf("Copy to chain next\n");

    //Copy to chain
    if( !( (i+1)%iMetropThin ))
      Memcpy(pEffectPars + (((i+1)/iMetropThin)-1)*nEffectPars, params, nEffectPars);

    if(i<(iIters-1)) pLogLikes[i+1]=pLogLikes[i];
  }
  
  //printf("Preparing to return...\n");
   
    //Prepare things to return

    SET_VECTOR_ELT(returnList, 0, effectPars);
    SET_VECTOR_ELT(returnList, 1, logLikes);
    SET_VECTOR_ELT(returnList, 2, predVals);

    UNPROTECT(6);
    
    //printf("Returning.\n");
	if(!iStorePred) Free(pPredVals); 
    return(returnList);
}
 
