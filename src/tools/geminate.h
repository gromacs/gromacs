#ifndef _GEMINATE_H
#define _GEMINATE_H

enum { gemNULL, gemNONE, gemDD, gemAD, gemAA, gemA4, gemNR};
static char *gemType[] = {NULL, "none", "dd", "ad", "aa", "a4", NULL};


/* This first part is derived from complex.c which I recieved from Omer Markowitch.
 * - Erik Marklund
 *
 * ------------- from complex.c ------------- */

#include <math.h>
/* definition of PI */
#ifndef PI
#define PI (acos(-1.0))
#endif

/* definition of the type complex */
typedef struct
{
  double r,i;
} complex;

/* Complexation of a paired number (r,i)                                     */
complex cmplx(double r, double i);

/* Complexation of a real number, x */
static complex _c(double x);

/* Real and Imaginary part of a complex number z -- Re (z) and Im (z)        */
static double Re(complex z);
static double Im(complex z);

/* Magnitude of a complex number z                                           */
static double cx_abs(complex z);

/* Addition of two complex numbers z1 and z2                                 */
static complex cxadd(complex z1, complex z2);

/* Addition of a complex number z1 and a real number r */
static complex cxradd(complex z, double r);

/* Subtraction of two complex numbers z1 and z2                              */
static complex cxsub(complex z1, complex z2);

/* Multiplication of two complex numbers z1 and z2                           */
static complex cxmul(complex z1, complex z2);

/* Square of a complex number z                                              */
static complex cxsq(complex z);

/* multiplication of a complex number z and a real number r */
static complex cxrmul(complex z, double r);

/* Division of two complex numbers z1 and z2                                 */
static complex cxdiv(complex z1, complex z2);

/* division of a complex z number by a real number x */
static complex cxrdiv(complex z, double r);

/* division of a real number r by a complex number x */
static complex rxcdiv(double r, complex z);

/* Integer power of a complex number z -- z^x                                */
static complex cxintpow(complex z, int x);

/* Exponential of a complex number-- exp (z)=|exp(z.r)|*{cos(z.i)+I*sin(z.i)}*/
static complex cxdexp(complex z);

/* Logarithm of a complex number -- log(z)=log|z|+I*Arg(z),                  */
/*                                  where -PI < Arg(z) < PI                  */ 
static complex cxlog(complex z);

/* Square root of a complex number z = |z| exp(I*the) -- z^(1/2)             */
/*                               z^(1/2)=|z|^(1/2)*[cos(the/2)+I*sin(the/2)] */
/*                               where 0 < the < 2*PI                        */
static complex cxdsqrt(complex z);

/* square root of a real number r  */
static complex cxrsqrt(double r);

/* Complex power of a complex number z1^z2                                   */
static complex cxdpow(complex z1, complex z2);

/* Print out a complex number z as z: z.r, z.i                               */
static void cxprintf(complex z);

/* ------------ end of complex.c ------------ */


/* This next part was derived from cerror.c and rerror.c,
 * also received from Omer Markovitch.
 * ------------- from [cr]error.c ------------- */

#ifndef sPI
#define sPI (sqrt(PI))
#endif
/************************************************************/
/* Real valued error function and related functions         */
/* x, y     : real variables                                */
/* erf(x)   : error function                                */
/* erfc(x)  : complementary error function                  */
/* omega(x) : exp(x*x)*erfc(x)                              */
/* W(x,y)   : exp(-x*x)*omega(x+y)=exp(2*x*y+y^2)*erfc(x+y) */
/************************************************************/

/*---------------------------------------------------------------------------*/
/* Utilzed the series approximation of erf(x)                                */
/* Relative error=|err(x)|/erf(x)<EPS                                        */
/* Handbook of Mathematical functions, Abramowitz, p 297                     */
/* Note: When x>=6 series sum deteriorates -> Asymptotic series used instead */
/*---------------------------------------------------------------------------*/
/*  This gives erfc(z) correctly only upto >10-15 */
static double _erf(double x);

/* Result --> Alex's code for erfc and experfc behaves better than this      */
/* complementray error function.                Returns 1.-erf(x)            */
static double _erfc(double x);

/* omega(x)=exp(x^2)erfc(x)                                                  */
static double omega(double x);

/* W(x,y)=exp(-x^2)*omega(x+y)=exp(2xy+y^2)*erfc(x+y)                        */
static double W(double x, double y);

/**************************************************************/
/* Complex error function and related functions               */
/* x, y     : real variables                                  */
/* z        : complex variable                                */
/* cerf(z)  : error function                                  */
/* comega(z): exp(z*z)*cerfc(z)                               */
/* W(x,z)   : exp(-x*x)*comega(x+z)=exp(2*x*z+z^2)*cerfc(x+z) */
/**************************************************************/
static complex cerf(complex z);

/*---------------------------------------------------------------------------*/
/* Utilzed the series approximation of erf(z=x+iy)                           */
/* Relative error=|err(z)|/|erf(z)|<EPS                                      */
/* Handbook of Mathematical functions, Abramowitz, p 299                     */
/* comega(z=x+iy)=cexp(z^2)*cerfc(z)                                         */
/*---------------------------------------------------------------------------*/
static complex comega(complex z);

/* W(x,z) exp(-x^2)*omega(x+z)                                               */
static complex cW(double x, complex z);

/* ------------ end of [cr]error.c ------------ */


/* ///////////////////////////////////////////////////////////////////////
 * Here follow routines and structs for reversible geminate recombination.
 */

typedef struct{
  size_t n;
  double *y;
  double tDelta;
  int nexp;
} balData;


typedef struct {
  /* Used in the rewritten version of Omer's gem. recomb. analysis */
  double ka, kd;                 /* -- || -- results[]  */
  double sigma;                  /* -- || -- sigma      */
  double D;                      /* The diffusion coefficient */
  double kD;                     /* Replaces kD in analyse_corr_gem3d() */

  /* The following members are for calcsquare() and takeAwayBallistic() */
  double tDelta;              /* Time between frames */
  double logAfterTime;        /* Time after which we do the lsq calculations on a logarithmic timescale. */
  double logDelta;
  double logPF;
  /* To get an equal number of points in the lin and log regime,
   * we'll use logDelta to calculate where to sample the ACF.
   * if i and j are indices in the linear and log regime, then:
   *   j = Cexp(A(i+nLin)),
   * where C = (nLin**2 / len) and A = log(len/nLin) / nLin.
   * This expands to j = (nLin**2 / len) * exp((i+nLin) * log(len/nLin) / nLin).
   * We'll save part of our brains and some computing time if we pre-calculate
   *  1) logDelta = log(len/nLin) / nLin
   *  2) logPF    = nLin**2 / len
   * and get j = logPF * exp((i+nLin) * logDelta). */
#define GETLOGINDEX(i,params) (params)->logPF * exp(((i)+(params)->nLin) * (params)->logDelta)
  int nLin;                 /* Number of timepoints in the linear regime */
  int len;                  /* Length of time and ct arrays */
  int nExpFit;              /* Number of exponentials to fit */       
  real ballistic;           /* Time before which the ballistic term should be fitted */
  bool bDt;                 /* TRUE =>  use time derivative at time 0
			     *          to find fastest component.
			     * FALSE => use coefficient in exponenetial
			     *          to find fastest component. */
} t_gemParams;


typedef struct{
  size_t n;         /* Number of data points (lin-log) */
  double *y;        /* The datapoints */
  double *ctTheory; /* The theoretical ct which will be built by gemFunc_f. */
  double *LinLog;
  double *time;
  double ka;
  double kd;
  double tDelta;    /* time difference between subsequent datapoints */
  size_t nData;     /* real size of the data */
  t_gemParams *params;
} gemFitData;

extern void takeAwayBallistic(double *ct, double *t,
		       int len, real tMax,
		       int nexp, bool bDerivative);


extern t_gemParams *init_gemParams(double sigma, double D,
				   real *t, double logAfterTime,
				   int len, real ballistic, int nBalExp, bool bDt);

/* Fit to geminate recombination model.
   Returns root mean square error of fit. */
extern real fitGemRecomb(double *ct, double *time, double **ctFit,
		  const int nData, t_gemParams *params);

#endif
