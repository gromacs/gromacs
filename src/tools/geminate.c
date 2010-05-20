#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBGSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_version.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "geminate.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
/* This first part is from complex.c which I recieved from Omer Markowitch.
 * - Erik Marklund
 *
 * ------------- from complex.c ------------- */

complex cmplx(double r, double i)
{
  complex value;
  value.r = r;
  value.i=i;
  return value;
}

static complex _c(double x)
{
  complex value;
  value.r=x;
  value.i=0.;
  return value;
}
static double Re(complex z) {return z.r;}
static double Im(complex z) {return z.i;}
static double cx_abs(complex z) { return (sqrt(z.r*z.r+z.i*z.i)); }
static complex cxadd(complex z1, complex z2)
{
  complex value;
  value.r=z1.r+z2.r;
  value.i=z1.i+z2.i;
  return value;
}
static complex cxradd(complex z, double r)
{
  complex value;
  value.r = z.r + r;
  value.i = z.i;
  return value;
}
static complex cxsub(complex z1, complex z2)
{
  complex value;
  value.r=z1.r-z2.r;
  value.i=z1.i-z2.i;
  return value;
}
static complex cxmul(complex z1, complex z2)
{
  complex value;
  value.r = z1.r*z2.r-z1.i*z2.i;
  value.i = z1.r*z2.i+z1.i*z2.r;
  return value;
}
static complex cxsq(complex z)
{
  complex value;
  value.r = z.r*z.r-z.i*z.i;
  value.i = z.r*z.i*2.;
  return value;
}
static complex cxrmul(complex z, double r)
{
  complex value;
  value.r = z.r*r;
  value.i= z.i*r;
  return value;
}
static complex cxdiv(complex z1, complex z2)
{
  complex value;
  double num;
  num = z2.r*z2.r+z2.i*z2.i;
  if(num == 0.) printf("ERROR in cxdiv function\n");
  value.r = (z1.r*z2.r+z1.i*z2.i)/num; value.i = (z1.i*z2.r-z1.r*z2.i)/num;
  return value;
}
static complex cxrdiv(complex z, double r)
{
  complex value;
  value.r = z.r/r;
  value.i = z.i/r;
  return value;
}
static complex rxcdiv(double r, complex z)
{
  complex value;
  double f;
  f = r/(z.r*z.r+z.i*z.i);
  value.r = f*z.r;
 value.i = -f*z.i;
 return value;
}
static complex cxintpow(complex z, int x)
{
  int i;
  complex value={1.,0.};
  if(x>0) {
    for(i=0; i < x; i++)
      value = cxmul(value, z);
    return value;
  } else if(x<0) {
    for(i=0; i > x; i--)
      value = cxdiv(value, z);
    return value;
  } else return value;
}
static complex cxdexp(complex z)
{
  complex value;
  double exp_z_r;
  exp_z_r = exp(z.r);
  value.r = exp_z_r*cos(z.i);
  value.i = exp_z_r*sin(z.i);
  return value;
}
static complex cxlog(complex z)
{
  complex value;
  double mag2;
  mag2 = z.r*z.r+z.i*z.i;
  if(mag2 < 0.)
    printf("ERROR in cxlog func\n");
  value.r = log(sqrt(mag2));
  if(z.r == 0.) {
    value.i = PI/2.;
    if(z.i <0.)
      value.i = -value.i;
  } else
    value.i = atan2(z.i, z.r);
  return value;
}
static complex cxdsqrt(complex z)
{
  complex value;
  double sq;
  sq = cx_abs(z);
  value.r=sqrt(fabs((sq+z.r)*0.5)); /* z'.r={|z|*[1+cos(the)]/2}^(1/2) */
  value.i=sqrt(fabs((sq-z.r)*0.5)); /* z'.i={|z|*[1-cos(the)]/2}^(1/2) */
  if(z.i < 0.)
    value.r = -value.r;
  return value;
}
static complex cxrsqrt(double r) {
  if (r < 0)
    return(cmplx(0, sqrt(-r)));
  else
    return(_c(sqrt(r)));
}
static complex cxdpow(complex z1, complex z2)
{
  complex value;
  value=cxdexp(cxmul(cxlog(z1),z2));
  return value;
}
static void cxprintf(complex z) { printf("z: %lg + %lg_i\n", z.r, z.i); }
/* ------------ end of complex.c ------------ */

/* This next part was derived from cubic.c, also received from Omer Markovitch.
 * ------------- from cubic.c ------------- */

/* Solver for a cubic equation: x^3-a*x^2+b*x-c=0                            */
static void solve(complex *al, complex *be, complex *gam,double a,double b,double c)
{
	double t1, t2, two_3, temp;
	complex ctemp, ct3;
	
	two_3=pow(2., 1./3.); t1 = -a*a+3.*b; t2 = 2.*a*a*a-9.*a*b+27.*c;
	temp = 4.*t1*t1*t1+t2*t2;

	ctemp = cmplx(temp,0.);	ctemp = cxadd(cmplx(t2,0.),cxdsqrt(ctemp));
	ct3 = cxdpow(ctemp, cmplx(1./3.,0.));

	ctemp = rxcdiv(-two_3*t1/3., ct3);
	ctemp = cxadd(ctemp, cxrdiv(ct3, 3.*two_3));
	
	*gam = cxadd(cmplx(a/3.,0.), ctemp);
	
	ctemp=cxmul(cxsq(*gam), cxsq(cxsub(*gam, cmplx(a,0.))));
	ctemp=cxadd(ctemp, cxmul(cmplx(-4.*c,0.), *gam));
	ctemp = cxdiv(cxdsqrt(ctemp), *gam);
	*al = cxrmul(cxsub(cxsub(cmplx(a,0.), *gam),ctemp),0.5);
	*be = cxrmul(cxadd(cxsub(cmplx(a,0.), *gam),ctemp),0.5);	
}

/* ------------ end of cubic.c ------------ */

/* This next part was derived from cerror.c and rerror.c, also received from Omer Markovitch.
 * ------------- from [cr]error.c ------------- */

static double _erf(double x)
{
  double n,sum,temp,exp2,xm,x2,x4,x6,x8,x10,x12;
  temp=x;
  sum=temp;
  xm=26.;
  x2=x*x;
  x4=x2*x2;
  x6=x4*x2;
  x8=x6*x2;
  x10=x8*x2;
  x12=x10*x2;
  exp2=exp(-x2);
  if(x<=xm){
    for(n=1.; n<=2000.; n+=1.){
      temp *= 2.*x2/(2.*n+1.);
      sum  += temp;
      if(fabs(temp/sum)<1.E-16)
	break;
    }
    if(n>=2000.)
      printf("In Erf calc - iteration exceeds %lg\n",n);
    sum *= 2./sPI*exp2;
  } else {
    /* from the asymptotic expansion of experfc(x) */
    sum = (1. - 0.5/x2 + 0.75/x4 - 1.875/x6 + 6.5625/x8 - 29.53125/x10 + 162.421875/x12) / sPI/x;
    sum *= exp2; /* now sum is erfc(x) */
    sum=-sum+1.;
  }
  return x>=0.0 ? sum : -sum;
}

static double _erfc(double x)
{
  double t,z,ans;
  z = fabs(x);
  t = 1.0/(1.0+0.5*z);
  
  ans = t * exp ( -z*z - 1.26551223 +
		  t*(1.00002368 +
		     t*(0.37409196 +
			t*(0.09678418 +
			   t*(-0.18628806 +
			      t*(0.27886807 +
				 t*(-1.13520398 +
				    t*(1.48851587 +
				       t*(-0.82215223 +
					  t*0.17087277)))))))));
	
  return  x>=0.0 ? ans : 2.0-ans;
}

static double omega(double x)
{
  double xm=26., ans;
  double xx=x*x, x4=xx*xx, x6=x4*xx, x8=x6*xx, x10=x8*xx, x12=x10*xx;
  if(x<=xm)
    ans = exp(xx)*_erfc(x);
  else {
    /* Asymptotic expansion */
    ans = (1. - 0.5/xx + 0.75/x4 - 1.875/x6 + 6.5625/x8 - 29.53125/x10 + 162.421875/x12) / sPI/x;
  }
  return ans;
}

static double W(double x, double y){ return(exp(-x*x)*omega(x+y)); }

static complex cerf(complex z)
{
  complex value;
  double x,y;
  double sumr,sumi,n,n2,f,temp,temp1;
  double x2,cos_2xy,sin_2xy,cosh_2xy,sinh_2xy,cosh_ny,sinh_ny;

  x=z.r;
  y=z.i;
  x2=x*x;
  sumr=0.;
  sumi=0.;
  cos_2xy=cos(2.*x*y);
  sin_2xy=sin(2.*x*y);
  cosh_2xy=cosh(2.*x*y);
  sinh_2xy=sinh(2.*x*y);

  for(n=1.0,temp=0.; n<=2000.; n+=1.0)
    {
      n2=n*n;
      cosh_ny=cosh(n*y);
      sinh_ny=sinh(n*y);
      f=exp(-n2/4.)/(n2+4.*x2);
      sumr += (2.*x - 2.*x*cosh_ny*cos_2xy+n*sinh_ny*sin_2xy)*f;
      sumi += (2.*x*cosh_ny*sin_2xy + n*sinh_ny*cos_2xy)*f;
      temp1 = sqrt(sumr*sumr+sumi*sumi);
      if(fabs((temp1-temp)/temp1)<1.E-16)
	break;	
      temp = temp1;
    }
  if(n==2000.)
    printf("iteration exceeds %lg\n",n);
  sumr*=2./PI;
  sumi*=2./PI;

  if(x!=0.)
    f=1./2./PI/x;
  else
    f=0.;
  value.r = _erf(x) + (f*(1.-cos_2xy) + sumr)*exp(-x2);
  value.i = (f*sin_2xy+sumi)*exp(-x2);
  return value;
}

static complex comega(complex z)
{
  complex value;
  double x,y;
  double sumr,sumi,n,n2,f,temp,temp1;
  double x2, y2,cos_2xy,sin_2xy,cosh_2xy,sinh_2xy,cosh_ny,sinh_ny,exp_y2;

  x=z.r;
  y=z.i;
  x2=x*x;
  y2=y*y;
  sumr=0.;
  sumi=0.;
  cos_2xy=cos(2.*x*y);
  sin_2xy=sin(2.*x*y);
  cosh_2xy=cosh(2.*x*y);
  sinh_2xy=sinh(2.*x*y);
  exp_y2=exp(-y2);

  for(n=1.0, temp=0.; n<=2000.; n+=1.0){
    n2 = n*n;
    cosh_ny = cosh(n*y);
    sinh_ny = sinh(n*y);
    f = exp(-n2/4.)/(n2+4.*x2);
    /* if(f<1.E-200) break; */
    sumr += (2.*x - 2.*x*cosh_ny*cos_2xy + n*sinh_ny*sin_2xy)*f;
    sumi += (2.*x*cosh_ny*sin_2xy + n*sinh_ny*cos_2xy)*f;
    temp1 = sqrt(sumr*sumr+sumi*sumi);
    if(fabs((temp1-temp)/temp1)<1.E-16)
      break;
    temp=temp1;
  }
  if(n==2000.)
    printf("iteration exceeds %lg\n",n);
  sumr*=2./PI;
  sumi*=2./PI;

  if(x!=0.)
    f=1./2./PI/x;
  else
    f=0.;
  value.r = omega(x)-(f*(1.-cos_2xy)+sumr);
  value.i = -(f*sin_2xy+sumi);
  value   = cxmul(value,cmplx(exp_y2*cos_2xy,exp_y2*sin_2xy));
  return (value);
}

static complex cW(double x, complex z){ return(cxrmul(comega(cxradd(z,x)),exp(-x*x))); }

/* ------------ end of [cr]error.c ------------ */


/* Changes the unit from square cm per s to square Ångström per ps,
 * since Omers code uses the latter units while g_mds outputs the former.
 * g_hbond expects a diffusion coefficent given in square cm per s. */
static double sqcm_per_s_to_sqA_per_ps (real D) {
  fprintf(stderr, "Diffusion coefficient is %f A^2/ps\n", D*1e4);
  return (double)(D*1e4);
}


static double eq10v2(double theoryCt[], double time[], int manytimes,
		     double ka, double kd, t_gemParams *params)
{
  //Finding the 3 roots
  double
    kD = params->kD,
    D  = params->D,
    r = params->sigma,

    a = (1.0 + ka/kD) * sqrt(D)/r,
    b = kd,
    c = kd * sqrt(D)/r;
  complex alpha, beta, gamma;
  solve(&alpha, &beta, &gamma, a, b, c);
  //Finding the 3 roots

  int i;
  complex c1, c2, c3, c4, oma, omb, omc, part1, part2, part3, part4;
  double tsqrt, sumimaginary=0.0;

  part1 = cxmul(alpha, cxmul(cxadd(beta,  gamma), cxsub(beta,  gamma))); //1(2+3)(2-3)
  part2 = cxmul(beta,  cxmul(cxadd(gamma, alpha), cxsub(gamma, alpha))); //2(3+1)(3-1)
  part3 = cxmul(gamma, cxmul(cxadd(alpha, beta) , cxsub(alpha, beta)));  //3(1+2)(1-2)
  part4 = cxmul(cxsub(gamma, alpha), cxmul(cxsub(alpha, beta), cxsub(beta, gamma))); //(3-1)(1-2)(2-3)

#ifdef HAVE_OPENMP
#pragma omp parallel for				\
  private(i, tsqrt, oma, omb, omc, c1, c2, c3, c4),	\
  reduction(+:sumimaginary),				\
  default(shared),					\
  schedule(guided)
  for (i=0; i<manytimes; i++){
    tsqrt = sqrt(time[i]);
    oma   = comega(cxrmul(alpha, tsqrt));
    c1    = cxmul(oma, cxdiv(part1, part4));
    omb   = comega(cxrmul(beta, tsqrt));
    c2    = cxmul(omb, cxdiv(part2, part4));
    omc   = comega(cxrmul(gamma, tsqrt));
    c3    = cxmul(omc, cxdiv(part3, part4));
    c4.r  = c1.r+c2.r+c3.r;
    c4.i  = c1.i+c2.i+c3.i;
    theoryCt[i]  = c4.r;
    sumimaginary += c4.i * c4.i;
  }
#endif

  return sumimaginary;

} //eq10v2


extern t_gemParams *init_gemParams(double sigma, double D,
				   real *t, double logAfterTime,
				   int len, real ballistic, int nBalExp, bool bDt)
{
  double tDelta;
  t_gemParams *p;
  snew(p,1);

  /* A few hardcoded things here. For now anyway. */
/*   p->ka_min   = 0; */
/*   p->ka_max   = 100; */
/*   p->dka      = 10; */
/*   p->kd_min   = 0; */
/*   p->kd_max   = 2; */
/*   p->dkd      = 0.2; */
  p->ka       = 0;
  p->kd       = 0;
/*   p->lsq      = -1; */
/*   p->lifetime = 0; */
  p->sigma    = sigma;
/*   p->lsq_old  = 99999; */
  p->D        = sqcm_per_s_to_sqA_per_ps(D);
  p->kD       = 4*acos(-1.0)*sigma*p->D;


  /* Parameters used by calcsquare(). Better to calculate them
   * here than in calcsquare every time it's called. */
  p->len = len;
  p->logAfterTime = logAfterTime;
  tDelta       = (t[len-1]-t[0]) / len;
  if (tDelta <= 0)
    gmx_fatal(FARGS, "Time between frames is non-positive!");
  else p->tDelta = tDelta;

  p->nExpFit      = nBalExp;
  p->nLin         = logAfterTime / tDelta;
/*   if (p->nLin <= 0) { */
/*     fprintf(stderr, "Number of data points in the linear regime is non-positive!\n"); */
/*     sfree(p); */
/*     return NULL; */
/*   } */
  /* We want the same number of data points in the log regime. Not crucial, but seems like a good idea. */
  p->logDelta = log(((float)len)/p->nLin) / p->nLin;
  p->logPF    = p->nLin*p->nLin/(float)len;
  
/*   p->logMult      = pow((float)len, 1.0/nLin);/\* pow(t[len-1]-t[0], 1.0/p->nLin); *\/ */
  p->ballistic    =  ballistic;
  p->bDt;
  return p;
}

#ifdef HAVE_LIBGSL
static double gemFunc_residual2(const gsl_vector *p, void *data)
{
  gemFitData *GD = (gemFitData *)data;
  int i,iLog,
    nLin=GD->params->nLin,
    nData=GD->nData;
  double r, residual2=0,
    *ctTheory=GD->ctTheory,
    *y=GD->y;

  eq10v2(GD->ctTheory, GD->time, GD->nData,
    	 gsl_vector_get(p, 0), gsl_vector_get(p, 1),
    	 GD->params);
  
  /* Removing a bunch of points from the log-part. */
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic)	\
  firstprivate(nLin, nData, ctTheory, y),	\
  private (i, iLog, r),				\
  reduction(+:residual2),			\
  default(shared)
#endif
  for (i=0; i<nLin; i++) {
    /* Linear part ----------*/
    r = ctTheory[i];
    residual2 += sqr(r-y[i]);
    /* Log part -------------*/
    iLog = GETLOGINDEX(i, GD->params);
    if (iLog >= nData)
      gmx_fatal(FARGS, "log index out of bounds: %i", iLog);
    r = ctTheory[iLog];
    residual2 += sqr(r-y[iLog]);

  }
  residual2 /= GD->n;
  return residual2;
}
#endif
 
extern real fitGemRecomb(double *ct, double *time, double **ctFit,
			const int nData, t_gemParams *params)
{
#ifndef HAVE_LIBGSL
  printf("Sorry, can't do reversible geminate recombination without gsl. "
	 "Recompile using --with-gsl.\n");
  return -1;
#else
  printf("Will fit ka and kd to the ACF according to the reversible geminate recombination model.\n");

  real   size, d2, tol=1e-10;
  int i,
    iter    = 0,
    status  = 0,
    maxiter = 1000;

  gsl_multimin_fminimizer *s;
  gsl_vector *x,*dx; /* parameters and initial step size */
  gsl_multimin_function fitFunc;
  size_t
    p = 2,              /* Number of parameters to fit. ka and kd.  */
    n = params->nLin*2;       /* Number of points in the reduced dataset  */

  if (params->D <= 0)
    {
      printf("Fitting of D is not implemented yet. It must be provided on the command line.\n");
      return -1;
    }
  /* nmsimplex2 had convergence problems prior to gsl v1.14,
   * but it's O(N) instead of O(N) operations, so let's use it if v >= 1.14 */
  gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
#ifdef GSL_MAJOR_VERSION
#ifdef GSL_MINOR_VERSION
  if ((GSL_MAJOR_VERSION == 1 && GSL_MINOR_VERSION >= 14) || 
      (GSL_MAJOR_VERSION > 1))
    T = gsl_multimin_fminimizer_nmsimplex2;
#endif
#endif
  
  if (nData<n) {
    fprintf(stderr, "Reduced data set larger than the complete data set!\n");
    n=nData;
  }
  gemFitData *GD;
  snew(GD,1);

  GD->n = n;
  GD->y = ct;
  GD->ctTheory=NULL;
  snew(GD->ctTheory, nData);
  GD->LinLog=NULL;
  snew(GD->LinLog, n);
  GD->time = time;
  GD->ka = 0;
  GD->kd = 0;
  GD->tDelta = time[1]-time[0];
  GD->nData = nData;
  GD->params = params;

  fitFunc.f = &gemFunc_residual2;
  fitFunc.n = 2;
  fitFunc.params = (void*)GD;

  x  = gsl_vector_alloc (fitFunc.n);
  dx = gsl_vector_alloc (fitFunc.n);
  gsl_vector_set (x,  0, 25);
  gsl_vector_set (x,  1, 0.5);
  gsl_vector_set (dx, i, 0.1);
  gsl_vector_set (dx, i, 0.01);
  
  
  s = gsl_multimin_fminimizer_alloc (T, fitFunc.n);
  gsl_multimin_fminimizer_set (s, &fitFunc, x, dx);
  gsl_vector_free (x);
  gsl_vector_free (dx);

  do  {
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);
    
    if (status != 0)
      gmx_fatal(FARGS,"Something went wrong in the iteration in minimizer %s:\n \"%s\"\n",
		gsl_multimin_fminimizer_name(s), gsl_strerror(status));
    
    d2     = gsl_multimin_fminimizer_minimum(s);
    size   = gsl_multimin_fminimizer_size(s);
    params->ka = gsl_vector_get (s->x, 0);
    params->kd = gsl_vector_get (s->x, 1);
    
    if (status)
      {
	fprintf(stderr, "%s\n", gsl_strerror(status));
	break;
      }

    status = gsl_multimin_test_size(size,tol);
    if (status == GSL_SUCCESS)
      printf("Converged to minimum at\n");
    printf ("iter %5d: ka = %2.5f  kd = %2.5f  f() = %7.3f  size = %.3f  chi2 = %2.5f\n",
	    iter,
	    params->ka,
	    params->kd,
	    s->fval, size, d2);
  }
  while ((status == GSL_CONTINUE) && (iter < maxiter));

  /*   /\* Calculate the theoretical ACF from the parameters one last time. *\/ */
  eq10v2(GD->ctTheory, time, nData, params->ka, params->kd, params);
  *ctFit = GD->ctTheory;

  sfree(GD);
  gsl_multimin_fminimizer_free (s);


  return d2;

#endif /* HAVE_LIBGSL */
}

#ifdef HAVE_LIBGSL
static int balFunc_f(const gsl_vector *x, void *data, gsl_vector *f)
{
  /* C + sum{ A_i * exp(-B_i * t) }*/

  balData *BD = (balData *)data;
  int n = BD->n;
  int i,j, nexp = BD->nexp;
  double
    *y     = BD->y,
    *A, *B, C,        /* There are the parameters to be optimized. */
    t, ct;

  snew(A, nexp);
  snew(B, nexp);
  
  for (i = 0; i<nexp; i++)
    {
      A[i] = gsl_vector_get(x, i*2);
      B[i] = gsl_vector_get(x, i*2+1);
    }
  C = gsl_vector_get(x, nexp*2);

  for (i=0; i<n; i++)
    {
      t = i*BD->tDelta;
      ct = 0;
      for (j=0; j<nexp; j++)
	ct += A[j] * exp(B[j] * t);
      ct += C;
      gsl_vector_set (f, i, ct - y[i]);
    }
  return GSL_SUCCESS;
}

/* The derivative stored in jacobian form (J)*/
static int balFunc_df(const gsl_vector *params, void *data, gsl_matrix *J)
{
  balData *BD = (balData*)data;
  size_t n = BD->n,
    i,j;
  double    *y     = BD->y,
    *A, *B, C, /* There are the parameters. */
    t;
  int nexp = BD->nexp;
  snew(A, nexp);
  snew(B, nexp);

  for (i=0; i<nexp; i++)
    {
      A[i] = gsl_vector_get(params, i*2);
      B[i] = gsl_vector_get(params, i*2+1);
    }
  C = gsl_vector_get(params, nexp*2);
  for (i=0; i<n; i++)
    {
      t = i*BD->tDelta;
      for (j=0; j<nexp; j++)
	{
	  gsl_matrix_set (J, i, j*2,   exp(B[j]*t));        /* df(t)/dA_j */
	  gsl_matrix_set (J, i, j*2+1, A[j]*t*exp(B[j]*t)); /* df(t)/dB_j */
	}
      gsl_matrix_set (J, i, nexp*2, 1); /* df(t)/dC */
    }
  return GSL_SUCCESS;
}

/* Calculation of the function and its derivative */
static int balFunc_fdf(const gsl_vector *params, void *data,
		       gsl_vector *f, gsl_matrix *J)
{
  balFunc_f(params, data, f);
  balFunc_df(params, data, J);
  return GSL_SUCCESS;
}
#endif /* HAVE_LIBGSL */

/* Removes the ballistic term from the beginning of the ACF,
 * just like in Omer's paper.
 */
extern void takeAwayBallistic(double *ct, double *t, int len, real tMax, int nexp, bool bDerivative)
{

  /* Use nonlinear regression with GSL instead.
   * Fit with 4 exponentials and one constant term,
   * subtract the fatest exponential. */

  /* /============  Set up the solver =============\ */

#ifndef HAVE_LIBGSL
  printf("Sorry, can't take away ballistic component without gsl. "
	 "Recompile using --with-gsl.\n");

#else
  int nData=0,i,status, iter=0;
  do {
    nData++;
  } while (t[nData]<tMax+t[0] && nData<len);

  /* Set the solver type. */
  const gsl_multifit_fdfsolver_type *T
    = gsl_multifit_fdfsolver_lmsder;

  const size_t p = nexp*2+1;              /* Number of parameters. */
  gsl_multifit_fdfsolver *s;              /* The solver itself. */
  gsl_multifit_function_fdf fitFunction;  /* The function to be fitted. */
  gsl_matrix *covar = gsl_matrix_alloc (p, p);  /* Covariance matrix for the parameters.
						 * We'll not use the result, though. */
  double *guess=NULL,                /* Initial guess. */
    *A,                              /* The fitted parameters. (A1, B1, A2, B2,... C) */
    a[2];
  const size_t n = nData;
  bool sorted;

  snew(guess, p);
  snew(A, p);

  /* Set up an initial gess for the parameters.
   * The solver is somewhat sensitive to the initial guess,
   * but this worked fine for a TIP3P box with -geminate dd
   * EDIT: In fact, this seems like a good starting pont for other watermodels too. */
  for (i=0; i<nexp; i++)
    {
      guess[i*2] = 0.1;
      guess[i*2+1] = -0.5 + (((double)i)/nexp - 0.5)*0.3;
    }
  guess[nexp * 2] = 0.01;

  gsl_vector_view theParams = gsl_vector_view_array(guess, p);
  balData *BD;
  snew(BD,1);
  BD->n     = n;
  BD->y     = ct;
  BD->tDelta = t[1]-t[0];
  BD->nexp = nexp;

  fitFunction.f      =  &balFunc_f;
  fitFunction.df     =  &balFunc_df;
  fitFunction.fdf    =  &balFunc_fdf;
  fitFunction.n      =  nData;
  fitFunction.p      =  p;
  fitFunction.params =  BD;

  s = gsl_multifit_fdfsolver_alloc (T, nData, p);
  if (s==NULL)
    gmx_fatal(FARGS, "Could not set up the nonlinear solver.");

  gsl_multifit_fdfsolver_set(s, &fitFunction, &theParams.vector);

  /* \=============================================/ */


  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      
      if (status)
	break;
      status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
  while (iter < 5000 && status == GSL_CONTINUE);
  if (iter == 5000)
    printf("The non-linear fitting did not converge in 5000 steps.\n"
	   "Check the quality of the fit!\n");
  else
    printf("Non-linear fitting of ballistic term converged in %d steps.\n\n", (int)iter);
  for (i=0; i<nexp; i++)
    printf("%c * exp(%c * t) + ", 'A'+(char)i*2, 'B'+(char)i*2);
  printf("%c\n", 'A'+(char)nexp*2);
  printf("Here are the actual numbers for A-%c:\n", 'A'+nexp*2);
  for (i=0; i<nexp; i++)
    {
      A[i*2] = gsl_vector_get(s->x, i*2);
      A[i*2+1] = gsl_vector_get(s->x, i*2+1);
      printf(" %g*exp(%g * x) +", A[i*2], A[i*2+1]);
    }
  A[i*2] = gsl_vector_get(s->x, i*2);          /* The last and constant term */
  printf(" %g\n", A[i*2]);

  /* Implement some check for parameter quality */
  for (i=0; i<nexp; i++)
    {
      if (A[i*2]<0 || A[i*2]>1)
	printf("WARNING: ----------------------------------\n"
	       " | A coefficient does not lie within [0,1].\n"
	       " | This may or may not be a problem.\n"
	       " | Double check the quality of the fit!\n");
      if (A[i*2+1]>0)
	printf("WARNING: ----------------------------------\n"
	       " | One factor in the exponent is positive.\n"
	       " | This could be a problem if the coefficient\n"
	       " | is large. Double check the quality of the fit!\n");
    }
  if (A[i*2]<0 || A[i*2]>1)
    printf("WARNING: ----------------------------------\n"
	   " | The constant term does not lie within [0,1].\n"
	   " | This may or may not be a problem.\n"
	   " | Double check the quality of the fit!\n");

  /* Sort the terms */
  sorted = (nexp > 1) ?  FALSE : TRUE;
  while (!sorted)
    {
      sorted = TRUE;
      for (i=0;i<nexp-1;i++)
	{
	  double ddt[2] = {
	    A[i*2] * A[i*2+1],
	    A[i*2+2] * A[i*2+3]
	  };
	  
	  if ((bDerivative && (ddt[0]<0 && ddt[1]<0 && ddt[0]>ddt[1])) || /* Compare derivative at t=0... */
	      (!bDerivative && (A[i*2+1] > A[i*2+3])))                    /* Or just the coefficient in the exponent */
	    {
	      sorted = FALSE;
	      a[0] = A[i*2];  /* coefficient */
	      a[1] = A[i*2+1]; /* parameter in the exponent */
	      
	      A[i*2] = A[i*2+2];
	      A[i*2+1] = A[i*2+3];
	      
	      A[i*2+2] = a[0];
	      A[i*2+3] = a[1];
	    }
	}
    }

  /* Subtract the fastest component */
  printf("Fastest component is %g * exp(%g * t)\n"
	 "Subtracting fastest component from ACF.\n", A[0], A[1]);

  for (i=0; i<len; i++)
    ct[i] = (ct[i] - A[0] * exp(A[1] * i*BD->tDelta)) / (1-A[0]);
  
  sfree(guess);
  sfree(A);

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
  fflush(stdout);
#endif /* HAVE_LIBGSL */
}
