#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"
#include "fatal.h"
#include "random.h"
#include "futil.h"
#include "physics.h"
#include "ehdata.h"

typedef struct {
  int  nener,nomega,nq;     /* Number of entries in the table      */
  real *ener;               /* Energy values                       */
  real **omega;             /* Omega values (energy dependent)     */
  real ***prob,***q;        /* Probability and Q                   */
} t_pq_inel;

typedef struct {
  int  nener,n2Ddata;       /* Number of entries in the table      */
  real *ener;               /* Energy values                       */
  real **prob,**data;       /* Probability and data (energy loss)  */
} t_p2Ddata;

static int realcomp(const void *a,const void *b)
{
  real *ra,*rb;
  
  ra = (real *)a;
  rb = (real *)b;
  if (ra < rb)
    return -1;
  else if (ra > rb)
    return 1;
  else 
    return 0;
}

static t_p2Ddata *read_p2Ddata(char *fn,int nener,int n2Ddata)
{
  FILE    *fp;
  t_p2Ddata *p2Ddata;
  int     i,j;
  double  e,p,o;
  
  fprintf(stdout,"Going to read %s\n",fn);
  fp = ffopen(fn,"r");

  /* Allocate memory and set constants */
  snew(p2Ddata,1);
  p2Ddata->nener  = nener;
  p2Ddata->n2Ddata = n2Ddata;
  snew(p2Ddata->ener,p2Ddata->nener);
  snew(p2Ddata->prob,p2Ddata->nener);
  snew(p2Ddata->data,p2Ddata->nener);
  
  /* Double loop to read data */
  for(i=0; (i<p2Ddata->nener); i++) {
    fprintf(stderr,"\rEnergy %d/%d",i+1,p2Ddata->nener);
    snew(p2Ddata->prob[i],p2Ddata->n2Ddata);
    snew(p2Ddata->data[i],p2Ddata->n2Ddata);
    
    for(j=0; (j<p2Ddata->n2Ddata); j++) {
      fscanf(fp,"%lf%lf%lf",&e,&p,&o);

      /* Consistency check */
      if (j==0)
	p2Ddata->ener[i] = e;
      else if (fabs(p2Ddata->ener[i]-e) > 1e-6*e)
	fatal_error(0,"Inconsistent energy %f i=%d j=%d",e,i,j);
      p2Ddata->prob[i][j] = p;
      p2Ddata->data[i][j] = o;
    }
    /* There is some noise on the data, take it away by sorting,
     * because otherwise binary search does not work.
     * This is equivalent to shifting in the data slightly along the X-axis
     * but better than linear search with the "real" data.
     */
    qsort(p2Ddata->data[i],p2Ddata->n2Ddata,sizeof(p2Ddata->data[0][0]),realcomp);
  }
  fprintf(stderr,"\n");
  
  fclose(fp);
  
  return p2Ddata;
}

static t_pq_inel *read_pq(char *fn)
{
  FILE      *fp;
  t_pq_inel *pq;
  int       i,j,k;
  double    e,p,o,t;
  
  fprintf(stdout,"Going to read %s\n",fn);
  fp = ffopen(fn,"r");
  
  /* Allocate memory and set constants */
  snew(pq,1);
  pq->nener  = 420;
  pq->nomega = 49;
  pq->nq = 101;
  snew(pq->ener,pq->nener);
  snew(pq->omega,pq->nener);
  snew(pq->prob,pq->nener);
  snew(pq->q,pq->nener);
  
  /* Triple loop to read data */
  for(i=0; (i<pq->nener); i++) {
    fprintf(stderr,"\rEnergy %d/%d",i+1,pq->nener);
    snew(pq->prob[i],pq->nomega);
    snew(pq->q[i],pq->nomega);
    snew(pq->omega[i],pq->nomega);
    
    for(j=0; (j<pq->nomega); j++) {
      snew(pq->prob[i][j],pq->nq);
      snew(pq->q[i][j],pq->nq);
      
      for(k=0; (k<pq->nq); k++) {
	fscanf(fp,"%lf%lf%lf%lf",&e,&o,&p,&t);
	
	/* Consistency check */
	if ((j == 0) && (k == 0)) 
	  pq->ener[i] = e;
	else if (fabs(pq->ener[i]-e) > 1e-6*e)
	  fatal_error(0,"Inconsistent energy %f i=%d j=%d k=%d",e,i,j,k);
	
	if (k == 0)
	  pq->omega[i][j] = o;
	else if (fabs(pq->omega[i][j]-o) > 1e-6*o)
	  fatal_error(0,"Inconsistent omega %f i=%d j=%d k=%d",o,i,j,k);
	
	pq->prob[i][j][k] = p;
	pq->q[i][j][k] = t;
      }
    }
  }
  fprintf(stderr,"\n");
  
  fclose(fp);
  
  return pq;
}

static int my_bsearch(real val,int ndata,real data[])
{
  int ilo,ihi,imed;

  ilo = 0; 
  ihi = ndata-1;
  while ((ihi - ilo) > 1) {
    imed = (ilo+ihi)/2;
    if (data[imed] > val) 
      ihi = imed;
    else
      ilo = imed;
  }
  /* Now val should be in between data[ilo] and data[ihi] */
  /* Decide which one is closest */
  if ((val-data[ilo]) > (data[ihi]-val))
    return ihi;
  else
    return ilo;
}

real get_omega(real ekin,int *seed,FILE *fp)
{
  static t_p2Ddata *p2Ddata = NULL;
  real r,ome;
  int  eindex,oindex;
  
  if (p2Ddata == NULL) 
    p2Ddata = read_p2Ddata("respd.dat",1019,100);
  
  /* Get energy index by binary search */
  eindex = my_bsearch(ekin,p2Ddata->nener,p2Ddata->ener);
  
  /* Start with random number */
  r = rando(seed);
  
  /* Do binary search in the energy table */
  oindex = my_bsearch(r,p2Ddata->n2Ddata,p2Ddata->prob[eindex]);
  
  ome = p2Ddata->data[eindex][oindex];
  
  if (fp) 
    fprintf(fp,"%8.3f  %8.3f\n",ome,r);
  
  return ome;
}

real get_theta_el(real ekin,int *seed,FILE *fp)
{
  static t_p2Ddata *p2Ddata = NULL;
  real r,theta;
  int  eindex,tindex;
  
  if (p2Ddata == NULL) 
    p2Ddata = read_p2Ddata("proeldds.dat",1200,101);
  
  /* Start with random number */
  r = rando(seed);
    
  /* Get energy index by binary search */
  eindex = my_bsearch(ekin,p2Ddata->nener,p2Ddata->ener);
  
  /* Do binary search in the energy table */
  tindex = my_bsearch(r,p2Ddata->n2Ddata,p2Ddata->prob[eindex]);
  
  theta = p2Ddata->data[eindex][tindex];
  
  if (fp) 
    fprintf(fp,"%8.3f  %8.3f\n",theta,r);
  
  return theta;
}

real get_q_inel(real ekin,real omega,int *seed,FILE *fp)
{
  static t_pq_inel *pq = NULL;
  int    eindex,oindex,tindex;
  real   r,the;
  
  if (pq == NULL)
    pq = read_pq("spresp-all.dat");

  /* Get energy index by binary search */
  eindex = my_bsearch(ekin,pq->nener,pq->ener);
  
  /* Do binary search in the energy table */
  oindex = my_bsearch(omega,pq->nomega,pq->omega[eindex]);
  
  /* Start with random number */
  r = rando(seed);
  
  tindex = my_bsearch(r,pq->nq,pq->prob[eindex][oindex]);
  
  the = pq->q[eindex][oindex][tindex];
  
  if (fp)
    fprintf(fp,"%8.3f  %8.3f\n",the,r);
    
  return the;
}


static void read_cross(char *fn,int n,real **ener,real **cross,real factor)
{
  FILE   *fp;
  double e,c;
  int    i;
  
  fprintf(stdout,"Going to read %s\n",fn);
  fp = ffopen(fn,"r");

  /* Allocate memory */
  snew(*cross,n);
  snew(*ener,n);
  for(i=0; (i<n); i++) {
    fscanf(fp,"%lf%lf",&e,&c);
    (*ener)[i] = e;
    (*cross)[i] = c*factor;
  }
  fclose(fp);
  
}

real cross_inel(real ekin)
{
  static real *ener  = NULL;
  static real *cross = NULL;
  static int  ninel = 520;
  int eindex;
  
  /* Read data at first call, convert A^2 to nm^2 */
  if (cross == NULL) 
    read_cross("totpineld.dat",ninel,&ener,&cross,0.01);
  
  /* Compute index with binary search */
  eindex = my_bsearch(ekin,ninel,ener);
  
  return cross[eindex];
}

real cross_el(real ekin)
{
  static real *ener  = NULL;
  static real *cross = NULL;
  static int  nel = 1200;
  int eindex;
  
  /* Read data at first call, convert A^2 to nm^2  */
  if (cross == NULL) {
    real rho    = 3.51; /* g/cm^3 */
    real factor = (rho/12.011)*AVOGADRO*1e-14*NANO;
    read_cross("totpeldds.dat",nel,&ener,&cross,factor);
  }
  /* Compute index with binary search */
  eindex = my_bsearch(ekin,nel,ener);
  
  return cross[eindex];
}

real band_ener(int *seed,FILE *fp)
{
  static real *ener  = NULL;
  static real *prob  = NULL;
  static int  nener = 500;
  int  eindex;
  real r;
  
  /* Read data at first call, misuse read_cross function */
  if (prob == NULL) 
    read_cross("realbd.dat",nener,&ener,&prob,1.0);
  
  r = rando(seed);
  
  eindex = my_bsearch(r,nener,prob);
  
  if (fp)
    fprintf(fp,"%8.3f  %8.3f\n",ener[eindex],r);
  
  return ener[eindex];
}

static void test_omega(FILE *fp,int *seed)
{
  int i;

  fprintf(fp,"Testing the energy loss tables\n");
  for(i=0; (i<1000); i++) {
    (void) get_omega(500*rando(seed),seed,fp);
  }
}

static void test_q_inel(FILE *fp,int *seed)
{
  int i;
  
  fprintf(fp,"Testing the energy/omega dependent inelastic scattering q tables\n");
  for(i=0; (i<1000); i++) {
    (void) get_q_inel(500*rando(seed),400*rando(seed),seed,fp);
  }
}

static void test_theta_el(FILE *fp,int *seed)
{
  int i;
  
  fprintf(fp,"Testing the energy dependent elastic scattering theta tables\n");
  for(i=0; (i<1000); i++) {
    (void) get_theta_el(500*rando(seed),seed,fp);
  }
}

static void test_sigma_el(FILE *fp)
{
  int  i;

  fprintf(fp,"Testing the elastic cross sections table\n");
  for(i=0; (i<500); i++) {
    fprintf(fp,"%3d  %8.3f\n",i,cross_el(i));
  }
}

static void test_sigma_inel(FILE *fp)
{
  int  i;

  fprintf(fp,"Testing the inelastic cross sections table\n");
  for(i=0; (i<500); i++) {
    fprintf(fp,"%3d  %8.3f\n",i,cross_inel(i));
  }
}

static void test_band_ener(int *seed,FILE *fp)
{
  int i;
  
  for(i=0; (i<500); i++) {
    band_ener(seed,fp);
  }
}

void test_tables(int *seed,char *fn)
{
  FILE *fp;
  
  fp = fopen(fn,"w");

  test_omega(fp,seed);
  test_q_inel(fp,seed);
  test_theta_el(fp,seed);
  test_sigma_el(fp);
  test_sigma_inel(fp);
  test_band_ener(seed,fp);
  
  fclose(fp);
}

