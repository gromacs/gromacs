#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "random.h"
#include "xvgr.h"
#include "vec.h"
#include "ionize.h"

typedef struct {
  real photo,coh,incoh,incoh_abs;
} t_cross;

typedef struct {
  char *name;
  int  nel;
  t_cross *cross;
} t_element;

void ionize(FILE *log,t_mdatoms *md,char **atomname[],real t,t_inputrec *ir)
{
  static FILE *xvg;
  static char *leg[] = { "Probability", "Number/Atom", "Total" };
  static bool bFirst = TRUE,bDoIt;
  static real pt,t0,imax,width,inv_nratoms,rho,nphot;
  static int  seed,total,*nelectrons;
  static real *cross_sec_atom;
  static int  Eindex=-1;
  static int  Energies[] = { 6, 8, 10, 15, 20 };
#define NENER asize(Energies)
  static t_cross   cross_sec_h[] = {
    { 2.63e-2,     1.01e-1,      5.49e-1,         7.12e-3 },
    { 9.79e-3,     6.18e-2,      5.83e-1,         9.60e-3 },
    { 4.55e-3,     4.16e-2,      5.99e-1,         1.19e-2 },
    { 1.12e-3,     1.96e-2,      6.09e-1,         1.73e-2 },
    { 4.16e-4,     1.13e-2,      6.07e-1,         2.23e-2 }
  };
  static t_cross   cross_sec_c[] = {
    { 2.04e+2,     5.88e+0,      2.29e+0,         3.06e-2 },
    { 8.22e+1,     4.22e+0,      2.56e+0,         4.38e-2 },
    { 4.03e+1,     3.26e+0,      2.74e+0,         5.72e-2 },
    { 1.09e+1,     1.97e+0,      3.04e+0,         9.15e-2 },
    { 4.29e+0,     1.30e+0,      3.20e+0,         1.24e-1 }
  };
  static t_cross   cross_sec_n[] = {
    { 4.04e+2,     8.99e+0,      2.49e+0,         3.43e-2 },
    { 1.65e+2,     6.29e+0,      2.86e+0,         5.01e-2 },
    { 8.15e+1,     4.76e+0,      3.10e+0,         6.57e-2 },
    { 2.24e+1,     2.82e+0,      3.46e+0,         1.05e-1 },
    { 8.87e+0,     1.88e+0,      3.65e+0,         1.43e-1 }
  };
  static t_cross   cross_sec_o[] = {
    { 7.18e+2,     1.33e+1,      2.66e+0,         3.75e-2 },
    { 2.96e+2,     9.21e+0,      3.14e+0,         5.62e-2 },
    { 1.47e+2,     6.85e+0,      3.44e+0,         7.43e-2 },
    { 4.09e+1,     3.97e+0,      3.87e+0,         1.18e-1 },
    { 1.63e+1,     2.64e+0,      4.10e+0,         1.61e-1 }
  };
  static t_cross   cross_sec_s[] = {
    { 8.78e+3,      5.54e+1,     3.98e+0,         5.42e-2 },
    { 3.83e+3,      4.29e+1,     4.71e+0,         8.38e-2 },
    { 1.99e+3,      3.36e+1,     5.32e+0,         1.16e-1 },
    { 5.96e+2,      1.97e+1,     6.36e+0,         2.00e-1 },
    { 2.49e+2,      1.29e+1,     6.94e+0,         2.80e-1 }
  };
  
  static t_element element[] = {
    { "H",   1, cross_sec_h },
    { "C",   6, cross_sec_c },
    { "N",   7, cross_sec_n },
    { "O",   8, cross_sec_o },
    { "S",  16, cross_sec_s }
  };
#define NELEM asize(element)
  real r;
  char *cc;
  int  i,j,dq,*elemcnt;
  
  if (bFirst) {
    /* Get parameters for gaussian photon pulse from inputrec */
    t0    = ir->userreal1;  /* Peak of the gaussian pulse            */
    nphot = ir->userreal2;  /* Intensity                             */
    width = ir->userreal3;  /* Width of the peak                     */
    rho   = ir->userreal4;  /* Diameter of the focal spot            */
    seed  = ir->userint1;   /* Random seed for stochastic ionization */
    imax  = (nphot/(M_PI*rho*rho))*1e-10*sqrt(1.0/M_PI)*(2.0/width);

    bDoIt = ((width > 0) && (imax > 0));
    
    if (bDoIt) {
      if (seed == 0)
	seed = 1993;
      
      for(Eindex=0; (Eindex < NENER) && (Energies[Eindex] != ir->userint2); 
	  Eindex++)
	;
      if (Eindex == NENER)
	fatal_error(0,"Energy level of %d keV not supported",ir->userint2);
      
      total = 0;
      inv_nratoms = 1.0/md->nr;
      snew(nelectrons,md->nr);
      snew(cross_sec_atom,md->nr);
      snew(elemcnt,NELEM);
      for(i=0; (i<md->nr); i++) {
	cc = *(atomname[i]);
	for(j=0; (j<NELEM); j++)
	  if (strstr(cc,element[j].name) != NULL) {
	    nelectrons[i]     = element[j].nel;
	    cross_sec_atom[i] = element[j].cross[Eindex].photo + 
	      element[j].cross[Eindex].incoh;
	    break;
	  }
	if (j == NELEM) {
	  j = 1; /* Carbon, hopefully */
	  nelectrons[i]     = element[j].nel;
	  cross_sec_atom[i] = element[j].cross[Eindex].photo + 
	    element[j].cross[Eindex].incoh;
	  fprintf(stderr,"Don't know number of electrons for %s, assuming %d\n",
		  cc,nelectrons[i]);
	}
	elemcnt[j]++;
      }
      xvg   = xvgropen("ionize.xvg","Ionization Events","Time (ps)","()");
      xvgr_legend(xvg,asize(leg),leg);
      
      fprintf(log,"Parameters for ionization events:\n");
      fprintf(log,"Imax = %g, t0 = %g, width = %g, seed = %d\n"
	      "# Photons = %g, rho = %g\n",
	      imax,t0,width,seed,nphot,rho);
      fprintf(log,"You have the following elements in your system:\n");
      for(j=0; (j<NELEM); j++)
	fprintf(log,"  %s: %d",element[j].name,elemcnt[j]);
      fprintf(log," atoms\n");
    }
    bFirst = FALSE;
  }
  if (!bDoIt)
    return;
    
  /* Calculate probability */
  pt = imax*ir->delta_t*exp(-sqr(2.0*(t-t0)/width));
  dq = 0;
  for(i=0; (i<md->nr); i++) {
    r = rando(&seed);
    if (r < pt*cross_sec_atom[i]) {
      if (nelectrons[i] >= 1) {
	md->chargeA[i] -= 1.0;
	md->chargeB[i] -= 1.0;
	dq ++;
	total ++;
      }
    }
  }
  fprintf(xvg,"%10g  %10g  %10g  %10g\n",
	  t,pt,dq*inv_nratoms,total*inv_nratoms);
  fflush(xvg);
}
