#include <math.h>
#include <string.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "random.h"
#include "physics.h"
#include "xvgr.h"
#include "vec.h"
#include "txtdump.h"
#include "ionize.h"
#include "names.h"
#include "futil.h"

typedef struct {
  real photo,coh,incoh,incoh_abs;
} t_cross;

static t_cross   cross_sec_h[] = {
  { 2.63e-2,     1.01e-1,      5.49e-1,         7.12e-3 },
  { 9.79e-3,     6.18e-2,      5.83e-1,         9.60e-3 },
  { 4.55e-3,     4.16e-2,      5.99e-1,         1.19e-2 },
  { 1.12e-3,     1.96e-2,      6.09e-1,         1.73e-2 },
  { 4.16e-4,     1.13e-2,      6.07e-1,         2.23e-2 }
};
static t_cross   cross_sec_he[] = {
  { 1.08e+0,     8.55e-1,      8.75e-1,         1.23e-2 },
  { 4.10e-1,     5.63e-1,      1.01e+0,         1.79e-2 },
  { 1.92e-1,     3.95e-1,      1.08e+0,         2.30e-2 },
  { 4.85e-2,     1.98e-1,      1.16e+0,         3.43e-2 },
  { 1.82e-2,     1.18e-1,      1.18e+0,         4.43e-2 }
};
static t_cross   cross_sec_li[] = {
  { 8.23e+0,     1.60e+0,      1.20e+0,         1.59e-2 },
  { 3.17e+0,     1.16e+0,      1.38e+0,         2.41e-2 },
  { 1.51e+0,     8.65e-1,      1.51e+0,         3.21e-2 },
  { 3.88e-1,     4.65e-1,      1.66e+0,         5.01e-2 },
  { 1.48e-1,     2.87e-1,      1.72e+0,         6.58e-2 }
  };
static t_cross   cross_sec_be[] = {
  { 3.34e+1,     2.47e+0,      1.60e+0,         2.04e-2 },
  { 1.31e+1,     1.88e+0,      1.78e+0,         3.01e-2 },
  { 6.28e+0,     1.46e+0,      1.92e+0,         4.02e-2 },
  { 1.64e+0,     8.39e-1,      2.14e+0,         6.46e-2 },
  { 6.29e-1,     5.32e-1,      2.23e+0,         8.63e-2 }
};
static t_cross   cross_sec_bo[] = {
  { 9.08e+1,     3.78e+0,      1.98e+0,         2.57e-2 },
  { 3.61e+1,     2.82e+0,      2.18e+0,         3.68e-2 },
  { 1.76e+1,     2.21e+0,      2.34e+0,         4.85e-2 },
  { 4.69e+0,     1.33e+0,      2.60e+0,         7.81e-2 },
  { 1.82e+0,     8.62e-1,      2.72e+0,         1.06e-1 }
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
static t_cross   cross_sec_f[] = {
  { 1.16e+3,     1.89e+1,      2.76e+0,         3.96e-2 },
  { 4.82e+2,     1.31e+1,      3.36e+0,         6.13e-2 },
  { 2.42e+2,     9.69e+0,      3.74e+0,         8.21e-2 },
  { 6.83e+1,     5.50e+0,      4.26e+0,         1.31e-1 },
  { 2.75e+1,     3.63e+0,      4.53e+0,         1.79e-1 }
};
static t_cross   cross_sec_ne[] = {
  { 1.78e+3,     2.60e+1,      2.78e+0,         4.04e-2 },
  { 7.51e+2,     1.82e+1,      3.48e+0,         6.47e-2 },
  { 3.80e+2,     1.34e+1,      3.96e+0,         8.85e-2 },
  { 1.08e+2,     7.50e+0,      4.61e+0,         1.44e-1 },
  { 4.39e+1,     4.90e+0,      4.94e+0,         1.96e-1 }
};
static t_cross   cross_sec_na[] = {
  { 2.61e+3,     3.18e+1,      2.95e+0,         4.19e-2 },
  { 1.11e+3,     2.27e+1,      3.69e+0,         6.79e-2 },
  { 5.63e+2,     1.69e+1,      4.22e+0,         9.43e-2 },
  { 1.63e+2,     9.47e+0,      4.97e+0,         1.56e-1 },
  { 6.65e+1,     6.18e+0,      5.35e+0,         2.13e-1 }
};
static t_cross   cross_sec_mg[] = {
  { 3.69e+3,     3.75e+1,      3.19e+0,         4.42e-2 },
  { 1.57e+3,     2.75e+1,      3.92e+0,         7.13e-2 },
  { 8.06e+2,     2.06e+1,      4.49e+0,         9.98e-2 },
  { 2.34e+2,     1.17e+1,      5.33e+0,         1.68e-1 },
  { 9.60e+1,     7.61e+0,      5.76e+0,         2.31e-1 }
};
static t_cross   cross_sec_al[] = {
  { 5.15e+3,     4.32e+1,      3.46e+0,         4.72e-2 },
  { 2.22e+3,     3.24e+1,      4.18e+0,         7.50e-2 },
  { 1.14e+3,     2.47e+1,      4.76e+0,         1.05e-1 },
  { 3.35e+2,     1.41e+1,      5.68e+0,         1.79e-1 },
  { 1.38e+2,     9.19e+0,      6.16e+0,         2.48e-1 }
};
static t_cross   cross_sec_si[] = {
  { 6.70e+3,     4.90e+1,      3.73e+0,         5.06e-2 },
  { 2.90e+3,     3.75e+1,      4.45e+0,         7.92e-2 },
  { 1.49e+3,     2.90e+1,      5.04e+0,         1.11e-1 },
  { 4.40e+2,     1.68e+1,      6.03e+0,         1.90e-1 },
  { 1.83e+2,     1.09e+1,      6.56e+0,         2.64e-1 }
};
static t_cross   cross_sec_p[] = {
  { 8.78e+3,     5.54e+1,      3.98e+0,         5.42e-2 },
  { 3.83e+3,     4.29e+1,      4.71e+0,         8.38e-2 },
  { 1.99e+3,     3.36e+1,      5.32e+0,         1.16e-1 },
  { 5.96e+2,     1.97e+1,      6.36e+0,         2.00e-1 },
  { 2.49e+2,     1.29e+1,      6.94e+0,         2.80e-1 }
};
static t_cross   cross_sec_s[] = {
  { 8.78e+3,      5.54e+1,     3.98e+0,         5.42e-2 },
  { 3.83e+3,      4.29e+1,     4.71e+0,         8.38e-2 },
  { 1.99e+3,      3.36e+1,     5.32e+0,         1.16e-1 },
  { 5.96e+2,      1.97e+1,     6.36e+0,         2.00e-1 },
  { 2.49e+2,      1.29e+1,     6.94e+0,         2.80e-1 }
};
static t_cross   cross_sec_cl[] = {
  { 1.40e+4,      7.06e+1,     4.47e+0,         6.17e-2 },
  { 6.20e+3,      5.47e+1,     5.27e+0,         9.40e-2 },
  { 3.26e+3,      4.37e+1,     5.91e+0,         1.29e-1 },
  { 9.88e+2,      2.65e+1,     7.04e+0,         2.21e-1 },
  { 4.17e+2,      1.75e+1,     7.71e+0,         3.12e-1 }
};
static t_cross   cross_sec_ar[] = {
  { 1.77e+4,      8.00e+1,     4.63e+0,         6.46e-2 },
  { 7.83e+3,      6.15e+1,     5.51e+0,         9.90e-2 },
  { 4.12e+3,      4.92e+1,     6.18e+0,         1.35e-1 },
  { 1.25e+3,      3.03e+1,     7.36e+0,         2.31e-1 },
  { 5.30e+2,      2.01e+1,     8.08e+0,         3.27e-1 }
};
static t_cross   cross_sec_k[] = {
  { 2.11e+4,      8.86e+1,     4.88e+0,         6.77e-2 },
  { 9.44e+3,      6.78e+1,     5.81e+0,         1.04e-1 },
  { 5.00e+3,      5.44e+1,     6.50e+0,         1.42e-1 },
  { 1.54e+3,      3.40e+1,     7.71e+0,         2.42e-1 },
  { 6.58e+2,      2.27e+1,     8.46e+0,         3.42e-1 }
};
static t_cross   cross_sec_ca[] = {
  { 2.54e+4,      9.77e+1,     5.14e+0,         7.08e-2 },
  { 1.15e+4,      7.44e+1,     6.12e+0,         1.09e-1 },
  { 6.10e+3,      5.96e+1,     6.84e+0,         1.49e-1 },
  { 1.89e+3,      3.77e+1,     8.06e+0,         2.52e-1 },
  { 8.09e+2,      2.54e+1,     8.84e+0,         3.57e-1 }
};
static t_cross   cross_sec_sc[] = {
  { 3.03e+4,      1.09e+2,     5.28e+0,         7.29e-2 },
  { 1.37e+4,      8.28e+1,     6.33e+0,         1.14e-1 },
  { 7.35e+3,      6.60e+1,     7.10e+0,         1.55e-1 },
  { 2.30e+3,      4.19e+1,     8.39e+0,         2.63e-1 },
  { 9.88e+2,      2.85e+1,     9.20e+0,         3.72e-1 }
};
static t_cross   cross_sec_ti[] = {
  { 3.53e+4,      1.23e+2,     5.39e+0,         7.47e-2 },
  { 1.63e+4,      9.23e+1,     6.50e+0,         1.17e-1 },
  { 8.82e+3,      7.32e+1,     7.33e+0,         1.61e-1 },
  { 2.79e+3,      4.65e+1,     8.70e+0,         2.74e-1 },
  { 1.21e+3,      3.18e+1,     9.56e+0,         3.86e-1 }
};
static t_cross   cross_sec_v[] = {
  { 4.16e+4,      1.38e+2,     5.45e+0,         7.58e-2 },
  { 1.91e+4,      1.03e+2,     6.63e+0,         1.20e-1 },
  { 1.03e+4,      8.13e+1,     7.53e+0,         1.66e-1 },
  { 3.27e+3,      5.14e+1,     8.99e+0,         2.84e-1 },
  { 1.42e+3,      3.54e+1,     9.89e+0,         4.01e-1 }
};
static t_cross   cross_sec_cr[] = {
  { 4.90e+4,      1.56e+2,     5.37e+0,         7.55e-2 },
  { 2.26e+4,      1.16e+2,     6.64e+0,         1.22e-1 },
  { 1.23e+4,      9.11e+1,     7.62e+0,         1.70e-1 },
  { 3.90e+3,      5.72e+1,     9.22e+0,         2.93e-1 },
  { 1.70e+3,      3.94e+1,     1.02e+1,         4.15e-1 }
};
static t_cross   cross_sec_mn[] = {
  { 6.75e+3,      1.72e+2,     5.50e+0,         7.70e-2 },
  { 2.58e+4,      1.28e+2,     6.79e+0,         1.24e-1 },
  { 1.40e+4,      1.00e+2,     7.81e+0,         1.74e-1 },
  { 4.49e+3,      6.27e+1,     9.50e+0,         3.03e-1 },
  { 1.96e+3,      4.33e+1,     1.05e+1,         4.30e-1 }
};
static t_cross   cross_sec_fe[] = {
  { 7.89e+3,      1.92e+2,     5.55e+0,         7.80e-2 },
  { 2.95e+4,      1.43e+2,     6.89e+0,         1.26e-1 },
  { 1.62e+4,      1.11e+2,     7.96e+0,         1.78e-1 },
  { 5.24e+3,      6.91e+1,     9.75e+0,         3.12e-1 },
  { 2.31e+3,      4.78e+1,     1.08e+1,         4.44e-1 }
};
typedef struct {
  char *name;
  int  nel;
  t_cross *cross;
} t_element;

static t_element element[] = {
  { NULL,  0, cross_sec_h }, /* Empty atom, do not use! */
  { "H",   1, cross_sec_h },
  { "He",  2, cross_sec_he },
  { "Li",  3, cross_sec_li },
  { "Be",  4, cross_sec_be },
  { "B",   5, cross_sec_bo },
  { "C",   6, cross_sec_c },
  { "N",   7, cross_sec_n },
  { "O",   8, cross_sec_o },
  { "F",   9, cross_sec_f },
  { "Ne", 10, cross_sec_ne },
  { "Na", 11, cross_sec_na },
  { "Mg", 12, cross_sec_mg },
  { "Al", 13, cross_sec_al },
  { "Si", 14, cross_sec_si },
  { "P",  15, cross_sec_p  },
  { "S",  16, cross_sec_s  },
  { "Cl", 17, cross_sec_cl },
  { "Ar", 18, cross_sec_ar },
  { "K",  19, cross_sec_k  },
  { "Ca", 20, cross_sec_ca },
  { "Sc", 21, cross_sec_sc },  
  { "Ti", 22, cross_sec_ti },
  { "V",  23, cross_sec_v  },
  { "Cr", 24, cross_sec_cr },
  { "Mn", 25, cross_sec_mn },
  { "Fe", 26, cross_sec_fe },
};
#define NELEM asize(element)


/* This structure holds the energy of a K electron, L electron and the
 * probability that an L electron is shot off in the photoelectric effect
 */
typedef struct {
  real E_K,E_L,Prob_L;
} t_recoil;

t_recoil recoil[NELEM] = {
  { 0.0136,0.0,   0.0  },
  { 0.025, 0.0,   0.0  },
  { 0.054, 0.0,   0.0  },
  { 0.117, 0.0,   0.0  },
  { 0.189, 0.0,   0.0  },
  { 0.285, 0.0,   0.0  },
  { 0.402, 0.0,   0.0  },
  { 0.533, 0.0,   0.0  },
  { 0.687, 0.0,   0.0  },
  { 0.868, 0.031, 0.08 },
  { 1.072, 0.041, 0.07 },
  { 1.305, 0.054, 0.07 },
  { 1.560, 0.077, 0.08 },
  { 1.839, 0.105, 0.08 },
  { 2.146, 0.133, 0.09 },
  { 2.472, 0.166, 0.09 },
  { 2.822, 0.212, 0.10 },
  { 3.203, 0.247, 0.10 },
  { 3.607, 0.298, 0.11 },
  { 4.038, 0.348, 0.11 },
  { 4.490, 0.404, 0.11 },
  { 4.966, 0.458, 0.12 },
  { 5.465, 0.516, 0.12 },
  { 5.989, 0.578, 0.11 },
  { 6.539, 0.645, 0.12 },
  { 7.112, 0.713, 0.12 }
};

double myexp(double x)
{
  if (x < -70)
    return 0.0;
  else
    return exp(x);
}

real ptheta_incoh(int Eindex,real theta) 
     /* theta should be in degrees */
{
  /* These numbers generated by fitting 5 gaussians to the real function
   * that describes the probability for theta.
   * We use symmetry in the gaussian (see 180-angle) therefore there
   * are fewer parameters (only 8 per energylevel).
   */
  static double ppp[5][8] = {
    { -0.00295169, 10.4847, 0.0341099, /*-*/43.1963, 
      -0.0164054,  30.2452, 71.0311,    2.50282 },
    { -0.00370852, 9.02037, 0.100559,  /*-*/42.9962,
      -0.0537891,  35.5077, 71.4305,    1.05515 },
    { -0.00427039, 7.86831, 0.118042,  /*-*/45.9846,
      -0.0634505,  38.6134, 70.3857,    0.240082 },
    { -0.00488796, 5.87988, 0.159574,  /*-*/51.5556,
      -0.0855767,  44.7307, 69.0251,   -0.414604 },
    { -0.00504604, 4.56299, 0.201064,  /*-*/54.8599,
      -0.107153,   48.7016, 68.8212,   -0.487699 }
  };
  double g1,g2,g3,g4,g5,ptheta;

  g1 = myexp(-0.5*sqr((theta-ppp[Eindex][7])/ppp[Eindex][1]));
  g2 = myexp(-0.5*sqr((theta-180+ppp[Eindex][7])/ppp[Eindex][1]));
  g3 = myexp(-0.5*sqr((theta-90)/ppp[Eindex][3]));
  g4 = myexp(-0.5*sqr((theta-ppp[Eindex][6])/ppp[Eindex][5]));
  g5 = myexp(-0.5*sqr((theta-180+ppp[Eindex][6])/ppp[Eindex][5]));

  ptheta = ppp[Eindex][0]*(g1+g2) + ppp[Eindex][2]*g3 + ppp[Eindex][4]*(g4+g5);

  return ptheta;
}

real rand_theta_incoh(int Eindex,int *seed) 
{
#define NINTP 450
#define prev (1-cur)
  static bool bFirst = TRUE;
  static real **intp;
  static int  i,j,cur=1;
  real theta,sum,rrr,dx;
  real g[5],y[2];
  
  dx = 90.0/(real)NINTP;
  if (bFirst) {
    /* Compute cumulative integrals of all probability distributions */
    snew(intp,5);
    for(i=0; (i<5); i++) {
      snew(intp[i],NINTP+1);
      y[prev]    = ptheta_incoh(i,0.0);
      /*sum        = y[prev];*/
      for(j=1; (j<=NINTP); j++) {
	y[cur]     = ptheta_incoh(i,j*dx);
	/*sum       += y[cur];*/
	intp[i][j] = intp[i][j-1] + (y[cur]+y[prev])*dx;
	cur        = prev;
      }
    }
    if (debug) {
      fprintf(debug,"Integrated probability functions for theta incoherent\n");
      for(j=0; (j<NINTP); j++) {
	fprintf(debug,"%10f",dx*j);
	for(i=0; (i<5); i++) 
	  fprintf(debug,"  %10f",intp[i][j]);
	fprintf(debug,"\n");
      }
    }
    bFirst = FALSE;
  }

  rrr = rando(seed);
  for(j=0; (j<NINTP) && (rrr > intp[Eindex][j]); j++)
    ;

  return (j-1+(rrr-intp[Eindex][j-1])/(intp[Eindex][j]-intp[Eindex][j-1]))*dx;
}

enum { ecollPHOTO, ecollINELASTIC, ecollNR };

void ionize(FILE *log,t_mdatoms *md,char **atomname[],real t,t_inputrec *ir,
	    rvec v[])
{
  static FILE *xvg,*ion;
  static char *leg[] = { "Probability", "Number/Atom", "Total" };
  static bool bFirst = TRUE,bDoIt;
  static bool bImpulse = TRUE;
  static real pt,t0,imax,width,inv_nratoms,rho,nphot;
  static int  seed,total,ephot,*nelectrons;
  static int  Eindex=-1;
  static int  Energies[] = { 6, 8, 10, 15, 20 };
#define NENER asize(Energies)
  real r,factor,ndv,E_lost,cross_atom,dvz,Kalpha,rrc;
  real incoh,incoh_abs,sigmaPincoh;
  rvec dv;
  bool bIonize;
  char *cc;
  int  i,j,k,m,dq,*elemcnt,nel_atom;
  
  if (bFirst) {
    /* Get parameters for gaussian photon pulse from inputrec */
    t0    = ir->userreal1;  /* Peak of the gaussian pulse            */
    nphot = ir->userreal2;  /* Intensity                             */
    width = ir->userreal3;  /* Width of the peak                     */
    rho   = ir->userreal4;  /* Diameter of the focal spot            */
    seed  = ir->userint1;   /* Random seed for stochastic ionization */
    ephot = ir->userint2;   /* Energy of the photons                 */
    
    imax  = (nphot/(M_PI*rho*rho))*1e-10*sqrt(1.0/M_PI)*(2.0/width);

    bDoIt = ((width > 0) && (imax > 0));
    
    if (bDoIt) {
      if (seed == 0)
	seed = 1993;
      
      for(Eindex=0; (Eindex < NENER) && (Energies[Eindex] != ephot); 
	  Eindex++)
	;
      if (Eindex == NENER)
	fatal_error(0,"Energy level of %d keV not supported",ephot);
      
      total = 0;
      inv_nratoms = 1.0/md->nr;
      snew(nelectrons,md->nr);
      snew(elemcnt,NELEM);
      for(i=0; (i<md->nr); i++) {
	cc = *(atomname[i]);
	for(j=1; (j<NELEM); j++)
	  if (strstr(cc,element[j].name) != NULL) {
	    nelectrons[i]     = element[j].nel;
	    break;
	  }
	if (j == NELEM) {
	  j = 6; /* Carbon, hopefully */
	  nelectrons[i]     = element[j].nel;
	  fprintf(stderr,
		  "Don't know number of electrons for %s, assuming %d\n",
		  cc,nelectrons[i]);
	}
	elemcnt[j]++;
      }
      xvg   = xvgropen("ionize.xvg","Ionization Events","Time (ps)","()");
      xvgr_legend(xvg,asize(leg),leg);
      ion   = ffopen("ionize.log","w");
            
      fprintf(log,"Parameters for ionization events:\n");
      fprintf(log,"Imax = %g, t0 = %g, width = %g, seed = %d\n"
	      "# Photons = %g, rho = %g, ephot = %d (keV), Impulse = %s\n",
	      imax,t0,width,seed,nphot,rho,ephot,yesno_names[bImpulse]);
      fprintf(log,"Electron_mass: %10.3e(keV) Atomic_mass: %10.3e(keV)\n"
	      "Speed_of_light: %10.3e(nm/ps)\n",
	      ELECTRONMASS_keV,ATOMICMASS_keV,SPEEDOFLIGHT);
      fprintf(log,"You have the following elements in your system:\n");
      for(j=1; (j<NELEM); j++)
	if (elemcnt[j] > 0)
	  fprintf(log,"  %s: %d",element[j].name,elemcnt[j]);
      fprintf(log," atoms\n");
    }
    bFirst = FALSE;
  }
  if (!bDoIt)
    return;

  /******************************************************
   *
   *    H E R E    S T A R T S   I O N I Z A T I O N
   *
   ******************************************************/
        
  /* Calculate probability */
  pt          = imax*ir->delta_t*exp(-sqr(2.0*(t-t0)/width));
  dq          = 0;
  
  /* Width of gaussian for probability of incoherent scattering */
  sigmaPincoh = 1/sqrt(44.0);
  
  /* Print time to the file that holds ionization events per atom */
  fprintf(ion,"%12.5f",t);
  
  /* Loop over atoms */
  for(i=0; (i<md->nr); i++) {
    /* Loop over collision types */
    for(k=0; (k<ecollNR); k++) {
      /* Check for electrons */
      nel_atom = nelectrons[i]/*-1*/;
      if (nel_atom >= 0) {
	/* Determine cross section */
	switch (k) {
	case ecollPHOTO: 
	  cross_atom = element[nel_atom].cross[Eindex].photo;
	  break;
	case ecollINELASTIC: 
	  cross_atom = element[nel_atom].cross[Eindex].incoh;
	  break;
	default:
	  fatal_error(0,"Vad?");
	}
	/* Check whether to ionize this guy */
	r = rando(&seed);
	if (r < pt*cross_atom) {
	  
	  clear_rvec(dv);
	  /* Energy of incoming photon, divided by mass of atom,
	   * yields velocity of atom
	   */
	  dv[ZZ] = (ephot*SPEEDOFLIGHT)/(ATOMICMASS_keV*md->massT[i]);
	  
	  switch (k) {
	  case ecollPHOTO: {
	    real theta,phi;
	    
	    /* Get parameters for photoelestic effect */
	    /* Note that in the article this called 2 theta */
	    theta = DEG2RAD*gauss(70.0,26.0,&seed);
	    phi   = 2*M_PI*rando(&seed);
	    
	    /* Now see which electron we kick off (K or L) */
	    rrc        = rando(&seed);
	    if (rrc < recoil[nel_atom].Prob_L)
	      E_lost = ephot-recoil[nel_atom].E_L;
	    else
	      E_lost = ephot-recoil[nel_atom].E_K;
	    if (E_lost < 0) {
	      E_lost  = 0.0;
	      bIonize = FALSE;
	    }
	    else {
	      /* Compute the components of the velocity vector */
	      factor = ((ELECTRONMASS_keV/(ATOMICMASS_keV*md->massT[i]))*
			(SPEEDOFLIGHT*sqrt(2*E_lost/ELECTRONMASS_keV)));
	      
	      /* Subtract momentum of recoiling electron */
	      dv[XX] -= factor*cos(phi)*sin(theta);
	      dv[YY] -= factor*sin(phi)*sin(theta);
	      dv[ZZ] -= factor*cos(theta);

	      if (debug)
		pr_rvecs(debug,0,"PHOTO",&dv,1);
	      	      
	      bIonize = TRUE;
	    }
	      
	    break;
	  }
	  case ecollINELASTIC: {
	    real theta,phi;

	    incoh      = element[nel_atom].cross[Eindex].incoh;
	    incoh_abs  = element[nel_atom].cross[Eindex].incoh_abs;
	    factor     = dv[ZZ]*((incoh-incoh_abs)/incoh);
	    theta      = DEG2RAD*rand_theta_incoh(Eindex,&seed);
	    phi        = 2*M_PI*rando(&seed);
	    
	    /* Subtract momentum of recoiling electron */
	    dv[XX] -= factor*cos(phi)*sin(theta);
	    dv[YY] -= factor*sin(phi)*sin(theta);
	    dv[ZZ] -= factor*cos(theta);

	    if (debug)
	      pr_rvecs(debug,0,"INELL",&dv,1);
	    
	    bIonize = TRUE;
	    
	    break;
	  }
	  }
	  
	  if (bIonize) {
	    /* First decrease the charge */
	    fprintf(ion,"  %d",i);
	    md->chargeA[i] -= 1.0;
	    md->chargeB[i] -= 1.0;
	    dq ++;
	    total ++;
	    
	    if (debug) {
	      fprintf(debug,"Random-dv[%3d] = %10.3e,%10.3e,%10.3e,"
		      " ephot = %d, Elost=%10.3e\n",
		      i,dv[XX],dv[YY],dv[ZZ],ephot,E_lost);
	    }
	    /* Now actually add the impulse to the velocities */
	    for(m=0; (m<DIM); m++)
	      v[i][m] += dv[m];
	  }
	}
      }
    }
  }
  fprintf(ion,"\n");
  fprintf(xvg,"%10g  %10g  %10g  %10g\n",
	  t,pt,dq*inv_nratoms,total*inv_nratoms);
  fflush(ion);
  fflush(xvg);
}
