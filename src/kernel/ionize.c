/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_ionize_c = "$Id$";
#include <string.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "random.h"
#include "physics.h"
#include "xvgr.h"
#include "vec.h"
#include "pbc.h"
#include "txtdump.h"
#include "ionize.h"
#include "names.h"
#include "futil.h"
#include "network.h"
#include "ion_data.h"

#define PREFIX "IONIZE: "

enum { eionCYL, eionSURF, eionGAUSS, eionNR };

enum { ecollPHOTO, ecollINELASTIC, ecollNR };

typedef struct {
  int  z,n,k;
  real fj,sigPh,sigIn,vAuger;
} t_cross_atom;

typedef struct {
  int nelec,maxelec,elec0,elmin_type;
} t_electron_db;

/* BEGIN GLOBAL VARIABLES */
static int   Energies[] = { 6, 8, 10, 12, 15, 20 };
static int   ionize_seed = 1993;
#define NENER asize(Energies)

static t_electron_db edb;

/* END GLOBAL VARIABLES */

void dump_ca(FILE *fp,t_cross_atom *ca,int i,char *file,int line)
{
  fprintf(fp,PREFIX"(line %d) atom %d, z = %d, n = %d, k = %d\n",
	  line,i,ca->z,ca->n,ca->k);
}

t_cross_atom *mk_cross_atom(FILE *log,t_mdatoms *md,
			    char **atomname[],int Eindex)
{
  int elem_index[] = { 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 5, 6 };
  t_cross_atom *ca;
  int  *elemcnt;
  char *cc;
  int  i,j;
  
  fprintf(log,PREFIX"Filling data structure for ionization\n");
  fprintf(log,PREFIX"Warning: all fj values set to 0.95 for now\n");
  snew(ca,md->nr);
  snew(elemcnt,NELEM+1);
  for(i=0; (i<md->nr); i++) {
    ca[i].n = 0;
    ca[i].k = 0;
    cc = *(atomname[i]);
    for(j=0; (j<NELEM); j++)
      if (strncmp(cc,element[j].name,strlen(element[j].name)) == 0) {
	ca[i].z = element[j].nel;
	break;
      }
    if (j == NELEM) 
      fatal_error(0,PREFIX"Don't know number of electrons for %s",
		  *atomname[i]);
    elemcnt[j]++;

    ca[i].sigPh = element[elem_index[ca[i].z]].cross[Eindex].photo;
    ca[i].sigIn = element[elem_index[ca[i].z]].cross[Eindex].incoh;
    ca[i].fj    = recoil[ca[i].z].Prob_K;
    switch (ca[i].z) {
    case 6:
      ca[i].vAuger  = 0.904;
      break;
    case 7:
      ca[i].vAuger  = 0.920;
      break;
    case 8:
      ca[i].vAuger  = 0.929;
      break;
    case 16:
    case 20:
      ca[i].vAuger = 1.0;
      break;
    default:
      ca[i].vAuger= -1;
    }
  }
  
  fprintf(log,PREFIX"You have the following elements in your system (%d atoms):\n"PREFIX,md->nr);
  for(j=0; (j<NELEM); j++)
    if (elemcnt[j] > 0)
      fprintf(log,"  %s: %d",element[j].name,elemcnt[j]);
  fprintf(log," atoms\n");

  sfree(elemcnt);
  
  return ca;
}

int number_K(t_cross_atom *ca)
{
  if (ca->z <= 2)
    return ca->z-ca->n;
  else
    return 2-ca->k;
}

int number_L(t_cross_atom *ca)
{
  return ca->k-2+ca->z-ca->n;
}

real xray_cross_section(int eColl,t_cross_atom *ca)
{
  real c=0;
  int  nK,nL;
  
  switch (eColl) {
  case ecollPHOTO:
    nK = number_K(ca);
    nL = number_L(ca);
    if (ca->z == 1)
      c = ca->sigPh;
    else if (ca->z == 2)
      c = ca->sigPh*0.5;
    else
      c = (nK*0.5*ca->fj + nL/(ca->z-2)*(1-ca->fj))*ca->sigPh;
    break;
  case ecollINELASTIC:
    c = (ca->z-ca->n)*ca->sigIn/ca->z;
    break;
  default:
    fatal_error(0,"No such collision type %d\n",eColl);
  }
  return c;
}

real prob_K(int eColl,t_cross_atom *ca)
{
  real Pl,Pk,P=0;
  
  if ((ca->z <= 2) || (ca->z == ca->n))
    return 0;

  switch (eColl) {
  case ecollPHOTO:
    Pl = (ca->k-2+ca->z-ca->n)*(1-ca->fj)/(ca->z-2);
    Pk = (2-ca->k)*ca->fj*0.5;
    P  = Pk/(Pl+Pk);
    break;
  case ecollINELASTIC:
    P = (2-ca->k)/(ca->z-ca->n);
    break;
  default:
    fatal_error(0,"No such collision type %d\n",eColl);
  } 
  return P;
}

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
  static double ppp[NENER][8] = {
    { -0.00295169, 10.4847, 0.0341099, /*-*/43.1963, 
      -0.0164054,  30.2452, 71.0311,    2.50282 },
    { -0.00370852, 9.02037, 0.100559,  /*-*/42.9962,
      -0.0537891,  35.5077, 71.4305,    1.05515 },
    { -0.00427039, 7.86831, 0.118042,  /*-*/45.9846,
      -0.0634505,  38.6134, 70.3857,    0.240082 },
    { -0.004514,   7.0728,  0.13464,  /*-*/48.213,
      -0.0723,     41.06,   69.38,     -0.02 },
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
  real g[NENER],y[2];
  
  dx = 90.0/(real)NINTP;
  if (bFirst) {
    /* Compute cumulative integrals of all probability distributions */
    snew(intp,NENER);
    for(i=0; (i<NENER); i++) {
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
	for(i=0; (i<NENER); i++) 
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

static void polar2cart(real phi,real theta,rvec v)
{
  v[XX] = cos(phi)*sin(theta);
  v[YY] = sin(phi)*sin(theta);
  v[ZZ] = cos(theta);
}

void rand_vector(rvec v,int *seed)
{
  real theta,phi;

  theta = 180.0*rando(seed);
  phi   = 360.0*rando(seed);
  polar2cart(phi,theta,v);
}

real electron_cross_section(FILE *fp,rvec v,real mass,int nelec)
{
  /* Compute cross section for electrons */
  real T,B,U,S,Q,R,N,t,u,lnt,sigma;
  real a0 = 0.05292; /* nm */
  
  /* Have to determine T (kinetic energy of electron) */
  T = 0.5*mass*iprod(v,v);
  
  /* R is the binding energy of the electron in hydrogen */
  R = 13.61*ELECTRONVOLT;
  
  /* Have to determine the binding energy B, differs per orbital of course */
  B = R;
  
  /* Have to determine the orbital kinetic energy U */
  U = R;
  
  /* Have to know number of electrons */
  N = nelec;
  
  /* Magic constant Q */
  Q = 1;
  
  /* Some help variables */
  t     = T/B;
  u     = U/B;
  S     = 4*M_PI*sqr(a0)*N*sqr(R/B);
  lnt   = log(t);
  
  /* Resulting variable */
  sigma = (S/(t+u+1))*( 0.5*Q*lnt*(1-1/sqr(t)) + (2-Q)*(1-1/t-lnt/(t+1)) ); 
  
  return sigma;
}

static bool analyze_electrons(FILE *fp,t_electron_db *edb,
			      int natom,char **atomname[])
{
  int  i,etp;
  char *cc;
 
  if (((cc = getenv("GENERATE_ELECTRONS")) != NULL) &&
      (sscanf(cc,"%d",&etp) == 1)) {
    for(i=0; (i<natom); i++) {
      if (strcmp(*atomname[i],"EL") == 0)
	break;
    }
    edb->elec0      = i;
    edb->nelec      = 0;
    edb->maxelec    = natom-i;
    edb->elmin_type = etp;
    fprintf(fp,PREFIX"There are %d possible electrons\n",edb->maxelec);
    
    return TRUE;
  }
  else {
    fprintf(fp,PREFIX"No electron features today.\n");
    return FALSE;
  }
}

void add_electron(FILE *fp,t_mdatoms *md,t_electron_db *edb,int ion,
		  rvec x[],rvec v[],rvec dv,real dt)
{
  int  m,ee;
  real nv;
  
  if (edb->nelec < edb->maxelec) {
    ee = edb->elec0+edb->nelec++;
    md->chargeA[ee] = md->chargeB[ee] = md->chargeT[ee] = -1;
    md->typeA[ee]   = md->typeB[ee]   = edb->elmin_type;

    /* Velocity! */
    svmul(-md->massA[ion]*md->invmass[ee],dv,v[ee]);
    /* Do a first step to prevent the electron from being on top of the 
     * nucleus, move it 0.05 nm from the nucleus 
     */
    nv = 1.0/norm(v[ee]);
    for(m=0; (m<DIM); m++) 
      x[ee][m] = x[ion][m] + v[ee][m]*nv*0.05;
  } 
  else
    fatal_error(0,PREFIX"No more particles to turn into electrons\n");
}

bool khole_decay(FILE *fp,t_cross_atom *ca,rvec x[],rvec v[],int ion,
		 int *seed,real dt,bool bElectron,
		 t_mdatoms *md,t_electron_db *edb)
{
  rvec dv;
  real ndv,factor;
  int  m;
  
  if ((ca->vAuger < 0) || (recoil[ca->z].tau == 0)) {
    dump_ca(stderr,ca,ion,__FILE__,__LINE__);
    exit(1);
  }
  if (rando(seed) < dt/recoil[ca->z].tau) {
    if (debug)
      fprintf(debug,"DECAY: Going to decay a k hole\n");
    ca->n++;
    ca->k--;
    /* Generate random vector */
    rand_vector(dv,seed);

    factor = ca->vAuger;
    if (debug)
      fprintf(debug,"DECAY: factor=%10g, dv = (%8.3f, %8.3f, %8.3f)\n",
	      factor,dv[XX],dv[YY],dv[ZZ]);
    svmul(factor,dv,dv);
    rvec_inc(v[ion],dv);

    /* Now put the electron in place */
    if (bElectron)    
      add_electron(fp,md,edb,ion,x,v,dv,dt);

    return TRUE;
  }
  else
    return FALSE;
}

real electron_atom_interactions(FILE *fp,t_mdatoms *md,t_inputrec *ir,
				int start,int end,
				rvec x[],rvec v[],rvec f[],matrix box)
{
  /* Compute what the name says... */
  int  i,j,m,elec1,e1;
  rvec dx;
  real mindist2,vc,vtot,fscal,fc,dx2,dx_1,qi,*q;
  
  mindist2 = sqr(0.05);
  vtot     = 0;
  
  if (edb.nelec > 0) {
    /* Do a search... */
    q = md->chargeT;
    if (ir->ePBC != epbcNONE) 
      init_pbc(box,FALSE);
    /* the end variable usually includes electrons */
    e1 = min(end,edb.elec0);
    for(i=start; (i<e1); i++) {
      elec1 = edb.elec0 + edb.nelec;
      qi = q[i]*ONE_4PI_EPS0;
      for(j=edb.elec0; (j<elec1); j++) {
	if (ir->ePBC == epbcNONE) 
	  rvec_sub(x[i],x[j],dx);
	else
	  pbc_dx(x[i],x[j],dx);
	dx2 = iprod(dx,dx);
	if (dx2 < mindist2) {
	  /* Model collision */
	}
	else {
	  /* Do normal coulomb */
	  dx_1  = invsqrt(dx2);
	  vc    = qi*q[j]*dx_1;
	  vtot += vc;
	  fscal = vc*dx_1*dx_1;
	  for(m=0; (m<DIM); m++) {
	    fc       = fscal*dx[m];
	    f[i][m] += fc;
	    f[j][m] -= fc;
	  }
	}
      }
    }
  }
  return vtot;
}

void ionize(FILE *fp,t_mdatoms *md,char **atomname[],real t,t_inputrec *ir,
	    rvec x[],rvec v[],int start,int end,matrix box,t_commrec *cr)
{
  static FILE  *xvg,*ion;
  static char  *leg[] = { "Probability", "Primary Ionization", "Integral over PI", "KHole-Decay", "Integral over KD" };
  static bool  bFirst = TRUE,bElectron = FALSE;
  static real  t0,imax,width,inv_nratoms,rho,nphot;
  static real  interval;
  static int   dq_tot,nkd_tot,ephot,mode;
  static t_cross_atom *ca;
  static int   Eindex=-1;
    
  real r,factor,ndv,E_lost=0,cross_atom,dvz,rrc;
  real pt,ptot,pphot,pcoll[ecollNR],tmax;
  real incoh,incoh_abs,sigmaPincoh,hboxx,hboxy,rho2;
  rvec dv,ddv;
  bool bIonize=FALSE,bKHole,bL,bDOIT;
  char *cc;
  int  i,j,k,kk,m,nK,nL,dq,nkh,nkdecay,elmin_type;
  int  *nionize,*nkhole,*ndecay,nbuf[2];
  
  if (bFirst) {
    /* Get parameters for gaussian photon pulse from inputrec */
    t0       = ir->userreal1;  /* Peak of the gaussian pulse            */
    nphot    = ir->userreal2;  /* Intensity                             */
    width    = ir->userreal3;  /* Width of the peak (in time)           */
    rho      = ir->userreal4;  /* Diameter of the focal spot (nm)       */
    ionize_seed = ir->userint1;   /* Random seed for stochastic ionization */
    ephot    = ir->userint2;   /* Energy of the photons                 */
    mode     = ir->userint3;   /* Mode of ionizing                      */
    interval = 0.001*ir->userint4;   /* Interval between pulses (ps)    */
     
    if ((width <= 0) || (nphot <= 0))
      fatal_error(0,"Your parameters for ionization are not set properly\n"
		  "width (userreal3) = %f,  nphot (userreal2) = %f",
		  width,nphot);
    
    if ((mode < 0) || (mode >= eionNR))
      fatal_error(0,"Ionization mode (userint3)"
		  " should be in the range 0 .. %d",eionNR-1);
    
    switch (mode) {
    case eionCYL:
      imax  = (nphot/(M_PI*sqr(rho/2)))*1e-10*1.0/(width*sqrt(2.0*M_PI));
      break;
    case eionSURF:
      imax  = (nphot/(M_PI*sqr(rho/2)))*1e-10*1.0/(width*sqrt(2.0*M_PI));
      break;
    }
    if (ionize_seed == 0)
      ionize_seed = make_seed();
    if (PAR(cr)) {
      for(i=0; (i<cr->nodeid); i++)
	ionize_seed = INT_MAX*rando(&ionize_seed);
      fprintf(fp,PREFIX"Modifying seed on parallel processor to %d\n",
	      ionize_seed);
    }
          
    for(Eindex=0; (Eindex < NENER) && (Energies[Eindex] != ephot); Eindex++)
      ;
    if (Eindex == NENER)
      fatal_error(0,PREFIX"Energy level of %d keV not supported",ephot);
    
    /* Initiate cross section data etc. */
    ca      = mk_cross_atom(fp,md,atomname,Eindex);
    
    dq_tot  = 0;
    nkd_tot = 0;
    inv_nratoms = 1.0/md->nr;

    xvg   = xvgropen("ionize.xvg","Ionization Events","Time (ps)","()");
    xvgr_legend(xvg,asize(leg),leg);
    ion   = ffopen("ionize.log","w");

    bElectron = analyze_electrons(fp,&edb,md->nr,atomname);
    
    fprintf(fp,PREFIX"Parameters for ionization events:\n");
    fprintf(fp,PREFIX"Imax = %g, t0 = %g, width = %g, seed = %d\n"
	    PREFIX"# Photons = %g, rho = %g, ephot = %d (keV), Electrons = %s\n",
	    imax,t0,width,ionize_seed,nphot,rho,ephot,yesno_names[bElectron]);
    fprintf(fp,PREFIX"Electron_mass: %10.3e(keV) Atomic_mass: %10.3e(keV)\n"
	    PREFIX"Speed_of_light: %10.3e(nm/ps)\n",
	    ELECTRONMASS_keV,ATOMICMASS_keV,SPEED_OF_LIGHT);
    fprintf(fp,PREFIX"Interval between shots: %g ps\n",interval);
    fprintf(fp,PREFIX"Eindex = %d\n",Eindex);
    fprintf(fp,PREFIX"Doing ionizations for atoms %d - %d\n",start,end);
    
    fflush(fp);

    bFirst = FALSE;
  }

  /******************************************************
   *
   *    H E R E    S T A R T S   I O N I Z A T I O N
   *
   ******************************************************/

  /* Calculate probability */
  tmax        = t0;
  if (interval > 0)
    while (t > (tmax+interval*0.5))
      tmax += interval;
  /*  End when t <= t0 + (N+0.5) interval */
  
  pt          = imax*ir->delta_t*exp(-0.5*sqr((t-tmax)/width));
  dq          = 0;
  nkdecay     = 0;

  hboxx       = 0.5*box[XX][XX];
  hboxy       = 0.5*box[YY][YY];
  rho2        = sqr(rho);
  
  /* Width of gaussian for probability of incoherent scattering */
  sigmaPincoh = 1/sqrt(44.0);

  /* Arrays for ionization statistics */
  snew(nionize,md->nr);
  snew(nkhole,md->nr);
  snew(ndecay,md->nr);
    
  /* Loop over atoms */
  for(i=start; (i<end); i++) {
    /* Loop over collision types */
    bKHole = FALSE;
    for(k=0; (k<ecollNR); k++) 
      /* Determine cross section for this collision type */
      pcoll[k]= pt*xray_cross_section(k,&(ca[i]));
    
    /* Total probability of ionisation */
    ptot = 1 - (1-pcoll[ecollPHOTO])*(1-pcoll[ecollINELASTIC]);
    if (debug && (i==0)) 
      fprintf(debug,PREFIX"Ptot = %g, t = %g\n",ptot,t);
    
    /* Check whether to ionize this guy */
    bDOIT = FALSE;
    switch (mode) {
    case eionCYL:
      bDOIT = (((rando(&ionize_seed) < ptot) && (ca[i].n < ca[i].z)) && 
	       ((sqr(x[i][XX] - hboxx) + sqr(x[i][YY] - hboxy)) < rho2));
      break;
    case eionSURF:
      bDOIT = FALSE;
      break;
    default:
      fatal_error(0,"Unknown ionization mode %d (%s, line %d)",mode,
		  __FILE__,__LINE__);
    }
      
    if (bDOIT) {
      clear_rvec(dv);
      
      /* The relative probability for a photoellastic event is given by: */
      pphot = pcoll[ecollPHOTO]/(pcoll[ecollPHOTO]+pcoll[ecollINELASTIC]);
      
      if (rando(&ionize_seed) < pphot) 
	k = ecollPHOTO;
      else
	k = ecollINELASTIC;
      
      /* If a random number is smaller than the probability for 
       * an L ionization than do that. Note that the probability
       * may be zero (H, He), but the < instead of <= covers that.
       */
      nK = number_K(&ca[i]);
      nL = number_L(&ca[i]);
      bL = (nK == 0) || ( (nL > 0) && (rando(&ionize_seed) > prob_K(k,&(ca[i]))));

      switch (k) {
      case ecollPHOTO: {
	/* Select which one to take by yet another random numer */
	real theta,phi;
	
	/* Get parameters for photoelestic effect */
	/* Note that in the article this is called 2 theta */
	theta = DEG2RAD*gauss(70.0,26.0,&ionize_seed);
	phi   = 2*M_PI*rando(&ionize_seed);
	
	if (bL)
	  E_lost = ephot-recoil[ca[i].z].E_L*(ca[i].n+1);
	else {
	  E_lost = ephot-recoil[ca[i].z].E_K;
	  if ((ca[i].z > 2) && (nL > 0))
	    bKHole = TRUE;
	}
	if (debug)
	  fprintf(debug,"i = %d, nK = %d, nL = %d, bL = %s, bKHole = %s\n",
		  i,nK,nL,BOOL(bL),BOOL(bKHole));
	if (E_lost < 0) {
	  E_lost  = 0.0;
	  bIonize = FALSE;
	  bKHole  = FALSE;
	}
	else {
	  /* Compute the components of the velocity vector */
	  factor = ((ELECTRONMASS_keV/(ATOMICMASS_keV*md->massT[i]))*
		    (SPEED_OF_LIGHT*sqrt(2*E_lost/ELECTRONMASS_keV)));
	  
	  /* Subtract momentum of recoiling electron */
	  polar2cart(phi,theta,ddv);
	  for(m=0; (m<DIM); m++)
	    dv[m] -= factor*ddv[m];
	
	  if (bElectron)
	    add_electron(fp,md,&edb,i,x,v,dv,ir->delta_t);
	  
	  if (debug)
	    pr_rvec(debug,0,"ELL",dv,DIM);
	  
	  bIonize = TRUE;
	}
	break;
      }
      case ecollINELASTIC: {
	real theta,phi,Ebind,Eelec;
	
	if (bL)
	  Ebind = (ca[i].n+1)*recoil[ca[i].z].E_L;
	else {
	  Ebind  = recoil[ca[i].z].E_K;
	  if ((ca[i].z > 2) && (nL > 0))
	    bKHole = TRUE;
	}
	theta      = DEG2RAD*rand_theta_incoh(Eindex,&ionize_seed);
	Eelec      = (sqr(ephot)/512)*(1-cos(2*theta));
	bIonize    = (Ebind <= Eelec);
	bKHole     = bKHole && bIonize;
	if (debug)
	  fprintf(debug,PREFIX"Ebind: %g, Eelectron: %g\n",Ebind,Eelec);
	if (!bIonize) {
	  /* Subtract momentum of recoiling photon */
	  /*phi     = 2*M_PI*rando(&ionize_seed);
 	    bKHole  = FALSE;  
	    factor  = ephot*438;  
	    dv[XX] -= factor*cos(phi)*sin(theta);
	    dv[YY] -= factor*sin(phi)*sin(theta);
	    dv[ZZ] -= factor*cos(theta);
	  */
	  if (debug)
	    pr_rvec(debug,0,"INELL",dv,DIM);
	}
	break;
      }
      default:
	fatal_error(0,"Ga direct naar de gevangenis. Ga niet langs start");
      }
      if (bIonize) {
	/* First increase the charge */
	if (ca[i].n < ca[i].z) {
	  md->chargeA[i] += 1.0;
	  md->chargeB[i] += 1.0;
	  ca[i].n++;
	  dq ++;
	}
	if (debug) {
	  fprintf(debug,"Random-dv[%3d] = %10.3e,%10.3e,%10.3e,"
		  " ephot = %d, Elost=%10.3e\n",
		  i,dv[XX],dv[YY],dv[ZZ],ephot,E_lost);
	}
      }
      /* Now actually add the impulse to the velocities */
      for(m=0; (m<DIM); m++)
	v[i][m] += dv[m];
      if (bKHole) {
	ca[i].k ++;
	nkhole[i]++;
      }
      else if (bIonize)
	nionize[i]++;
    }
    
    /* Now check old event: Loop over k holes! */
    nkh = ca[i].k;
    for (kk = 0; (kk < nkh); kk++) 
      if (khole_decay(fp,&(ca[i]),x,v,i,&ionize_seed,ir->delta_t,
		      bElectron,md,&edb)) {
	nkdecay ++;
	ndecay[i]++;
      }
    
    if (debug && (ca[i].n > 0))
      dump_ca(debug,&(ca[i]),i,__FILE__,__LINE__);
  }

  /* Sum events for statistics if necessary */
  if (PAR(cr)) {
    gmx_sumi(md->nr,nionize,cr);
    gmx_sumi(md->nr,nkhole,cr);
    gmx_sumi(md->nr,ndecay,cr);
    nbuf[0] = dq; nbuf[1] = nkdecay;
    gmx_sumi(2,nbuf,cr);
    dq = nbuf[0]; nkdecay = nbuf[1];
  }
  /* Now sum global events on this timestep to cumulative numbers */
  dq_tot  += dq;
  nkd_tot += nkdecay;
  
  /* Printing time */
  if (MASTER(cr)) {
    /* Print data to the file that holds ionization events per atom */
    fprintf(ion,"%12.8f",t);
    for(i=0; (i<md->nr); i++) {
      if (nionize[i])
	fprintf(ion,"  I:%d",i+1);
      if (nkhole[i])
	fprintf(ion,"  K:%d",i+1);
      if (ndecay[i])
	fprintf(ion,"  D:%d",i+1);
    }
    fprintf(ion,"\n");
    if (debug)
      fflush(ion);
  
    /* Print statictics to file */
    fprintf(xvg,"%10.5f  %10.3e  %6d  %6d  %6d  %6d",
	    t,pt,dq,dq_tot,nkdecay,nkd_tot);
    fprintf(xvg,"\n");
    if (debug)
      fflush(xvg);
  }
  sfree(nionize);
  sfree(nkhole);
  sfree(ndecay);
}

