/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_update_c = "$Id$";

#include <stdio.h>
#include <math.h>

#include "assert.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "nrnb.h"
#include "physics.h"
#include "macros.h"
#include "vveclib.h"
#include "vec.h"
#include "main.h"
#include "confio.h"
#include "update.h"
#include "random.h"
#include "futil.h"
#include "mshift.h"
#include "tgroup.h"
#include "force.h"
#include "names.h"
#include "txtdump.h"
#include "mdrun.h"
#include "copyrite.h"
#include "constr.h"
#include "edsam.h"
#include "pull.h"

typedef struct {
  real gdt;
  real eph;
  real emh;
  real em;
  real b;
  real c;
  real d;
} t_sdconst;

/* constants for a (bad) random number generator */
const unsigned long im = 0xffff;
const unsigned long ia = 1093;
const unsigned long ic = 18257;
const real inv_im      = 1.0/(0xffff);

static t_sdconst *sdc;

static void do_update_md(int start,int homenr,double dt,
			 rvec lamb[],t_grp_acc *gstat,t_grp_tcstat *tcstat,
			 rvec accel[],ivec nFreeze[],real invmass[],
			 unsigned short ptype[],unsigned short cFREEZE[],
			 unsigned short cACC[],unsigned short cTC[],
			 rvec x[],rvec xprime[],rvec v[],rvec vold[],
			 rvec f[],matrix M,bool bExtended)
{
  double imass,w_dt;
  int    gf,ga,gt;
  rvec   vrel;
  real   vn,vv,va,vb,vnrel;
  real   lg,xi,uold;
  int    n,d;

  if(bExtended) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parinello-Rahman coupling
     */
    for (n=start; n<start+homenr; n++) {  
      imass = invmass[n];
      gf   = cFREEZE[n];
      ga   = cACC[n];
      gt   = cTC[n];
      xi   = tcstat[gt].xi;
      
      rvec_sub(v[n],gstat[ga].uold,vrel);

      for (d=0; d<DIM; d++) {
	lg             = lamb[gt][d]; 
	vold[n][d]     = v[n][d];
	
	if ((ptype[n] != eptDummy) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
	  
	  vnrel= lg*(vrel[d] + dt*(imass*f[n][d]-xi*vrel[d]-iprod(M[d],vrel)));  
	  /* do not scale the mean velocities u */
	  vn             = gstat[ga].uold[d] + accel[ga][d]*dt + vnrel; 
	  v[n][d]        = vn;
	  xprime[n][d]   = x[n][d]+vn*dt;
	} else
	  xprime[n][d]   = x[n][d];
      }
    }
    
  } else {
    /* Classic version of update, used with berendsen coupling */
    for (n=start; n<start+homenr; n++) {  
      w_dt = invmass[n]*dt;
      gf   = cFREEZE[n];
      ga   = cACC[n];
      gt   = cTC[n];
      
      for (d=0; d<DIM; d++) {
	vn             = v[n][d];
	lg             = lamb[gt][d];
	vold[n][d]     = vn;
	
	if ((ptype[n] != eptDummy) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
	  vv             = lg*(vn + f[n][d]*w_dt);
          
	  /* do not scale the mean velocities u */
	  uold           = gstat[ga].uold[d];
	  va             = vv + accel[ga][d]*dt;
	  vb             = va + (1.0-lg)*uold;
	  v[n][d]        = vb;
	  xprime[n][d]   = x[n][d]+vb*dt;
	} else {
	  v[n][d]        = 0.0;
	  xprime[n][d]   = x[n][d];
	}
      }
    }
  }
}

static void do_update_visc(int start,int homenr,double dt,
			   rvec lamb[],real invmass[],t_grp_tcstat *tcstat,
			   unsigned short ptype[],unsigned short cTC[],
			   rvec x[],rvec xprime[],rvec v[],rvec vold[],
			   rvec f[],matrix M,matrix box,real
			   cos_accel,real vcos,bool bExtended)
{
  double imass,w_dt;
  int    gt;
  real   vn,vc;
  real   lg,xi,vv;
  real   fac,cosz;
  rvec   vrel;
  int    n,d;
  
  fac = 2*M_PI/(box[ZZ][ZZ]);

  if(bExtended) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parinello-Rahman coupling
     */
    for (n=start; n<start+homenr; n++) {  
      imass = invmass[n];
      gt   = cTC[n];
      cosz = cos(fac*x[n][ZZ]);
      
      copy_rvec(v[n],vold[n]);
      copy_rvec(v[n],vrel);
      
      vc            = cosz*vcos;
      vrel[XX]     -= vc;
      xi           = tcstat[gt].xi;
      
      for (d=0; d<DIM; d++) {
	vn             = v[n][d];
	lg             = lamb[gt][d];
	
	if ((ptype[n] != eptDummy) && (ptype[n] != eptShell)) {
	  vn              = lg*(vrel[d] + dt*(imass*f[n][d]-xi*vrel[d]-iprod(M[d],vrel)));
	  if (d == XX) 
	    vn           += vc + dt*cosz*cos_accel;
	  
	  v[n][d]        = vn;
	  xprime[n][d]   = x[n][d]+vn*dt;
	} else
	  xprime[n][d]   = x[n][d];
      }
    }
    
  } else {
    /* Classic version of update, used with berendsen coupling */
    for (n=start; n<start+homenr; n++) {  
      w_dt = invmass[n]*dt;
      gt   = cTC[n];
      cosz = cos(fac*x[n][ZZ]);
      
      for (d=0; d<DIM; d++) {
	vn             = v[n][d];
	lg             = lamb[gt][d];
	vold[n][d]     = vn;
	
	if ((ptype[n] != eptDummy) && (ptype[n] != eptShell)) {
	  if (d == XX) {
	    vc           = cosz*vcos;
	    /* Do not scale the cosine velocity profile */
	    vv           = vc + lg*(vn - vc + f[n][d]*w_dt);
	    /* Add the cosine accelaration profile */
	    vv          += dt*cosz*cos_accel;
	  } else
	    vv           = lg*(vn + f[n][d]*w_dt);
	  
	  v[n][d]        = vv;
	  xprime[n][d]   = x[n][d]+vv*dt;
	} else {
	  v[n][d]        = 0.0;
	  xprime[n][d]   = x[n][d];
	}
      }
    }
  }
}

static real fgauss(unsigned long *jran)
{
  static real sqrt3 = 1.7320508075688772;
  real jr;

  *jran = (*jran*ia+ic) & im;
  jr = (real)*jran;
  *jran = (*jran*ia+ic) & im;
  jr += (real)*jran;
  *jran = (*jran*ia+ic) & im;
  jr += (real)*jran;
  *jran = (*jran*ia+ic) & im;
  jr += (real)*jran;
  
  return sqrt3*(jr*inv_im-2);
}

void init_sd_consts(int ngtc,real tau_t[],real dt)
{
  int  n;
  real y;

  snew(sdc,ngtc);

  for(n=0; n<ngtc; n++) {
    sdc[n].gdt = dt/tau_t[n];
    sdc[n].eph = exp(sdc[n].gdt/2);
    sdc[n].emh = exp(-sdc[n].gdt/2);
    sdc[n].em  = exp(-sdc[n].gdt);
    if (sdc[n].gdt >= 0.1) {
      sdc[n].b = sdc[n].gdt*(sqr(sdc[n].eph)-1) - 4*sqr(sdc[n].eph-1);
      sdc[n].c = sdc[n].gdt - 3 + 4*sdc[n].emh - sdc[n].em;
      sdc[n].d = 2 - sdc[n].eph - sdc[n].emh;
    } else {
      y = sdc[n].gdt/2;
      /* Seventh order expansions for small y */
      sdc[n].b = y*y*y*y*(1/3.0+y*(1/3.0+y*(17/90.0+y*7/9.0)));
      sdc[n].c = y*y*y*(2/3.0+y*(-1/2.0+y*(7/30.0+y*(-1/12.0+y*31/1260.0))));
      sdc[n].d = y*y*(-1+y*y*(-1/12.0-y*y/360.0));
    }
    if (debug)
      fprintf(debug,"SD const tc-grp %d: b %g  c %g  d %g\n",
	      n,sdc[n].b,sdc[n].c,sdc[n].d);
  }
}

static void do_update_sd(int start,int homenr,
			 rvec accel[],ivec nFreeze[],
			 real invmass[],unsigned short ptype[],
			 unsigned short cFREEZE[],unsigned short cACC[],
			 unsigned short cTC[],real SAfactor,
			 rvec x[],rvec xprime[],rvec v[],rvec vold[],rvec f[],
			 int ngtc,real tau_t[],real ref_t[],
			 int *seed, bool bFirstHalf)
{
  typedef struct {
    real V;
    real X;
    real Yv;
    real Yx;
  } t_sd_sigmas;

  static bool bFirst = TRUE;
  static t_sd_sigmas *sig=NULL;
  static rvec *X,*V;
  real   kT;
  int    gf,ga,gt;
  real   vn=0,Vmh,Xmh;
  real   ism;
  int    n,d;
  unsigned long  jran;
  
  if (sig == NULL) {
    snew(sig,ngtc);
    snew(X,homenr);
    snew(V,homenr);
  }
  
  if (bFirstHalf) {
    for(n=0; n<ngtc; n++) {
      kT = BOLTZ*SAfactor*ref_t[n];
      /* The mass is encounted for later, since this differs per atom */
      sig[n].V  = sqrt(kT*(1-sdc[n].em));
      sig[n].X  = sqrt(kT*sqr(tau_t[n])*sdc[n].c);
      sig[n].Yv = sqrt(kT*sdc[n].b/sdc[n].c);
      sig[n].Yx = sqrt(kT*sqr(tau_t[n])*sdc[n].b/(1-sdc[n].em));
    }
  }

  jran = (unsigned long)((real)im*rando(seed));

  for (n=start; n<start+homenr; n++) {  
    ism = sqrt(invmass[n]);
    gf  = cFREEZE[n];
    ga  = cACC[n];
    gt  = cTC[n];
    
    for (d=0; d<DIM; d++) {
      if (bFirstHalf) {
	vn             = v[n][d];
	vold[n][d]     = vn;
      }
      if ((ptype[n] != eptDummy) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
	if (bFirstHalf) {

	  if (bFirst)
	    X[n-start][d] = ism*sig[gt].X*fgauss(&jran);
	  
	  Vmh = X[n-start][d]*sdc[gt].d/(tau_t[gt]*sdc[gt].c) 
	    + ism*sig[gt].Yv*fgauss(&jran);
	  V[n-start][d] = ism*sig[gt].V*fgauss(&jran);
	  
	  v[n][d] = vn*sdc[gt].em 
	    + (invmass[n]*f[n][d] + accel[ga][d])*tau_t[gt]*(1 - sdc[gt].em)
	    + V[n-start][d] - sdc[gt].em*Vmh;
	    
	  xprime[n][d] = x[n][d] + v[n][d]*tau_t[gt]*(sdc[gt].eph - sdc[gt].emh); 
  
	} else {
	  
	  /* Correct the velocties for the constraints */
	  v[n][d] = 
	    (xprime[n][d] - x[n][d])/(tau_t[gt]*(sdc[gt].eph - sdc[gt].emh));  

	  Xmh = V[n-start][d]*tau_t[gt]*sdc[gt].d/(sdc[gt].em-1) 
	    + ism*sig[gt].Yx*fgauss(&jran);
	  X[n-start][d] = ism*sig[gt].X*fgauss(&jran);
	  
	  xprime[n][d] += X[n-start][d] - Xmh;
	  
	}
      } else {
	if (bFirstHalf) {
	  v[n][d]        = 0.0;
	  xprime[n][d]   = x[n][d];
	}
      }
    }
  }
  
  bFirst = FALSE;
}

static void do_update_bd(int start,int homenr,double dt,
			 ivec nFreeze[],
			 real invmass[],unsigned short ptype[],
			 unsigned short cFREEZE[],unsigned short cTC[],
			 rvec x[],rvec xprime[],rvec v[],rvec vold[],
			 rvec f[],real temp,real fr,
			 int ngtc,real tau_t[],real ref_t[],
			 int *seed)
{
  int    gf,gt;
  real   vn;
  static real *rf=NULL;
  real   rfac=0,invfr=0;
  int    n,d;
  unsigned long  jran;
  
  if (rf == NULL)
    snew(rf,ngtc);

  if (fr) {
    rfac  = sqrt(2.0*BOLTZ*temp/(fr*dt));
    invfr = 1.0/fr;
  } else
    for(n=0; n<ngtc; n++)
      rf[n]  = sqrt(2.0*BOLTZ*ref_t[n]);
  
  jran = (unsigned long)((real)im*rando(seed));

  for (n=start; (n<start+homenr); n++) {  
    gf = cFREEZE[n];
    gt = cTC[n];
    for (d=0; (d<DIM); d++) {
      vold[n][d]     = v[n][d];
      if ((ptype[n]!=eptDummy) && (ptype[n]!=eptShell) && !nFreeze[gf][d]) {
	if (fr)
	  vn         = invfr*f[n][d] + rfac*fgauss(&jran);
	else
	  /* NOTE: invmass = 1/(mass*fric_const) */
	  vn         = invmass[n]*f[n][d]*dt 
	               + sqrt(invmass[n])*rf[gt]*fgauss(&jran);

	v[n][d]      = vn;
	xprime[n][d] = x[n][d]+vn*dt;
      } else {
	v[n][d]      = 0.0;
	xprime[n][d] = x[n][d];
      }
    }
  }
}

static void shake_calc_vir(FILE *log,int nxf,rvec x[],rvec f[],tensor vir,
                           t_commrec *cr)
{
  int    i,m,n;
  matrix dvir;
  
  clear_mat(dvir);
  for(i=0; (i<nxf); i++) {
    for(m=0; (m<DIM); m++)
      for(n=0; (n<DIM); n++)
        dvir[m][n]+=x[i][m]*f[i][n];
  }
  
  for(m=0; (m<DIM); m++)
    for(n=0; (n<DIM); n++)
      vir[m][n]-=0.5*dvir[m][n];
}

static void dump_it_all(FILE *fp,char *title,
		 int natoms,rvec x[],rvec xp[],rvec v[],
		 rvec vold[],rvec f[])
{
#ifdef DEBUG
  fprintf(fp,"%s\n",title);
  pr_rvecs(fp,0,"x",x,natoms);
  pr_rvecs(fp,0,"xp",xp,natoms);
  pr_rvecs(fp,0,"v",v,natoms);
  pr_rvecs(fp,0,"vold",vold,natoms);
  pr_rvecs(fp,0,"f",f,natoms);
#endif
}

void calc_ke_part(bool bFirstStep,bool bSD,int start,int homenr,
		  rvec vold[],rvec v[],rvec vt[],
		  t_grpopts *opts,t_mdatoms *md,t_groups *grps,
		  t_nrnb *nrnb,real lambda,real *dvdlambda)
{
  int          g,d,n,ga,gt;
  rvec         v_corrt;
  real         hsqrt2,fac,hm,vvt,vct;
  t_grp_tcstat *tcstat=grps->tcstat;
  t_grp_acc    *grpstat=grps->grpstat;
  real         dvdl;

  /* group velocities are calculated in update_grps and
   * accumulated in acumulate_groups.
   * Now the partial global and groups ekin.
   */
  for (g=0; (g<opts->ngtc); g++)
    clear_mat(grps->tcstat[g].ekin); 
    
  if (bFirstStep) {
    for(n=start; (n<start+homenr); n++) {
      copy_rvec(v[n],vold[n]);
    }
    for (g=0; (g<opts->ngacc); g++) {
      for(d=0; (d<DIM); d++)
	grps->grpstat[g].ut[d]=grps->grpstat[g].u[d];
    }
  }
  else {
    for (g=0; (g<opts->ngacc); g++) { 
      for(d=0; (d<DIM); d++)
	grps->grpstat[g].ut[d]=0.5*(grps->grpstat[g].u[d]+
				    grps->grpstat[g].uold[d]);
    }
  }

  hsqrt2 = 0.5*sqrt(2.0);
  dvdl = 0;
  for(n=start; (n<start+homenr); n++) {  
    ga   = md->cACC[n];
    gt   = md->cTC[n];
    hm   = 0.5*md->massT[n];

    if (bSD)
      /* Scale up the average to compensate for the friction */
      fac = (0.5 - hsqrt2)*sdc[gt].em + hsqrt2;
    else
      fac = 0.5;

    for(d=0; (d<DIM); d++) {
      vvt        = fac*(v[n][d]+vold[n][d]);
      vt[n][d]   = vvt;
      vct        = vvt - grpstat[ga].ut[d];
      v_corrt[d] = vct;
    }
    for(d=0; (d<DIM); d++) {
      tcstat[gt].ekin[XX][d]+=hm*v_corrt[XX]*v_corrt[d];
      tcstat[gt].ekin[YY][d]+=hm*v_corrt[YY]*v_corrt[d];
      tcstat[gt].ekin[ZZ][d]+=hm*v_corrt[ZZ]*v_corrt[d];
    }
    if (dvdlambda!=NULL && md->bPerturbed[n]) {
      dvdl-=0.5*(md->massB[n]-md->massA[n])*iprod(v_corrt,v_corrt);
    }
  }
  if(dvdlambda!=NULL)
    *dvdlambda += dvdl;
  
#ifdef DEBUG
  fprintf(stdlog,"ekin: U=(%12e,%12e,%12e)\n",
	  grpstat[0].ut[XX],grpstat[0].ut[YY],grpstat[0].ut[ZZ]);
  fprintf(stdlog,"ekin: %12e\n",trace(tcstat[0].ekin));
#endif

  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

void calc_ke_part_visc(bool bFirstStep,int start,int homenr,
		       matrix box,rvec x[],
		       rvec vold[],rvec v[],rvec vt[],
		       t_grpopts *opts,t_mdatoms *md,t_groups *grps,
		       t_nrnb *nrnb,real lambda,real *dvdlambda)
{
  int          g,d,n,gt;
  rvec         v_corrt;
  real         hm,vvt;
  t_grp_tcstat *tcstat=grps->tcstat;
  t_cos_acc    *cosacc=&(grps->cosacc);
  real         dvdl;
  real         fac,cosz;
  double       mvcos;

  for (g=0; g<opts->ngtc; g++)
    clear_mat(grps->tcstat[g].ekin); 
    
  if (bFirstStep)
    for(n=start; n<start+homenr; n++)
      copy_rvec(v[n],vold[n]);

  fac = 2*M_PI/box[ZZ][ZZ];
  mvcos = 0;
  dvdl = 0;
  for(n=start; n<start+homenr; n++) {  
    gt   = md->cTC[n];
    hm   = 0.5*md->massT[n];
    
    for(d=0; d<DIM; d++) {
      vvt        = 0.5*(v[n][d]+vold[n][d]);
      vt[n][d]   = vvt;
      v_corrt[d] = vvt;
    }
    cosz         = cos(fac*x[n][ZZ]);
    /* Subtract the profile for the kinetic energy */
    v_corrt[XX] -= cosz*cosacc->vcos;
    /* Calculate the amplitude of the new velocity profile */
    mvcos       += 2*cosz*md->massT[n]*v[n][XX];

    for(d=0; d<DIM; d++) {
      tcstat[gt].ekin[XX][d]+=hm*v_corrt[XX]*v_corrt[d];
      tcstat[gt].ekin[YY][d]+=hm*v_corrt[YY]*v_corrt[d];
      tcstat[gt].ekin[ZZ][d]+=hm*v_corrt[ZZ]*v_corrt[d];
    }
    if (dvdlambda!=NULL && md->bPerturbed[n]) {
      dvdl-=0.5*(md->massB[n]-md->massA[n])*iprod(v_corrt,v_corrt);
    }
  }
  if(dvdlambda!=NULL)
    *dvdlambda += dvdl;
  cosacc->mvcos = mvcos;

  inc_nrnb(nrnb,eNR_EKIN,homenr);
}



void update(int          natoms, 	/* number of atoms in simulation */
	    int      	 start,
	    int          homenr,	/* number of home particles 	*/
	    int          step,
	    real         lambda,
	    real         *dvdlambda,    /* FEP stuff */
	    t_parm       *parm,         /* input record and box stuff	*/
	    real         SAfactor,      /* simulated annealing factor   */
	    t_mdatoms    *md,
	    rvec         x[],	        /* coordinates of home particles */
	    t_graph      *graph,	
	    rvec         force[], 	/* forces on home particles 	*/
	    rvec         delta_f[],
	    rvec         vold[],	/* Old velocities		   */
	    rvec         vt[], 		/* velocities at whole timestep  */
	    rvec         v[],  	        /* velocity at next halfstep   	*/
	    t_topology   *top,
	    t_groups     *grps,
	    tensor       vir_part,
	    t_commrec    *cr,
	    t_nrnb       *nrnb,
	    bool         bTYZ,
	    bool         bDoUpdate,
	    t_edsamyn    *edyn,
	    t_pull       *pulldata,
	    bool         bConstrain,
	    bool         bNEMD)
{
  static bool      bFirst=TRUE;
  static rvec      *xprime,*x_unc=NULL;
  static int       ngtc,ngacc;
  static rvec      *lamb;
  static t_edpar   edpar;
  static bool      bHaveConstr,bExtended;
  double           dt;
  real             dt_1,dt_2,mdt_2;
  int              i,n,m,g,vol;
  matrix           M;
  t_inputrec       *ir=&(parm->ir);
  
  if (bFirst) {
    bHaveConstr = init_constraints(stdlog,top,&(parm->ir),md,start,homenr,
				    ir->eI!=eiSteep && bConstrain);
    bHaveConstr = bHaveConstr || pulldata->bPull;
    bExtended   = (ir->etc==etcNOSEHOOVER) || (ir->epc==epcPARINELLORAHMAN);
    
    if (edyn->bEdsam) 
      init_edsam(stdlog,top,md,start,homenr,x,parm->box,
		 edyn,&edpar);
    
    /* Allocate memory for xold, original atomic positions
     * and for xprime.
     */
    snew(xprime,natoms);
    snew(x_unc,homenr);
    /* Copy the pointer to the external acceleration in the opts */
    ngacc=ir->opts.ngacc;    
    ngtc=ir->opts.ngtc;    
       
    snew(lamb,ngtc);

    /* done with initializing */
    bFirst=FALSE;
  }
  
  dt   = ir->delta_t;
  dt_1 = 1.0/dt;
  dt_2 = 1.0/(dt*dt);
  vol  = det(parm->box);

  for(i=0; i<ngtc; i++) {
    real l=grps->tcstat[i].lambda;
    
    if (bTYZ)
      lamb[i][XX]=1;
    else
      lamb[i][XX]=l;
    lamb[i][YY]=l;
    lamb[i][ZZ]=l;
  }
  if (bDoUpdate) {
    /* update mean velocities */
    for (g=0; g<ngacc; g++) {
      copy_rvec(grps->grpstat[g].u,grps->grpstat[g].uold);
      clear_rvec(grps->grpstat[g].u);
    }
    clear_mat(M);

    if(ir->epc == epcPARINELLORAHMAN)
      parinellorahman_pcoupl(&(parm->ir),step,parm->pres,parm->box,parm->boxv,M);
    /* Now do the actual update of velocities and positions */
    where();
    dump_it_all(stdlog,"Before update",natoms,x,xprime,v,vold,force);
    if (ir->eI==eiMD) {
      if (grps->cosacc.cos_accel == 0)
	/* use normal version of update */
	do_update_md(start,homenr,dt,lamb,grps->grpstat,grps->tcstat,
		     ir->opts.acc,ir->opts.nFreeze,md->invmass,md->ptype,
		     md->cFREEZE,md->cACC,md->cTC,x,xprime,v,vold,force,M,
		     bExtended);
      else
	do_update_visc(start,homenr,dt,lamb,md->invmass,grps->tcstat,
		       md->ptype,md->cTC,x,xprime,v,vold,force,M,
		       parm->box,grps->cosacc.cos_accel,grps->cosacc.vcos,bExtended);
    } else if (ir->eI==eiSD) {
      /* The SD update is done in 2 parts, because an extra constraint step
       * is needed 
       */
      do_update_sd(start,homenr,
		   ir->opts.acc,ir->opts.nFreeze,
		   md->invmass,md->ptype,
		   md->cFREEZE,md->cACC,md->cTC,SAfactor,
		   x,xprime,v,vold,force,
		   ir->opts.ngtc,ir->opts.tau_t,ir->opts.ref_t,
		   &ir->ld_seed,TRUE);
      if (bHaveConstr && bConstrain) {
	for(n=start; n<start+homenr; n++)
	  copy_rvec(xprime[n],x_unc[n-start]);
	/* Constrain the coordinates xprime */
	constrain(stdlog,top,ir,step,md,start,homenr,x,xprime,NULL,
		  parm->box,lambda,dvdlambda,nrnb,TRUE);

	for(n=start; n<start+homenr; n++) {
	  /* A correction factor eph is needed for the SD constraint force */
	  mdt_2 = dt_2*md->massT[n]*sdc[md->cTC[n]].eph;
	  for(i=0; i<DIM; i++)
	    delta_f[n][i] = (xprime[n][i] - x_unc[n-start][i])*mdt_2;
	}
      }
      do_update_sd(start,homenr,
		   ir->opts.acc,ir->opts.nFreeze,
		   md->invmass,md->ptype,
		   md->cFREEZE,md->cACC,md->cTC,SAfactor,
		   x,xprime,v,vold,force,
		   ir->opts.ngtc,ir->opts.tau_t,ir->opts.ref_t,
		   &ir->ld_seed,FALSE);
    } else if (ir->eI==eiBD) 
      do_update_bd(start,homenr,dt,
		   ir->opts.nFreeze,md->invmass,md->ptype,
		   md->cFREEZE,md->cTC,
		   x,xprime,v,vold,force,
		   ir->bd_temp,ir->bd_fric,
		   ir->opts.ngtc,ir->opts.tau_t,ir->opts.ref_t,
		   &ir->ld_seed);
    else
      fatal_error(0,"Don't know how to update coordinates");
    
    where();
    inc_nrnb(nrnb, bExtended ? eNR_EXTUPDATE : eNR_UPDATE,homenr);
    dump_it_all(stdlog,"After update",natoms,x,xprime,v,vold,force);
  }
  else {
    /* If we're not updating we're doing shakefirst!
     * In this case the extra coordinates are passed in v array
     */
    for(n=start; (n<start+homenr); n++) {
      copy_rvec(v[n],xprime[n]);
    }
  }
  
  /* 
   *  Steps (7C, 8C)
   *  APPLY CONSTRAINTS:
   *  BLOCK SHAKE 
   */
  
  /* When doing PR pressure coupling we have to constrain the
   * bonds in each iteration. If we are only using Nose-Hoover tcoupling
   * it is enough to do this once though, since the relative velocities 
   * after this will be normal to the bond vector
   */
  if (bHaveConstr && bConstrain) {
    if (ir->eI != eiSD)
      /* Copy Unconstrained X to temp array */
      for(n=start; n<start+homenr; n++)
	copy_rvec(xprime[n],x_unc[n-start]);
    
    /* Constrain the coordinates xprime */
    constrain(stdlog,top,ir,step,md,start,homenr,x,xprime,NULL,
	      parm->box,lambda,dvdlambda,nrnb,TRUE);
    
    where();
    
    dump_it_all(stdlog,"After Shake",natoms,x,xprime,v,vold,force);
    
    /* apply Essential Dynamics constraints when required */
    if (edyn->bEdsam)
      do_edsam(stdlog,top,ir,step,md,start,homenr,xprime,x,
	       x_unc,force,parm->box,edyn,&edpar,bDoUpdate);
    
    /* apply pull constraints when required. Act on xprime, the SHAKED
       coordinates.  Don't do anything to f */
    if (pulldata->bPull && pulldata->runtype != eAfm && 
	pulldata->runtype != eUmbrella &&
	pulldata->runtype != eTest) 
      pull(pulldata,xprime,force,parm->box,top,dt,step,homenr,md); 
    
    where();      
    
    if (bDoUpdate) {
      if (ir->eI != eiSD) {
	/* The constraint virial and the velocities are incorrect for BD */
	for(n=start; n<start+homenr; n++) {
	  mdt_2 = dt_2*md->massT[n];
	  for(i=0; i<DIM; i++) {
	    delta_f[n][i] = (xprime[n][i] - x_unc[n-start][i])*mdt_2;
	    /* recalculate the velocities from the old and new positions */
	    v[n][i]       = (xprime[n][i] - x[n][i])*dt_1;
	  }
	}
	where();
      }
      
      inc_nrnb(nrnb,eNR_SHAKE_V,homenr);
      dump_it_all(stdlog,"After Shake-V",natoms,x,xprime,v,vold,force);
      where();
      
      /* Calculate virial due to constraints (for this node) */
      calc_vir(stdlog,homenr,&(x[start]),&(delta_f[start]),vir_part,cr);
      inc_nrnb(nrnb,eNR_SHAKE_VIR,homenr);
      where();
    }  
  }
  
  /* We must always unshift here, also if we did not shake
   * x was shifted in do_force */
  where();
  if ((graph->nnodes > 0) && bDoUpdate && (ir->ePBC != epbcNONE)) {
    unshift_x(graph,parm->box,x,xprime);
    if(TRICLINIC(parm->box))
      inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
    else
      inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);	  
    for(n=start; (n<graph->start); n++)
      copy_rvec(xprime[n],x[n]);
    for(n=graph->start+graph->nnodes; (n<start+homenr); n++)
      copy_rvec(xprime[n],x[n]);
  }
  else {
    for(n=start; (n<start+homenr); n++)
      copy_rvec(xprime[n],x[n]);
  }
  dump_it_all(stdlog,"After unshift",natoms,x,xprime,v,vold,force);
  where();
  
  if (bDoUpdate) {
    update_grps(start,homenr,grps,&(ir->opts),v,md,bNEMD);
    if (ir->epc == epcBERENDSEN)
      berendsen_pcoupl(ir,step,parm->pres,parm->box,start,homenr,x,md->cFREEZE,nrnb,
		       ir->opts.nFreeze);
    else if (ir->epc == epcPARINELLORAHMAN) {
      /* The box velocities were updated in do_pr_pcoupl in the update
       * iteration, but we dont change the box vectors until we get here
       * since we need to be able to shift/unshift above.
       */
      for(i=0;i<DIM;i++)
	for(m=0;m<=i;m++)
	  parm->box[i][m] += dt * parm->boxv[i][m];
    }
    if (ir->epc != epcNO)
      correct_box(parm->box);
    where();
    /* (un)shifting should NOT be done after this,
     * since the box vectors might have changed
     */
  }
}

  
void correct_ekin(FILE *log,int start,int end,rvec v[],rvec vcm,real mass[],
		  real tmass,tensor ekin)
{
  /* 
   * This is a debugging routine. It should not be called for production code
   *
   * The kinetic energy should calculated according to:
   *   Ekin = 1/2 m (v-vcm)^2
   * However the correction is not always applied, since vcm may not be
   * known in time and we compute
   *   Ekin' = 1/2 m v^2 instead
   * This can be corrected afterwards by computing
   *   Ekin = Ekin' + 1/2 m ( -2 v vcm + vcm^2)
   * or in hsorthand:
   *   Ekin = Ekin' - m v vcm + 1/2 m vcm^2
   */
  int    i,j,k;
  real   m,tm;
  rvec   hvcm,mv;
  tensor dekin;

  /* Local particles */  
  clear_rvec(mv);
  
  /* Processor dependent part. */
  tm = 0;
  for(i=start; (i<end); i++) {
    m      = mass[i];
    tm    += m;
    for(j=0; (j<DIM); j++)
      mv[j] += m*v[i][j];
  }
  /* Shortcut */ 
  svmul(1/tmass,vcm,vcm); 
  svmul(0.5,vcm,hvcm);
  clear_mat(dekin);
  for(j=0; (j<DIM); j++)
    for(k=0; (k<DIM); k++)
      dekin[j][k] += vcm[k]*(tm*hvcm[j]-mv[j]);

  pr_rvecs(log,0,"dekin",dekin,DIM);
  pr_rvecs(log,0," ekin", ekin,DIM);
  fprintf(log,"dekin = %g, ekin = %g  vcm = (%8.4f %8.4f %8.4f)\n",
	  trace(dekin),trace(ekin),vcm[XX],vcm[YY],vcm[ZZ]);
  fprintf(log,"mv = (%8.4f %8.4f %8.4f)\n",
	  mv[XX],mv[YY],mv[ZZ]);
}
