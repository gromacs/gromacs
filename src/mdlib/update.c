/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <math.h>

#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "nrnb.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "main.h"
#include "confio.h"
#include "update.h"
#include "gmx_random.h"
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
#include "disre.h"
#include "orires.h"
#include "gmx_wallcycle.h"

typedef struct {
  double gdt;
  double eph;
  double emh;
  double em;
  double b;
  double c;
  double d;
} gmx_sd_const_t;

typedef struct {
  real V;
  real X;
  real Yv;
  real Yx;
} gmx_sd_sigma_t;

typedef struct gmx_stochd {
  /* The random state */
  gmx_rng_t gaussrand;
  /* BD stuff */
  real *bd_rf;
  /* SD stuff */
  gmx_sd_const_t *sdc;
  gmx_sd_sigma_t *sdsig;
  rvec *sd_V;
  int  sd_V_nalloc;
} t_gmx_stochd;

static void do_update_md(int start,int homenr,double dt,
                         t_grp_tcstat *tcstat,t_grp_acc *gstat,real nh_xi[],
                         rvec accel[],ivec nFreeze[],real invmass[],
                         unsigned short ptype[],unsigned short cFREEZE[],
                         unsigned short cACC[],unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],
                         rvec f[],matrix M,bool bExtended)
{
  double imass,w_dt;
  int    gf=0,ga=0,gt=0;
  rvec   vrel;
  real   vn,vv,va,vb,vnrel;
  real   lg,xi,u;
  int    n,d;

  if(bExtended) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parrinello-Rahman coupling
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */
    for(n=start; n<start+homenr; n++) {
      imass = invmass[n];
      if (cFREEZE)
	gf   = cFREEZE[n];
      if (cACC)
	ga   = cACC[n];
      if (cTC)
	gt   = cTC[n];
      lg   = tcstat[gt].lambda;
      xi   = nh_xi[gt];

      rvec_sub(v[n],gstat[ga].u,vrel);

      for(d=0; d<DIM; d++) {
        if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
          vnrel = lg*(vrel[d] + dt*(imass*f[n][d] - 0.5*xi*vrel[d]
				    - iprod(M[d],vrel)))/(1 + 0.5*xi*dt);  
          /* do not scale the mean velocities u */
          vn             = gstat[ga].u[d] + accel[ga][d]*dt + vnrel; 
          v[n][d]        = vn;
          xprime[n][d]   = x[n][d]+vn*dt;
        } else
          xprime[n][d]   = x[n][d];
      }
    }

  } else {
    /* Classic version of update, used with berendsen coupling */
    for(n=start; n<start+homenr; n++) {
      w_dt = invmass[n]*dt;
      if (cFREEZE)
	gf   = cFREEZE[n];
      if (cACC)
	ga   = cACC[n];
      if (cTC)
	gt   = cTC[n];
      lg   = tcstat[gt].lambda;

      for(d=0; d<DIM; d++) {
        vn             = v[n][d];

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
          vv             = lg*(vn + f[n][d]*w_dt);

          /* do not scale the mean velocities u */
          u              = gstat[ga].u[d];
          va             = vv + accel[ga][d]*dt;
          vb             = va + (1.0-lg)*u;
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
			   t_grp_tcstat *tcstat,real invmass[],real nh_xi[],
                           unsigned short ptype[],unsigned short cTC[],
                           rvec x[],rvec xprime[],rvec v[],
                           rvec f[],matrix M,matrix box,real
                           cos_accel,real vcos,bool bExtended)
{
  double imass,w_dt;
  int    gt=0;
  real   vn,vc;
  real   lg,xi,vv;
  real   fac,cosz;
  rvec   vrel;
  int    n,d;

  fac = 2*M_PI/(box[ZZ][ZZ]);

  if(bExtended) {
    /* Update with coupling to extended ensembles, used for
     * Nose-Hoover and Parrinello-Rahman coupling
     */
    for(n=start; n<start+homenr; n++) {
      imass = invmass[n];
      if (cTC)
	gt   = cTC[n];
      lg   = tcstat[gt].lambda;
      cosz = cos(fac*x[n][ZZ]);

      copy_rvec(v[n],vrel);

      vc            = cosz*vcos;
      vrel[XX]     -= vc;
      xi           = nh_xi[gt];

      for(d=0; d<DIM; d++) {
        vn             = v[n][d];

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
          vn    = lg*(vrel[d] + dt*(imass*f[n][d] - 0.5*xi*vrel[d]
				    - iprod(M[d],vrel)))/(1 + 0.5*xi*dt);
          if(d == XX)
            vn += vc + dt*cosz*cos_accel;

          v[n][d]        = vn;
          xprime[n][d]   = x[n][d]+vn*dt;
        } else
          xprime[n][d]   = x[n][d];
      }
    }

  } else {
    /* Classic version of update, used with berendsen coupling */
    for(n=start; n<start+homenr; n++) {
      w_dt = invmass[n]*dt;
      if (cTC)
        gt   = cTC[n];
      lg   = tcstat[gt].lambda;
      cosz = cos(fac*x[n][ZZ]);

      for(d=0; d<DIM; d++) {
        vn             = v[n][d];

        if((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
          if(d == XX) {
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

gmx_stochd_t init_stochd(FILE *fplog,t_inputrec *ir)
{
  t_gmx_stochd *sd;
  gmx_sd_const_t *sdc;
  int  ngtc,n;
  real y;

  snew(sd,1);

  /* Initiate random number generator for langevin type dynamics,
   * for BD, SD or velocity rescaling temperature coupling.
   */
  sd->gaussrand = gmx_rng_init(ir->ld_seed);

  ngtc = ir->opts.ngtc;

  if (ir->eI == eiBD) {
    snew(sd->bd_rf,ngtc);
  } else if (EI_SD(ir->eI)) {
    snew(sd->sdc,ngtc);
    snew(sd->sdsig,ngtc);
    
    sdc = sd->sdc;
    for(n=0; n<ngtc; n++) {
      sdc[n].gdt = ir->delta_t/ir->opts.tau_t[n];
      sdc[n].eph = exp(sdc[n].gdt/2);
      sdc[n].emh = exp(-sdc[n].gdt/2);
      sdc[n].em  = exp(-sdc[n].gdt);
      if (sdc[n].gdt >= 0.05) {
	sdc[n].b = sdc[n].gdt*(sdc[n].eph*sdc[n].eph - 1) 
	  - 4*(sdc[n].eph - 1)*(sdc[n].eph - 1);
	sdc[n].c = sdc[n].gdt - 3 + 4*sdc[n].emh - sdc[n].em;
	sdc[n].d = 2 - sdc[n].eph - sdc[n].emh;
      } else {
	y = sdc[n].gdt/2;
	/* Seventh order expansions for small y */
	sdc[n].b = y*y*y*y*(1/3.0+y*(1/3.0+y*(17/90.0+y*7/9.0)));
	sdc[n].c = y*y*y*(2/3.0+y*(-1/2.0+y*(7/30.0+y*(-1/12.0+y*31/1260.0))));
	sdc[n].d = y*y*(-1+y*y*(-1/12.0-y*y/360.0));
      }
      if(debug)
	fprintf(debug,"SD const tc-grp %d: b %g  c %g  d %g\n",
		n,sdc[n].b,sdc[n].c,sdc[n].d);
    }
  }

  return sd;
}

void get_stochd_state(t_gmx_stochd *sd,t_state *state)
{
  gmx_rng_get_state(sd->gaussrand,state->ld_rng,state->ld_rngi);
}

void set_stochd_state(t_gmx_stochd *sd,t_state *state)
{
  gmx_rng_set_state(sd->gaussrand,state->ld_rng,state->ld_rngi[0]);
}

static void do_update_sd1(t_gmx_stochd *sd,
			  int start,int homenr,double dt,
			  rvec accel[],ivec nFreeze[],
			  real invmass[],unsigned short ptype[],
			  unsigned short cFREEZE[],unsigned short cACC[],
			  unsigned short cTC[],
			  rvec x[],rvec xprime[],rvec v[],rvec f[],
			  rvec sd_X[],
			  int ngtc,real tau_t[],real ref_t[])
{
  gmx_sd_const_t *sdc;
  gmx_sd_sigma_t *sig;
  gmx_rng_t gaussrand;
  real   kT;
  int    gf=0,ga=0,gt=0;
  real   ism,sd_V;
  int    n,d;

  sdc = sd->sdc;
  sig = sd->sdsig;
  if (homenr > sd->sd_V_nalloc) {
    sd->sd_V_nalloc = over_alloc_dd(homenr);
    srenew(sd->sd_V,sd->sd_V_nalloc);
  }
  gaussrand = sd->gaussrand;
  
  for(n=0; n<ngtc; n++) {
    kT = BOLTZ*ref_t[n];
    /* The mass is encounted for later, since this differs per atom */
    sig[n].V  = sqrt(2*kT*(1 - sdc[n].em));
  }

  for(n=start; n<start+homenr; n++) {
    ism = sqrt(invmass[n]);
    if (cFREEZE)
      gf  = cFREEZE[n];
    if (cACC)
      ga  = cACC[n];
    if (cTC)
      gt  = cTC[n];

    for(d=0; d<DIM; d++) {
      if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
	sd_V = ism*sig[gt].V*gmx_rng_gaussian_table(gaussrand);
	
	v[n][d] = v[n][d]*sdc[gt].em 
	  + (invmass[n]*f[n][d] + accel[ga][d])*tau_t[gt]*(1 - sdc[gt].em)
	  + sd_V;

	xprime[n][d] = x[n][d] + v[n][d]*dt;
      } else {
	v[n][d]      = 0.0;
	xprime[n][d] = x[n][d];
      }
    }
  }
}

static void do_update_sd2(t_gmx_stochd *sd,bool bInitStep,
			  int start,int homenr,
			  rvec accel[],ivec nFreeze[],
			  real invmass[],unsigned short ptype[],
			  unsigned short cFREEZE[],unsigned short cACC[],
			  unsigned short cTC[],
			  rvec x[],rvec xprime[],rvec v[],rvec f[],
			  rvec sd_X[],
			  int ngtc,real tau_t[],real ref_t[],
			  bool bFirstHalf)
{
  gmx_sd_const_t *sdc;
  gmx_sd_sigma_t *sig;
  /* The random part of the velocity update, generated in the first
   * half of the update, needs to be remembered for the second half.
   */
  rvec *sd_V;
  gmx_rng_t gaussrand;
  real   kT;
  int    gf=0,ga=0,gt=0;
  real   vn=0,Vmh,Xmh;
  real   ism;
  int    n,d;

  sdc = sd->sdc;
  sig = sd->sdsig;
  if (homenr > sd->sd_V_nalloc) {
    sd->sd_V_nalloc = over_alloc_dd(homenr);
    srenew(sd->sd_V,sd->sd_V_nalloc);
  }
  sd_V = sd->sd_V;
  gaussrand = sd->gaussrand;

  if(bFirstHalf) {
    for(n=0; n<ngtc; n++) {
      kT = BOLTZ*ref_t[n];
      /* The mass is encounted for later, since this differs per atom */
      sig[n].V  = sqrt(kT*(1-sdc[n].em));
      sig[n].X  = sqrt(kT*sqr(tau_t[n])*sdc[n].c);
      sig[n].Yv = sqrt(kT*sdc[n].b/sdc[n].c);
      sig[n].Yx = sqrt(kT*sqr(tau_t[n])*sdc[n].b/(1-sdc[n].em));
    }
  }

  for(n=start; n<start+homenr; n++) {
    ism = sqrt(invmass[n]);
    if (cFREEZE)
      gf  = cFREEZE[n];
    if (cACC)
      ga  = cACC[n];
    if (cTC)
      gt  = cTC[n];

    for(d=0; d<DIM; d++) {
      if(bFirstHalf) {
        vn             = v[n][d];
      }
      if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) {
        if (bFirstHalf) {

          if (bInitStep)
            sd_X[n][d] = ism*sig[gt].X*gmx_rng_gaussian_table(gaussrand);

          Vmh = sd_X[n][d]*sdc[gt].d/(tau_t[gt]*sdc[gt].c) 
                + ism*sig[gt].Yv*gmx_rng_gaussian_table(gaussrand);
          sd_V[n-start][d] = ism*sig[gt].V*gmx_rng_gaussian_table(gaussrand);

          v[n][d] = vn*sdc[gt].em 
                    + (invmass[n]*f[n][d] + accel[ga][d])*tau_t[gt]*(1 - sdc[gt].em)
                    + sd_V[n-start][d] - sdc[gt].em*Vmh;

          xprime[n][d] = x[n][d] + v[n][d]*tau_t[gt]*(sdc[gt].eph - sdc[gt].emh); 

        } else {

          /* Correct the velocities for the constraints.
	   * This operation introduces some inaccuracy,
	   * since the velocity is determined from differences in coordinates.
	   */
          v[n][d] = 
          (xprime[n][d] - x[n][d])/(tau_t[gt]*(sdc[gt].eph - sdc[gt].emh));  

          Xmh = sd_V[n-start][d]*tau_t[gt]*sdc[gt].d/(sdc[gt].em-1) 
                + ism*sig[gt].Yx*gmx_rng_gaussian_table(gaussrand);
          sd_X[n][d] = ism*sig[gt].X*gmx_rng_gaussian_table(gaussrand);

          xprime[n][d] += sd_X[n][d] - Xmh;

        }
      } else {
        if(bFirstHalf) {
          v[n][d]        = 0.0;
          xprime[n][d]   = x[n][d];
        }
      }
    }
  }
}

static void do_update_bd(int start,int homenr,double dt,
                         ivec nFreeze[],
                         real invmass[],unsigned short ptype[],
                         unsigned short cFREEZE[],unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],
                         rvec f[],real friction_coefficient,
                         int ngtc,real tau_t[],real ref_t[],
			 real *rf,gmx_rng_t gaussrand)
{
  int    gf=0,gt=0;
  real   vn;
  real   invfr=0;
  int    n,d;

  if (friction_coefficient != 0) {
    invfr = 1.0/friction_coefficient;
    for(n=0; n<ngtc; n++)
      rf[n] = sqrt(2.0*BOLTZ*ref_t[n]/(friction_coefficient*dt));
  } else
    for(n=0; n<ngtc; n++)
      rf[n] = sqrt(2.0*BOLTZ*ref_t[n]);

  for(n=start; (n<start+homenr); n++) {
    if (cFREEZE)
      gf = cFREEZE[n];
    if (cTC)
      gt = cTC[n];
    for(d=0; (d<DIM); d++) {
      if((ptype[n]!=eptVSite) && (ptype[n]!=eptShell) && !nFreeze[gf][d]) {
        if (friction_coefficient != 0)
          vn = invfr*f[n][d] + rf[gt]*gmx_rng_gaussian_table(gaussrand);
        else
          /* NOTE: invmass = 1/(mass*friction_constant*dt) */
          vn = invmass[n]*f[n][d]*dt 
	    + sqrt(invmass[n])*rf[gt]*gmx_rng_gaussian_table(gaussrand);

        v[n][d]      = vn;
        xprime[n][d] = x[n][d]+vn*dt;
      } else {
        v[n][d]      = 0.0;
        xprime[n][d] = x[n][d];
      }
    }
  }
}

static void dump_it_all(FILE *fp,char *title,
                        int natoms,rvec x[],rvec xp[],rvec v[],rvec f[])
{
#ifdef DEBUG
  if (fp) {
    fprintf(fp,"%s\n",title);
    pr_rvecs(fp,0,"x",x,natoms);
    pr_rvecs(fp,0,"xp",xp,natoms);
    pr_rvecs(fp,0,"v",v,natoms);
    pr_rvecs(fp,0,"f",f,natoms);
  }
#endif
}

void calc_ke_part(rvec v[],t_grpopts *opts,t_mdatoms *md,t_groups *grps,
                  t_nrnb *nrnb,real lambda)
{
  int          start=md->start,homenr=md->homenr;
  int          g,d,n,ga=0,gt=0;
  rvec         v_corrt;
  real         hm;
  t_grp_tcstat *tcstat=grps->tcstat;
  t_grp_acc    *grpstat=grps->grpstat;
  real         dekindl;

  /* group velocities are calculated in update_grps and
   * accumulated in acumulate_groups.
   * Now the partial global and groups ekin.
   */
  for(g=0; (g<opts->ngtc); g++) {
    copy_mat(grps->tcstat[g].ekinh,grps->tcstat[g].ekinh_old);
    clear_mat(grps->tcstat[g].ekinh);
  }
  grps->dekindl_old = grps->dekindl;

  dekindl = 0;
  for(n=start; (n<start+homenr); n++) {
    if (md->cACC)
      ga = md->cACC[n];
    if (md->cTC)
      gt = md->cTC[n];
    hm   = 0.5*md->massT[n];

    for(d=0; (d<DIM); d++) {
      v_corrt[d] = v[n][d] - grpstat[ga].u[d];
    }
    for(d=0; (d<DIM); d++) {
      tcstat[gt].ekinh[XX][d]+=hm*v_corrt[XX]*v_corrt[d];
      tcstat[gt].ekinh[YY][d]+=hm*v_corrt[YY]*v_corrt[d];
      tcstat[gt].ekinh[ZZ][d]+=hm*v_corrt[ZZ]*v_corrt[d];
    }
    if (md->nMassPerturbed && md->bPerturbed[n])
      dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt,v_corrt);
  }
  grps->dekindl = dekindl;

  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

void calc_ke_part_visc(matrix box,rvec x[],rvec v[],
                       t_grpopts *opts,t_mdatoms *md,t_groups *grps,
                       t_nrnb *nrnb,real lambda)
{
  int          start=md->start,homenr=md->homenr;
  int          g,d,n,gt=0;
  rvec         v_corrt;
  real         hm;
  t_grp_tcstat *tcstat=grps->tcstat;
  t_cos_acc    *cosacc=&(grps->cosacc);
  real         dekindl;
  real         fac,cosz;
  double       mvcos;

  for(g=0; g<opts->ngtc; g++) {
    copy_mat(grps->tcstat[g].ekinh,grps->tcstat[g].ekinh_old);
    clear_mat(grps->tcstat[g].ekinh);
  }
  grps->dekindl_old = grps->dekindl;

  fac = 2*M_PI/box[ZZ][ZZ];
  mvcos = 0;
  dekindl = 0;
  for(n=start; n<start+homenr; n++) {
    if (md->cTC)
      gt = md->cTC[n];
    hm   = 0.5*md->massT[n];

    /* Note that the times of x and v differ by half a step */
    cosz         = cos(fac*x[n][ZZ]);
    /* Calculate the amplitude of the new velocity profile */
    mvcos       += 2*cosz*md->massT[n]*v[n][XX];

    copy_rvec(v[n],v_corrt);
    /* Subtract the profile for the kinetic energy */
    v_corrt[XX] -= cosz*cosacc->vcos;
    for(d=0; (d<DIM); d++) {
      tcstat[gt].ekinh[XX][d]+=hm*v_corrt[XX]*v_corrt[d];
      tcstat[gt].ekinh[YY][d]+=hm*v_corrt[YY]*v_corrt[d];
      tcstat[gt].ekinh[ZZ][d]+=hm*v_corrt[ZZ]*v_corrt[d];
    }
    if(md->nPerturbed && md->bPerturbed[n])
      dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt,v_corrt);
  }
  grps->dekindl = dekindl;
  cosacc->mvcos = mvcos;

  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

/* Static variables for the deform algorithm */
static int    deformref_step;
static matrix deformref_box;

void set_deform_reference_box(int step,matrix box)
{
  deformref_step = step;
  copy_mat(box,deformref_box);
}

static void deform(int start,int homenr,rvec x[],matrix box,matrix *scale_tot,
		   const t_inputrec *ir,int step)
{
  matrix new,invbox,mu;
  real   elapsed_time;
  int    i,j;  

  elapsed_time = (step + 1 - deformref_step)*ir->delta_t;
  copy_mat(box,new);
  for(i=0; i<DIM; i++) {
    for(j=0; j<DIM; j++) {
      if (ir->deform[i][j] != 0) {
	new[i][j] = deformref_box[i][j] + elapsed_time*ir->deform[i][j];
      }
    }
  }
  /* We correct the off-diagonal elements,
   * which can grow indefinitely during shearing,
   * so the shifts do not get messed up.
   */
  for(i=1; i<DIM; i++) {
    for(j=i-1; j>=0; j--) {
      while (new[i][j] - box[i][j] > 0.5*new[j][j]) {
	rvec_dec(new[i],new[j]);
      }
      while (new[i][j] - box[i][j] < -0.5*new[j][j]) {
	rvec_inc(new[i],new[j]);
      }
    }
  }
  m_inv_ur0(box,invbox);
  copy_mat(new,box);
  mmul_ur0(box,invbox,mu);
  
  for(i=start; i<start+homenr; i++) {
    x[i][XX] = mu[XX][XX]*x[i][XX]+mu[YY][XX]*x[i][YY]+mu[ZZ][XX]*x[i][ZZ];
    x[i][YY] = mu[YY][YY]*x[i][YY]+mu[ZZ][YY]*x[i][ZZ];
    x[i][ZZ] = mu[ZZ][ZZ]*x[i][ZZ];
  }
  if (*scale_tot) {
    /* The transposes of the scaling matrices are stored,
     * so we need to do matrix multiplication in the inverse order.
     */
    mmul_ur0(*scale_tot,mu,*scale_tot);
  }
}

void update(FILE         *fplog,
	    int          step,
            real         *dvdlambda,    /* FEP stuff */
            t_inputrec   *inputrec,      /* input record and box stuff	*/
            t_mdatoms    *md,
	    t_state      *state,
            t_graph      *graph,  
            rvec         force[],        /* forces on home particles */
	    rvec         xprime[],       /* buffer for x for update  */
	    t_fcdata     *fcd,
            t_topology   *top,
            t_groups     *grps,
            tensor       vir_part,
	    matrix       *scale_tot,
            t_commrec    *cr,
            t_nrnb       *nrnb,
	    gmx_wallcycle_t wcycle,
	    gmx_stochd_t sd,
	    gmx_constr_t constr,
	    bool         bHaveConstr,
            bool         bNEMD,
	    bool         bInitStep)
{
  bool             bExtended,bLastStep,bLog=FALSE,bEner=FALSE;
  double           dt;
  real             dt_1;
  int              start,homenr,i,n,m,g;
  matrix           pcoupl_mu,M;
  tensor           vir_con;
  
  start  = md->start;
  homenr = md->homenr;
  
  /* We need to update the NMR restraint history when time averaging is used */
  if (state->flags & (1<<estDISRE_RM3TAV)) {
    update_disres_history(fcd,&state->hist);
  }
  if (state->flags & (1<<estORIRE_DTAV)) {
    update_orires_history(fcd,&state->hist);
  }

  bExtended =
    (inputrec->etc == etcNOSEHOOVER) ||
    (inputrec->epc == epcPARRINELLORAHMAN);

  dt   = inputrec->delta_t;
  dt_1 = 1.0/dt;

  clear_mat(pcoupl_mu);
  for(i=0; i<DIM; i++) {
    pcoupl_mu[i][i] = 1.0;
  }
  clear_mat(M);

  /* We can always tcoupl, even if we did not sum the energies
   * the previous step, since grps->tcstat[i].Th is only updated
   * when the energies have been summed.
   */
  switch (inputrec->etc) {
  case etcNO:
    break;
  case etcBERENDSEN:
    berendsen_tcoupl(&(inputrec->opts),grps,inputrec->delta_t);
    break;
  case etcNOSEHOOVER:
    nosehoover_tcoupl(&(inputrec->opts),grps,inputrec->delta_t,
		      state->nosehoover_xi,state->therm_integral);
    break;
  case etcVRESCALE:
    vrescale_tcoupl(&(inputrec->opts),grps,inputrec->delta_t,
		    state->therm_integral,sd->gaussrand);
    break;
  }
  /* We can always pcoupl, even if we did not sum the energies
   * the previous step, since state->pres_prev is only updated
   * when the energies have been summed.
   */
  if (inputrec->epc == epcBERENDSEN && !bInitStep) {
    berendsen_pcoupl(fplog,step,inputrec,state->pres_prev,state->box,
		     pcoupl_mu);
  }
  if (inputrec->epc == epcPARRINELLORAHMAN) {
    parrinellorahman_pcoupl(fplog,step,inputrec,state->pres_prev,
			    state->box,state->box_rel,state->boxv,
			    M,scale_tot,bInitStep);
  }
  
  /* Now do the actual update of velocities and positions */
  where();
  dump_it_all(fplog,"Before update",
	      state->natoms,state->x,xprime,state->v,force);
  if (inputrec->eI == eiMD) {
    if (grps->cosacc.cos_accel == 0) {
      /* use normal version of update */
      do_update_md(start,homenr,dt,
		   grps->tcstat,grps->grpstat,state->nosehoover_xi,
		   inputrec->opts.acc,inputrec->opts.nFreeze,md->invmass,md->ptype,
		   md->cFREEZE,md->cACC,md->cTC,
		   state->x,xprime,state->v,force,M,
		   bExtended);
    } else {
      do_update_visc(start,homenr,dt,
		     grps->tcstat,md->invmass,state->nosehoover_xi,
		     md->ptype,md->cTC,state->x,xprime,state->v,force,M,
		     state->box,grps->cosacc.cos_accel,grps->cosacc.vcos,bExtended);
    }
  } else if (inputrec->eI == eiSD1) {
    do_update_sd1(sd,start,homenr,dt,
		  inputrec->opts.acc,inputrec->opts.nFreeze,
		  md->invmass,md->ptype,
		  md->cFREEZE,md->cACC,md->cTC,
		  state->x,xprime,state->v,force,state->sd_X,
		  inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t);
  } else if (inputrec->eI == eiSD2) {
    /* The SD update is done in 2 parts, because an extra constraint step
     * is needed 
     */
    do_update_sd2(sd,bInitStep,start,homenr,
		  inputrec->opts.acc,inputrec->opts.nFreeze,
		  md->invmass,md->ptype,
		  md->cFREEZE,md->cACC,md->cTC,
		  state->x,xprime,state->v,force,state->sd_X,
		  inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
		  TRUE);
  } else if (inputrec->eI == eiBD) {
    do_update_bd(start,homenr,dt,
		 inputrec->opts.nFreeze,md->invmass,md->ptype,
		 md->cFREEZE,md->cTC,
		 state->x,xprime,state->v,force,
		 inputrec->bd_fric,
		 inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
		 sd->bd_rf,sd->gaussrand);
  } else {
    gmx_fatal(FARGS,"Don't know how to update coordinates");
  }
  where();
  inc_nrnb(nrnb, bExtended ? eNR_EXTUPDATE : eNR_UPDATE, homenr);
  dump_it_all(fplog,"After update",
	      state->natoms,state->x,xprime,state->v,force);

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
  if (bHaveConstr) {
    bLastStep = (step == inputrec->init_step+inputrec->nsteps);
    bLog  = (do_per_step(step,inputrec->nstlog) || bLastStep || (step < 0));
    bEner = (do_per_step(step,inputrec->nstenergy) || bLastStep);
    if (constr) {
      /* Constrain the coordinates xprime */
      wallcycle_start(wcycle,ewcCONSTR);
      constrain(NULL,bLog,bEner,constr,top,
		inputrec,cr,step,1,md,
		state->x,xprime,NULL,
		state->box,state->lambda,dvdlambda,
		state->v,&vir_con,nrnb,econqCoord);
      wallcycle_stop(wcycle,ewcCONSTR);
    }
    where();

    dump_it_all(fplog,"After Shake",
		state->natoms,state->x,xprime,state->v,force);

    if (inputrec->eI == eiSD2) {
      /* A correction factor eph is needed for the SD constraint force */
      /* Here we can, unfortunately, not have proper corrections
       * for different friction constants, so we use the first one.
       */
      for(i=0; i<DIM; i++)
	for(m=0; m<DIM; m++)
	  vir_part[i][m] += sd->sdc[0].eph*vir_con[i][m];
    } else {
      m_add(vir_part,vir_con,vir_part);
    }
    if (debug)
      pr_rvecs(debug,0,"constraint virial",vir_part,DIM);
    where();
  }
  
  where();
  if (inputrec->eI == eiSD2) {
    /* The second part of the SD integration */
    do_update_sd2(sd,FALSE,start,homenr,
		  inputrec->opts.acc,inputrec->opts.nFreeze,
		  md->invmass,md->ptype,
		  md->cFREEZE,md->cACC,md->cTC,
		  state->x,xprime,state->v,force,state->sd_X,
		  inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
		  FALSE);
    inc_nrnb(nrnb, eNR_UPDATE, homenr);
    
    if (constr) {
      /* Constrain the coordinates xprime */
      wallcycle_start(wcycle,ewcCONSTR);
      constrain(NULL,bLog,bEner,constr,top,
		inputrec,cr,step,1,md,
		state->x,xprime,NULL,
		state->box,state->lambda,dvdlambda,
		NULL,NULL,nrnb,econqCoord);
      wallcycle_stop(wcycle,ewcCONSTR);
    }
  }
  
  /* We must always unshift here, also if we did not shake
   * x was shifted in do_force */
  
  if (graph && (graph->nnodes > 0)) {
    unshift_x(graph,state->box,state->x,xprime);
    if (TRICLINIC(state->box))
      inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
    else
      inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);    
    for(n=start; (n<graph->start); n++)
      copy_rvec(xprime[n],state->x[n]);
    for(n=graph->start+graph->nnodes; (n<start+homenr); n++)
      copy_rvec(xprime[n],state->x[n]);
  } else {
    for(n=start; (n<start+homenr); n++)
      copy_rvec(xprime[n],state->x[n]);
  }
  dump_it_all(fplog,"After unshift",
	      state->natoms,state->x,xprime,state->v,force);
  where();

  update_grps(start,homenr,grps,&(inputrec->opts),state->v,md,state->lambda,
	      bNEMD);

  if (inputrec->epc == epcBERENDSEN) {
    berendsen_pscale(inputrec,pcoupl_mu,state->box,state->box_rel,
		     start,homenr,state->x,md->cFREEZE,nrnb);
    if (scale_tot) {
      /* The transposes of the scaling matrices are stored,
       * therefore we need to reverse the order in the multiplication.
       */
      mmul_ur0(*scale_tot,pcoupl_mu,*scale_tot);
    }
  } else if (inputrec->epc == epcPARRINELLORAHMAN) {
    /* The box velocities were updated in do_pr_pcoupl in the update
     * iteration, but we dont change the box vectors until we get here
     * since we need to be able to shift/unshift above.
     */
    for(i=0;i<DIM;i++)
      for(m=0;m<=i;m++)
	state->box[i][m] += dt*state->boxv[i][m];

    preserve_box_shape(inputrec,state->box_rel,state->box);
  }
  if (DEFORM(*inputrec)) {
    deform(start,homenr,state->x,state->box,scale_tot,inputrec,step);
  }
  where();
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
