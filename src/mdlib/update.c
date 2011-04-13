/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

/*For debugging, start at v(-dt/2) for velolcity verlet -- uncomment next line */
/*#define STARTFROMDT2*/

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

typedef struct {
  /* The random state */
  gmx_rng_t gaussrand;
  /* BD stuff */
  real *bd_rf;
  /* SD stuff */
  gmx_sd_const_t *sdc;
  gmx_sd_sigma_t *sdsig;
  rvec *sd_V;
  int  sd_V_nalloc;
} gmx_stochd_t;

typedef struct gmx_update
{
    gmx_stochd_t *sd;
    rvec *xp;
    int  xp_nalloc;
    /* Variables for the deform algorithm */
    gmx_large_int_t deformref_step;
    matrix     deformref_box;
} t_gmx_update;


void store_rvec(rvec *from, rvec *to, int n) {
  int i;
  for (i=0;i<n;i++) {
    copy_rvec(from[i],to[i]);
  }
}

static void do_update_md(int start,int nrend,double dt,
                         t_grp_tcstat *tcstat,t_grp_acc *gstat,double nh_vxi[],
                         rvec accel[],ivec nFreeze[],real invmass[],
                         unsigned short ptype[],unsigned short cFREEZE[],
                         unsigned short cACC[],unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],
                         rvec f[],matrix M,
                         gmx_bool bNH,gmx_bool bPR)
{
  double imass,w_dt;
  int    gf=0,ga=0,gt=0;
  rvec   vrel;
  real   vn,vv,va,vb,vnrel;
  real   lg,vxi=0,u;
  int    n,d;

  if (bNH || bPR) 
  {
      /* Update with coupling to extended ensembles, used for                     
       * Nose-Hoover and Parrinello-Rahman coupling                               
       * Nose-Hoover uses the reversible leap-frog integrator from                
       * Holian et al. Phys Rev E 52(3) : 2338, 1995                              
       */
      for(n=start; n<nrend; n++) 
      {
          imass = invmass[n];
          if (cFREEZE) 
          {
              gf   = cFREEZE[n];
          }
          if (cACC) 
          {
              ga   = cACC[n];
          }
          if (cTC)
          {
              gt   = cTC[n];
          }
          lg   = tcstat[gt].lambda;
          if (bNH) {
              vxi   = nh_vxi[gt];
          }
          rvec_sub(v[n],gstat[ga].u,vrel);
          
          for(d=0; d<DIM; d++) 
          {
              if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) 
              {
                  vnrel = (lg*vrel[d] + dt*(imass*f[n][d] - 0.5*vxi*vrel[d]
                                            - iprod(M[d],vrel)))/(1 + 0.5*vxi*dt);  
                  /* do not scale the mean velocities u */
                  vn             = gstat[ga].u[d] + accel[ga][d]*dt + vnrel; 
                  v[n][d]        = vn;
                  xprime[n][d]   = x[n][d]+vn*dt;
              } 
              else 
              {
                  v[n][d]        = 0.0;
                  xprime[n][d]   = x[n][d];
              }
          }
      }
  } 
  else 
  {
      /* Classic version of update, used with berendsen coupling */
      for(n=start; n<nrend; n++) 
      {
          w_dt = invmass[n]*dt;
          if (cFREEZE) 
          {
              gf   = cFREEZE[n];
          }
          if (cACC)
          {
              ga   = cACC[n];
          }
          if (cTC) 
          {
              gt   = cTC[n];
          }
          lg   = tcstat[gt].lambda;
          
          for(d=0; d<DIM; d++) 
          {
              vn             = v[n][d];
              if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) 
              {
                  vv             = lg*vn + f[n][d]*w_dt;
                  
                  /* do not scale the mean velocities u */
                  u              = gstat[ga].u[d];
                  va             = vv + accel[ga][d]*dt;
                  vb             = va + (1.0-lg)*u;
                  v[n][d]        = vb;
                  xprime[n][d]   = x[n][d]+vb*dt;
              } 
              else 
              {
                  v[n][d]        = 0.0;
                  xprime[n][d]   = x[n][d];
              }
          }
      }
  }
}

static void do_update_vv_vel(int start,int nrend,double dt,
                             t_grp_tcstat *tcstat,t_grp_acc *gstat,
                             rvec accel[],ivec nFreeze[],real invmass[],
                             unsigned short ptype[],
                             unsigned short cFREEZE[],unsigned short cACC[],
                             rvec v[],rvec f[],
                             gmx_bool bExtended, real veta, real alpha)
{
    double imass,w_dt;
    int    gf=0,ga=0,gt=0;
    rvec   vrel;
    real   u,vn,vv,va,vb,vnrel;
    int    n,d;
    double g,mv1,mv2;
    
    if (bExtended)
    {
        g        = 0.25*dt*veta*alpha;
        mv1      = exp(-g);
        mv2      = series_sinhx(g);
    }
    else 
    {
        mv1      = 1.0;
        mv2      = 1.0;
    }
    for(n=start; n<nrend; n++) 
    {
        w_dt = invmass[n]*dt;
        if (cFREEZE)
        {
            gf   = cFREEZE[n];
        }
        if (cACC)
        {
            ga   = cACC[n];
        }
        
        for(d=0; d<DIM; d++) 
        {
            if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) 
            {
                v[n][d]             = mv1*(mv1*v[n][d] + 0.5*(w_dt*mv2*f[n][d]))+0.5*accel[ga][d]*dt;
            } 
            else 
            {
                v[n][d]        = 0.0;
            }
        }
    }
} /* do_update_vv_vel */

static void do_update_vv_pos(int start,int nrend,double dt,
                             t_grp_tcstat *tcstat,t_grp_acc *gstat,
                             rvec accel[],ivec nFreeze[],real invmass[],
                             unsigned short ptype[],
                             unsigned short cFREEZE[],
                             rvec x[],rvec xprime[],rvec v[],
                             rvec f[],gmx_bool bExtended, real veta, real alpha)
{
  double imass,w_dt;
  int    gf=0;
  int    n,d;
  double g,mr1,mr2;

  if (bExtended) {
      g        = 0.5*dt*veta;
      mr1      = exp(g);
      mr2      = series_sinhx(g);
  }
  else 
  {
      mr1      = 1.0;
      mr2      = 1.0;
  }
  
  for(n=start; n<nrend; n++) {
      if (cFREEZE)
      {
          gf   = cFREEZE[n];
      }
      
      for(d=0; d<DIM; d++) 
      {
          if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) 
          {
              xprime[n][d]   = mr1*(mr1*x[n][d]+mr2*dt*v[n][d]);
          } 
          else 
          {
              xprime[n][d]   = x[n][d];
          }
      }
  }
}/* do_update_vv_pos */

static void do_update_visc(int start,int nrend,double dt,
                           t_grp_tcstat *tcstat,real invmass[],double nh_vxi[],
                           unsigned short ptype[],unsigned short cTC[],
                           rvec x[],rvec xprime[],rvec v[],
                           rvec f[],matrix M,matrix box,real
                           cos_accel,real vcos,
                           gmx_bool bNH,gmx_bool bPR)
{
    double imass,w_dt;
    int    gt=0;
    real   vn,vc;
    real   lg,vxi=0,vv;
    real   fac,cosz;
    rvec   vrel;
    int    n,d;
    
    fac = 2*M_PI/(box[ZZ][ZZ]);
    
    if (bNH || bPR) {
        /* Update with coupling to extended ensembles, used for
         * Nose-Hoover and Parrinello-Rahman coupling
         */
        for(n=start; n<nrend; n++) {
            imass = invmass[n];
            if (cTC) 
            {
                gt   = cTC[n];
            }
            lg   = tcstat[gt].lambda;
            cosz = cos(fac*x[n][ZZ]);
            
            copy_rvec(v[n],vrel);
            
            vc            = cosz*vcos;
            vrel[XX]     -= vc;
            if (bNH) 
            {
                vxi        = nh_vxi[gt];
            }
            for(d=0; d<DIM; d++) 
            {
                vn             = v[n][d];
                
                if((ptype[n] != eptVSite) && (ptype[n] != eptShell)) 
                {
                    vn  = (lg*vrel[d] + dt*(imass*f[n][d] - 0.5*vxi*vrel[d]
                                            - iprod(M[d],vrel)))/(1 + 0.5*vxi*dt);
                    if(d == XX) 
                    {
                        vn += vc + dt*cosz*cos_accel;
                    }
                    v[n][d]        = vn;
                    xprime[n][d]   = x[n][d]+vn*dt;
                } 
                else 
                {
                    xprime[n][d]   = x[n][d];
                }
            }
        }
    } 
    else 
    {
        /* Classic version of update, used with berendsen coupling */
        for(n=start; n<nrend; n++) 
        {
            w_dt = invmass[n]*dt;
            if (cTC) 
            {
                gt   = cTC[n];
            }
            lg   = tcstat[gt].lambda;
            cosz = cos(fac*x[n][ZZ]);
            
            for(d=0; d<DIM; d++) 
            {
                vn             = v[n][d];
                
                if((ptype[n] != eptVSite) && (ptype[n] != eptShell)) 
                {
                    if(d == XX) 
                    {
                        vc           = cosz*vcos;
                        /* Do not scale the cosine velocity profile */
                        vv           = vc + lg*(vn - vc + f[n][d]*w_dt);
                        /* Add the cosine accelaration profile */
                        vv          += dt*cosz*cos_accel;
                    } 
                    else 
                    {
                        vv           = lg*(vn + f[n][d]*w_dt);
                    }
                    v[n][d]        = vv;
                    xprime[n][d]   = x[n][d]+vv*dt;
                } 
                else 
                {
                    v[n][d]        = 0.0;
                    xprime[n][d]   = x[n][d];
                }
            }
        }
    }
}

static gmx_stochd_t *init_stochd(FILE *fplog,t_inputrec *ir)
{
    gmx_stochd_t *sd;
    gmx_sd_const_t *sdc;
    int  ngtc,n;
    real y;
    
    snew(sd,1);

    /* Initiate random number generator for langevin type dynamics,
     * for BD, SD or velocity rescaling temperature coupling.
     */
    sd->gaussrand = gmx_rng_init(ir->ld_seed);

    ngtc = ir->opts.ngtc;

    if (ir->eI == eiBD)
    {
        snew(sd->bd_rf,ngtc);
    }
    else if (EI_SD(ir->eI))
    {
        snew(sd->sdc,ngtc);
        snew(sd->sdsig,ngtc);
    
        sdc = sd->sdc;
        for(n=0; n<ngtc; n++)
        {
            if (ir->opts.tau_t[n] > 0)
            {
                sdc[n].gdt = ir->delta_t/ir->opts.tau_t[n];
                sdc[n].eph = exp(sdc[n].gdt/2);
                sdc[n].emh = exp(-sdc[n].gdt/2);
                sdc[n].em  = exp(-sdc[n].gdt);
            }
            else
            {
                /* No friction and noise on this group */
                sdc[n].gdt = 0;
                sdc[n].eph = 1;
                sdc[n].emh = 1;
                sdc[n].em  = 1;
            }
            if (sdc[n].gdt >= 0.05)
            {
                sdc[n].b = sdc[n].gdt*(sdc[n].eph*sdc[n].eph - 1) 
                    - 4*(sdc[n].eph - 1)*(sdc[n].eph - 1);
                sdc[n].c = sdc[n].gdt - 3 + 4*sdc[n].emh - sdc[n].em;
                sdc[n].d = 2 - sdc[n].eph - sdc[n].emh;
            }
            else
            {
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

void get_stochd_state(gmx_update_t upd,t_state *state)
{
    gmx_rng_get_state(upd->sd->gaussrand,state->ld_rng,state->ld_rngi);
}

void set_stochd_state(gmx_update_t upd,t_state *state)
{
    gmx_rng_set_state(upd->sd->gaussrand,state->ld_rng,state->ld_rngi[0]);
}

gmx_update_t init_update(FILE *fplog,t_inputrec *ir)
{
    t_gmx_update *upd;
    
    snew(upd,1);
    
    if (ir->eI == eiBD || EI_SD(ir->eI) || ir->etc == etcVRESCALE)
    {
        upd->sd = init_stochd(fplog,ir);
    }

    upd->xp = NULL;
    upd->xp_nalloc = 0;

    return upd;
}

static void do_update_sd1(gmx_stochd_t *sd,
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
  if (homenr > sd->sd_V_nalloc) 
  {
      sd->sd_V_nalloc = over_alloc_dd(homenr);
      srenew(sd->sd_V,sd->sd_V_nalloc);
  }
  gaussrand = sd->gaussrand;
  
  for(n=0; n<ngtc; n++) 
  {
      kT = BOLTZ*ref_t[n];
      /* The mass is encounted for later, since this differs per atom */
      sig[n].V  = sqrt(2*kT*(1 - sdc[n].em));
  }
  
  for(n=start; n<start+homenr; n++) 
  {
      ism = sqrt(invmass[n]);
      if (cFREEZE) 
      {
          gf  = cFREEZE[n];
      }
      if (cACC)
      {
          ga  = cACC[n];
      }
      if (cTC)
      {
          gt  = cTC[n];
      }
      
      for(d=0; d<DIM; d++) 
      {
          if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) 
          {
              sd_V = ism*sig[gt].V*gmx_rng_gaussian_table(gaussrand);
              
              v[n][d] = v[n][d]*sdc[gt].em 
                  + (invmass[n]*f[n][d] + accel[ga][d])*tau_t[gt]*(1 - sdc[gt].em)
                  + sd_V;
              
              xprime[n][d] = x[n][d] + v[n][d]*dt;
          } 
          else 
          {
              v[n][d]      = 0.0;
              xprime[n][d] = x[n][d];
          }
      }
  }
}

static void do_update_sd2(gmx_stochd_t *sd,gmx_bool bInitStep,
                          int start,int homenr,
                          rvec accel[],ivec nFreeze[],
                          real invmass[],unsigned short ptype[],
                          unsigned short cFREEZE[],unsigned short cACC[],
                          unsigned short cTC[],
                          rvec x[],rvec xprime[],rvec v[],rvec f[],
                          rvec sd_X[],
                          int ngtc,real tau_t[],real ref_t[],
                          gmx_bool bFirstHalf)
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
  if (homenr > sd->sd_V_nalloc) 
  {
      sd->sd_V_nalloc = over_alloc_dd(homenr);
      srenew(sd->sd_V,sd->sd_V_nalloc);
  }
  sd_V = sd->sd_V;
  gaussrand = sd->gaussrand;

  if (bFirstHalf) 
  {
      for (n=0; n<ngtc; n++) 
      {
          kT = BOLTZ*ref_t[n];
          /* The mass is encounted for later, since this differs per atom */
          sig[n].V  = sqrt(kT*(1-sdc[n].em));
          sig[n].X  = sqrt(kT*sqr(tau_t[n])*sdc[n].c);
          sig[n].Yv = sqrt(kT*sdc[n].b/sdc[n].c);
          sig[n].Yx = sqrt(kT*sqr(tau_t[n])*sdc[n].b/(1-sdc[n].em));
      }
  }
  
  for (n=start; n<start+homenr; n++) 
  {
      ism = sqrt(invmass[n]);
      if (cFREEZE) 
      {
          gf  = cFREEZE[n];
      }
      if (cACC) 
      {
          ga  = cACC[n];
      }
      if (cTC)
      {
          gt  = cTC[n];
      }
      
      for(d=0; d<DIM; d++) 
      {
          if (bFirstHalf) 
          {
              vn             = v[n][d];
          }
          if((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d]) 
          {
              if (bFirstHalf) 
              {
                  if (bInitStep) 
                  {
                      sd_X[n][d] = ism*sig[gt].X*gmx_rng_gaussian_table(gaussrand);
                  }
                  Vmh = sd_X[n][d]*sdc[gt].d/(tau_t[gt]*sdc[gt].c) 
                      + ism*sig[gt].Yv*gmx_rng_gaussian_table(gaussrand);
                  sd_V[n-start][d] = ism*sig[gt].V*gmx_rng_gaussian_table(gaussrand);
                  
                  v[n][d] = vn*sdc[gt].em 
                      + (invmass[n]*f[n][d] + accel[ga][d])*tau_t[gt]*(1 - sdc[gt].em)
                      + sd_V[n-start][d] - sdc[gt].em*Vmh;
                  
                  xprime[n][d] = x[n][d] + v[n][d]*tau_t[gt]*(sdc[gt].eph - sdc[gt].emh); 
              } 
              else 
              {
                  
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
          } 
          else 
          {
              if (bFirstHalf) 
              {
                  v[n][d]        = 0.0;
                  xprime[n][d]   = x[n][d];
              }
          }
      }
  }
}

static void do_update_bd(int start,int nrend,double dt,
                         ivec nFreeze[],
                         real invmass[],unsigned short ptype[],
                         unsigned short cFREEZE[],unsigned short cTC[],
                         rvec x[],rvec xprime[],rvec v[],
                         rvec f[],real friction_coefficient,
                         int ngtc,real tau_t[],real ref_t[],
                         real *rf,gmx_rng_t gaussrand)
{
    /* note -- these appear to be full step velocities . . .  */
    int    gf=0,gt=0;
    real   vn;
    real   invfr=0;
    int    n,d;
    
    if (friction_coefficient != 0) 
    {
        invfr = 1.0/friction_coefficient;
        for(n=0; n<ngtc; n++) 
        {
            rf[n] = sqrt(2.0*BOLTZ*ref_t[n]/(friction_coefficient*dt));
        } 
    }
    else 
    {
        for(n=0; n<ngtc; n++) 
        {
            rf[n] = sqrt(2.0*BOLTZ*ref_t[n]);
        }
    }
    for(n=start; (n<nrend); n++) 
    {
        if (cFREEZE) 
        {
            gf = cFREEZE[n];
        }
        if (cTC)
        {
            gt = cTC[n];
        }
        for(d=0; (d<DIM); d++) 
        {
            if((ptype[n]!=eptVSite) && (ptype[n]!=eptShell) && !nFreeze[gf][d]) 
            {
                if (friction_coefficient != 0) {
                    vn = invfr*f[n][d] + rf[gt]*gmx_rng_gaussian_table(gaussrand);
                } 
                else 
                {
                    /* NOTE: invmass = 1/(mass*friction_constant*dt) */
                    vn = invmass[n]*f[n][d]*dt 
                        + sqrt(invmass[n])*rf[gt]*gmx_rng_gaussian_table(gaussrand);
                }

                v[n][d]      = vn;
                xprime[n][d] = x[n][d]+vn*dt;
            }
            else 
            {
                v[n][d]      = 0.0;
                xprime[n][d] = x[n][d];
            }
        }
    }
}

static void dump_it_all(FILE *fp,const char *title,
                        int natoms,rvec x[],rvec xp[],rvec v[],rvec f[])
{
#ifdef DEBUG
  if (fp) 
  {
    fprintf(fp,"%s\n",title);
    pr_rvecs(fp,0,"x",x,natoms);
    pr_rvecs(fp,0,"xp",xp,natoms);
    pr_rvecs(fp,0,"v",v,natoms);
    pr_rvecs(fp,0,"f",f,natoms);
  }
#endif
}

static void calc_ke_part_normal(rvec v[], t_grpopts *opts,t_mdatoms *md,
                                gmx_ekindata_t *ekind,t_nrnb *nrnb,gmx_bool bEkinAveVel, 
                                gmx_bool bSaveEkinOld)
{
  int          start=md->start,homenr=md->homenr;
  int          g,d,n,m,ga=0,gt=0;
  rvec         v_corrt;
  real         hm;
  t_grp_tcstat *tcstat=ekind->tcstat;
  t_grp_acc    *grpstat=ekind->grpstat;
  real         dekindl;

  /* three main: VV with AveVel, vv with AveEkin, leap with AveEkin.  Leap with AveVel is also
     an option, but not supported now.  Additionally, if we are doing iterations.  
     bEkinAveVel: If TRUE, we sum into ekin, if FALSE, into ekinh.
     bSavEkinOld: If TRUE (in the case of iteration = bIterate is TRUE), we don't copy over the ekinh_old.  
     If FALSE, we overrwrite it.
  */

  /* group velocities are calculated in update_ekindata and
   * accumulated in acumulate_groups.
   * Now the partial global and groups ekin.
   */
  for(g=0; (g<opts->ngtc); g++) 
  {
      
      if (!bSaveEkinOld) {
          copy_mat(tcstat[g].ekinh,tcstat[g].ekinh_old);
      } 
      if(bEkinAveVel) {
          clear_mat(tcstat[g].ekinf);
      } else {
          clear_mat(tcstat[g].ekinh);
      }
      if (bEkinAveVel) {
          tcstat[g].ekinscalef_nhc = 1.0;   /* need to clear this -- logic is complicated! */
      }
  }
  ekind->dekindl_old = ekind->dekindl;
  
  dekindl = 0;
  for(n=start; (n<start+homenr); n++) 
  {
      if (md->cACC)
      {
          ga = md->cACC[n];
      }
      if (md->cTC)
      {
          gt = md->cTC[n];
      }
      hm   = 0.5*md->massT[n];
      
      for(d=0; (d<DIM); d++) 
      {
          v_corrt[d]  = v[n][d]  - grpstat[ga].u[d];
      }
      for(d=0; (d<DIM); d++) 
      {
          for (m=0;(m<DIM); m++) 
          {
              /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
              if (bEkinAveVel) 
              {
                  tcstat[gt].ekinf[m][d]+=hm*v_corrt[m]*v_corrt[d];
              } 
              else 
              {
                  tcstat[gt].ekinh[m][d]+=hm*v_corrt[m]*v_corrt[d]; 
              }
          }
      }
      if (md->nMassPerturbed && md->bPerturbed[n]) 
      {
          dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt,v_corrt);
      }
  }
  ekind->dekindl = dekindl;
  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

static void calc_ke_part_visc(matrix box,rvec x[],rvec v[],
                              t_grpopts *opts,t_mdatoms *md,
                              gmx_ekindata_t *ekind,
                              t_nrnb *nrnb, gmx_bool bEkinAveVel, gmx_bool bSaveEkinOld)
{
  int          start=md->start,homenr=md->homenr;
  int          g,d,n,m,gt=0;
  rvec         v_corrt;
  real         hm;
  t_grp_tcstat *tcstat=ekind->tcstat;
  t_cos_acc    *cosacc=&(ekind->cosacc);
  real         dekindl;
  real         fac,cosz;
  double       mvcos;

  for(g=0; g<opts->ngtc; g++) 
  {
      copy_mat(ekind->tcstat[g].ekinh,ekind->tcstat[g].ekinh_old);
      clear_mat(ekind->tcstat[g].ekinh);
  }
  ekind->dekindl_old = ekind->dekindl;

  fac = 2*M_PI/box[ZZ][ZZ];
  mvcos = 0;
  dekindl = 0;
  for(n=start; n<start+homenr; n++) 
  {
      if (md->cTC) 
      {
          gt = md->cTC[n];
      }
      hm   = 0.5*md->massT[n];
      
      /* Note that the times of x and v differ by half a step */
      /* MRS -- would have to be changed for VV */
      cosz         = cos(fac*x[n][ZZ]);
      /* Calculate the amplitude of the new velocity profile */
      mvcos       += 2*cosz*md->massT[n]*v[n][XX];
      
      copy_rvec(v[n],v_corrt);
      /* Subtract the profile for the kinetic energy */
      v_corrt[XX] -= cosz*cosacc->vcos;
      for (d=0; (d<DIM); d++) 
      {
          for (m=0; (m<DIM); m++) 
          {
              /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
              if (bEkinAveVel) 
              {
                  tcstat[gt].ekinf[m][d]+=hm*v_corrt[m]*v_corrt[d];
              } 
              else 
              {
                  tcstat[gt].ekinh[m][d]+=hm*v_corrt[m]*v_corrt[d];
              }
          }
      }
      if(md->nPerturbed && md->bPerturbed[n]) 
      {
          dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt,v_corrt);
      }
  }
  ekind->dekindl = dekindl;
  cosacc->mvcos = mvcos;
  
  inc_nrnb(nrnb,eNR_EKIN,homenr);
}

void calc_ke_part(t_state *state,t_grpopts *opts,t_mdatoms *md,
                  gmx_ekindata_t *ekind,t_nrnb *nrnb, gmx_bool bEkinAveVel, gmx_bool bSaveEkinOld)
{
    if (ekind->cosacc.cos_accel == 0)
    {
        calc_ke_part_normal(state->v,opts,md,ekind,nrnb,bEkinAveVel,bSaveEkinOld);
    }
    else
    {
        calc_ke_part_visc(state->box,state->x,state->v,opts,md,ekind,nrnb,bEkinAveVel,bSaveEkinOld);
    }
}

void init_ekinstate(ekinstate_t *ekinstate,const t_inputrec *ir)
{
    ekinstate->ekin_n = ir->opts.ngtc;
    snew(ekinstate->ekinh,ekinstate->ekin_n);
    snew(ekinstate->ekinf,ekinstate->ekin_n);
    snew(ekinstate->ekinh_old,ekinstate->ekin_n);
    snew(ekinstate->ekinscalef_nhc,ekinstate->ekin_n);
    snew(ekinstate->ekinscaleh_nhc,ekinstate->ekin_n);
    snew(ekinstate->vscale_nhc,ekinstate->ekin_n);
    ekinstate->dekindl = 0;
    ekinstate->mvcos   = 0;
}

void update_ekinstate(ekinstate_t *ekinstate,gmx_ekindata_t *ekind)
{
  int i;
  
  for(i=0;i<ekinstate->ekin_n;i++) 
  {
      copy_mat(ekind->tcstat[i].ekinh,ekinstate->ekinh[i]);
      copy_mat(ekind->tcstat[i].ekinf,ekinstate->ekinf[i]); 
      copy_mat(ekind->tcstat[i].ekinh_old,ekinstate->ekinh_old[i]); 
      ekinstate->ekinscalef_nhc[i] = ekind->tcstat[i].ekinscalef_nhc;
      ekinstate->ekinscaleh_nhc[i] = ekind->tcstat[i].ekinscaleh_nhc;
      ekinstate->vscale_nhc[i] = ekind->tcstat[i].vscale_nhc;
  }

  copy_mat(ekind->ekin,ekinstate->ekin_total);
  ekinstate->dekindl = ekind->dekindl;
  ekinstate->mvcos = ekind->cosacc.mvcos;
  
}

void restore_ekinstate_from_state(t_commrec *cr,
                                  gmx_ekindata_t *ekind,ekinstate_t *ekinstate)
{
  int i,n;

  if (MASTER(cr)) 
  {
      for(i=0;i<ekinstate->ekin_n;i++) 
      {
          copy_mat(ekinstate->ekinh[i],ekind->tcstat[i].ekinh);
          copy_mat(ekinstate->ekinf[i],ekind->tcstat[i].ekinf);
          copy_mat(ekinstate->ekinh_old[i],ekind->tcstat[i].ekinh_old);
          ekind->tcstat[i].ekinscalef_nhc = ekinstate->ekinscalef_nhc[i];
          ekind->tcstat[i].ekinscaleh_nhc = ekinstate->ekinscaleh_nhc[i];
          ekind->tcstat[i].vscale_nhc = ekinstate->vscale_nhc[i];
      }
      
      copy_mat(ekinstate->ekin_total,ekind->ekin);

      ekind->dekindl = ekinstate->dekindl;
      ekind->cosacc.mvcos = ekinstate->mvcos;
      n = ekinstate->ekin_n;
  }
 
  if (PAR(cr)) 
  {
      gmx_bcast(sizeof(n),&n,cr);
      for(i=0;i<n;i++) 
      {
          gmx_bcast(DIM*DIM*sizeof(ekind->tcstat[i].ekinh[0][0]),
                    ekind->tcstat[i].ekinh[0],cr);
          gmx_bcast(DIM*DIM*sizeof(ekind->tcstat[i].ekinf[0][0]),
                    ekind->tcstat[i].ekinf[0],cr);
          gmx_bcast(DIM*DIM*sizeof(ekind->tcstat[i].ekinh_old[0][0]),
                    ekind->tcstat[i].ekinh_old[0],cr);

          gmx_bcast(sizeof(ekind->tcstat[i].ekinscalef_nhc),
                    &(ekind->tcstat[i].ekinscalef_nhc),cr);
          gmx_bcast(sizeof(ekind->tcstat[i].ekinscaleh_nhc),
                    &(ekind->tcstat[i].ekinscaleh_nhc),cr);
          gmx_bcast(sizeof(ekind->tcstat[i].vscale_nhc),
                    &(ekind->tcstat[i].vscale_nhc),cr);
      }
      gmx_bcast(DIM*DIM*sizeof(ekind->ekin[0][0]),
                ekind->ekin[0],cr);

      gmx_bcast(sizeof(ekind->dekindl),&ekind->dekindl,cr);
      gmx_bcast(sizeof(ekind->cosacc.mvcos),&ekind->cosacc.mvcos,cr);
  }
}

void set_deform_reference_box(gmx_update_t upd,gmx_large_int_t step,matrix box)
{
    upd->deformref_step = step;
    copy_mat(box,upd->deformref_box);
}

static void deform(gmx_update_t upd,
                   int start,int homenr,rvec x[],matrix box,matrix *scale_tot,
                   const t_inputrec *ir,gmx_large_int_t step)
{
    matrix bnew,invbox,mu;
    real   elapsed_time;
    int    i,j;  
    
    elapsed_time = (step + 1 - upd->deformref_step)*ir->delta_t;
    copy_mat(box,bnew);
    for(i=0; i<DIM; i++)
    {
        for(j=0; j<DIM; j++)
        {
            if (ir->deform[i][j] != 0)
            {
                bnew[i][j] =
                    upd->deformref_box[i][j] + elapsed_time*ir->deform[i][j];
            }
        }
    }
    /* We correct the off-diagonal elements,
     * which can grow indefinitely during shearing,
     * so the shifts do not get messed up.
     */
    for(i=1; i<DIM; i++)
    {
        for(j=i-1; j>=0; j--)
        {
            while (bnew[i][j] - box[i][j] > 0.5*bnew[j][j])
            {
                rvec_dec(bnew[i],bnew[j]);
            }
            while (bnew[i][j] - box[i][j] < -0.5*bnew[j][j])
            {
                rvec_inc(bnew[i],bnew[j]);
            }
        }
    }
    m_inv_ur0(box,invbox);
    copy_mat(bnew,box);
    mmul_ur0(box,invbox,mu);
  
    for(i=start; i<start+homenr; i++)
    {
        x[i][XX] = mu[XX][XX]*x[i][XX]+mu[YY][XX]*x[i][YY]+mu[ZZ][XX]*x[i][ZZ];
        x[i][YY] = mu[YY][YY]*x[i][YY]+mu[ZZ][YY]*x[i][ZZ];
        x[i][ZZ] = mu[ZZ][ZZ]*x[i][ZZ];
    }
    if (*scale_tot)
    {
        /* The transposes of the scaling matrices are stored,
         * so we need to do matrix multiplication in the inverse order.
         */
        mmul_ur0(*scale_tot,mu,*scale_tot);
    }
}

static void combine_forces(int nstlist,
                           gmx_constr_t constr,
                           t_inputrec *ir,t_mdatoms *md,t_idef *idef,
                           t_commrec *cr,gmx_large_int_t step,t_state *state,
                           int start,int nrend,
                           rvec f[],rvec f_lr[],
                           t_nrnb *nrnb)
{
    int  i,d,nm1;

    /* f contains the short-range forces + the long range forces
     * which are stored separately in f_lr.
     */

    if (constr != NULL && !(ir->eConstrAlg == econtSHAKE && ir->epc == epcNO))
    {
        /* We need to constrain the LR forces separately,
         * because due to the different pre-factor for the SR and LR
         * forces in the update algorithm, we can not determine
         * the constraint force for the coordinate constraining.
         * Constrain only the additional LR part of the force.
         */
        /* MRS -- need to make sure this works with trotter integration -- the constraint calls may not be right.*/
        constrain(NULL,FALSE,FALSE,constr,idef,ir,NULL,cr,step,0,md,
                  state->x,f_lr,f_lr,state->box,state->lambda,NULL,
                  NULL,NULL,nrnb,econqForce,ir->epc==epcMTTK,state->veta,state->veta);
    }
    
    /* Add nstlist-1 times the LR force to the sum of both forces
     * and store the result in forces_lr.
     */
    nm1 = nstlist - 1;
    for(i=start; i<nrend; i++)
    {
        for(d=0; d<DIM; d++)
        {
            f_lr[i][d] = f[i][d] + nm1*f_lr[i][d];
        }
    }
}

void update_tcouple(FILE         *fplog,
                    gmx_large_int_t   step,
                    t_inputrec   *inputrec,   
                    t_state      *state,
                    gmx_ekindata_t *ekind,
                    gmx_wallcycle_t wcycle,
                    gmx_update_t upd,
                    t_extmass    *MassQ,
                    t_mdatoms  *md)
    
{
    gmx_bool   bTCouple=FALSE;
    real   dttc;
    int    i,start,end,homenr;
    
    /* if using vv, we do this elsewhere in the code */
    if (inputrec->etc != etcNO &&
        !(IR_NVT_TROTTER(inputrec) || IR_NPT_TROTTER(inputrec)))
    {
        /* We should only couple after a step where energies were determined */
        bTCouple = (inputrec->nsttcouple == 1 ||
                    do_per_step(step+inputrec->nsttcouple-1,
                                inputrec->nsttcouple));
    }
    
    if (bTCouple)
    {
        dttc = inputrec->nsttcouple*inputrec->delta_t;

        switch (inputrec->etc) 
        {
        case etcNO:
            break;
        case etcBERENDSEN:
            berendsen_tcoupl(inputrec,ekind,dttc);
            break;
        case etcNOSEHOOVER:
            nosehoover_tcoupl(&(inputrec->opts),ekind,dttc,
                              state->nosehoover_xi,state->nosehoover_vxi,MassQ);
            break;
        case etcVRESCALE:
            vrescale_tcoupl(inputrec,ekind,dttc,
                            state->therm_integral,upd->sd->gaussrand);
            break;
        }
        /* rescale in place here */
        if (EI_VV(inputrec->eI))
        {
            rescale_velocities(ekind,md,md->start,md->start+md->homenr,state->v);
        }
    }
    else 
    {
        /* Set the T scaling lambda to 1 to have no scaling */
        for(i=0; (i<inputrec->opts.ngtc); i++)
        {
            ekind->tcstat[i].lambda = 1.0;
        }
    }
}

void update_pcouple(FILE         *fplog,
                    gmx_large_int_t   step,
                    t_inputrec   *inputrec,   
                    t_state      *state,
                    matrix       pcoupl_mu,
                    matrix       M,
                    gmx_wallcycle_t wcycle,
                    gmx_update_t upd,
                    gmx_bool         bInitStep)
{
    gmx_bool   bPCouple=FALSE;
    real   dtpc=0;
    int    i;
    
    /* if using vv, we do this elsewhere in the code */
    if (inputrec->epc != epcNO &&
        !(IR_NVT_TROTTER(inputrec) || IR_NPT_TROTTER(inputrec)))
    {
        /* We should only couple after a step where energies were determined */
        bPCouple = (inputrec->nstpcouple == 1 ||
                    do_per_step(step+inputrec->nstpcouple-1,
                                inputrec->nstpcouple));
    }
    
    clear_mat(pcoupl_mu);
    for(i=0; i<DIM; i++)
    {
        pcoupl_mu[i][i] = 1.0;
    }
    
    clear_mat(M);
     
    if (bPCouple)
    {
        dtpc = inputrec->nstpcouple*inputrec->delta_t;

        switch (inputrec->epc) 
        {
            /* We can always pcoupl, even if we did not sum the energies
             * the previous step, since state->pres_prev is only updated
             * when the energies have been summed.
             */
        case (epcNO):
            break;
        case (epcBERENDSEN):
            if (!bInitStep) 
            {
                berendsen_pcoupl(fplog,step,inputrec,dtpc,state->pres_prev,state->box,
                                 pcoupl_mu);
            }
            break;
        case (epcPARRINELLORAHMAN):
            parrinellorahman_pcoupl(fplog,step,inputrec,dtpc,state->pres_prev,
                                    state->box,state->box_rel,state->boxv,
                                    M,pcoupl_mu,bInitStep);
            break;
        default:
            break;
        }  
    } 
}

static rvec *get_xprime(const t_state *state,gmx_update_t upd)
{
    if (state->nalloc > upd->xp_nalloc)
    {
        upd->xp_nalloc = state->nalloc;
        srenew(upd->xp,upd->xp_nalloc);
    }
 
    return upd->xp;
}

void update_constraints(FILE         *fplog,
                        gmx_large_int_t   step,
                        real         *dvdlambda,    /* FEP stuff */
                        t_inputrec   *inputrec,      /* input record and box stuff	*/
                        gmx_ekindata_t *ekind,
                        t_mdatoms    *md,
                        t_state      *state,
                        t_graph      *graph,  
                        rvec         force[],        /* forces on home particles */
                        t_idef       *idef,
                        tensor       vir_part,
                        tensor       vir,            /* tensors for virial and ekin, needed for computing */
                        t_commrec    *cr,
                        t_nrnb       *nrnb,
                        gmx_wallcycle_t wcycle,
                        gmx_update_t upd,
                        gmx_constr_t constr,
                        gmx_bool         bInitStep,
                        gmx_bool         bFirstHalf,
                        gmx_bool         bCalcVir,
                        real         vetanew)
{
    gmx_bool             bExtended,bTrotter,bLastStep,bLog=FALSE,bEner=FALSE,bDoConstr=FALSE;
    double           dt;
    real             dt_1;
    int              start,homenr,nrend,i,n,m,g,d;
    tensor           vir_con;
    rvec             *vbuf,*xprime=NULL;
    
    if (constr) {bDoConstr=TRUE;}
    if (bFirstHalf && !EI_VV(inputrec->eI)) {bDoConstr=FALSE;} 

    /* for now, SD update is here -- though it really seems like it 
       should be reformulated as a velocity verlet method, since it has two parts */
    
    start  = md->start;
    homenr = md->homenr;
    nrend = start+homenr;
    
    dt   = inputrec->delta_t;
    dt_1 = 1.0/dt;
    
    /* 
     *  Steps (7C, 8C)
     *  APPLY CONSTRAINTS:
     *  BLOCK SHAKE 

     * When doing PR pressure coupling we have to constrain the
     * bonds in each iteration. If we are only using Nose-Hoover tcoupling
     * it is enough to do this once though, since the relative velocities 
     * after this will be normal to the bond vector
     */
    
    if (bDoConstr) 
    {
        /* clear out constraints before applying */
        clear_mat(vir_part);

        xprime = get_xprime(state,upd);

        bLastStep = (step == inputrec->init_step+inputrec->nsteps);
        bLog  = (do_per_step(step,inputrec->nstlog) || bLastStep || (step < 0));
        bEner = (do_per_step(step,inputrec->nstenergy) || bLastStep);
        /* Constrain the coordinates xprime */
        wallcycle_start(wcycle,ewcCONSTR);
        if (EI_VV(inputrec->eI) && bFirstHalf) 
        {
            constrain(NULL,bLog,bEner,constr,idef,
                      inputrec,ekind,cr,step,1,md,
                      state->x,state->v,state->v,
                      state->box,state->lambda,dvdlambda,
                      NULL,bCalcVir ? &vir_con : NULL,nrnb,econqVeloc,
                      inputrec->epc==epcMTTK,state->veta,vetanew);
        } 
        else 
        {
            constrain(NULL,bLog,bEner,constr,idef,
                      inputrec,ekind,cr,step,1,md,
                      state->x,xprime,NULL,
                      state->box,state->lambda,dvdlambda,
                      state->v,bCalcVir ? &vir_con : NULL ,nrnb,econqCoord,
                      inputrec->epc==epcMTTK,state->veta,state->veta);
        }
        wallcycle_stop(wcycle,ewcCONSTR);
        
        where();
    
        dump_it_all(fplog,"After Shake",
                    state->natoms,state->x,xprime,state->v,force);
        
        if (bCalcVir) 
        {
            if (inputrec->eI == eiSD2) 
            {
                /* A correction factor eph is needed for the SD constraint force */
                /* Here we can, unfortunately, not have proper corrections
                 * for different friction constants, so we use the first one.
                 */
                for(i=0; i<DIM; i++) 
                {
                    for(m=0; m<DIM; m++)
                    {
                        vir_part[i][m] += upd->sd->sdc[0].eph*vir_con[i][m];
                    }
                }
            } 
            else 
            {
                m_add(vir_part,vir_con,vir_part);
            }
            if (debug) 
            {
                pr_rvecs(debug,0,"constraint virial",vir_part,DIM);
            }
        }
    }
    
    where();
    if ((inputrec->eI == eiSD2) && !(bFirstHalf))
    {
        xprime = get_xprime(state,upd);

        /* The second part of the SD integration */
        do_update_sd2(upd->sd,FALSE,start,homenr,
                      inputrec->opts.acc,inputrec->opts.nFreeze,
                      md->invmass,md->ptype,
                      md->cFREEZE,md->cACC,md->cTC,
                      state->x,xprime,state->v,force,state->sd_X,
                      inputrec->opts.ngtc,inputrec->opts.tau_t,
                      inputrec->opts.ref_t,FALSE);
        inc_nrnb(nrnb, eNR_UPDATE, homenr);
        
        if (bDoConstr) 
        {
            /* Constrain the coordinates xprime */
            wallcycle_start(wcycle,ewcCONSTR);
            constrain(NULL,bLog,bEner,constr,idef,
                      inputrec,NULL,cr,step,1,md,
                      state->x,xprime,NULL,
                      state->box,state->lambda,dvdlambda,
                      NULL,NULL,nrnb,econqCoord,FALSE,0,0);
            wallcycle_stop(wcycle,ewcCONSTR);
        }
    }    

    /* We must always unshift after updating coordinates; if we did not shake
       x was shifted in do_force */
    
    if (!(bFirstHalf)) /* in the first half of vv, no shift. */
    {
        if (graph && (graph->nnodes > 0)) 
        {
            unshift_x(graph,state->box,state->x,upd->xp);
            if (TRICLINIC(state->box)) 
            {
                inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
            }
            else
            {
                inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);    
            }
            copy_rvecn(upd->xp,state->x,start,graph->start);
            copy_rvecn(upd->xp,state->x,graph->start+graph->nnodes,nrend);
        }
        else 
        {
            copy_rvecn(upd->xp,state->x,start,nrend);
        }
        
        dump_it_all(fplog,"After unshift",
                    state->natoms,state->x,upd->xp,state->v,force);
    }
/* ############# END the update of velocities and positions ######### */
}

void update_box(FILE         *fplog,
                gmx_large_int_t   step,
                t_inputrec   *inputrec,      /* input record and box stuff	*/
                t_mdatoms    *md,
                t_state      *state,
                t_graph      *graph,
                rvec         force[],        /* forces on home particles */
                matrix       *scale_tot,
                matrix       pcoupl_mu,
                t_nrnb       *nrnb,
                gmx_wallcycle_t wcycle,
                gmx_update_t upd,
                gmx_bool         bInitStep,
                gmx_bool         bFirstHalf)
{
    gmx_bool             bExtended,bTrotter,bLastStep,bLog=FALSE,bEner=FALSE;
    double           dt;
    real             dt_1;
    int              start,homenr,nrend,i,n,m,g;
    tensor           vir_con;
    
    start  = md->start;
    homenr = md->homenr;
    nrend = start+homenr;
    
    bExtended =
        (inputrec->etc == etcNOSEHOOVER) ||
        (inputrec->epc == epcPARRINELLORAHMAN) ||
        (inputrec->epc == epcMTTK);
    
    dt = inputrec->delta_t;
    
    where();

    /* now update boxes */
    switch (inputrec->epc) {
    case (epcNO):
        break;
    case (epcBERENDSEN):
        berendsen_pscale(inputrec,pcoupl_mu,state->box,state->box_rel,
                         start,homenr,state->x,md->cFREEZE,nrnb);
        break;
    case (epcPARRINELLORAHMAN): 
        /* The box velocities were updated in do_pr_pcoupl in the update
         * iteration, but we dont change the box vectors until we get here
         * since we need to be able to shift/unshift above.
         */
        for(i=0;i<DIM;i++)
        {
            for(m=0;m<=i;m++)
            {
                state->box[i][m] += dt*state->boxv[i][m];
            }
        }
        preserve_box_shape(inputrec,state->box_rel,state->box);
        
        /* Scale the coordinates */
        for(n=start; (n<start+homenr); n++) 
        {
            tmvmul_ur0(pcoupl_mu,state->x[n],state->x[n]);
        }
        break;
    case (epcMTTK):
        switch (inputrec->epct) 
        {
        case (epctISOTROPIC):
            /* DIM * eta = ln V.  so DIM*eta_new = DIM*eta_old + DIM*dt*veta => 
               ln V_new = ln V_old + 3*dt*veta => V_new = V_old*exp(3*dt*veta) => 
               Side length scales as exp(veta*dt) */
            
            msmul(state->box,exp(state->veta*dt),state->box);
            
            /* Relate veta to boxv.  veta = d(eta)/dT = (1/DIM)*1/V dV/dT.
o               If we assume isotropic scaling, and box length scaling
               factor L, then V = L^DIM (det(M)).  So dV/dt = DIM
               L^(DIM-1) dL/dt det(M), and veta = (1/L) dL/dt.  The
               determinant of B is L^DIM det(M), and the determinant
               of dB/dt is (dL/dT)^DIM det (M).  veta will be
               (det(dB/dT)/det(B))^(1/3).  Then since M =
               B_new*(vol_new)^(1/3), dB/dT_new = (veta_new)*B(new). */
            
            msmul(state->box,state->veta,state->boxv);
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }
    
    if ((!(IR_NPT_TROTTER(inputrec))) && scale_tot) 
    {
        /* The transposes of the scaling matrices are stored,
         * therefore we need to reverse the order in the multiplication.
         */
        mmul_ur0(*scale_tot,pcoupl_mu,*scale_tot);
    }

    if (DEFORM(*inputrec)) 
    {
        deform(upd,start,homenr,state->x,state->box,scale_tot,inputrec,step);
    }
    where();
    dump_it_all(fplog,"After update",
                state->natoms,state->x,upd->xp,state->v,force);
}

void update_coords(FILE         *fplog,
                   gmx_large_int_t   step,
                   t_inputrec   *inputrec,      /* input record and box stuff	*/
                   t_mdatoms    *md,
                   t_state      *state,
                   rvec         *f,        /* forces on home particles */
                   gmx_bool         bDoLR,
                   rvec         *f_lr,
                   t_fcdata     *fcd,
                   gmx_ekindata_t *ekind,
                   matrix       M,
                   gmx_wallcycle_t wcycle,
                   gmx_update_t upd,
                   gmx_bool         bInitStep,
                   int          UpdatePart,
                   t_commrec    *cr,  /* these shouldn't be here -- need to think about it */
                   t_nrnb       *nrnb,
                   gmx_constr_t constr,
                   t_idef       *idef)
{
    gmx_bool             bExtended,bNH,bPR,bTrotter,bLastStep,bLog=FALSE,bEner=FALSE;
    double           dt,alpha;
    real             *imass,*imassin;
    rvec             *force;
    real             dt_1;
    int              start,homenr,nrend,i,j,d,n,m,g;
    int              blen0,blen1,iatom,jatom,nshake,nsettle,nconstr,nexpand;
    int              *icom = NULL;
    tensor           vir_con;
    rvec             *vcom,*xcom,*vall,*xall,*xin,*vin,*forcein,*fall,*xpall,*xprimein,*xprime;
    

    /* Running the velocity half does nothing except for velocity verlet */
    if ((UpdatePart == etrtVELOCITY1 || UpdatePart == etrtVELOCITY2) &&
        !EI_VV(inputrec->eI))
    {
        gmx_incons("update_coords called for velocity without VV integrator");
    }

    start  = md->start;
    homenr = md->homenr;
    nrend = start+homenr;
    
    xprime = get_xprime(state,upd);
    
    dt   = inputrec->delta_t;
    dt_1 = 1.0/dt;

    /* We need to update the NMR restraint history when time averaging is used */
    if (state->flags & (1<<estDISRE_RM3TAV)) 
    {
        update_disres_history(fcd,&state->hist);
    }
    if (state->flags & (1<<estORIRE_DTAV)) 
    {
        update_orires_history(fcd,&state->hist);
    }
    
    bNH = inputrec->etc == etcNOSEHOOVER;
    bPR = ((inputrec->epc == epcPARRINELLORAHMAN) || (inputrec->epc == epcMTTK)); 

    bExtended = bNH || bPR;
    
    if (bDoLR && inputrec->nstlist > 1 && !EI_VV(inputrec->eI))  /* get this working with VV? */
    {
        /* Store the total force + nstlist-1 times the LR force
         * in forces_lr, so it can be used in a normal update algorithm
         * to produce twin time stepping.
         */
        /* is this correct in the new construction? MRS */
        combine_forces(inputrec->nstlist,constr,inputrec,md,idef,cr,step,state,
                       start,nrend,f,f_lr,nrnb);
        force = f_lr;
    }
    else
    {
        force = f;
    }

    /* ############# START The update of velocities and positions ######### */
    where();
    dump_it_all(fplog,"Before update",
                state->natoms,state->x,xprime,state->v,force);
    
    switch (inputrec->eI) {
    case (eiMD):
        if (ekind->cosacc.cos_accel == 0) {
            /* use normal version of update */
            do_update_md(start,nrend,dt,
                         ekind->tcstat,ekind->grpstat,state->nosehoover_vxi,
                         inputrec->opts.acc,inputrec->opts.nFreeze,md->invmass,md->ptype,
                         md->cFREEZE,md->cACC,md->cTC,
                         state->x,xprime,state->v,force,M,
                         bNH,bPR);
        } 
        else 
        {
            do_update_visc(start,nrend,dt,
                           ekind->tcstat,md->invmass,state->nosehoover_vxi,
                           md->ptype,md->cTC,state->x,xprime,state->v,force,M,
                           state->box,ekind->cosacc.cos_accel,ekind->cosacc.vcos,bNH,bPR);
        }
        break;
    case (eiSD1):
        do_update_sd1(upd->sd,start,homenr,dt,
                      inputrec->opts.acc,inputrec->opts.nFreeze,
                      md->invmass,md->ptype,
                      md->cFREEZE,md->cACC,md->cTC,
                      state->x,xprime,state->v,force,state->sd_X,
                      inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t);
        break;
    case (eiSD2):
        /* The SD update is done in 2 parts, because an extra constraint step
         * is needed 
         */
        do_update_sd2(upd->sd,bInitStep,start,homenr,
                      inputrec->opts.acc,inputrec->opts.nFreeze,
                      md->invmass,md->ptype,
                      md->cFREEZE,md->cACC,md->cTC,
                      state->x,xprime,state->v,force,state->sd_X,
                      inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
                      TRUE);
        break;
    case (eiBD):
        do_update_bd(start,nrend,dt,
                     inputrec->opts.nFreeze,md->invmass,md->ptype,
                     md->cFREEZE,md->cTC,
                     state->x,xprime,state->v,force,
                     inputrec->bd_fric,
                     inputrec->opts.ngtc,inputrec->opts.tau_t,inputrec->opts.ref_t,
                     upd->sd->bd_rf,upd->sd->gaussrand);
        break;
    case (eiVV):
    case (eiVVAK):
        alpha = 1.0 + DIM/((double)inputrec->opts.nrdf[0]); /* assuming barostat coupled to group 0. */
        switch (UpdatePart) {
        case etrtVELOCITY1:
        case etrtVELOCITY2:
            do_update_vv_vel(start,nrend,dt,
                             ekind->tcstat,ekind->grpstat,
                             inputrec->opts.acc,inputrec->opts.nFreeze,
                             md->invmass,md->ptype,
                             md->cFREEZE,md->cACC,
                             state->v,force,
                             bExtended,state->veta,alpha);  
            break;
        case etrtPOSITION:
            do_update_vv_pos(start,nrend,dt,
                             ekind->tcstat,ekind->grpstat,
                             inputrec->opts.acc,inputrec->opts.nFreeze,
                             md->invmass,md->ptype,md->cFREEZE,
                             state->x,xprime,state->v,force,
                             bExtended,state->veta,alpha);
            break;
        }
        break;
    default:
        gmx_fatal(FARGS,"Don't know how to update coordinates");
        break;
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
  for(i=start; (i<end); i++) 
  {
    m      = mass[i];
    tm    += m;
    for(j=0; (j<DIM); j++) 
    {
        mv[j] += m*v[i][j];
    }
  }
  /* Shortcut */ 
  svmul(1/tmass,vcm,vcm); 
  svmul(0.5,vcm,hvcm);
  clear_mat(dekin);
  for(j=0; (j<DIM); j++) 
  {
      for(k=0; (k<DIM); k++) 
      {
          dekin[j][k] += vcm[k]*(tm*hvcm[j]-mv[j]);
      }
  }
  pr_rvecs(log,0,"dekin",dekin,DIM);
  pr_rvecs(log,0," ekin", ekin,DIM);
  fprintf(log,"dekin = %g, ekin = %g  vcm = (%8.4f %8.4f %8.4f)\n",
          trace(dekin),trace(ekin),vcm[XX],vcm[YY],vcm[ZZ]);
  fprintf(log,"mv = (%8.4f %8.4f %8.4f)\n",
          mv[XX],mv[YY],mv[ZZ]);
}
