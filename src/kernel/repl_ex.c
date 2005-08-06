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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "repl_ex.h"
#include "network.h"
#include "random.h"
#include "smalloc.h"
#include "physics.h"
#include "copyrite.h"
#include "macros.h"
#include "vec.h"
#include "names.h"

gmx_repl_ex_t *init_replica_exchange(FILE *fplog,
				     const t_commrec *mcr,
				     const t_state *state,
				     const t_inputrec *ir,
				     int nst,int init_seed,
				     bool bNPT)
{
  real temp;
  int  i,j,k;
  gmx_repl_ex_t *re;

  fprintf(fplog,"\nInitializing Replica Exchange\n");
  please_cite(fplog,"Hukushima96a");
  if (bNPT) {
    fprintf(fplog,"Repl  Using Constant Pressure REMD.\n");
    please_cite(fplog,"Okabe2001a");
  }

  if (mcr == NULL || mcr->nnodes == 1)
    gmx_fatal(FARGS,"Nothing to exchange with only one replica");

  temp = ir->opts.ref_t[0];
  for(i=1; (i<ir->opts.ngtc); i++) {
    if (ir->opts.ref_t[i] != temp) {
      fprintf(fplog,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
      fprintf(stderr,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
    }
  }
  
  check_multi_int(fplog,mcr,state->natoms,"the number of atoms");
  check_multi_int(fplog,mcr,ir->eI,"the integrator");
  check_multi_int(fplog,mcr,ir->nsteps,"nsteps");
  check_multi_int(fplog,mcr,ir->init_step,"init_step");
  check_multi_int(fplog,mcr,ir->opts.ngtc,
		  "the number of temperature coupling groups");
  

  snew(re,1);

  re->repl  = mcr->nodeid;
  re->nrepl = mcr->nnodes;
  re->bNPT  = bNPT;
  
  fprintf(fplog,"Repl  There are %d replicas:\n",re->nrepl);

  snew(re->temp,re->nrepl);
  re->temp[re->repl] = temp;
  gmx_sum(re->nrepl,re->temp,mcr);
  snew(re->ind,re->nrepl);
  /* Make an index for increasing temperature order */
  for(i=0; i<re->nrepl; i++)
    re->ind[i] = i;
  for(i=0; i<re->nrepl; i++) {
    for(j=i+1; j<re->nrepl; j++) {
      if (re->temp[re->ind[j]] < re->temp[re->ind[i]]) {
	k = re->ind[i];
	re->ind[i] = re->ind[j];
	re->ind[j] = k;
      } else if (re->temp[re->ind[j]] == re->temp[re->ind[i]]) {
	gmx_fatal(FARGS,"Two replicas have identical temperatures");
      }
    }
  }
  fprintf(fplog,"Repl ");
  for(i=0; i<re->nrepl; i++)
    fprintf(fplog," %3d  ",re->ind[i]);
  fprintf(fplog,"\nRepl ");
  for(i=0; i<re->nrepl; i++)
    fprintf(fplog," %5.1f",re->temp[re->ind[i]]);
  fprintf(fplog,"\nRepl  ");

  re->nst = nst;
  if (init_seed == -1) {
    if (MASTER(mcr))
      re->seed = make_seed();
    else
      re->seed = 0;
    gmx_sumi(1,&(re->seed),mcr);
  } else {
    re->seed = init_seed;
  }
  fprintf(fplog,"\nRepl  random seed: %d\n",re->seed);

  re->nattempt[0] = 0;
  re->nattempt[1] = 0;
  snew(re->prob_sum,re->nrepl);
  snew(re->nexchange,re->nrepl);

  fprintf(fplog,
	  "Repl  below: x=exchange, pr=probability\n");

  return re;
}

static void exchange_reals(const t_commrec *mcr,int a,int b,real *v,int n)
{
  real *buf;
  int  i;
  
  if (v) {
    if (mcr->nodeid == a) {
      snew(buf,n);
      gmx_rxs(b,buf,n*sizeof(real));
      gmx_txs(b,v,n*sizeof(real));
      for(i=0; i<n; i++)
	v[i] = buf[i];
      sfree(buf);
    } else if (mcr->nodeid == b) {
      gmx_txs(a,v,n*sizeof(real));
      gmx_rxs(a,v,n*sizeof(real));
    }
  }
}

static void exchange_rvecs(const t_commrec *mcr,int a,int b,rvec *v,int n)
{
  rvec *buf;
  int  i;
  
  if (v) {
    if (mcr->nodeid == a) {
      snew(buf,n);
      gmx_rxs(b,buf[0],n*sizeof(rvec));
      gmx_txs(b,v[0],n*sizeof(rvec));
      for(i=0; i<n; i++)
	copy_rvec(buf[i],v[i]);
      sfree(buf);
    } else if (mcr->nodeid == b) {
      gmx_txs(a,v[0],n*sizeof(rvec));
      gmx_rxs(a,v[0],n*sizeof(rvec));
    }
  }
}

static void exchange_state(const t_commrec *mcr,int a,int b,t_state *state)
{
  /* When t_state changes, this code should be updated. */
  exchange_rvecs(mcr,a,b,state->box,DIM);
  exchange_rvecs(mcr,a,b,state->boxv,DIM);
  exchange_rvecs(mcr,a,b,state->pcoupl_mu,DIM);
  exchange_reals(mcr,a,b,state->nosehoover_xi,state->ngtc);
  exchange_rvecs(mcr,a,b,state->x,state->natoms);
  exchange_rvecs(mcr,a,b,state->v,state->natoms);
  exchange_rvecs(mcr,a,b,state->sd_X,state->natoms);
}

static void scale_velocities(t_state *state,real fac)
{
  int i;

  if (state->v)
    for(i=0; i<state->natoms; i++)
      svmul(fac,state->v[i],state->v[i]);
}

static void print_ind(FILE *fplog,char *leg,int n,int *ind,bool *bEx)
{
  int i;

  fprintf(fplog,"Repl %2s %2d",leg,ind[0]);
  for(i=1; i<n; i++) {
    fprintf(fplog," %c %2d",(bEx!=0 && bEx[i]) ? 'x' : ' ',ind[i]);
  }
  fprintf(fplog,"\n");
}

static void print_prob(FILE *fplog,char *leg,int n,real *prob)
{
  int  i;
  char buf[8];
  
  fprintf(fplog,"Repl %2s ",leg);
  for(i=1; i<n; i++) {
    if (prob[i] >= 0) {
      sprintf(buf,"%4.2f",prob[i]);
      fprintf(fplog,"  %3s",buf[0]=='1' ? "1.0" : buf+1);
    } else {
      fprintf(fplog,"     ");
    }
  }
  fprintf(fplog,"\n");
}

static void print_count(FILE *fplog,char *leg,int n,int *count)
{
  int i;

  fprintf(fplog,"Repl %2s ",leg);
  for(i=1; i<n; i++) {
    fprintf(fplog," %4d",count[i]);
  }
  fprintf(fplog,"\n");
}

bool replica_exchange(FILE *fplog,const t_commrec *mcr,gmx_repl_ex_t *re,
		      t_state *state,real epot,int step,real time,
		      real pres)
{
  int  m,i,a,b;
  real *Epot,*prob,ediff,delta,ddd,*Pres,*Vol,betaA,betaB;
  bool *bEx,bExchanged;

  fprintf(fplog,"Replica exchange at step %d time %g\n",step,time);
  snew(Epot,re->nrepl);
  snew(Pres,re->nrepl);
  snew(Vol,re->nrepl);
  Epot[re->repl] = epot;
  Pres[re->repl] = pres;
  Vol[re->repl]  = det(state->box);
  gmx_sum(re->nrepl,Epot,mcr);
  gmx_sum(re->nrepl,Vol,mcr);
  gmx_sum(re->nrepl,Pres,mcr);

  snew(bEx,re->nrepl);
  snew(prob,re->nrepl);

  bExchanged = FALSE;
  m = (step / re->nst) % 2;
  for(i=1; i<re->nrepl; i++) {
    a = re->ind[i-1];
    b = re->ind[i];
    if (i % 2 == m) {
      /* Use equations from:
       * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
       */
      ediff = Epot[b] - Epot[a];
      betaA = 1.0/(re->temp[a]*BOLTZ);
      betaB = 1.0/(re->temp[b]*BOLTZ);
      delta = (betaA-betaB)*ediff;
      if (re->bNPT) {
	ddd = (betaA*Pres[a]-betaB*Pres[b])*(Vol[b]-Vol[a]);
        fprintf(fplog,"i = %d deltaE = %12.5e deltaP = %12.5e delta = %12.5e\n",
	        i,delta,ddd,delta+ddd);
	delta += ddd;
      }
      if (delta <= 0) {
	prob[i] = 1;
	bEx[i] = TRUE;
      } else {
        if (delta > 100)
          prob[i] = 0;
        else
	  prob[i] = exp(-delta);
	bEx[i] = (rando(&(re->seed)) < prob[i]);
      }
      re->prob_sum[i] += prob[i];    
      if (bEx[i]) {
	exchange_state(mcr,a,b,state);
	if (a == re->repl) {
	  scale_velocities(state,sqrt(re->temp[a]/re->temp[b]));
	  bExchanged = TRUE;
	} else if (b == re->repl) {
	  scale_velocities(state,sqrt(re->temp[b]/re->temp[a]));
	  bExchanged = TRUE;
	}
	re->nexchange[i]++;
      }
    } else {
      prob[i] = -1;
      bEx[i] = FALSE;
    }
  }
  print_ind(fplog,"ex",re->nrepl,re->ind,bEx);
  print_prob(fplog,"pr",re->nrepl,prob);
  fprintf(fplog,"\n");

  sfree(bEx);
  sfree(prob);
  sfree(Epot);
  sfree(Pres);
  sfree(Vol);
  
  re->nattempt[m]++;

  return bExchanged;
}

void print_replica_exchange_statistics(FILE *fplog,gmx_repl_ex_t *re)
{
  real *prob;
  int  i;
  
  fprintf(fplog,"\nReplica exchange statistics\n");
  fprintf(fplog,"Repl  %d attempts, %d odd, %d even\n",
	  re->nattempt[0]+re->nattempt[1],re->nattempt[1],re->nattempt[0]);

  snew(prob,re->nrepl);
  
  fprintf(fplog,"Repl  average probabilities:\n");
  for(i=1; i<re->nrepl; i++) {
    if (re->nattempt[i%2] == 0)
      prob[i] = 0;
    else
      prob[i] =  re->prob_sum[i]/re->nattempt[i%2];
  }
  print_ind(fplog,"",re->nrepl,re->ind,NULL);
  print_prob(fplog,"",re->nrepl,prob);

  fprintf(fplog,"Repl  number of exchanges:\n");
  print_ind(fplog,"",re->nrepl,re->ind,NULL);
  print_count(fplog,"",re->nrepl,re->nexchange);
  
  fprintf(fplog,"Repl  average number of exchanges:\n");
  for(i=1; i<re->nrepl; i++) {
    if (re->nattempt[i%2] == 0)
      prob[i] = 0;
    else
      prob[i] =  ((real)re->nexchange[i])/re->nattempt[i%2];
  }
  print_ind(fplog,"",re->nrepl,re->ind,NULL);
  print_prob(fplog,"",re->nrepl,prob);

  sfree(prob);
  
  fprintf(fplog,"\n");
}
