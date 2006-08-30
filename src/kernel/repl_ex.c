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
#include "mvdata.h"
#include "domdec.h"

gmx_repl_ex_t *init_replica_exchange(FILE *fplog,
				     const gmx_multisim_t *ms,
				     const t_state *state,
				     const t_inputrec *ir,
				     int nst,int init_seed)
{
  real temp,pres;
  int  i,j,k;
  gmx_repl_ex_t *re;

  fprintf(fplog,"\nInitializing Replica Exchange\n");
  please_cite(fplog,"Hukushima96a");
  if (ir->epc != epcNO) {
    fprintf(fplog,"Repl  Using Constant Pressure REMD.\n");
    please_cite(fplog,"Okabe2001a");
  }

  if (ms == NULL || ms->nsim == 1)
    gmx_fatal(FARGS,"Nothing to exchange with only one replica, maybe you forgot to set the -multi option of mdrun?");

  temp = ir->opts.ref_t[0];
  for(i=1; (i<ir->opts.ngtc); i++) {
    if (ir->opts.ref_t[i] != temp) {
      fprintf(fplog,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
      fprintf(stderr,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
    }
  }
  
  check_multi_int(fplog,ms,state->natoms,"the number of atoms");
  check_multi_int(fplog,ms,ir->eI,"the integrator");
  check_multi_int(fplog,ms,ir->nsteps,"nsteps");
  check_multi_int(fplog,ms,ir->init_step,"init_step");
  check_multi_int(fplog,ms,ir->etc,"the temperature coupling");
  check_multi_int(fplog,ms,ir->opts.ngtc,
		  "the number of temperature coupling groups");
  check_multi_int(fplog,ms,ir->epc,"the pressure coupling");

  snew(re,1);

  re->repl     = ms->sim;
  re->nrepl    = ms->nsim;
  re->bNPT     = (ir->epc != epcNO);
  
  fprintf(fplog,"Repl  There are %d replicas:\n",re->nrepl);

  snew(re->temp,re->nrepl);
  re->temp[re->repl] = temp;
  gmx_sum_sim(re->nrepl,re->temp,ms);
  if (re->bNPT) {
    snew(re->pres,re->nrepl);
    if (ir->epct == epctSURFACETENSION) {
      pres = ir->ref_p[ZZ][ZZ];
    } else {
      pres = 0;
      j = 0; 
      for(i=0; i<DIM; i++)
	if (ir->compress[i][i] != 0) {
	  pres += ir->ref_p[i][i];
	  j++;
	}
      pres /= j;
    }
    re->pres[re->repl] = pres;
    gmx_sum_sim(re->nrepl,re->pres,ms);  }
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
  fprintf(fplog,"Repl   ");
  for(i=0; i<re->nrepl; i++)
    fprintf(fplog," %3d  ",re->ind[i]);
  fprintf(fplog,"\nRepl  T");
  for(i=0; i<re->nrepl; i++)
    fprintf(fplog," %5.1f",re->temp[re->ind[i]]);
  if (re->bNPT) {
    fprintf(fplog,"\nRepl  p");
    for(i=0; i<re->nrepl; i++)
    {
      fprintf(fplog," %5.2f",re->pres[re->ind[i]]);
    }

    for(i=0; i<re->nrepl; i++)
    {
      if ((i > 0) && (re->pres[re->ind[i]] < re->pres[re->ind[i-1]]))
      {
        gmx_fatal(FARGS,"The reference pressure decreases with increasing temperature");
      }
    }
  }
  fprintf(fplog,"\nRepl  ");
  
  re->nst = nst;
  if (init_seed == -1) {
    if (MASTERSIM(ms))
      re->seed = make_seed();
    else
      re->seed = 0;
    gmx_sumi_sim(1,&(re->seed),ms);
  } else {
    re->seed = init_seed;
  }
  fprintf(fplog,"\nRepl  exchange interval: %d\n",re->nst);
  fprintf(fplog,"\nRepl  random seed: %d\n",re->seed);

  re->nattempt[0] = 0;
  re->nattempt[1] = 0;
  snew(re->prob_sum,re->nrepl);
  snew(re->nexchange,re->nrepl);

  fprintf(fplog,
	  "Repl  below: x=exchange, pr=probability\n");

  return re;
}

static void exchange_reals(const gmx_multisim_t *ms,int b,real *v,int n)
{
  real *buf;
  int  i;

  if (v) {
    snew(buf,n);
#ifdef GMX_MPI
    /*
    MPI_Sendrecv(v,  n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
		 buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
		 ms->mpi_comm_masters,MPI_STATUS_IGNORE);
    */
    {
      MPI_Request mpi_req;

      MPI_Isend(v,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
		ms->mpi_comm_masters,&mpi_req);
      MPI_Recv(buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
	       ms->mpi_comm_masters,MPI_STATUS_IGNORE);
      MPI_Wait(&mpi_req,MPI_STATUS_IGNORE);
    }
#endif
    for(i=0; i<n; i++)
      v[i] = buf[i];
    sfree(buf);
  }
}

static void exchange_rvecs(const gmx_multisim_t *ms,int b,rvec *v,int n)
{
  rvec *buf;
  int  i;
  
  if (v) {
    snew(buf,n);
#ifdef GMX_MPI
    /*
    MPI_Sendrecv(v[0],  n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
		 buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
		 ms->mpi_comm_masters,MPI_STATUS_IGNORE);
    */
    {
      MPI_Request mpi_req;

      MPI_Isend(v[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
		ms->mpi_comm_masters,&mpi_req);
      MPI_Recv(buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
	       ms->mpi_comm_masters,MPI_STATUS_IGNORE);
      MPI_Wait(&mpi_req,MPI_STATUS_IGNORE);
    }
#endif
    for(i=0; i<n; i++)
      copy_rvec(buf[i],v[i]);
    sfree(buf);
  }
}

static void exchange_state(const gmx_multisim_t *ms,int b,t_state *state)
{
  /* When t_state changes, this code should be updated. */
  exchange_rvecs(ms,b,state->box,DIM);
  exchange_rvecs(ms,b,state->boxv,DIM);
  exchange_rvecs(ms,b,state->pcoupl_mu,DIM);
  exchange_reals(ms,b,state->nosehoover_xi,state->ngtc);
  exchange_rvecs(ms,b,state->x,state->natoms);
  exchange_rvecs(ms,b,state->v,state->natoms);
  exchange_rvecs(ms,b,state->sd_X,state->natoms);
}

static void scale_velocities(t_state *state,real fac)
{
  int i;

  if (state->v)
    for(i=0; i<state->natoms; i++)
      svmul(fac,state->v[i],state->v[i]);
}

static void broadcast_reals(const t_commrec *cr,int n,real *v)
{
#ifdef GMX_MPI
  if (v)
    MPI_Bcast(v,n*sizeof(real),MPI_BYTE,MASTERRANK(cr),
	      cr->mpi_comm_mysim);
#endif
}

static void broadcast_rvecs(const t_commrec *cr,int n,rvec *v)
{
#ifdef GMX_MPI
  if (v)
    MPI_Bcast(v[0],n*sizeof(rvec),MPI_BYTE,MASTERRANK(cr),
	      cr->mpi_comm_mysim);
#endif
}

static void broadcast_matrix(const t_commrec *cr,matrix m)
{
#ifdef GMX_MPI
  MPI_Bcast(m[0],sizeof(matrix),MPI_BYTE,MASTERRANK(cr),
	    cr->mpi_comm_mysim);
#endif
}

static void pd_collect_state(const t_commrec *cr,t_nsborder *nsb,
			     t_state *state)
{
  int shift;
  
  if (debug)
    fprintf(debug,"Collecting state before exchange\n");
  
  shift = cr->nnodes - cr->npmenodes - 1;
  move_rvecs(cr,FALSE,FALSE,cr->left,cr->right,
	     state->x,NULL,shift,nsb,NULL);
  if (state->v)
    move_rvecs(cr,FALSE,FALSE,cr->left,cr->right,
	       state->v,NULL,shift,nsb,NULL);
  if (state->sd_X)
    move_rvecs(cr,FALSE,FALSE,cr->left,cr->right,
	       state->sd_X,NULL,shift,nsb,NULL);
}

void pd_distribute_state(const t_commrec *cr,t_state *state)
{
  /* This code does seem to distribute the state,
     but also seems to cause memory problems,
     as after move_cgcm the cg's on the non-master nodes are corrupted
  if (MASTER(cr))
    mv_state(cr,cr->right,state);
  else
    ld_state(cr,cr->left,state);
  if (cr->nodeid != cr->nnodes-1)
    mv_state(cr,cr->right,state);
  */

  /* All array dimensions are already known on all nodes */
  broadcast_matrix(cr,state->box);
  broadcast_matrix(cr,state->box);
  broadcast_matrix(cr,state->pcoupl_mu);
  broadcast_reals(cr,state->ngtc,state->nosehoover_xi);
  broadcast_rvecs(cr,state->natoms,state->x);
  broadcast_rvecs(cr,state->natoms,state->v);
  broadcast_rvecs(cr,state->natoms,state->sd_X);
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

static int get_replica_exchange(FILE *fplog,const gmx_multisim_t *ms,
				gmx_repl_ex_t *re,real epot,real vol,
				int step,real time)
{
  int  m,i,a,b;
  real *Epot,*prob,ediff,delta,dpV,*Vol,betaA,betaB;
  bool *bEx,bPrint;
  int  exchange;

  fprintf(fplog,"Replica exchange at step %d time %g\n",step,time);
  
  snew(Epot,re->nrepl);
  snew(Vol,re->nrepl);
  Epot[re->repl] = epot;
  Vol[re->repl]  = vol;
  gmx_sum_sim(re->nrepl,Epot,ms);
  gmx_sum_sim(re->nrepl,Vol,ms);

  snew(bEx,re->nrepl);
  snew(prob,re->nrepl);

  exchange = -1;
  m = (step / re->nst) % 2;
  for(i=1; i<re->nrepl; i++) {
    a = re->ind[i-1];
    b = re->ind[i];
    bPrint = (re->repl==a || re->repl==b);
    if (i % 2 == m) {
      /* Use equations from:
       * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
       */
      ediff = Epot[b] - Epot[a];
      betaA = 1.0/(re->temp[a]*BOLTZ);
      betaB = 1.0/(re->temp[b]*BOLTZ);
      delta = (betaA-betaB)*ediff;
      if (bPrint)
	fprintf(fplog,"Repl %d <-> %d  dE = %10.3e",a,b,delta);
      if (re->bNPT) {
	dpV = (betaA*re->pres[a]-betaB*re->pres[b])*(Vol[b]-Vol[a])/PRESFAC;
	if (bPrint)
	  fprintf(fplog,"  dpV = %10.3e  d = %10.3e",dpV,delta + dpV);
	delta += dpV;
      }
      if (bPrint)
	fprintf(fplog,"\n");
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
	if (a == re->repl) {
	  exchange = b;
	} else if (b == re->repl) {
	  exchange = a;
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
  sfree(Vol);
  
  re->nattempt[m]++;

  return exchange;
}

static void write_debug_x(t_state *state)
{
  int i;

  if (debug) {
    for(i=0; i<state->natoms; i+=10)
      fprintf(debug,"dx %5d %10.5f %10.5f %10.5f\n",i,state->x[i][XX],state->x[i][YY],state->x[i][ZZ]);
  }
}

bool replica_exchange(FILE *fplog,const t_commrec *cr,gmx_repl_ex_t *re,
		      t_state *state,real epot,t_nsborder *nsb,
		      t_block *cgs,t_state *state_local,
		      int step,real time)
{
  gmx_multisim_t *ms;
  int  exchange=-1,shift;
  bool bExchanged=FALSE;

  ms = cr->ms;

  if (MASTER(cr)) {
    exchange = get_replica_exchange(fplog,ms,re,epot,det(state->box),
				    step,time);
    bExchanged = (exchange >= 0);
  }
      
  if (PAR(cr)) {
#ifdef GMX_MPI
    MPI_Bcast(&bExchanged,sizeof(bool),MPI_BYTE,MASTERRANK(cr),
	      cr->mpi_comm_mysim);
#endif
  }
  
  if (bExchanged) {
    if (PAR(cr)) {
      if (DOMAINDECOMP(cr))
	dd_collect_state(cr->dd,cgs,state_local,state);
      else
	pd_collect_state(cr,nsb,state);
    }
    if (MASTER(cr)) {
      if (debug)
	fprintf(debug,"Exchanging %d with %d\n",ms->sim,exchange);
      exchange_state(ms,exchange,state);
      scale_velocities(state,sqrt(re->temp[ms->sim]/re->temp[exchange]));
    }
  }
  
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
