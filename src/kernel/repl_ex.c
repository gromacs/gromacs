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
#include "partdec.h"

typedef struct gmx_repl_ex
{
    int  repl;
    int  nrepl;
    real temp;
    int  type;
    real *q;
    gmx_bool bNPT;
    real *pres;
    int  *ind;
    int  nst;
    int  seed;
    int  nattempt[2];
    real *prob_sum;
    int  *nexchange;
} t_gmx_repl_ex;

enum { ereTEMP, ereLAMBDA, ereNR };
const char *erename[ereNR] = { "temperature", "lambda" };

static void repl_quantity(FILE *fplog,const gmx_multisim_t *ms,
                          struct gmx_repl_ex *re,int ere,real q)
{
    real *qall;
    gmx_bool bDiff;
    int  s;

    snew(qall,ms->nsim);
    qall[re->repl] = q;
    gmx_sum_sim(ms->nsim,qall,ms);

    bDiff = FALSE;
    for(s=1; s<ms->nsim; s++)
    {
        if (qall[s] != qall[0])
        {
            bDiff = TRUE;
        }
    }
    if (bDiff)
    {
        if (re->type >= 0 && re->type < ereNR)
        {
            gmx_fatal(FARGS,"For replica exchange both %s and %s differ",
                      erename[re->type],erename[ere]);
        }
        /* Set the replica exchange type and quantities */
        re->type = ere;
        snew(re->q,re->nrepl);
        for(s=0; s<ms->nsim; s++)
        {
            re->q[s] = qall[s];
        }
    }

    sfree(qall);
}

gmx_repl_ex_t init_replica_exchange(FILE *fplog,
                                    const gmx_multisim_t *ms,
                                    const t_state *state,
                                    const t_inputrec *ir,
                                    int nst,int init_seed)
{
    real temp,pres;
    int  i,j,k;
    struct gmx_repl_ex *re;

    fprintf(fplog,"\nInitializing Replica Exchange\n");

    if (ms == NULL || ms->nsim == 1)
    {
        gmx_fatal(FARGS,"Nothing to exchange with only one replica, maybe you forgot to set the -multi option of mdrun?");
    }

    snew(re,1);

    re->repl     = ms->sim;
    re->nrepl    = ms->nsim;

    fprintf(fplog,"Repl  There are %d replicas:\n",re->nrepl);

    check_multi_int(fplog,ms,state->natoms,"the number of atoms");
    check_multi_int(fplog,ms,ir->eI,"the integrator");
    check_multi_large_int(fplog,ms,ir->init_step+ir->nsteps,"init_step+nsteps");
    check_multi_large_int(fplog,ms,(ir->init_step+nst-1)/nst,
                          "first exchange step: init_step/-replex");
    check_multi_int(fplog,ms,ir->etc,"the temperature coupling");
    check_multi_int(fplog,ms,ir->opts.ngtc,
                    "the number of temperature coupling groups");
    check_multi_int(fplog,ms,ir->epc,"the pressure coupling");
    check_multi_int(fplog,ms,ir->efep,"free energy");

    re->temp = ir->opts.ref_t[0];
    for(i=1; (i<ir->opts.ngtc); i++)
    {
        if (ir->opts.ref_t[i] != re->temp)
        {
            fprintf(fplog,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
            fprintf(stderr,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
        }
    }

    re->type = -1;
    for(i=0; i<ereNR; i++)
    {
        switch (i)
        {
        case ereTEMP:
            repl_quantity(fplog,ms,re,i,re->temp);
            break;
        case ereLAMBDA:
            if (ir->efep != efepNO)
            {
                repl_quantity(fplog,ms,re,i,ir->init_lambda);
            }
            break;
        default:
            gmx_incons("Unknown replica exchange quantity");
        }
    }
    if (re->type == -1)
    {
        gmx_fatal(FARGS,"The properties of the %d systems are all the same, there is nothing to exchange",re->nrepl);
    }

    switch (re->type)
    {
    case ereTEMP:
        please_cite(fplog,"Sugita1999a");
        if (ir->epc != epcNO)
        {
            re->bNPT = TRUE;
            fprintf(fplog,"Repl  Using Constant Pressure REMD.\n");
            please_cite(fplog,"Okabe2001a");
        }
        if (ir->etc == etcBERENDSEN)
        {
            gmx_fatal(FARGS,"REMD with the %s thermostat does not produce correct potential energy distributions, consider using the %s thermostat instead",
                      ETCOUPLTYPE(ir->etc),ETCOUPLTYPE(etcVRESCALE));
        }
        break;
    case ereLAMBDA:
        if (ir->delta_lambda != 0)
        {
            gmx_fatal(FARGS,"delta_lambda is not zero");
        }
        break;
    }

    if (re->bNPT)
    {
        snew(re->pres,re->nrepl);
        if (ir->epct == epctSURFACETENSION)
        {
            pres = ir->ref_p[ZZ][ZZ];
        }
        else
        {
            pres = 0;
            j = 0;
            for(i=0; i<DIM; i++)
            {
                if (ir->compress[i][i] != 0)
                {
                    pres += ir->ref_p[i][i];
                    j++;
                }
            }
            pres /= j;
        }
        re->pres[re->repl] = pres;
        gmx_sum_sim(re->nrepl,re->pres,ms);
    }

    snew(re->ind,re->nrepl);
    /* Make an index for increasing temperature order */
    for(i=0; i<re->nrepl; i++)
    {
        re->ind[i] = i;
    }
    for(i=0; i<re->nrepl; i++)
    {
        for(j=i+1; j<re->nrepl; j++)
        {
            if (re->q[re->ind[j]] < re->q[re->ind[i]])
            {
                k = re->ind[i];
                re->ind[i] = re->ind[j];
                re->ind[j] = k;
            }
            else if (re->q[re->ind[j]] == re->q[re->ind[i]])
            {
                gmx_fatal(FARGS,"Two replicas have identical %ss",erename[re->type]);
            }
        }
    }
    fprintf(fplog,"Repl   ");
    for(i=0; i<re->nrepl; i++)
    {
        fprintf(fplog," %3d  ",re->ind[i]);
    }
    switch (re->type)
    {
    case ereTEMP:
        fprintf(fplog,"\nRepl  T");
        for(i=0; i<re->nrepl; i++)
        {
            fprintf(fplog," %5.1f",re->q[re->ind[i]]);
        }
        break;
    case ereLAMBDA:
        fprintf(fplog,"\nRepl  l");
        for(i=0; i<re->nrepl; i++)
        {
            fprintf(fplog," %5.3f",re->q[re->ind[i]]);
        }
        break;
    default:
        gmx_incons("Unknown replica exchange quantity");
    }
    if (re->bNPT)
    {
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
    if (init_seed == -1)
    {
        if (MASTERSIM(ms))
        {
            re->seed = make_seed();
        }
        else
        {
            re->seed = 0;
        }
        gmx_sumi_sim(1,&(re->seed),ms);
    }
    else
    {
        re->seed = init_seed;
    }
    fprintf(fplog,"\nRepl  exchange interval: %d\n",re->nst);
    fprintf(fplog,"\nRepl  random seed: %d\n",re->seed);

    re->nattempt[0] = 0;
    re->nattempt[1] = 0;
    snew(re->prob_sum,re->nrepl);
    snew(re->nexchange,re->nrepl);

    fprintf(fplog,"Repl  below: x=exchange, pr=probability\n");

    return re;
}

static void exchange_reals(const gmx_multisim_t *ms,int b,real *v,int n)
{
    real *buf;
    int  i;

    if (v)
    {
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
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_doubles(const gmx_multisim_t *ms,int b,double *v,int n)
{
    double *buf;
    int  i;

    if (v)
    {
        snew(buf,n);
#ifdef GMX_MPI
        /*
          MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
          buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
          ms->mpi_comm_masters,MPI_STATUS_IGNORE);
        */
        {
            MPI_Request mpi_req;

            MPI_Isend(v,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
                      ms->mpi_comm_masters,&mpi_req);
            MPI_Recv(buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
                     ms->mpi_comm_masters,MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req,MPI_STATUS_IGNORE);
        }
#endif
        for(i=0; i<n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_rvecs(const gmx_multisim_t *ms,int b,rvec *v,int n)
{
    rvec *buf;
    int  i;
  
    if (v)
    {
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
        {
            copy_rvec(buf[i],v[i]);
        }
        sfree(buf);
    }
}

static void exchange_state(const gmx_multisim_t *ms,int b,t_state *state)
{
    /* When t_state changes, this code should be updated. */
    int ngtc,nnhpres;
    ngtc = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    exchange_rvecs(ms,b,state->box,DIM);
    exchange_rvecs(ms,b,state->box_rel,DIM);
    exchange_rvecs(ms,b,state->boxv,DIM);
    exchange_reals(ms,b,&(state->veta),1);
    exchange_reals(ms,b,&(state->vol0),1);
    exchange_rvecs(ms,b,state->svir_prev,DIM);
    exchange_rvecs(ms,b,state->fvir_prev,DIM);
    exchange_rvecs(ms,b,state->pres_prev,DIM);
    exchange_doubles(ms,b,state->nosehoover_xi,ngtc);
    exchange_doubles(ms,b,state->nosehoover_vxi,ngtc);
    exchange_doubles(ms,b,state->nhpres_xi,nnhpres);
    exchange_doubles(ms,b,state->nhpres_vxi,nnhpres);
    exchange_doubles(ms,b,state->therm_integral,state->ngtc);
    exchange_rvecs(ms,b,state->x,state->natoms);
    exchange_rvecs(ms,b,state->v,state->natoms);
    exchange_rvecs(ms,b,state->sd_X,state->natoms);
}

static void copy_rvecs(rvec *s,rvec *d,int n)
{
    int i;

    if (d != NULL)
    {
        for(i=0; i<n; i++)
        {
            copy_rvec(s[i],d[i]);
        }
    }
}

static void copy_doubles(const double *s,double *d,int n)
{
    int i;

    if (d != NULL)
    {
        for(i=0; i<n; i++)
        {
            d[i] = s[i];
        }
    }
}

#define scopy_rvecs(v,n)   copy_rvecs(state->v,state_local->v,n);
#define scopy_doubles(v,n) copy_doubles(state->v,state_local->v,n);

static void copy_state_nonatomdata(t_state *state,t_state *state_local)
{
    /* When t_state changes, this code should be updated. */
    int ngtc,nnhpres;
    ngtc = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    scopy_rvecs(box,DIM);
    scopy_rvecs(box_rel,DIM);
    scopy_rvecs(boxv,DIM);
    state_local->veta = state->veta;
    state_local->vol0 = state->vol0;
    scopy_rvecs(svir_prev,DIM);
    scopy_rvecs(fvir_prev,DIM);
    scopy_rvecs(pres_prev,DIM);
    scopy_doubles(nosehoover_xi,ngtc);
    scopy_doubles(nosehoover_vxi,ngtc);
    scopy_doubles(nhpres_xi,nnhpres);
    scopy_doubles(nhpres_vxi,nnhpres);
    scopy_doubles(therm_integral,state->ngtc);
    scopy_rvecs(x,state->natoms);
    scopy_rvecs(v,state->natoms);
    scopy_rvecs(sd_X,state->natoms);
}

static void scale_velocities(t_state *state,real fac)
{
    int i;

    if (state->v)
    {
        for(i=0; i<state->natoms; i++)
        {
            svmul(fac,state->v[i],state->v[i]);
        }
    }
}

static void pd_collect_state(const t_commrec *cr,t_state *state)
{
    int shift;
  
    if (debug)
    {
        fprintf(debug,"Collecting state before exchange\n");
    }
    shift = cr->nnodes - cr->npmenodes - 1;
    move_rvecs(cr,FALSE,FALSE,GMX_LEFT,GMX_RIGHT,state->x,NULL,shift,NULL);
    if (state->v)
    {
        move_rvecs(cr,FALSE,FALSE,GMX_LEFT,GMX_RIGHT,state->v,NULL,shift,NULL);
    }
    if (state->sd_X)
    {
        move_rvecs(cr,FALSE,FALSE,GMX_LEFT,GMX_RIGHT,state->sd_X,NULL,shift,NULL);
    }
}

static void print_ind(FILE *fplog,const char *leg,int n,int *ind,gmx_bool *bEx)
{
    int i;

    fprintf(fplog,"Repl %2s %2d",leg,ind[0]);
    for(i=1; i<n; i++)
    {
        fprintf(fplog," %c %2d",(bEx!=0 && bEx[i]) ? 'x' : ' ',ind[i]);
    }
    fprintf(fplog,"\n");
}

static void print_prob(FILE *fplog,const char *leg,int n,real *prob)
{
    int  i;
    char buf[8];
  
    fprintf(fplog,"Repl %2s ",leg);
    for(i=1; i<n; i++)
    {
        if (prob[i] >= 0)
        {
            sprintf(buf,"%4.2f",prob[i]);
            fprintf(fplog,"  %3s",buf[0]=='1' ? "1.0" : buf+1);
        }
        else
        {
            fprintf(fplog,"     ");
        }
    }
    fprintf(fplog,"\n");
}

static void print_count(FILE *fplog,const char *leg,int n,int *count)
{
    int i;

    fprintf(fplog,"Repl %2s ",leg);
    for(i=1; i<n; i++)
    {
        fprintf(fplog," %4d",count[i]);
    }
    fprintf(fplog,"\n");
}

static int get_replica_exchange(FILE *fplog,const gmx_multisim_t *ms,
                                struct gmx_repl_ex *re,real *ener,real vol,
                                gmx_large_int_t step,real time)
{
    int  m,i,a,b;
    real *Epot=NULL,*Vol=NULL,*dvdl=NULL,*prob;
    real ediff=0,delta=0,dpV=0,betaA=0,betaB=0;
    gmx_bool *bEx,bPrint;
    int  exchange;

    fprintf(fplog,"Replica exchange at step " gmx_large_int_pfmt " time %g\n",step,time);
  
    switch (re->type)
    {
    case ereTEMP:
        snew(Epot,re->nrepl);
        snew(Vol,re->nrepl);
        Epot[re->repl] = ener[F_EPOT];
        Vol[re->repl]  = vol;
        gmx_sum_sim(re->nrepl,Epot,ms);
        gmx_sum_sim(re->nrepl,Vol,ms);
        break;
    case ereLAMBDA:
        snew(dvdl,re->nrepl);
        dvdl[re->repl] = ener[F_DVDL];
        gmx_sum_sim(re->nrepl,dvdl,ms);
        break;
    }

    snew(bEx,re->nrepl);
    snew(prob,re->nrepl);

    exchange = -1;
    m = (step / re->nst) % 2;
    for(i=1; i<re->nrepl; i++)
    {
        a = re->ind[i-1];
        b = re->ind[i];
        bPrint = (re->repl==a || re->repl==b);
        if (i % 2 == m)
        {
            switch (re->type)
            {
            case ereTEMP:
                /* Use equations from:
                 * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
                 */
                ediff = Epot[b] - Epot[a];
                betaA = 1.0/(re->q[a]*BOLTZ);
                betaB = 1.0/(re->q[b]*BOLTZ);
                delta = (betaA - betaB)*ediff;
                break;
            case ereLAMBDA:
                /* Here we exchange based on a linear extrapolation of dV/dlambda.
                 * We would like to have the real energies
                 * from foreign lambda calculations.
                 */
                ediff = (dvdl[a] - dvdl[b])*(re->q[b] - re->q[a]);
                delta = ediff/(BOLTZ*re->temp);
                break;
            default:
                gmx_incons("Unknown replica exchange quantity");
            }
            if (bPrint)
            {
                fprintf(fplog,"Repl %d <-> %d  dE = %10.3e",a,b,delta);
            }
            if (re->bNPT)
            {
                dpV = (betaA*re->pres[a]-betaB*re->pres[b])*(Vol[b]-Vol[a])/PRESFAC;
                if (bPrint)
                {
                    fprintf(fplog,"  dpV = %10.3e  d = %10.3e",dpV,delta + dpV);
                }
                delta += dpV;
            }
            if (bPrint)
            {
                fprintf(fplog,"\n");
            }
            if (delta <= 0)
            {
                prob[i] = 1;
                bEx[i] = TRUE;
            }
            else
            {
                if (delta > 100)
                {
                    prob[i] = 0;
                }
                else
                {
                    prob[i] = exp(-delta);
                }
                bEx[i] = (rando(&(re->seed)) < prob[i]);
            }
            re->prob_sum[i] += prob[i];    
            if (bEx[i])
            {
                if (a == re->repl)
                {
                    exchange = b;
                }
                else if (b == re->repl)
                {
                    exchange = a;
                }
                re->nexchange[i]++;
            }
        }
        else
        {
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
    sfree(dvdl);
  
    re->nattempt[m]++;

    return exchange;
}

static void write_debug_x(t_state *state)
{
    int i;

    if (debug)
    {
        for(i=0; i<state->natoms; i+=10)
        {
            fprintf(debug,"dx %5d %10.5f %10.5f %10.5f\n",i,state->x[i][XX],state->x[i][YY],state->x[i][ZZ]);
        }
    }
}

gmx_bool replica_exchange(FILE *fplog,const t_commrec *cr,struct gmx_repl_ex *re,
                          t_state *state,real *ener,
                          t_state *state_local,
                          gmx_large_int_t step,real time)
{
    gmx_multisim_t *ms;
    int  exchange=-1,shift;
    gmx_bool bExchanged=FALSE;
    
    ms = cr->ms;
  
    if (MASTER(cr))
    {
        exchange = get_replica_exchange(fplog,ms,re,ener,det(state->box),
                                        step,time);
        bExchanged = (exchange >= 0);
    }
    
    if (PAR(cr))
    {
#ifdef GMX_MPI
        MPI_Bcast(&bExchanged,sizeof(gmx_bool),MPI_BYTE,MASTERRANK(cr),
                  cr->mpi_comm_mygroup);
#endif
    }
    
    if (bExchanged)
    {
        /* Exchange the states */

        if (PAR(cr))
        {
            /* Collect the global state on the master node */
            if (DOMAINDECOMP(cr))
            {
                dd_collect_state(cr->dd,state_local,state);
            }
            else
            {
                pd_collect_state(cr,state);
            }
        }
        
        if (MASTER(cr))
        {
            /* Exchange the global states between the master nodes */
            if (debug)
            {
                fprintf(debug,"Exchanging %d with %d\n",ms->sim,exchange);
            }
            exchange_state(ms,exchange,state);
            
            if (re->type == ereTEMP)
            {
                scale_velocities(state,sqrt(re->q[ms->sim]/re->q[exchange]));
            }
        }

        /* With domain decomposition the global state is distributed later */
        if (!DOMAINDECOMP(cr))
        {
            /* Copy the global state to the local state data structure */
            copy_state_nonatomdata(state,state_local);
            
            if (PAR(cr))
            {
                bcast_state(cr,state,FALSE);
            }
        }
    }
        
    return bExchanged;
}

void print_replica_exchange_statistics(FILE *fplog,struct gmx_repl_ex *re)
{
    real *prob;
    int  i;
  
    fprintf(fplog,"\nReplica exchange statistics\n");
    fprintf(fplog,"Repl  %d attempts, %d odd, %d even\n",
            re->nattempt[0]+re->nattempt[1],re->nattempt[1],re->nattempt[0]);

    snew(prob,re->nrepl);
  
    fprintf(fplog,"Repl  average probabilities:\n");
    for(i=1; i<re->nrepl; i++)
    {
        if (re->nattempt[i%2] == 0)
        {
            prob[i] = 0;
        }
        else
        {
            prob[i] =  re->prob_sum[i]/re->nattempt[i%2];
        }
    }
    print_ind(fplog,"",re->nrepl,re->ind,NULL);
    print_prob(fplog,"",re->nrepl,prob);

    fprintf(fplog,"Repl  number of exchanges:\n");
    print_ind(fplog,"",re->nrepl,re->ind,NULL);
    print_count(fplog,"",re->nrepl,re->nexchange);
  
    fprintf(fplog,"Repl  average number of exchanges:\n");
    for(i=1; i<re->nrepl; i++)
    {
        if (re->nattempt[i%2] == 0)
        {
            prob[i] = 0;
        }
        else
        {
            prob[i] =  ((real)re->nexchange[i])/re->nattempt[i%2];
        }
    }
    print_ind(fplog,"",re->nrepl,re->ind,NULL);
    print_prob(fplog,"",re->nrepl,prob);

    sfree(prob);
  
    fprintf(fplog,"\n");
}
