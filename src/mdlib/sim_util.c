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

#ifdef GMX_CRAY_XT3
#include<catamount/dclock.h>
#endif


#include <stdio.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <math.h>
#include "typedefs.h"
#include "string2.h"
#include "gmxfio.h"
#include "smalloc.h"
#include "names.h"
#include "confio.h"
#include "mvdata.h"
#include "txtdump.h"
#include "pbc.h"
#include "chargegroup.h"
#include "vec.h"
#include "time.h"
#include "nrnb.h"
#include "mshift.h"
#include "mdrun.h"
#include "update.h"
#include "physics.h"
#include "main.h"
#include "mdatoms.h"
#include "force.h"
#include "bondf.h"
#include "pme.h"
#include "pppm.h"
#include "disre.h"
#include "orires.h"
#include "network.h"
#include "calcmu.h"
#include "constr.h"
#include "xvgr.h"
#include "trnio.h"
#include "xtcio.h"
#include "copyrite.h"
#include "gmx_random.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "gmx_wallcycle.h"
#include "genborn.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#include "qmmm.h"

#if 0
typedef struct gmx_timeprint {
    
} t_gmx_timeprint;
#endif

/* Portable version of ctime_r implemented in src/gmxlib/string2.c, but we do not want it declared in public installed headers */
char *
gmx_ctime_r(const time_t *clock,char *buf, int n);


double
gmx_gettime()
{
#ifdef HAVE_GETTIMEOFDAY
	struct timeval t;
	double seconds;
	
	gettimeofday(&t,NULL);
	
	seconds = (double) t.tv_sec + 1e-6*(double)t.tv_usec;
	
	return seconds;
#else
	double  seconds;
	
	seconds = time(NULL);
	
	return seconds;
#endif
}


#define difftime(end,start) ((double)(end)-(double)(start))

void print_time(FILE *out,gmx_runtime_t *runtime,gmx_large_int_t step,   
                t_inputrec *ir, t_commrec *cr)
{
    time_t finish;
    char   timebuf[STRLEN];
    double dt;
    char buf[48];
    
#ifndef GMX_THREADS
    if (!PAR(cr))
#endif
    {
        fprintf(out,"\r");
    }
    fprintf(out,"step %s",gmx_step_str(step,buf));
    if ((step >= ir->nstlist))
    {
        if ((ir->nstlist == 0) || ((step % ir->nstlist) == 0))
        {
            /* We have done a full cycle let's update time_per_step */
            runtime->last = gmx_gettime();
            dt = difftime(runtime->last,runtime->real);
            runtime->time_per_step = dt/(step - ir->init_step + 1);
        }
        dt = (ir->nsteps + ir->init_step - step)*runtime->time_per_step;
        
        if (ir->nsteps >= 0)
        {
            if (dt >= 300)
            {    
                finish = (time_t) (runtime->last + dt);
                gmx_ctime_r(&finish,timebuf,STRLEN);
                sprintf(buf,"%s",timebuf);
                buf[strlen(buf)-1]='\0';
                fprintf(out,", will finish %s",buf);
            }
            else
                fprintf(out,", remaining runtime: %5d s          ",(int)dt);
        }
        else
        {
            fprintf(out," performance: %.1f ns/day    ",
                    ir->delta_t/1000*24*60*60/runtime->time_per_step);
        }
    }
#ifndef GMX_THREADS
    if (PAR(cr))
    {
        fprintf(out,"\n");
    }
#endif

    fflush(out);
}

#ifdef NO_CLOCK 
#define clock() -1
#endif

static double set_proctime(gmx_runtime_t *runtime)
{
    double diff;
#ifdef GMX_CRAY_XT3
    double prev;

    prev = runtime->proc;
    runtime->proc = dclock();
    
    diff = runtime->proc - prev;
#else
    clock_t prev;

    prev = runtime->proc;
    runtime->proc = clock();

    diff = (double)(runtime->proc - prev)/(double)CLOCKS_PER_SEC;
#endif
    if (diff < 0)
    {
        /* The counter has probably looped, ignore this data */
        diff = 0;
    }

    return diff;
}

void runtime_start(gmx_runtime_t *runtime)
{
    runtime->real = gmx_gettime();
    runtime->proc          = 0;
    set_proctime(runtime);
    runtime->realtime      = 0;
    runtime->proctime      = 0;
    runtime->last          = 0;
    runtime->time_per_step = 0;
}

void runtime_end(gmx_runtime_t *runtime)
{
    double now;
    
    now = gmx_gettime();
    
    runtime->proctime += set_proctime(runtime);
    runtime->realtime  = now - runtime->real;
    runtime->real      = now;
}

void runtime_upd_proc(gmx_runtime_t *runtime)
{
    runtime->proctime += set_proctime(runtime);
}

void print_date_and_time(FILE *fplog,int nodeid,const char *title,
                         const gmx_runtime_t *runtime)
{
    int i;
    char timebuf[STRLEN];
    char time_string[STRLEN];
    time_t tmptime;

    if (fplog)
    {
        if (runtime != NULL)
        {
            tmptime = (time_t) runtime->real;
            gmx_ctime_r(&tmptime,timebuf,STRLEN);
        }
        else
        {
            tmptime = (time_t) gmx_gettime();
            gmx_ctime_r(&tmptime,timebuf,STRLEN);
        }
        for(i=0; timebuf[i]>=' '; i++)
        {
            time_string[i]=timebuf[i];
        }
        time_string[i]='\0';

        fprintf(fplog,"%s on node %d %s\n",title,nodeid,time_string);
    }
}

static void sum_forces(int start,int end,rvec f[],rvec flr[])
{
  int i;
  
  if (gmx_debug_at) {
    pr_rvecs(debug,0,"fsr",f+start,end-start);
    pr_rvecs(debug,0,"flr",flr+start,end-start);
  }
  for(i=start; (i<end); i++)
    rvec_inc(f[i],flr[i]);
}

/* 
 * calc_f_el calculates forces due to an electric field.
 *
 * force is kJ mol^-1 nm^-1 = e * kJ mol^-1 nm^-1 / e 
 *
 * Et[] contains the parameters for the time dependent 
 * part of the field (not yet used). 
 * Ex[] contains the parameters for
 * the spatial dependent part of the field. You can have cool periodic
 * fields in principle, but only a constant field is supported
 * now. 
 * The function should return the energy due to the electric field
 * (if any) but for now returns 0.
 *
 * WARNING:
 * There can be problems with the virial.
 * Since the field is not self-consistent this is unavoidable.
 * For neutral molecules the virial is correct within this approximation.
 * For neutral systems with many charged molecules the error is small.
 * But for systems with a net charge or a few charged molecules
 * the error can be significant when the field is high.
 * Solution: implement a self-consitent electric field into PME.
 */
static void calc_f_el(FILE *fp,int  start,int homenr,
                      real charge[],rvec x[],rvec f[],
                      t_cosines Ex[],t_cosines Et[],double t)
{
    rvec Ext;
    real t0;
    int  i,m;
    
    for(m=0; (m<DIM); m++)
    {
        if (Et[m].n > 0)
        {
            if (Et[m].n == 3)
            {
                t0 = Et[m].a[1];
                Ext[m] = cos(Et[m].a[0]*(t-t0))*exp(-sqr(t-t0)/(2.0*sqr(Et[m].a[2])));
            }
            else
            {
                Ext[m] = cos(Et[m].a[0]*t);
            }
        }
        else
        {
            Ext[m] = 1.0;
        }
        if (Ex[m].n > 0)
        {
            /* Convert the field strength from V/nm to MD-units */
            Ext[m] *= Ex[m].a[0]*FIELDFAC;
            for(i=start; (i<start+homenr); i++)
                f[i][m] += charge[i]*Ext[m];
        }
        else
        {
            Ext[m] = 0;
        }
    }
    if (fp != NULL)
    {
        fprintf(fp,"%10g  %10g  %10g  %10g #FIELD\n",t,
                Ext[XX]/FIELDFAC,Ext[YY]/FIELDFAC,Ext[ZZ]/FIELDFAC);
    }
}

static void calc_virial(FILE *fplog,int start,int homenr,rvec x[],rvec f[],
			tensor vir_part,t_graph *graph,matrix box,
			t_nrnb *nrnb,const t_forcerec *fr,int ePBC)
{
  int i,j;
  tensor virtest;

  /* The short-range virial from surrounding boxes */
  clear_mat(vir_part);
  calc_vir(fplog,SHIFTS,fr->shift_vec,fr->fshift,vir_part,ePBC==epbcSCREW,box);
  inc_nrnb(nrnb,eNR_VIRIAL,SHIFTS);
  
  /* Calculate partial virial, for local atoms only, based on short range. 
   * Total virial is computed in global_stat, called from do_md 
   */
  f_calc_vir(fplog,start,start+homenr,x,f,vir_part,graph,box);
  inc_nrnb(nrnb,eNR_VIRIAL,homenr);

  /* Add position restraint contribution */
  for(i=0; i<DIM; i++) {
    vir_part[i][i] += fr->vir_diag_posres[i];
  }

  /* Add wall contribution */
  for(i=0; i<DIM; i++) {
    vir_part[i][ZZ] += fr->vir_wall_z[i];
  }

  if (debug)
    pr_rvecs(debug,0,"vir_part",vir_part,DIM);
}

static void print_large_forces(FILE *fp,t_mdatoms *md,t_commrec *cr,
			       gmx_large_int_t step,real pforce,rvec *x,rvec *f)
{
  int  i;
  real pf2,fn2;
  char buf[STEPSTRSIZE];

  pf2 = sqr(pforce);
  for(i=md->start; i<md->start+md->homenr; i++) {
    fn2 = norm2(f[i]);
    /* We also catch NAN, if the compiler does not optimize this away. */
    if (fn2 >= pf2 || fn2 != fn2) {
      fprintf(fp,"step %s  atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
	      gmx_step_str(step,buf),
	      ddglatnr(cr->dd,i),x[i][XX],x[i][YY],x[i][ZZ],sqrt(fn2));
    }
  }
}

void do_force(FILE *fplog,t_commrec *cr,
              t_inputrec *inputrec,
              gmx_large_int_t step,t_nrnb *nrnb,gmx_wallcycle_t wcycle,
              gmx_localtop_t *top,
              gmx_mtop_t *mtop,
              gmx_groups_t *groups,
              matrix box,rvec x[],history_t *hist,
              rvec f[],
              tensor vir_force,
              t_mdatoms *mdatoms,
              gmx_enerdata_t *enerd,t_fcdata *fcd,
              real *lambda,t_graph *graph,
              t_forcerec *fr,gmx_vsite_t *vsite,rvec mu_tot,
              double t,FILE *field,gmx_edsam_t ed,
              gmx_bool bBornRadii,
              int flags)
{
    int    cg0,cg1,i,j;
    int    start,homenr;
    double mu[2*DIM]; 
    gmx_bool   bSepDVDL,bStateChanged,bNS,bFillGrid,bCalcCGCM,bBS;
    gmx_bool   bDoLongRange,bDoForces,bSepLRF;
    matrix boxs;
    real   e,v,dvdlambda[efptNR];
    real   dvdl_dum,lambda_dum;
    t_pbc  pbc;
    float  cycles_ppdpme,cycles_pme,cycles_seppme,cycles_force;
  
    start  = mdatoms->start;
    homenr = mdatoms->homenr;

    bSepDVDL = (fr->bSepDVDL && do_per_step(step,inputrec->nstlog));

    clear_mat(vir_force);

    if (PARTDECOMP(cr))
    {
        pd_cg_range(cr,&cg0,&cg1);
    }
    else
    {
        cg0 = 0;
        if (DOMAINDECOMP(cr))
        {
            cg1 = cr->dd->ncg_tot;
        }
        else
        {
            cg1 = top->cgs.nr;
        }
        if (fr->n_tpi > 0)
        {
            cg1--;
        }
    }

    bStateChanged = (flags & GMX_FORCE_STATECHANGED);
    bNS           = (flags & GMX_FORCE_NS) && (fr->bAllvsAll==FALSE); 
    bFillGrid     = (bNS && bStateChanged);
    bCalcCGCM     = (bFillGrid && !DOMAINDECOMP(cr));
    bDoLongRange  = (fr->bTwinRange && bNS && (flags & GMX_FORCE_DOLR));
    bDoForces     = (flags & GMX_FORCE_FORCES);
    bSepLRF       = (bDoLongRange && bDoForces && (flags & GMX_FORCE_SEPLRF));

    if (bStateChanged)
    {
        update_forcerec(fplog,fr,box);
        
        /* Calculate total (local) dipole moment in a temporary common array. 
         * This makes it possible to sum them over nodes faster.
         */
        calc_mu(start,homenr,
                x,mdatoms->chargeA,mdatoms->chargeB,mdatoms->nChargePerturbed,
                mu,mu+DIM);
    }
  
  if (fr->ePBC != epbcNONE) { 
    /* Compute shift vectors every step,
     * because of pressure coupling or box deformation!
     */
    if ((flags & GMX_FORCE_DYNAMICBOX) && bStateChanged)
      calc_shifts(box,fr->shift_vec);
    
    if (bCalcCGCM) { 
      put_charge_groups_in_box(fplog,cg0,cg1,fr->ePBC,box,
			       &(top->cgs),x,fr->cg_cm);
      inc_nrnb(nrnb,eNR_CGCM,homenr);
      inc_nrnb(nrnb,eNR_RESETX,cg1-cg0);
    } 
    else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph) {
      unshift_self(graph,box,x);
    }
  } 
  else if (bCalcCGCM) {
    calc_cgcm(fplog,cg0,cg1,&(top->cgs),x,fr->cg_cm);
    inc_nrnb(nrnb,eNR_CGCM,homenr);
  }
  
  if (bCalcCGCM) {
    if (PAR(cr)) {
      move_cgcm(fplog,cr,fr->cg_cm);
    }
    if (gmx_debug_at)
      pr_rvecs(debug,0,"cgcm",fr->cg_cm,top->cgs.nr);
  }

#ifdef GMX_MPI
  if (!(cr->duty & DUTY_PME)) {
    /* Send particle coordinates to the pme nodes.
     * Since this is only implemented for domain decomposition
     * and domain decomposition does not use the graph,
     * we do not need to worry about shifting.
     */    

    wallcycle_start(wcycle,ewcPP_PMESENDX);
    GMX_MPE_LOG(ev_send_coordinates_start);

    bBS = (inputrec->nwall == 2);
    if (bBS) {
      copy_mat(box,boxs);
      svmul(inputrec->wall_ewald_zfac,boxs[ZZ],boxs[ZZ]);
    }

    gmx_pme_send_x(cr,bBS ? boxs : box,x,
                   mdatoms->nChargePerturbed,lambda[efptCOUL],
                   ( flags & GMX_FORCE_VIRIAL),step);

    GMX_MPE_LOG(ev_send_coordinates_finish);
    wallcycle_stop(wcycle,ewcPP_PMESENDX);
  }
#endif /* GMX_MPI */

    /* Communicate coordinates and sum dipole if necessary */
    if (PAR(cr))
    {
        wallcycle_start(wcycle,ewcMOVEX);
        if (DOMAINDECOMP(cr))
        {
            dd_move_x(cr->dd,box,x);
        }
        else
        {
            move_x(fplog,cr,GMX_LEFT,GMX_RIGHT,x,nrnb);
        }
        /* When we don't need the total dipole we sum it in global_stat */
        if (bStateChanged && NEED_MUTOT(*inputrec))
        {
            gmx_sumd(2*DIM,mu,cr);
        }
        wallcycle_stop(wcycle,ewcMOVEX);
    }
    if (bStateChanged)
    {
        for(i=0; i<2; i++)
        {
            for(j=0;j<DIM;j++)
            {
                fr->mu_tot[i][j] = mu[i*DIM + j];
            }
        }
    }
    if (fr->efep == efepNO)
    {
        copy_rvec(fr->mu_tot[0],mu_tot);
    }
    else
    {
        for(j=0; j<DIM; j++)
        {
            mu_tot[j] =
                (1.0 - lambda[efptCOUL])*fr->mu_tot[0][j] + lambda[efptCOUL]*fr->mu_tot[1][j];
        }
    }

    /* Reset energies */
    reset_enerdata(&(inputrec->opts),fr,bNS,enerd,MASTER(cr));
    clear_rvecs(SHIFTS,fr->fshift);

    if (bNS)
    {
        wallcycle_start(wcycle,ewcNS);
        
        if (graph && bStateChanged)
        {
            /* Calculate intramolecular shift vectors to make molecules whole */
            mk_mshift(fplog,graph,fr->ePBC,box,x);
        }

        /* Reset long range forces if necessary */
        if (fr->bTwinRange)
        {
            /* Reset the (long-range) forces if necessary */
            clear_rvecs(fr->natoms_force_constr,bSepLRF ? fr->f_twin : f);
        }

        /* Do the actual neighbour searching and if twin range electrostatics
         * also do the calculation of long range forces and energies.
         */
        for (i=0;i<efptNR;i++) {dvdlambda[i] = 0;}
        ns(fplog,fr,x,box,
           groups,&(inputrec->opts),top,mdatoms,
           cr,nrnb,lambda,dvdlambda,&enerd->grpp,bFillGrid,
           bDoLongRange,bDoForces,bSepLRF ? fr->f_twin : f);
        if (bSepDVDL)
        {
            fprintf(fplog,sepdvdlformat,"LR non-bonded",0.0,dvdlambda);
        }
        enerd->dvdl_lin[efptVDW] += dvdlambda[efptVDW];
        enerd->dvdl_lin[efptCOUL] += dvdlambda[efptCOUL];
        
        wallcycle_stop(wcycle,ewcNS);
    }
	
    if (inputrec->implicit_solvent && bNS) 
    {
        make_gb_nblist(cr,inputrec->gb_algorithm,inputrec->rlist,
                       x,box,fr,&top->idef,graph,fr->born);
    }
	
    if (DOMAINDECOMP(cr))
    {
        if (!(cr->duty & DUTY_PME))
        {
            wallcycle_start(wcycle,ewcPPDURINGPME);
            dd_force_flop_start(cr->dd,nrnb);
        }
    }
	
    /* Start the force cycle counter.
     * This counter is stopped in do_forcelow_level.
     * No parallel communication should occur while this counter is running,
     * since that will interfere with the dynamic load balancing.
     */
    wallcycle_start(wcycle,ewcFORCE);
    
    if (bDoForces)
    {
        /* Reset forces for which the virial is calculated separately:
         * PME/Ewald forces if necessary */
        if (fr->bF_NoVirSum) 
        {
            if (flags & GMX_FORCE_VIRIAL)
            {
                fr->f_novirsum = fr->f_novirsum_alloc;
                GMX_BARRIER(cr->mpi_comm_mygroup);
                if (fr->bDomDec)
                {
                    clear_rvecs(fr->f_novirsum_n,fr->f_novirsum);
                }
                else
                {
                    clear_rvecs(homenr,fr->f_novirsum+start);
                }
                GMX_BARRIER(cr->mpi_comm_mygroup);
            }
            else
            {
                /* We are not calculating the pressure so we do not need
                 * a separate array for forces that do not contribute
                 * to the pressure.
                 */
                fr->f_novirsum = f;
            }
        }

        if (bSepLRF)
        {
            /* Add the long range forces to the short range forces */
            for(i=0; i<fr->natoms_force_constr; i++)
            {
                copy_rvec(fr->f_twin[i],f[i]);
            }
        }
        else if (!(fr->bTwinRange && bNS))
        {
            /* Clear the short-range forces */
            clear_rvecs(fr->natoms_force_constr,f);
        }

        clear_rvec(fr->vir_diag_posres);

        GMX_BARRIER(cr->mpi_comm_mygroup);
    }
    if (inputrec->ePull == epullCONSTRAINT)
    {
        clear_pull_forces(inputrec->pull);
    }

    /* update QMMMrec, if necessary */
    if(fr->bQMMM)
    {
        update_QMMMrec(cr,fr,x,mdatoms,box,top);
    }

    if ((flags & GMX_FORCE_BONDED) && top->idef.il[F_POSRES].nr > 0)
    {
        /* Position restraints always require full pbc */
        set_pbc(&pbc,inputrec->ePBC,box);
        v = posres(top->idef.il[F_POSRES].nr,top->idef.il[F_POSRES].iatoms,
                   top->idef.iparams_posres,
                   (const rvec*)x,fr->f_novirsum,fr->vir_diag_posres,
                   inputrec->ePBC==epbcNONE ? NULL : &pbc,lambda[efptRESTRAINT],&(dvdlambda[efptRESTRAINT]),
                   fr->rc_scaling,fr->ePBC,fr->posres_com,fr->posres_comB);
        if (bSepDVDL)
        {
            fprintf(fplog,sepdvdlformat,
                    interaction_function[F_POSRES].longname,v,dvdlambda);
        }
        enerd->term[F_POSRES] += v;
        /* This linear lambda dependence assumption is only correct
         * when only k depends on lambda,
         * not when the reference position depends on lambda.
         * grompp checks for this.  (verify this is still the case?)
         */
        enerd->dvdl_nonlin[efptRESTRAINT] += dvdlambda[efptRESTRAINT]; /* if just the force constant changes, this is linear, 
                                                                     but we can't be sure w/o additional checking that is
                                                                     hard to do at this level of code. Otherwise, 
                                                                     the dvdl is not differentiable */
         inc_nrnb(nrnb,eNR_POSRES,top->idef.il[F_POSRES].nr/2);
        if ((inputrec->fepvals->n_lambda > 0) && (flags & GMX_FORCE_DHDL))
        {
            for(i=0; i<enerd->n_lambda; i++)
            {
                lambda_dum = (i==0 ? lambda[efptRESTRAINT] : inputrec->fepvals->all_lambda[efptRESTRAINT][i-1]);
                v = posres(top->idef.il[F_POSRES].nr,top->idef.il[F_POSRES].iatoms,
                           top->idef.iparams_posres,
                           (const rvec*)x,NULL,NULL,
                           inputrec->ePBC==epbcNONE ? NULL : &pbc,lambda_dum,&dvdl_dum,
                           fr->rc_scaling,fr->ePBC,fr->posres_com,fr->posres_comB);
                enerd->enerpart_lambda[i] += v;
            }
        } 

   }

    /* Compute the bonded and non-bonded energies and optionally forces */    
    do_force_lowlevel(fplog,step,fr,inputrec,&(top->idef),
                      cr,nrnb,wcycle,mdatoms,&(inputrec->opts),
                      x,hist,f,enerd,fcd,mtop,top,fr->born,
                      &(top->atomtypes),bBornRadii,box,
                      inputrec->fepvals,lambda,graph,&(top->excls),fr->mu_tot,
                      flags,&cycles_pme);
    
    cycles_force = wallcycle_stop(wcycle,ewcFORCE);
    GMX_BARRIER(cr->mpi_comm_mygroup);
    
    if (ed)
    {
        do_flood(fplog,cr,x,f,ed,box,step);
    }
	
    if (DOMAINDECOMP(cr))
    {
        dd_force_flop_stop(cr->dd,nrnb);
        if (wcycle)
        {
            dd_cycles_add(cr->dd,cycles_force-cycles_pme,ddCyclF);
        }
    }
    
    if (bDoForces)
    {
        if (IR_ELEC_FIELD(*inputrec))
        {
            /* Compute forces due to electric field */
            calc_f_el(MASTER(cr) ? field : NULL,
                      start,homenr,mdatoms->chargeA,x,fr->f_novirsum,
                      inputrec->ex,inputrec->et,t);
        }
        
        /* Communicate the forces */
        if (PAR(cr))
        {
            wallcycle_start(wcycle,ewcMOVEF);
            if (DOMAINDECOMP(cr))
            {
                dd_move_f(cr->dd,f,fr->fshift);
                /* Do we need to communicate the separate force array
                 * for terms that do not contribute to the single sum virial?
                 * Position restraints and electric fields do not introduce
                 * inter-cg forces, only full electrostatics methods do.
                 * When we do not calculate the virial, fr->f_novirsum = f,
                 * so we have already communicated these forces.
                 */
                if (EEL_FULL(fr->eeltype) && cr->dd->n_intercg_excl &&
                    (flags & GMX_FORCE_VIRIAL))
                {
                    dd_move_f(cr->dd,fr->f_novirsum,NULL);
                }
                if (bSepLRF)
                {
                    /* We should not update the shift forces here,
                     * since f_twin is already included in f.
                     */
                    dd_move_f(cr->dd,fr->f_twin,NULL);
                }
            }
            else
            {
                pd_move_f(cr,f,nrnb);
                if (bSepLRF)
                {
                    pd_move_f(cr,fr->f_twin,nrnb);
                }
            }
            wallcycle_stop(wcycle,ewcMOVEF);
        }

        /* If we have NoVirSum forces, but we do not calculate the virial,
         * we sum fr->f_novirum=f later.
         */
        if (vsite && !(fr->bF_NoVirSum && !(flags & GMX_FORCE_VIRIAL)))
        {
            wallcycle_start(wcycle,ewcVSITESPREAD);
            spread_vsite_f(fplog,vsite,x,f,fr->fshift,nrnb,
                           &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
            wallcycle_stop(wcycle,ewcVSITESPREAD);

            if (bSepLRF)
            {
                wallcycle_start(wcycle,ewcVSITESPREAD);
                spread_vsite_f(fplog,vsite,x,fr->f_twin,NULL,
                               nrnb,
                               &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
                wallcycle_stop(wcycle,ewcVSITESPREAD);
            }
        }
        
        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Calculation of the virial must be done after vsites! */
            calc_virial(fplog,mdatoms->start,mdatoms->homenr,x,f,
                        vir_force,graph,box,nrnb,fr,inputrec->ePBC);
        }
    }

    if (inputrec->ePull == epullUMBRELLA || inputrec->ePull == epullCONST_F)
    {
        /* Calculate the center of mass forces, this requires communication,
         * which is why pull_potential is called close to other communication.
         * The virial contribution is calculated directly,
         * which is why we call pull_potential after calc_virial.
         */
        set_pbc(&pbc,inputrec->ePBC,box);
        dvdlambda[efptRESTRAINT] = 0; 
        enerd->term[F_COM_PULL] =
            pull_potential(inputrec->ePull,inputrec->pull,mdatoms,&pbc,
                           cr,t,lambda[efptRESTRAINT],x,f,vir_force,&(dvdlambda[efptRESTRAINT]));
        if (bSepDVDL)
        {
            fprintf(fplog,sepdvdlformat,"Com pull",enerd->term[F_COM_PULL],dvdlambda[efptRESTRAINT]);
        }
        enerd->dvdl_lin[efptRESTRAINT] += dvdlambda[efptRESTRAINT];
    }

    if (PAR(cr) && !(cr->duty & DUTY_PME))
    {
        cycles_ppdpme = wallcycle_stop(wcycle,ewcPPDURINGPME);
        dd_cycles_add(cr->dd,cycles_ppdpme,ddCyclPPduringPME);

        /* In case of node-splitting, the PP nodes receive the long-range 
         * forces, virial and energy from the PME nodes here.
         */    
        wallcycle_start(wcycle,ewcPP_PMEWAITRECVF);
        dvdlambda[efptCOUL] = 0;
        gmx_pme_receive_f(cr,fr->f_novirsum,fr->vir_el_recip,&e,&dvdlambda[efptCOUL],
                          &cycles_seppme);
        if (bSepDVDL)
        {
            fprintf(fplog,sepdvdlformat,"PME mesh",e,dvdlambda[efptCOUL]);
        }
        enerd->term[F_COUL_RECIP] += e;
        enerd->dvdl_lin[efptCOUL] += dvdlambda[efptCOUL];
        if (wcycle)
        {
            dd_cycles_add(cr->dd,cycles_seppme,ddCyclPME);
        }
        wallcycle_stop(wcycle,ewcPP_PMEWAITRECVF);
    }

    if (bDoForces && fr->bF_NoVirSum)
    {
        if (vsite)
        {
            /* Spread the mesh force on virtual sites to the other particles... 
             * This is parallellized. MPI communication is performed
             * if the constructing atoms aren't local.
             */
            wallcycle_start(wcycle,ewcVSITESPREAD);
            spread_vsite_f(fplog,vsite,x,fr->f_novirsum,NULL,nrnb,
                           &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
            wallcycle_stop(wcycle,ewcVSITESPREAD);
        }
        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Now add the forces, this is local */
            if (fr->bDomDec)
            {
                sum_forces(0,fr->f_novirsum_n,f,fr->f_novirsum);
            }
            else
            {
                sum_forces(start,start+homenr,f,fr->f_novirsum);
            }
            if (EEL_FULL(fr->eeltype))
            {
                /* Add the mesh contribution to the virial */
                m_add(vir_force,fr->vir_el_recip,vir_force);
            }
            if (debug)
            {
                pr_rvecs(debug,0,"vir_force",vir_force,DIM);
            }
        }
    }
    
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(inputrec->opts),enerd);
    
    if (fr->print_force >= 0 && bDoForces)
    {
        print_large_forces(stderr,mdatoms,cr,step,fr->print_force,x,f);
    }
}

void do_constrain_first(FILE *fplog,gmx_constr_t constr,
                        t_inputrec *ir,t_mdatoms *md,
                        t_state *state,rvec *f,
                        t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
                        t_forcerec *fr, gmx_localtop_t *top, tensor shake_vir)
{
    int    i,m,start,end;
    gmx_large_int_t step;
    double mass,tmass,vcm[4];
    real   dt=ir->delta_t;
    real   dvdl_dum;
    rvec   *savex;
    
    snew(savex,state->natoms);

    start = md->start;
    end   = md->homenr + start;
    
    if (debug)
        fprintf(debug,"vcm: start=%d, homenr=%d, end=%d\n",
                start,md->homenr,end);
    /* Do a first constrain to reset particles... */
    step = ir->init_step;
    if (fplog)
    {
        char buf[STEPSTRSIZE];
        fprintf(fplog,"\nConstraining the starting coordinates (step %s)\n",
                gmx_step_str(step,buf));
    }
    dvdl_dum = 0;
    
    /* constrain the current position */
    constrain(NULL,TRUE,FALSE,constr,&(top->idef),
              ir,NULL,cr,step,0,md,
              state->x,state->x,NULL,
              state->box,state->lambda[efptBONDED],&dvdl_dum,
              NULL,NULL,nrnb,econqCoord,ir->epc==epcMTTK,state->veta,state->veta);
    if (EI_VV(ir->eI)) 
    {
        /* constrain the inital velocity, and save it */
        /* also may be useful if we need the ekin from the halfstep for velocity verlet */
        /* might not yet treat veta correctly */
        constrain(NULL,TRUE,FALSE,constr,&(top->idef),
                  ir,NULL,cr,step,0,md,
                  state->x,state->v,state->v,
                  state->box,state->lambda[efptBONDED],&dvdl_dum,
                  NULL,NULL,nrnb,econqVeloc,ir->epc==epcMTTK,state->veta,state->veta);
    }
    /* constrain the inital velocities at t-dt/2 */
    if (EI_STATE_VELOCITY(ir->eI) && ir->eI!=eiVV)
    {
        for(i=start; (i<end); i++) 
        {
            for(m=0; (m<DIM); m++) 
            {
                /* Reverse the velocity */
                state->v[i][m] = -state->v[i][m];
                /* Store the position at t-dt in buf */
                savex[i][m] = state->x[i][m] + dt*state->v[i][m];
            }
        }
    /* Shake the positions at t=-dt with the positions at t=0                        
     * as reference coordinates.                                                     
         */
        if (fplog)
        {
            char buf[STEPSTRSIZE];
            fprintf(fplog,"\nConstraining the coordinates at t0-dt (step %s)\n",
                    gmx_step_str(step,buf));
        }
        dvdl_dum = 0;
        constrain(NULL,TRUE,FALSE,constr,&(top->idef),
                  ir,NULL,cr,step,-1,md,
                  state->x,savex,NULL,
                  state->box,state->lambda[efptBONDED],&dvdl_dum,
                  state->v,NULL,nrnb,econqCoord,ir->epc==epcMTTK,state->veta,state->veta);
        
        for(i=start; i<end; i++) {
            for(m=0; m<DIM; m++) {
                /* Re-reverse the velocities */
                state->v[i][m] = -state->v[i][m];
            }
        }
    }
    
    for(m=0; (m<4); m++)
        vcm[m] = 0;
    for(i=start; i<end; i++) {
        mass = md->massT[i];
        for(m=0; m<DIM; m++) {
            vcm[m] += state->v[i][m]*mass;
        }
        vcm[3] += mass;
    }
    
    if (ir->nstcomm != 0 || debug) {
        /* Compute the global sum of vcm */
        if (debug)
            fprintf(debug,"vcm: %8.3f  %8.3f  %8.3f,"
                    " total mass = %12.5e\n",vcm[XX],vcm[YY],vcm[ZZ],vcm[3]);
        if (PAR(cr))
            gmx_sumd(4,vcm,cr);
        tmass = vcm[3];
        for(m=0; (m<DIM); m++)
            vcm[m] /= tmass;
        if (debug) 
            fprintf(debug,"vcm: %8.3f  %8.3f  %8.3f,"
                    " total mass = %12.5e\n",vcm[XX],vcm[YY],vcm[ZZ],tmass);
        if (ir->nstcomm != 0) {
            /* Now we have the velocity of center of mass, let's remove it */
            for(i=start; (i<end); i++) {
                for(m=0; (m<DIM); m++)
                    state->v[i][m] -= vcm[m];
            }

        }
    }
    sfree(savex);
}

void calc_enervirdiff(FILE *fplog,int eDispCorr,t_forcerec *fr)
{
  double eners[2],virs[2],enersum,virsum,y0,f,g,h;
  double r0,r1,r,rc3,rc9,ea,eb,ec,pa,pb,pc,pd;
  double invscale,invscale2,invscale3;
  int    ri0,ri1,ri,i,offstart,offset;
  real   scale,*vdwtab; 

  fr->enershiftsix = 0;
  fr->enershifttwelve = 0;
  fr->enerdiffsix = 0;
  fr->enerdifftwelve = 0;
  fr->virdiffsix = 0;
  fr->virdifftwelve = 0;

  if (eDispCorr != edispcNO) {
    for(i=0; i<2; i++) {
      eners[i] = 0;
      virs[i]  = 0;
    }
    if ((fr->vdwtype == evdwSWITCH) || (fr->vdwtype == evdwSHIFT)) {
      if (fr->rvdw_switch == 0)
	gmx_fatal(FARGS,
		  "With dispersion correction rvdw-switch can not be zero "
		  "for vdw-type = %s",evdw_names[fr->vdwtype]);

      scale  = fr->nblists[0].tab.scale;
      vdwtab = fr->nblists[0].vdwtab;

      /* Round the cut-offs to exact table values for precision */
      ri0 = floor(fr->rvdw_switch*scale);
      ri1 = ceil(fr->rvdw*scale);
      r0  = ri0/scale;
      r1  = ri1/scale;
      rc3 = r0*r0*r0;
      rc9  = rc3*rc3*rc3;

      if (fr->vdwtype == evdwSHIFT) {
	/* Determine the constant energy shift below rvdw_switch */
	fr->enershiftsix    = (real)(-1.0/(rc3*rc3)) - vdwtab[8*ri0];
	fr->enershifttwelve = (real)( 1.0/(rc9*rc3)) - vdwtab[8*ri0 + 4];
      }
      /* Add the constant part from 0 to rvdw_switch.
       * This integration from 0 to rvdw_switch overcounts the number
       * of interactions by 1, as it also counts the self interaction.
       * We will correct for this later.
       */
      eners[0] += 4.0*M_PI*fr->enershiftsix*rc3/3.0;
      eners[1] += 4.0*M_PI*fr->enershifttwelve*rc3/3.0;
      
      invscale = 1.0/(scale);  
      invscale2 = invscale*invscale;
      invscale3 = invscale*invscale2;

      /* following summation derived from cubic spline definition,
	Numerical Recipies in C, second edition, p. 113-116.  Exact
	for the cubic spline.  We first calculate the negative of
	the energy from rvdw to rvdw_switch, assuming that g(r)=1,
	and then add the more standard, abrupt cutoff correction to
	that result, yielding the long-range correction for a
	switched function.  We perform both the pressure and energy
	loops at the same time for simplicity, as the computational
	cost is low. */
      
      for (i=0;i<2;i++) {
        enersum = 0.0; virsum = 0.0;
        if (i==0)
	  offstart = 0;
	else
	  offstart = 4;
	for (ri=ri0; ri<ri1; ri++) {
          r = ri*invscale;
          ea = invscale3;
          eb = 2.0*invscale2*r;
          ec = invscale*r*r;
          
          pa = invscale3;
          pb = 3.0*invscale2*r;
          pc = 3.0*invscale*r*r;
          pd = r*r*r;
          
          /* this "8" is from the packing in the vdwtab array - perhaps
	    should be #define'ed? */
          offset = 8*ri + offstart;
          y0 = vdwtab[offset];
          f = vdwtab[offset+1];
          g = vdwtab[offset+2];
          h = vdwtab[offset+3];
	  
          enersum += y0*(ea/3 + eb/2 + ec) + f*(ea/4 + eb/3 + ec/2)+
            g*(ea/5 + eb/4 + ec/3) + h*(ea/6 + eb/5 + ec/4);  
          virsum  +=  f*(pa/4 + pb/3 + pc/2 + pd) + 
            2*g*(pa/5 + pb/4 + pc/3 + pd/2) + 3*h*(pa/6 + pb/5 + pc/4 + pd/3);
	  
        }
        enersum *= 4.0*M_PI;
        virsum  *= 4.0*M_PI; 
        eners[i] -= enersum;
        virs[i]  -= virsum;
      }

      /* now add the correction for rvdw_switch to infinity */
      eners[0] += -4.0*M_PI/(3.0*rc3);
      eners[1] +=  4.0*M_PI/(9.0*rc9);
      virs[0]  +=  8.0*M_PI/rc3;
      virs[1]  += -16.0*M_PI/(3.0*rc9);
    } 
    else if ((fr->vdwtype == evdwCUT) || (fr->vdwtype == evdwUSER)) {
      if (fr->vdwtype == evdwUSER && fplog)
	fprintf(fplog,
		"WARNING: using dispersion correction with user tables\n");
      rc3  = fr->rvdw*fr->rvdw*fr->rvdw;
      rc9  = rc3*rc3*rc3;
      eners[0] += -4.0*M_PI/(3.0*rc3);
      eners[1] +=  4.0*M_PI/(9.0*rc9);
      virs[0]  +=  8.0*M_PI/rc3;
      virs[1]  += -16.0*M_PI/(3.0*rc9);
    } else {
      gmx_fatal(FARGS,
		"Dispersion correction is not implemented for vdw-type = %s",
		evdw_names[fr->vdwtype]);
    }
    fr->enerdiffsix    = eners[0];
    fr->enerdifftwelve = eners[1];
    /* The 0.5 is due to the Gromacs definition of the virial */
    fr->virdiffsix     = 0.5*virs[0];
    fr->virdifftwelve  = 0.5*virs[1];
  }
}

void calc_dispcorr(FILE *fplog,t_inputrec *ir,t_forcerec *fr,
                   gmx_large_int_t step,int natoms,
                   matrix box,real lambda,tensor pres,tensor virial,
                   real *prescorr, real *enercorr, real *dvdlcorr)
{
    gmx_bool bCorrAll,bCorrPres;
    real dvdlambda,invvol,dens,ninter,avcsix,avctwelve,enerdiff,svir=0,spres=0;
    int  m;
    
    *prescorr = 0;
    *enercorr = 0;
    *dvdlcorr = 0;
    
    clear_mat(virial);
    clear_mat(pres);
    
    if (ir->eDispCorr != edispcNO) {
        bCorrAll  = (ir->eDispCorr == edispcAllEner ||
                     ir->eDispCorr == edispcAllEnerPres);
        bCorrPres = (ir->eDispCorr == edispcEnerPres ||
                     ir->eDispCorr == edispcAllEnerPres);
        
        invvol = 1/det(box);
        if (fr->n_tpi) 
        {
            /* Only correct for the interactions with the inserted molecule */
            dens = (natoms - fr->n_tpi)*invvol;
            ninter = fr->n_tpi;
        } 
        else 
        {
            dens = natoms*invvol;
            ninter = 0.5*natoms;
        }
        
        if (ir->efep == efepNO) 
        {
            avcsix    = fr->avcsix[0];
            avctwelve = fr->avctwelve[0];
        } 
        else 
        {
            avcsix    = (1 - lambda)*fr->avcsix[0]    + lambda*fr->avcsix[1];
            avctwelve = (1 - lambda)*fr->avctwelve[0] + lambda*fr->avctwelve[1];
        }
        
        enerdiff = ninter*(dens*fr->enerdiffsix - fr->enershiftsix);
        *enercorr += avcsix*enerdiff;
        dvdlambda = 0.0;
        if (ir->efep != efepNO) 
        {
            dvdlambda += (fr->avcsix[1] - fr->avcsix[0])*enerdiff;
        }
        if (bCorrAll) 
        {
            enerdiff = ninter*(dens*fr->enerdifftwelve - fr->enershifttwelve);
            *enercorr += avctwelve*enerdiff;
            if (fr->efep != efepNO) 
            {
                dvdlambda += (fr->avctwelve[1] - fr->avctwelve[0])*enerdiff;
            }
        }
        
        if (bCorrPres) 
        {
            svir = ninter*dens*avcsix*fr->virdiffsix/3.0;
            if (ir->eDispCorr == edispcAllEnerPres)
            {
                svir += ninter*dens*avctwelve*fr->virdifftwelve/3.0;
            }
            /* The factor 2 is because of the Gromacs virial definition */
            spres = -2.0*invvol*svir*PRESFAC;
            
            for(m=0; m<DIM; m++) {
                virial[m][m] += svir;
                pres[m][m] += spres;
            }
            *prescorr += spres;
        }
        
        /* Can't currently control when it prints, for now, just print when degugging */
        if (debug)
        {
            if (bCorrAll) {
                fprintf(debug,"Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
                        avcsix,avctwelve);
            }
            if (bCorrPres) 
            {
                fprintf(debug,
                        "Long Range LJ corr.: Epot %10g, Pres: %10g, Vir: %10g\n",
                        *enercorr,spres,svir);
            }
            else
            {
                fprintf(debug,"Long Range LJ corr.: Epot %10g\n",*enercorr);
            }
        }
        
        if (fr->bSepDVDL && do_per_step(step,ir->nstlog))
        {
            fprintf(fplog,sepdvdlformat,"Dispersion correction",
                    *enercorr,dvdlambda);
        }
        if (fr->efep != efepNO) 
        {
            *dvdlcorr += dvdlambda;
        }
    }
}

void do_pbc_first(FILE *fplog,matrix box,t_forcerec *fr,
		  t_graph *graph,rvec x[])
{
  if (fplog)
    fprintf(fplog,"Removing pbc first time\n");
  calc_shifts(box,fr->shift_vec);
  if (graph) {
    mk_mshift(fplog,graph,fr->ePBC,box,x);
    if (gmx_debug_at)
      p_graph(debug,"do_pbc_first 1",graph);
    shift_self(graph,box,x);
    /* By doing an extra mk_mshift the molecules that are broken
     * because they were e.g. imported from another software
     * will be made whole again. Such are the healing powers
     * of GROMACS.
     */
    mk_mshift(fplog,graph,fr->ePBC,box,x);
    if (gmx_debug_at)
      p_graph(debug,"do_pbc_first 2",graph);
  }
  if (fplog)
    fprintf(fplog,"Done rmpbc\n");
}

static void low_do_pbc_mtop(FILE *fplog,int ePBC,matrix box,
			    gmx_mtop_t *mtop,rvec x[],
			    gmx_bool bFirst)
{
  t_graph *graph;
  int mb,as,mol;
  gmx_molblock_t *molb;

  if (bFirst && fplog)
    fprintf(fplog,"Removing pbc first time\n");

  snew(graph,1);
  as = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molb = &mtop->molblock[mb];
    if (molb->natoms_mol == 1 || 
	(!bFirst && mtop->moltype[molb->type].cgs.nr == 1)) {
      /* Just one atom or charge group in the molecule, no PBC required */
      as += molb->nmol*molb->natoms_mol;
    } else {
      /* Pass NULL iso fplog to avoid graph prints for each molecule type */
      mk_graph_ilist(NULL,mtop->moltype[molb->type].ilist,
		     0,molb->natoms_mol,FALSE,FALSE,graph);
      
      for(mol=0; mol<molb->nmol; mol++) {
	mk_mshift(fplog,graph,ePBC,box,x+as);
	
	shift_self(graph,box,x+as);
	/* The molecule is whole now.
	 * We don't need the second mk_mshift call as in do_pbc_first,
	 * since we no longer need this graph.
	 */
	
	as += molb->natoms_mol;
      }
      done_graph(graph);
    }
  }
  sfree(graph);
}

void do_pbc_first_mtop(FILE *fplog,int ePBC,matrix box,
		       gmx_mtop_t *mtop,rvec x[])
{
  low_do_pbc_mtop(fplog,ePBC,box,mtop,x,TRUE);
}

void do_pbc_mtop(FILE *fplog,int ePBC,matrix box,
		 gmx_mtop_t *mtop,rvec x[])
{
  low_do_pbc_mtop(fplog,ePBC,box,mtop,x,FALSE);
}

void finish_run(FILE *fplog,t_commrec *cr,const char *confout,
                t_inputrec *inputrec,
                t_nrnb nrnb[],gmx_wallcycle_t wcycle,
                gmx_runtime_t *runtime,
                gmx_bool bWriteStat)
{
  int    i,j;
  t_nrnb *nrnb_tot=NULL;
  real   delta_t;
  double nbfs,mflop;
  double cycles[ewcNR];

  wallcycle_sum(cr,wcycle,cycles);

  if (cr->nnodes > 1) {
    if (SIMMASTER(cr))
      snew(nrnb_tot,1);
#ifdef GMX_MPI
    MPI_Reduce(nrnb->n,nrnb_tot->n,eNRNB,MPI_DOUBLE,MPI_SUM,
               MASTERRANK(cr),cr->mpi_comm_mysim);
#endif  
  } else {
    nrnb_tot = nrnb;
  }
    
  if (SIMMASTER(cr)) {
    print_flop(fplog,nrnb_tot,&nbfs,&mflop);
    if (cr->nnodes > 1) {
      sfree(nrnb_tot);
    }
  }

  if ((cr->duty & DUTY_PP) && DOMAINDECOMP(cr)) {
    print_dd_statistics(cr,inputrec,fplog);
  }

#ifdef GMX_MPI
    if (PARTDECOMP(cr))
    {
        if (MASTER(cr))
        {
            t_nrnb     *nrnb_all;
            int        s;
            MPI_Status stat;

            snew(nrnb_all,cr->nnodes);
            nrnb_all[0] = *nrnb;
            for(s=1; s<cr->nnodes; s++)
            {
                MPI_Recv(nrnb_all[s].n,eNRNB,MPI_DOUBLE,s,0,
                         cr->mpi_comm_mysim,&stat);
            }
            pr_load(fplog,cr,nrnb_all);
            sfree(nrnb_all);
        }
        else
        {
            MPI_Send(nrnb->n,eNRNB,MPI_DOUBLE,MASTERRANK(cr),0,
                     cr->mpi_comm_mysim);
        }
    }
#endif  

  if (SIMMASTER(cr)) {
    wallcycle_print(fplog,cr->nnodes,cr->npmenodes,runtime->realtime,
                    wcycle,cycles);

    if (EI_DYNAMICS(inputrec->eI)) {
      delta_t = inputrec->delta_t;
    } else {
      delta_t = 0;
    }
    
    if (fplog) {
        print_perf(fplog,runtime->proctime,runtime->realtime,
                   cr->nnodes-cr->npmenodes,
                   runtime->nsteps_done,delta_t,nbfs,mflop);
    }
    if (bWriteStat) {
        print_perf(stderr,runtime->proctime,runtime->realtime,
                   cr->nnodes-cr->npmenodes,
                   runtime->nsteps_done,delta_t,nbfs,mflop);
    }

    /*
    runtime=inputrec->nsteps*inputrec->delta_t;
    if (bWriteStat) {
      if (cr->nnodes == 1)
	fprintf(stderr,"\n\n");
      print_perf(stderr,nodetime,realtime,runtime,&ntot,
		 cr->nnodes-cr->npmenodes,FALSE);
    }
    wallcycle_print(fplog,cr->nnodes,cr->npmenodes,realtime,wcycle,cycles);
    print_perf(fplog,nodetime,realtime,runtime,&ntot,cr->nnodes-cr->npmenodes,
	       TRUE);
    if (PARTDECOMP(cr))
      pr_load(fplog,cr,nrnb_all);
    if (cr->nnodes > 1)
      sfree(nrnb_all);
    */
  }
}

extern void initialize_lambdas(FILE *fplog,int efep,t_lambda *fep,int *fep_state,real *lambda,double *lam0)
{
    int i;
    
    if (efep==efepNO) {
        for (i=0;i<efptNR;i++)  {
            lambda[i] = 0.0;
            if (lam0) 
            {
                lam0[i] = 0.0;
            }
        }
        return;
    }
    
    for (i=0;i<efptNR;i++) 
    {
        *fep_state = fep->init_fep_state;
        /* overwrite lambda state with init_lambda for now for backwards compatibility */
        if (fep->init_lambda>=0) /* if it's -1, it was never initializd */
        {
            lambda[i] = fep->init_lambda;
            if (lam0) {
                lam0[i] = lambda[i];
            }
        }
        else 
        {
            lambda[i] = fep->all_lambda[i][*fep_state];
            if (lam0) {
                lam0[i] = lambda[i];
            }
        }
    } 
    
    /* Send to the log the information on the current lambdas */
    if (fplog != NULL) 
    {
        fprintf(fplog,"Initial lambda vector:[ ");
        for (i=0;i<efptNR;i++) 
        {
            fprintf(fplog,"%10.4f ",lambda[i]);
        }
        fprintf(fplog,"]\n");
    }
    return;
}

    
void init_md(FILE *fplog,
             t_commrec *cr,t_inputrec *ir,const output_env_t oenv,
             double *t,double *t0,
             real *lambda, int *fep_state, double *lam0,
             t_nrnb *nrnb,gmx_mtop_t *mtop,
             gmx_update_t *upd,
             int nfile,const t_filenm fnm[],
             gmx_mdoutf_t **outf,t_mdebin **mdebin,
             tensor force_vir,tensor shake_vir,rvec mu_tot,
             gmx_bool *bSimAnn,t_vcm **vcm, t_state *state, unsigned long Flags)
{
    int  i,j,n;
    real tmpt,mod;
	
    /* Initial values */
    *t = *t0       = ir->init_t;

    *bSimAnn=FALSE;
    for(i=0;i<ir->opts.ngtc;i++)
    {
        /* set bSimAnn if any group is being annealed */
        if(ir->opts.annealing[i]!=eannNO)
        {
            *bSimAnn = TRUE;
        }
    }
    if (*bSimAnn)
    {
        update_annealing_target_temp(&(ir->opts),ir->init_t);
    }

    /* Initialize lambda variables */
    initialize_lambdas(fplog,(ir->efep != efepNO),ir->fepvals,fep_state,lambda,lam0);
    
    if (upd)
    {
        *upd = init_update(fplog,ir);
    }

    
    if (vcm != NULL)
    {
        *vcm = init_vcm(fplog,&mtop->groups,ir);
    }
    
    if (EI_DYNAMICS(ir->eI) && !(Flags & MD_APPENDFILES))
    {
        if (ir->etc == etcBERENDSEN)
        {
            please_cite(fplog,"Berendsen84a");
        }
        if (ir->etc == etcVRESCALE)
        {
            please_cite(fplog,"Bussi2007a");
        }
    }
    
    init_nrnb(nrnb);
    
    if (nfile != -1)
    {
        *outf = init_mdoutf(nfile,fnm,Flags,cr,ir,oenv);

        *mdebin = init_mdebin((Flags & MD_APPENDFILES) ? NULL : (*outf)->fp_ene,
                              mtop,ir, (*outf)->fp_dhdl);
    }
    
    /* Initiate variables */  
    clear_mat(force_vir);
    clear_mat(shake_vir);
    clear_rvec(mu_tot);
    
    debug_gmx();
}


void GenerateGibbsProbabilities(real *scaled_lamee, real *p_k, real *pks, int minfep, int maxfep) {
  
    int ifep; 
    real denom, maxlamee;
    
    *pks = 0.0;
    maxlamee = scaled_lamee[minfep];
    /* find the maximum value */
    for (ifep=minfep;ifep<=maxfep;ifep++) 
    {
        if (scaled_lamee[ifep]>maxlamee) 
        {
            maxlamee = scaled_lamee[ifep];
        }
    }
    /* find the denominator */
    for (ifep=minfep;ifep<=maxfep;ifep++) 
    {
        *pks += exp(scaled_lamee[ifep]-maxlamee);
    }  
    /*numerators*/
    for (ifep=minfep;ifep<=maxfep;ifep++) 
    {
        p_k[ifep] = exp(scaled_lamee[ifep]-maxlamee) / *pks;
    }
}

real do_logsum(int N, real *a_n) {
    
    //     RETURN VALUE
    // log(\sum_{i=0}^(N-1) exp[a_n])
    real maxarg;
    real sum;
    int i;
    real logsum;
    //     compute maximum argument to exp(.)

    maxarg = a_n[0];
    for(i=1;i<N;i++) 
    {
        maxarg = max(maxarg,a_n[i]);
    }
  
    // compute sum of exp(a_n - maxarg)
    sum = 0.0;
    for (i=0;i<N;i++) 
    {
        sum = sum + exp(a_n[i] - maxarg);
    }
  
    //     compute log sum
    logsum = log(sum) + maxarg;
    return logsum;
} 

int FindMinimum(real *min_metric, int N) {
    
    real min_val;
    int min_nval,nval;
    
    min_nval = 0;
    min_val = min_metric[0];
    
    for (nval=0; nval<N; nval++) 
    {
        if (min_metric[nval] < min_val) 
        {
            min_val = min_metric[nval];
            min_nval = nval;
        }
    }
    return min_nval;
}

static gmx_bool CheckHistogramRatios(int nhisto, real *histo, real ratio, int flatcriteria) 
{
    
    int i;
    real nmean,maxval,minval;
    gmx_bool bIfFlat;
    
    nmean = 0;
    minval = histo[0];
    maxval = histo[0];
    for (i=0;i<nhisto;i++) 
    {
        if (histo[i] > maxval) 
        {
            maxval = histo[i];
        }
        
        if (histo[i] < minval) 
        {
            minval = histo[i];
        }
        nmean += histo[i];
    }

    if (nmean == 0) 
    {
        /* no samples! is bad!*/
        bIfFlat = FALSE;
        return bIfFlat;
    }
    nmean /= (real)nhisto;
    
    bIfFlat = FALSE;
    if (flatcriteria==flatbyEXTREMES) 
    {
        if (minval/maxval > ratio)
        {
            bIfFlat = TRUE;
        }
    }
    if (flatcriteria==flatbyAVERAGE) 
    {
        bIfFlat = TRUE;
        for (i=0;i<nhisto;i++) 
        {
            /* make sure that all points are in the ratio < x <  1/ratio range  */
            if (!((histo[i]/nmean < 1.0/ratio) && (histo[i]/nmean > ratio)))
            {
                bIfFlat = FALSE;
                break;
            }
        }
    }
    return bIfFlat;
}

static gmx_bool CheckIfDoneEquilibrating(t_lambda *fep, df_history_t *dfhist, gmx_large_int_t step) 
{

    int i,nlim,totalsamples;
    gmx_bool bDoneEquilibrating,bIfFlat;
    
    nlim = fep->n_lambda; /* for simplicity */

    /* assume we have equilibrated the weights, then check to see if any of the conditions are not met */

    /* calculate the total number of samples */
    switch (fep->elmceq) {
        
    case elmceqNO:
        /* We have not equilibrated, and won't, ever. */
        return FALSE;
    case elmceqYES:
        /* we have equilibrated -- we're done */
        return TRUE;
    case elmceqSTEPS:
        /* first, check if we are equilibrating by steps, if we're still under */
        if (step < fep->equil_steps)
        {
            bDoneEquilibrating = FALSE;
        }
        break;
    case elmceqSAMPLES:
        totalsamples = 0;
        for (i=0;i<nlim;i++) 
        {
            totalsamples += dfhist->n_at_lam[i];
        }
        if (totalsamples < fep->equil_samples) 
        {
            bDoneEquilibrating = FALSE;
        }
        break;
    case elmceqNUMATLAM:
        for (i=0;i<nlim;i++) 
        {
            if (dfhist->n_at_lam[i] < fep->equil_n_at_lam) /* we are still doing the initial sweep, so we're definitely not
                                                              done equilibrating*/
            {
                bDoneEquilibrating  = FALSE;
                break;
            }
        }
        break;
    case elmceqWLDELTA:
        if (EWL(fep->elamstats)) /* This check should be in readir as well, but 
                                    just to be sure */
        {
            if (dfhist->wl_delta > fep->equil_wl_delta) 
            {
                bDoneEquilibrating = FALSE;
            }
        }
        break;
    case elmceqRATIO:
        
        /* we can use the flatness as a judge of good weights, as long as
           we're not doing minvar, or Wang-Landau. 
           But turn off for now until we figure out exactly how we do this.
        */
        
        if (!(EWL(fep->elamstats) || fep->elamstats==elamstatsMINVAR)) 
        { 
            /* we want to use flatness -avoiding- the forced-through samples.  Plus, we need to convert to 
               floats for this histogram function. */

            real *modhisto;
            snew(modhisto,fep->n_lambda);
            for (i=0;i<fep->n_lambda;i++) {
                modhisto[i] = 1.0*(dfhist->n_at_lam[i]-fep->lmc_forced_nstart);
            }
            bIfFlat = CheckHistogramRatios(fep->n_lambda,modhisto,fep->equil_ratio,flatbyEXTREMES);
            sfree(modhisto);
            if (!bIfFlat) 
            {
                bDoneEquilibrating = FALSE;
            }
        }
    default:
        bDoneEquilibrating = TRUE;
    }
    /* one last case to go though, if we are doing slow growth to get initial values, we haven't finished equilibrating */
    
    if (fep->lmc_forced_nstart > 0) 
    {
        for (i=0;i<nlim;i++) 
        {
            if (dfhist->n_at_lam[i] < fep->lmc_forced_nstart) /* we are still doing the initial sweep, so we're definitely not
                                                                 done equilibrating*/
            { 
                bDoneEquilibrating = FALSE;
                break;
            }
        }
    }
    return bDoneEquilibrating;
}

static gmx_bool UpdateWeights(t_lambda *fep, df_history_t *dfhist, int fep_state, real *scaled_lamee, real *weighted_lamee, gmx_large_int_t step) 
{
    real maxdiff = 0.000000001;
    gmx_bool bSufficientSamples;
    int i, k, n, nz, indexi, indexk, min_n, max_n, nlam, nlim, totali;
    int n0,np1,nm1,nval,min_nvalm,min_nvalp,maxc;
    real chi_m1_0,chi_p1_0,chi_m2_0,chi_p2_0,chi_p1_m1,chi_p2_m1,chi_m1_p1,chi_m2_p1;
    real omega_m1_0,omega_p1_m1,omega_m1_p1,omega_p1_0,clam_osum;
    real de,de_function,dr,denom,maxdr,pks=0;
    real min_val,cnval,zero_sum_weights;
    real *omegam_array, *weightsm_array, *omegap_array, *weightsp_array, *varm_array, *varp_array, *dwp_array, *dwm_array;    
    real clam_varm, clam_varp, clam_weightsm, clam_weightsp, clam_minvar;
    real *lam_weights, *lam_minvar_corr, *lam_variance, *lam_dg, *p_k;
    int *nonzero;

    /* if we have equilibrated the weights, exit now */
    if (dfhist->bEquil) 
    {
        return FALSE;
    }
    
    /* for simplicity */
    nlim = fep->n_lambda;

    if (CheckIfDoneEquilibrating(fep,dfhist,step))
    {
        dfhist->bEquil = TRUE; 
        /* zero out the visited states so we know how many equilibrated states we have
           from here on out.*/
        for (i=0;i<nlim;i++)
        {
            dfhist->n_at_lam[i] = 0;
        }
        return TRUE;
    }
    
    /* If we reached this far, we have not equilibrated yet, keep on
       going resetting the weights */
    
    if (EWL(fep->elamstats))
    {
        if (fep->elamstats==elamstatsWL) 
        {
            dfhist->sum_weights[fep_state] -= dfhist->wl_delta; 
            dfhist->wl_histo[fep_state] += 1.0;
        } 
        else if (fep->elamstats==elamstatsGWL) 
        {
            snew(p_k,nlim);
            GenerateGibbsProbabilities(weighted_lamee,p_k,&pks,0,nlim-1);
            for (i=0;i<nlim;i++) 
            {
                dfhist->sum_weights[i] -= dfhist->wl_delta*p_k[i];
                dfhist->wl_histo[i] += p_k[i];
            }
            sfree(p_k);
        }
        
        zero_sum_weights =  dfhist->sum_weights[0];
        for (i=0;i<nlim;i++) 
        {
            dfhist->sum_weights[i] -= zero_sum_weights;
        }
    }
    
    if (fep->elamstats==elamstatsBARKER || fep->elamstats==elamstatsMETROPOLIS || fep->elamstats==elamstatsMINVAR) {
        
        maxc = 2*fep->c_range+1;
        
        snew(lam_dg,nlim);
        snew(lam_variance,nlim);
        
        snew(omegap_array,maxc);
        snew(weightsp_array,maxc);
        snew(varp_array,maxc);    
        snew(dwp_array,maxc);
        
        snew(omegam_array,maxc);
        snew(weightsm_array,maxc);
        snew(varm_array,maxc);
        snew(dwm_array,maxc);
        
        /* unpack the current lambdas -- we will only update 2 of these */
        
        for (i=0;i<nlim-1;i++) 
        { /* only through the second to last */
            lam_dg[i] = dfhist->sum_dg[i+1] - dfhist->sum_dg[i]; 
            lam_variance[i] = pow(dfhist->sum_variance[i+1],2) - pow(dfhist->sum_variance[i],2); 
        }
        
        /* accumulate running averages */
        for (nval = 0; nval<maxc; nval++) 
        {
            /* constants for later use */
            cnval = (real)(nval-fep->c_range); 
            /* actually, should be able to rewrite it w/o exponential, for better numerical stability */
            if (fep_state > 0) 
            {
                de = exp(cnval - (scaled_lamee[fep_state]-scaled_lamee[fep_state-1]));
                if (fep->elamstats==elamstatsBARKER || fep->elamstats==elamstatsMINVAR) 
                {
                    de_function = 1.0/(1.0+de);
                } 
                else if (fep->elamstats==elamstatsMETROPOLIS) 
                {
                    if (de < 1.0) 
                    {
                        de_function = 1.0;
                    } 
                    else 
                    {
                        de_function = 1.0/de;
                    }
                }
                dfhist->accum_m[fep_state][nval] += de_function;
                dfhist->accum_m2[fep_state][nval] += de_function*de_function;
            }
            
            if (fep_state < nlim-1) 
            {
                de = exp(-cnval + (scaled_lamee[fep_state+1]-scaled_lamee[fep_state]));
                if (fep->elamstats==elamstatsBARKER || fep->elamstats==elamstatsMINVAR) 
                {
                    de_function = 1.0/(1.0+de);
                } 
                else if (fep->elamstats==elamstatsMETROPOLIS) 
                {
                    if (de < 1.0) 
                    {
                        de_function = 1.0;
                    } 
                    else 
                    {
                        de_function = 1.0/de;
                    }
                }
                dfhist->accum_p[fep_state][nval] += de_function;
                dfhist->accum_p2[fep_state][nval] += de_function*de_function;
            }
            
            /* Metropolis transition and Barker transition (unoptimized Bennett) acceptance weight determination */
            
            n0  = dfhist->n_at_lam[fep_state];
            if (fep_state > 0) {nm1 = dfhist->n_at_lam[fep_state-1];} else {nm1 = 0;}     
            if (fep_state < nlim-1) {np1 = dfhist->n_at_lam[fep_state+1];} else {np1 = 0;}
            
            if (n0 > 0) 
            {
                chi_m1_0 = dfhist->accum_m[fep_state][nval]/n0; 
                chi_p1_0 = dfhist->accum_p[fep_state][nval]/n0;
                chi_m2_0 = dfhist->accum_m2[fep_state][nval]/n0; 
                chi_p2_0 = dfhist->accum_p2[fep_state][nval]/n0;
            }
            
            if ((fep_state > 0 ) && (nm1 > 0)) 
            {    
                chi_p1_m1 = dfhist->accum_p[fep_state-1][nval]/nm1;
                chi_p2_m1 = dfhist->accum_p2[fep_state-1][nval]/nm1;
            }
            
            if ((fep_state < nlim-1) && (np1 > 0)) 
            {
                chi_m1_p1 = dfhist->accum_m[fep_state+1][nval]/np1;	
                chi_m2_p1 = dfhist->accum_m2[fep_state+1][nval]/np1;	
            }
            
            omega_m1_0 = 0;
            omega_p1_0 = 0;
            clam_weightsm = 0;
            clam_weightsp = 0;
            clam_varm = 0;
            clam_varp = 0;
            
            if (fep_state > 0) 
            {
                omega_m1_0 = chi_m2_0/pow(chi_m1_0,2) - 1.0;
                if (nm1 > 0) 
                {
                    omega_p1_m1 = chi_p2_m1/pow(chi_p1_m1,2) - 1.0;
                    clam_weightsm = (log(chi_m1_0) - log(chi_p1_m1)) + cnval;
                    clam_varm = (1.0/n0)*(omega_m1_0) + (1.0/nm1)*(omega_p1_m1);
                } 
            }
            
            if (fep_state < nlim-1) 
            {  
                omega_p1_0 = chi_p2_0/pow(chi_p1_0,2) - 1.0;
                if (np1 > 0) 
                {
                    omega_m1_p1 = chi_m2_p1/pow(chi_m1_p1,2) - 1.0;
                    clam_weightsp = (log(chi_m1_p1) - log(chi_p1_0)) + cnval;
                    clam_varp = (1.0/np1)*(omega_m1_p1) + (1.0/n0)*(omega_p1_0);
                }
            }
            
            omegam_array[nval]             = omega_m1_0;
            weightsm_array[nval]           = clam_weightsm;
            varm_array[nval]               = clam_varm;
            if (nm1 > 0) 
            {
                dwm_array[nval]  = fabs( (cnval + log((1.0*n0)/nm1)) - lam_dg[fep_state-1] );
            } 
            else 
            {
                dwm_array[nval]  = fabs( cnval - lam_dg[fep_state-1] );
            }
            
            omegap_array[nval]             = omega_p1_0;
            weightsp_array[nval]           = clam_weightsp;
            varp_array[nval]               = clam_varp;
            if (np1 > 0) 
            {
                dwp_array[nval]  = fabs( (cnval + log((1.0*np1)/n0)) - lam_dg[fep_state] );
            } 
            else 
            {
                dwp_array[nval]  = fabs( cnval - lam_dg[fep_state] );
            }
            
        }
        
        /* find the C's closest to the old weights value */ 
        
        min_nvalm = FindMinimum(dwm_array,maxc);
        omega_m1_0    = omegam_array[min_nvalm];
        clam_weightsm = weightsm_array[min_nvalm];
        clam_varm     = varm_array[min_nvalm];
        
        min_nvalp = FindMinimum(dwp_array,maxc);
        omega_p1_0    = omegap_array[min_nvalp];
        clam_weightsp = weightsp_array[min_nvalp];
        clam_varp     = varp_array[min_nvalp];
        
        clam_osum = omega_m1_0 + omega_p1_0;
        clam_minvar = 0;
        if (clam_osum > 0) 
        {
            clam_minvar = 0.5*log(clam_osum);
        }
        
        if (fep_state > 0) 
        {
            lam_dg[fep_state-1] = clam_weightsm; 
            lam_variance[fep_state-1] = clam_varm;
        } 
        
        if (fep_state < nlim-1) 
        {
            lam_dg[fep_state] = clam_weightsp;
            lam_variance[fep_state] = clam_varp;
        }
        
        if (fep->elamstats==elamstatsMINVAR) 
        {
            bSufficientSamples = TRUE;
            /* make sure they are all past a threshold */
            for (i=0;i<nlim;i++) 
            {
                if (dfhist->n_at_lam[i] < fep->minvarmin) 
                {
                    bSufficientSamples = FALSE;
                } 
            }
            if (bSufficientSamples) 
            {
                /* MRS -- check this logic */
                dfhist->sum_minvar[fep_state] = clam_minvar; 
                if (fep_state==0) 
                {
                    for (i=0;i<nlim;i++) 
                    {
                        dfhist->sum_minvar[i]+=(fep->minvar_const-clam_minvar);
                    }
                    fep->minvar_const = clam_minvar;
                    dfhist->sum_minvar[fep_state] = 0.0;
                } 
                else 
                {
                    dfhist->sum_minvar[fep_state] -= fep->minvar_const;
                }
            }
        } 
        
        /* we need to rezero minvar now, since it could change at fep_state = 0 */
        dfhist->sum_dg[0] = 0.0;
        dfhist->sum_variance[0] = 0.0;
        dfhist->sum_weights[0] = dfhist->sum_dg[0] + dfhist->sum_minvar[0]; /* should be zero */
        
        for (i=1;i<nlim;i++) 
        {
            dfhist->sum_dg[i] = lam_dg[i-1] + dfhist->sum_dg[i-1];
            dfhist->sum_variance[i] = sqrt(lam_variance[i-1] + pow(dfhist->sum_variance[i-1],2)); 
            dfhist->sum_weights[i] = dfhist->sum_dg[i] + dfhist->sum_minvar[i]; 
        }
        
        sfree(lam_dg);
        sfree(lam_variance);
        
        sfree(omegam_array);
        sfree(weightsm_array);
        sfree(varm_array);
        sfree(dwm_array);
        
        sfree(omegap_array);    
        sfree(weightsp_array);
        sfree(varp_array);    
        sfree(dwp_array);
    }
    return FALSE;
}

static int ChooseNewLambda(FILE *log, t_inputrec *ir, df_history_t *dfhist, int fep_state, real *weighted_lamee, real *p_k, gmx_rng_t rng) 
{
    /* Choose new lambda value, and update transition matrix */
    
    int i,ifep,jfep,minfep,maxfep,nlim,lamnew,lamtrial,starting_fep_state;
    real r1,r2,pks,de_old,de_new,de,trialprob,tprob=0;
    real **Tij;
    real *propose,*accept,*remainder;
    real sum,pnorm;
    t_lambda *fep;

    fep = ir->fepvals;
    nlim = fep->n_lambda;
    starting_fep_state = fep_state;
    lamnew = fep_state; /* so that there is a default setting -- stays the same */ 
     
    if (!EWL(fep->elamstats))   /* ignore equilibrating the weights if using WL */
    {
        if ((fep->lmc_forced_nstart > 0) && (dfhist->n_at_lam[nlim-1] <= fep->lmc_forced_nstart))
        {
            /* Use a marching method to run through the lambdas and get preliminary free energy data, 
               before starting 'free' sampling.  We start free sampling when we have enough at each lambda */
            
            /* if we have enough at this lambda, move on to the next one */
            
            if (dfhist->n_at_lam[fep_state] == fep->lmc_forced_nstart)
            {
                lamnew = fep_state+1;
                if (lamnew == nlim)  /* whoops, stepped too far! */
                {
                    lamnew -= 1;
                }
            } 
            else
            {
                lamnew = fep_state;
            }
            return lamnew;
        }
    }

    snew(propose,fep->n_lambda);
    snew(accept,fep->n_lambda);
    snew(remainder,fep->n_lambda);
    
    for (i=0;i<fep->lmc_repeats;i++) 
    {
        
        for(ifep=0;ifep<nlim;ifep++) 
        {
            propose[ifep] = 0;
            accept[ifep] = 0;
        }
        
        if ((fep->elmcmove==elmcmoveGIBBS) || (fep->elmcmove==elmcmoveMETGIBBS)) 
        {
            /* use the Gibbs sampler, with variable locality */
            if (fep->gibbsdeltalam < 0) 
            {
                minfep = 0; 
                maxfep = nlim-1;
            } 
            else 
            {
                minfep = fep_state - fep->gibbsdeltalam;
                maxfep = fep_state + fep->gibbsdeltalam;
                if (minfep < 0) 
                { 
                    minfep = 0;
                }
                if (maxfep > nlim-1) 
                {
                    maxfep = nlim-1;
                } 
            }
            
            GenerateGibbsProbabilities(weighted_lamee,p_k,&pks,minfep,maxfep); 
            
            if (fep->elmcmove == elmcmoveGIBBS) 
            {
                for (ifep=minfep;ifep<=maxfep;ifep++) 
                {
                    propose[ifep] = p_k[ifep];
                    accept[ifep] = 1.0;
                }
                /* Gibbs sampling */
                r1 = gmx_rng_uniform_real(rng);
                for (lamnew=minfep;lamnew<=maxfep;lamnew++) 
                {
                    if (r1 <= p_k[lamnew]) 
                    {
                        break;
                    }
                    r1 -= p_k[lamnew];
                }
            } 
            else if (fep->elmcmove==elmcmoveMETGIBBS) 
            {
                
                /* Metropolized Gibbs sampling */
                for (ifep=minfep;ifep<=maxfep;ifep++) 
                {
                    remainder[ifep] = 1 - p_k[ifep];
                }
                
                /* find the proposal probabilities */

                if (remainder[fep_state] == 0) {
                    /* only the current state has any probability */
                    /* we have to stay at the current state */
                    lamnew=fep_state;
                } else {
                    for (ifep=minfep;ifep<=maxfep;ifep++) 
                    {
                        if (ifep != fep_state) 
                        {
                            propose[ifep] = p_k[ifep]/remainder[fep_state];
                        }
                        else 
                        {
                            propose[ifep] = 0;
                        }
                    }

                    r1 = gmx_rng_uniform_real(rng);
                    for (lamtrial=minfep;lamtrial<=maxfep;lamtrial++) 
                    {
                        pnorm = p_k[lamtrial]/remainder[fep_state];
                        if (lamtrial!=fep_state) 
                        {
                            if (r1 <= pnorm) 
                            {
                                break;
                            }
                            r1 -= pnorm;
                        }
                    }
                
                    /* we have now selected lamtrial according to p(lamtrial)/1-p(fep_state) */
                    tprob = 1.0;
                    /* trial probability is min{1,\frac{1 - p(old)}{1-p(new)} MRS 1/8/2008 */
                    trialprob = (remainder[fep_state])/(remainder[lamtrial]);
                    if (trialprob < tprob) 
                    {
                        tprob = trialprob;
                    }
                    r2 = gmx_rng_uniform_real(rng);
                    if (r2 < tprob) 
                    {
                        lamnew = lamtrial;
                    } 
                    else 
                    {
                        lamnew = fep_state;
                    }
                }

                /* now figure out the acceptance probability for each */
                for (ifep=minfep;ifep<=maxfep;ifep++) 
                {
                    tprob = 1.0;
                    if (remainder[ifep] != 0) {
                        trialprob = (remainder[fep_state])/(remainder[ifep]);
                    }
                    if (trialprob < tprob) 
                    {
                        tprob = trialprob;
                    }
                    /* probability for fep_state=0, but that's fine, it's never proposed! */
                    accept[ifep] = tprob;
                }
            }
            
            if (lamnew > maxfep) 
            {
                /* if its greater than maxfep, then something went wrong -- probably underflow in the calculation
                   of sum weights */
                fprintf(log,"Denominator %3d%17.10e\n",0,pks);	      	    
                fprintf(log,"  i               dE        numerator          weights\n");
                for (ifep=minfep;ifep<=maxfep;ifep++) 
                {
                    fprintf(log,"%3d,%17.10e%17.10e%17.10e\n",ifep,weighted_lamee[ifep],p_k[ifep],dfhist->sum_weights[ifep]);
                }
                gmx_fatal(FARGS,"Something wrong in choosing new lambda state with a Gibbs move -- probably underflow in weight determination");
            }
        } 
        else if ((fep->elmcmove==elmcmoveMETROPOLIS) || (fep->elmcmove==elmcmoveBARKER)) 
        {
            /* use the metropolis sampler with trial +/- 1 */
            r1 = gmx_rng_uniform_real(rng);
            if (r1 < 0.5) 
            {
                if (fep_state == 0) {
                    lamtrial = fep_state;
                }  
                else 
                {
                    lamtrial = fep_state-1;
                }
            } 
            else 
            {
                if (fep_state == nlim-1) {
                    lamtrial = fep_state;
                }  
                else 
                {
                    lamtrial = fep_state+1;
                }
            }
            
            de = weighted_lamee[lamtrial] - weighted_lamee[fep_state];
            if (fep->elmcmove==elmcmoveMETROPOLIS) 
            {
                tprob = 1.0;
                trialprob = exp(de);
                if (trialprob < tprob) 
                {
                    tprob = trialprob;
                }
                propose[fep_state] = 0; 
                propose[lamtrial] = 1.0; /* note that this overwrites the above line if fep_state = ntrial, which only occurs at the ends */
                accept[fep_state] = 1.0; /* doesn't actually matter, never proposed unless fep_state = ntrial, in which case it's 1.0 anyway */
                accept[lamtrial] = tprob;
                
            } 
            else if (fep->elmcmove==elmcmoveBARKER) 
            {
                tprob = 1.0/(1.0+exp(-de));
                
                propose[fep_state] = (1-tprob); 
                propose[lamtrial] += tprob; /* we add, to account for the fact that at the end, they might be the same point */ 
                accept[fep_state] = 1.0;
                accept[lamtrial] = 1.0;
            }
            
            r2 = gmx_rng_uniform_real(rng);
            if (r2 < tprob) {
                lamnew = lamtrial;
            } else {
                lamnew = fep_state;
            }
        }
        
        for (ifep=0;ifep<nlim;ifep++) 
        {
            dfhist->Tij[fep_state][ifep] += propose[ifep]*accept[ifep];
            dfhist->Tij[fep_state][fep_state] += propose[ifep]*(1.0-accept[ifep]);
        }
        fep_state = lamnew;
    }
    
    dfhist->Tij_empirical[starting_fep_state][lamnew] += 1.0;  
    
    sfree(propose);
    sfree(accept);
    sfree(remainder);
    
    return lamnew;
}

/* print out the weights to the log, along with current state */
extern void PrintFreeEnergyInfoToFile(FILE *outfile, t_lambda *fep, df_history_t *dfhist, 
                                      int nlam, int frequency, gmx_large_int_t step)
{
    int nlim,i,ifep,jfep;
    real dw,dg,dv,dm,Tprint;
    const char *print_names[efptNR] = {" FEPL","MassL","CoulL"," VdwL","BondL","RestT"};
    nlim = fep->n_lambda;

    if (mod(step,frequency)==0)
    {
        fprintf(outfile,"             MC-lambda information\n");
        fprintf(outfile,"  N");
        for (i=0;i<efptNR;i++) 
        {
            if (fep->separate_dvdl[i]) 
            {
                fprintf(outfile,"%7s",print_names[i]);
            }
        }
        fprintf(outfile,"    Count   ");
        if (fep->elamstats==elamstatsMINVAR) 
        {
            fprintf(outfile,"W(in kT)   G(in kT)  dG(in kT)  dV(in kT)\n");
        } 
        else 
        {
            fprintf(outfile,"G(in kT)  dG(in kT)\n");	    
        }
        for (ifep=0;ifep<nlim;ifep++) 
        {
            if (ifep==nlim-1) 
            {
                dw=0.0;
                dg=0.0;
                dv=0.0;
                dm=0.0;
            } 
            else 
            {
                dw = dfhist->sum_weights[ifep+1] - dfhist->sum_weights[ifep];
                dg = dfhist->sum_dg[ifep+1] - dfhist->sum_dg[ifep];
                dv = sqrt(pow(dfhist->sum_variance[ifep+1],2) - pow(dfhist->sum_variance[ifep],2));
                dm = dfhist->sum_minvar[ifep+1] - dfhist->sum_minvar[ifep];

            }
            fprintf(outfile,"%3d",(ifep+1));
            for (i=0;i<efptNR;i++) 
            {
                if (fep->separate_dvdl[i]) 
                {
                    fprintf(outfile,"%7.3f",fep->all_lambda[i][ifep]);
                }
            }
            if (EWL(fep->elamstats) && (!(dfhist->bEquil)))  /* if performing WL and still haven't equilibrated */
            {
                fprintf(outfile,"%9.3f",dfhist->wl_histo[ifep]);
            } 
            else   /* we have equilibrated weights */
            {
                fprintf(outfile,"%9d",dfhist->n_at_lam[ifep]);	      
            }
            if (fep->elamstats==elamstatsMINVAR) 
            {
                fprintf(outfile,"%11.5f%11.5f%11.5f%11.5f",dfhist->sum_weights[ifep],dfhist->sum_dg[ifep],dg,dv);
            } 
            else 
            {
                fprintf(outfile,"%11.5f%11.5f",dfhist->sum_weights[ifep],dw);
            }
            if (ifep == nlam) {
                fprintf(outfile," <<\n");
            } 
            else 
            {
                fprintf(outfile,"   \n");
            }
        }
        fprintf(outfile,"\n");
        
        if ((mod(step,fep->nstTij)==0) && (fep->nstTij > 0) && (step > 0))
        {            
            fprintf(outfile,"                     Transition Matrix\n");
            for (ifep=0;ifep<nlim;ifep++) 
            {
                fprintf(outfile,"%12d",(ifep+1));
            } 
            fprintf(outfile,"\n");
            for (ifep=0;ifep<nlim;ifep++) 
            {
                for (jfep=0;jfep<nlim;jfep++) 
                {
                    if (dfhist->n_at_lam[ifep] > 0) 
                    {
                        if (fep->bSymmetrizedTMatrix) 
                        {
                            Tprint = (dfhist->Tij[ifep][jfep]+dfhist->Tij[jfep][ifep])/(dfhist->n_at_lam[ifep]+dfhist->n_at_lam[jfep]);
                        } else {
                            Tprint = (dfhist->Tij[ifep][jfep])/(dfhist->n_at_lam[ifep]);
                        }
                    } 
                    else 
                    {
                        Tprint = 0.0;
                    }
                    fprintf(outfile,"%12.8f",Tprint);
                }
                fprintf(outfile,"%3d\n",(ifep+1));
             }
            
            fprintf(outfile,"                  Empirical Transition Matrix\n");
            for (ifep=0;ifep<nlim;ifep++) 
            {
                fprintf(outfile,"%12d",(ifep+1));
            }
            fprintf(outfile,"\n");	  
            for (ifep=0;ifep<nlim;ifep++) 
            {
                for (jfep=0;jfep<nlim;jfep++) 
                {
                    if (dfhist->n_at_lam[ifep] > 0) 
                    {
                        if (fep->bSymmetrizedTMatrix) 
                        {
                            Tprint = (dfhist->Tij_empirical[ifep][jfep]+dfhist->Tij_empirical[jfep][ifep])/(dfhist->n_at_lam[ifep]+dfhist->n_at_lam[jfep]);
                        } else {
                            Tprint = dfhist->Tij_empirical[ifep][jfep]/(dfhist->n_at_lam[ifep]);
                        }
                    } 
                    else 
                    {
                        Tprint = 0.0;
                    }
                    fprintf(outfile,"%12.8f",Tprint);
                }
                fprintf(outfile,"%3d\n",(ifep+1));
            }
        }
	}
}

extern void get_mc_state(gmx_rng_t rng,t_state *state)
{
    gmx_rng_get_state(rng,state->mc_rng,state->mc_rngi);
}

extern void set_mc_state(gmx_rng_t rng,t_state *state)
{
    gmx_rng_set_state(rng,state->mc_rng,state->mc_rngi[0]);
}

extern int ExpandedEnsembleDynamics(FILE *log,t_inputrec *ir, gmx_enerdata_t *enerd, 
                                     int nlam, df_history_t *dfhist, gmx_large_int_t step, gmx_rng_t mcrng)
{ 
    real *pfep_lamee,*p_k, *scaled_lamee, *weighted_lamee;
    int i,nlim,lamnew;
    real mckt,maxscaled=0,maxweighted=0;
    t_lambda *fep;
    gmx_bool bIfReset,bDoneEquilibrating=FALSE;

    fep = ir->fepvals;
	nlim = fep->n_lambda;

    snew(scaled_lamee,nlim);
    snew(weighted_lamee,nlim);
    snew(pfep_lamee,nlim);
    snew(p_k,nlim);

    if ((int)step == fep->nstfep)  /* this is the first step that this routine is visited.  Should make this more robust */
    {
        dfhist->wl_delta = fep->init_wl_delta;  /* MRS -- this would fit better somewhere else? */
        for (i=0;i<nlim;i++) {
            dfhist->sum_weights[i] = fep->init_lambda_weights[i];
            dfhist->sum_dg[i] = fep->init_lambda_weights[i];
        }
    } 
    
	/* update the count at the current lambda*/
	dfhist->n_at_lam[nlam]++;	
    
    /* need to calculate the PV term somewhere, but not needed here? Not until there's a lambda state that's 
       pressure controlled.*/
    /*
      pVTerm = 0;
      where does this PV term go?
      for (i=0;i<nlim;i++) 
      {	
      fep_lamee[i] += pVTerm; 
      }
    */
    
	/* set some constants */
	mckt = BOLTZ*fep->mc_temp; /* currently set to the system reft unless otherwise defined */

	/* determine the minimum value to avoid overflow.  Probably a better way to do this */
	/* we don't need to include the pressure term, since the volume is the same between the two.
	   is there some term we are neglecting, however? */
    
	for (i=0;i<nlim;i++) {
        scaled_lamee[i] = (enerd->enerpart_lambda[i+1]-enerd->enerpart_lambda[0])/mckt;
        
        /* save these energies for printing, so they don't get overwritten by the next step */
        /* they aren't overwritten in the non-free energy case, but we always print with these
           for simplicity */
        
        pfep_lamee[i] = scaled_lamee[i];
        
        weighted_lamee[i] = dfhist->sum_weights[i] - scaled_lamee[i]; 
        if (i==0) 
        {
            maxscaled = scaled_lamee[i];
            maxweighted = weighted_lamee[i]; 
        } 
        else 
        {
            if (scaled_lamee[i] > maxscaled) 
            {
                maxscaled = scaled_lamee[i];
            }
            if (weighted_lamee[i] > maxweighted) 
            {
                maxweighted = weighted_lamee[i];
            }
        }
	}
    
	for (i=0;i<nlim;i++) 
    {
        scaled_lamee[i] -= maxscaled;
        weighted_lamee[i] -= maxweighted;
	}
	
	/* update weights - we decide whether or not to actually this inside */

	bDoneEquilibrating = UpdateWeights(fep,dfhist,nlam,scaled_lamee,weighted_lamee,step);
    if (bDoneEquilibrating) 
    {
        if (log) {
            fprintf(log,"\nStep %d: Weights have equilibrated, using criteria: %s\n",(int)step,elmceq_names[fep->elmceq]);
        }
    }
    
    lamnew = ChooseNewLambda(log,ir,dfhist,nlam,weighted_lamee,p_k,mcrng);
    
	/* required for serial tempering? */
	/*fep->opts.ref_t[0]          = TemperatureBase*fep->temperature_lambdas[lamnew]; */

	/* now check on the Wang-Landau updating critera */
	
	if (EWL(fep->elamstats)) 
    {
        bIfReset = CheckHistogramRatios(fep->n_lambda,dfhist->wl_histo,fep->wl_ratio,flatbyAVERAGE);
        
        if (bIfReset) 
        {
            for (i=0;i<nlim;i++) 
            {
                dfhist->wl_histo[i] = 0;
            }
            dfhist->wl_delta *= fep->wl_scale;
            if (log) {
                fprintf(log,"\nStep %d: Wang-Landau weight is now %14.8f\n",(int)step,dfhist->wl_delta);
            }
        }
	}
    sfree(scaled_lamee);
    sfree(weighted_lamee);
    sfree(p_k);

    return lamnew;
}


