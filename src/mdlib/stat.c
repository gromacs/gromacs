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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdio.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "gmx_fatal.h"
#include "network.h"
#include "txtdump.h"
#include "names.h"
#include "physics.h"
#include "vec.h"
#include "maths.h"
#include "mvdata.h"
#include "main.h"
#include "force.h"
#include "vcm.h"
#include "smalloc.h"
#include "futil.h"
#include "network.h"
#include "rbin.h"
#include "tgroup.h"
#include "xtcio.h"
#include "gmxfio.h"
#include "trnio.h"
#include "statutil.h"
#include "domdec.h"
#include "partdec.h"
#include "constr.h"
#include "checkpoint.h"
#include "mdrun.h"

typedef struct gmx_global_stat
{
    t_bin *rb;
    int   *itc0;
    int   *itc1;
} t_gmx_global_stat;

gmx_global_stat_t global_stat_init(t_inputrec *ir)
{
    gmx_global_stat_t gs;

    snew(gs,1);
    
    gs->rb = mk_bin();
    snew(gs->itc0,ir->opts.ngtc);
    snew(gs->itc1,ir->opts.ngtc);

    return gs;
}

void global_stat_destroy(gmx_global_stat_t gs)
{
    destroy_bin(gs->rb);
    sfree(gs->itc0);
    sfree(gs->itc1);
    sfree(gs);
}

void global_stat(FILE *fplog,gmx_global_stat_t gs,
		 t_commrec *cr,gmx_enerdata_t *enerd,
		 tensor fvir,tensor svir,rvec mu_tot,
		 t_inputrec *inputrec,
		 gmx_ekindata_t *ekind,bool bSumEkinhOld,
		 gmx_constr_t constr,
		 t_vcm *vcm,int *nabnsb,
		 real *chkpt,real *terminate,
		 gmx_mtop_t *top_global, t_state *state_local)
{
  t_bin  *rb;
  int    *itc0,*itc1;
  int    ie,ifv,isv,irmsd=0,imu=0;
  int    idedl=0,idvdll=0,idvdlnl=0,iepl=0,icm=0,imass=0,ica=0,inb=0;
  int    ibnsb=-1,ichkpt=-1,iterminate;
  int    icj=-1,ici=-1,icx=-1;
  int    inn[egNR];
  int    j;
  real   *rmsd_data,rbnsb;
  double nb;
  
  rb   = gs->rb;
  itc0 = gs->itc0;
  itc1 = gs->itc1;

  reset_bin(rb);

  /* This routine copies all the data to be summed to one big buffer
   * using the t_bin struct. 
   */
  where();
  ie  = add_binr(rb,F_NRE,enerd->term);
  where();
  ifv = add_binr(rb,DIM*DIM,fvir[0]);
  where();
  isv = add_binr(rb,DIM*DIM,svir[0]);
  where();
  if (constr) {
    rmsd_data = constr_rmsd_data(constr);
    if (rmsd_data)
      irmsd = add_binr(rb,inputrec->eI==eiSD2 ? 3 : 2,rmsd_data);
  } else {
    rmsd_data = NULL;
  }
  if (!NEED_MUTOT(*inputrec)) {
    imu = add_binr(rb,DIM,mu_tot);
    where();
  }
  if (ekind) {
    for(j=0; (j<inputrec->opts.ngtc); j++) {
      if (bSumEkinhOld) {
	itc0[j]=add_binr(rb,DIM*DIM,ekind->tcstat[j].ekinh_old[0]);
      }
      itc1[j]=add_binr(rb,DIM*DIM,ekind->tcstat[j].ekinh[0]);
    }
    where();
    idedl = add_binr(rb,1,&(ekind->dekindl));
    where();
    ica   = add_binr(rb,1,&(ekind->cosacc.mvcos));
    where();
  }
  for(j=0; (j<egNR); j++)
    inn[j]=add_binr(rb,enerd->grpp.nener,enerd->grpp.ener[j]);
  where();
  if (inputrec->efep != efepNO) {
    idvdll  = add_bind(rb,1,&enerd->dvdl_lin);
    idvdlnl = add_bind(rb,1,&enerd->dvdl_nonlin);
    if (enerd->n_lambda > 0) {
      iepl = add_bind(rb,enerd->n_lambda,enerd->enerpart_lambda);
    }
  }
  if (vcm) {
    icm   = add_binr(rb,DIM*vcm->nr,vcm->group_p[0]);
    where();
    imass = add_binr(rb,vcm->nr,vcm->group_mass);
    where();
    if (vcm->mode == ecmANGULAR) {
      icj   = add_binr(rb,DIM*vcm->nr,vcm->group_j[0]);
      where();
      icx   = add_binr(rb,DIM*vcm->nr,vcm->group_x[0]);
      where();
      ici   = add_binr(rb,DIM*DIM*vcm->nr,vcm->group_i[0][0]);
      where();
    }
  }
  if (DOMAINDECOMP(cr)) {
    nb = cr->dd->nbonded_local;
    inb = add_bind(rb,1,&nb);
  }
  where();
  if (nabnsb) {
    rbnsb = *nabnsb;
    ibnsb = add_binr(rb,1,&rbnsb);
  }
  if (chkpt)
    ichkpt   = add_binr(rb,1,chkpt);
  iterminate = add_binr(rb,1,terminate);
  
  /* Global sum it all */
  if (debug)
    fprintf(debug,"Summing %d energies\n",rb->maxreal);
  sum_bin(rb,cr);
  where();
  
  /* Extract all the data locally */
  extract_binr(rb,ie  ,F_NRE,enerd->term);
  extract_binr(rb,ifv ,DIM*DIM,fvir[0]);
  extract_binr(rb,isv ,DIM*DIM,svir[0]);
  if (rmsd_data)
    extract_binr(rb,irmsd,inputrec->eI==eiSD2 ? 3 : 2,rmsd_data);
  if (!NEED_MUTOT(*inputrec))
    extract_binr(rb,imu,DIM,mu_tot);
  if (ekind) {
    for(j=0; (j<inputrec->opts.ngtc); j++) {
      if (bSumEkinhOld)
	extract_binr(rb,itc0[j],DIM*DIM,ekind->tcstat[j].ekinh_old[0]);
      extract_binr(rb,itc1[j],DIM*DIM,ekind->tcstat[j].ekinh[0]);
    }
    extract_binr(rb,idedl,1,&(ekind->dekindl));
    extract_binr(rb,ica,1,&(ekind->cosacc.mvcos));
    where();
  }
  for(j=0; (j<egNR); j++)
    extract_binr(rb,inn[j],enerd->grpp.nener,enerd->grpp.ener[j]);
  if (inputrec->efep != efepNO) {
    extract_bind(rb,idvdll ,1,&enerd->dvdl_lin);
    extract_bind(rb,idvdlnl,1,&enerd->dvdl_nonlin);
    if (enerd->n_lambda > 0) {
      extract_bind(rb,iepl,enerd->n_lambda,enerd->enerpart_lambda);
    }
  }
  if (vcm) {
    extract_binr(rb,icm,DIM*vcm->nr,vcm->group_p[0]);
    where();
    extract_binr(rb,imass,vcm->nr,vcm->group_mass);
    where();
    if (vcm->mode == ecmANGULAR) {
      extract_binr(rb,icj,DIM*vcm->nr,vcm->group_j[0]);
      where();
      extract_binr(rb,icx,DIM*vcm->nr,vcm->group_x[0]);
      where();
      extract_binr(rb,ici,DIM*DIM*vcm->nr,vcm->group_i[0][0]);
      where();
    }
  }
  if (DOMAINDECOMP(cr)) {
    extract_bind(rb,inb,1,&nb);
    if ((int)(nb + 0.5) != cr->dd->nbonded_global)
      dd_print_missing_interactions(fplog,cr,(int)(nb + 0.5),top_global,state_local);
  }
  where();
  if (nabnsb) {
    extract_binr(rb,ibnsb,1,&rbnsb);
    *nabnsb = (int)(rbnsb + 0.5);
  }
  where();
  if (chkpt)
    extract_binr(rb,ichkpt,1,chkpt);
  extract_binr(rb,iterminate,1,terminate);
  where();

  /* Small hack for temp only */
  enerd->term[F_TEMP] /= (cr->nnodes - cr->npmenodes);
}

int do_per_step(gmx_step_t step,gmx_step_t nstep)
{
  if (nstep != 0) 
    return ((step % nstep)==0); 
  else 
    return 0;
}

static void moveit(t_commrec *cr,
		   int left,int right,const char *s,rvec xx[])
{
  if (!xx) 
    return;

  move_rvecs(cr,FALSE,FALSE,left,right,
	     xx,NULL,(cr->nnodes-cr->npmenodes)-1,NULL);
}

void write_traj(FILE *fplog,t_commrec *cr,
                int fp_trn,bool bX,bool bV,bool bF,
                int fp_xtc,bool bXTC,int xtc_prec,
                char *fn_cpt,bool bCPT,
                gmx_mtop_t *top_global,
                int eIntegrator,int simulation_part,gmx_step_t step,double t,
                t_state *state_local,t_state *state_global,
                rvec *f_local,rvec *f_global,
                int *n_xtc,rvec **x_xtc)
{
    int     i,j;
    gmx_groups_t *groups;
    rvec    *xxtc;
    
#define MX(xvf) moveit(cr,GMX_LEFT,GMX_RIGHT,#xvf,xvf)
    
    if (DOMAINDECOMP(cr))
    {
        if (bCPT)
        {
            dd_collect_state(cr->dd,state_local,state_global);
        }
        else
        {
            if (bX || bXTC)
            {
                dd_collect_vec(cr->dd,state_local,state_local->x,
                               state_global->x);
            }
            if (bV)
            {
                dd_collect_vec(cr->dd,state_local,state_local->v,
                               state_global->v);
            }
        }
        if (bF)
        {
            dd_collect_vec(cr->dd,state_local,f_local,f_global);
        }
    }
    else
    {
        if (bCPT)
        {
            /* All pointers in state_local are equal to state_global,
             * but we need to copy the non-pointer entries.
             */
            state_global->lambda = state_local->lambda;
            copy_mat(state_local->box,state_global->box);
            copy_mat(state_local->boxv,state_global->boxv);
            copy_mat(state_local->pres_prev,state_global->pres_prev);
        }
        if (cr->nnodes > 1)
        {
            /* Particle decomposition, collect the data on the master node */
            if (bCPT)
            {
                if (state_local->flags & (1<<estX))   MX(state_global->x);
                if (state_local->flags & (1<<estV))   MX(state_global->v);
                if (state_local->flags & (1<<estSDX)) MX(state_global->sd_X);
                if (state_global->nrngi > 1) {
                    if (state_local->flags & (1<<estLD_RNG)) {
#ifdef GMX_MPI
                        MPI_Gather(state_local->ld_rng ,
                                   state_local->nrng*sizeof(state_local->ld_rng[0]),MPI_BYTE,
                                   state_global->ld_rng,
                                   state_local->nrng*sizeof(state_local->ld_rng[0]),MPI_BYTE,
                                   MASTERRANK(cr),cr->mpi_comm_mygroup);
#endif
                    }
                    if (state_local->flags & (1<<estLD_RNGI))
                    {
#ifdef GMX_MPI
                        MPI_Gather(state_local->ld_rngi,
                                   sizeof(state_local->ld_rngi[0]),MPI_BYTE,
                                   state_global->ld_rngi,
                                   sizeof(state_local->ld_rngi[0]),MPI_BYTE,
                                   MASTERRANK(cr),cr->mpi_comm_mygroup);
#endif
                    }
                }
            }
            else
            {
                if (bX || bXTC) MX(state_global->x);
                if (bV)         MX(state_global->v);
            }
            if (bF)         MX(f_global);
        }
    }
    
    if (MASTER(cr)) {
        if (bCPT) {
            write_checkpoint(fn_cpt,fplog,cr,eIntegrator,simulation_part,step,t,state_global);
        }
        
        if (bX || bV || bF) {
            fwrite_trn(fp_trn,step,t,state_local->lambda,
                       state_local->box,top_global->natoms,
                       bX ? state_global->x : NULL,
                       bV ? state_global->v : NULL,
                       bF ? f_global : NULL);
            if(gmx_fio_flush(fp_trn) != 0)
            {
                gmx_file("Cannot write trajectory; maybe you are out of quota?");
            }
        }      
        if (bXTC) {
            groups = &top_global->groups;
            if (*n_xtc == -1)
            {
                *n_xtc = 0;
                for(i=0; (i<top_global->natoms); i++)
                {
                    if (ggrpnr(groups,egcXTC,i) == 0)
                    {
                        (*n_xtc)++;
                    }
                }
                if (*n_xtc != top_global->natoms)
                {
                    snew(*x_xtc,*n_xtc);
                }
            }
            if (*n_xtc == top_global->natoms)
            {
                xxtc = state_global->x;
            }
            else
            {
                xxtc = *x_xtc;
                j = 0;
                for(i=0; (i<top_global->natoms); i++)
                {
                    if (ggrpnr(groups,egcXTC,i) == 0)
                    {
                        copy_rvec(state_global->x[i],xxtc[j++]);
                    }
                }
            }
            if (write_xtc(fp_xtc,*n_xtc,step,t,
                          state_local->box,xxtc,xtc_prec) == 0)
            {
                gmx_fatal(FARGS,"XTC error - maybe you are out of quota?");
            }
        }
    }
}

