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
#include <math.h>
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
#include "xvgr.h"

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

static int filter_enerdterm(real *afrom, gmx_bool bToBuffer, real *ato,
                            gmx_bool bTemp, gmx_bool bPres, gmx_bool bEner) {
    int i,to,from;

    from = 0;
    to   = 0;
    for (i=0;i<F_NRE;i++)
    {
        if (bToBuffer)
        {
            from = i;
        }
        else
        {
            to = i;
        }
        switch (i) {
        case F_EKIN:
        case F_TEMP:
        case F_DKDL:
            if (bTemp)
            {
                ato[to++] = afrom[from++];
            }
            break;
        case F_PRES:    
        case F_PDISPCORR:
        case F_VTEMP:
            if (bPres)
            {
                ato[to++] = afrom[from++];
            }
            break;
        default:
            if (bEner)
            {
                ato[to++] = afrom[from++];
            }
            break;
        }
    }

    return to;
}

void global_stat(FILE *fplog,gmx_global_stat_t gs,
                 t_commrec *cr,gmx_enerdata_t *enerd,
                 tensor fvir,tensor svir,rvec mu_tot,
                 t_inputrec *inputrec,
                 gmx_ekindata_t *ekind,gmx_constr_t constr,
                 t_vcm *vcm,
                 int nsig,real *sig,
                 gmx_mtop_t *top_global, t_state *state_local, 
                 gmx_bool bSumEkinhOld, int flags)
/* instead of current system, gmx_booleans for summing virial, kinetic energy, and other terms */
{
  t_bin  *rb;
  int    *itc0,*itc1;
  int    ie=0,ifv=0,isv=0,irmsd=0,imu=0;
  int    idedl=0,idvdll=0,idvdlnl=0,iepl=0,icm=0,imass=0,ica=0,inb=0;
  int    isig=-1;
  int    icj=-1,ici=-1,icx=-1;
  int    inn[egNR];
  real   copyenerd[F_NRE];
  int    nener,j;
  real   *rmsd_data=NULL;
  double nb;
  gmx_bool   bVV,bTemp,bEner,bPres,bConstrVir,bEkinAveVel,bFirstIterate,bReadEkin;

  bVV           = EI_VV(inputrec->eI);
  bTemp         = flags & CGLO_TEMPERATURE;
  bEner         = flags & CGLO_ENERGY;
  bPres         = (flags & CGLO_PRESSURE); 
  bConstrVir    = (flags & CGLO_CONSTRAINT);
  bFirstIterate = (flags & CGLO_FIRSTITERATE);
  bEkinAveVel   = (inputrec->eI==eiVV || (inputrec->eI==eiVVAK && bPres));
  bReadEkin     = (flags & CGLO_READEKIN);

  rb   = gs->rb;
  itc0 = gs->itc0;
  itc1 = gs->itc1;
  

  reset_bin(rb);
  /* This routine copies all the data to be summed to one big buffer
   * using the t_bin struct. 
   */

  /* First, we neeed to identify which enerd->term should be
     communicated.  Temperature and pressure terms should only be
     communicated and summed when they need to be, to avoid repeating
     the sums and overcounting. */

  nener = filter_enerdterm(enerd->term,TRUE,copyenerd,bTemp,bPres,bEner);
  
  /* First, the data that needs to be communicated with velocity verlet every time
     This is just the constraint virial.*/
  if (bConstrVir) {
      isv = add_binr(rb,DIM*DIM,svir[0]);
      where();
  }
  
/* We need the force virial and the kinetic energy for the first time through with velocity verlet */
  if (bTemp || !bVV)
  {
      if (ekind) 
      {
          for(j=0; (j<inputrec->opts.ngtc); j++) 
          {
              if (bSumEkinhOld) 
              {
                  itc0[j]=add_binr(rb,DIM*DIM,ekind->tcstat[j].ekinh_old[0]);
              }
              if (bEkinAveVel && !bReadEkin) 
              {
                  itc1[j]=add_binr(rb,DIM*DIM,ekind->tcstat[j].ekinf[0]);
              } 
              else if (!bReadEkin)
              {
                  itc1[j]=add_binr(rb,DIM*DIM,ekind->tcstat[j].ekinh[0]);
              }
          }
          /* these probably need to be put into one of these categories */
          where();
          idedl = add_binr(rb,1,&(ekind->dekindl));
          where();
          ica   = add_binr(rb,1,&(ekind->cosacc.mvcos));
          where();
      }  
  }      
  where();
  
  if ((bPres || !bVV) && bFirstIterate)
  {
      ifv = add_binr(rb,DIM*DIM,fvir[0]);
  }


  if (bEner) 
  { 
      where();
      if (bFirstIterate) 
      {
          ie  = add_binr(rb,nener,copyenerd);
      }
      where();
      if (constr) 
      {
          rmsd_data = constr_rmsd_data(constr);
          if (rmsd_data) 
          {
              irmsd = add_binr(rb,inputrec->eI==eiSD2 ? 3 : 2,rmsd_data);
          }
      } 
      if (!NEED_MUTOT(*inputrec)) 
      {
          imu = add_binr(rb,DIM,mu_tot);
          where();
      }
      
      if (bFirstIterate) 
      {
          for(j=0; (j<egNR); j++)
          {
              inn[j]=add_binr(rb,enerd->grpp.nener,enerd->grpp.ener[j]);
          }
          where();
          if (inputrec->efep != efepNO) 
          {
              idvdll  = add_bind(rb,1,&enerd->dvdl_lin);
              idvdlnl = add_bind(rb,1,&enerd->dvdl_nonlin);
              if (enerd->n_lambda > 0) 
              {
                  iepl = add_bind(rb,enerd->n_lambda,enerd->enerpart_lambda);
              }
          }
      }

      if (vcm) 
      {
          icm   = add_binr(rb,DIM*vcm->nr,vcm->group_p[0]);
          where();
          imass = add_binr(rb,vcm->nr,vcm->group_mass);
          where();
          if (vcm->mode == ecmANGULAR) 
          {
              icj   = add_binr(rb,DIM*vcm->nr,vcm->group_j[0]);
              where();
              icx   = add_binr(rb,DIM*vcm->nr,vcm->group_x[0]);
              where();
              ici   = add_binr(rb,DIM*DIM*vcm->nr,vcm->group_i[0][0]);
              where();
          }
      }
  }
  if (DOMAINDECOMP(cr)) 
  {
      nb = cr->dd->nbonded_local;
      inb = add_bind(rb,1,&nb);
      }
  where();
  if (nsig > 0) 
  {
      isig = add_binr(rb,nsig,sig);
  }

  /* Global sum it all */
  if (debug)
  {
      fprintf(debug,"Summing %d energies\n",rb->maxreal);
  }
  sum_bin(rb,cr);
  where();

  /* Extract all the data locally */

  if (bConstrVir) 
  {
      extract_binr(rb,isv ,DIM*DIM,svir[0]);
  }

  /* We need the force virial and the kinetic energy for the first time through with velocity verlet */
  if (bTemp || !bVV)
  {
      if (ekind) 
      {
          for(j=0; (j<inputrec->opts.ngtc); j++) 
          {
              if (bSumEkinhOld)
              {
                  extract_binr(rb,itc0[j],DIM*DIM,ekind->tcstat[j].ekinh_old[0]);
              }
              if (bEkinAveVel && !bReadEkin) {
                  extract_binr(rb,itc1[j],DIM*DIM,ekind->tcstat[j].ekinf[0]);
              }
              else if (!bReadEkin)
              {
                  extract_binr(rb,itc1[j],DIM*DIM,ekind->tcstat[j].ekinh[0]);              
              }
          }
          extract_binr(rb,idedl,1,&(ekind->dekindl));
          extract_binr(rb,ica,1,&(ekind->cosacc.mvcos));
          where();
      }
  }
  if ((bPres || !bVV) && bFirstIterate)
  {
      extract_binr(rb,ifv ,DIM*DIM,fvir[0]);
  }

  if (bEner) 
  {
      if (bFirstIterate) 
      {
          extract_binr(rb,ie,nener,copyenerd);
          if (rmsd_data) 
          {
              extract_binr(rb,irmsd,inputrec->eI==eiSD2 ? 3 : 2,rmsd_data);
          }
          if (!NEED_MUTOT(*inputrec))
          {
              extract_binr(rb,imu,DIM,mu_tot);
          }

          for(j=0; (j<egNR); j++)
          {
              extract_binr(rb,inn[j],enerd->grpp.nener,enerd->grpp.ener[j]);
          }
          if (inputrec->efep != efepNO) 
          {
              extract_bind(rb,idvdll ,1,&enerd->dvdl_lin);
              extract_bind(rb,idvdlnl,1,&enerd->dvdl_nonlin);
              if (enerd->n_lambda > 0) 
              {
                  extract_bind(rb,iepl,enerd->n_lambda,enerd->enerpart_lambda);
              }
          }
          /* should this be here, or with ekin?*/
          if (vcm) 
          {
              extract_binr(rb,icm,DIM*vcm->nr,vcm->group_p[0]);
              where();
              extract_binr(rb,imass,vcm->nr,vcm->group_mass);
              where();
              if (vcm->mode == ecmANGULAR) 
              {
                  extract_binr(rb,icj,DIM*vcm->nr,vcm->group_j[0]);
                  where();
                  extract_binr(rb,icx,DIM*vcm->nr,vcm->group_x[0]);
                  where();
                  extract_binr(rb,ici,DIM*DIM*vcm->nr,vcm->group_i[0][0]);
                  where();
              }
          }
          if (DOMAINDECOMP(cr)) 
          {
              extract_bind(rb,inb,1,&nb);
              if ((int)(nb + 0.5) != cr->dd->nbonded_global) 
              {
                  dd_print_missing_interactions(fplog,cr,(int)(nb + 0.5),top_global,state_local);
              }
          }
          where();

          filter_enerdterm(copyenerd,FALSE,enerd->term,bTemp,bPres,bEner);    
/* Small hack for temp only - not entirely clear if still needed?*/
          /* enerd->term[F_TEMP] /= (cr->nnodes - cr->npmenodes); */
      }
  }

  if (nsig > 0) 
  {
      extract_binr(rb,isig,nsig,sig);
  }
  where();
}

int do_per_step(gmx_large_int_t step,gmx_large_int_t nstep)
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

gmx_mdoutf_t *init_mdoutf(int nfile,const t_filenm fnm[],int mdrun_flags,
                          const t_commrec *cr,const t_inputrec *ir,
                          const output_env_t oenv)
{
    gmx_mdoutf_t *of;
    char filemode[3];
    gmx_bool bAppendFiles;

    snew(of,1);

    of->fp_trn   = NULL;
    of->fp_ene   = NULL;
    of->fp_xtc   = NULL;
    of->fp_dhdl  = NULL;
    of->fp_field = NULL;
    
    of->eIntegrator     = ir->eI;
    of->simulation_part = ir->simulation_part;


    bAppendFiles = (mdrun_flags & MD_APPENDFILES);

    of->bKeepAndNumCPT = (mdrun_flags & MD_KEEPANDNUMCPT);

    sprintf(filemode, bAppendFiles ? "a+" : "w+");
        
    if (MASTER(cr))
    {
        if ((EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
#ifndef GMX_FAHCORE
            &&
            !(EI_DYNAMICS(ir->eI) &&
              ir->nstxout == 0 &&
              ir->nstvout == 0 &&
              ir->nstfout == 0)
#endif
	    )
        {
            of->fp_trn = open_trn(ftp2fn(efTRN,nfile,fnm), filemode);
        }

        if (EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
        {
            of->fp_ene = open_enx(ftp2fn(efEDR,nfile,fnm), filemode);
        }
        of->fn_cpt = opt2fn("-cpo",nfile,fnm);
        
        if (ir->efep != efepNO && ir->nstdhdl > 0 &&
            (ir->separate_dhdl_file == sepdhdlfileYES ) && 
            EI_DYNAMICS(ir->eI))
        {
            if (bAppendFiles)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-dhdl",nfile,fnm),filemode);
            }
            else
            {
                of->fp_dhdl = open_dhdl(opt2fn("-dhdl",nfile,fnm),ir,oenv);
            }
        }
        
        if (opt2bSet("-field",nfile,fnm) &&
            (ir->ex[XX].n || ir->ex[YY].n || ir->ex[ZZ].n))
        {
            if (bAppendFiles)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-field",nfile,fnm),
                                            filemode);
            }
            else
            {				  
                of->fp_field = xvgropen(opt2fn("-field",nfile,fnm),
                                        "Applied electric field","Time (ps)",
                                        "E (V/nm)",oenv);
            }
        }
    }
    if (EI_DYNAMICS(ir->eI) &&
        ir->nstxtcout > 0)
    {
        of->fp_xtc = open_xtc(ftp2fn(efXTC,nfile,fnm), filemode, cr);
        of->xtc_prec = ir->xtcprec;
    }
    return of;
}

void done_mdoutf(gmx_mdoutf_t *of)
{
    if (of->fp_ene != NULL)
    {
        close_enx(of->fp_ene);
    }
    if (of->fp_xtc)
    {
        close_xtc(of->fp_xtc);
    }
    if (of->fp_trn)
    {
        close_trn(of->fp_trn);
    }
    if (of->fp_dhdl != NULL)
    {
        gmx_fio_fclose(of->fp_dhdl);
    }
    if (of->fp_field != NULL)
    {
        gmx_fio_fclose(of->fp_field);
    }

    sfree(of);
}

// This is used for use with write_buf inside write_traj, it copies the 
// original state_local (old_sl) and copies it into new_sl
static int copy_state_local(t_state *new_sl,t_state *old_sl) 
{
	int *cg_gl_new = new_sl->cg_gl;
	rvec *x_new = new_sl->x;

	memcpy(new_sl,old_sl,sizeof(t_state));
	new_sl->cg_gl = cg_gl_new;
	new_sl->x = x_new;
        srenew(new_sl->cg_gl,old_sl->cg_gl_nalloc);
        srenew(new_sl->x,old_sl->nalloc);
	memcpy(new_sl->cg_gl, old_sl->cg_gl, sizeof(int) * old_sl->cg_gl_nalloc);
	memcpy(new_sl->x, old_sl->x, sizeof(rvec) * old_sl->natoms);
	return 0;
}



void write_traj(FILE *fplog,t_commrec *cr,
                gmx_mdoutf_t *of,
                int mdof_flags,
                gmx_mtop_t *top_global,
                gmx_large_int_t step,double t,
                t_state *state_local,t_state *state_global,
                rvec *f_local,rvec *f_global,
                int *n_xtc,rvec **x_xtc,
                t_inputrec *ir, gmx_bool bLastStep,
                t_write_buffer* write_buf,
                gmx_wallcycle_t wcycle)
{
    int     i,j;
    gmx_groups_t *groups;
    rvec    *xxtc = NULL;
    rvec *local_v;
    rvec *global_v;
    
    int bufferStep = 0;
    gmx_bool bBuffer = cr->nionodes > 1; // Used to determine if buffers will be used
    gmx_bool writeXTCNow = (mdof_flags & MDOF_XTC);

    if (bBuffer)// If buffering will be used
    {
        bufferStep = (step/ir->nstxtcout - (int)ceil((double)write_buf->step_after_checkpoint/ir->nstxtcout)) % cr->nionodes;// bufferStep = step/(how often to write) - (round up) step_at_checkpoint/(how often to write)  MOD (how often we actually do write)
        writeXTCNow = ((mdof_flags & MDOF_XTC) && bufferStep == cr->nionodes-1)   //write XTC in this step and buffer is full
                || (ir->nstxtcout>0 &&  bufferStep < cr->nionodes-1 && (bLastStep || (mdof_flags & MDOF_CPT) || (mdof_flags & MDOF_X)));// XTC is written AND we haven't just written because buffer was full AND its the last step OR its a checkpoint  OR write uncompressed X

        if ((mdof_flags & MDOF_CPT) || (mdof_flags & MDOF_X))
        {
            write_buf->step_after_checkpoint = step + 1;
        }
    }

#define MX(xvf) moveit(cr,GMX_LEFT,GMX_RIGHT,#xvf,xvf)


    /* MRS -- defining these variables is to manage the difference
     * between half step and full step velocities, but there must be a better way . . . */

    local_v  = state_local->v;
    global_v = state_global->v;

    if (DOMAINDECOMP(cr))
    {
        if (mdof_flags & MDOF_CPT)

        {
            dd_collect_state(cr->dd,state_local,state_global,wcycle);
        }
        else
        {
            //Collect X if writing X. Also Collect if writing XTC and not buffering
            if ((mdof_flags & MDOF_X) || ((mdof_flags & MDOF_XTC) && !bBuffer))
            {
                dd_collect_vec(cr->dd,state_local,state_local->x,
                               state_global->x,wcycle);
            }
            if (mdof_flags & MDOF_V)
            {
                dd_collect_vec(cr->dd,state_local,local_v,
                               global_v,wcycle);
            }
        }
        if (mdof_flags & MDOF_F)
        {
            dd_collect_vec(cr->dd,state_local,f_local,f_global,wcycle);
        }
        //Could be optimized by not collecting all coordinates but only those in the xtc selection.
        if (bBuffer) {
            if (mdof_flags & MDOF_XTC)
            {
                //This block of code copies the current dd and state_local to buffers to prepare for writing later.
                //The last frame being buffered (then writeXTCNow is TRUE) is always collected on the MASTER
                //We always remember the time on the master in case we write a checkpoint next before writing the next XTC frame
                if (MASTER(cr) || bufferStep == cr->dd->iorank)
                {
                    write_buf->step=step;
                    write_buf->t=t;
                }
                copy_dd(write_buf->dd[bufferStep],cr->dd);
                copy_state_local(write_buf->state_local[bufferStep],state_local);
            }

            if (writeXTCNow)
            {
                for (i = 0; i <= bufferStep; i++)//Collect each buffered frame to one of the IO nodes. The data is collected to the node with rank write_buf->dd[i]->masterrank.
                {
                    write_buf->dd[i]->masterrank = cr->dd->iorank2ddrank[i];
                    if (!(i==bufferStep && ((mdof_flags & MDOF_CPT) || (mdof_flags & MDOF_X))))
                    {
                        dd_collect_vec(write_buf->dd[i],write_buf->state_local[i],write_buf->state_local[i]->x,state_global->x,wcycle);
                    }
                }
            }
        }
    }
    else
    {
        if (mdof_flags & MDOF_CPT)
        {
            /* All pointers in state_local are equal to state_global,
             * but we need to copy the non-pointer entries.
             */
            state_global->lambda = state_local->lambda;
            state_global->veta = state_local->veta;
            state_global->vol0 = state_local->vol0;
            copy_mat(state_local->box,state_global->box);
            copy_mat(state_local->boxv,state_global->boxv);
            copy_mat(state_local->svir_prev,state_global->svir_prev);
            copy_mat(state_local->fvir_prev,state_global->fvir_prev);
            copy_mat(state_local->pres_prev,state_global->pres_prev);
        }
        if (cr->nnodes > 1)
        {
            /* Particle decomposition, collect the data on the master node */
            if (mdof_flags & MDOF_CPT)
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
                if (mdof_flags & (MDOF_X | MDOF_XTC)) MX(state_global->x);
                if (mdof_flags & MDOF_V)              MX(global_v);
            }
            if (mdof_flags & MDOF_F) MX(f_global);
         }
     }     

    /* The order of write_checkpoint and write_xtc/fwrite_trn is crucial, because the position of the trajectories is stored in the checkpoint.
     * The checkpoint is written before the current step because the current step is written to the trajectory when appending from the checkpoint
     * This is because write_traj is not called at the end of the step but in the middle of the step. It is called before the (coordinate) update and before the step number is increased.
     * Mainly the forces have been computed (to allow to write both coordinates and forces).
     *
     * If we buffer we need to call write_traj after write_xtc. Otherwise the xtc position stored in the checkpoint would point to the position before any of the buffered steps.
     * In that case the position is calculated to be correct (not including the current frame) even though we have already written all frames.
     */
     if (!bBuffer)
	 {
		if (mdof_flags & MDOF_CPT)
		{
		 write_checkpoint(of->fn_cpt,of->bKeepAndNumCPT,
						  fplog,cr,of->eIntegrator,
						  of->simulation_part,step,t,state_global,wcycle);
		}
	 }

     if (writeXTCNow && IONODE(cr)) {  //this is an IO node (we have to call write_traj on all IO nodes!)

		gmx_bool bWrite = MASTER(cr) || cr->dd->iorank<bufferStep;  //this node is actually writing
		int write_step;
		real write_t;

		if (bWrite) { // If this node is one of the writing nodes
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
						copy_rvec (state_global->x[i], xxtc[j++]);

					}
				}
			}
		}
		if (bBuffer)
		{
			write_step = write_buf->step;
			write_t = write_buf->t;
		}
		else
		{
			write_step = step;
			write_t = t;
		}
		if (write_xtc(of->fp_xtc,*n_xtc,write_step,write_t,
				  state_local->box,xxtc,of->xtc_prec,bWrite,wcycle) == 0)//If it is NOT ACTUALLY being written
		{
			gmx_fatal(FARGS,"XTC error - maybe you are out of quota?");
		}
		gmx_fio_check_file_position(of->fp_xtc);
     }
     if (bBuffer)
     {
        if (mdof_flags & MDOF_CPT)
		{
		 write_checkpoint(of->fn_cpt,of->bKeepAndNumCPT,
						  fplog,cr,of->eIntegrator,
						  of->simulation_part,step,t,state_global,wcycle);
		}
     }

     if (MASTER(cr))
     {

         if (mdof_flags & (MDOF_X | MDOF_V | MDOF_F))
         {
            fwrite_trn(of->fp_trn,step,t,state_local->lambda,
                       state_local->box,top_global->natoms,
                       (mdof_flags & MDOF_X) ? state_global->x : NULL,
                       (mdof_flags & MDOF_V) ? global_v : NULL,
                       (mdof_flags & MDOF_F) ? f_global : NULL);
            if (gmx_fio_flush(of->fp_trn) != 0)
            {
                gmx_file("Cannot write trajectory; maybe you are out of quota?");
            }
            gmx_fio_check_file_position(of->fp_trn);
        }
     }

}

