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

#ifdef GMX_CRAY_XT3
#include<catamount/dclock.h>
#endif


#include <stdio.h>
#include <time.h>
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

#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "gmx_wallcycle.h"

#ifdef GMX_MPI
#include "mpi.h"
#endif
#include "qmmm.h"


#define difftime(end,start) ((double)(end)-(double)(start))

void print_time(FILE *out,time_t start,int step,t_inputrec *ir)
{
  static real time_per_step;
  static time_t end;
  time_t finish;
  double dt;
  char buf[48];

  if (!gmx_parallel_env)
    fprintf(out,"\r");
  fprintf(out,"step %d",step);
  if ((step >= ir->nstlist)) {
    if ((ir->nstlist == 0) || ((step % ir->nstlist) == 0)) {
      /* We have done a full cycle let's update time_per_step */
      end=time(NULL);
      dt=difftime(end,start);
      time_per_step=dt/(step - ir->init_step + 1);
    }
    dt=(ir->nsteps + ir->init_step - step)*time_per_step;

    if (dt >= 300) {    
      finish = end+(time_t)dt;
      sprintf(buf,"%s",ctime(&finish));
      buf[strlen(buf)-1]='\0';
      fprintf(out,", will finish %s",buf);
    }
    else
      fprintf(out,", remaining runtime: %5d s          ",(int)dt);
  }
  if (gmx_parallel_env)
    fprintf(out,"\n");

  fflush(out);
}

time_t print_date_and_time(FILE *fplog,int nodeid,char *title)
{
  int i;
  char *ts,time_string[STRLEN];
  time_t now;

  now=time(NULL);
  ts=ctime(&now);
  for (i=0; ts[i]>=' '; i++) time_string[i]=ts[i];
  time_string[i]='\0';
  if (fplog)
    fprintf(fplog,"%s on node %d %s\n",title,nodeid,time_string);

  return now;
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

static void reset_energies(t_grpopts *opts,
			   t_forcerec *fr,bool bNS,
			   gmx_enerdata_t *enerd,
			   bool bMaster)
{
  bool bKeepLR;
  int  i,j;
  
  /* First reset all energy components, except for the long range terms
   * on the master at non neighbor search steps, since the long range
   * terms have already been summed at the last neighbor search step.
   */
  bKeepLR = (fr->bTwinRange && !bNS);
  for(i=0; (i<egNR); i++) {
    if (!(bKeepLR && bMaster && (i == egCOULLR || i == egLJLR))) {
      for(j=0; (j<enerd->grpp.nener); j++)
	enerd->grpp.ener[i][j] = 0.0;
    }
  }
  if (!bKeepLR) {
    enerd->dvdl_lr = 0.0;
  }

  /* Normal potential energy components */
  for(i=0; (i<=F_EPOT); i++) {
    enerd->term[i] = 0.0;
  }
  /* Initialize the dVdlambda term with the long range contribution */
  enerd->term[F_DVDL]     = enerd->dvdl_lr;
  enerd->term[F_DKDL]     = 0.0;
  enerd->term[F_DGDL_CON] = 0.0;
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
		      t_cosines Ex[],t_cosines Et[],real t)
{
  rvec Ext;
  real t0;
  int  i,m;
  
  for(m=0; (m<DIM); m++) {
    if (Et[m].n) {
      if (Et[m].n == 3) {
	t0 = Et[m].a[1];
	Ext[m] = cos(Et[m].a[0]*(t-t0))*exp(-sqr(t-t0)/(2.0*sqr(Et[m].a[2])));
      }
      else
	Ext[m] = cos(Et[m].a[0]*t);
    }
    else
      Ext[m] = 1.0;
    if (Ex[m].n) {
      /* Convert the field strength from V/nm to MD-units */
      Ext[m] *= Ex[m].a[0]*FIELDFAC;
      for(i=start; (i<start+homenr); i++)
	f[i][m] += charge[i]*Ext[m];
    }
    else
      Ext[m] = 0;
  }
  if (fp) 
    fprintf(fp,"%10g  %10g  %10g  %10g #FIELD\n",t,
	    Ext[XX]/FIELDFAC,Ext[YY]/FIELDFAC,Ext[ZZ]/FIELDFAC);
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

static real sum_v(int n,real v[])
{
  real t;
  int  i;
  
  t = 0.0;
  for(i=0; (i<n); i++)
    t = t + v[i];
    
  return t;
}

static void sum_epot(t_grpopts *opts,gmx_enerdata_t *enerd)
{
  gmx_grppairener_t *grpp;
  real *epot;
  int i;
  
  grpp = &enerd->grpp;
  epot = enerd->term;

  /* Accumulate energies */
  epot[F_COUL_SR]  = sum_v(grpp->nener,grpp->ener[egCOULSR]);
  epot[F_LJ]       = sum_v(grpp->nener,grpp->ener[egLJSR]);
  epot[F_LJ14]     = sum_v(grpp->nener,grpp->ener[egLJ14]);
  epot[F_COUL14]   = sum_v(grpp->nener,grpp->ener[egCOUL14]);
  epot[F_COUL_LR]  = sum_v(grpp->nener,grpp->ener[egCOULLR]);
  epot[F_LJ_LR]    = sum_v(grpp->nener,grpp->ener[egLJLR]);
/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
  epot[F_BHAM]     = sum_v(grpp->nener,grpp->ener[egBHAMSR]);
  epot[F_BHAM_LR]  = sum_v(grpp->nener,grpp->ener[egBHAMLR]);

  epot[F_EPOT] = 0;
  for(i=0; (i<F_EPOT); i++)
    if (i != F_DISRESVIOL && i != F_ORIRESDEV && i != F_DIHRESVIOL)
      epot[F_EPOT] += epot[i];
}

static void print_large_forces(FILE *fp,t_mdatoms *md,t_commrec *cr,
			       int step,real pforce,rvec *x,rvec *f)
{
  int  i;
  real pf2,fn2;

  pf2 = sqr(pforce);
  for(i=md->start; i<md->start+md->homenr; i++) {
    fn2 = norm2(f[i]);
    if (fn2 >= pf2) {
      fprintf(fp,"step %d  atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
	      step,ddglatnr(cr->dd,i),x[i][XX],x[i][YY],x[i][ZZ],sqrt(fn2));
    }
  }
}

void do_force(FILE *fplog,t_commrec *cr,
	      t_inputrec *inputrec,
	      int step,t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	      gmx_localtop_t *top,
	      gmx_groups_t *groups,
	      matrix box,rvec x[],history_t *hist,
	      rvec f[],rvec buf[],
	      tensor vir_force,
	      t_mdatoms *mdatoms,
	      gmx_enerdata_t *enerd,t_fcdata *fcd,
	      real lambda,t_graph *graph,
	      t_forcerec *fr,gmx_vsite_t *vsite,rvec mu_tot,
	      real t,FILE *field,gmx_edsam_t ed,
	      int flags)
{
  static rvec box_size;
  int    cg0,cg1,i,j;
  int    start,homenr;
  static double mu[2*DIM]; 
  rvec   mu_tot_AB[2];
  bool   bSepDVDL,bStateChanged,bNS,bFillGrid,bCalcCGCM,bBS,bDoForces;
  matrix boxs;
  real   e,v,dvdl;
  t_pbc  pbc;
  float  cycles_ppdpme,cycles_pme,cycles_force;
  
  start  = mdatoms->start;
  homenr = mdatoms->homenr;

  bSepDVDL = (fr->bSepDVDL && do_per_step(step,inputrec->nstlog));

  clear_mat(vir_force);
  
  if (PARTDECOMP(cr)) {
    pd_cg_range(cr,&cg0,&cg1);
  } else {
    cg0 = 0;
    if (DOMAINDECOMP(cr))
      cg1 = cr->dd->ncg_tot;
    else
      cg1 = top->cgs.nr;
    if (fr->n_tpi > 0)
      cg1--;
  }

  bStateChanged = (flags & GMX_FORCE_STATECHANGED);
  bNS           = (flags & GMX_FORCE_NS);
  bFillGrid     = (bNS && bStateChanged);
  bCalcCGCM     = (bFillGrid && !DOMAINDECOMP(cr));
  bDoForces     = (flags & GMX_FORCE_FORCES);

  if (bStateChanged) {
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
    if (DYNAMIC_BOX(*inputrec) && bStateChanged)
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

    gmx_pme_send_x(cr,bBS ? boxs : box,x,mdatoms->nChargePerturbed,lambda);

    GMX_MPE_LOG(ev_send_coordinates_finish);
    wallcycle_stop(wcycle,ewcPP_PMESENDX);
  }
#endif /* GMX_MPI */

  /* Communicate coordinates and sum dipole if necessary */
  if (PAR(cr)) {
    wallcycle_start(wcycle,ewcMOVEX);
    if (DOMAINDECOMP(cr)) {
      dd_move_x(cr->dd,box,x,buf);
    } else {
      move_x(fplog,cr,GMX_LEFT,GMX_RIGHT,x,nrnb);
    }
    /* When we don't need the total dipole we sum it in global_stat */
    if (NEED_MUTOT(*inputrec))
      gmx_sumd(2*DIM,mu,cr);
    wallcycle_stop(wcycle,ewcMOVEX);
  }
  for(i=0; i<2; i++)
    for(j=0;j<DIM;j++)
      mu_tot_AB[i][j] = mu[i*DIM + j];
  if (fr->efep == efepNO)
    copy_rvec(mu_tot_AB[0],mu_tot);
  else
    for(j=0; j<DIM; j++)
      mu_tot[j] = (1.0 - lambda)*mu_tot_AB[0][j] + lambda*mu_tot_AB[1][j];

  /* Reset energies */
  reset_energies(&(inputrec->opts),fr,bNS,enerd,MASTER(cr));    
  if (bNS) {
    wallcycle_start(wcycle,ewcNS);
    
    if (graph && bStateChanged)
      /* Calculate intramolecular shift vectors to make molecules whole */
      mk_mshift(fplog,graph,fr->ePBC,box,x);

    /* Reset long range forces if necessary */
    if (fr->bTwinRange) {
      clear_rvecs(fr->f_twin_n,fr->f_twin);
      clear_rvecs(SHIFTS,fr->fshift_twin);
    }
    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    dvdl = 0; 
    ns(fplog,fr,x,f,box,groups,&(inputrec->opts),top,mdatoms,
       cr,nrnb,step,lambda,&dvdl,&enerd->grpp,bFillGrid,bDoForces);
    if (bSepDVDL)
      fprintf(fplog,sepdvdlformat,"LR non-bonded",0,dvdl);
    enerd->dvdl_lr       = dvdl;
    enerd->term[F_DVDL] += dvdl;

    wallcycle_stop(wcycle,ewcNS);
  }
  
  if (DOMAINDECOMP(cr)) {
    if (!(cr->duty & DUTY_PME)) {
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

  if (bDoForces) {
      /* Reset PME/Ewald forces if necessary */
    if (fr->bF_NoVirSum) 
    {
      GMX_BARRIER(cr->mpi_comm_mygroup);
      if (fr->bDomDec)
	clear_rvecs(fr->f_novirsum_n,fr->f_novirsum);
      else
	clear_rvecs(homenr,fr->f_novirsum+start);
      GMX_BARRIER(cr->mpi_comm_mygroup);
    }
    /* Copy long range forces into normal buffers */
    if (fr->bTwinRange) {
      for(i=0; i<fr->f_twin_n; i++)
	copy_rvec(fr->f_twin[i],f[i]);
      for(i=0; i<SHIFTS; i++)
	copy_rvec(fr->fshift_twin[i],fr->fshift[i]);
    } 
    else {
      if (DOMAINDECOMP(cr))
	clear_rvecs(cr->dd->nat_tot,f);
      else
	clear_rvecs(mdatoms->nr,f);
      clear_rvecs(SHIFTS,fr->fshift);
    }
    clear_rvec(fr->vir_diag_posres);
    GMX_BARRIER(cr->mpi_comm_mygroup);
  }
  if (inputrec->ePull == epullCONSTRAINT)
    clear_pull_forces(inputrec->pull);

  /* update QMMMrec, if necessary */
  if(fr->bQMMM)
    update_QMMMrec(cr,fr,x,mdatoms,box,top);

  if ((flags & GMX_FORCE_BONDED) && top->idef.il[F_POSRES].nr > 0) {
    /* Position restraints always require full pbc */
    set_pbc(&pbc,inputrec->ePBC,box);
    v = posres(top->idef.il[F_POSRES].nr,top->idef.il[F_POSRES].iatoms,
	       top->idef.iparams_posres,
	       (const rvec*)x,fr->f_novirsum,fr->vir_diag_posres,
	       inputrec->ePBC==epbcNONE ? NULL : &pbc,lambda,&dvdl,
	       fr->rc_scaling,fr->ePBC,fr->posres_com,fr->posres_comB);
    if (bSepDVDL) {
      fprintf(fplog,sepdvdlformat,
	      interaction_function[F_POSRES].longname,v,dvdl);
    }
    enerd->term[F_POSRES] += v;
    enerd->term[F_DVDL]   += dvdl;
    inc_nrnb(nrnb,eNR_POSRES,top->idef.il[F_POSRES].nr/2);
  }
  /* Compute the bonded and non-bonded forces */    
  do_force_lowlevel(fplog,step,fr,inputrec,&(top->idef),
		    cr,nrnb,wcycle,mdatoms,&(inputrec->opts),
		    x,hist,f,enerd,fcd,box,lambda,graph,&(top->excls),mu_tot_AB,
		    flags,&cycles_force);
  GMX_BARRIER(cr->mpi_comm_mygroup);

  if (ed) {
    do_flood(fplog,cr,x,f,ed,box,step);
  }
	
  if (DOMAINDECOMP(cr)) {
    dd_force_flop_stop(cr->dd,nrnb);
    if (wcycle)
      dd_cycles_add(cr->dd,cycles_force,ddCyclF);
  }
  
  if (bDoForces) {
    /* Compute forces due to electric field */
    calc_f_el(MASTER(cr) ? field : NULL,
	      start,homenr,mdatoms->chargeA,x,f,inputrec->ex,inputrec->et,t);
    
    /* When using PME/Ewald we compute the long range virial there.
     * otherwise we do it based on long range forces from twin range
     * cut-off based calculation (or not at all).
     */
    
    /* Communicate the forces */
    if (PAR(cr)) {
      wallcycle_start(wcycle,ewcMOVEF);
      if (DOMAINDECOMP(cr)) {
	dd_move_f(cr->dd,f,buf,fr->fshift);
	/* Position restraint do not introduce inter-cg forces */
	if (EEL_FULL(fr->eeltype) && cr->dd->n_intercg_excl)
	  dd_move_f(cr->dd,fr->f_novirsum,buf,NULL);
      } else {
	move_f(fplog,cr,GMX_LEFT,GMX_RIGHT,f,buf,nrnb);
      }
      wallcycle_stop(wcycle,ewcMOVEF);
    }
  }

  if (bDoForces) {
    if (vsite) {
      wallcycle_start(wcycle,ewcVSITESPREAD);
      spread_vsite_f(fplog,vsite,x,f,fr->fshift,nrnb,
		     &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
      wallcycle_stop(wcycle,ewcVSITESPREAD);
    }
    
    /* Calculation of the virial must be done after vsites! */
    calc_virial(fplog,mdatoms->start,mdatoms->homenr,x,f,
		vir_force,graph,box,nrnb,fr,inputrec->ePBC);
  }

  if (inputrec->ePull == epullUMBRELLA || inputrec->ePull == epullCONST_F) {
    /* Calculate the center of mass forces, this requires communication,
     * which is why pull_potential is called close to other communication.
     * The virial contribution is calculated directly,
     * which is why we call pull_potential after calc_virial.
     */
    set_pbc(&pbc,inputrec->ePBC,box);
    dvdl = 0; 
    enerd->term[F_COM_PULL] =
      pull_potential(inputrec->ePull,inputrec->pull,mdatoms,&pbc,
		     cr,t,lambda,x,f,vir_force,&dvdl);
    if (bSepDVDL)
      fprintf(fplog,sepdvdlformat,"Com pull",enerd->term[F_COM_PULL],dvdl);
    enerd->term[F_DVDL] += dvdl;
  }

  if (!(cr->duty & DUTY_PME)) {
    cycles_ppdpme = wallcycle_stop(wcycle,ewcPPDURINGPME);
    dd_cycles_add(cr->dd,cycles_ppdpme,ddCyclPPduringPME);
  }

#ifdef GMX_MPI
  if (PAR(cr) && !(cr->duty & DUTY_PME)) {
    /* In case of node-splitting, the PP nodes receive the long-range 
     * forces, virial and energy from the PME nodes here.
     */    
    wallcycle_start(wcycle,ewcPP_PMEWAITRECVF);
    dvdl = 0;
    gmx_pme_receive_f(cr,fr->f_novirsum,fr->vir_el_recip,&e,&dvdl,
		      &cycles_pme);
    if (bSepDVDL)
      fprintf(fplog,sepdvdlformat,"PME mesh",e,dvdl);
    enerd->term[F_COUL_RECIP] += e;
    enerd->term[F_DVDL] += dvdl;
    if (wcycle)
      dd_cycles_add(cr->dd,cycles_pme,ddCyclPME);
    wallcycle_stop(wcycle,ewcPP_PMEWAITRECVF);
  }
#endif

  if (bDoForces && fr->bF_NoVirSum) {
    if (vsite) {
      /* Spread the mesh force on virtual sites to the other particles... 
       * This is parallellized. MPI communication is performed
       * if the constructing atoms aren't local.
       */
      wallcycle_start(wcycle,ewcVSITESPREAD);
      spread_vsite_f(fplog,vsite,x,fr->f_novirsum,NULL,nrnb,
		     &top->idef,fr->ePBC,fr->bMolPBC,graph,box,cr);
      wallcycle_stop(wcycle,ewcVSITESPREAD);
    }
    /* Now add the forces, this is local */
    if (fr->bDomDec) {
      sum_forces(0,fr->f_novirsum_n,f,fr->f_novirsum);
    } else {
      sum_forces(start,start+homenr,f,fr->f_novirsum);
    }
    if (EEL_FULL(fr->eeltype)) {
      /* Add the mesh contribution to the virial */
      m_add(vir_force,fr->vir_el_recip,vir_force);
    }
    if (debug)
      pr_rvecs(debug,0,"vir_force",vir_force,DIM);
  }

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(inputrec->opts),enerd);

  if (fr->print_force >= 0 && bDoForces)
    print_large_forces(stderr,mdatoms,cr,step,fr->print_force,x,f);
}

#ifdef NO_CLOCK 
#define clock() -1
#endif
static double runtime=0;
#ifdef GMX_CRAY_XT3
static double cprev;

void start_time(void)
{
  cprev   = dclock();
  runtime = 0.0;
}

void update_time(void)
{
  double c;

  c        = dclock();
  runtime += (c-cprev);
  cprev    = c;
}
#else
static clock_t cprev;

void start_time(void)
{
  cprev   = clock();
  runtime = 0.0;
}

void update_time(void)
{
  clock_t c;
  
  c        = clock();
  runtime += (c-cprev)/(double)CLOCKS_PER_SEC;
  cprev    = c;
}
#endif
double node_time(void)
{
  return runtime;
}

void do_shakefirst(FILE *fplog,gmx_constr_t constr,
		   t_inputrec *inputrec,t_mdatoms *md,
		   t_state *state,rvec buf[],rvec f[],
		   t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
		   t_forcerec *fr,t_idef *idef)
{
  int    i,m,start,end,step;
  double mass,tmass,vcm[4];
  real   dt=inputrec->delta_t;
  real   dvdlambda;

  start = md->start;
  end   = md->homenr + start;
  if (debug)
    fprintf(debug,"vcm: start=%d, homenr=%d, end=%d\n",
	    start,md->homenr,end);
  /* Do a first SHAKE to reset particles... */
  step = inputrec->init_step;
  if (fplog)
    fprintf(fplog,"\nConstraining the starting coordinates (step %d)\n",step);
  dvdlambda = 0;
  constrain(NULL,TRUE,FALSE,constr,idef,
	    inputrec,cr,step,0,md,
	    state->x,state->x,NULL,
	    state->box,state->lambda,&dvdlambda,
	    NULL,NULL,nrnb,econqCoord);

  if (EI_STATE_VELOCITY(inputrec->eI)) {
    for(i=start; (i<end); i++) {
      for(m=0; (m<DIM); m++) {
	/* Reverse the velocity */
	state->v[i][m] = -state->v[i][m];
	/* Store the position at t-dt in buf */
	buf[i][m] = state->x[i][m] + dt*state->v[i][m];
      }
    }
    
    /* Shake the positions at t=-dt with the positions at t=0
     * as reference coordinates.
     */
    if (fplog)
      fprintf(fplog,"\nConstraining the coordinates at t0-dt (step %d)\n",
	      step);
    dvdlambda = 0;
    constrain(NULL,TRUE,FALSE,constr,idef,
	      inputrec,cr,step,-1,md,
	      state->x,buf,NULL,
	      state->box,state->lambda,&dvdlambda,
	      state->v,NULL,nrnb,econqCoord);
    
    for(m=0; (m<4); m++)
      vcm[m] = 0;
    for(i=start; i<end; i++) {
      mass = md->massT[i];
      for(m=0; m<DIM; m++) {
	/* Re-reverse the velocities */
	state->v[i][m] = -state->v[i][m];
	vcm[m] += state->v[i][m]*mass;
      }
      vcm[3] += mass;
    }
    
    if (inputrec->nstcomm != 0 || debug) {
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
      if (inputrec->nstcomm != 0) {
	/* Now we have the velocity of center of mass, let's remove it */
	for(i=start; (i<end); i++) {
	  for(m=0; (m<DIM); m++)
	    state->v[i][m] -= vcm[m];
	}
      }
    }
  }
}

static void calc_enervirdiff(FILE *fplog,int eDispCorr,t_forcerec *fr)
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

void calc_dispcorr(FILE *fplog,t_inputrec *ir,t_forcerec *fr,int step,
		   int natoms,matrix box,real lambda,
		   tensor pres,tensor virial,real ener[])
{
  static bool bFirst=TRUE;
  bool bCorrAll,bCorrPres;
  real dvdlambda,invvol,dens,ninter,avcsix,avctwelve,enerdiff,svir=0,spres=0;
  int  m;

  ener[F_DISPCORR] = 0.0;

  if (ir->eDispCorr != edispcNO) {
    bCorrAll  = (ir->eDispCorr == edispcAllEner ||
		 ir->eDispCorr == edispcAllEnerPres);
    bCorrPres = (ir->eDispCorr == edispcEnerPres ||
		 ir->eDispCorr == edispcAllEnerPres);

    if (bFirst)
      calc_enervirdiff(fplog,ir->eDispCorr,fr);
    
    invvol = 1/det(box);
    if (fr->n_tpi) {
      /* Only correct for the interactions with the inserted molecule */
      dens = (natoms - fr->n_tpi)*invvol;
      ninter = fr->n_tpi;
    } else {
      dens = natoms*invvol;
      ninter = 0.5*natoms;
    }

    if (ir->efep == efepNO) {
      avcsix    = fr->avcsix[0];
      avctwelve = fr->avctwelve[0];
    } else {
      avcsix    = (1 - lambda)*fr->avcsix[0]    + lambda*fr->avcsix[1];
      avctwelve = (1 - lambda)*fr->avctwelve[0] + lambda*fr->avctwelve[1];
    }
    
    enerdiff = ninter*(dens*fr->enerdiffsix - fr->enershiftsix);
    ener[F_DISPCORR] += avcsix*enerdiff;
    dvdlambda = 0.0;
    if (ir->efep != efepNO)
      dvdlambda += (fr->avcsix[1] - fr->avcsix[0])*enerdiff;

    if (bCorrAll) {
      enerdiff = ninter*(dens*fr->enerdifftwelve - fr->enershifttwelve);
      ener[F_DISPCORR] += avctwelve*enerdiff;
      if (fr->efep != efepNO)
	dvdlambda += (fr->avctwelve[1] - fr->avctwelve[0])*enerdiff;
    }

    if (bCorrPres) {
      svir = ninter*dens*avcsix*fr->virdiffsix/3.0;
      if (ir->eDispCorr == edispcAllEnerPres)
	svir += ninter*dens*avctwelve*fr->virdifftwelve/3.0;
      
      /* The factor 2 is because of the Gromacs virial definition */
      spres = -2.0*invvol*svir*PRESFAC;
      
      for(m=0; m<DIM; m++) {
	virial[m][m] += svir;
	pres[m][m] += spres;
      }
      ener[F_PRES] += spres;
    }
    
    if (bFirst && fplog) {
      if (bCorrAll)
	fprintf(fplog,"Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
		avcsix,avctwelve);
      else
	fprintf(fplog,"Long Range LJ corr.: <C6> %10.4e\n",avcsix);
      
      if (bCorrPres)
	fprintf(fplog,
		"Long Range LJ corr.: Epot %10g, Pres: %10g, Vir: %10g\n",
                ener[F_DISPCORR],spres,svir);
      else
	fprintf(fplog,"Long Range LJ corr.: Epot %10g\n",ener[F_DISPCORR]);
    }
    bFirst = FALSE;

    if (fr->bSepDVDL && do_per_step(step,ir->nstlog))
      fprintf(fplog,sepdvdlformat,"Dispersion correction",
	      ener[F_DISPCORR],dvdlambda);
    
    ener[F_EPOT] += ener[F_DISPCORR];
    if (fr->efep != efepNO)
      ener[F_DVDL] += dvdlambda;
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
			    bool bFirst)
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

void finish_run(FILE *fplog,t_commrec *cr,char *confout,
		t_inputrec *inputrec,
		t_nrnb nrnb[],gmx_wallcycle_t wcycle,
		double nodetime,double realtime,int nsteps_done,
		bool bWriteStat)
{
  int    i,j;
  t_nrnb *nrnb_all=NULL,ntot;
  real   delta_t;
  double nbfs,mflop;
  double cycles[ewcNR];
#ifdef GMX_MPI
  int    sender;
  double nrnb_buf[4];
  MPI_Status status;
#endif

  wallcycle_sum(cr,wcycle,cycles);

  if (cr->nnodes > 1) {
    if (SIMMASTER(cr))
      snew(nrnb_all,cr->nnodes);
#ifdef GMX_MPI
    MPI_Gather(nrnb,sizeof(t_nrnb),MPI_BYTE,
	       nrnb_all,sizeof(t_nrnb),MPI_BYTE,
	       0,cr->mpi_comm_mysim);
#endif  
  } else {
    nrnb_all = nrnb;
  }
    
  if (SIMMASTER(cr)) {
    for(i=0; (i<eNRNB); i++)
      ntot.n[i]=0;
    for(i=0; (i<cr->nnodes); i++)
      for(j=0; (j<eNRNB); j++)
	ntot.n[j] += nrnb_all[i].n[j];

    print_flop(fplog,&ntot,&nbfs,&mflop);
    if (nrnb_all) {
      sfree(nrnb_all);
    }
  }

  if ((cr->duty & DUTY_PP) && DOMAINDECOMP(cr)) {
    print_dd_statistics(cr,inputrec,fplog);
  }

  if (SIMMASTER(cr)) {
    if (PARTDECOMP(cr)) {
      pr_load(fplog,cr,nrnb_all);
    }

    wallcycle_print(fplog,cr->nnodes,cr->npmenodes,realtime,wcycle,cycles);

    if (EI_DYNAMICS(inputrec->eI)) {
      delta_t = inputrec->delta_t;
    } else {
      delta_t = 0;
    }
    
    if (fplog) {
      print_perf(fplog,nodetime,realtime,cr->nnodes-cr->npmenodes,
		 nsteps_done,delta_t,nbfs,mflop);
    }
    if (bWriteStat) {
      print_perf(stderr,nodetime,realtime,cr->nnodes-cr->npmenodes,
		 nsteps_done,delta_t,nbfs,mflop);
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

void init_md(FILE *fplog,
			 t_commrec *cr,t_inputrec *ir,real *t,real *t0,
			 real *lambda,real *lam0,
			 t_nrnb *nrnb,gmx_mtop_t *mtop,
			 gmx_stochd_t *sd,
			 int nfile,t_filenm fnm[],
			 int *fp_trn,int *fp_xtc,int *fp_ene,char **fn_cpt,
			 FILE **fp_dgdl,FILE **fp_field,
			 t_mdebin **mdebin,
			 tensor force_vir,tensor shake_vir,rvec mu_tot,
			 bool *bNEMD,bool *bSimAnn,t_vcm **vcm, unsigned long Flags)
{
  int  i,j,n;
  real tmpt,mod;
  char filemode[2];

  sprintf(filemode, (Flags & MD_APPENDFILES) ? "a" : "w");
	
  /* Initial values */
  *t = *t0       = ir->init_t;
  if (ir->efep != efepNO) {
    *lam0 = ir->init_lambda;
    *lambda = *lam0 + ir->init_step*ir->delta_lambda;
  }
  else {
    *lambda = *lam0   = 0.0;
  } 
 
  *bSimAnn=FALSE;
  for(i=0;i<ir->opts.ngtc;i++) {
    /* set bSimAnn if any group is being annealed */
    if(ir->opts.annealing[i]!=eannNO)
      *bSimAnn = TRUE;
  }
  if (*bSimAnn) {
    update_annealing_target_temp(&(ir->opts),ir->init_t);
  }

  *bNEMD = (ir->opts.ngacc > 1) || (norm(ir->opts.acc[0]) > 0);
  
  if (sd && (ir->eI == eiBD || EI_SD(ir->eI) || ir->etc == etcVRESCALE)) {
    *sd = init_stochd(fplog,ir);
  }

  if (vcm) {
    *vcm = init_vcm(fplog,&mtop->groups,ir);
  }
   
  if (EI_DYNAMICS(ir->eI) && !(Flags & MD_APPENDFILES)) {
    if (ir->etc == etcBERENDSEN) {
      please_cite(fplog,"Berendsen84a");
    }
    if (ir->etc == etcVRESCALE) {
      please_cite(fplog,"Bussi2007a");
    }
  }
 
  init_nrnb(nrnb);

  if (nfile != -1)
  {
	  *fp_trn = -1;
	  *fp_ene = -1;
	  *fp_xtc = -1;
	  
	  if (MASTER(cr)) 
	  {
		  *fp_trn = open_trn(ftp2fn(efTRN,nfile,fnm), filemode);
		  if (ir->nstxtcout > 0)
		  {
			  *fp_xtc = open_xtc(ftp2fn(efXTC,nfile,fnm), filemode);
		  }
		  *fp_ene = open_enx(ftp2fn(efENX,nfile,fnm), filemode);
		  *fn_cpt = opt2fn("-cpo",nfile,fnm);
		  
		  if ((fp_dgdl != NULL) && ir->efep!=efepNO)
		  {
			  if(Flags & MD_APPENDFILES)
			  {
				  *fp_dgdl= gmx_fio_fopen(opt2fn("-dgdl",nfile,fnm),filemode);
			  }
			  else
			  {
				  *fp_dgdl = xvgropen(opt2fn("-dgdl",nfile,fnm),
									  "dG/d\\8l\\4","Time (ps)",
									  "dG/d\\8l\\4 (kJ mol\\S-1\\N [\\8l\\4]\\S-1\\N)");
			  }
		  }
		  
		  if ((fp_field != NULL) && (ir->ex[XX].n || ir->ex[YY].n ||ir->ex[ZZ].n))
		  {
			  if(Flags & MD_APPENDFILES)
			  {
				  *fp_field = gmx_fio_fopen(opt2fn("-field",nfile,fnm),filemode);
			  }
			  else
			  {				  
				  *fp_field = xvgropen(opt2fn("-field",nfile,fnm),
									   "Applied electric field","Time (ps)",
									   "E (V/nm)");
			  }
		  }
	  }
	  *mdebin = init_mdebin( (Flags & MD_APPENDFILES) ? -1 : *fp_ene,mtop,ir);
  }

  /* Initiate variables */  
  clear_mat(force_vir);
  clear_mat(shake_vir);
  clear_rvec(mu_tot);

  debug_gmx();
}

