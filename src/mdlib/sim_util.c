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
  fprintf(out,"step %d",step - ir->init_step);
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

static void reset_energies(t_grpopts *opts,t_groups *grp,
			   t_forcerec *fr,bool bNS,real epot[])
{
  int   i,j;
  
  /* First reset all energy components but the Long Range, except in
   * some special cases.
   */
  for(i=0; (i<egNR); i++)
    if (((i != egCOULLR) && (i != egLJLR)) ||
	(fr->bTwinRange && bNS) || (!fr->bTwinRange))
      for(j=0; (j<grp->estat.nn); j++)
	grp->estat.ee[i][j]=0.0;
  
  /* Normal potential energy components */
  for(i=0; (i<=F_EPOT); i++)
    epot[i] = 0.0;
  epot[F_DVDL]    = 0.0;
  epot[F_DVDLKIN] = 0.0;
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

void do_force(FILE *fplog,t_commrec *cr,
	      t_inputrec *inputrec,
	      int step,t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	      t_topology *top,t_groups *grps,
	      matrix box,rvec x[],rvec f[],rvec buf[],
	      t_mdatoms *mdatoms,real ener[],t_fcdata *fcd,
	      real lambda,t_graph *graph,
	      bool bStateChanged,bool bNS,bool bNBFonly,bool bDoForces,
	      t_forcerec *fr,rvec mu_tot,
	      bool bGatherOnly,real t,FILE *field,t_edsamyn *edyn)
{
  static rvec box_size;
  static real dvdl_lr = 0;
  int    cg0,cg1,i,j;
  int    start,homenr;
  static double mu[2*DIM]; 
  rvec   mu_tot_AB[2];
  bool   bFillGrid,bCalcCGCM;
  real   e,d;
  float  pme_cycles;
  
  start  = mdatoms->start;
  homenr = mdatoms->homenr;
  
  if (PARTDECOMP(cr)) {
    pd_cg_range(cr,&cg0,&cg1);
  } else {
    cg0 = 0;
    if (DOMAINDECOMP(cr))
      cg1 = cr->dd->ncg_tot;
    else
      cg1 = top->blocks[ebCGS].nr;
  }

  bFillGrid = (bNS && bStateChanged);
  bCalcCGCM = (bFillGrid && !DOMAINDECOMP(cr));

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
      put_charge_groups_in_box(fplog,cg0,cg1,box,
			       &(top->blocks[ebCGS]),x,fr->cg_cm);
      inc_nrnb(nrnb,eNR_CGCM,homenr);
      inc_nrnb(nrnb,eNR_RESETX,cg1-cg0);
    } 
    else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph) {
      unshift_self(graph,box,x);
    }
  } 
  else if (bCalcCGCM) {
    calc_cgcm(fplog,cg0,cg1,&(top->blocks[ebCGS]),x,fr->cg_cm);
    inc_nrnb(nrnb,eNR_CGCM,homenr);
  }
  
  if (bCalcCGCM) {
    if (PAR(cr)) {
      move_cgcm(fplog,cr,fr->cg_cm);
    }
    if (gmx_debug_at)
      pr_rvecs(debug,0,"cgcm",fr->cg_cm,top->blocks[ebCGS].nr);
  }

#ifdef GMX_MPI
  /* send particle coordinates and charges to the pme nodes if necessary */
  if (!(cr->duty & DUTY_PME))
  { 
    GMX_MPE_LOG(ev_send_coordinates_start);
    
    if (graph) {
      GMX_MPE_LOG(ev_shift_start);
      shift_self(graph,box,x);
      GMX_MPE_LOG(ev_shift_finish);
    }

    gmx_pme_send_x_q(cr,box,x,NULL,NULL,
		     mdatoms->nChargePerturbed,lambda,step>=inputrec->nsteps);

    if (graph) {
      GMX_MPE_LOG(ev_shift_start);
      unshift_self(graph,box,x);
      GMX_MPE_LOG(ev_shift_finish);
    }    

    GMX_MPE_LOG(ev_send_coordinates_finish);
  }
#endif /* GMX_MPI */

  /* Communicate coordinates and sum dipole if necessary */
  if (PAR(cr)) {
    wallcycle_start(wcycle,ewcMOVEX);
    if (DOMAINDECOMP(cr)) {
      dd_move_x(cr->dd,box,x,buf);
    } else {
      move_x(fplog,cr,cr->left,cr->right,x,nrnb);
    }
    /* When we don't need the total dipole we sum it in global_stat */
    if (NEED_MUTOT(*inputrec))
      gmx_sumd(2*DIM,mu,cr);
    wallcycle_stop(wcycle,ewcMOVEX);
    if (DOMAINDECOMP(cr) && wcycle)
      dd_cycles_add(cr->dd,wallcycle_lastcycle(wcycle,ewcMOVEX),ddCyclMoveX);
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
  reset_energies(&(inputrec->opts),grps,fr,bNS,ener);    
  if (bNS) {
    wallcycle_start(wcycle,ewcNS);
    
    if (graph && bStateChanged)
      /* Calculate intramolecular shift vectors to make molecules whole */
      mk_mshift(fplog,graph,box,x);

    /* Reset long range forces if necessary */
    if (fr->bTwinRange) {
      clear_rvecs(fr->f_twin_n,fr->f_twin);
      clear_rvecs(SHIFTS,fr->fshift_twin);
    }
    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    dvdl_lr = 0; 

    ns(fplog,fr,x,f,box,grps,&(inputrec->opts),top,mdatoms,
       cr,nrnb,step,lambda,&dvdl_lr,bFillGrid,bDoForces);

    wallcycle_stop(wcycle,ewcNS);
  }
  
  wallcycle_start(wcycle,ewcFORCE);

  if (bDoForces) {
      /* Reset PME/Ewald forces if necessary */
    if (EEL_FULL(fr->eeltype)) 
    {
      GMX_BARRIER(cr->mpi_comm_mygroup);
      if (fr->bDomDec)
	clear_rvecs(fr->f_el_recip_n,fr->f_el_recip);
      else
	clear_rvecs(homenr,fr->f_el_recip+start);
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
	clear_rvecs(cr->dd->nat_tot_vsite,f);
      else
	clear_rvecs(top->atoms.nr,f);
      clear_rvecs(SHIFTS,fr->fshift);
    }
    GMX_BARRIER(cr->mpi_comm_mygroup);
  }

  /* update QMMMrec, if necessary */
  if(fr->bQMMM)
    update_QMMMrec(cr,fr,x,mdatoms,box,top);

  /* Compute the forces */    
  force(fplog,step,fr,inputrec,&(top->idef),cr,nrnb,wcycle,grps,mdatoms,
	top->atoms.grps[egcENER].nr,&(inputrec->opts),
	x,f,ener,fcd,box,lambda,graph,&(top->blocks[ebEXCLS]),
	bNBFonly,bDoForces,mu_tot_AB,bGatherOnly,edyn);
  GMX_BARRIER(cr->mpi_comm_mygroup);
	
  /* Take long range contribution to free energy into account */
  ener[F_DVDL] += dvdl_lr;

  wallcycle_stop(wcycle,ewcFORCE);
  if (DOMAINDECOMP(cr) && wcycle)
    dd_cycles_add(cr->dd,wallcycle_lastcycle(wcycle,ewcFORCE),ddCyclF);
  
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
	if (EEL_FULL(fr->eeltype) && cr->dd->n_intercg_excl)
	  dd_move_f(cr->dd,fr->f_el_recip,buf,NULL);
      } else {
	move_f(fplog,cr,cr->left,cr->right,f,buf,nrnb);
      }
      wallcycle_stop(wcycle,ewcMOVEF);
      if (DOMAINDECOMP(cr) && wcycle)
	dd_cycles_add(cr->dd,wallcycle_lastcycle(wcycle,ewcMOVEF),ddCyclMoveF);
      /* In case of node-splitting, the PP nodes receive the long-range 
       * forces, virial and energy from the PME nodes here 
       */    
#ifdef GMX_MPI       
      if (!(cr->duty & DUTY_PME))
      {
	d = 0;
	gmx_pme_receive_f(cr,fr->f_el_recip,fr->vir_el_recip,&e,&d,
			  &pme_cycles);
	if (fr->bSepDVDL && do_per_step(step,inputrec->nstlog))
	  fprintf(fplog,sepdvdlformat,"PME mesh",e,d);
        ener[F_COUL_RECIP] += e;
	ener[F_DVDL] += d;
	if (wcycle)
	  dd_cycles_add(cr->dd,pme_cycles,ddCyclPME);
      }
#endif
    }
  }
}

void sum_lrforces(rvec f[],t_forcerec *fr,int start,int homenr)
{
  /* Now add the forces from the PME calculation. Since this only produces
   * forces on the local atoms, this can be safely done after the
   * communication step.
   */
  if (EEL_FULL(fr->eeltype)) {
    if (fr->bDomDec)
      sum_forces(0,fr->f_el_recip_n,f,fr->f_el_recip);
    else
      sum_forces(start,start+homenr,f,fr->f_el_recip);
  }
}

void calc_virial(FILE *fplog,int start,int homenr,rvec x[],rvec f[],
		 tensor vir_part,tensor vir_el_recip,
		 t_graph *graph,matrix box,
		 t_nrnb *nrnb,const t_forcerec *fr)
{
  int i,j;
  tensor virtest;

  /* The short-range virial from surrounding boxes */
  clear_mat(vir_part);
  calc_vir(fplog,SHIFTS,fr->shift_vec,fr->fshift,vir_part);
  inc_nrnb(nrnb,eNR_VIRIAL,SHIFTS);
  
  /* Calculate partial virial, for local atoms only, based on short range. 
   * Total virial is computed in global_stat, called from do_md 
   */
  f_calc_vir(fplog,start,start+homenr,x,f,vir_part,graph,box);
  inc_nrnb(nrnb,eNR_VIRIAL,homenr);

  /* Add up virial if necessary */  
  if (EEL_FULL(fr->eeltype) && (fr->eeltype != eelPPPM)) {
    
    for(i=0; (i<DIM); i++) 
      for(j=0; (j<DIM); j++) 
	vir_part[i][j] += vir_el_recip[i][j];
  }
  if (debug)
    pr_rvecs(debug,0,"vir_part",vir_part,DIM);
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

void do_shakefirst(FILE *fplog,real ener[],
		   t_inputrec *inputrec,t_mdatoms *md,
		   t_state *state,rvec buf[],rvec f[],
		   t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
		   t_groups *grps,t_forcerec *fr,t_topology *top,
		   t_edsamyn *edyn,t_pull *pulldata)
{
  int    i,m,start,end,step;
  tensor shake_vir,pres;
  double mass,tmass,vcm[4];
  real   dt=inputrec->delta_t;
  real   dt_1;
  rvec   *xptr;

  if (count_constraints(top,cr)) {
    start = md->start;
    end   = md->homenr + start;
    if (debug)
      fprintf(debug,"vcm: start=%d, homenr=%d, end=%d\n",
	      start,md->homenr,end);
    /* Do a first SHAKE to reset particles... */
    step = -2;
    fprintf(fplog,"\nConstraining the starting coordinates (step %d)\n",step);
    clear_mat(shake_vir);
    clear_mat(pres);
    update(step,&ener[F_DVDL],inputrec,md,state,graph,
	   NULL,top,grps,shake_vir,cr,nrnb,
	   edyn,pulldata,FALSE,FALSE,FALSE,state->x,pres);
    /* Compute coordinates at t=-dt, store them in buf */
    /* for(i=0; (i<nsb->natoms); i++) {*/
    for(i=start; (i<end); i++) {
      for(m=0; (m<DIM); m++) {
	buf[i][m] = state->x[i][m] - dt*state->v[i][m];
      }
    }
    
    /* Shake the positions at t=-dt with the positions at t=0
     * as reference coordinates.
     */
    step = -1;
    fprintf(fplog,"\nConstraining the coordinates at t0-dt (step %d)\n",step);
    clear_mat(shake_vir);
    clear_mat(pres);
    update(step,&ener[F_DVDL],inputrec,md,state,graph,
	   NULL,top,grps,shake_vir,cr,nrnb,
	   edyn,pulldata,FALSE,FALSE,FALSE,buf,pres);

    /* Compute the velocities at t=-dt/2 using the coordinates at
     * t=-dt and t=0
     * Compute velocity of center of mass and total mass
     */
    for(m=0; (m<4); m++)
      vcm[m] = 0;
    dt_1=1.0/dt;
    for(i=start; (i<end); i++) {
      /*for(i=0; (i<nsb->natoms); i++) {*/
      mass = md->massT[i];
      for(m=0; (m<DIM); m++) {
	state->v[i][m] = (state->x[i][m] - buf[i][m])*dt_1;
	vcm[m] += state->v[i][m]*mass;
      }
      vcm[3] += mass;
    }
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
      if (fr->vdwtype == evdwUSER)
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

  bCorrAll  = (ir->eDispCorr == edispcAllEner ||
	       ir->eDispCorr == edispcAllEnerPres);
  bCorrPres = (ir->eDispCorr == edispcEnerPres ||
	       ir->eDispCorr == edispcAllEnerPres);
  ener[F_DISPCORR] = 0.0;
  ener[F_PRES]     = trace(pres)/3.0;
  dvdlambda        = 0.0;

  if (ir->eDispCorr != edispcNO) {
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
    
    if (bFirst) {
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
      bFirst = FALSE;
    }
  } 

  if (fr->bSepDVDL && do_per_step(step,ir->nstlog))
    fprintf(fplog,sepdvdlformat,"Dispersion correction",
	    ener[F_DISPCORR],dvdlambda);
  
  ener[F_EPOT] += ener[F_DISPCORR];
  if (fr->efep != efepNO)
    ener[F_DVDL] += dvdlambda;
}


void do_pbc_first(FILE *fplog,matrix box,t_forcerec *fr,
		  t_graph *graph,rvec x[])
{
  fprintf(fplog,"Removing pbc first time\n");
  calc_shifts(box,fr->shift_vec);
  if (graph) {
    mk_mshift(fplog,graph,box,x);
    if (gmx_debug_at)
      p_graph(debug,"do_pbc_first 1",graph);
    shift_self(graph,box,x);
    /* By doing an extra mk_mshift the molecules that are broken
     * because they were e.g. imported from another software
     * will be made whole again. Such are the healing powers
     * of GROMACS.
     */
    mk_mshift(fplog,graph,box,x);
    if (gmx_debug_at)
      p_graph(debug,"do_pbc_first 2",graph);
  }
  fprintf(fplog,"Done rmpbc\n");
}

void finish_run(FILE *fplog,t_commrec *cr,char *confout,
		t_topology *top,t_inputrec *inputrec,
		t_nrnb nrnb[],gmx_wallcycle_t wcycle,
		double nodetime,double realtime,int step,
		bool bWriteStat)
{
  int    i,j;
  t_nrnb *nrnb_all=NULL,ntot;
  real   runtime;
  double cycles[ewcNR];
#ifdef GMX_MPI
  int    sender;
  double nrnb_buf[4];
  MPI_Status status;
#endif

  wallcycle_sum(cr,wcycle,cycles);

  if (cr->nnodes > 1) {
    if (MASTER(cr))
      snew(nrnb_all,cr->nnodes);
#ifdef GMX_MPI
    MPI_Gather(nrnb,sizeof(t_nrnb),MPI_BYTE,
	       nrnb_all,sizeof(t_nrnb),MPI_BYTE,
	       0,cr->mpi_comm_mysim);
#endif  
  } else {
    nrnb_all = nrnb;
  }
    
  if (MASTER(cr)) {
    for(i=0; (i<eNRNB); i++)
      ntot.n[i]=0;
    for(i=0; (i<cr->nnodes); i++)
      for(j=0; (j<eNRNB); j++)
	ntot.n[j] += nrnb_all[i].n[j];
  }
  if (MASTER(cr)) {
    runtime=inputrec->nsteps*inputrec->delta_t;
    if (bWriteStat) {
      fprintf(stderr,"\n\n");
      wallcycle_print(stderr,cr->nnodes,cr->npmenodes,realtime,wcycle,cycles);
      print_perf(stderr,nodetime,realtime,runtime,&ntot,
		 cr->nnodes-cr->npmenodes);
    }
    wallcycle_print(fplog,cr->nnodes,cr->npmenodes,realtime,wcycle,cycles);
    print_perf(fplog,nodetime,realtime,runtime,&ntot,cr->nnodes-cr->npmenodes);
    if (PARTDECOMP(cr))
      pr_load(fplog,cr,nrnb_all);
    if (cr->nnodes > 1)
      sfree(nrnb_all);
  }
}

void init_md(t_commrec *cr,t_inputrec *ir,real *t,real *t0,
	     real *lambda,real *lam0,
	     t_nrnb *nrnb,t_topology *top,
	     int nfile,t_filenm fnm[],
	     int *fp_trn,int *fp_xtc,int *fp_ene,
	     FILE **fp_dgdl,FILE **fp_field,t_mdebin **mdebin,t_groups *grps,
	     tensor force_vir,tensor shake_vir,rvec mu_tot,
	     bool *bNEMD,bool *bSimAnn,t_vcm **vcm)
{
  int  i,j,n;
  real tmpt,mod;

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
  if(*bSimAnn) 
    update_annealing_target_temp(&(ir->opts),ir->init_t); 

  init_nrnb(nrnb);
  
  if (nfile != -1) {
    if (MASTER(cr)) {
      *fp_trn = open_trn(ftp2fn(efTRN,nfile,fnm),"w");
      if (ir->nstxtcout > 0)
	*fp_xtc = open_xtc(ftp2fn(efXTC,nfile,fnm),"w");
      *fp_ene = open_enx(ftp2fn(efENX,nfile,fnm),"w");
      if ((fp_dgdl != NULL) && ir->efep!=efepNO)
	*fp_dgdl =
	  xvgropen(opt2fn("-dgdl",nfile,fnm),
		   "dG/d\\8l\\4","Time (ps)",
		   "dG/d\\8l\\4 (kJ mol\\S-1\\N [\\8l\\4]\\S-1\\N)");
      if ((fp_field != NULL) && (ir->ex[XX].n || ir->ex[YY].n ||ir->ex[ZZ].n))
	*fp_field = xvgropen(opt2fn("-field",nfile,fnm),
			     "Applied electric field","Time (ps)",
			     "E (V/nm)");
    } else
      *fp_ene = -1;

    *mdebin = init_mdebin(*fp_ene,grps,&(top->atoms),&(top->idef),ir,cr);
  }
  
  /* Initiate variables */  
  clear_mat(force_vir);
  clear_mat(shake_vir);
  clear_rvec(mu_tot);

  *vcm = init_vcm(stdlog,&top->atoms,ir->nstcomm,ir->comm_mode);
    
  debug_gmx();

  *bNEMD = (ir->opts.ngacc > 1) || (norm(ir->opts.acc[0]) > 0);

  if (ir->eI == eiSD)
    init_sd_consts(ir->opts.ngtc,ir->opts.tau_t,ir->delta_t);

}

