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

#define difftime(end,start) ((double)(end)-(double)(start))

void print_time(FILE *out,time_t start,int step,t_inputrec *ir)
{
  static real time_per_step;
  static time_t end;
  time_t finish;
  double dt;
  char buf[48];

  fprintf(out,"\rstep %d",step - ir->init_step);
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
      fprintf(out,", will finish at %s",buf);
    }
    else
      fprintf(out,", remaining runtime: %5d s               ",(int)dt);
  }
  if (gmx_parallel_env)
    fprintf(out,"\n");

  fflush(out);
}

time_t print_date_and_time(FILE *log,int nodeid,char *title)
{
  int i;
  char *ts,time_string[STRLEN];
  time_t now;

  now=time(NULL);
  ts=ctime(&now);
  for (i=0; ts[i]>=' '; i++) time_string[i]=ts[i];
  time_string[i]='\0';
  fprintf(log,"%s on node %d %s\n",title,nodeid,time_string);
  return now;
}

static void pr_commrec(FILE *log,t_commrec *cr)
{
  fprintf(log,"commrec: nodeid=%d, nnodes=%d, left=%d, right=%d, threadid=%d, nthreads=%d\n",
	  cr->nodeid,cr->nnodes,cr->left,cr->right,cr->threadid,cr->nthreads);
}

static void sum_forces(int start,int end,rvec f[],rvec flr[])
{
  int i;
  
  if (debug) {
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
    if (((i != egLR) && (i != egLJLR)) ||
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

void do_force(FILE *log,t_commrec *cr,t_commrec *mcr,
	      t_parm *parm,t_nsborder *nsb,tensor vir_part,tensor pme_vir,
	      int step,t_nrnb *nrnb,t_topology *top,t_groups *grps,
	      matrix box,rvec x[],rvec f[],rvec buf[],
	      t_mdatoms *mdatoms,real ener[],t_fcdata *fcd,bool bVerbose,
	      real lambda,t_graph *graph,
	      bool bNS,bool bNBFonly,t_forcerec *fr, rvec mu_tot,
	      bool bGatherOnly,real t,FILE *field)
{
  static rvec box_size;
  static real dvdl_lr = 0;
  int    cg0,cg1,i,j;
  int    start,homenr;
  static real mu_and_q[DIM+1]; 
  real   qsum;
  
  start  = START(nsb);
  homenr = HOMENR(nsb);
  cg0    = CG0(nsb);
  cg1    = CG1(nsb);
  
  update_forcerec(log,fr,box);

  /* Calculate total (local) dipole moment in a temporary common array. 
   * This makes it possible to sum them over nodes faster.
   */
  calc_mu_and_q(nsb,x,mdatoms->chargeT,mu_and_q,mu_and_q+DIM);

  if (fr->ePBC != epbcNONE) { 
    /* Compute shift vectors every step, because of pressure coupling! */
    if (parm->ir.epc != epcNO)
      calc_shifts(box,box_size,fr->shift_vec);
    
    if (bNS) { 
      put_charge_groups_in_box(log,cg0,cg1,box,box_size,
			       &(top->blocks[ebCGS]),x,fr->cg_cm);
      inc_nrnb(nrnb,eNR_RESETX,homenr);
    } 
    else if ((parm->ir.eI==eiSteep || parm->ir.eI==eiCG) && graph)
      unshift_self(graph,box,x);
    
  }
  else if (bNS)
    calc_cgcm(log,cg0,cg1,&(top->blocks[ebCGS]),x,fr->cg_cm);
 
  if (bNS) {
    inc_nrnb(nrnb,eNR_CGCM,cg1-cg0);
    if (PAR(cr))
      move_cgcm(log,cr,fr->cg_cm,nsb->workload);
    if (debug)
      pr_rvecs(debug,0,"cgcm",fr->cg_cm,nsb->cgtotal);
  }
  
  /* Communicate coordinates and sum dipole and net charge if necessary */
  if (PAR(cr)) {
    move_x(log,cr->left,cr->right,x,nsb,nrnb);
    gmx_sum(DIM+1,mu_and_q,cr);
  }
  for(i=0;i<DIM;i++)
    mu_tot[i]=mu_and_q[i];
  qsum=mu_and_q[DIM];
  
  /* Reset energies */
  reset_energies(&(parm->ir.opts),grps,fr,bNS,ener);    
  if (bNS) {
    if (fr->ePBC == epbcXYZ)
      /* Calculate intramolecular shift vectors to make molecules whole */
      mk_mshift(log,graph,box,x);
	       
    /* Reset long range forces if necessary */
    if (fr->bTwinRange) {
      clear_rvecs(nsb->natoms,fr->f_twin);
      clear_rvecs(SHIFTS,fr->fshift_twin);
    }
    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    dvdl_lr = 0; 

    ns(log,fr,x,f,box,grps,&(parm->ir.opts),top,mdatoms,
       cr,nrnb,nsb,step,lambda,&dvdl_lr);
  }
  /* Reset PME/Ewald forces if necessary */
  if (EEL_LR(fr->eeltype)) 
    clear_rvecs(homenr,fr->f_pme+start);
    
  /* Copy long range forces into normal buffers */
  if (fr->bTwinRange) {
    for(i=0; i<nsb->natoms; i++)
      copy_rvec(fr->f_twin[i],f[i]);
    for(i=0; i<SHIFTS; i++)
      copy_rvec(fr->fshift_twin[i],fr->fshift[i]);
  } 
  else {
    clear_rvecs(nsb->natoms,f);
    clear_rvecs(SHIFTS,fr->fshift);
  }
  
  /* Compute the forces */    
  force(log,step,fr,&(parm->ir),&(top->idef),nsb,cr,mcr,nrnb,grps,mdatoms,
	top->atoms.grps[egcENER].nr,&(parm->ir.opts),
	x,f,ener,fcd,bVerbose,box,lambda,graph,&(top->atoms.excl),
	bNBFonly,pme_vir,mu_tot,qsum,bGatherOnly);
	
  /* Take long range contribution to free energy into account */
  ener[F_DVDL] += dvdl_lr;
  
#ifdef DEBUG
  if (bNS)
    print_nrnb(log,nrnb);
#endif

  /* The short-range virial from surrounding boxes */
  clear_mat(vir_part);
  calc_vir(log,SHIFTS,fr->shift_vec,fr->fshift,vir_part);
  inc_nrnb(nrnb,eNR_VIRIAL,SHIFTS);

  if (debug) 
    pr_rvecs(debug,0,"vir_shifts",vir_part,DIM);

  /* Compute forces due to electric field */
  calc_f_el(MASTER(cr) ? field : NULL,
	    start,homenr,mdatoms->chargeT,x,f,parm->ir.ex,parm->ir.et,t);

  /* When using PME/Ewald we compute the long range virial (pme_vir) there.
   * otherwise we do it based on long range forces from twin range
   * cut-off based calculation (or not at all).
   */
  
  /* Communicate the forces */
  if (PAR(cr))
    move_f(log,cr->left,cr->right,f,buf,nsb,nrnb);
}

void sum_lrforces(rvec f[],t_forcerec *fr,int start,int homenr)
{
  /* Now add the forces from the PME calculation. Since this only produces
   * forces on the local atoms, this can be safely done after the
   * communication step.
   */
  if (EEL_LR(fr->eeltype))
    sum_forces(start,start+homenr,f,fr->f_pme);
}

void calc_virial(FILE *log,int start,int homenr,rvec x[],rvec f[],
		 tensor vir_part,tensor pme_vir,
		 t_graph *graph,matrix box,
		 t_nrnb *nrnb,t_forcerec *fr,bool bTweak)
{
  int i,j;
  tensor virtest;
  
  /* Now it is time for the short range virial. At this timepoint vir_part
   * already contains the virial from surrounding boxes.
   * Calculate partial virial, for local atoms only, based on short range. 
   * Total virial is computed in global_stat, called from do_md 
   */
  f_calc_vir(log,start,start+homenr,x,f,vir_part,graph,box);
  inc_nrnb(nrnb,eNR_VIRIAL,homenr);

  /* Add up the long range forces if necessary */
  /* if (!bTweak) {
    sum_lrforces(f,fr,start,homenr);
    }*/
  
  /* Add up virial if necessary */  
  if (EEL_LR(fr->eeltype) && (fr->eeltype != eelPPPM)) {
    if (debug && bTweak) {
      clear_mat(virtest);
      f_calc_vir(log,start,start+homenr,x,fr->f_pme,virtest,graph,box);
      pr_rvecs(debug,0,"virtest",virtest,DIM);
      pr_rvecs(debug,0,"pme_vir",pme_vir,DIM);
    }    
    /* PPPM virial sucks */
    if (!bTweak)
      for(i=0; (i<DIM); i++) 
	for(j=0; (j<DIM); j++) 
	  vir_part[i][j]+=pme_vir[i][j];
  }
  if (debug)
    pr_rvecs(debug,0,"vir_part",vir_part,DIM);
}

#ifdef NO_CLOCK 
#define clock() -1
#endif
static double runtime=0;
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

double node_time(void)
{
  return runtime;
}

void do_shakefirst(FILE *log,bool bTYZ,real ener[],
		   t_parm *parm,t_nsborder *nsb,t_mdatoms *md,
		   t_state *state,rvec vold[],rvec buf[],rvec f[],
		   t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
		   t_groups *grps,t_forcerec *fr,t_topology *top,
		   t_edsamyn *edyn,t_pull *pulldata)
{
  int    i,m,start,homenr,end,step;
  tensor shake_vir;
  double mass,tmass,vcm[4];
  real   dt=parm->ir.delta_t;
  real   dt_1;
  rvec   *xptr;

  if (count_constraints(top,cr)) {
    start  = START(nsb);
    homenr = HOMENR(nsb);
    end    = start+homenr;
    if (debug)
      fprintf(debug,"vcm: start=%d, homenr=%d, end=%d\n",start,homenr,end);
    /* Do a first SHAKE to reset particles... */
    step = -2;
    fprintf(log,"\nConstraining the starting coordinates (step %d)\n",step);
    clear_mat(shake_vir);
    update(nsb->natoms,start,homenr,step,&ener[F_DVDL],
	   parm,md,state,graph,
	   NULL,NULL,vold,top,grps,shake_vir,cr,nrnb,bTYZ,
	   edyn,pulldata,FALSE,FALSE,FALSE,state->x);
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
    fprintf(log,"\nConstraining the coordinates at t0-dt (step %d)\n",step);
    clear_mat(shake_vir);
    update(nsb->natoms,start,homenr,step,&ener[F_DVDL],
	   parm,md,state,graph,
	   NULL,NULL,vold,top,grps,shake_vir,cr,nrnb,bTYZ,
	   edyn,pulldata,FALSE,FALSE,FALSE,buf);

    /* Compute the velocities at t=-dt/2 using the coordinates at
     * t=-dt and t=0
     * Compute velocity of center of mass and total mass
     */
    for(m=0; (m<4); m++)
      vcm[m] = 0;
    dt_1=1.0/dt;
    for(i=start; (i<end); i++) {
      /*for(i=0; (i<nsb->natoms); i++) {*/
      mass = md->massA[i];
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
    /* Now we have the velocity of center of mass, let's remove it */
    for(i=start; (i<end); i++) {
      for(m=0; (m<DIM); m++)
	state->v[i][m] -= vcm[m];
    }
  }
}

static void calc_enervirdiff(FILE *log,int eDispCorr,t_forcerec *fr)
{
  double eners[2],virs[2],enersum,virsum,y0,f,g,h;
  double r0,r1,r,rc3,rc9,ea,eb,ec,pa,pb,pc,pd;
  double scale,scale2,scale3;
  real   pot,vir;
  int    ri0,ri1,ri,i,offstart,offset;
  real   *vdwtab; 

  fr->enerdiffsix = 0;
  fr->enerdifftwelve = 0;
  fr->virdiffsix = 0;
  fr->virdifftwelve = 0;

  if (eDispCorr != edispcNO) {
    for(i=0; i<2; i++) {
      eners[i] = 0;
      virs[i]  = 0;
    }
    if (fr->vdwtype == evdwSWITCH || fr->vdwtype == evdwSHIFT
	|| fr->vdwtype == evdwUSER) {
      if (fr->rvdw_switch == 0)
	fatal_error(0,"With dispersion correction rvdw-switch can not be zero "
		    "for vdw-type = %s",evdw_names[fr->vdwtype]);

      /* Round the cut-offs to exact table values for precision */
      ri0 = floor(fr->rvdw_switch*fr->tabscale);
      ri1 = ceil(fr->rvdw*fr->tabscale);
      r0  = ri0/fr->tabscale;
      r1  = ri1/fr->tabscale;

      vdwtab = fr->vdwtab;

      if (fr->vdwtype == evdwSHIFT || fr->vdwtype == evdwUSER) {
	/* Add the constant part from 0 to rvdw_switch */
	rc3  = r0*r0*r0;
	rc9  = rc3*rc3*rc3;
	for(i=0; i<2; i++) {
	  offset = 8*ri0 + (i==0 ? 0 : 4);
          y0 = vdwtab[offset];
          f  = vdwtab[offset+1];
	  if (i == 0) {
	    pot = -1.0/(rc3*rc3);
	    vir =  6.0/(rc3*rc3);
	    if (fr->vdwtype == evdwUSER)
	      fprintf(log,
		      "WARNING: using dispersion correction with user tables\n"
		      "         assuming that the missing dispersion from 0\n"
		      "         to rvdw-switch = %f is constant at %f nm^-6\n",
		      fr->rvdw_switch,pot-y0);
	  } else {
	    pot =   1.0/(rc9*rc3);
	    vir = -12.0/(rc9*rc3);
	  }
	  eners[i] += 4.0*M_PI*(pot - y0)*rc3/3.0;
	  virs[i]  += 4.0*M_PI*(vir - f )*rc3/3.0;
	}
      }

      scale = 1.0/(fr->tabscale);  
      scale2 = scale*scale;
      scale3 = scale*scale2;

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
          r = ri*scale;
          ea = scale3;
          eb = 2.0*scale2*r;
          ec = scale*r*r;
          
          pa = scale3;
          pb = 3.0*scale2*r;
          pc = 3.0*scale*r*r;
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
      rc3  = r0*r0*r0;
      rc9  = rc3*rc3*rc3;
      eners[0] += -4.0*M_PI/(3.0*rc3);
      eners[1] +=  4.0*M_PI/(9.0*rc9);
      virs[0]  +=  8.0*M_PI/rc3;
      virs[1]  += -16.0*M_PI/(3.0*rc9);
    } else if (fr->vdwtype == evdwCUT) {
      rc3  = fr->rvdw*fr->rvdw*fr->rvdw;
      rc9  = rc3*rc3*rc3;
      eners[0] += -4.0*M_PI/(3.0*rc3);
      eners[1] +=  4.0*M_PI/(9.0*rc9);
      virs[0]  +=  8.0*M_PI/rc3;
      virs[1]  += -16.0*M_PI/(3.0*rc9);
    } else {
      fatal_error(0,
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

void calc_dispcorr(FILE *log,int eDispCorr,t_forcerec *fr,int natoms,
		   matrix box,tensor pres,tensor virial,real ener[])
{
  /* modified for switched VdW corrections Michael R. Shirts 2/21/03 */
  static bool bFirst=TRUE;
  real invvol,dens,svir,spres;
  int  m;
  
  ener[F_DISPCORR] = 0.0;
  ener[F_PRES]     = trace(pres)/3.0;
  
  if (eDispCorr != edispcNO) {
    if (bFirst)
      calc_enervirdiff(log,eDispCorr,fr);
    
    invvol = 1/det(box);
    dens = natoms*invvol;
    
    ener[F_DISPCORR] = 0.5*natoms*dens*fr->avcsix*fr->enerdiffsix;
    if (eDispCorr == edispcAllEner || eDispCorr == edispcAllEnerPres)
      ener[F_DISPCORR] += 0.5*natoms*dens*fr->avctwelve*fr->enerdifftwelve;
    
    svir = 0;
    if (eDispCorr == edispcEnerPres || eDispCorr == edispcAllEnerPres) {
      svir += 0.5*natoms*dens*fr->avcsix*fr->virdiffsix/3.0;
      if (eDispCorr == edispcAllEnerPres)
	svir += 0.5*natoms*dens*fr->avctwelve*fr->virdifftwelve/3.0;
    }
    /* The factor 2 is because of the Gromacs virial definition */
    spres = -2.0*invvol*svir*PRESFAC;
    
    for(m=0; m<DIM; m++) {
      virial[m][m] += svir;
      pres[m][m] += spres;
    }
    ener[F_PRES] += spres;

    if (bFirst) {
      if (eDispCorr == edispcEner || eDispCorr == edispcAllEner)
        fprintf(log,"Long Range LJ corr. to Epot: %10g\n",ener[F_DISPCORR]);
      else if (eDispCorr == edispcEnerPres || eDispCorr == edispcAllEnerPres)
        fprintf(log,
		"Long Range LJ corr. to Epot: %10g, Pres: %10g, Vir: %10g\n",
                ener[F_DISPCORR],spres,svir);
      bFirst = FALSE;
    }
  } 
  
  ener[F_EPOT] += ener[F_DISPCORR];
  ener[F_ETOT] += ener[F_DISPCORR];
}


void do_pbc_first(FILE *log,matrix box,rvec box_size,t_forcerec *fr,
		  t_graph *graph,rvec x[])
{
  fprintf(log,"Removing pbc first time\n");
  calc_shifts(box,box_size,fr->shift_vec);
  if (graph) {
    mk_mshift(log,graph,box,x);
    if (getenv ("NOPBC") == NULL)
      shift_self(graph,box,x);
    else
      fprintf(log,"Not doing first shift_self\n");
  }
  fprintf(log,"Done rmpbc\n");
}

void set_pot_bools(t_inputrec *ir,t_topology *top,
		   bool *bLR,bool *bLJLR,bool *bBHAM,bool *b14)
{
  *bLR   = (ir->rcoulomb > ir->rlist) || EEL_LR(ir->coulombtype);
  *bLJLR = (ir->rvdw > ir->rlist);
  *bBHAM = (top->idef.functype[0]==F_BHAM);
  *b14   = (top->idef.il[F_LJ14].nr > 0);
}

void finish_run(FILE *log,t_commrec *cr,char *confout,
		t_nsborder *nsb,t_topology *top,t_parm *parm,
		t_nrnb nrnb[],double nodetime,double realtime,int step,
		bool bWriteStat)
{
  int    i,j;
  t_nrnb ntot;
  real   runtime;
  for(i=0; (i<eNRNB); i++)
    ntot.n[i]=0;
  for(i=0; (i<nsb->nnodes); i++)
    for(j=0; (j<eNRNB); j++)
      ntot.n[j]+=nrnb[i].n[j];
  runtime=0;
  if (bWriteStat) {
    runtime=parm->ir.nsteps*parm->ir.delta_t;
    if (MASTER(cr)) {
      fprintf(stderr,"\n\n");
      print_perf(stderr,nodetime,realtime,runtime,&ntot,nsb->nnodes);
    }
    else
      print_nrnb(log,&(nrnb[nsb->nodeid]));
  }

  if (MASTER(cr)) {
    print_perf(log,nodetime,realtime,runtime,&ntot,nsb->nnodes);
    if (nsb->nnodes > 1)
      pr_load(log,nsb->nnodes,nrnb);
  }
}

void init_md(t_commrec *cr,t_inputrec *ir,tensor box,real *t,real *t0,
	     real *lambda,real *lam0,
	     t_nrnb *mynrnb,bool *bTYZ,t_topology *top,
	     int nfile,t_filenm fnm[],char **traj,
	     char **xtc_traj,int *fp_ene,
	     FILE **fp_dgdl,FILE **fp_field,t_mdebin **mdebin,t_groups *grps,
	     tensor force_vir,tensor pme_vir,
	     tensor shake_vir,t_mdatoms *mdatoms,rvec mu_tot,
	     bool *bNEMD,bool *bSimAnn,t_vcm **vcm,t_nsborder *nsb)
{
  bool bBHAM,b14,bLR,bLJLR;
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

  init_nrnb(mynrnb);
  
  /* Check Environment variables & other booleans */
  *bTYZ=getenv("TYZ") != NULL;
  set_pot_bools(ir,top,&bLR,&bLJLR,&bBHAM,&b14);
  
  if (nfile != -1) {
    /* Filenames */
    *traj     = ftp2fn(efTRN,nfile,fnm);
    *xtc_traj = ftp2fn(efXTC,nfile,fnm);
    
    if (MASTER(cr)) {
      *fp_ene = open_enx(ftp2fn(efENX,nfile,fnm),"w");
      if ((fp_dgdl != NULL) && ir->efep!=efepNO)
	*fp_dgdl =
	  xvgropen(opt2fn("-dgdl",nfile,fnm),
		   "dG/d\\8l\\4","Time (ps)",
		   "dG/d\\8l\\4 (kJ mol\\S-1\\N nm\\S-2\\N \\8l\\4\\S-1\\N)");
      if ((fp_field != NULL) && (ir->ex[XX].n || ir->ex[YY].n ||ir->ex[ZZ].n))
	*fp_field = xvgropen(opt2fn("-field",nfile,fnm),
			     "Applied electric field","Time (ps)",
			     "E (V/nm)");
    } else
      *fp_ene = -1;

    *mdebin = init_mdebin(*fp_ene,grps,&(top->atoms),&(top->idef),
			  bLR,bLJLR,bBHAM,b14,ir->efep!=efepNO,ir->epc,
			  ir->eDispCorr,(TRICLINIC(ir->compress) || TRICLINIC(box)),
			  ir->etc,cr);
  }
  
  /* Initiate variables */  
  clear_mat(force_vir);
  clear_mat(pme_vir);
  clear_mat(shake_vir);
  clear_rvec(mu_tot);
  
  /* Set initial values for invmass etc. */
  update_mdatoms(mdatoms,*lambda,TRUE);

  *vcm = init_vcm(stdlog,top,cr,mdatoms,START(nsb),HOMENR(nsb),ir->nstcomm);
    
  debug_gmx();

  *bNEMD = (ir->opts.ngacc > 1) || (norm(ir->opts.acc[0]) > 0);

  if (ir->eI == eiSD)
    init_sd_consts(ir->opts.ngtc,ir->opts.tau_t,ir->delta_t);

}

