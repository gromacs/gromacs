/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_sim_util_c = "$Id$";

#include <stdio.h>
#include <time.h>
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

#define difftime(end,start) ((double)(end)-(double)(start))

void print_time(FILE *out,time_t start,int step,t_inputrec *ir)
{
  static real time_per_step;
  static time_t end;
  time_t finish;
  double dt;
  char buf[48];

  fprintf(out,"\rstep %d",step);
  if ((step >= ir->nstlist)) {
    if ((ir->nstlist == 0) || ((step % ir->nstlist) == 0)) {
      /* We have done a full cycle let's update time_per_step */
      end=time(NULL);
      dt=difftime(end,start);
      time_per_step=dt/step;
    }
    dt=(ir->nsteps-step)*time_per_step;

    if (dt >= 300) {    
      finish = end+(time_t)dt;
      sprintf(buf,"%s",ctime(&finish));
      buf[strlen(buf)-1]='\0';
      fprintf(out,", will finish at %s",buf);
    }
    else
      fprintf(out,", remaining runtime: %5d s    ",(int)dt);
  }
  fflush(out);
}

time_t print_date_and_time(FILE *log,int pid,char *title)
{
  int i;
  char *ts,time_string[STRLEN];
  time_t now;

  now=time(NULL);
  ts=ctime(&now);
  for (i=0; ts[i]>=' '; i++) time_string[i]=ts[i];
  time_string[i]='\0';
  fprintf(log,"%s on processor %d %s\n",title,pid,time_string);
  return now;
}

static void pr_commrec(FILE *log,t_commrec *cr)
{
  fprintf(log,"commrec: pid=%d, nprocs=%d, left=%d, right=%d\n",
	  cr->pid,cr->nprocs,cr->left,cr->right);
}

static void sum_forces(int start,int end,rvec f[],rvec flr[])
{
  int i;
  
  for(i=start; (i<end); i++)
    rvec_inc(f[i],flr[i]);
}

static void reset_forces(bool bNS,rvec f[],t_forcerec *fr,int natoms)
{
  int i;
  
  if (fr->bTwinRange) {
    if (bNS) {
      clear_rvecs(natoms,fr->flr);
      clear_rvecs(SHIFTS,fr->fshift_lr);
    }
    else {
      for(i=0; (i<natoms); i++)
	copy_rvec(fr->flr[i],f[i]);
      for(i=0; (i<SHIFTS); i++)
	copy_rvec(fr->fshift_lr[i],fr->fshift[i]);
    } 
  }
  else {
    if (fr->eeltype == eelPPPM) 
      clear_rvecs(natoms,fr->flr);
    clear_rvecs(natoms,f);
    clear_rvecs(SHIFTS,fr->fshift);
  }
}

static void reset_energies(t_grpopts *opts,t_groups *grp,
			   t_forcerec *fr,bool bNS,real epot[])
{
  const real zero=0.0;
  int   i,j;
  
  /* First reset all energy components but the Long Range */
  for(i=0; (i<egNR); i++)
    if (i != egLR)
      for(j=0; (j<grp->estat.nn); j++)
	grp->estat.ee[i][j]=0.0;

  /* If method of long range is twin range and we do neighboursearching,
   * or if the method is PPPM
   */
  if ((fr->bTwinRange && bNS) || (fr->eeltype == eelPPPM)) {
    for(i=0; (i<grp->estat.nn); i++) 
      grp->estat.ee[egLR][i]=0.0;
  }
	
  /* Normal potential energy components */
  for(i=0; (i<=F_EPOT); i++)
    epot[i] = zero;
  epot[F_DVDL]    = zero;
  epot[F_DVDLKIN] = zero;

}

void do_force(FILE *log,t_commrec *cr,
	      t_parm *parm,t_nsborder *nsb,tensor vir_part,
	      int step,t_nrnb *nrnb,t_topology *top,t_groups *grps,
	      rvec x[],rvec v[],rvec f[],rvec buf[],
	      t_mdatoms *mdatoms,real ener[],bool bVerbose,
	      real lambda,t_graph *graph,
	      bool bNS,bool bMolEpot,t_forcerec *fr)
{
  static rvec box_size;
  int    pid,cg0,cg1;
  int    start,homenr;
  
  pid    = cr->pid;
  start  = START(nsb);
  homenr = HOMENR(nsb);
  cg0    = (pid == 0) ? 0 : nsb->cgload[pid-1];
  cg1    = nsb->cgload[pid];
  
  where();
  update_forcerec(log,fr,parm->box);
  where();
  
  /* Compute shift vectors every step, because of pressure coupling! */
  if (parm->ir.epc != epcNO)
    calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  where();
  
  if (bNS) {
    put_atoms_in_box(log,cg0,cg1,FALSE,
		     parm->box,box_size,&(top->blocks[ebCGS]),x,
		     fr->shift_vec,fr->cg_cm);
    inc_nrnb(nrnb,eNR_RESETX,homenr);
    inc_nrnb(nrnb,eNR_CGCM,cg1-cg0);

    where();    
    if (PAR(cr))
      move_cgcm(log,cr,fr->cg_cm,nsb->cgload);
#ifdef DEBUG
    pr_rvecs(log,0,"cgcm",fr->cg_cm,nsb->cgtotal);
#endif
  }
  where();
  if (PAR(cr)) 
    move_x(log,cr->left,cr->right,x,nsb,nrnb);
  where();

  /* Reset the forces and the enrgies */
  reset_forces(bNS,f,fr,nsb->natoms);
  reset_energies(&(parm->ir.opts),grps,fr,bNS,ener);    
  where();
  
  if (bNS) {
    /* Calculate intramolecular shift vectors to make molecules whole again */
    mk_mshift(log,graph,parm->box,x);

    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    ns(log,fr,x,f,parm->box,grps,&(parm->ir.opts),top,mdatoms,
       cr,nrnb,nsb,step);
  }
  
  /* Compute the forces */    
  force(log,step,fr,&(parm->ir),&(top->idef),nsb,cr,nrnb,grps,mdatoms,
	top->atoms.grps[egcENER].nr,&(parm->ir.opts),
	x,f,vir_part,ener,bVerbose,parm->box,lambda,graph,&(top->atoms.excl));
  where();
#ifdef DEBUG
  if (bNS)
    print_nrnb(log,nrnb);
#endif

  if (fr->eeltype != eelPPPM) {
    /* This summing need actually not be done for all coordinates... */
    if (fr->bTwinRange)
      sum_forces(0,nsb->natoms,f,fr->flr);
    if (PAR(cr)) 
      move_f(log,cr->left,cr->right,f,buf,nsb,nrnb);
  
    /* Calculate virial */
    f_calc_vir(log,start,start+homenr,x,f,vir_part,cr,graph,fr->shift_vec);
    inc_nrnb(nrnb,eNR_VIRIAL,homenr);
  }
  else {
    /* If we do PPPM the long range forces should not be taken into account
     * for computation of the virial. Rather a special formula
     * due to Neumann is used.
     */
    f_calc_vir(log,0,nsb->natoms,x,f,vir_part,cr,graph,fr->shift_vec);
    inc_nrnb(nrnb,eNR_VIRIAL,nsb->natoms);
    sum_forces(0,nsb->natoms,f,fr->flr);
    if (PAR(cr)) 
      move_f(log,cr->left,cr->right,f,buf,nsb,nrnb);
  }
  where(); 
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

double cpu_time(void)
{
  return runtime;
}

void do_shakefirst(FILE *log,bool bTYZ,real lambda,real ener[],
		   t_parm *parm,t_nsborder *nsb,t_mdatoms *md,
		   rvec x[],rvec vold[],rvec buf[],rvec f[],
		   rvec v[],t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
		   t_groups *grps,t_forcerec *fr,t_topology *top,
		   t_edsamyn *edyn)
{
  int    i,m;
  tensor shake_vir;
  real   dt=parm->ir.delta_t;
  real   dt_1;

  if ((top->idef.il[F_SHAKE].nr > 0)   ||
      (top->idef.il[F_SETTLE].nr > 0)) {
    /* Do a first SHAKE to reset particles... */
    clear_mat(shake_vir);
    update(nsb->natoms,START(nsb),HOMENR(nsb),
	   -1,lambda,&ener[F_DVDL],
	   &(parm->ir),FALSE,md,x,graph,
	   fr->shift_vec,NULL,NULL,vold,x,NULL,parm->pres,parm->box,
	   top,grps,shake_vir,cr,nrnb,bTYZ,FALSE,edyn);
    
    /* Compute coordinates at t=-dt, store them in buf */
    for(i=0; (i<nsb->natoms); i++) {
      for(m=0; (m<DIM); m++) {
	f[i][m]=x[i][m];
	buf[i][m]=x[i][m]-dt*v[i][m];
      }
    }
    
    /* Shake the positions at t=-dt with the positions at t=0
     * as reference coordinates.
     */
    clear_mat(shake_vir);
    update(nsb->natoms,START(nsb),HOMENR(nsb),
	   0,lambda,&ener[F_DVDL],&(parm->ir),FALSE,md,f,graph,
	   fr->shift_vec,NULL,NULL,vold,buf,NULL,parm->pres,parm->box,
	   top,grps,shake_vir,cr,nrnb,bTYZ,FALSE,edyn);
    
    /* Compute the velocities at t=-dt/2 using the coordinates at
     * t=-dt and t=0
     */
    dt_1=1.0/dt;
    for(i=0; (i<nsb->natoms); i++) {
      for(m=0; (m<DIM); m++)
	v[i][m]=(x[i][m]-f[i][m])*dt_1;
    }
  
    /* Shake the positions at t=-dt with the positions at t=0
     * as reference coordinates.
     */
    clear_mat(shake_vir);
    update(nsb->natoms,START(nsb),HOMENR(nsb),
	   0,lambda,&ener[F_DVDL],&(parm->ir),FALSE,
	   md,f,graph,
	   fr->shift_vec,NULL,NULL,vold,buf,NULL,parm->pres,parm->box,
	   top,grps,shake_vir,cr,nrnb,bTYZ,FALSE,edyn);
    
    /* Compute the velocities at t=-dt/2 using the coordinates at
     * t=-dt and t=0
     */
    dt_1=1.0/dt;
    for(i=0; (i<nsb->natoms); i++) {
      for(m=0; (m<DIM); m++)
	v[i][m]=(x[i][m]-f[i][m])*dt_1;
    }
  }
}

void calc_ljcorr(FILE *log,bool bLJcorr,t_forcerec *fr,int natoms,
		 matrix box,tensor pres,tensor virial,real ener[])
{
  static bool bFirst=TRUE;
  real vol,rc3,spres,svir;
  int  m;
  
  if (bLJcorr) {
    vol           = det(box);
    rc3           = fr->rshort*fr->rshort*fr->rshort;
    ener[F_LJLR]  = -2.0*natoms*natoms*M_PI*fr->avcsix/(3.0*vol*rc3);
    spres         = 2.0*ener[F_LJLR]*PRESFAC/vol;
    svir          = -6.0*ener[F_LJLR];
    ener[F_PRES]  = trace(pres)/3.0+spres;
    for(m=0; (m<DIM); m++) {
      pres[m][m]    += spres;
      virial[m][m]  += svir;
    }
    if (bFirst) {
      fprintf(log,"Long Range LJ corrections: Epot=%10g, Pres=%10g, Vir=%10g\n",
	      ener[F_LJLR],spres,svir);
      bFirst = FALSE;
    }
  }
  else {
    ener[F_LJLR]  = 0.0;
    ener[F_PRES]  = trace(pres)/3.0;
  }
  ener[F_EPOT] += ener[F_LJLR];
  ener[F_ETOT] += ener[F_LJLR];
}

