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
static char *SRCID_stat_c = "$Id$";

#include <string.h>
#include <stdio.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "fatal.h"
#include "network.h"
#include "txtdump.h"
#include "names.h"
#include "physics.h"
#include "vec.h"
#include "maths.h"
#include "mvdata.h"
#include "main.h"
#include "force.h"
#include "nrnb.h"
#include "smalloc.h"
#include "futil.h"
#include "network.h"
#include "rbin.h"
#include "tgroup.h"
#include "xtcio.h"
#include "gmxfio.h"
#include "trnio.h"
#include "statutil.h"

void global_stat(FILE *log,
		 t_commrec *cr,real ener[],
		 tensor fvir,tensor svir,
		 t_grpopts *opts,t_groups *grps,
		 t_nrnb *mynrnb,t_nrnb nrnb[],
		 rvec vcm,rvec mu_tot)
{
  static t_bin *rb=NULL; 
  static int   *itc;
  int    imu,ie,ifv,isv,icm,in[MAXPROC],inn[egNR];
  int    j;
  
  if (rb==NULL) {
    rb=mk_bin();
    snew(itc,opts->ngtc);
  }
  else
    reset_bin(rb);
  
  /* Reset nrnb stuff */
  for(j=0; (j<cr->nprocs); j++)
    init_nrnb(&(nrnb[j]));
  cp_nrnb(&(nrnb[cr->pid]),mynrnb);
  
  /* This routine copies all the data to be summed to one big buffer
   * using the t_bin struct. 
   */
  where();
  ie  = add_binr(log,rb,F_NRE,ener);
  where();
  ifv = add_binr(log,rb,DIM*DIM,fvir[0]);
  where();
  isv = add_binr(log,rb,DIM*DIM,svir[0]);
  where();
  for(j=0; (j<cr->nprocs); j++)
    in[j] = add_bind(log,rb,eNRNB,nrnb[j].n);
  where();
  for(j=0; (j<opts->ngtc); j++) 
    itc[j]=add_binr(log,rb,DIM*DIM,grps->tcstat[j].ekin[0]);
  where();
  for(j=0; (j<egNR); j++)
    inn[j]=add_binr(log,rb,grps->estat.nn,grps->estat.ee[j]);
  where();
  icm = add_binr(log,rb,DIM,vcm);
  where();
  imu = add_binr(log,rb,DIM,mu_tot);
  
  /* Global sum it all */
  sum_bin(rb,cr);
  where();
  
  /* Extract all the data locally */
  extract_binr(rb,ie  ,F_NRE,ener);
  extract_binr(rb,ifv ,DIM*DIM,fvir[0]);
  extract_binr(rb,isv ,DIM*DIM,svir[0]);
  for(j=0; (j<cr->nprocs); j++)
    extract_bind(rb,in[j],eNRNB,nrnb[j].n);
  for(j=0; (j<opts->ngtc); j++) 
    extract_binr(rb,itc[j],DIM*DIM,grps->tcstat[j].ekin[0]);
  for(j=0; (j<egNR); j++)
    extract_binr(rb,inn[j],grps->estat.nn,grps->estat.ee[j]);
  extract_binr(rb,icm,DIM,vcm);
  where();
  extract_binr(rb,imu,DIM,mu_tot);
  where();
  
  /* Small hack for temp only */
  ener[F_TEMP]/=cr->nprocs;
}

int do_per_step(int step,int nstep)
{
  if (nstep != 0) 
    return ((step % nstep)==0); 
  else 
    return 0;
}

int do_any_io(int step, t_inputrec *ir)
{
  return (do_per_step(step,ir->nstxout)  ||
	  do_per_step(step,ir->nstprint) ||
	  do_per_step(step,ir->nstfout)  ||
	  do_per_step(step,ir->nstvout));
}

static void moveit(FILE *log,
		   int left,int right,char *s,rvec xx[],t_nsborder *nsb)
{
  rvec  *temp;
  int   i,m,bP,start,homenr;
    
  if (!xx) 
    return;

  start=nsb->index[nsb->pid];
  homenr=nsb->homenr[nsb->pid];
#ifdef DEBUG
  fprintf(log,"Moving %s for trajectory file, start=%d, homenr=%d\n",
	  s,start,homenr);
#endif
  snew(temp,homenr);
  for(i=0; (i<homenr); i++)
    copy_rvec(xx[start+i],temp[i]);

  move_rvecs(log,FALSE,FALSE,left,right,xx,NULL,nsb->nprocs-1,nsb,NULL);
  
  for(i=0; (i<homenr); i++) {
    bP=0;
    for(m=0; (m<DIM); m++)
      if (xx[start+i][m] != temp[i][m])
	bP=1;
    if (bP)
      fprintf(log,"%s[%5d] before: (%8.3f,%8.3f,%8.3f)"
	      " After: (%8.3f,%8.3f,%8.3f)\n",
	      s,start+i,temp[i][XX],temp[i][YY],temp[i][ZZ],
	      xx[start+i][XX],xx[start+i][YY],xx[start+i][ZZ]);
  }
  sfree(temp);
}

void write_traj(FILE *log,t_commrec *cr,
		char *traj,t_nsborder *nsb,
		int step,real t,real lambda,t_nrnb nrnb[],
		int natoms,rvec *xx,rvec *vv,rvec *ff,matrix box)
{
  static int fp=-1;
  
  if ((fp == -1) && MASTER(cr)) {
#ifdef DEBUG
    fprintf(log,"Going to open trajectory file: %s\n",traj);
#endif
    fp = open_trn(traj,"w");
  }
  
#define MX(xvf) moveit(log,cr->left,cr->right,#xvf,xvf,nsb)
  if (cr->nprocs > 1) {
    MX(xx);
    MX(vv);
    MX(ff);
  }
  if ((xx || vv || ff) && MASTER(cr)) {
    fwrite_trn(fp,step,t,lambda,box,natoms,xx,vv,ff);
    fio_flush(fp);
  }
}

/* XDR stuff for compressed trajectories */
static int xd;

void write_xtc_traj(FILE *log,t_commrec *cr,
		    char *xtc_traj,t_nsborder *nsb,t_mdatoms *md,
		    int step,real t,rvec *xx,matrix box,real prec)
{
  static bool bFirst=TRUE;
  static rvec *x_sel;
  static int  natoms;
  int    i,j;
  
  if ((bFirst) && MASTER(cr)) {
#ifdef DEBUG
    fprintf(log,"Going to open compressed trajectory file: %s\n",xtc_traj);
#endif
    xd = open_xtc(xtc_traj,"w");
    
    /* Count the number of atoms in the selection */
    natoms=0;
    for(i=0; (i<md->nr); i++)
      if (md->cXTC[i] == 0)
	natoms++;
    fprintf(log,"There are %d atoms in your output selection\n",natoms);
    if (natoms != md->nr)
      snew(x_sel,natoms);
    
    bFirst=FALSE;
  }
  
  if (cr->nprocs > 1) {
    MX(xx);
  }
  
  if ((xx) && MASTER(cr)) {
    if (natoms == md->nr)
      x_sel = xx;
    else {
      /* We need to copy everything into a temp array */
      for(i=j=0; (i<md->nr); i++) {
	if (md->cXTC[i] == 0) {
	  copy_rvec(xx[i],x_sel[j]);
	  j++;
	}
      }
    }
    if (write_xtc(xd,natoms,step,t,box,x_sel,prec) == 0) {
      fprintf(stderr,"XTC error. Quitting %s\n",Program());
      exit(1);
    }
  }
}

void close_xtc_traj(void)
{
  close_xtc(xd);
}



