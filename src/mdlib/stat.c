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

void global_stat(FILE *fplog,
		 t_commrec *cr,real ener[],
		 tensor fvir,tensor svir,rvec mu_tot,
		 t_inputrec *inputrec,
		 t_groups *grps,bool bSumEkinhOld,
		 gmx_constr_t constr,
		 t_vcm *vcm,int *nabnsb,real *terminate)
{
  static t_bin *rb=NULL; 
  static int   *itc0,*itc1;
  int    ie,ifv,isv,irmsd=0,imu=0,idedl,icm=0,imass=0,ica,inb=0;
  int    ibnsb=-1,iterminate;
  int    icj=-1,ici=-1,icx=-1;
  int    inn[egNR];
  int    j;
  real   *rmsd_data,rbnsb;
  double nb;
  
  if (rb==NULL) {
    rb=mk_bin();
    snew(itc0,inputrec->opts.ngtc);
    snew(itc1,inputrec->opts.ngtc);
  }
  else
    reset_bin(rb);

  /* This routine copies all the data to be summed to one big buffer
   * using the t_bin struct. 
   */
  where();
  ie  = add_binr(rb,F_NRE,ener);
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
  for(j=0; (j<inputrec->opts.ngtc); j++) {
    if (bSumEkinhOld)
      itc0[j]=add_binr(rb,DIM*DIM,grps->tcstat[j].ekinh_old[0]);
    itc1[j]=add_binr(rb,DIM*DIM,grps->tcstat[j].ekinh[0]);
  }
  where();
  idedl = add_binr(rb,1,&(grps->dekindl));
  where();
  for(j=0; (j<egNR); j++)
    inn[j]=add_binr(rb,grps->estat.nn,grps->estat.ee[j]);
  where();
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
  ica   = add_binr(rb,1,&(grps->cosacc.mvcos));
  where();
  if (DOMAINDECOMP(cr)) {
    nb = cr->dd->nbonded_local;
    inb = add_bind(rb,1,&nb);
  }
  where();
  if (nabnsb) {
    rbnsb = *nabnsb;
    ibnsb = add_binr(rb,1,&rbnsb);
  }
  iterminate = add_binr(rb,1,terminate);
  
  /* Global sum it all */
  if (debug)
    fprintf(debug,"Summing %d energies\n",rb->maxreal);
  sum_bin(rb,cr);
  where();
  
  /* Extract all the data locally */
  extract_binr(rb,ie  ,F_NRE,ener);
  extract_binr(rb,ifv ,DIM*DIM,fvir[0]);
  extract_binr(rb,isv ,DIM*DIM,svir[0]);
  if (rmsd_data)
    extract_binr(rb,irmsd,inputrec->eI==eiSD2 ? 3 : 2,rmsd_data);
  if (!NEED_MUTOT(*inputrec))
    extract_binr(rb,imu,DIM,mu_tot);
  for(j=0; (j<inputrec->opts.ngtc); j++) {
    if (bSumEkinhOld)
      extract_binr(rb,itc0[j],DIM*DIM,grps->tcstat[j].ekinh_old[0]);
    extract_binr(rb,itc1[j],DIM*DIM,grps->tcstat[j].ekinh[0]);
  }
  extract_binr(rb,idedl,1,&(grps->dekindl));
  for(j=0; (j<egNR); j++)
    extract_binr(rb,inn[j],grps->estat.nn,grps->estat.ee[j]);
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
  extract_binr(rb,ica,1,&(grps->cosacc.mvcos));
  where();
  if (DOMAINDECOMP(cr)) {
    extract_bind(rb,inb,1,&nb);
    if ((int)(nb + 0.5) != cr->dd->nbonded_global)
      dd_print_missing_interactions(fplog,cr,(int)(nb + 0.5));
  }
  where();
  if (nabnsb) {
    extract_binr(rb,ibnsb,1,&rbnsb);
    *nabnsb = (int)(rbnsb + 0.5);
  }
  where();
  extract_binr(rb,iterminate,1,terminate);
  where();

  /* Small hack for temp only */
  ener[F_TEMP] /= (cr->nnodes - cr->npmenodes);
}

int do_per_step(int step,int nstep)
{
  if (nstep != 0) 
    return ((step % nstep)==0); 
  else 
    return 0;
}

static void moveit(t_commrec *cr,
		   int left,int right,char *s,rvec xx[])
{
  if (!xx) 
    return;

  move_rvecs(cr,FALSE,FALSE,left,right,
	     xx,NULL,(cr->nnodes-cr->npmenodes)-1,NULL);
}

void write_traj(t_commrec *cr,
		int fp_trn,bool bX,bool bV,bool bF,
		int fp_xtc,bool bXTC,int xtc_prec,
		t_topology *top_global,
		int step,real t,
		t_state *state_local,t_state *state_global,
		rvec *f_local,rvec *f_global)
{
  static int  nxtc=-1;
  static rvec *xxtc=NULL;
  t_atoms *atoms;
  t_block *cgs;
  int     i,j;

#define MX(xvf) moveit(cr,GMX_LEFT,GMX_RIGHT,#xvf,xvf)

  if (DOMAINDECOMP(cr)) {
    cgs = &top_global->cgs;
    if (bX || bXTC) dd_collect_vec(cr->dd,cgs,state_local,
				   state_local->x,state_global->x);
    if (bV)         dd_collect_vec(cr->dd,cgs,state_local,
				   state_local->v,state_global->v);
    if (bF)         dd_collect_vec(cr->dd,cgs,state_local,
				   f_local,f_global);
  } else if (cr->nnodes > 1) {
    if (bX || bXTC) MX(state_global->x);
    if (bV)         MX(state_global->v);
    if (bF)         MX(f_global);
  }

  if (MASTER(cr)) {
   atoms = &top_global->atoms;
   if (bX || bV || bF) {
      fwrite_trn(fp_trn,step,t,state_local->lambda,
		 state_local->box,atoms->nr,
		 bX ? state_global->x : NULL,
		 bV ? state_global->v : NULL,
		 bF ? f_global : NULL);
      gmx_fio_flush(fp_trn);
    }
    if (bXTC) {
      if (nxtc == -1) {
	nxtc = 0;
	for(i=0; (i<atoms->nr); i++)
	  if (atoms->atom[i].grpnr[egcXTC] == 0)
	    nxtc++;
	if (nxtc != atoms->nr)
	  snew(xxtc,nxtc);
      }
      if (nxtc == atoms->nr) {
	xxtc = state_global->x;
      } else {
	j = 0;
	for(i=0; (i<atoms->nr); i++)
	  if (atoms->atom[i].grpnr[egcXTC] == 0)
	    copy_rvec(state_global->x[i],xxtc[j++]);
      }
      if (write_xtc(fp_xtc,nxtc,step,t,state_local->box,xxtc,xtc_prec) == 0)
	gmx_fatal(FARGS,"XTC error");
    }
  }
}
