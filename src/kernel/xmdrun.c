/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Glycine aRginine prOline Methionine Alanine Cystine Serine
 */
static char *SRCID_xmdrun_c = "$Id$";
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "sysstuff.h"
#include "string2.h"
#include "nrnb.h"
#include "network.h"
#include "confio.h"
#include "binio.h"
#include "copyrite.h"
#include "smalloc.h"
#include "main.h"
#include "pbc.h"
#include "force.h"
#include "macros.h"
#include "names.h"
#include "mdrun.h"
#include "fatal.h"
#include "txtdump.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "filenm.h"
#include "statutil.h"
#include "tgroup.h"
#include "vcm.h"
#include "trnio.h"
#include "ebin.h"
#include "mdebin.h"
#include "disre.h"
#include "dummies.h"
#include "init_sh.h"
#include "do_gct.h"
#include "physics.h"
#include "sim_util.h"
#include "block_tx.h"
#include "rdgroup.h"
#include "glaasje.h"
#include "edsam.h"
#include "calcmu.h"
#include "ionize.h" 

static t_commrec *cr_msim;

t_commrec *init_msim(t_commrec *cr,int nfile,t_filenm fnm[])
{
  t_commrec *cr_new;
  int  i,ftp;
  char *buf;
  
  cr_msim = cr;
  snew(cr_new,1);
  cr_new->nodeid = 0;
  cr_new->nnodes = 1;
  cr_new->left   = cr->left;
  cr_new->right  = cr->right;
  
  /* Patch file names (except log which has been done already) */
  for(i=0; (i<nfile); i++) {
    /* Because of possible multiple extensions per type we must look 
     * at the actual file name 
     */
    ftp = fn2ftp(fnm[i].fn);
    if (ftp != efLOG) {
#ifdef DEBUGPAR
      fprintf(stderr,"Old file name: %s",fnm[i].fn);
#endif
      buf = par_fn(fnm[i].fn,ftp,cr);
      sfree(fnm[i].fn);
      fnm[i].fn = strdup(buf);
#ifdef DEBUGPAR
      fprintf(stderr,", new: %s\n",fnm[i].fn);
#endif
    }
  }
  
  return cr_new;
}

real mol_dipole(int k0,int k1,atom_id ma[],rvec x[],real q[])
{
  int  k,kk,m;
  rvec mu;
  
  clear_rvec(mu);
  for(k=k0; (k<k1); k++) {
    kk = ma[k];
    for(m=0; (m<DIM); m++) 
      mu[m] += q[kk]*x[kk][m];
  }
  return norm(mu);  /* Dipole moment of this molecule in e nm */
}

real calc_mu_aver(t_commrec *cr,t_nsborder *nsb,rvec x[],real q[],rvec mu,
		  t_topology *top,t_mdatoms *md,int gnx,atom_id grpindex[])
{
  int     i;
  real    mu_ave;
  t_block *mols;
    
  mols = &(top->blocks[ebMOLS]);

  /* I guess we have to parallelise this one! */

  if (gnx > 0) {
    mu_ave = 0.0;
    for(i=0; (i<gnx); i++) {
      int gi = grpindex[i];
      mu_ave += mol_dipole(mols->index[gi],mols->index[gi+1],mols->a,x,q);
    }
    
    return(mu_ave/gnx);
  }
  else
    return 0;
}

#define XMDRUN
#include "../mdlib/md.c"
#include "../kernel/mdrun.c"
