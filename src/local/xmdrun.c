/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
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
  cr_new->pid    = 0;
  cr_new->nprocs = 1;
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
  int     i,start,end;
  real    mu_ave;
  t_atom  *atom;
  t_block *mols;
  
  start = START(nsb);
  end   = start + HOMENR(nsb);  
  
  atom = top->atoms.atom;
  mols = &(top->blocks[ebMOLS]);
  /*
  clear_rvec(mu);
  for(i=start; (i<end); i++)
    for(m=0; (m<DIM); m++)
      mu[m] += q[i]*x[i][m];
  if (PAR(cr)) {
    gmx_sum(DIM,mu,cr);
  }
  */
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
