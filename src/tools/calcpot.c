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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_calcpot_c = "$Id$";

#include "vec.h"
#include "calcpot.h"
#include "nrnb.h"
#include "mdebin.h"
#include "mshift.h"
#include "smalloc.h"
#include "force.h"
#include "main.h"
#include "filenm.h"
#include "fatal.h"
#include "mdrun.h"
#include "ns.h"
#include "txtdump.h"

static void c_tabpot(real tabscale,   real VFtab[],
		     int  nri,        int  iinr[],
		     int  shift[],
		     int  jindex[],   int  jjnr[],
		     real pos[],      
		     real facel,      real charge[],
		     real pot[],      real shiftvec[])
{
  /* Local variables */
  const real nul = 0.000000;

  /* Table stuff */
  const real one = 1.000000;
  const real two = 2.000000;
  real r1,r1t,fijC,eps,eps2,Y,F,Fp,Geps,Heps2,VV,FF;
  int  n0,n1,nnn;

  /* General and coulomb stuff */
  int  ii,k,n,jnr,ii3,nj0,nj1,is3,j3,ggid;
  real fxJ,fyJ,fzJ,fxO,fyO,fzO;
  real ixO,iyO,izO,dxO,dyO,dzO;
  real txO,tyO,tzO,vcO,fsO,qO,rsqO,rinv1O,rinv2O;
  real qqO,qj;
  real jx,jy,jz,shX,shY,shZ,poti;

  /* Outer loop (over i particles) starts here */
  for(n=0; (n<nri); n++) {

    /* Unpack shift vector */
    is3               = 3*shift[n];
    shX               = shiftvec[is3];
    shY               = shiftvec[is3+1];
    shZ               = shiftvec[is3+2];

    /* Unpack I particle */
    ii                = iinr[n];
    ii3               = 3*ii;

    /* Charge of i particle(s) divided by 4 pi eps0 */
    qO                = facel*charge[ii];

    /* Bounds for the innerloop */
    nj0               = jindex[n];
    nj1               = jindex[n+1];

    /* Compute shifted i position */
    ixO               = shX + pos[ii3];
    iyO               = shY + pos[ii3+1];
    izO               = shZ + pos[ii3+2];
    poti              = nul;
    
    /* Inner loop (over j-particles) starts right here */
    for(k=nj0; (k<nj1); k++) {

      /* Unpack neighbourlist */
      jnr               = jjnr[k];
      j3                = 3*jnr;
      qj                = facel*charge[jnr];
      jx                = pos[j3];
      jy                = pos[j3+1];
      jz                = pos[j3+2];

      /* First one is for oxygen, with LJ */
      dxO               = ixO - jx;
      dyO               = iyO - jy;
      dzO               = izO - jz;
      rsqO              = dxO*dxO + dyO*dyO + dzO*dzO;

      /* Doing fast invsqrt */
      rinv1O            = invsqrt(rsqO);

      /* O block */
      r1                = one/rinv1O;
      r1t               = r1*tabscale;
      n0                = r1t;
      n1                = 12*n0;
      eps               = r1t-n0;
      eps2              = eps*eps;
      nnn               = n1;
      Y                 = VFtab[nnn];
      F                 = VFtab[nnn+1];
      Geps              = eps*VFtab[nnn+2];
      Heps2             = eps2*VFtab[nnn+3];
      Fp                = F+Geps+Heps2;
      VV                = Y+eps*Fp;

      pot[jnr]         += VV*qO;
      poti             += VV*qj;
      
    }
    pot[ii] += poti;
  }
}

static void low_calc_pot(FILE *log,int nl_type,t_forcerec *fr,
			 rvec x[],t_mdatoms *mdatoms,real pot[])
{
  t_nblist *nlist;
  
  nlist = &fr->nlist_sr[nl_type];
  
  c_tabpot(fr->tabscale,fr->coulvdwtab,nlist->nri,nlist->iinr,
	   nlist->shift,nlist->jindex,nlist->jjnr,
	   x[0],fr->epsfac,mdatoms->chargeA,pot,fr->shift_vec[0]);

  fprintf(log,"There were %d interactions\n",nlist->nrj);
}

void calc_pot(FILE *logf,t_nsborder *nsb,t_commrec *cr,t_groups *grps,
	      t_parm *parm,t_topology *top,rvec x[],t_forcerec *fr,
	      t_mdatoms *mdatoms,real pot[])
{
  static bool        bFirst=TRUE;
  static t_nrnb      nrnb;
  static rvec        *f;
  real        lam=0,dum=0;
  rvec        box_size;
  int         i,m;

  /* Calc the force */
  fprintf(stderr,"Doing single force calculation...\n");

  if (bFirst) {
    snew(f,   nsb->natoms);
    
    bFirst = FALSE;
  }
  /* Reset long range forces if necessary */
  if (fr->bTwinRange) {
    clear_rvecs(nsb->natoms,fr->flr);
    clear_rvecs(SHIFTS,fr->fshift_lr);
  }
  if (parm->ir.epc != epcNO)
      calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  put_charge_groups_in_box(stdlog,0,top->blocks[ebCGS].nr,FALSE,
			   parm->box,box_size,&(top->blocks[ebCGS]),x,
			   fr->cg_cm);
  /* mk_mshift(stdlog,graph,parm->box,x);*/
  /* Do the actual neighbour searching and if twin range electrostatics
   * also do the calculation of long range forces and energies.
   */
  
  ns(logf,fr,x,f,parm->box,grps,&(parm->ir.opts),top,mdatoms,cr,
     &nrnb,nsb,0,lam,&dum);
  for(m=0; (m<DIM); m++)
    box_size[m] = parm->box[m][m];
  for(i=0; (i<mdatoms->nr); i++)
    pot[i] = 0;
  if (debug) {
    pr_rvecs(debug,0,"x",x,mdatoms->nr);
    pr_rvecs(debug,0,"cgcm",fr->cg_cm,top->blocks[ebCGS].nr);
  }
  /* electrostatics from any atom to atoms without LJ */
  low_calc_pot(logf,eNL_QQ,fr,x,mdatoms,pot);
    /* electrostatics from any atom to atoms without charge */
  low_calc_pot(logf,eNL_VDW,fr,x,mdatoms,pot);
  /* electrostatics from any atom to atoms with LJ */
  low_calc_pot(logf,eNL_VDWQQ,fr,x,mdatoms,pot);
}

void init_calcpot(int nfile,t_filenm fnm[],t_topology *top,
		  rvec **x,t_parm *parm,t_commrec *cr,
		  t_graph **graph,t_mdatoms **mdatoms,
		  t_nsborder *nsb,t_groups *grps,
		  t_forcerec **fr,real **pot,
		  matrix box)
{
  real     t,t0,lam,lam0,SAfac;
  bool     bTYZ,bNEMD;
  char     *traj,*xtc_traj;
  rvec     *v,mutot;
  t_nrnb   nrnb;
  t_mdebin *mdebin;
  t_vcm    *vcm=NULL;
  int      fp_ene,m;
  rvec     box_size;
  tensor   force_vir,shake_vir;
  
  /* Initiate */
  cr->nnodes = 1; cr->nodeid    = 0; cr->left   = 0; cr->right  = 1;
  cr->nthreads = 1 ; cr->threadid = 0;
  open_log(ftp2fn(efLOG,nfile,fnm),cr);

  init_nrnb(&nrnb);
  init_single(stdlog,parm,ftp2fn(efTPX,nfile,fnm),top,x,&v,mdatoms,nsb);
  init_md(cr,&(parm->ir),parm->box,&t,&t0,&lam,&lam0,&SAfac,
	  &nrnb,&bTYZ,top,-1,NULL,&traj,&xtc_traj,&fp_ene,NULL,
	  &mdebin,grps,force_vir,shake_vir,*mdatoms,mutot,&bNEMD,&vcm);
  init_groups(stdlog,*mdatoms,&(parm->ir.opts),grps);  

  /* Calculate intramolecular shift vectors to make molecules whole again */
  *graph = mk_graph(&(top->idef),top->atoms.nr,FALSE,FALSE);
  mk_mshift(stdlog,*graph,parm->box,*x);
  
  /* Turn off twin range if appropriate */
  parm->ir.rvdw  = parm->ir.rcoulomb;
  parm->ir.rlist = parm->ir.rcoulomb;
  fprintf(stderr,"Using a coulomb cut-off of %g nm\n",parm->ir.rcoulomb); 
  
  /* Turn off free energy computation */
  parm->ir.efep = 0;

  /* Set vanderwaals to shift, to force tables */
  parm->ir.vdwtype     = evdwSHIFT;
  parm->ir.rvdw_switch = 0.0;
    
  /* Initiate forcerecord */
  *fr = mk_forcerec();
  init_forcerec(stdlog,*fr,&(parm->ir),top,cr,*mdatoms,
		nsb,parm->box,FALSE);

  /* Remove periodicity */  
  for(m=0; (m<DIM); m++)
    box_size[m] = parm->box[m][m];
  if (parm->ir.ePBC != epbcNONE)
    do_pbc_first(stdlog,parm,box_size,*fr,*graph,*x);

  copy_mat(parm->box,box);
      
  snew(*pot,nsb->natoms);
  
  sfree(v);
}
