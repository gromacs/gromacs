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
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_runner_c = "$Id$";
#include <string.h>
#include <time.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "string2.h"
#include "network.h"
#include "confio.h"
#include "copyrite.h"
#include "smalloc.h"
#include "main.h"
#include "pbc.h"
#include "force.h"
#include "macros.h"
#include "names.h"
#include "fatal.h"
#include "txtdump.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "statutil.h"
#include "tgroup.h"
#include "nrnb.h"
#include "disre.h"
#include "mdatoms.h"
#include "mdrun.h"
#include "callf77.h"
#include "pppm.h"
#include "pme.h"
#include "xvgr.h"
#include "dummies.h"

void get_cmparm(t_inputrec *ir,int step,bool *bStopCM,bool *bStopRot)
{
  if (ir->nstcomm == 0) {
    *bStopCM  = FALSE;
    *bStopRot = FALSE;
  } else if (ir->nstcomm > 0) {
    *bStopCM  = do_per_step(step,ir->nstcomm);
    *bStopRot = FALSE;
  } else {
    *bStopCM  = FALSE;
    *bStopRot = do_per_step(step,-ir->nstcomm);
  }
}

void set_pot_bools(t_inputrec *ir,t_topology *top,
		   bool *bLR,bool *bLJLR,bool *bBHAM,bool *b14)
{
  *bLR   = (ir->rcoulomb > ir->rlist) || EEL_LR(ir->coulombtype);
  *bLJLR = (ir->rvdw > ir->rlist);
  *bBHAM = (top->idef.functype[0]==F_BHAM);
  *b14   = (top->idef.il[F_LJ14].nr > 0);
}




void finish_run(FILE *log,t_commrec *cr,
		char *confout,
                t_nsborder *nsb,
		t_topology *top,
		t_parm *parm,
		t_nrnb nrnb[],
		double nodetime,double realtime,int step,
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

void mdrunner(t_commrec *cr,int nfile,t_filenm fnm[],bool bVerbose,
	      bool bCompact,int nDlb,bool bNM,int nstepout,t_edsamyn *edyn,
	      unsigned long Flags)
{
  double     nodetime=0,realtime;
  t_parm     *parm;
  rvec       *buf,*f,*vold,*v,*vt,*x,box_size;
  real       tmpr1,tmpr2;
  real       *ener;
  t_nrnb     *nrnb;
  t_nsborder *nsb;
  t_topology *top;
  t_groups   *grps;
  t_graph    *graph;
  t_mdatoms  *mdatoms;
  t_forcerec *fr;
  time_t     start_t=0;
  bool       bDummies,bParDummies;
  t_comm_dummies dummycomm;
  int        i,m;
    
  /* Initiate everything (snew sets to zero!) */
  snew(ener,F_NRE);
  snew(nsb,1);
  snew(top,1);
  snew(grps,1);
  snew(parm,1);
  snew(nrnb,cr->nnodes);
  
  if (bVerbose && MASTER(cr)) 
    fprintf(stderr,"Getting Loaded...\n");

  if (PAR(cr)) {
    /* The master thread on the master node reads from disk, then passes everything
     * around the ring, and finally frees the stuff
     */
    if (MASTER(cr)) 
      distribute_parts(cr->left,cr->right,cr->nodeid,cr->nnodes,parm,
		       ftp2fn(efTPX,nfile,fnm),nDlb);
    
    /* Every node (including the master) reads the data from the ring */
    init_parts(stdlog,cr,
	       parm,top,&x,&v,&mdatoms,nsb,
	       MASTER(cr) ? LIST_SCALARS | LIST_PARM : 0, 
	       &bParDummies,&dummycomm);
  } else {
    /* Read it up... */
    init_single(stdlog,parm,ftp2fn(efTPX,nfile,fnm),top,&x,&v,&mdatoms,nsb);
    bParDummies=FALSE;
  }
  snew(buf,nsb->natoms);
  snew(f,nsb->natoms);
  snew(vt,nsb->natoms);
  snew(vold,nsb->natoms);

  if (bVerbose && MASTER(cr))
    fprintf(stderr,"Loaded with Money\n\n");
  
  /* Index numbers for parallellism... */
  nsb->nodeid      = cr->nodeid;
  top->idef.nodeid = cr->nodeid;

  /* Group stuff (energies etc) */
  init_groups(stdlog,mdatoms,&(parm->ir.opts),grps);
  /* Copy the cos acceleration to the groups struct */
  grps->cosacc.cos_accel = parm->ir.cos_accel;
  
  /* Periodicity stuff */  
  graph=mk_graph(&(top->idef),top->atoms.nr,FALSE,FALSE);
  if (debug)
    p_graph(debug,"Initial graph",graph);
  
  /* Distance Restraints */
  init_disres(stdlog,top->idef.il[F_DISRES].nr,&(parm->ir));

  /* check if there are dummies */
  bDummies=FALSE;
  for(i=0; (i<F_NRE) && !bDummies; i++)
    bDummies = ((interaction_function[i].flags & IF_DUMMY) && 
		(top->idef.il[i].nr > 0));

  /* Initiate forcerecord */
  fr = mk_forcerec();
  init_forcerec(stdlog,fr,&(parm->ir),top,cr,mdatoms,nsb,parm->box,FALSE,
		opt2fn("-table",nfile,fnm),FALSE);
  fr->bSepDVDL = ((Flags & MD_SEPDVDL) == MD_SEPDVDL);
  /* Initiate box */
  for(m=0; (m<DIM); m++)
    box_size[m]=parm->box[m][m];
    
  /* Initiate PPPM if necessary */
  if (fr->eeltype == eelPPPM)
    init_pppm(stdlog,cr,nsb,FALSE,TRUE,box_size,getenv("GMXGHAT"),&parm->ir);
  if (fr->eeltype == eelPME)
    init_pme(stdlog,cr,parm->ir.nkx,parm->ir.nky,parm->ir.nkz,parm->ir.pme_order,
	     HOMENR(nsb),parm->ir.bOptFFT);
  
  /* Now do whatever the user wants us to do (how flexible...) */
  if (bNM) {
    start_t=do_nm(stdlog,cr,nfile,fnm,
		  bVerbose,bCompact,nstepout,parm,grps,
		  top,ener,x,vold,v,vt,f,buf,
		  mdatoms,nsb,nrnb,graph,edyn,fr,box_size);
  }
  else {
    switch (parm->ir.eI) {
    case eiMD:
    case eiSD:
    case eiBD:
      start_t=do_md(stdlog,cr,nfile,fnm,
		    bVerbose,bCompact,bDummies,
		    bParDummies ? &dummycomm : NULL,
		    nstepout,parm,grps,top,ener,x,vold,v,vt,f,buf,
		    mdatoms,nsb,nrnb,graph,edyn,fr,box_size,Flags);
      break;
    case eiCG:
      start_t=do_cg(stdlog,nfile,fnm,parm,top,grps,nsb,
		    x,f,buf,mdatoms,parm->ekin,ener,
		    nrnb,bVerbose,bDummies,
		    bParDummies ? &dummycomm : NULL,
		    cr,graph,fr,box_size);
      break;
    case eiSteep:
      start_t=do_steep(stdlog,nfile,fnm,parm,top,grps,nsb,
		       x,f,buf,mdatoms,parm->ekin,ener,
		       nrnb,bVerbose,bDummies,
		       bParDummies ? &dummycomm : NULL,
		       cr,graph,fr,box_size);
      break;
    default:
      fatal_error(0,"Invalid integrator (%d)...\n",parm->ir.eI);
    }

    /* Some timing stats */  
    if (MASTER(cr)) {
      realtime=difftime(time(NULL),start_t);
      if ((nodetime=node_time()) == 0)
	nodetime=realtime;
    }
    else 
      realtime=0;
    md2atoms(mdatoms,&(top->atoms),TRUE);
    /* Finish up, write some stuff */
    { 
      char *gro=ftp2fn(efSTO,nfile,fnm);
      /* if rerunMD, don't write last frame again */
      finish_run(stdlog,cr,gro,nsb,top,parm,
		 nrnb,nodetime,realtime,parm->ir.nsteps,
		 parm->ir.eI==eiMD || parm->ir.eI==eiSD || parm->ir.eI==eiBD);
    }
  }
  /* Does what it says */  
  print_date_and_time(stdlog,cr->nodeid,"Finished mdrun");

  if (MASTER(cr)) {
    thanx(stderr);
  }
}


void init_md(t_commrec *cr,t_inputrec *ir,tensor box,real *t,real *t0,
	     real *lambda,real *lam0,real *SAfactor,
	     t_nrnb *mynrnb,bool *bTYZ,t_topology *top,
	     int nfile,t_filenm fnm[],char **traj,char **xtc_traj,int *fp_ene,
	     FILE **fp_dgdl,t_mdebin **mdebin,t_groups *grps,
	     tensor force_vir,tensor pme_vir,
	     tensor shake_vir,t_mdatoms *mdatoms,rvec mu_tot,
	     bool *bNEMD,t_vcm **vcm,t_nsborder *nsb)
{
  bool bBHAM,b14,bLR,bLJLR;
  int  i;
  
  /* Initial values */
  *t = *t0       = ir->init_t;
  if (ir->efep != efepNO) {
    *lambda = *lam0 = ir->init_lambda;
  }
  else {
    *lambda = *lam0   = 0.0;
  } 
  if (ir->bSimAnn) {
    *SAfactor = 1.0 - *t0/ir->zero_temp_time;
    if (*SAfactor < 0) 
      *SAfactor = 0;
  } else
    *SAfactor     = 1.0;
    
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
    } else
      *fp_ene = -1;

    *mdebin = init_mdebin(*fp_ene,grps,&(top->atoms),&(top->idef),
			  bLR,bLJLR,bBHAM,b14,ir->efep!=efepNO,ir->epc,
			  ir->eDispCorr,(TRICLINIC(ir->compress) || TRICLINIC(box)),
			  (ir->etc==etcNOSEHOOVER),cr);
  }
  
  /* Initiate variables */  
  clear_mat(force_vir);
  clear_mat(pme_vir);
  clear_mat(shake_vir);
  clear_rvec(mu_tot);
  
  /* Set initial values for invmass etc. */
  init_mdatoms(mdatoms,*lambda,TRUE);

  *vcm = init_vcm(stdlog,top,mdatoms,
		  START(nsb),HOMENR(nsb),ir->nstcomm);
    
  debug_gmx();

  *bNEMD = (ir->opts.ngacc > 1) || (norm(ir->opts.acc[0]) > 0);

  if (ir->eI == eiSD)
    init_sd_consts(ir->opts.ngtc,ir->opts.tau_t,ir->delta_t);

}

void do_pbc_first(FILE *log,t_parm *parm,rvec box_size,t_forcerec *fr,
		  t_graph *graph,rvec x[])
{
  fprintf(log,"Removing pbc first time\n");
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  mk_mshift(log,graph,parm->box,x);
  if (getenv ("NOPBC") == NULL)
    shift_self(graph,parm->box,x);
  else
    fprintf(log,"Not doing first shift_self\n");
  fprintf(log,"Done rmpbc\n");
}

