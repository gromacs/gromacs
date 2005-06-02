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
 * Good gRace! Old Maple Actually Chews Slate
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "sysstuff.h"
#include "string2.h"
#include "led.h"
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
#include "stat.h"
#include "fatal.h"
#include "txtdump.h"
#include "futil.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "filenm.h"
#include "statutil.h"
#include "tgroup.h"
#include "mdrun.h"
#include "vcm.h"
#include "ebin.h"
#include "mdebin.h"
#include "vsite.h"
#include "mdrun.h"
#include "physics.h"
#include "glaasje.h"
    

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph)
{
  t_forcerec *fr;
  t_mdebin   *mdebin;
  FILE       *ene;
  int        step;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bTYZ,bLR,bBHAM,b14;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  char       *traj,*enerfile;
  char       *xtc_traj; /* compressed trajectory filename */
  int        i,m;
  rvec       vcm,box_size;
  rvec       *xx,*vv,*ff;
  
  /* Initial values */
  t0           = parm->ir.init_t;
  if (parm->ir.bPert) {
    lam0         = parm->ir.init_lambda;
    lambda       = lam0;
  }
  else {
    lam0   = 0.0;
    lambda = 0.0;
  } 
  SAfactor     = 1.0;
  
  /* Check Environment variables */
  bTYZ=getenv("TYZ") != NULL;
  
  init_nrnb(&mynrnb);
  
  fr=mk_forcerec();
  init_forcerec(log,fr,&(parm->ir),&(top->blocks[ebMOLS]),cr,
		&(top->blocks[ebCGS]),&(top->idef),mdatoms,parm->box,FALSE);
  for(m=0; (m<DIM); m++)
    box_size[m]=parm->box[m][m];
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  
  fprintf(log,"Removing pbc first time\n");
  mk_mshift(log,graph,parm->box,x);
  shift_self(graph,fr->shift_vec,x);
  fprintf(log,"Done rmpbc\n");

  traj=ftp2fn(efTRJ,nfile,fnm);
  xtc_traj=ftp2fn(efXTC,nfile,fnm);
  where();
  
  if (MASTER(cr)) 
    ene=ftp2FILE(efENE,nfile,fnm,"w");
  else
    ene=NULL;
  where();
  
  bLR=(parm->ir.rlong > parm->ir.rshort);
  bBHAM=(top->idef.functype[0]==F_BHAM);
  b14=(top->idef.il[F_LJ14].nr > 0);
  mdebin=init_mdebin(ene,grps,&(top->atoms),bLR,bBHAM,b14);
  where();
  
  clear_rvec(vcm);
  
  /* Set initial values for invmass etc. */
  init_mdatoms(mdatoms,lambda,TRUE);
  where();
  
  if (parm->ir.bShakeFirst) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,mdatoms,x,vold,buf,f,v,
		  graph,cr,nrnb,grps,fr,top);
  where();
    
  /* Compute initial EKin for all.. */
  calc_ke_part(TRUE,0,top->atoms.nr,
	       vold,v,vt,&(parm->ir.opts),
	       mdatoms,grps,&mynrnb,
	       lambda,&ener[F_DVDLKIN]);
	       
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
  where();
  
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->pid,"Started mdrun");
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n\n",*(top->name),
	    parm->ir.nsteps,parm->ir.nsteps*parm->ir.delta_t);
  }
  /***********************************************************
   *
   *             Loop over MD steps 
   *
   ************************************************************/
  for (step=0; (step<parm->ir.nsteps); step++) {
    /* Stop Center of Mass motion */
    bStopCM=do_per_step(step,parm->ir.nstcomm);
    
    /* Determine whether or not to do Neighbour Searching */
    bNS=((parm->ir.nstlist && ((step % parm->ir.nstlist)==0)) || (step==0));
    
    /* Construct virtual sites */
    construct_vsites(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
    /* Set values for invmass etc. */
    if (parm->ir.bPert)
      init_mdatoms(mdatoms,lambda,FALSE);
    
    /*init_forcerec(log,fr,&(parm->ir),&(top->blocks[ebMOLS]),cr,
      &(top->blocks[ebCGS]),&(top->idef),mdatoms,parm->box,FALSE);
      */	  
    clear_mat(force_vir);
    do_force(log,cr,parm,nsb,force_vir,step,&mynrnb,
	     top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	     lambda,graph,TRUE,bNS,FALSE,fr);
	         
    do_glas(log,START(nsb),HOMENR(nsb),x,f,fr,mdatoms,top->idef.atnr,
	    &parm->ir);
	         
    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. We must output all of this before
     * the update.
     */
    t        = t0   + step*parm->ir.delta_t;
    if (parm->ir.bPert)
      lambda   = lam0 + step*parm->ir.delta_lambda;
    SAfactor = 1.0  - step*parm->ir.delta_t*parm->ir.cooling_rate;
    
    /* Spread the force on virtual sites to the other particles... */
    spread_vsite_f(log,x,f,&mynrnb,&top->idef);
    
    if (do_per_step(step,parm->ir.nstxout)) xx=x; else xx=NULL;
    if (do_per_step(step,parm->ir.nstvout)) vv=v; else vv=NULL;
    if (do_per_step(step,parm->ir.nstfout)) ff=f; else ff=NULL;
    write_traj(log,cr,traj,
	       nsb,step,t,lambda,nrnb,nsb->natoms,xx,vv,ff,parm->box);
    where();

    if (do_per_step(step,parm->ir.nstxtcout)) {
      write_xtc_traj(log,cr,xtc_traj,nsb,
		     step,t,nsb->natoms,x,parm->box,parm->ir.xtcprec);
      where();
    }
    
    where();
    clear_mat(shake_vir);
    update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
	   &(parm->ir),FALSE,mdatoms,x,graph,
	   fr->shift_vec,f,buf,vold,v,vt,parm->pres,parm->box,
	   top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE);
    if (PAR(cr)) 
      accumulate_u(cr,&(parm->ir.opts),grps);
      
    where();
    /* Calculate partial Kinetic Energy (for this processor) 
     * per group!
     */
    calc_ke_part(FALSE,START(nsb),HOMENR(nsb),
		 vold,v,vt,&(parm->ir.opts),
		 mdatoms,grps,&mynrnb,
		 lambda,&ener[F_DVDL]);
    where();
    if (bStopCM)
      calc_vcm(log,HOMENR(nsb),START(nsb),mdatoms->massT,v,vcm);
    
    if (PAR(cr)) {
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm);
      if (!bNS)
	for(i=0; (i<grps->estat.nn); i++)
	  grps->estat.ee[egLR][i] /= cr->nprocs;
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
    if (bStopCM) {
      check_cm(log,vcm,mdatoms->tmass);
      do_stopcm(log,HOMENR(nsb),START(nsb),v,vcm,mdatoms->tmass);
      inc_nrnb(&mynrnb,eNR_STOPCM,HOMENR(nsb));
    }
    
    /* Add force and shake contribution to the virial */
    m_add(force_vir,shake_vir,parm->vir);
    
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    
    /* Sum the kinetic energies of the groups & calc temp */
    ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
    ener[F_EKIN]=trace(parm->ekin);
    ener[F_ETOT]=ener[F_EPOT]+ener[F_EKIN];
    
    /* Calculate Temperature coupling parameters lambda */
    tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
    
    /* Calculate pressure ! */
    calc_pres(parm->box,parm->ekin,parm->vir,parm->pres);
    
    /* Calculate long range corrections to pressure and energy */
    calc_ljcorr(log,parm->ir.userint1,
		fr,mdatoms->nr,parm->box,parm->pres,ener);
    
    upd_mdebin(mdebin,mdatoms->tmass,step,ener,parm->box,shake_vir,
	       force_vir,parm->vir,parm->pres,grps);
	       
    where();
    if ((MASTER(cr) && do_per_step(step,parm->ir.nstprint))) {
      print_ebin(ene,log,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    }
    if (bVerbose)
      fflush(log);
    
    if ((step % 10) == 0)
      update_time();
    if (MASTER(cr) && bVerbose && ((step % stepout)==0))
      print_time(stderr,start_t,step,&parm->ir);
  }
  if (MASTER(cr)) {
    if (parm->ir.nstprint > 1)
      print_ebin(ene,log,step-1,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    
    print_ebin(NULL,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(NULL,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));
  }
  if (ene)
    ffclose(ene);
  
  /* Construct virtual sites, for last output frame */
  construct_vsites(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
  /*free_nslist(log);*/
  
  return start_t;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The mdrun program performs Molecular Dynamics simulations.",
    "It reads the binary topology (.tpb) file and distributes the",
    "topology over processors if needed. The coordinates are passed",
    "around, so that computations can begin.",
    "First a neighbourlist is made, then the forces are computed.",
    "The forces are globally summed, and the velocities and",
    "positions are updated. If necessary shake is performed to constrain",
    "bond lengths and/or bond angles.",
    "Temperature and Pressure can be controlled using weak coupling to a",
    "bath.[PAR]",
    "A number of environment variables can be set to influence the behaviour",
    "of the mdrun program. Most of these are for debugging purposes, but they",
    "sometimes come in handy when porting the software to an",
    "unsupported platform as well. These environment variables", 
    "are listed elsewhere in the manual.[PAR]",
    "The mdrun produces three output file, plus one log file per processor.",
    "The first file is the trajectory, containing coordinates, velocities",
    "etc. The second file contains the coordinates and velocities at the end",
    "of the run plus the computational box. The third file contains energies.",
    "In the log file from processor 0 the energies, temperature, etc. are printed.[PAR]",
    "When run on a parallel computer or with PVM on a cluster of workstations",
    "the [BB]-np[bb] option must be given to indicate the number of",
    "processors. Note that at current PVM does work, but it may launch",
    "multiple processes on a single processor, which is not very sensible."
  };
  char         *lognm=NULL;
  t_commrec    *cr;
  static t_filenm fnm[] = {
    { efTPB, NULL, NULL,      ffREAD },
    { efTRJ, "-o", NULL,      ffWRITE },
    { efXTC, "-x", NULL,      ffOPTWR },
    { efGRO, "-c", "confout", ffWRITE },
    { efENE, "-e", "ener",    ffWRITE },
    { efLOG, "-g", "md",      ffWRITE },
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bVerbose=FALSE,bCompact=TRUE;
  static int  nprocs=1,nDLB=0,nstepout=10;
  static t_pargs pa[] = {
    { "-np",      FALSE, etINT, &nprocs,
      "Number of processors, must be the same as used for grompp. THIS SHOULD BE THE FIRST ARGUMENT ON THE COMMAND LINE FOR MPI" },
    { "-v",       FALSE, etBOOL,&bVerbose, "Verbose mode" },
    { "-compact", FALSE, etBOOL,&bCompact,
      "Write a compact log file, i.e. do not write full virial and energy group matrix (these are also in the energy file, so this is redundant) " },
    { "-dlb",     FALSE, etINT, &nDLB,
      "Use dynamic load balancing every ... step. BUGGY do not use" },
    { "-stepout", FALSE, etINT, &nstepout,
      "Frequency of writing the remaining runtime" }
  };

  get_pargs(argc,argv,asize(pa),pa);
  cr = init_par(nprocs,argv);
  bVerbose = bVerbose && MASTER(cr);
  
  if (MASTER(cr)) {
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,PCA_KEEP_ARGS | PCA_BE_NICE ,
		      NFILE,fnm,0,NULL,asize(desc),desc,0,NULL);
    print_pargs(stderr,asize(pa),pa);
  }	      
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,eCITEGMX);
  }

  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,FALSE,nstepout);
  
  exit(0);
  
  return 0;
}

