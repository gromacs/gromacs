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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_runner_c = "$Id$";

#include <string.h>
#include <time.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "string2.h"
#include "led.h"
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
#include "lutab.h"
#include "mdrun.h"
#include "callf77.h"
#include "pppm.h"

bool optRerunMDset (int nfile, t_filenm fnm[])
{
  return opt2bSet("-rerun",nfile,fnm);
}

void finish_run(FILE *log,t_commrec *cr,
		char *confout,
                t_nsborder *nsb,
		t_topology *top,
		t_parm *parm,
		rvec x[],rvec v[],
		t_nrnb nrnb[],
		double cputime,double realtime,int step,
		bool bWriteStat)
{
  int    i,j;
  t_nrnb ntot;
  
  for(i=0; (i<eNRNB); i++)
    ntot.n[i]=0;
  for(i=0; (i<nsb->nprocs); i++)
    for(j=0; (j<eNRNB); j++)
      ntot.n[j]+=nrnb[i].n[j];
  
  if (bWriteStat) {
    if (MASTER(cr)) {
      fprintf(stderr,"\n\n");
      print_perf(stderr,cputime,realtime,&ntot,nsb->nprocs);
    }
    else
      print_nrnb(log,&(nrnb[nsb->pid]));
  }

  if (MASTER(cr)) {
    print_perf(log,cputime,realtime,&ntot,nsb->nprocs);
    if (nsb->nprocs > 1)
      pr_load(log,nsb->nprocs,nrnb);
    
    fprintf(stderr,"writing final coordinates...\n");
    write_sto_conf(confout,*top->name, &(top->atoms),x,v,parm->box);
  }
}

void mdrunner(t_commrec *cr,int nfile,t_filenm fnm[],bool bVerbose,
	      bool bCompact,int nDlb,bool bNM,int nstepout,t_edsamyn *edyn)
{
  double     cputime,realtime;
  t_parm     *parm;
  rvec       *buf,*f,*vold,*v,*vt,*x,box_size;
  real       *ener;
  t_nrnb     *nrnb;
  t_nsborder *nsb;
  t_topology *top;
  t_groups   *grps;
  t_graph    *graph;
  t_mdatoms  *mdatoms;
  t_forcerec *fr;
  time_t     start_t=0;
  bool       bDummies;
  int        i,m;
    
  /* Initiate everything (snew sets to zero!) */
  snew(ener,F_NRE);
  snew(nsb,1);
  snew(top,1);
  snew(grps,1);
  snew(parm,1);
  snew(nrnb,cr->nprocs);

  /* Initiate invsqrt routines */
#ifdef USEF77
#ifdef FINVSQRT
  fillbuf();
#endif
#endif
#ifdef CINVSQRT
  init_lookup_table(stdlog);
#endif
  
  if (bVerbose && MASTER(cr)) 
    fprintf(stderr,"Getting Loaded...\n");

  if (PAR(cr)) {
    /* The master processor reads from disk, then passes everything
     * around the ring, and finally frees the stuff
     */
    if (MASTER(cr)) 
      distribute_parts(cr->left,cr->right,cr->pid,cr->nprocs,parm,
		       ftp2fn(efTPX,nfile,fnm),nDlb);
    
    /* Every processor (including the master) reads the data from the ring */
    init_parts(stdlog,cr,
	       parm,top,&x,&v,&mdatoms,nsb,
	       MASTER(cr) ? LIST_SCALARS | LIST_PARM : 0);
	       
  }
  else {
    /* Read it up... */
    init_single(stdlog,parm,ftp2fn(efTPX,nfile,fnm),top,&x,&v,&mdatoms,nsb);
  }
  snew(buf,nsb->natoms);
  snew(f,nsb->natoms);
  snew(vt,nsb->natoms);
  snew(vold,nsb->natoms);

  if (bVerbose && MASTER(cr))
    fprintf(stderr,"Loaded with Money\n\n");

  /* Index numbers for parallellism... */
  nsb->pid      = cr->pid;
  top->idef.pid = cr->pid;

  /* Group stuff (energies etc) */
  init_groups(stdlog,mdatoms,&(parm->ir.opts),grps);
  
  /* Periodicity stuff */  
  graph=mk_graph(&(top->idef),top->atoms.nr,FALSE);
#ifdef DEBUG
  p_graph(stdlog,"graph",graph);
#endif
  
  /* Distance Restraints */
  init_disres(stdlog,top->idef.il[F_DISRES].nr,&(parm->ir));

  /* check if there are dummies */
  bDummies=FALSE;
  for(i=0; (i<F_NRE) && !bDummies; i++)
    bDummies = ((interaction_function[i].flags & IF_DUMMY) && 
		(top->idef.il[i].nr > 0));

  /* Initiate forcerecord */
  fr=mk_forcerec();
  init_forcerec(stdlog,fr,&(parm->ir),&(top->blocks[ebMOLS]),cr,
		&(top->blocks[ebCGS]),&(top->idef),mdatoms,parm->box,FALSE);
  /* Initiate box */
  for(m=0; (m<DIM); m++)
    box_size[m]=parm->box[m][m];
    
  /* Initiate PPPM if necessary */
  if (fr->eeltype == eelPPPM)
    init_pppm(stdlog,cr,nsb,FALSE,TRUE,box_size,ftp2fn(efHAT,nfile,fnm),&parm->ir);
		
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
    case eiLD:
      start_t=do_md(stdlog,cr,nfile,fnm,
		    bVerbose,bCompact,bDummies,nstepout,parm,grps,
		    top,ener,x,vold,v,vt,f,buf,
		    mdatoms,nsb,nrnb,graph,edyn,fr,box_size);
      break;
    case eiCG:
      start_t=do_cg(stdlog,nfile,fnm,parm,top,grps,nsb,
		    x,f,buf,mdatoms,parm->ekin,ener,
		    nrnb,bVerbose,cr,graph,fr,box_size);
      break;
    case eiSteep:
      start_t=do_steep(stdlog,nfile,fnm,parm,top,grps,nsb,
		       x,f,buf,mdatoms,parm->ekin,ener,
		       nrnb,bVerbose,bDummies,cr,graph,fr,box_size);
      break;
    default:
      fatal_error(0,"Invalid integrator (%d)...\n",parm->ir.eI);
    }

    /* Some timing stats */  
    if (MASTER(cr)) {
      print_time(stderr,start_t,parm->ir.nsteps,&parm->ir);
      realtime=difftime(time(NULL),start_t);
      if ((cputime=cpu_time()) == 0)
	cputime=realtime;
    }
    else 
      realtime=0;
    
    if (((parm->ir.eI==eiMD) || (parm->ir.eI==eiLD)) &&
	!optRerunMDset(nfile,fnm)) {
      real t      = parm->ir.init_t+parm->ir.nsteps*parm->ir.delta_t;
      real lambda = parm->ir.init_lambda+parm->ir.nsteps*parm->ir.delta_lambda;
      int  step   = parm->ir.nsteps;
      
      write_traj(stdlog,cr,opt2fn("-o",nfile,fnm),
		 nsb,step,t,lambda,nrnb,top->atoms.nr,
		 x,v,f,parm->box);
      if (parm->ir.nstxtcout != 0) { 
	if (do_per_step(step,parm->ir.nstxtcout))
	  write_xtc_traj(stdlog,cr,opt2fn("-x",nfile,fnm),nsb,mdatoms,
			 step,t,x,parm->box,parm->ir.xtcprec);
	if (MASTER(cr))
	  close_xtc_traj();
      }
    }
    
    md2atoms(mdatoms,&(top->atoms),TRUE);
    
    /* Finish up, write some stuff */
    { 
      char *gro=ftp2fn(efSTO,nfile,fnm);
      
      /* if rerunMD, don't write last frame again */
      finish_run(stdlog,cr,gro,nsb,top,parm,
		 x,v,nrnb,cputime,realtime,parm->ir.nsteps,
		 (parm->ir.eI==eiMD) || (parm->ir.eI==eiLD));
    }
  }
  
  /* Does what it says */  
  print_date_and_time(stdlog,cr->pid,"Finished mdrun");

  if (MASTER(cr)) {
    thanx(stderr);
  }
}

