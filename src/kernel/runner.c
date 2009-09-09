/*
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

bool optRerunMDset (int nfile, t_filenm fnm[])
{
  return opt2bSet("-rerun",nfile,fnm);
}

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
		const char *confout,
                t_nsborder *nsb,
		t_topology *top,
		t_parm *parm,
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
  }
}

void mdrunner(t_commrec *cr,int nfile,t_filenm fnm[],bool bVerbose,
	      bool bCompact,int nDlb,bool bNM,int nstepout,t_edsamyn *edyn)
{
  double     cputime=0,realtime;
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
  init_lookup_table(stdlog);
  
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
  } else {
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
		    nrnb,bVerbose,bDummies,cr,graph,fr,box_size);
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
    
    md2atoms(mdatoms,&(top->atoms),TRUE);
    
    /* Finish up, write some stuff */
    { 
      const char *gro=ftp2fn(efSTO,nfile,fnm);
      
      /* if rerunMD, don't write last frame again */
      finish_run(stdlog,cr,gro,nsb,top,parm,
		 nrnb,cputime,realtime,parm->ir.nsteps,
		 (parm->ir.eI==eiMD) || (parm->ir.eI==eiLD));
    }
  }
  
  /* Does what it says */  
  print_date_and_time(stdlog,cr->pid,"Finished mdrun");

  if (MASTER(cr)) {
    thanx(stderr);
  }
}

void init_md(t_commrec *cr,t_inputrec *ir,real *t,real *t0,
	     real *lambda,real *lam0,real *SAfactor,
	     t_nrnb *mynrnb,bool *bTYZ,t_topology *top,
	     int nfile,const t_filenm fnm[],char **traj,char **xtc_traj,
             ener_file_t *fp_ene, t_mdebin **mdebin,t_groups *grps,rvec vcm,
	     tensor force_vir,tensor shake_vir,t_mdatoms *mdatoms)
{
  bool bBHAM,b14,bLR,bLJLR;
  
  /* Initial values */
  *t = *t0       = ir->init_t;
  if (ir->bPert) {
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
  
  /* Filenames */
  *traj     = ftp2fn(efTRN,nfile,fnm);
  *xtc_traj = ftp2fn(efXTC,nfile,fnm);

  if (MASTER(cr)) {
    *fp_ene = open_enx(ftp2fn(efENX,nfile,fnm),"w");
    *mdebin = init_mdebin(*fp_ene,grps,&(top->atoms),&(top->idef),
			  bLR,bLJLR,bBHAM,b14,ir->bPert,ir->epc,
			  ir->bDispCorr);
  }
  else {
    *fp_ene = NULL;
    *mdebin = NULL;
  }
  
  /* Initiate variables */  
  clear_rvec(vcm);
  clear_mat(force_vir);
  clear_mat(shake_vir);
  
  /* Set initial values for invmass etc. */
  init_mdatoms(mdatoms,*lambda,TRUE);
  
  where();
}

void do_pbc_first(FILE *log,t_parm *parm,rvec box_size,t_forcerec *fr,
		  t_graph *graph,rvec x[])
{
  fprintf(log,"Removing pbc first time\n");
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  mk_mshift(log,graph,parm->box,x);
  shift_self(graph,fr->shift_vec,x);
  fprintf(log,"Done rmpbc\n");
}

