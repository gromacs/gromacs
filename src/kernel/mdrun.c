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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_mdrun_c = "$Id$";

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
#include "pdebug.h"
#include "mdrun.h"
#include "vcm.h"
#include "ebin.h"
#include "mdebin.h"
#include "dummies.h"
#include "mdrun.h"
#include "physics.h"
#include "pppm.h"
#include "edsam.h"
#include "calcmu.h"

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn)
{
  t_forcerec *fr;
  t_mdebin   *mdebin;
  FILE       *ene,*fmu;
  int        step;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bStopRot,bTYZ,bLR,bBHAM,b14,bRerunMD,bNotLastFrame,bMU;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  char       *traj,*enerfile;
  char       *xtc_traj; /* compressed trajectory filename */
  int        i,m;
  rvec       vcm,box_size,mu_tot;
  rvec       *xx,*vv,*ff;
  int        natoms,status;
  
  /* check if "-rerun" is used. */
  bRerunMD = optRerunMDset(nfile,fnm);
  
  /* check if "-mu" is used */
  bMU = opt2bSet("-mu",nfile,fnm);
  if ((bMU) && MASTER(cr))
    fmu = opt2FILE("-mu",nfile,fnm,"w");
  
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
  if (parm->ir.bSimAnn) {
    SAfactor = 1.0  - t0/parm->ir.zero_temp_time;
    if (SAfactor < 0) 
      SAfactor = 0;
  } else
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
		  graph,cr,&mynrnb,grps,fr,top,edyn);
  where();
    
  /* Compute initial EKin for all.. */
  clear_mat(force_vir);
  clear_mat(shake_vir);
  calc_ke_part(TRUE,START(nsb),HOMENR(nsb),
	       vold,v,vt,&(parm->ir.opts),
	       mdatoms,grps,&mynrnb,
	       lambda,&ener[F_DVDLKIN]);
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm);
  clear_rvec(vcm);
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
  where();
  
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->pid,"Started mdrun");
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    if (bRerunMD) 
      fprintf(stderr,"starting md rerun '%s', reading coordinates from input trajectory '%s'\n\n",*(top->name),opt2fn("-rerun",nfile,fnm));
    else
      fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n\n",
	      *(top->name),parm->ir.nsteps,parm->ir.nsteps*parm->ir.delta_t);
  }
  
  /* Initiate PPPM if necessary */
  if (fr->eeltype == eelPPPM) {
    bool bGGhat = ! fexist(ftp2fn(efHAT,nfile,fnm));
    (void)do_pppm(log,FALSE,bGGhat,ftp2fn(efHAT,nfile,fnm),
		  &parm->ir,top->atoms.nr,x,f,mdatoms->chargeT,box_size,
		  fr->phi,cr,&mynrnb);
  }
  
  /***********************************************************
   *
   *             Loop over MD steps 
   *
   ************************************************************/
  
  /* if rerunMD then read coordinates from input trajectory,
   * set velocities to zero */
  if (bRerunMD) {
    int i;

    natoms=read_first_x(&status,opt2fn("-rerun",nfile,fnm),&t,&x,parm->box);
    for(i=0;(i<natoms);i++) {
      v[i][XX]=0;
      v[i][YY]=0;
      v[i][ZZ]=0;
    }
    bNotLastFrame = (natoms!=0);
  } 

  
  /* loop over MD steps or if rerunMD to end of 
   * input trajectory */
  for(step=0; ((!bRerunMD && (step<parm->ir.nsteps)) || 
	       (bRerunMD && bNotLastFrame)); step++) {
    
    /* for rerun MD always do Neighbour Searching */
    if (bRerunMD) 
      bNS = TRUE;
    else {
      /* Stop Center of Mass motion */
      if (parm->ir.nstcomm == 0) {
	bStopCM=FALSE;
	bStopRot=FALSE;
      } else if (parm->ir.nstcomm > 0) {
	bStopCM=do_per_step(step,parm->ir.nstcomm);
	bStopRot=FALSE;
      } else {
	bStopCM=FALSE;
	bStopRot=do_per_step(step,-parm->ir.nstcomm);
      }
      /* Determine whether or not to do Neighbour Searching */
      bNS=((parm->ir.nstlist && ((step % parm->ir.nstlist)==0)) || (step==0));
    }
    
    /* Construct dummy particles */
    construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
    /* Set values for invmass etc. */
    init_mdatoms(mdatoms,lambda,FALSE);

    /* Calculate total dipole moment if necessary */    
    if (bMU) {
      calc_mu(cr,nsb,x,mdatoms->chargeT,mu_tot);
      if (MASTER(cr))
	write_mu(fmu,mu_tot,parm->box);
    }
    
    /* Calc forces and virial */
    clear_mat(force_vir);
    do_force(log,cr,parm,nsb,force_vir,step,&mynrnb,
	     top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	     lambda,graph,bNS,FALSE,fr);
#ifdef DEBUG
    pr_rvecs(log,0,"force_vir",force_vir,DIM);
#endif     
    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */
    if (!bRerunMD)
      t        = t0   + step*parm->ir.delta_t;
    
    if (parm->ir.bPert)
      lambda   = lam0 + step*parm->ir.delta_lambda;
    if (parm->ir.bSimAnn) {
      SAfactor = 1.0  - t/parm->ir.zero_temp_time;
      if (SAfactor < 0) 
	SAfactor = 0;
    }
     
    /* Spread the force on dummy particle to the other particles... */
    spread_dummy_f(log,x,f,&mynrnb,&top->idef);
    
    if (do_per_step(step,parm->ir.nstxout)) xx=x; else xx=NULL;
    if (do_per_step(step,parm->ir.nstvout)) vv=v; else vv=NULL;
    if (do_per_step(step,parm->ir.nstfout)) ff=f; else ff=NULL;
    
    write_traj(log,cr,traj,
	       nsb,step,t,lambda,nrnb,nsb->natoms,xx,vv,ff,parm->box);
    where();
    
    /* for rerunMD, certain things don't have to be done */
    if (!bRerunMD) {
      if (do_per_step(step,parm->ir.nstxtcout)) {
	write_xtc_traj(log,cr,xtc_traj,nsb,mdatoms,
		       step,t,x,parm->box,parm->ir.xtcprec);
	where();
      }
      
      where();
      
      clear_mat(shake_vir);
      update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
	     &(parm->ir),FALSE,mdatoms,x,graph,
      fr->shift_vec,f,buf,vold,v,vt,parm->pres,parm->box,
      top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn);
#ifdef DEBUG
      pr_rvecs(log,0,"shake_vir",shake_vir,DIM);
#endif
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
    }
    
    if (PAR(cr)) {
      /* Globally (over all CPUs) sum energy, virial etc. */
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm);
      if (!bNS) 
	for(i=0; (i<grps->estat.nn); i++)
	  grps->estat.ee[egLR][i] /= cr->nprocs;
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
    if (!bRerunMD) {
      /* Do center of mass motion removal */
      if (bStopCM) {
	check_cm(log,vcm,mdatoms->tmass);
	do_stopcm(log,HOMENR(nsb),START(nsb),v,vcm,mdatoms->tmass);
	inc_nrnb(&mynrnb,eNR_STOPCM,HOMENR(nsb));
      }
      
      /* Do fit to remove overall rotation */
      if (bStopRot)
	do_stoprot(log,top->atoms.nr,box_size,x,mdatoms->massT);
    
      /* Add force and shake contribution to the virial */
      m_add(force_vir,shake_vir,parm->vir);
    }
    
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    
    if (!bRerunMD) {
      /* Sum the kinetic energies of the groups & calc temp */
      ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
      ener[F_EKIN]=trace(parm->ekin);
      ener[F_ETOT]=ener[F_EPOT]+ener[F_EKIN];
      
      /* Calculate Temperature coupling parameters lambda */
      tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
    
      /* Calculate pressure ! */
      calc_pres(parm->box,parm->ekin,parm->vir,parm->pres);
    }
    
    /* Calculate long range corrections to pressure and energy */
    calc_ljcorr(log,parm->ir.bLJcorr,
		fr,mdatoms->nr,parm->box,parm->pres,ener);
    
    upd_mdebin(mdebin,mdatoms->tmass,step,ener,parm->box,shake_vir,
	       force_vir,parm->vir,parm->pres,grps);
    
    where();
    
    if ( MASTER(cr) && (do_per_step(step,parm->ir.nstprint)) ) {
      print_ebin(ene,log,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    }
    if (bVerbose)
      fflush(log);
    
    /* Time for performance */
    if ((step % 10) == 0)
      update_time();
      
    /* Remaining runtime */
    if (MASTER(cr) && bVerbose && ((step % stepout)==0))
      print_time(stderr,start_t,step,&parm->ir);
    
    /* if rerunMD read next frame from input trajectory */
    if (bRerunMD) 
      bNotLastFrame = read_next_x(status,&t,natoms,x,parm->box);
  }
  /* End of main MD loop */
  
  if (bMU && MASTER(cr))
    fclose(fmu);
    
  if (MASTER(cr)) {
    /* if RerunMD don't print energies of last frame again */
    if ( (parm->ir.nstprint > 1) && !bRerunMD )
      print_ebin(ene,log,step-1,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    
    print_ebin(NULL,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(NULL,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));
  }
  if (ene)
    ffclose(ene);
  
  /* Construct dummy particles, for last output frame */
  construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
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
    "multiple processes on a single processor, which is not very sensible.[PAR]",
    "ED (essential dynamics) sampling is switched on by using the [BB]-ei[bb] flag followed by an .edi",
    "file. The .edi file can be produced using options in the essdyn menu",
    "of the WHAT IF program. mdrun produces a .edo file that contains",
    "projections of positions, velocities and forces onto selected",
    "eigenvectors.[PAR]",
    "With [BB]-rerun[bb] an input trajectory can be given for which ",
    "forces and energies will be (re)calculated.[PAR]",
    "When the [BB]-mu[bb] option is given the total dipole of the simulation",
    "box will be written to a file at every step of the simulation.",
    "This file is in binary and can only be read by [TT]g_dipoles[tt]."
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
    { efTRX, "-rerun", "rerun", ffOPTRD },
    /* function "optRerunMDset" (in runner.c) checks if -rerun is specified */
    { efHAT, "-hat","ghat",   ffOPTRD },
    { efEDI, "-ei", "sam",    ffOPTRD },
    { efEDO, "-eo", "sam",    ffOPTWR },
    { efDAT, "-mu", "dipole", ffOPTWR }
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
  t_edsamyn edyn;
  
  get_pargs(&argc,argv,asize(pa),pa,TRUE);
  cr = init_par(nprocs,argv);
  bVerbose = bVerbose && MASTER(cr);
  edyn.bEdsam=FALSE;
  
  if (MASTER(cr)) {
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,
		      PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_NOGET_PARGS,
		      TRUE,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  }	      
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Berendsen95a");
  }
  
  if (opt2bSet("-ei",NFILE,fnm)) 
    ed_open(NFILE,fnm,&edyn);
    
  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,FALSE,nstepout,&edyn);
  
  exit(0);
  
  return 0;
}

