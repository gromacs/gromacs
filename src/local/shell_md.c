/*
 *       %Z% %M% %I% %G%
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 * Copyright (c) 1990-1995,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * S  C  A  M  O  R  G
 */

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
#include "ebin.h"
#include "mdebin.h"
#include "disre.h"
#include "dummies.h"
#include "init_sh.h"
#include "do_gct.h"
#include "physics.h"
#include "block_tx.h"
#include "rdgroup.h"
#include "edsam.h"
#include "calcmu.h"

static bool      bMultiSim = FALSE;
static t_commrec *cr_msim;

static void do_1pos(rvec xnew,rvec xold,rvec f,real k_1,real step)
{
  real xo,yo,zo;
  real dx,dy,dz,dx2;
  
  xo=xold[XX];
  yo=xold[YY];
  zo=xold[ZZ];
  
  dx=f[XX]*k_1;
  dy=f[YY]*k_1;
  dz=f[ZZ]*k_1;
  
  xnew[XX]=xo+dx*step;
  xnew[YY]=yo+dy*step;
  xnew[ZZ]=zo+dz*step;
}

static void shell_pos_sd(FILE *log,real step,rvec xold[],rvec xnew[],rvec f[],
			 int ns,t_shell s[])
{
  int  i,shell;
  real k_1;
  
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    k_1   = s[i].k_1;
    do_1pos(xnew[shell],xold[shell],f[shell],k_1,step);
  }
}

static void predict_shells(FILE *log,rvec x[],rvec v[],real dt,int ns,t_shell s[])
{
  int  i,m;
  real hdt,tdt;
  
  hdt = 0.5*dt;
  tdt = dt/3.0;
  
  for(i=0; (i<ns); i++) {
    switch (s[i].nnucl) {
    case 1:
      for(m=0; (m<DIM); m++)
	x[s[i].shell][m]+=v[s[i].nucl1][m]*dt;
      break;
    case 2:
      for(m=0; (m<DIM); m++)
	x[s[i].shell][m]+=(v[s[i].nucl1][m]+v[s[i].nucl2][m])*hdt;
      break;
    case 3:
      for(m=0; (m<DIM); m++)
	x[s[i].shell][m]+=(v[s[i].nucl1][m]+v[s[i].nucl2][m]+v[s[i].nucl3][m])*tdt;
      break;
    default:
      fatal_error(0,"Shell %d has %d nuclei!",i,s[i].nnucl);
    }
  }
}

static void print_epot(char *which,
		       int mdstep,int count,real step,real epot,real df)
{
  fprintf(stderr,"MDStep=%5d/%2d lamb: %6g, E-Pot: %12.8e",
	  mdstep,count,step,epot);
  
  if (count != 0)
    fprintf(stderr,", rmsF: %12.8e %s\n",df,which);
  else
    fprintf(stderr,"\n");
}


static real rms_force(rvec f[],int ns,t_shell s[])
{
  int  i,shell;
  real df2;
  
  df2=0.0;
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    df2  += iprod(f[shell],f[shell]);
  }
  return sqrt(df2/ns);
}

static void dump_shells(FILE *fp,rvec x[],rvec f[],real ftol,int ns,t_shell s[])
{
  int  i,shell;
  real ft2,ff2;
  
  ft2 = sqr(ftol);
  
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    ff2   = iprod(f[shell],f[shell]);
    if (ff2 > ft2)
      fprintf(fp,"SHELL %5d, force %10.5f  %10.5f  %10.5f\n",
	      shell,f[shell][XX],f[shell][YY],f[shell][ZZ]);
    if (ff2 > 100*ft2)
      pr_rvecs(fp,0,"SHELL-X",x+(5*(shell / 5)),5);
  }
}

static int relax_shells(FILE *ene,FILE *log,t_commrec *cr,bool bVerbose,
			int mdstep,t_parm *parm,bool bDoNS,bool bStopCM,
			t_topology *top,real ener[],
			rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
			rvec buf[],t_mdatoms *md,t_nsborder *nsb,t_nrnb *nrnb,
			t_graph *graph,t_groups *grps,tensor vir_part,
			int nshell,t_shell shells[],t_forcerec *fr,
			char *traj,real t,real lambda,
			int natoms,matrix box,t_mdebin *mdebin)
{
  static bool bFirst=TRUE;
  static rvec *pos[2],*force[2];
  real   Epot[2],df[2];
#define NEPOT asize(Epot)
  real   ftol,step;
  real   step0;
  bool   bDone,bMinSet;
  int    g;
  int    number_steps;
  int    count=0;
  int    i,start=START(nsb),homenr=HOMENR(nsb),end=START(nsb)+HOMENR(nsb);
  int    Min=0;
#define  Try (1-Min)             /* At start Try = 1 */

  if (bFirst) {
    /* Allocate local arrays */
    for(i=0; (i<2); i++) {
      snew(pos[i],nsb->natoms);
      snew(force[i],nsb->natoms);
    }
    bFirst=FALSE;
  }
  
  ftol         = parm->ir.em_tol;
  number_steps = parm->ir.niter;
  step0        = 1.0;

  /* Do a prediction of the shell positions */
  predict_shells(log,x,v,parm->ir.delta_t,nshell,shells);
    
  /* Calculate the forces first time around */
  do_force(log,cr,parm,nsb,vir_part,mdstep,nrnb,
	   top,grps,x,v,force[Min],buf,md,ener,bVerbose && !PAR(cr),
	   lambda,graph,bDoNS,FALSE,fr);
  df[Min]=rms_force(force[Min],nshell,shells);
  
  /* Copy x to pos[Min] & pos[Try]: during minimization only the
   * shell positions are updated, therefore the other particles must
   * be set here.
   */
  memcpy(pos[Min],x,nsb->natoms*sizeof(x[0]));
  memcpy(pos[Try],x,nsb->natoms*sizeof(x[0]));

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  Epot[Min]=ener[F_EPOT];
  if (PAR(cr))
    gmx_sum(NEPOT,Epot,cr);
  
  step=step0;
  if (bVerbose && MASTER(cr))
    print_epot("",mdstep,0,step,Epot[Min],df[Min]);

  if (debug) {
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EKIN].longname, ener[F_EKIN]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EPOT].longname, ener[F_EPOT]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_ETOT].longname, ener[F_ETOT]);
  }

  if (debug)
    fprintf(debug,"SHELLSTEP %d\n",mdstep);
    
  bDone=((nshell == 0) || (df[Min] < ftol));
  bMinSet=FALSE;
  for(count=1; 
      !(bDone || ((number_steps > 0) && (count>=number_steps))); count++) {
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try],force[Min],nshell,shells); 

#ifdef DEBUG
    if (debug) {
      pr_rvecs(debug,0,"pos[Try] b4 do_force",&(pos[Try][start]),homenr);
      pr_rvecs(debug,0,"pos[Min] b4 do_force",&(pos[Min][start]),homenr);
    }
#endif
    /* Try the new positions */
    do_force(log,cr,parm,nsb,vir_part,1,nrnb,
	     top,grps,pos[Try],v,force[Try],buf,md,ener,bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr);
    df[Try]=rms_force(force[Try],nshell,shells);
#ifdef DEBUG
    if (debug) 
      pr_rvecs(debug,0,"F na do_force",&(force[Try][start]),homenr);
#endif
    if (debug)
      dump_shells(debug,pos[Try],force[Try],ftol,nshell,shells);
      
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    Epot[Try]=ener[F_EPOT];
    if (PAR(cr)) 
      gmx_sum(NEPOT,Epot,cr);
    
    if (bVerbose && MASTER(cr))
      print_epot("",mdstep,count,step,Epot[Try],df[Try]);

    bDone=(df[Try] < ftol);
    
    if ((Epot[Try] < Epot[Min])) {
      Min  = Try;
      step = step0;
    }
    else
      step *= 0.2;
  }
  if (MASTER(cr) && !bDone) 
    fprintf(stderr,"EM did not converge in %d steps\n",number_steps);
  
  /* Parallelise this one! */
  memcpy(x,pos[Min],nsb->natoms*sizeof(x[0]));
  memcpy(f,force[Min],nsb->natoms*sizeof(f[0]));
  /*
  pr_rvecs(log,0,"X na do_relax",&(x[start]),homenr);
  pr_rvecs(log,0,"F na do_relax",&(f[start]),homenr);
  */

  return count; 
}

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,bool bDummies,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *md,
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,t_forcerec *fr,rvec box_size)
{
  t_mdebin   *mdebin;
  FILE       *ene;
  int        fp_ene,step,nre,k,n,count;
  double     tcount;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bTYZ,bLR,bBHAM,b14;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  real       etmp[F_NRE];
  real       mu_aver,mu_aver2;
  char       *traj;
  char       *xtc_traj; /* compressed trajectory filename */
  int        nDLB;
  int        i,m;
  rvec       vcm,mu_tot;
  t_coupl_rec *tcr;
  rvec       *xx,*vv,*ff;  
  bool       bTCR;

  /* Shell stuff */
  int         nshell;
  t_shell     *shells;
    
  where();

  /* Initial values */
  t = t0           = parm->ir.init_t;
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
  tcount       = 0;
  
  where();

  /* Check Environment variables */
  bTCR=ftp2bSet(efGCT,nfile,fnm);
  if (MASTER(cr)) {
    if (bTCR)
      fprintf(stderr,"Will do General Coupling Theory!\n");
    else
      fprintf(stderr,"*NO* General Coupling Theory ! ? !\n");
  }

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
  
  /* Check Environment variables */
  bTYZ=getenv("TYZ") != NULL;
  
  init_nrnb(&mynrnb);
  
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  
  fprintf(log,"Removing pbc first time\n");
  mk_mshift(log,graph,parm->box,x);
  shift_self(graph,fr->shift_vec,x);
  fprintf(log,"Done rmpbc\n");
  
  traj     = ftp2fn(efTRR,nfile,fnm);
  xtc_traj = ftp2fn(efXTC,nfile,fnm);
  
  bLR    = (parm->ir.rlong > parm->ir.rshort);
  bBHAM  = (top->idef.functype[0]==F_BHAM);
  b14    = (top->idef.il[F_LJ14].nr > 0);

  if (MASTER(cr)) {
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w");
    mdebin=init_mdebin(fp_ene,grps,&(top->atoms),bLR,bBHAM,b14);
  } else {
    fp_ene = -1;
    mdebin = NULL;
  }
  
  /*  init_dummies(log,md,START(nsb),HOMENR(nsb)); */
  where();
  
  clear_rvec(vcm);
  
  if (!parm->ir.bUncStart) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,md,x,vold,buf,f,v,
		  graph,cr,nrnb,grps,fr,top,edyn);
  
  /* Compute initial EKin for all.. */
  calc_ke_part(TRUE,0,top->atoms.nr,
               vold,v,vt,&(parm->ir.opts),
               md,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,
	 0,parm->ir.ntcmemory);
  where();
  
  shells = init_shells(log,START(nsb),HOMENR(nsb),&top->idef,md,&nshell);
  where();
  
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->pid,"Started mdrun");
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    printf("starting mdrun '%s'\n%d steps, %8.1f ps.\n\n",*(top->name),
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

    /* Construct dummy particles */
    construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
    /* Calculate total dipole moment if necessary */
    calc_mu(nsb,x,md->chargeT,mu_tot);

    /* Now is the time to relax the shells */
    count=relax_shells(ene,log,cr,bVerbose,step,parm,
		       bNS,bStopCM,top,ener,
		       x,vold,v,vt,f,buf,md,nsb,&mynrnb,
		       graph,grps,force_vir,nshell,shells,fr,
		       traj,0,0,nsb->natoms,parm->box,mdebin);
    tcount+=count;

    if (bTCR && MASTER(cr) && (step == 0)) 
      tcr=init_coupling(log,nfile,fnm,cr,fr,md,&(top->idef));

    /* Spread the force on dummy particle to the other particles... */
    spread_dummy_f(log,x,f,&mynrnb,&top->idef);

    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. 
     * We must output all of this before the update.
     */
    t=t0+step*parm->ir.delta_t;
    if (parm->ir.bPert)
      lambda   = lam0 + step*parm->ir.delta_lambda;
    if (parm->ir.bSimAnn) {
      SAfactor = 1.0  - t/parm->ir.zero_temp_time;
      if (SAfactor < 0) 
	SAfactor = 0;
    }
    
    if (do_per_step(step,parm->ir.nstxout)) xx=x; else xx=NULL;
    if (do_per_step(step,parm->ir.nstvout)) vv=v; else vv=NULL;
    if (do_per_step(step,parm->ir.nstfout)) ff=f; else ff=NULL;
    write_traj(log,cr,traj,
               nsb,step,t,lambda,&mynrnb,nsb->natoms,xx,vv,ff,parm->box);
    where();

    if (do_per_step(step,parm->ir.nstxtcout)) {
      write_xtc_traj(log,cr,xtc_traj,nsb,md,
                     step,t,x,parm->box,parm->ir.xtcprec);
      where();
    }

    where();
    clear_mat(shake_vir);
    update(nsb->natoms,START(nsb),HOMENR(nsb),
	   0,lambda,&ener[F_DVDL],&(parm->ir),FALSE,
           md,x,graph,
           fr->shift_vec,f,buf,vold,v,vt,parm->pres,parm->box,
           top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn);
    if (PAR(cr)) 
      accumulate_u(cr,&(parm->ir.opts),grps);
 
    /* Calculate partial Kinetic Energy (for this processor) 
     * per group!
     */
    calc_ke_part(FALSE,START(nsb),HOMENR(nsb),
                 vold,v,vt,&(parm->ir.opts),
                 md,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
    where();
    if (bStopCM)
      calc_vcm(log,HOMENR(nsb),START(nsb),md->massT,v,vcm);
    
    /* Copy the partial virial to the global virial (to be summed) */
    if (PAR(cr)) {
      global_stat(log,cr,ener,force_vir,shake_vir,
                  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot);
      if (!bNS)
        for(i=0; (i<grps->estat.nn); i++)
          grps->estat.ee[egLR][i] /= cr->nprocs;
    }
    else 
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
    if (bStopCM) {
      check_cm(log,vcm,md->tmass);
      do_stopcm(log,HOMENR(nsb),START(nsb),v,vcm,md->tmass);
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
#ifdef DEBUG
    fprintf(stderr,"Ekin 1: %14.10e", ener[F_EKIN]);
#endif
    /* Calculate Temperature coupling parameters lambda */
    tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,
	   step,parm->ir.ntcmemory);
    
    /* Calculate pressure ! */
    calc_pres(parm->box,parm->ekin,parm->vir,parm->pres,
	      EEL_LR(fr->eeltype) ? ener[F_LR] : 0.0);

    /* Calculate long range corrections to pressure and energy */
    if (bTCR)
      set_avcsix(log,fr,&top->idef,md);
    calc_ljcorr(log,parm->ir.bLJcorr,
		fr,md->nr,parm->box,parm->pres,parm->vir,ener);
    
    if (bTCR)
      do_coupling(log,nfile,fnm,tcr,t,step,ener,fr,
		  &(parm->ir),MASTER(cr) || bMultiSim,md,&(top->idef),mu_aver,
		  top->blocks[ebMOLS].nr,bMultiSim ? cr_msim : cr,
		  parm->box,parm->vir);

    upd_mdebin(mdebin,md->tmass,step,ener,parm->box,shake_vir,
               force_vir,parm->vir,parm->pres,grps,mu_tot);
               
    where();
    if ((MASTER(cr) && do_per_step(step,parm->ir.nstprint))) {
      print_ebin(fp_ene,log,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    }
    if (bVerbose)
      fflush(log);
      
    if (MASTER(cr) && bVerbose && ((step % stepout)==0))
      print_time(stderr,start_t,step,&parm->ir);
  }
  t=t0+step*parm->ir.delta_t;
  lambda=lam0+step*parm->ir.delta_lambda;
  
  if (MASTER(cr)) {
    if (parm->ir.nstprint > 1)
      print_ebin(fp_ene,log,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));

    print_ebin(-1,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(-1,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));    
  }
  fprintf(log,"Average number of force evaluations per MD step: %.2f\n",
	  tcount/parm->ir.nsteps);

  return start_t;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "shell_md is the polarizable MD program"
  };
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

  int          i;
  t_edsamyn edyn;

  t_commrec    *cr;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL,      ffREAD },
    { efTRR, "-o", NULL,      ffWRITE },
    { efXTC, "-x", NULL,      ffOPTWR },
    { efSTO, "-c", "confout", ffWRITE },
    { efGCT, "-j", "wham",    ffOPTRD },
    { efGCT, "-jo", "bam",    ffOPTRD },
    { efHAT, "-hat","ghat",    ffOPTRD },
    { efENX, "-e", "ener",    ffWRITE },
    { efLOG, "-g", "md",      ffWRITE },
    { efNDX, "-n",  "mols",    ffREAD }
  };
#define NFILE asize(fnm)

  stdlog = stderr;
  
  cr = init_par(&argc,argv);
  bVerbose = bVerbose && MASTER(cr);
  edyn.bEdsam = FALSE;

  if (MASTER(cr)) 
    CopyRight(stderr,argv[0]);

    parse_common_args(&argc,argv,
		      PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | 
		      (MASTER(cr) ? 0 : PCA_QUIET),
		      TRUE,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  	      
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Berendsen95a");
  }
  where();
  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,FALSE,nstepout,&edyn);
  
  exit(0);
  
  return 0;
}
