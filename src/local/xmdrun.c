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
 * S  C  A  M  O  R  G
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
#include "ebin.h"
#include "mdebin.h"
#include "disre.h"
#include "dummies.h"
#include "init_sh.h"
#include "do_gct.h"
#include "physics.h"
#include "block_tx.h"
#include "rdgroup.h"
#include "glaasje.h"
#include "edsam.h"
#include "calcmu.h"
#include "ionize.h"

void Etcoupl(bool bTC,t_grpopts *opts,t_groups *grps,real dt,real SAfactor,
	     int step,int nmem)
{
  static real *Told=NULL;
  int    i;
  real   T,reft,lll;

  if (!Told) {
    snew(Told,opts->ngtc);
    for(i=0; (i<opts->ngtc); i++) 
      Told[i]=opts->ref_t[i]*SAfactor;
  }
  
  for(i=0; (i<opts->ngtc); i++) {
    reft=opts->ref_t[i]*SAfactor;
    if (reft < 0)
      reft=0;
    
    Told[i] = run_aver(Told[i],grps->tcstat[i].T,step,nmem);
    T       = Told[i];
    
    if ((bTC) && (T != 0.0)) {
      lll=sqrt(1.0 + (dt/opts->tau_t[i])*(reft/T-1.0));
      grps->tcstat[i].lambda=max(min(lll,1.25),0.8);
    }
    else
      grps->tcstat[i].lambda=1.0;
#ifdef DEBUGTC
    fprintf(stdlog,"group %d: T: %g, Lambda: %g\n",
	    i,T,grps->tcstat[i].lambda);
#endif
  }
}

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

static void shell_pos_sd(FILE *log,real step,
			 rvec xold[],rvec xnew[],rvec f[],
			 int nas,t_atomshell as[],
			 int nbs,t_bondshell bs[])
{
  int  i,shell;
  real k_1;
  
  for(i=0; (i<nas); i++) {
    shell = as[i].shell;
    k_1   = as[i].k_1;
    do_1pos(xnew[shell],xold[shell],f[shell],k_1,step);
  }
  for(i=0; (i<nbs); i++) {
    shell = bs[i].shell;
    k_1   = bs[i].k_1;
    do_1pos(xnew[shell],xold[shell],f[shell],k_1,step);
  }
}


static void shell_pos_sd2(FILE *log,real step,
			  rvec xold[],rvec xnew[],rvec fmin[],rvec ftry[],
			  int nas,t_atomshell as[],
			  int nbs,t_bondshell bs[])
{
  int  i,m,shell;
  real k_1;
  
  for(i=0; (i<nas); i++) {
    shell = as[i].shell;
    k_1   = as[i].k_1;
    rvec_inc(fmin[shell],ftry[shell]);
    do_1pos(xnew[shell],xold[shell],fmin[shell],k_1,step);
  }
  for(i=0; (i<nbs); i++) {
    shell = bs[i].shell;
    k_1   = bs[i].k_1;
    rvec_inc(fmin[shell],ftry[shell]);
    do_1pos(xnew[shell],xold[shell],fmin[shell],k_1,step);
  }
}


static void predict_shells(FILE *log,rvec x[],rvec v[],real dt,
			   int nas,t_atomshell as[],
			   int nbs,t_bondshell bs[])
{
  int  i,m;
  real hdt;
  
  for(i=0; (i<nas); i++)
    for(m=0; (m<DIM); m++)
      x[as[i].shell][m]+=v[as[i].nucl1][m]*dt;
  
  hdt=0.5*dt;
  for(i=0; (i<nbs); i++) 
    for(m=0; (m<DIM); m++)
      x[bs[i].shell][m]+=(v[bs[i].nucl1][m]+v[bs[i].nucl2][m])*hdt;
}

static void print_epot(char *which,
		       int mdstep,int count,real step,real epot,real df)
{
  fprintf(stderr,"MDStep=%5d/%2d lamb: %10.5f, E-Pot: %12.8e",
	  mdstep,count,step,epot);
  
  if (count != 0)
    fprintf(stderr,", rmsF: %12.8e %s\n",df,which);
  else
    fprintf(stderr,"\n");
}


static real rms_force(rvec f[],
		      int nas,t_atomshell as[],
		      int nbs,t_bondshell bs[])
{
  int  i,shell;
  real df2;
  
  df2=0.0;
  for(i=0; (i<nas); i++) {
    shell = as[i].shell;
    df2  += iprod(f[shell],f[shell]);
  }
  for(i=0; (i<nbs); i++) {
    shell = bs[i].shell;
    df2  += iprod(f[shell],f[shell]);
  }
  return sqrt(df2/(nas+nbs));
}

static int relax_shells(FILE *ene,FILE *log,t_commrec *cr,
			bool bVerbose,int mdstep,
			t_parm *parm,bool bDoNS,bool bStopCM,
			t_topology *top,real ener[],
			rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
			rvec buf[],t_mdatoms *md,
			t_nsborder *nsb,t_nrnb *nrnb,
			t_graph *graph,
			t_groups *grps,tensor vir_part,
			int nas,t_atomshell as[],
			int nbs,t_bondshell bs[],
			t_forcerec *fr,
			char *traj,
			real t,real lambda,
			int natoms,matrix box,t_mdebin *mdebin)
{
  static bool bFirst=TRUE;
  static rvec *pos[3],*force[3];
  static real step;
#define NEPOT asize(Epot)
  real   ftol,s1,s2,eold,step0;
  rvec   EEE,abc,rmsF;
  matrix SSS = { { 0, 0, 1 }, { 1, 1, 1}, { 0, 1, 0 }};
  matrix Sinv;
  bool   bDone,bMinSet;
  int    g;
  int    number_steps;
  int    count=0;
  int    i,start=START(nsb),homenr=HOMENR(nsb),end=START(nsb)+HOMENR(nsb);
  int    Min=0;
  /* #define  Try (1-Min)  */           /* At start Try = 1 */
#define  Try1 ((Min+1) % 3)
#define  Try2 ((Min+2) % 3)

  if (bFirst) {
    /* Allocate local arrays */
    for(i=0; (i<3); i++) {
      snew(pos[i],nsb->natoms);
      snew(force[i],nsb->natoms);
    }
    bFirst=FALSE;
  }
  
  ftol  = parm->ir.em_tol;
  number_steps=parm->ir.userint4;
  step0 = parm->ir.userreal1;   
  step  = step0;
  
  /* Do a prediction of the shell positions */
  predict_shells(log,x,v,parm->ir.delta_t,nas,as,nbs,bs);
    
  /* Calculate the forces first time around */
  do_force(log,cr,parm,nsb,vir_part,mdstep,nrnb,
	   top,grps,x,v,force[Min],buf,md,ener,bVerbose && !PAR(cr),
	   lambda,graph,bDoNS,FALSE,fr);
  rmsF[0]=rms_force(force[Min],nas,as,nbs,bs);
  
  /* Copy x to pos[Min] & pos[Try1]: during minimization only the
   * shell positions are updated, therefore the other particles must
   * be set here.
   */
  memcpy(pos[Min],x,nsb->natoms*sizeof(x[0]));
  memcpy(pos[Try1],x,nsb->natoms*sizeof(x[0]));
  memcpy(pos[Try2],x,nsb->natoms*sizeof(x[0]));
  for(i=0; (i<nsb->natoms); i++) {
    clear_rvec(force[Try1][i]);
    clear_rvec(force[Try2][i]);
  }
  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  EEE[0] = ener[F_EPOT];
  
#ifdef DEBUG
  if (bVerbose && MASTER(cr))
    print_epot("",mdstep,0,step,Epot[Min],df[Min]);
  fprintf(stderr,"%17s: %14.10e\n",
	  interaction_function[F_EKIN].longname, ener[F_EKIN]);
  fprintf(stderr,"%17s: %14.10e\n",
	  interaction_function[F_EPOT].longname, ener[F_EPOT]);
  fprintf(stderr,"%17s: %14.10e\n",
      interaction_function[F_ETOT].longname, ener[F_ETOT]);
#endif
      
  bDone=((nas == 0) && (nbs == 0) || (rmsF[0] < ftol));
  bMinSet=FALSE;
  
  /****************************************************** 
   *  Start the shell-relaxation loop 
   ******************************************************/
  for(count=1; 
      !(bDone || ((number_steps > 0) && (count>=number_steps))); ) {
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try1],force[Min],nas,as,nbs,bs); 

#ifdef DEBUGHARD
    pr_rvecs(log,0,"pos[Try1] b4 do_force",&(pos[Try1][start]),homenr);
    pr_rvecs(log,0,"pos[Try2] b4 do_force",&(pos[Try2][start]),homenr);
    pr_rvecs(log,0,"pos[Min] b4 do_force",&(pos[Min][start]),homenr);
#endif
    
    /* Try the new positions */
    do_force(log,cr,parm,nsb,vir_part,1,nrnb,
	     top,grps,pos[Try1],v,force[Try1],buf,md,ener,bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr);
    count++;
    rmsF[1]=rms_force(force[Try1],nas,as,nbs,bs);
#ifdef DEBUGHARD
    pr_rvecs(log,0,"F na do_force",&(force[Try1][start]),homenr);
#endif

    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    EEE[1] = ener[F_EPOT];
    
    if (bVerbose && MASTER(cr))
      print_epot("",mdstep,count,step,EEE[1],rmsF[1]);
    bDone=(rmsF[1] < ftol);
    
    /* NOW! Do line mimization */
    EEE[2] = rmsF[0];
    
    m_inv(SSS,Sinv);
    mvmul(Sinv,EEE,abc);
	
    if ((abc[0] == 0) || (abc[0]*abc[1] > 0)) {
      pr_rvecs(log,0,"SSS",SSS,DIM);
      pr_rvecs(log,0,"Sinv",Sinv,DIM);
      pr_rvec(log,0,"EEE",EEE,DIM);
      pr_rvec(log,0,"abc",abc,DIM);
      fatal_error(0,"Relaxing shells: line minimization failed. Check log");
    }
    /* We know and checked that the solution for step must be > 0 because 
     * the new positions are in the direction of the forces
     * i.e. the first derivative of the energy is < 0 at step = 0
     */
    step = min(2.0*step0,-abc[1]/(2*abc[0]));
	
    /* New positions at half the step size, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try2],force[Min],nas,as,nbs,bs); 
      
    /* Try the new positions */
    do_force(log,cr,parm,nsb,vir_part,1,nrnb,
	     top,grps,pos[Try2],v,force[Try2],buf,md,ener,
	     bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr);
    count++;
      
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    
    eold    = EEE[0];
    EEE[0]  = ener[F_EPOT];
    rmsF[0] = rms_force(force[Try2],nas,as,nbs,bs);
    bDone   = (rmsF[0] < ftol);
      
    if (bVerbose && MASTER(cr)) {
      print_epot("",mdstep,count,step,EEE[0],rmsF[0]);
      fprintf(stderr,"DE1=%10g, DE0= %10g\n",eold-EEE[1],eold-EEE[0]);
    }
    step = step0;
    Min  = Try2;
  }
  if (MASTER(cr) && !bDone) 
    fprintf(stderr,"EM did not converge in %d steps\n",number_steps);
  
  /* Parallelise this one! */
  memcpy(x,pos[Min],nsb->natoms*sizeof(x[0]));
  memcpy(f,force[Min],nsb->natoms*sizeof(f[0]));
#ifdef DEBUGHARD
  pr_rvecs(log,0,"X na do_relax",&(x[start]),homenr);
  pr_rvecs(log,0,"F na do_relax",&(f[start]),homenr);
#endif

  return count; 
}

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *md,
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn)
{
  t_forcerec *fr;
  t_mdebin   *mdebin;
  FILE       *ene,*fmu;
  int        step,nre,k,n,count;
  double     tcount;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bTYZ,bLR,bBHAM,b14,bMU;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  real       etmp[F_NRE];
  char       strbuf[256];
  char       *traj;
  char       *xtc_traj; /* compressed trajectory filename */
  int        nDLB;
  int        i,m;
  rvec       box_size,vcm,mu_tot;
  t_coupl_rec *tcr;
  rvec       *xx,*vv,*ff;  
  bool       bTCR;
  real       mu_aver,mu_aver2;
  int        gnx;
  atom_id    *grpindex;
  char       *grpname;

  /* Shell stuff */
  int         nas;
  int         nbs;
  t_atomshell *as;
  t_bondshell *bs;
    
  /* Initial values */
  t = t0       = parm->ir.init_t;
  if (parm->ir.bSimAnn) {
    SAfactor = 1.0  - t0/parm->ir.zero_temp_time;
    if (SAfactor < 0) 
      SAfactor = 0;
  } else
    SAfactor     = 1.0;
  lam0         = parm->ir.init_lambda;
  nre          = F_NRE;
  tcount       = 0;
  
  where();

  /* check if "-mu" is used */
  bMU = opt2bSet("-mu",nfile,fnm);
  if ((bMU) && MASTER(cr))
    fmu = opt2FILE("-mu",nfile,fnm,"w");

  rd_index(ftp2fn(efNDX,nfile,fnm),1,&gnx,&grpindex,&grpname);

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
    fprintf(log,"Doing Free Energy Perturbation calculation. Lam0 = %g\n",
	    lam0);
  }
  else {
    lam0   = 0.0;
    lambda = 0.0;
  } 
  
  /* Check Environment variables */
  bTYZ=getenv("TYZ") != NULL;
  
  init_nrnb(&mynrnb);
  
  fr       = mk_forcerec();
  init_forcerec(log,fr,&(parm->ir),&(top->blocks[ebMOLS]),cr,
		&(top->blocks[ebCGS]),&(top->idef),md,parm->box,FALSE);
  for(m=0; (m<DIM); m++)
    box_size[m]=parm->box[m][m];
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  
  fprintf(log,"Removing pbc first time\n");
  mk_mshift(log,graph,parm->box,x);
  shift_self(graph,fr->shift_vec,x);
  fprintf(log,"Done rmpbc\n");
  
  traj     = ftp2fn(efTRJ,nfile,fnm);
  xtc_traj = ftp2fn(efXTC,nfile,fnm);
  
  if (MASTER(cr))
    ene=ftp2FILE(efENE,nfile,fnm,"w");
  else
    ene=NULL;
  bLR    = ((parm->ir.rlong > parm->ir.rshort) && 
	    (parm->ir.eeltype == eelTWIN));
  bBHAM  = (top->idef.functype[0]==F_BHAM);
  b14    = (top->idef.il[F_LJ14].nr > 0);
  mdebin = init_mdebin(ene,grps,&(top->atoms),bLR,bBHAM,b14);
  
  /*  init_dummies(log,md,START(nsb),HOMENR(nsb)); */
  where();
  
  clear_rvec(vcm);
  
  if (parm->ir.bShakeFirst) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,md,x,vold,buf,f,v,
		  graph,cr,nrnb,grps,fr,top,edyn);
  
  /* Compute initial EKin for all.. */
  calc_ke_part(TRUE,0,top->atoms.nr,
               vold,v,vt,&(parm->ir.opts),
               md,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
  /* Compute initial EKin for all.. dit zat er in gmx151 nog wel in???? 
  clear_mat(force_vir);
  clear_mat(shake_vir);
  calc_ke_part(TRUE,0,top->atoms.nr,
               vold,v,vt,&(parm->ir.opts),
               md,grps,&mynrnb,
	       lambda,&ener[F_DVDLKIN]);
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm);
  clear_rvec(vcm);
  */

  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  Etcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,
	  0,parm->ir.ntcmemory);
  where();
  
  init_shells(log,START(nsb),HOMENR(nsb),&top->idef,md,&nas,&as,&nbs,&bs);
  where();
  
  /* Write start time and temperature */
  sprintf(strbuf,"Starting %s",Program());
  start_t=print_date_and_time(log,cr->pid,strbuf);
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    printf("starting %s '%s'\n%d steps, %8.1f ps.\n\n",Program(),
	   *(top->name),parm->ir.nsteps,parm->ir.nsteps*parm->ir.delta_t);
  }

  /***********************************************************
   *
   *             Loop over MD steps 
   *
   ************************************************************/
  for (step=0; (step<parm->ir.nsteps); step++) {
    /* Stop Center of Mass motion */
    bStopCM=do_per_step(step,parm->ir.nstcomm);

    ionize(log,mdatoms,top->atoms.atomname,t,&parm->ir);

    /* Determine whether or not to do Neighbour Searching */
    bNS=((parm->ir.nstlist && ((step % parm->ir.nstlist)==0)) || (step==0));

    /* Construct dummy particles */
    construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);

    /* ??? Or should this be after relax_shells ??? */
    /* Calculate total dipole moment of the simulation box and the average
       dipole moment of the molecules */    

    mu_aver=calc_mu(cr,nsb,x,md->chargeA,mu_tot,top,md,gnx,grpindex);
    
    if (bMU && MASTER(cr))
      write_mu(fmu,mu_tot,parm->box);
    
    /* Now is the time to relax the shells */
    count=relax_shells(ene,log,cr,bVerbose,step,parm,
		       bNS,bStopCM,top,ener,
		       x,vold,v,vt,f,
		       buf,md,nsb,&mynrnb,
		       graph,grps,force_vir,
		       nas,as,nbs,bs,fr,
		       traj,t,lambda,nsb->natoms,parm->box,mdebin);
    tcount+=count;

    do_glas(log,START(nsb),HOMENR(nsb),x,f,fr,md,top->idef.atnr,
	    &parm->ir,ener);
	         
    if (bTCR && MASTER(cr) && (step == 0)) 
      tcr=init_coupling(log,nfile,fnm,cr,fr,md,&(top->idef));

    /* Spread the force on dummy particle to the other particles... */
    spread_dummy_f(log,x,f,&mynrnb,&top->idef);

    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. 
     * We must output all of this before the update.
     */
    t        = t0   + step*parm->ir.delta_t;
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
                  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm);
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
    Etcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,
	    step,parm->ir.ntcmemory);
    
    /* Calculate pressure ! */
    calc_pres(parm->box,parm->ekin,parm->vir,parm->pres);

    /* Calculate long range corrections to pressure and energy */
    calc_ljcorr(log,parm->ir.bLJcorr,
		fr,md->nr,parm->box,parm->pres,ener);
    
    if (bTCR)
      do_coupling(log,tcr,t,step,ener,fr,
		  &(parm->ir),MASTER(cr),md,&(top->idef),mu_aver);
    
    upd_mdebin(mdebin,md->tmass,step,ener,parm->box,shake_vir,
               force_vir,parm->vir,parm->pres,grps);
               
    where();
    if ((MASTER(cr) && do_per_step(step,parm->ir.nstprint))) {
      print_ebin(ene,log,step,t,lambda,SAfactor,eprNORMAL,bCompact,
                 mdebin,grps,&(top->atoms));
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
      print_ebin(ene,log,step-1,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    
    print_ebin(NULL,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(NULL,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));
  }
  fprintf(log,"Average number of force evaluations per MD step: %.2f\n",
	  tcount/parm->ir.nsteps);

  return start_t;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "xmdrun is the experimental MD program. New features are tested in this",
    "program before being implemented in the default mdrun. Currently under",
    "investigation are: polarizibility, glass and Free energy perturbation."
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

  t_commrec    *cr;
  t_filenm fnm[] = {
    { efTPB, NULL, NULL,      ffREAD },
    { efTRJ, "-o", NULL,      ffWRITE },
    { efXTC, "-x", NULL,      ffOPTWR },
    { efGRO, "-c", "confout", ffWRITE },
    { efGCT, "-j", "wham",    FALSE },
    { efGCT, "-jo", "bam",    FALSE },
    { efENE, "-e", "ener",    ffWRITE },
    { efLOG, "-g", "md",      ffWRITE },
    { efDAT, "-mu", "dipole", ffOPTWR },
    { efNDX, "-n", "mols",    ffREAD }
  };
#define NFILE asize(fnm)
  t_edsamyn edyn;
  
  stdlog = stderr;
  
  get_pargs(&argc,argv,asize(pa),pa,TRUE);
  cr = init_par(nprocs,argv);
  bVerbose = bVerbose && MASTER(cr);
  edyn.bEdsam = FALSE;
  
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

  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,FALSE,nstepout,&edyn);
  
  exit(0);
  
  return 0;
}
