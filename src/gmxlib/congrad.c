/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
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
 * Great Red Oystrich Makes All Chemists Sane
 */
#include <string.h>
#include <time.h>
#include <math.h>
#include "sysstuff.h"
#include "string2.h"
#include "led.h"
#include "network.h"
#include "confio.h"
#include "binio.h"
#include "copyrite.h"
#include "smalloc.h"
#include "nrnb.h"
#include "main.h"
#include "pbc.h"
#include "force.h"
#include "macros.h"
#include "random.h"
#include "names.h"
#include "stat.h"
#include "fatal.h"
#include "txtdump.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "statutil.h"
#include "tgroup.h"
#include "mdebin.h"
#include "mdrun.h"
#include "congrad.h"

static void sp_header(FILE *out,real epot,real ftol)
{
  fprintf(out,"Conjugate gradients:\n");
  fprintf(out,"   Tolerance         = %12.5e\n",ftol);
  fprintf(out,"   Starting Energy   = %20.15e\n",epot);
}

time_t do_cg(FILE *log,int nfile,t_filenm fnm[],
		t_parm *parm,t_topology *top,
		t_groups *grps,t_nsborder *nsb,
		rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		tensor ekin,real ener[],
		t_nrnb nrnb[],
		bool bVerbose,t_commrec *cr,t_graph *graph)
{
  t_forcerec *fr;
  static char *CG="Conjugate Gradients";
  real   step0,lambda,ftol,fmax,testf,zet,w,smin;
  rvec   *p,*f,*xprime;
  rvec   *xx,*ff,box_size;
  real   EpotA=0.0,EpotB=0.0,a=0.0,b,beta=0.0,gpa,gpb;
  real   vcm[4],fnorm,pnorm,fnorm_old;
  FILE       *ene;
  t_mdebin   *mdebin;
  t_nrnb mynrnb;
  bool   bNS=TRUE,bDone,bLR,bBHAM,b14,bRand,brerun;
  time_t start_t;
  tensor force_vir,shake_vir;
  int    number_steps,naccept=0,nstcg=parm->ir.userint1;
  int    count=0;
  int    i,m,start,end,niti;

  /* Initiate some variables */
  if (parm->ir.bPert)
    lambda       = parm->ir.init_lambda;
  else 
    lambda = 0.0;
  fr=mk_forcerec();
  init_forcerec(log,fr,&(parm->ir),&(top->blocks[ebMOLS]),cr,
		&(top->blocks[ebCGS]),&(top->idef),mdatoms,parm->box,FALSE);
  for(m=0; (m<DIM); m++)
    box_size[m]=parm->box[m][m];
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  
  vcm[0]=vcm[1]=vcm[2]=vcm[3]=0.0;
  
  /* Print to log file */
  start_t=print_date_and_time(log,cr->pid,"Started Conjugate Gradients");
  
  /* p is the search direction, f the force, xprime the new positions */
  snew(p,nsb->natoms);
  snew(f,nsb->natoms);
  snew(xprime,nsb->natoms);

  start=nsb->index[cr->pid];
  end=nsb->homenr[cr->pid]-start;

  /* Open the energy file */  
  if (MASTER(cr))
    ene=ftp2FILE(efENE,nfile,fnm,"w");
  else
    ene=NULL;
    
  /* Set some booleans for the epot routines */
  bLR=(parm->ir.rlong > parm->ir.rshort);   /* Long Range Coulomb   ? */
  bBHAM=(top->idef.functype[0]==F_BHAM);    /* Use buckingham       ? */
  b14=(top->idef.il[F_LJ14].nr > 0);        /* Use 1-4 interactions ? */

  /* Init bin for energy stuff */
  mdebin=init_mdebin(ene,grps,&(top->atoms),bLR,bBHAM,b14);
  
  /* Clear some matrix variables */
  clear_mat(force_vir);
  clear_mat(shake_vir);
  
  /* Set variables for stepsize (in nm). This is the largest 
   * step that we are going to make in any direction.
   */
  step0=parm->ir.em_stepsize;
  
  /* Tolerance for convergence */
  ftol=parm->ir.em_tol;
  
  /* Max number of steps */
  number_steps=parm->ir.nsteps;

  /* Call the force routine and some auxiliary (neighboursearching etc.) */
  do_force(log,cr,
	   parm,nsb,force_vir,
	   0,&(nrnb[cr->pid]),top,grps,
	   x,buf,f,buf,
	   mdatoms,ener,bVerbose && !(PAR(cr)),
	   lambda,graph,bNS,FALSE,fr);
  unshift_self(graph,fr->shift_vec,x);
  where();

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  where();
  
  /* Clear var. */
  clear_mat(force_vir);
  where();
  
  /* Communicat energies etc. */
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm);
  where();
  
  /* Copy stuff to the energy bin for easy printing etc. */
  upd_mdebin(mdebin,mdatoms->tmass,count,ener,parm->box,shake_vir,
	     force_vir,parm->vir,parm->pres,grps);
  where();
  	
  /* Print only if we are the master processor */
  if (MASTER(cr))
    print_ebin(ene,log,count,count,lambda,0.0,eprNORMAL,TRUE,
	       mdebin,grps,&(top->atoms));
  where();
  
  /* This is the starting energy */
  EpotA=ener[F_EPOT];

  if (MASTER(cr)) {
    /* Print to the screen */
    start_t=print_date_and_time(log,cr->pid,"Started EM");
    sp_header(stderr,EpotA,ftol);
    sp_header(log,EpotA,ftol);
  }

  /* normalising step size, this saves a few hundred steps in the
   * beginning of the run.
   */
  fnorm=f_norm(log,cr->left,cr->right,nsb->nprocs,start,end,f);
  fnorm_old=fnorm;
  
  /* Print stepsize */
  if (MASTER(cr)) {
    fprintf(stderr,"   F-Norm            = %12.5e\n",fnorm);
    fprintf(stderr,"\n");
  }
  
  /* Here starts the loop, count is the counter for the number of steps
   * bDone is a BOOLEAN variable, which will be set TRUE when
   * the minimization has converged.
   */
  for(count=1,bDone=FALSE; 
      !(bDone || ((number_steps > 0) && (count==number_steps))); count++) {

    /* start conjugate gradient, determine search interval a,b */
    gpa=0.0;
    for(i=start; i<end; i++) {
      for(m=0;m<DIM;m++){
	p[i][m]=f[i][m]+ beta*p[i][m];
	gpa=gpa-p[i][m]*f[i][m];
      }
    }
    pnorm=f_norm(log,cr->left,cr->right,nsb->nprocs,start,end,p);

    a=0.0;
    b=step0/pnorm;
    niti=0;

    /* search a,b iteratively, if necessary */
    brerun=TRUE;
    while (brerun) {
      for (i=start;i<end;i++) {
	for(m=0;m<DIM;m++) {
	  xprime[i][m]=x[i][m] + b*p[i][m];
	}
      }
      bNS=((parm->ir.nstlist && ((count % parm->ir.nstlist)==0)) || (count==0)); 
      /* Calc force & energy on new trial position  */
      do_force(log,cr,parm,nsb,force_vir,
	       count,&(nrnb[cr->pid]),top,grps,xprime,buf,f,
	       buf,mdatoms,ener,bVerbose && !(PAR(cr)),
	       lambda,graph,bNS,FALSE,fr);
      unshift_self(graph,fr->shift_vec,x);
      bNS=FALSE;
      gpb=0.0;
      for(i=start;i<end;i++) {
	for(m=0;m<DIM;m++) {
	  gpb=gpb-p[i][m]*f[i][m];
	}
      } 
      
      /* Sum the potential energy terms from group contributions */
      sum_epot(&(parm->ir.opts),grps,ener);
      
      /* Clear stuff again */
      clear_mat(force_vir);
      clear_mat(shake_vir);
      
      /* Communicate stuff when parallel */
      if (PAR(cr)) 
	global_stat(log,cr,ener,force_vir,shake_vir,
		    &(parm->ir.opts),grps,&mynrnb,nrnb,vcm);

      EpotB=ener[F_EPOT];
      
      if ((gpb >= 0.0) || (EpotB >= EpotA))
	brerun=FALSE;
      else {
	a=b;
	EpotA=EpotB;
	gpa=gpb;
	b+=b;
      }
     niti++;
    }

    /* find stepsize smin in interval a-b */
    zet = 3.0 * (EpotA-EpotB) / (b-a) + gpa + gpb;
    w = zet*zet - gpa*gpb;
    if (w < 0.0) {
      fprintf(stderr,"Negative w: %20.12e\n",w);
      fprintf(stderr,"z= %20.12e\n",zet);
      fprintf(stderr,"gpa= %20.12e, gpb= %20.12e\n",gpa,gpb);
      fprintf(stderr,"a= %20.12e, b= %20.12e\n",a,b);
      fprintf(stderr,"EpotA= %20.12e, EpotB= %20.12e\n",EpotA,EpotB);
      fatal_error(0,"Negative number for sqrt encountered (%f)",w);
    }      
    w = sqrt(w);
    smin = b - ((gpb+w-zet)*(b-a))/((gpb-gpa)+2.0*w);

    /* new positions */
    for (i=start;i<end;i++){
      for(m=0;m<DIM;m++){
	xprime[i][m]=x[i][m] + smin*p[i][m];
      }
    }
    /* new energy, forces */
    do_force(log,cr,parm,nsb,force_vir,
	     count,&(nrnb[cr->pid]),top,grps,xprime,buf,f,
	     buf,mdatoms,ener,bVerbose && !(PAR(cr)),
	     lambda,graph,bNS,FALSE,fr);
    unshift_self(graph,fr->shift_vec,x);
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener); 
    fnorm=f_norm(log,cr->left,cr->right,nsb->nprocs,start,end,f);

    /* Clear stuff again */
      clear_mat(force_vir);
      clear_mat(shake_vir);
      
      /* Communicate stuff when parallel */
      if (PAR(cr)) 
	global_stat(log,cr,ener,force_vir,shake_vir,
		    &(parm->ir.opts),grps,&mynrnb,nrnb,vcm);

    EpotA=ener[F_EPOT];

    /* new search direction */
    /* beta = 0 means steepest descents */
    if (nstcg && ((count % nstcg)==0)) 
      beta = 0.0;
    else
      beta=fnorm*fnorm/(fnorm_old*fnorm_old);

    /* update x, fnorm_old */
    for (i=start;i<end;i++)
      copy_rvec(xprime[i],x[i]);
    fnorm_old=fnorm;

    /* Test whether the convergence criterion is met */
    fmax=f_max(log,cr->left,cr->right,nsb->nprocs,start,end,f);     

    /* Print it if necessary */
    if (bVerbose && MASTER(cr)) {
      fprintf(stderr,"\rStep %d, E-Pot = %16.10e, F-max = %12.5e\n",
	      count,EpotA,fmax);
      /* Store the new (lower) energies */
      upd_mdebin(mdebin,mdatoms->tmass,count,ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps);
      /* Print the energies allways when we should be verbose */
      if (MASTER(cr))
	print_ebin(ene,log,count,count,lambda,0.0,eprNORMAL,TRUE,
		   mdebin,grps,&(top->atoms));
    }
    
    /* Stop when the maximum force lies below tolerance */
    bDone=(fmax < ftol);

    /*if (MASTER(cr))
      fprintf(stderr,"\n");*/
	
  } /* End of the loop */

  /* Print some shit... */
  if (MASTER(cr)) {
    fprintf(stderr,"\nwriting lowest energy coords to traj...\n");
    xx=x;
    ff=f;
    write_traj(log,cr,ftp2fn(efTRJ,nfile,fnm),
	       nsb,count,(real) count,
	       lambda,nrnb,nsb->natoms,xx,xx,ff,parm->box);
    fmax=f_max(log,cr->left,cr->right,nsb->nprocs,start,end,f);
    fprintf(stderr,"Maximum force: %12.5e\n",fmax);
    if (bDone) {
      fprintf(stderr,"\n%s converged to %8.6f \n",CG,ftol);
      fprintf(log,"%s converged to %8.6f \n",CG,ftol);
    }
    else {
      fprintf(stderr,"\n%s did not converge in %d steps\n",CG,number_steps);
      fprintf(log,"%s did not converge in %d steps\n",CG,number_steps);
    }
    fprintf(stderr,"  Function value at minimum = %12.4e\n",EpotA);
    fprintf(log,"  Function value at minimum = %12.4e\n",EpotA);
  }
  
  /* To print the actual number of steps we needed somewhere */
  parm->ir.nsteps=count;

  return start_t;
} /* That's all folks */

