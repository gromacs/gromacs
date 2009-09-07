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
 * Great Red Oystrich Makes All Chemists Sane
 */

#define FORCE_CRIT

#include <string.h>
#include <time.h>
#include <math.h>
#include "sysstuff.h"
#include "string2.h"
#include "network.h"
#include "confio.h"
#include "copyrite.h"
#include "smalloc.h"
#include "nrnb.h"
#include "main.h"
#include "pbc.h"
#include "force.h"
#include "macros.h"
#include "random.h"
#include "names.h"
#include "fatal.h"
#include "txtdump.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "enxio.h"
#include "tgroup.h"
#include "mdebin.h"
#include "mdrun.h"
#include "pppm.h"
#include "dummies.h"

#ifdef FORCE_CRIT
static void sp_header(FILE *out,real epot,real fsqrt,real step,real ftol)
{
  fprintf(out,"STEEPEST DESCENTS:\n");
  fprintf(out,"   Stepsize          = %12.5e\n",step);
  fprintf(out,"   Tolerance         = %12.5e\n",ftol);
  fprintf(out,"   Starting rmsF     = %30.20e\n",fsqrt);
  fprintf(out,"   Starting Energy   = %30.20e\n",epot);
}
#else
static void sp_header(FILE *out,real epot,real step,real ftol)
{
  fprintf(out,"STEEPEST DESCENTS:\n");
  fprintf(out,"   Stepsize          = %12.5e\n",step);
  fprintf(out,"   Tolerance         = %12.5e\n",ftol);
  fprintf(out,"   Starting Energy   = %30.20e\n",epot);
}
#endif

real f_max(FILE *log,
	   int left,int right,int nprocs,
	   int start,int end,rvec grad[])
{
  real fmax,fmax_0,fam;
  int  i,m;

  /* This routine finds the largest force component (abs value)
   * and returns it. On parallel machines the global max is taken.
   */
  fmax = fabs(grad[start][0]);
  
  for(i=start; (i<end); i++)
    for(m=0; (m<DIM); m++) {
      fam=fabs(grad[i][m]);
      fmax=max(fmax,fam);
    }
  if (nprocs == 1)
    return fmax;

  for(i=0; (i<nprocs-1); i++) {
    gmx_tx(left,(void *)&fmax,sizeof(fmax));
    gmx_rx(right,(void *)&fmax_0,sizeof(fmax_0));
    gmx_wait(left,right);
    if (fmax_0 > fmax)
      fmax=fmax_0;
  }

  return fmax;
}

real f_norm(FILE *log,
	   int left,int right,int nprocs,
	   int start,int end,rvec grad[])
{
  real fnorm;
  int  i,m;

  /* This routine finds the norm of the force
   * and returns it. Parallel machines not supported.
   */
  fnorm = 0;
  
  for(i=start; (i<end); i++) 
    for(m=0; (m<DIM); m++) { 
      fnorm=fnorm+grad[i][m]*grad[i][m]; 
    } 
  fnorm=sqrt(fnorm); 
  
  if (nprocs == 1) 
    return fnorm; 
  
  fatal_error(0,"This version of Steepest Descents cannot be run in parallel"); 
  return fnorm; 
} 

real f_sqrt(FILE *log, 
 	    int left,int right,int nprocs, 
 	    int start,int end,rvec f[]) 
{ 
  int i,m; 
  real fsqr; 
  
  
  /* this routine calculates the route mean square force  */
  
  fsqr=0; 
  
  /* calculate the sum of the square force of this processor  */
  for(i=start;(i<end);i++) 
    for(m=0;(m<DIM);m++)  
      fsqr += sqr(f[i][m]); 
  
  
  if (nprocs == 1) 
    return sqrt(fsqr)/(end - start); 
  
  
  fatal_error(0,"This version of f_sqrt cannot run in parrallel"); 
 

  return ( sqrt(fsqr) );
} 

static void do_step(int start,int end,rvec x[],rvec f[], 
 		    bool bRandom,real step) 
{ 
  static int seed=1993; 
  int    i,m; 
  real   r; 
  
  if (bRandom) { 
    fprintf(stderr,"\rRandom step\n"); 
    for(i=start; (i<end); i++) { 
      for(m=0; (m<DIM); m++) { 
 	r=rando(&seed); 
 	x[i][m]=x[i][m]+step*r; 
      } 
    } 
  } 
  else { 
    /* New positions to try  */
    for(i=start; (i<end);i++) 
      for(m=0; (m<DIM); m++) 
 	x[i][m] = x[i][m]+step*f[i][m]; 
  } 
} 

time_t do_steep(FILE *log,int nfile,t_filenm fnm[], 
 		t_parm *parm,t_topology *top, 
 		t_groups *grps,t_nsborder *nsb, 
 		rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms, 
 		tensor ekin,real ener[],t_nrnb nrnb[], 
 		bool bVerbose,bool bDummies,t_commrec *cr,t_graph *graph,
		t_forcerec *fr,rvec box_size) 
{ 
  static char *SD="STEEPEST DESCENTS"; 
  real   step,lambda,ftol,fmax; 
  rvec   *pos[2],*force[2]; 
  rvec   *xx,*ff; 
#ifdef FORCE_CRIT 
  real   Fsqrt[2]; 
#endif 
  real   Epot[2]; 
  real   vcm[4],fnorm,ustep; 
  int        fp_ene; 
  t_mdebin   *mdebin; 
  t_nrnb mynrnb; 
  bool   bNS=TRUE,bDone,bLR,bLJLR,bBHAM,b14; 
  time_t start_t; 
  tensor force_vir,shake_vir; 
  rvec   mu_tot;
  int    number_steps;
  int    count=0; 
  int    i,m,start,end; 
  int    Min=0; 
  int    steps_accepted=0; 
#define  TRY (1-Min) 
  
  /* Initiate some variables  */
  if (parm->ir.bPert)
    lambda       = parm->ir.init_lambda;
  else 
    lambda = 0.0;

  clear_rvec(mu_tot);
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE); 
  
  vcm[0]=vcm[1]=vcm[2]=vcm[3]=0.0; 
  
  /* Print to log file  */
  print_date_and_time(log,cr->pid,"Started Steepest Descents"); 
  
  /* We need two coordinate arrays and two force arrays  */
  for(i=0; (i<2); i++) { 
    snew(pos[i],nsb->natoms); 
    snew(force[i],nsb->natoms); 
  } 

  start=nsb->index[cr->pid]; 
  end=nsb->homenr[cr->pid]-start; 
  
  /* Open the enrgy file */   
  if (MASTER(cr)) 
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w"); 
  else 
    fp_ene=-1; 
  
  /* Set some booleans for the epot routines  */
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);
  
  /* Init bin for energy stuff  */
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),bLR,bLJLR,
		     bBHAM,b14,parm->ir.bPert,parm->ir.epc,
		     parm->ir.bDispCorr); 
  
  /* Clear some matrix variables  */
  clear_mat(force_vir); 
  clear_mat(shake_vir); 
  
  /* Copy coord vectors to our temp array  */
  for(i=0; (i<nsb->natoms); i++) { 
    copy_rvec(x[i],pos[Min][i]); 
    copy_rvec(x[i],pos[TRY][i]); 
  } 
    
  /* Set variables for stepsize (in nm). This is the largest  
   * step that we are going to make in any direction. 
   */
  step=ustep=parm->ir.em_stepsize; 
  
  /* Tolerance for conversion  */
  ftol=parm->ir.em_tol; 
  
  /* Max number of steps  */
  number_steps=parm->ir.nsteps; 
  
  if (MASTER(cr)) { 
    /* Print to the screen  */
    print_date_and_time(log,cr->pid,"Started EM"); 
    fprintf(stderr,"STEEPEST DESCENTS:\n");
    fprintf(log,"STEEPEST DESCENTS:\n");
    fprintf(stderr,"   Tolerance         = %12.5g\n",ftol); 
    fprintf(log,"   Tolerance         = %12.5g\n",ftol); 
  } 
    
  /* Here starts the loop, count is the counter for the number of steps 
   * bDone is a BOOLEAN variable, which will be set TRUE when 
   * the minimization has converged. 
   */
  for(count=0,bDone=FALSE;  
      !(bDone || ((number_steps > 0) && (count==number_steps))); count++) { 
    
    /* set new coordinates, except for first step */
    if (count>0)
      for(i=start; (i<end); i++)  
	for(m=0;(m<DIM);m++) 
	  pos[TRY][i][m] = pos[Min][i][m] + step * force[Min][i][m]; 
    
    if (bDummies) {
      /* Construct dummy particles */
      shift_self(graph,fr->shift_vec,pos[TRY]);
      construct_dummies(log,pos[TRY],&(nrnb[cr->pid]),1,NULL,&top->idef);
      unshift_self(graph,fr->shift_vec,pos[TRY]);
    }
    
    /* Calc force & energy on new positions  */
    do_force(log,cr,parm,nsb,force_vir, 
 	     count,&(nrnb[cr->pid]),top,grps,pos[TRY],buf,force[TRY],buf,
	     mdatoms,ener,bVerbose && !(PAR(cr)), 
 	     lambda,graph,bNS,FALSE,fr); 
    unshift_self(graph,fr->shift_vec,pos[TRY]);

    /* Spread the force on dummy particle to the other particles... */
    spread_dummy_f(log,pos[TRY],force[TRY],&(nrnb[cr->pid]),&top->idef);
    
    /* Sum the potential energy terms from group contributions  */
    sum_epot(&(parm->ir.opts),grps,ener); 
    
    /* Clear stuff again  */
    clear_mat(force_vir); 
    clear_mat(shake_vir); 
    
    /* Communicat stuff when parallel  */
    if (PAR(cr))  
      global_stat(log,cr,ener,force_vir,shake_vir, 
 		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot); 
    
    /* This is the new energy  */
#ifdef FORCE_CRIT 
    Fsqrt[TRY]=f_sqrt(log,cr->left,cr->right,nsb->nprocs,start,end,force[TRY]); 
#endif 
    Epot[TRY]=ener[F_EPOT]; 
    
    /* Print it if necessary  */
#ifdef FORCE_CRIT 
    if ((bVerbose || (Fsqrt[TRY] < Fsqrt[Min])) && MASTER(cr)) { 
      fprintf(stderr,
	      "\rStep = %5d, Dx = %12.5e, Epot = %12.5e rmsF = %12.5e\n", 
 	      count,step,Epot[TRY],Fsqrt[TRY]); 
      /* Store the new (lower) energies  */
      upd_mdebin(mdebin,mdatoms->tmass,count,ener,parm->box,shake_vir, 
 		 force_vir,parm->vir,parm->pres,grps,mu_tot); 
      /* Print the energies allways when we should be verbose  */
      if (MASTER(cr)) 
 	print_ebin(fp_ene,TRUE,FALSE,log,count,count,lambda,
		   0.0,eprNORMAL,TRUE,mdebin,grps,&(top->atoms)); 
    } 
#else 
    if ((bVerbose || (Epot[TRY] < Epot[Min])) && MASTER(cr)) { 
      fprintf(stderr,"\rStep = %5d, Dx = %12.5e, E-Pot = %30.20e\n", 
 	      count,step,Epot[TRY]); 
      
      /* Store the new (lower) energies  */
      upd_mdebin(mdebin,mdatoms->tmass,count,ener,parm->box,shake_vir, 
 		 force_vir,parm->vir,parm->pres,grps); 
      /* Print the energies allways when we should be verbose  */
      if (MASTER(cr)) 
 	print_ebin(fp_ene,log,count,count,lambda,0.0,eprNORMAL,TRUE, 
 		   mdebin,grps,&(top->atoms)); 
    } 
#endif 
    
    /* Now if the new energy is smaller than the previous...  
     * or if this is the first step!
     * or if we did random steps! 
     */
    
#ifdef FORCE_CRIT 
    if ( (count==0) || (Fsqrt[TRY] < Fsqrt[Min]) ) { 
#else 
    if ( (count==0) || (Epot[TRY] < Epot[Min]) ) { 
#endif 
      steps_accepted++; 
      if (do_per_step(steps_accepted,parm->ir.nstfout)) 
	ff=force[TRY];  
      else 
	ff=NULL;    
      if (do_per_step(steps_accepted,parm->ir.nstxout)) {
	xx=pos[TRY];   
	write_traj(log,cr,ftp2fn(efTRN,nfile,fnm), 
		   nsb,count,(real) count, 
		   lambda,nrnb,nsb->natoms,xx,NULL,ff,parm->box); 
      } else 
	xx=NULL; 
      
      /* Test whether the convergence criterion is met...  */
#ifdef FORCE_CRIT 
      bDone=(Fsqrt[TRY] < ftol);
#else 
      /* Stop when the difference between two tries is less than 
       * the tolerance times the energy. 
       */
      bDone=(fabs(Epot[Min]-Epot[TRY]) < ftol*fabs(Epot[Min])); 
#endif 
      
      /* Copy the arrays for force, positions and energy  */
      /* The 'Min' array always holds the coords and forces of the minimal 
	 sampled energy  */
      Min = TRY; 
      
      /* increase step size  */
      if (count>0)
	ustep *= 1.2; 
      
      /*if (MASTER(cr)) 
	fprintf(stderr,"\n"); */
    } 
    else { 
      /* If energy is not smaller make the step smaller...  */
      if (ustep > 1.0e-14)
	ustep *= 0.5;
      else
	ustep = 1.0e-14;
    }
     
    /* Determine new step  */
    fnorm = f_norm(log,cr->left,cr->right,nsb->nprocs,start,end,force[Min]); 
    step=ustep/fnorm;
    
  } /* End of the loop  */
    
  /* Print some shit...  */
  if (MASTER(cr)) { 
    fprintf(stderr,"\nwriting lowest energy coordinates.\n"); 
    xx=pos[Min]; 
    ff=force[Min]; 
    write_traj(log,cr,ftp2fn(efTRN,nfile,fnm), 
 	       nsb,count,(real) count, 
 	       lambda,nrnb,nsb->natoms,xx,NULL,ff,parm->box); 
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top->name, &(top->atoms),xx,NULL,parm->box);
    
    fmax=f_max(log,cr->left,cr->right,nsb->nprocs,start,end,force[Min]); 
    fprintf(stderr,"Maximum force: %12.5e\n",fmax); 
    if (bDone) { 
      fprintf(stderr,"\n%s converged to %8.6f \n",SD,ftol); 
      fprintf(log,"%s converged to %8.6f \n",SD,ftol); 
    } 
    else { 
      fprintf(stderr,"\n%s did not converge in %d steps\n",SD,number_steps); 
      fprintf(log,"%s did not converge in %d steps\n",SD,number_steps); 
    } 
#ifdef FORCE_CRIT 
    fprintf(stderr,"  Minimum Root Mean Square Force  = %12.4e\n",Fsqrt[Min]); 
    fprintf(log,"  Minimum Root Mean Square Force = %12.4e\n",Fsqrt[Min]); 
#endif
    fprintf(stderr,"  Function value at minimum = %12.4e\n",Epot[Min]); 
    fprintf(log,"  Function value at minimum = %12.4e\n",Epot[Min]); 
  }
  if (MASTER(cr))
    close_enx(fp_ene);
    
  /* Put the coordinates back in the x array (otherwise the whole
   * minimization would be in vain)
   */
  for(i=0; (i<nsb->natoms); i++)
    copy_rvec(pos[Min][i],x[i]);
  
  /* To print the actual number of steps we needed somewhere */
  parm->ir.nsteps=count;
  
  return start_t;
} /* That's all folks */

