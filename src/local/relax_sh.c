#include <stdio.h>
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "fatal.h"
#include "vec.h"
#include "txtdump.h"
#include "mdrun.h"
#include "init_sh.h"
#include "mdatoms.h"
#include "network.h"

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
    if (debug) {
      pr_rvec(debug,0,"fshell",f[shell],DIM);
      pr_rvec(debug,0,"xold",xold[shell],DIM);
      pr_rvec(debug,0,"xnew",xnew[shell],DIM);
    }
  }
}

static void predict_shells(FILE *log,rvec x[],rvec v[],real dt,int ns,t_shell s[],
			   real mass[],bool bInit)
{
  int  i,m,s1,n1,n2,n3;
  real dt_1,dt_2,dt_3,fudge,tm,m1,m2,m3;
  rvec *ptr;
  
  /* We introduce a fudge factor for performance reasons: with this choice
   * the initial force on the shells is about a factor of two lower than without
   */
  fudge = 1.0; /*/sqrt(3.0);*/
  dt_1 = dt*fudge;
  /*dt_2 = (dt/2.0)*fudge;
    dt_3 = (dt/3.0)*fudge;*/
    
  if (bInit) {
    ptr  = x;
    dt_1 = 1;
  }
  else {
    ptr  = v;
    dt_1 = fudge/dt;
  }
    
  for(i=0; (i<ns); i++) {
    s1 = s[i].shell;
    if (bInit)
      clear_rvec(x[s1]);
    switch (s[i].nnucl) {
    case 1:
      n1 = s[i].nucl1;
      for(m=0; (m<DIM); m++)
	x[s1][m]+=ptr[n1][m]*dt_1;
      break;
    case 2:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      m1 = mass[n1];
      m2 = mass[n2];
      tm = dt_1/(m1+m2);
      for(m=0; (m<DIM); m++)
	x[s1][m]+=(m1*ptr[n1][m]+m2*ptr[n2][m])*tm;
      break;
    case 3:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      n3 = s[i].nucl3;
      m1 = mass[n1];
      m2 = mass[n2];
      m3 = mass[n3];
      tm = dt_1/(m1+m2+m3);
      for(m=0; (m<DIM); m++)
	x[s1][m]+=(m1*ptr[n1][m]+m2*ptr[n2][m]+m3*ptr[n3][m])*tm;
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

static void check_pbc(FILE *fp,rvec x[],int shell)
{
  int m;
  
  for(m=0; (m<DIM); m++)
    if (fabs(x[shell][m]-x[shell-3][m]) > 0.3)
      pr_rvecs(fp,0,"SHELL-X",x+(4*(shell / 4)),4);
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
      fprintf(fp,"SHELL %5d, force %10.5f  %10.5f  %10.5f, |f| %10.5f\n",
	      shell,f[shell][XX],f[shell][YY],f[shell][ZZ],sqrt(ff2));
    check_pbc(fp,x,shell);
  }
}

int relax_shells(FILE *log,t_commrec *cr,bool bVerbose,
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
  tensor my_vir[2];
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
  /*predict_shells(log,x,v,parm->ir.delta_t,nshell,shells,md->massT,(mdstep == 0));*/
    
  /* Calculate the forces first time around */
  clear_mat(my_vir[Min]);
  do_force(log,cr,parm,nsb,my_vir[Min],mdstep,nrnb,
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
  for(count=1; (step > 0.5) &&
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
    clear_mat(my_vir[Try]);
    do_force(log,cr,parm,nsb,my_vir[Try],1,nrnb,
	     top,grps,pos[Try],v,force[Try],buf,md,ener,bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr);
    df[Try]=rms_force(force[Try],nshell,shells);

    if (debug) 
      pr_rvecs(debug,0,"F na do_force",&(force[Try][start]),homenr);

    if (debug) {
      fprintf(debug,"SHELL ITER %d\n",count);
      dump_shells(debug,pos[Try],force[Try],ftol,nshell,shells);
    }
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    Epot[Try]=ener[F_EPOT];
    if (PAR(cr)) 
      gmx_sum(NEPOT,Epot,cr);
    
    if (bVerbose && MASTER(cr))
      print_epot("",mdstep,count,step,Epot[Try],df[Try]);

    bDone=(df[Try] < ftol);
    
    if ((Epot[Try] < Epot[Min]) /* || (df[Try] < df[Min])*/) {
      Min  = Try;
      step = step0;
    }
    else
      step *= 0.8;
  }
  if (MASTER(cr) && !bDone) 
    fprintf(stderr,"EM did not converge in %d steps\n",number_steps);
  
  /* Parallelise this one! */
  memcpy(x,pos[Min],nsb->natoms*sizeof(x[0]));
  memcpy(f,force[Min],nsb->natoms*sizeof(f[0]));
  copy_mat(my_vir[Min],vir_part);
  
  return count; 
}

int relax_shells2(FILE *log,t_commrec *cr,
		  bool bVerbose,int mdstep,
		  t_parm *parm,bool bDoNS,bool bStopCM,
		  t_topology *top,real ener[],
		  rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		  rvec buf[],t_mdatoms *md,
		  t_nsborder *nsb,t_nrnb *nrnb,
		  t_graph *graph,
		  t_groups *grps,tensor vir_part,
		  int nshell,t_shell shells[],
		  t_forcerec *fr,char *traj,
		  real t,real lambda,
		  int natoms,matrix box,t_mdebin *mdebin)
{
  static bool bFirst=TRUE;
  static rvec *pos[3],*force[3];
  static real step;
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
  
  ftol         = parm->ir.em_tol;
  number_steps = parm->ir.niter;
  step0        = 1.0;   
  step         = step0;
  
  /* Do a prediction of the shell positions */
  predict_shells(log,x,v,parm->ir.delta_t,nshell,shells,md->massT,(mdstep == 0));
    
  /* Calculate the forces first time around */
  do_force(log,cr,parm,nsb,vir_part,mdstep,nrnb,
	   top,grps,x,v,force[Min],buf,md,ener,bVerbose && !PAR(cr),
	   lambda,graph,bDoNS,FALSE,fr);

  if (nshell) 
    rmsF[0]=rms_force(force[Min],nshell,shells);
  
  /* Copy x to pos[Min] & pos[Try1]: during minimization only the
   * shell positions are updated, therefore the other particles must
   * be set here.
   */
  memcpy(pos[Min], x,nsb->natoms*sizeof(x[0]));
  memcpy(pos[Try1],x,nsb->natoms*sizeof(x[0]));
  memcpy(pos[Try2],x,nsb->natoms*sizeof(x[0]));
  for(i=0; (i<nsb->natoms); i++) {
    clear_rvec(force[Try1][i]);
    clear_rvec(force[Try2][i]);
  }
  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  EEE[0] = ener[F_EPOT];
  
  if (bVerbose && MASTER(cr))
    print_epot("",mdstep,0,step,EEE[0],rmsF[Min]);
  if (debug) {
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EKIN].longname, ener[F_EKIN]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EPOT].longname, ener[F_EPOT]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_ETOT].longname, ener[F_ETOT]);
  }
      
  bDone=((nshell == 0) || (rmsF[0] < ftol));
  bMinSet=FALSE;
  
  /****************************************************** 
   *  Start the shell-relaxation loop 
   ******************************************************/
  for(count=1; 
      !(bDone || ((number_steps > 0) && (count>=number_steps))); ) {
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try1],force[Min],nshell,shells); 

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
    rmsF[Try1]=rms_force(force[Try1],nshell,shells);
#ifdef DEBUGHARD
    pr_rvecs(log,0,"F na do_force",&(force[Try1][start]),homenr);
#endif

    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    EEE[1] = ener[F_EPOT];
    
    if (bVerbose && MASTER(cr))
      print_epot("",mdstep,count,step,EEE[1],rmsF[Try1]);
    bDone=(rmsF[Try1] < ftol);
    
    /* NOW! Do line mimization *
     * Assume the energy is a quadratic function of the stepsize x:
     * a x^2 + b x + c = E
     * We try to solve the following equations:
     * a x^2 + b x + c = E[0] (x = 0)
     * a x^2 + b x + c = E[1] (x = 1)
     * 2 a x + b       = F[1] (x = 1)
     * This we can write as a matrix equation:
     *
     * ( 0  0  1 ) a   E[0]
     * ( 1  1  1 ) b = E[1]
     * ( 0  2  1 ) c   F[1]
     */
    EEE[2] = rmsF[Min]*sqrt(nshell*1.0);
    
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
    shell_pos_sd(log,step,pos[Min],pos[Try2],force[Min],nshell,shells); 
      
    /* Try the new positions */
    do_force(log,cr,parm,nsb,vir_part,1,nrnb,
	     top,grps,pos[Try2],v,force[Try2],buf,md,ener,
	     bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr);
    count++;
      
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    
    eold      = EEE[0];
    EEE[0]  = ener[F_EPOT];
    rmsF[Min] = rms_force(force[Try2],nshell,shells);
    bDone     = (rmsF[Min] < ftol);
      
    if (bVerbose && MASTER(cr)) {
      print_epot("",mdstep,count,step,EEE[0],rmsF[Min]);
      fprintf(stderr,"DE1=%10g, DE0= %10g\n",eold-EEE[1],eold-EEE[0]);
    }
    step = step0;
    Min  = Try2;
  }
  if (MASTER(cr) && !bDone) 
    fprintf(stderr,"EM did not converge in %d steps\n",number_steps);
  
  /* Parallelise this one! */
  memcpy(x,pos[Min],nsb->natoms*sizeof(x[Min]));
  memcpy(f,force[Min],nsb->natoms*sizeof(f[Min]));
#ifdef DEBUGHARD
  pr_rvecs(log,0,"X na do_relax",&(x[start]),homenr);
  pr_rvecs(log,0,"F na do_relax",&(f[start]),homenr);
#endif

  return count; 
}

