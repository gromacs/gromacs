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
 * Great Red Owns Many ACres of Sand 
 */
#include "sysstuff.h"
#include "typedefs.h"
#include "futil.h"
#include "mvdata.h"
#include "rwtop.h"
#include "network.h"
#include "fatal.h"
#include "smalloc.h"
#include "nsb.h"
#include "statutil.h"
#include "statusio.h"
#include "main.h"
#include "parm.h"
#include "txtdump.h"
#include "macros.h"
#include "tconf.h"

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bMaster,bool bPar,bool bVerbose,int stepout,
	     int start,int homenr,
	     t_inputrec *ir,t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],t_energy stat_e[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],real mass[],
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,
	     t_energy fvir[9],t_energy svir[9],
	     real tmass)
{
  return 0;
}

static void printit(FILE *log,t_topology *top,t_nsborder *nsb,
		     t_parm *parm,rvec x[],rvec v[],t_commrec *cr)
{
  int indent=0;
  
  indent=pr_title(log,indent,"Test Top");
  pr_inputrec(log,indent,"ir",&(parm->ir));
  pr_rvecs(log,indent,"box",parm->box,DIM);
  pr_rvecs(log,indent,"vir",parm->vir,DIM);
  pr_rvecs(log,indent,"pres",parm->pres,DIM);
  pr_rvecs(log,indent,"x",x,nsb->natoms);
  pr_rvecs(log,indent,"v",v,nsb->natoms);
  pr_rvecs(log,indent,"f",NULL,nsb->natoms);
  pr_energies(log,indent,"e",NULL,0);
  pr_top(log,indent,"topology",top);
  
  print_nsb(log,nsb);
  fprintf(log,"Been there, Done it\n");
}

static void moveit(t_commrec *cr,char *tpbfile)
{
  FILE         *fp;
  t_statheader sh;
  t_parm       *parm;
  t_nsborder   nsb;
  t_topology   top;
  t_grps       grps;
  
  char         *szVer;
  int          natoms,nre,step;
  real         t,lambda;
  rvec         *x,*v;
  real         *mass;

  fp=ffopen(tpbfile,"r");
  printf("opened %s\n",tpbfile);
  szVer=rd_header(fp,&sh);
  printf("status version:\"%s\" done\n",szVer); fflush(stdout);
  snew(parm,1);
  snew(x,sh.natoms);
  snew(v,sh.natoms);
  snew(mass,sh.natoms);
  szVer=rd_hstatus(fp,&sh,&step,&t,&lambda,&parm->ir,parm->box,
                   NULL,NULL,&natoms,x,v,NULL,&nre,NULL,&top);
  printf("status version:\"%s\" done\n",szVer); fflush(stdout);
  calc_nsb(&(top.blocks[ebCGS]),cr->nprocs,&nsb);
  
  mv_data(cr->left,cr->right,parm,&nsb,&top,x,v);

  init_groups(stdlog,&(top.atoms),&(parm->ir.opts),&grps);
  pr_groups(stdlog,&(top.atoms),&grps,&(parm->ir.opts),10,FALSE);

  printit(stdlog,&top,&nsb,parm,x,v,cr);
}

static void loadit(t_commrec *cr)
{
  t_parm     parm;
  t_topology top;
  t_nsborder nsb;
  rvec       *x;
  rvec       *v;

  ld_data(cr->left,cr->right,&parm,&nsb,&top,&x,&v);
  
  printit(stdlog,&top,&nsb,&parm,x,v,cr);
}

int main(int argc,char *argv[])
{
  t_commrec *cr;
  t_filenm fnm[] = {
    { efTPB, NULL, NULL,      FALSE }
  };
#define NFILE asize(fnm)

  cr=init_par(&argc,argv);
  if (cr->nprocs < 2) 
    fatal_error(0,"Not enough processors");
    
  fprintf(stdlog,"I am CPU %d\n",cr->pid);

  if (cr->pid == 0) {
    CopyRight(stdout,argv[0]);
    parse_common_args(&argc,argv,0,NFILE,fnm,TRUE,NULL);
    moveit(cr,ftp2fn(efTPB,NFILE,fnm));
  }
  else if (cr->pid == 1) {
    loadit(cr);
  }
  return 0;
}
