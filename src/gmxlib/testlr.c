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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "vec.h"
#include "physics.h"
#include "macros.h"
#include "names.h"
#include "smalloc.h"
#include "tpxio.h"
#include "statutil.h"
#include "writeps.h"
#include "copyrite.h"
#include "xvgr.h"
#include "minvert.h"
#include "pppm.h"
#include "readinp.h"
#include "main.h"
#include "force.h"
#include "nrnb.h"
#include "shift_util.h"
#include "ewald_util.h"
#include "mshift.h"
#include "poisson.h"
#include "mdatoms.h"

void pr_f(char *fn,int natoms,rvec f[])
{
  FILE *fp;
  int  i;
  
  fp=ffopen(fn,"w");
  for(i=0; (i<natoms); i++)
    fprintf(fp,"  %12.5e\n  %12.5e\n  %12.5e\n",f[i][XX],f[i][YY],f[i][ZZ]);
  fclose(fp);
}

void test_pppm(FILE *log,       bool bVerbose,
	       bool bGenerGhat, char *ghatfn,
	       t_atoms *atoms,  t_inputrec *ir,
	       rvec x[],        rvec f[],
	       real charge[],   rvec box,
	       real phi[],      real phi_s[],
	       int nmol,        t_commrec *cr,
	       bool bOld,       t_block *cgs)
{
  char       buf[256];
  real       ener;
  int        i;
  t_nrnb     nrnb;
  t_nsborder nsb;
  
  init_nrnb(&nrnb);
  calc_nsb(cgs,1,&nsb,0);
  
  /* First time only setup is done! */
  init_pppm(log,cr,&nsb,bVerbose,bOld,box,ghatfn,ir);
  
  ener = do_pppm(log,bVerbose,x,f,charge,box,phi,cr,&nsb,&nrnb);
  fprintf(log,"Vpppm = %g\n",ener);
  
  sprintf(buf,"PPPM-%d.pdb",ir->nkx);
  write_pqr(buf,atoms,x,phi,0);
  
  pr_f("pppm-force",atoms->nr,f);
  
  calc_ener(log,buf,FALSE,nmol,atoms->nr,phi,charge,&atoms->excl);
  
  for(i=0; (i<atoms->nr); i++) 
    phi[i]+=phi_s[i];
  sprintf(buf,"PPPM-%d+SR",ir->nkx);
  calc_ener(log,buf,FALSE,nmol,atoms->nr,phi,charge,&atoms->excl);
  strcat(buf,".pdb");
  write_pqr(buf,atoms,x,phi,0);
}

void test_poisson(FILE *log,       bool bVerbose,
		  t_atoms *atoms,  t_inputrec *ir,
		  rvec x[],        rvec f[],
		  real charge[],   rvec box,
		  real phi[],      real phi_s[],
		  int nmol,        t_commrec *cr,
		  bool bFour,      rvec f_four[],
		  real phi_f[],    bool bOld)
{
  char buf[256];
  real ener;
  rvec beta;
  int  i,nit;
  t_nrnb nrnb;
  
  init_nrnb(&nrnb);
  
  /* First time only setup is done! */
  if (bFour) {
    for(i=0; (i<atoms->nr); i++)
      phi_f[i] -= phi_s[i];
    ener = do_optimize_poisson(log,bVerbose,ir,atoms->nr,x,f,charge,
			       box,phi,cr,&nrnb,f_four,phi_f,beta,bOld);
    for(i=0; (i<atoms->nr); i++)
      phi_f[i] += phi_s[i];
    nit = 0;
  }
  else {
    ener = do_poisson(log,bVerbose,ir,atoms->nr,x,f,charge,box,phi,
		      cr,&nrnb,&nit,bOld);
  }
    
  fprintf(log,"Vpoisson = %g, nit = %d\n",ener,nit);
  
  sprintf(buf,"POISSON-%d.pdb",ir->nkx);
  write_pqr(buf,atoms,x,phi,0);
  
  pr_f("poisson-force",atoms->nr,f);
  
  calc_ener(log,buf,FALSE,nmol,atoms->nr,phi,charge,&atoms->excl);
  
  for(i=0; (i<atoms->nr); i++) 
    phi[i]+=phi_s[i];
  sprintf(buf,"POISSON-%d+SR",ir->nkx);
  calc_ener(log,buf,FALSE,nmol,atoms->nr,phi,charge,&atoms->excl);
  strcat(buf,".pdb");
  write_pqr(buf,atoms,x,phi,0);
}

void test_four(FILE *log,int NFILE,t_filenm fnm[],t_atoms *atoms,
	       t_inputrec *ir,rvec x[],rvec f[],rvec box,real charge[],
	       real phi_f[],real phi_s[],int nmol,t_commrec *cr,
	       bool bOld,bool bOldEwald)
{
  int  i;
  real energy;

  if (bOldEwald)  
    energy = do_ewald(log,ir,atoms->nr,x,f,charge,box,phi_f,cr,bOld);
  else
    energy = do_ewald_new(log,ir,atoms->nr,x,f,charge,box,phi_f,cr,bOld);
  
  /*symmetrize_phi(log,atoms->nr,phi_f,bVerbose);*/
    
  /* Plot the long range fourier solution in a matrix */    
  plot_phi(opt2fn("-elf",NFILE,fnm),box,atoms->nr,x,phi_f);
  print_phi(opt2fn("-of",NFILE,fnm),atoms->nr,x,phi_f);
  calc_ener(log,"Fourier",FALSE,nmol,atoms->nr,phi_f,charge,&atoms->excl);
  write_pqr(opt2fn("-pf",NFILE,fnm),atoms,x,phi_f,0.0/*1.5*box[XX]*/);
  pr_f("four-force",atoms->nr,f);
  
  /* Calc and plot the total potential */
  for(i=0; (i<atoms->nr); i++) {
    phi_f[i]+=phi_s[i];
    /*clear_rvec(f[i]);*/
  }
  calc_ener(log,"Fourier+SR",FALSE,nmol,atoms->nr,phi_f,charge,&atoms->excl);
}

static void print_opts(FILE *fp,t_inputrec *ir,bool bFour)
{
  fprintf(fp,"Ewald solution: %s\n",bool_names[bFour]);
  fprintf(fp,"r1:       %10.3e\n",ir->rcoulomb_switch);
  fprintf(fp,"rc:       %10.3e\n",ir->rcoulomb);
  if (bFour)
    fprintf(fp,"KVectors: %10d  %10d  %10d\n",ir->nkx,ir->nky,ir->nkz);
  fprintf(fp,"\n");
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "testlr tests the PPPM and Ewald method for the",
    "long range electrostatics problem."
  };
  static t_filenm  fnm[] = {
    { efTPX, NULL,   NULL,       ffREAD },
    { efHAT, "-g",   "ghat",     ffOPTRD },
    { efOUT, "-o",   "rho",      ffOPTWR },
    { efOUT, "-op",  "lr-pb",    ffOPTWR },
    { efOUT, "-of",  "lr-four",  ffOPTWR },
    { efOUT, "-opt", "tot-pb",   ffOPTWR },
    { efOUT, "-oft", "tot-four", ffOPTWR },
    { efOUT, "-fin", "lr-four",  ffOPTWR },
    { efEPS, "-es",  "sr",       ffOPTWR },
    { efEPS, "-elf", "lr-four",  ffOPTWR },
    { efEPS, "-etf", "tot-four", ffOPTWR },
    { efEPS, "-qr",  "qk-real",  ffOPTWR },
    { efEPS, "-qi",  "qk-im",    ffOPTWR },
    { efEPS, "-elp", "lr-pb",    ffOPTWR },
    { efEPS, "-etp", "tot-pb",   ffOPTWR },
    { efEPS, "-rho", "rho",      ffOPTWR },
    { efEPS, "-qq",  "charge",   ffOPTWR },
    { efXVG, "-gt",  "gk-tab",   ffOPTWR },
    { efXVG, "-fcorr","fcorr",   ffWRITE },
    { efXVG, "-pcorr","pcorr",   ffWRITE },
    { efXVG, "-ftotcorr","ftotcorr",   ffWRITE },
    { efXVG, "-ptotcorr","ptotcorr",   ffWRITE },
    { efLOG, "-l",   "fptest",   ffWRITE },
    { efXVG, "-gr",  "spread",   ffOPTWR },
    { efPDB, "-pf",  "pqr-four", ffOPTWR },
    { efPDB, "-phitot", "pppm-phitot", ffOPTWR }
  };
#define NFILE asize(fnm)
  FILE         *log;
  t_topology   top;
  t_tpxheader  stath;
  t_inputrec   ir;
  t_block      *excl;
  t_forcerec   *fr;
  t_commrec    *cr;
  t_mdatoms    *mdatoms;
  t_graph      *graph;
  int          i,step,nre,natoms,nmol;
  rvec         *x,*f_sr,*f_excl,*f_four,*f_pppm,*f_pois,box_size,hbox;
  matrix       box;
  real         t,lambda,vsr,*charge,*phi_f,*phi_pois,*phi_s,*phi_p3m,*rho;
  
  static bool bFour=FALSE,bVerbose=FALSE,bGGhat=FALSE,bPPPM=TRUE,
    bPoisson=FALSE,bOld=FALSE,bOldEwald=TRUE;
  static int nprocs = 1;
  static t_pargs pa[] = {
    { "-np",     FALSE, etINT,  &nprocs,  "Do it in parallel" },
    { "-ewald",  FALSE, etBOOL, &bFour,   "Do an Ewald solution"},
    { "-pppm",   FALSE, etBOOL, &bPPPM,   "Do a PPPM solution" },
    { "-poisson",FALSE, etBOOL, &bPoisson,"Do a Poisson solution" },
    {    "-v",   FALSE, etBOOL, &bVerbose,"Verbose on"},
    { "-ghat",   FALSE, etBOOL, &bGGhat,  "Generate Ghat function"},
    { "-old",    FALSE, etBOOL, &bOld,    "Use old function types"},
    { "-oldewald",FALSE,etBOOL, &bOldEwald,"Use old Ewald code"}
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  if (nprocs > 1) {
    cr = init_par(&argc,argv);
    open_log(ftp2fn(efLOG,NFILE,fnm),cr);
    log = stdlog;
  }
  else {
    cr     = init_par(&argc,argv);
    log    = ftp2FILE(efLOG,NFILE,fnm,"w");
    stdlog = log;  }
  

  /* Read topology and coordinates */
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&stath,FALSE);
  snew(x,stath.natoms);
  snew(f_sr,stath.natoms);
  snew(f_excl,stath.natoms);
  snew(f_four,stath.natoms);
  snew(f_pppm,stath.natoms);
  snew(f_pois,stath.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	   box,&natoms,x,NULL,NULL,&top);
  excl=&(top.atoms.excl);
  nmol=top.blocks[ebMOLS].nr;

  /* Allocate space for potential, charges and rho (charge density) */
  snew(charge,stath.natoms);
  snew(phi_f,stath.natoms);
  snew(phi_p3m,stath.natoms);
  snew(phi_pois,stath.natoms);
  snew(phi_s,stath.natoms);
  snew(rho,stath.natoms);
  
  /* Set the charges */
  for(i=0; (i<natoms); i++)
    charge[i]=top.atoms.atom[i].q;

  /* Make a simple box vector instead of tensor */
  for(i=0; (i<DIM); i++) 
    box_size[i]=box[i][i];
  
  /* Set some constants */
  fr      = mk_forcerec();
  mdatoms = atoms2md(&(top.atoms),FALSE,FALSE);
  
  set_LRconsts(log,ir.rcoulomb_switch,ir.rcoulomb,box_size,fr);
  init_forcerec(log,fr,&ir,&(top.blocks[ebMOLS]),cr,
		&(top.blocks[ebCGS]),&(top.idef),mdatoms,box,FALSE);
  calc_shifts(box,box_size,fr->shift_vec,FALSE);

  /* Periodicity stuff */  
  graph = mk_graph(&(top.idef),top.atoms.nr,FALSE,FALSE);
  shift_self(graph,fr->shift_vec,x);

  calc_LRcorrections(log,0,natoms,ir.rcoulomb_switch,
		     ir.rcoulomb,charge,excl,x,f_excl,bOld);
  pr_f("f_excl.dat",natoms,f_excl);
  
  /* Compute the short range potential */
  put_atoms_in_box(natoms,box,x);
  vsr=phi_sr(log,natoms,x,charge,ir.rcoulomb,
	     ir.rcoulomb_switch,box_size,phi_s,excl,f_sr,bOld); 
  pr_f("f_sr.dat",natoms,f_sr);
  
  /* Plot the short range potential in a matrix */    
  calc_ener(log,"Short Range",TRUE,nmol,natoms,phi_s,charge,excl);
  
  
  if (bFour)   
    test_four(log,NFILE,fnm,&(top.atoms),&ir,x,f_four,box_size,charge,phi_f,
	      phi_s,nmol,cr,bOld,bOldEwald);
  
  if (bPPPM) 
    test_pppm(log,bVerbose,bGGhat,opt2fn("-g",NFILE,fnm),
	      &(top.atoms),&ir,x,f_pppm,charge,box_size,phi_p3m,phi_s,nmol,
	      cr,bOld,&(top.blocks[ebCGS]));
  
  if (bPoisson)
    test_poisson(log,bVerbose,
		 &(top.atoms),&ir,x,f_pois,charge,box_size,phi_pois,
		 phi_s,nmol,cr,bFour,f_four,phi_f,bOld);
	        
  if (bPPPM && bFour) 
    analyse_diff(log,"PPPM",
		 top.atoms.nr,f_four,f_pppm,phi_f,phi_p3m,phi_s,
		 opt2fn("-fcorr",NFILE,fnm),
		 opt2fn("-pcorr",NFILE,fnm),
		 opt2fn("-ftotcorr",NFILE,fnm),
		 opt2fn("-ptotcorr",NFILE,fnm));
  
  if (bPoisson && bFour) 
    analyse_diff(log,"Poisson",
		 top.atoms.nr,f_four,f_pois,phi_f,phi_pois,phi_s,
		 opt2fn("-fcorr",NFILE,fnm),
		 opt2fn("-pcorr",NFILE,fnm),
		 opt2fn("-ftotcorr",NFILE,fnm),
		 opt2fn("-ptotcorr",NFILE,fnm));
  
  fclose(log);
  
  thanx(stderr);
  
  return 0;
}

