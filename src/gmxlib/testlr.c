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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_testlr_c = "$Id$";

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
#include "assert.h"
#include "copyrite.h"
#include "xvgr.h"
#include "minvert.h"
#include "pppm.h"
#include "readinp.h"
#include "main.h"
#include "force.h"
#include "nrnb.h"
#include "lrutil.h"

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
	       int nmol,        t_commrec *cr)
{
  char buf[256];
  real ener;
  int  i;
  t_nrnb nrnb;
  
  init_nrnb(&nrnb);
  
  /* First time only setup is done! */
  ener = do_pppm(log,bVerbose,bGenerGhat,ghatfn,ir,atoms->nr,  
		 x,f,charge,box,phi,cr,&nrnb);
  ener = do_pppm(log,bVerbose,bGenerGhat,NULL,NULL,atoms->nr,  
		 x,f,charge,box,phi,cr,&nrnb);
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

void test_four(FILE *log,int NFILE,t_filenm fnm[],t_atoms *atoms,
	       t_inputrec *ir,rvec x[],rvec f[],rvec box,real charge[],
	       real phi_f[],real phi_s[],int nmol,t_commrec *cr)
{
  int  i;
  real energy;
  
  energy = do_ewald(log,ir,atoms->nr,x,f,charge,box,phi_f,cr);
  
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

void analyse_diff(FILE *log,int natom,rvec ffour[],rvec fpppm[],
		  real phi_f[],real phi_p[],char *fcorr,char *pcorr)
{
  int  i,m;
  FILE *fp;
  char buf[120];
  real f2sum=0,p2sum=0;
  real df,fmax,dp,pmax;
  
  fmax = fabs(ffour[0][0]-fpppm[0][0]);
  pmax = fabs(phi_f[0] - phi_p[0]);
  
  for(i=0; (i<natom); i++) {
    dp     = fabs(phi_f[i] - phi_p[i]);
    pmax   = max(dp,pmax);
    p2sum += dp*dp;
    for(m=0; (m<DIM); m++) {
      df     = fabs(ffour[i][m] - fpppm[i][m]);
      fmax   = max(df,fmax);
      f2sum += df*df;
    }
  }
  
  fprintf(log,"\n********************************\nERROR ANALYSIS\n");
  fprintf(log,"%-10s%12s%12s\n","Error:","Max Abs","RMS");
  fprintf(log,"%-10s  %10.3f  %10.3f\n",
	  "Force",fmax,sqrt(f2sum/(3.0*natom)));
  fprintf(log,"%-10s  %10.3f  %10.3f\n",
	  "Potential",pmax,sqrt(p2sum/(natom)));
  
  fp = xvgropen(fcorr,"Force Correlation","Four-Force","PPPM-Force");
  for(i=0; (i<natom); i++) {
    for(m=0; (m<DIM); m++)
      fprintf(fp,"%10.3f  %10.3f\n",ffour[i][m],fpppm[i][m]);
  }
  ffclose(fp);
  xvgr_file(fcorr,NULL);
  
  fp = xvgropen(pcorr,"Potential Correlation","Four-Pot","PPPM-Pot");
  for(i=0; (i<natom); i++) {
    fprintf(fp,"%10.3f  %10.3f\n",phi_f[i],phi_p[i]);
  }
  ffclose(fp);
  xvgr_file(pcorr,NULL);
}

static void print_opts(FILE *fp,t_inputrec *ir,bool bFour)
{
  fprintf(fp,"Ewald solution: %s\n",bool_names[bFour]);
  fprintf(fp,"r1:       %10.3e\n",ir->rshort);
  fprintf(fp,"rc:       %10.3e\n",ir->rlong);
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
  int          i,step,nre,natoms,nmol;
  rvec         *x,*f_four,*f_pppm,box_size,hbox;
  matrix       box;
  real         t,lambda,vsr,*charge,*phi_f,*phi_s,*phi_p3m,*rho;
  
  static bool bFour=FALSE,bVerbose=FALSE,bGGhat=FALSE,bPPPM=TRUE;
  static int nprocs = 1;
  static t_pargs pa[] = {
    { "-np",     FALSE, etINT,  &nprocs,  "Do it in parallel" },
    { "-ewald",  FALSE, etBOOL, &bFour,   "Do an Ewald solution"},
    { "-pppm",   FALSE, etBOOL, &bPPPM,   "Do a PPPM solution" },
    {    "-v",   FALSE, etBOOL, &bVerbose,"Verbose on"},
    { "-ghat",   FALSE, etBOOL, &bGGhat,  "Generate Ghat function"}
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  if (nprocs > 1) {
    cr = init_par(argv);
    open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  }
  else {
    cr = NULL;
    log=ftp2FILE(efLOG,NFILE,fnm,"w");
  }
  
#ifdef CINVSQRT
  init_lookup_table(log);
#endif

  /* Read topology and coordinates */
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&stath);
  snew(x,stath.natoms);
  snew(f_four,stath.natoms);
  snew(f_pppm,stath.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	   box,&natoms,x,NULL,NULL,&top);
  excl=&(top.atoms.excl);
  nmol=top.blocks[ebMOLS].nr;
  
  /* Allocate space for potential, charges and rho (charge density) */
  snew(charge,stath.natoms);
  snew(phi_f,stath.natoms);
  snew(phi_p3m,stath.natoms);
  snew(phi_s,stath.natoms);
  snew(rho,stath.natoms);
  
  /* Set the charges */
  for(i=0; (i<natoms); i++)
    charge[i]=top.atoms.atom[i].q;

  /* Make a simple box vector instead of tensor */
  for(i=0; (i<DIM); i++) 
    box_size[i]=box[i][i];
  
  /* Set some constants */
  fr = mk_forcerec();
  set_LRconsts(log,ir.rshort,ir.rlong,box_size,fr);

  /* Compute the short range potential */
  vsr=phi_sr(natoms,x,charge,ir.rlong,ir.rshort,box_size,phi_s,excl);  
  
  /* Plot the short range potential in a matrix */    
  calc_ener(log,"Short Range",TRUE,nmol,natoms,phi_s,charge,excl);

  if (bFour)   
    test_four(log,NFILE,fnm,&(top.atoms),&ir,x,f_four,box_size,charge,phi_f,
	      phi_s,nmol,cr);
  
  if (bPPPM) 
    test_pppm(log,bVerbose,bGGhat,opt2fn("-g",NFILE,fnm),
	      &(top.atoms),&ir,x,f_pppm,charge,box_size,phi_p3m,phi_s,nmol,cr);
  
  if (bPPPM && bFour) 
    analyse_diff(log,top.atoms.nr,f_four,f_pppm,phi_f,phi_p3m,
		 opt2fn("-fcorr",NFILE,fnm),
		 opt2fn("-pcorr",NFILE,fnm));
  
  fclose(log);
  
  thanx(stdout);
  
  return 0;
}

