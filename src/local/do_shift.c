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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_do_shift_c = "$Id$";

#include "sysstuff.h"
#include "typedefs.h"
#include "string2.h"
#include "strdb.h"
#include "macros.h"
#include "smalloc.h"
#include "mshift.h"
#include "statutil.h"
#include "copyrite.h"
#include "confio.h"
#include "fatal.h"
#include "xvgr.h"
#include "gstat.h"
#include "rdgroup.h"

void cat(FILE *out,char *fn,real t)
{
  FILE *in;
  char *ptr,buf[256];
  int    anr,rnr;
  char   anm[24],rnm[24];
  double f1,f2,f3,f4,f5,f6;
  
  in=ffopen(fn,"r");
  while ((ptr=fgets2(buf,255,in)) != NULL) {
    sscanf(buf,"%d%d%s%s%lf%lf%lf%lf%lf%lf",
	   &anr,&rnr,rnm,anm,&f1,&f2,&f3,&f4,&f5,&f6);
    fprintf(out,"%8g  %10g  %10g  %10g  %10g  %10g  %10g  %s%d-%s%d\n",
	    t,f6,f1,f2,f3,f4,f5,rnm,rnr,anm,anr);
  }
  if ((int)strlen(buf) > 0) 
    fprintf(out,"%s\n",buf);
  fclose(in);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "do_shift reads a trajectory file and computes the chemical",
    "shift for each time frame (or every [BB]dt[bb] ps) by",
    "calling the 'total' program. If you do not have the total program,",
    "get it. do_shift assumes that the total executable is in",
    "/home/mdgroup/total/total. If that is not the case, then you should",
    "set an environment variable [BB]TOTAL[bb] as in: [PAR]",
    "[TT]setenv TOTAL /usr/local/bin/total[tt][PAR]",
    "where the right hand side should point to the total executable.[PAR]",
    "Output is printed in files shift.out where t is the time of the frame.[PAR]",
    "The program also needs an input file called [BB]random.dat[bb] which",
    "contains the random coil chemical shifts of all protons."
  };
  static real dt=0.0;
  t_pargs pa[] = {
    { "-dt", FALSE, etREAL, &dt, "Time interval between frames." }
  };
  static char *bugs[] = {
    "The program is very slow"
  };
  static     char *OXYGEN="O";
  FILE       *out,*tot;
  t_topology *top;
  t_atoms    *atoms;
  int        status,nres;
  real       t,nt;
  int        i,natoms,nframe=0;
  matrix     box;
  int        gnx;
  char       *grpnm,*randf;
  atom_id    *index;
  rvec       *x,*x_s;
  char       pdbfile[L_tmpnam],tmpfile[L_tmpnam];
  char       total[256],*dptr;
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD },
    { efTPB, NULL,   NULL,     ffREAD },
    { efNDX, NULL,   NULL,     ffREAD },
    { efOUT, "-o",   "shift",  ffWRITE },
    { efDAT, "-d",   "random", ffREAD }
  };
  char *leg[] = { "shift","ring","anisCO","anisCN","sigmaE","sum" };
#define NFILE asize(fnm)

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,
		    asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
		    
  top=read_top(ftp2fn(efTPB,NFILE,fnm));
  atoms=&(top->atoms);
  nres=atoms->nres;
  for(i=0; (i<atoms->nr); i++)
    if ((strcmp(*atoms->atomname[i],"O1") == 0) ||
	(strcmp(*atoms->atomname[i],"O2") == 0) ||
	(strcmp(*atoms->atomname[i],"OXT") == 0) ||
	(strcmp(*atoms->atomname[i],"OT") == 0))
      atoms->atomname[i]=&OXYGEN;
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);
  
  snew(x_s,atoms->nr);

  (void) tmpnam(pdbfile);
  (void) tmpnam(tmpfile);
  if ((dptr=getenv("TOTAL")) == NULL)
    dptr="/home/mdgroup/total/total";
  sprintf(total,"%s > /dev/null",dptr);
  fprintf(stderr,"total cmd='%s'\n",total);
  randf=ftp2fn(efDAT,NFILE,fnm);
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if (natoms != atoms->nr) 
    fatal_error(0,"Trajectory does not match topology!");
  out=ftp2FILE(efOUT,NFILE,fnm,"w");
  xvgr_legend(out,asize(leg),leg);
  nt=t;
  do {
    if (t >= nt) {
      fprintf(stderr,"\rt=%.1f",t);
      rm_pbc(&(top->idef),top->atoms.nr,box,x,x_s);
      write_pdb_conf_indexed(pdbfile,atoms,x_s,box,gnx,index);
      
      tot=popen(total,"w");
      fprintf(tot,"%s\n%s\n3\nn\n%s\nn\n",pdbfile,tmpfile,randf);
      pclose(tot);
      cat(out,tmpfile,t);
      
      remove(pdbfile);
      remove(tmpfile);
      nt+=dt;
      nframe++;
    }
  } while(read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
