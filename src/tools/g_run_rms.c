/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_g_run_rms_c = "$Id$";

#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pbc.h"
#include "macros.h"
#include "gstat.h"
#include "xvgr.h"
#include "tpxio.h"

real calc_ener(int natoms,real w_rms[],rvec x[],rvec xp[])
{
  real e,tmas;
  int  j,m;
  
  /*calculate energy */
  e=0;
  tmas=0;
  for(j=0;(j<natoms);j++) {
    tmas+=w_rms[j];
    for(m=0;(m<3);m++) 
      e+=w_rms[j]*(x[j][m]-xp[j][m])*(x[j][m]-xp[j][m]);
  }
  /*return energy*/
  return(sqrt(e/tmas));
}

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_run_rms computes the root mean square deviation (RMSD) of a structure",
    "from a trajectory x(t), with respect to a reference structure from the",
    "same trajectory, but at a specified time before the current time",
    "x(t-dt) by LSQ fitting the structures on top of each other[PAR]",
    "This tells you something about the mobility of structure as a function",
    "of time"
  };
  static int run_time = 5;
  t_pargs pa[] = {
    { "-t", FALSE, etINT, &run_time,
      "Time interval dt between reference and actual structure" }
  };
  int          step,nre,natom,i,j,m,teller=0;
  real         t,lambda,*w_rls,*w_rms,tmas;
  int          status;
  t_tpxheader  header;
  t_inputrec   ir;
  t_topology   top;
  matrix       box;
  rvec         **x,*xp,*v,xcm,*temp;
  FILE         *fp;
  
#define RLS 0
#define RMS 1
  int          isize[2];
  atom_id      *index[2];
  char         *grpnames[2];
  t_filenm     fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffREAD },
    { efXVG, NULL, "runrms", ffWRITE }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  snew(x,run_time+1);

  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&header);
  snew(xp,header.natoms);
  for(i=0;(i<run_time+1);i++)
    snew(x[i],header.natoms);
  snew(w_rls,header.natoms);
  snew(w_rms,header.natoms);
  snew(v,header.natoms);
  
  read_tpx(ftp2fn(efTPX,NFILE,fnm),
	   &step,&t,&lambda,&ir,
	   box,&natom,xp,NULL,NULL,&top);

  /*set box type*/
  init_pbc(box,FALSE);
 
  fprintf(stderr,"select group for root least squares fit\n");
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize[RLS],&index[RLS],&grpnames[RLS]);
  for(i=0;(i<isize[RLS]);i++)
    w_rls[index[RLS][i]]=top.atoms.atom[index[RLS][i]].m;
  
  fprintf(stderr,"select group for root mean square calculation\n");
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize[RMS],&index[RMS],&grpnames[RMS]);
  for(i=0;(i<isize[RMS]);i++)
    w_rms[index[RMS][i]]=top.atoms.atom[index[RMS][i]].m;
    
  fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Running RMSD","Time (ps)","RMSD (nm)");
  
  rm_pbc(&(top.idef),top.atoms.nr,box,xp,xp);
  
  /*calculate mass*/
  tmas=0;
  for(i=0;(i<header.natoms);i++)
    tmas+=w_rls[i];
  
  /*set center of mass to zero for reference coordinates*/
  
  /*do a first step*/
  natom=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x[0],box);
  for(i=0; (i<run_time); i++) {
    read_next_x(status,&t,natom,x[i],box);
    rm_pbc(&(top.idef),top.atoms.nr,box,x[i],x[i]);
    for(m=0;(m<3);m++)
      xcm[m]=0;
    for(j=0;(j<header.natoms);j++)
      for(m=0;(m<3);m++) {
	xcm[m]+=x[i][j][m]*w_rls[j];
      }
    for(m=0;(m<3);m++)
      xcm[m]/=tmas;
    for(j=0;(j<header.natoms);j++) {
      for(m=0;(m<3);m++)
	x[i][j][m]-=xcm[m];
    }
  }
 
  read_next_x(status,&t,natom,x[run_time],box);  
  do {
    if ((teller++ % 10) == 0)
      fprintf(stderr,"\rFrame: %d",teller);
    rm_pbc(&(top.idef),top.atoms.nr,box,x[run_time],x[run_time]);
    for(m=0;(m<3);m++)
      xcm[m]=0;
    for(j=0;(j<header.natoms);j++)
      for(m=0;(m<3);m++) {
	xcm[m]+=x[run_time][j][m]*w_rls[j];
      }
    for(m=0;(m<3);m++)
      xcm[m]/=tmas;
    for(j=0;(j<header.natoms);j++) {
      for(m=0;(m<3);m++)
	x[run_time][j][m]-=xcm[m];
    }
    
    /*calculate energy of root_least_squares*/
    do_fit(header.natoms,w_rls,x[0],x[run_time]);
    fprintf(fp,"%8.4f %8.4f\n",
	    t,calc_ener(header.natoms,w_rms,x[0],x[run_time]));
    
    /*swap coordinates*/
    temp=x[0];
    for(i=0;(i<run_time);i++)
      x[i]=x[i+1];
    x[run_time]=temp;
  } while ((read_next_x(status,&t,natom,x[run_time],box)));
  fprintf(stderr,"\n");
  fclose(fp);

  close_trj(status);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
  
  thanx(stdout);
  
  return 0;
}
