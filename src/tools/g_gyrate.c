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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_g_gyrate_c = "$Id$";

#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "txtdump.h"

real calc_gyro(rvec x[],int gnx,atom_id index[],t_atom atom[],real tm,
	       rvec gvec,rvec d,bool bQ,bool bRot)
{
  int    i,ii,m;
  real   gyro,dx2,m0;
  matrix trans;
  rvec   comp;

  if (bRot) {
    principal_comp(gnx,index,atom,x,trans,d);
    for(m=0; (m<DIM); m++)
      d[m]=sqrt(d[m]/tm);
#ifdef DEBUG
    pr_rvecs(stderr,0,"trans",trans,DIM);
#endif
    /* rotate_atoms(gnx,index,x,trans); */
  }
  clear_rvec(comp);
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    if (bQ)
      m0=fabs(atom[ii].q);
    else
      m0=atom[ii].m;
    for(m=0; (m<DIM); m++) {
      dx2=x[ii][m]*x[ii][m];
      comp[m]+=dx2*m0;
    }
  }
  gyro=comp[XX]+comp[YY]+comp[ZZ];
  
  for(m=0; (m<DIM); m++)
    gvec[m]=sqrt((gyro-comp[m])/tm);
    
  return sqrt(gyro/tm);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_gyrate computes the radius of gyration of a group of atoms",
    "and the radii of gyration about the x, y and z axes,"
    "as a function of time. The atoms are explicitly mass weighted."
  };
  static bool bQ=FALSE,bRot=FALSE;
  t_pargs pa[] = {
    { "-q", FALSE, etBOOL, &bQ,
      "Use absolute value of the charge of an atom as weighting factor instead of mass" },
    { "-p", FALSE, etBOOL, &bRot,
      "Calculate the radii of gyration about the principal axes." }
  };
  FILE       *out;
  int        status;
  t_topology *top;
  rvec       *x,*x_s;
  rvec       xcm,gvec;
  rvec       d;         /* eigenvalues of inertia tensor */
  matrix     box;
  real       t,tm,gyro;
  int        natoms;
  char       *grpname;
  int        i,j,gnx;
  atom_id    *index;
  char       *leg[] = { "Rg", "RgX", "RgY", "RgZ" }; 
#define NLEG asize(leg) 
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL, ffREAD }, 
    { efTPX, NULL, NULL, ffREAD }, 
    { efXVG, NULL, NULL, ffWRITE }, 
    { efNDX, NULL, NULL, ffOPTRD } 
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 
  clear_rvec(d); 

  for(i=1; (i<argc); i++) { 
    if (strcmp(argv[i],"-q") == 0) { 
      bQ=TRUE; 
      fprintf(stderr,"Will print radius normalised by charge\n"); 
    } 
    else if (strcmp(argv[i],"-r") == 0) { 
      bRot=TRUE; 
      fprintf(stderr,"Will rotate system along principal axes\n"); 
    }
  } 
  top=read_top(ftp2fn(efTPX,NFILE,fnm)); 
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 
  snew(x_s,natoms); 

  j=0; 
  if (bQ) 
    out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
		 "Radius of Charge","Time (ps)","Rg (nm)"); 
  else 
    out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
		 "Radius of gyration","Time (ps)","Rg (nm)"); 
  if (bRot) 
    fprintf(out,"@ subtitle \"Axes are principal component axes\"\n");
  xvgr_legend(out,NLEG,leg);
  do {
    rm_pbc(&(top->idef),top->atoms.nr,box,x,x_s);
    tm=sub_xcm(x_s,gnx,index,top->atoms.atom,xcm,bQ);
    gyro=calc_gyro(x_s,gnx,index,top->atoms.atom,tm,gvec,d,bQ,bRot);    

    if (bRot) {
      fprintf(out,"%10g  %10g  %10g  %10g  %10g\n",
	      t,gyro,d[XX],d[YY],d[ZZ]); }
    else {
      fprintf(out,"%10g  %10g  %10g  %10g  %10g\n",
	      t,gyro,gvec[XX],gvec[YY],gvec[ZZ]); }
    j++;
  } while(read_next_x(status,&t,natoms,x,box));
  close_trj(status);
  
  fclose(out);

  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stdout);
  
  return 0;
}
