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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_nmeig_c = "$Id$";

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

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_nmeig calculates the eigenvectors/values of a (Hessian) matrix."
  };
  static bool bM=FALSE;
  static int  begin=1,end=100;
  t_pargs pa[] = {
    { "-m",  FALSE, etBOOL, &bM,
      "Divide elements of Hessian by product of sqrt(mass) of involved atoms prior to diagonalization. This should be used for 'Normal Modes' analyses" },
    { "-first", FALSE, etINT, &begin,     
      "first eigenvector to write away" },
    { "-last",  FALSE, etINT, &end, 
      "last eigenvector to write away" }
  };
  FILE       *out,*vec;
  int        status;
  t_topology *top;
  rvec       *x;
  matrix     box;
  real       t,*hess,*rdum1,*rdum2,rdum;
  int        natoms,ndim,count;
  char       *grpname,*vecnm;
  int        i,j,k,l,gnx;
  bool       bSuck,bEnd;
  atom_id    *index;
  t_filenm fnm[] = { 
    { efMTX, "-f", "hessian", ffREAD }, 
    { efTPX, NULL, NULL, ffREAD },
    { efXVG, NULL, "eigval", ffWRITE },
    { efVEC, "-k", "eigvec", ffWRITE }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 
  bEnd=opt2parg_bSet("-last",asize(pa),pa);

  out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
                 "Eigenvalues","Eigenvector index","Eigenvalue"); 

  if (bM)
    fprintf(out,"@ subtitle \"of mass weighted Hessian matrix\"\n");
  else 
    fprintf(out,"@ subtitle \"of Hessian matrix\"\n");

  top=read_top(ftp2fn(efTPX,NFILE,fnm)); 

  vecnm=ftp2fn(efVEC,NFILE,fnm);
  vec=(FILE *)ffopen(vecnm,"w");

  /*open Hessian matrix and read first 'line' */

  natoms=read_first_x(&status,ftp2fn(efMTX,NFILE,fnm),&t,&x,box); 
  ndim=natoms*DIM;

  fprintf(stderr,"Dimensionality of matrix: %d\n",ndim);

  /* store two dimensional Hessian matrix in one dimensional array 
   * to facilitate communication with fortran */

  snew(hess,ndim*ndim);

  snew(rdum1,ndim);
  snew(rdum2,ndim);

  for (i=0; (i<natoms); i++) {
    for (j=0; (j<DIM); j++) {
      hess[i*DIM+j]=0.0-x[i][j];
    }
  }

  /* read in rest of Hessian */
  
  for (i=1; (i<ndim); i++) {
    if (read_next_x(status,&t,natoms,x,box)) {
      for (j=0; (j<natoms); j++) {
	for (k=0; (k<DIM); k++) {
	  hess[i*ndim+j*DIM+k]=0.0-x[j][k];
	}
      }
    }
    else {
      fatal_error(0,"Premature end of file, hessian matrix");
    }
  }
  close_trj(status);


  /* check if matrix is approximately symmetric */
  for (i=0; (i<ndim); i++) {
    for (j=0; (j<i); j++) {
      if (hess[i*ndim+j] != 0.0 && hess[j*ndim+i] != 0.0)
	rdum=(hess[i*ndim+j]-hess[j*ndim+i])/(hess[i*ndim+j]+hess[j*ndim+i]);
      else {
	if (fabs(hess[i*ndim+j] - hess[j*ndim+i]) > 1.0e-5)
	  rdum=1.0;      
	else
	  rdum=0.0;
      }
      if (abs(rdum)>1.0e-5) {
	fprintf(stderr,"possible non-symmetric matrix:\n");
	fprintf(stderr,"x[%d][%d]=%e,x[%d][%d]=%e\n",i,j,\
		hess[i*ndim+j],j,i,hess[j*ndim+i]);
      }
      /* make sure that it is symmetric 
       * this is not very nice, but we did warn, did't we? */
      hess[j*ndim+i]=hess[i*ndim+j];
    }
  }
	  
  /* divide elements hess[i][j] 
     by sqrt(mas[i])*sqrt(mas[j]) when required */

  if (bM) {
    for (i=0; (i<natoms); i++) {
      for (j=0; (j<DIM); j++) {
	for (k=0; (k<natoms); k++) {
	  for (l=0; (l<DIM); l++) {
	    hess[(i*DIM+j)*ndim+k*DIM+l]/=
	      (sqrt(top->atoms.atom[i].m)*sqrt(top->atoms.atom[k].m));
	  }
	}
      }
    }
  }

  /* call diagonalization routine. Tested only fortran double precision */

  fprintf(stderr,"Diagonalizing...\n");
  fflush(stderr);

#ifdef USEF77
  CALLF77(fql77)(&ndim,hess,rdum1,rdum2,&ndim);
#else
  /*ql77(int n,real **x,real *d,real *e,int nmax)*/
  /*fprintf(stderr,"Calling ql77...\n");
    ql77 (ndim,hess,rdum1,rdum2,ndim);*/
  fatal_error(0,"C version of ql77 buggy. Use f77. Sorry.");
#endif
  
  /* check the output, first 6 eigenvalues should be quite small */

  bSuck=FALSE;
  for (i=0; (i<6); i++) {
    if (abs(rdum1[i]) > 1.0e-3) 
      bSuck=TRUE;
  }
  if (bSuck) {
    fprintf(stderr,"\nOne of the first eigenvalues has a non-zero value.\n");
    fprintf(stderr,"This could mean that the reference structure was not\n");
    fprintf(stderr,"properly energy minimised\n");
  }

  /* now write the output */

  fprintf (stderr,"Writing eigenvalues...\n");

  for (i=0; (i<ndim); i++) {
    fprintf (out,"%10d %15.6e\n",i+1,rdum1[i]);
  }
  fclose(out);  

  fprintf (stderr,"Writing eigenvectors %d to %d...\n",begin,end);

  count=0;
  for (i=begin-1; (i<end); i++) {
    for (j=0; (j<ndim); j++) {
      count++;
      fprintf (vec,"%12.5e",hess[i*ndim+j]);
      if ((count % 6) == 0)
	fprintf (vec,"\n");
    }
    if ((count % 6) != 0) 
      fprintf (vec,"\n");
    count=0;
  }

  fclose(vec);  

  fprintf(stderr,"\n");

  thanx(stdout);
  
  return 0;
}
