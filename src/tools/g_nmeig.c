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
 * GRowing Old MAkes el Chrono Sweat
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
#include "ql77.h"
#include "eigio.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_nmeig calculates the eigenvectors/values of a (Hessian) matrix,",
    "which can be calculated with [TT]nmrun[tt].",
    "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
    "The structure is written first with t=0. The eigenvectors",
    "are written as frames with the eigenvector number as timestamp.",
    "The eigenvectors can be analyzed with [TT]g_anaeig[tt].",
    "An ensemble of structures can be generated from the eigenvectors with",
    "[TT]g_nmens[tt]."
  };
  static bool bM=TRUE;
  static int  begin=1,end=100;
  t_pargs pa[] = {
    { "-m",  FALSE, etBOOL, {&bM},
      "Divide elements of Hessian by product of sqrt(mass) of involved "
      "atoms prior to diagonalization. This should be used for 'Normal Modes' "
      "analysis" },
    { "-first", FALSE, etINT, {&begin},     
      "First eigenvector to write away" },
    { "-last",  FALSE, etINT, {&end}, 
      "Last eigenvector to write away" }
  };
  FILE       *out;
  int        status,trjout;
  t_topology top;
  t_inputrec ir;
  rvec       *top_x,*x;
  matrix     box;
  real       t,*hess,*eigv,rdum,mass_fac;
  int        natoms,ndim,count;
  char       *grpname,title[256];
  int        i,j,k,l,d,gnx;
  bool       bSuck;
  atom_id    *index;
  t_filenm fnm[] = { 
    { efMTX, "-f", "hessian",  ffREAD  }, 
    { efTPS, NULL, NULL,       ffREAD  },
    { efXVG, NULL, "eigenval", ffWRITE },
    { efTRN, "-v", "eigenvec", ffWRITE }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_SET_NPRI,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&top_x,NULL,box,bM);

  /*open Hessian matrix and read first 'line' */

  natoms=read_first_x(&status,ftp2fn(efMTX,NFILE,fnm),&t,&x,box); 
  ndim=natoms*DIM;

  fprintf(stderr,"Dimensionality of matrix: %d\n",ndim);

  /* store two dimensional Hessian matrix in one dimensional array 
   * to facilitate communication with fortran */

  snew(hess,ndim*ndim);

  snew(eigv,ndim);

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
      if (fabs(rdum)>1.0e-5) {
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
	  mass_fac=invsqrt(top.atoms.atom[i].m*top.atoms.atom[k].m);
	  for (l=0; (l<DIM); l++)
	    hess[(i*DIM+j)*ndim+k*DIM+l]*=mass_fac;
	}
      }
    }
  }

  /* call diagonalization routine. Tested only fortran double precision */

  fprintf(stderr,"\nDiagonalizing...\n");
  fflush(stderr);

  ql77 (ndim,hess,eigv);
  
  /* check the output, first 6 eigenvalues should be quite small */

  bSuck=FALSE;
  for (i=0; (i<6); i++) {
    if (fabs(eigv[i]) > 1.0e-3) 
      bSuck=TRUE;
  }
  if (bSuck) {
    fprintf(stderr,"\nOne of the first eigenvalues has a non-zero value.\n");
    fprintf(stderr,"This could mean that the reference structure was not\n");
    fprintf(stderr,"properly energy minimised\n");
  }

  /* now write the output */
  fprintf (stderr,"Writing eigenvalues...\n");
  out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
	       "Eigenvalues","Eigenvector index","Eigenvalue");
  if (bM)
    fprintf(out,"@ subtitle \"of mass weighted Hessian matrix\"\n");
  else 
    fprintf(out,"@ subtitle \"of Hessian matrix\"\n");
  for (i=0; i<ndim; i++)
    fprintf (out,"%6d %15g\n",i+1,eigv[i]);
  fclose(out);
  
  if (begin<1)
    begin=1;
  if (end>ndim)
    end=ndim;

  write_eigenvectors(opt2fn("-v",NFILE,fnm),natoms,hess,FALSE,begin,end,
		     eWXR_NO,NULL,FALSE,top_x,bM);

  thanx(stderr);
  
  return 0;
}
