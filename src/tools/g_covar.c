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
#include "confio.h"
#include "trnio.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "txtdump.h"
#include "callf77.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_covar[tt] calculates and diagonalizes the (mass weighted)",
    "covariance matrix.",
    "All structures are fitted to the structure in the run input file.[PAR]",
    "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
    "When the same atoms are used for the fit and the covariance analysis,",
    "the reference structure is written first with t=-1.",
    "The average structure is written with t=0, the eigenvectors",
    "are written as frames with the eigenvector number as timestamp.",
    "The eigenvectors can be analyzed with [TT]g_anaeig[tt]."
  };
  static bool bM=TRUE;
  static int  begin=1,end=-1;
  t_pargs pa[] = {
    { "-m",  FALSE, etBOOL, &bM,
      "Calculate mass weighted eigenvectors"},
    { "-first", FALSE, etINT, &begin,     
      "first eigenvector to write away" },
    { "-last",  FALSE, etINT, &end, 
      "last eigenvector to write away (-1 is till the last)" }
  };
  FILE       *out;
  int        status,trjout;
  t_tpxheader  header;
  t_inputrec   ir;
  t_topology   top;
  t_idef       *idef;
  rvec       *x,*xread,*xref,*xav;
  matrix     box,zerobox;
  real       t,*mat,dev,*rdum1,*rdum2,rdum,avmass,inv_nframes;
  real       xid,*mass,*w_rls,lambda;
  int        ntopatoms,step;
  int        natoms,nat,ndim,count,nframes;
  char       *grpname,*infile;
  int        i,j,k,l,d,d2,nfit;
  atom_id    *index,*all_at,*ifit;
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL, ffREAD }, 
    { efTPX, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, "eigval", ffWRITE },
    { efTRN, "-v", "eigvec", ffWRITE },
    { efSTO, "-av", "average.pdb", ffWRITE }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  clear_mat(zerobox);
  /*
    if (bM)
    fprintf(out,"@ subtitle \"of mass weighted Hessian matrix\"\n");
    else 
    fprintf(out,"@ subtitle \"of Hessian matrix\"\n");
    */
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&header);

  ntopatoms=header.natoms;
  snew(xref,ntopatoms);

  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,box,
	   &ntopatoms,xref,NULL,NULL,&top);

  printf("\nChoose a group for the least squares fit\n"); 
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	    &nfit,&ifit,&grpname);
  printf("\nChoose a group for the covariance analysis\n"); 
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	    &natoms,&index,&grpname);

  if (nfit < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");

  snew(w_rls,ntopatoms);
  for(i=0; (i<nfit); i++)
    w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;

  snew(mass,natoms);
  for(i=0; (i<natoms); i++)
    if (bM)
      mass[i]=top.atoms.atom[index[i]].m;
    else
      mass[i]=1.0;

  /* Prepare reference frame */
  reset_x(nfit,ifit,nfit,ifit,xref,w_rls);

  infile=opt2fn("-f",NFILE,fnm);

  snew(x,natoms);
  snew(xav,natoms);
  ndim=natoms*DIM;
  snew(mat,ndim*ndim);

  fprintf(stderr,"Constructing covariance matrix (%dx%d)...\n",ndim,ndim);
  
  nframes=0;
  nat=read_first_x(&status,infile,&t,&xread,box);
  snew(all_at,nat);
  for(i=0; (i<nat); i++)
    all_at[i]=i;
  do {
    nframes++;
    /* calculate x: a fitted struture of the selected atoms */ 
    reset_x(nfit,ifit,nat,all_at,xread,w_rls);
    do_fit(nat,w_rls,xref,xread);
    for (i=0; i<natoms; i++)
      copy_rvec(xread[index[i]],x[i]);

    for (i=0; i<natoms; i++) {
      /* calculate average structure */
      rvec_inc(xav[i],x[i]);
      /* calculate cross product matrix */
      for (d=0; d<DIM; d++) {
	k=ndim*(DIM*i+d);
	xid=x[i][d];
	for (j=0; j<natoms; j++) {
	  l=k+DIM*j;
	  for(d2=0; d2<DIM; d2++)
	    mat[l+d2] += (double)(xid*x[j][d2]);
	}
      }
    }
  } while (read_next_x(status,&t,nat,xread,box));
  close_trj(status);

  /* calculate the mass-weighted covariance matrix */
  inv_nframes=1.0/nframes;
  avmass=0;
  for (i=0; i<natoms; i++) {
    svmul(inv_nframes,xav[i],xav[i]);
    avmass+=mass[i];
    mass[i]=sqrt(mass[i]);
    copy_rvec(xav[i],xread[index[i]]);
  }
  write_sto_conf_indexed(opt2fn("-av",NFILE,fnm),"Average structure",
			   &(top.atoms),xread,NULL,NULL,natoms,index);
  avmass=natoms/avmass;
  for (i=0; i<natoms; i++) 
    for (j=0; j<natoms; j++) 
      for (d=0; d<DIM; d++) {
	k=ndim*(DIM*i+d)+DIM*j;
	for (d2=0; d2<DIM; d2++)
	  mat[k+d2]=(mat[k+d2]*inv_nframes-xav[i][d]*xav[j][d2])
	    *mass[i]*mass[j]*avmass;
      }

  /* check if matrix is approximately symmetric */
  for (i=0; (i<ndim); i++) {
    for (j=0; (j<i); j++) {
      if (mat[i*ndim+j] != 0.0 && mat[j*ndim+i] != 0.0)
	rdum=(mat[i*ndim+j]-mat[j*ndim+i])/(mat[i*ndim+j]+mat[j*ndim+i]);
      else {
	if (fabs(mat[i*ndim+j] - mat[j*ndim+i]) > 1.0e-5)
	  rdum=1.0;      
	else
	  rdum=0.0;
      }
      if (abs(rdum)>1.0e-5) {
	fprintf(stderr,"possible non-symmetric matrix:\n");
	fprintf(stderr,"x[%d][%d]=%e,x[%d][%d]=%e\n",i,j,\
		mat[i*ndim+j],j,i,mat[j*ndim+i]);
      }
      /* make sure that it is symmetric 
       * this is not very nice, but we did warn, did't we? */
      mat[j*ndim+i]=mat[i*ndim+j];
    }
  }
	  
  /* call diagonalization routine. Tested only fortran double precision */

  fprintf(stderr,"Diagonalizing...\n");
  fflush(stderr);

  snew(rdum1,ndim);
  snew(rdum2,ndim);
#ifdef USEF77
  fql77(&ndim,mat,rdum1,rdum2,&ndim);
#else
  /*ql77(int n,real **x,real *d,real *e,int nmax)*/
  /*fprintf(stderr,"Calling ql77...\n");
    ql77 (ndim,mat,rdum1,rdum2,ndim);*/
  fatal_error(0,"C version of ql77 buggy. Use f77. Sorry.");
#endif
  
  /* now write the output */

  fprintf (stderr,"Writing eigenvalues to %s\n",ftp2fn(efXVG,NFILE,fnm));

  out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
                 "Eigenvalues","Eigenvector index","Eigenvalue");  

  for (i=0; (i<ndim); i++) {
    fprintf (out,"%10d %15.6e\n",i+1,rdum1[ndim-1-i]);
  }
  fclose(out);  

  if (end==-1)
    end=ndim;
  fprintf (stderr,
	   "Writing %saverage structure\nand eigenvectors %d to %d to %s\n",
	   (nfit==natoms) ? "reference and " : "",
	   begin,end,opt2fn("-v",NFILE,fnm));

  trjout = open_tpx(opt2fn("-v",NFILE,fnm),"w");
  if (nfit==natoms)
    fwrite_trn(trjout,-1,-1,0,zerobox,natoms,xref,NULL,NULL);
  fwrite_trn(trjout,0,0,0,zerobox,natoms,xav,NULL,NULL);
  for(i=begin; i<=end; i++) {
    for (j=0; j<natoms; j++)
      for(d=0; d<DIM; d++)
	x[j][d]=mat[(ndim-i)*ndim+DIM*j+d];
    fwrite_trn(trjout,i,(real)i,0,zerobox,natoms,x,NULL,NULL);
  }
  close_trn(trjout);

  thanx(stdout);
  
  return 0;
}
  
