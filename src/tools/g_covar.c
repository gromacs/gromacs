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
#include "do_fit.h"
#include "rmpbc.h"
#include "txtdump.h"
#include "eigio.h"
#include "ql77.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_covar[tt] calculates and diagonalizes the (mass-weighted)",
    "covariance matrix.",
    "All structures are fitted to the structure in the structure file.",
    "When this is not a run input file periodicity will not be taken into",
    "account. When the fit and analysis groups are identical and the analysis",
    "is non mass-weighted, the fit will also be non mass-weighted.[PAR]",
    "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
    "When the same atoms are used for the fit and the covariance analysis,",
    "the reference structure is written first with t=-1.",
    "The average structure is written with t=0, the eigenvectors",
    "are written as frames with the eigenvector number as timestamp.",
    "The eigenvectors can be analyzed with [TT]g_anaeig[tt]."
  };
  static bool bFit=TRUE,bM=FALSE;
  static int  end=-1;
  t_pargs pa[] = {
    { "-fit",  FALSE, etBOOL, {&bFit},
      "Fit to a reference structure"},
    { "-mwa",  FALSE, etBOOL, {&bM},
      "Mass-weighted covariance analysis"},
    { "-last",  FALSE, etINT, {&end}, 
      "Last eigenvector to write away (-1 is till the last)" }
  };
  FILE       *out;
  int        status,trjout;
  t_topology top;
  t_atoms    *atoms;  
  rvec       *x,*xread,*xref,*xav;
  matrix     box,zerobox;
  real       t,*mat,dev,trace,sum,*eigval,inv_nframes;
  real       xj,*sqrtm,*w_rls=NULL;
  int        ntopatoms,step;
  int        natoms,nat,ndim,count,nframes;
  int        WriteXref;
  char       *grpname,*infile,str[STRLEN];
  int        i,j,k,l,d,dj,nfit;
  atom_id    *index,*all_at,*ifit;
  bool       bDiffMass1,bDiffMass2;
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL, ffREAD }, 
    { efTPS, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, "eigenval", ffWRITE },
    { efTRN, "-v", "eigenvec", ffWRITE },
    { efSTO, "-av", "average.pdb", ffWRITE }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_SET_NPRI,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  clear_mat(zerobox);

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),str,&top,&xref,NULL,box,TRUE);
  atoms=&top.atoms;

  if (bFit) {
    printf("\nChoose a group for the least squares fit\n"); 
    get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	      &nfit,&ifit,&grpname);
    if (nfit < 3) 
      fatal_error(0,"Need >= 3 points to fit!\n");
  } else
    nfit=0;
  printf("\nChoose a group for the covariance analysis\n"); 
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	    &natoms,&index,&grpname);

  bDiffMass1=FALSE;
  if (bFit) {
    snew(w_rls,atoms->nr);
    for(i=0; (i<nfit); i++) {
      w_rls[ifit[i]]=atoms->atom[ifit[i]].m;
      if (i)
        bDiffMass1 = bDiffMass1 || (w_rls[ifit[i]]!=w_rls[ifit[i-1]]);
    }
  }
  bDiffMass2=FALSE;
  snew(sqrtm,natoms);
  for(i=0; (i<natoms); i++)
    if (bM) {
      sqrtm[i]=sqrt(atoms->atom[index[i]].m);
      if (i)
	bDiffMass2 = bDiffMass2 || (sqrtm[i]!=sqrtm[i-1]);
    }
    else
      sqrtm[i]=1.0;
  
  if (bFit && bDiffMass1 && !bDiffMass2) {
    bDiffMass1 = natoms != nfit;
    i=0;
    for (i=0; (i<natoms) && !bDiffMass1; i++)
      bDiffMass1 = index[i] != ifit[i];
    if (!bDiffMass1) {
      fprintf(stderr,"\n"
	      "Note: the fit and analysis group are identical,\n"
	      "      while the fit is mass weighted and the analysis is not.\n"
	      "      Making the fit non mass weighted.\n\n");
      for(i=0; (i<nfit); i++)
	w_rls[ifit[i]]=1.0;
    }
  }

  snew(all_at,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    all_at[i]=i;

  /* Prepare reference frame */
  rm_pbc(&(top.idef),atoms->nr,box,xref,xref);
  if (bFit)
    reset_x(nfit,ifit,atoms->nr,all_at,xref,w_rls);

  infile=opt2fn("-f",NFILE,fnm);

  snew(x,natoms);
  snew(xav,natoms);
  ndim=natoms*DIM;
  snew(mat,ndim*ndim);

  fprintf(stderr,"Constructing covariance matrix (%dx%d)...\n",ndim,ndim);

  nframes=0;
  nat=read_first_x(&status,infile,&t,&xread,box);
  do {
    nframes++;
    /* calculate x: a fitted struture of the selected atoms */
    rm_pbc(&(top.idef),nat,box,xread,xread);
    if (bFit) {
      reset_x(nfit,ifit,nat,all_at,xread,w_rls);
      do_fit(nat,w_rls,xref,xread);
    }
    for (i=0; i<natoms; i++)
      copy_rvec(xread[index[i]],x[i]);

    for (j=0; j<natoms; j++) {
      /* calculate average structure */
      rvec_inc(xav[j],x[j]);
      /* calculate cross product matrix */
      for (dj=0; dj<DIM; dj++) {
	k=ndim*(DIM*j+dj);
	xj=x[j][dj];
	for (i=j; i<natoms; i++) {
	  l=k+DIM*i;
	  for(d=0; d<DIM; d++)
	    mat[l+d] += x[i][d]*xj;
	}
      }
    }
  } while (read_next_x(status,&t,nat,xread,box));
  close_trj(status);

  /* calculate the mass-weighted covariance matrix */
  inv_nframes=1.0/nframes;
  for (i=0; i<natoms; i++) {
    svmul(inv_nframes,xav[i],xav[i]);
    copy_rvec(xav[i],xread[index[i]]);
  }

  write_sto_conf_indexed(opt2fn("-av",NFILE,fnm),"Average structure",
			 atoms,xread,NULL,NULL,natoms,index);
  
  for (j=0; j<natoms; j++) 
    for (dj=0; dj<DIM; dj++) 
      for (i=j; i<natoms; i++) { 
	k=ndim*(DIM*j+dj)+DIM*i;
	for (d=0; d<DIM; d++) {
	  mat[k+d]=(mat[k+d]*inv_nframes-xav[i][d]*xav[j][dj])
	    *sqrtm[i]*sqrtm[j];
	}
      }
  for (j=0; j<ndim; j++) 
    for (i=j; i<ndim; i++)
      mat[ndim*i+j]=mat[ndim*j+i];
  
  trace=0;
  for(i=0; i<ndim; i++)
    trace+=mat[i*ndim+i];
  fprintf(stderr,"\nTrace of the covariance matrix: %g (%snm^2)\n",
	  trace,bM ? "amu " : "");

  if (debug) {
    fprintf(stderr,"Dumping the covariance matrix to covmat.dat\n"); 
    out = ffopen("covmat.dat","w");
    for (j=0; j<ndim; j++) {
      for (i=0; i<ndim; i++)
	fprintf(out," %g",mat[ndim*j+i]);
      fprintf(out,"\n");
    }
  }

  /* call diagonalization routine */

  fprintf(stderr,"\nDiagonalizing...\n");
  fflush(stderr);

  snew(eigval,ndim);
  ql77 (ndim,mat,eigval);
  
  /* now write the output */

  sum=0;
  for(i=0; i<ndim; i++)
    sum+=eigval[i];
  fprintf(stderr,"\nSum of the eigenvalues: %g (%snm^2)\n",
	  sum,bM ? "amu " : "");
  if (fabs(trace-sum)>0.01*trace)
    fprintf(stderr,"\nWARNING: eigenvalue sum deviates from the trace of the covariance matrix\n");
  
  fprintf (stderr,"\nWriting eigenvalues to %s\n",ftp2fn(efXVG,NFILE,fnm));

  sprintf(str,"(%snm\\S2\\N)",bM ? "amu " : "");
  out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
	       "Eigenvalues of the covariance matrix",
	       "Eigenvector index",str);  

  for (i=0; (i<ndim); i++) {
    fprintf (out,"%10d %g\n",i+1,eigval[ndim-1-i]);
  }
  fclose(out);  

  if (end==-1) {
    if (nframes-1 < ndim)
      end=nframes-1;
    else
      end=ndim;
  }
  if (bFit) {
    /* misuse lambda: 0/1 mass weighted analysis no/yes */
    if (nfit==natoms) {
      WriteXref = eWXR_YES;
      for(i=0; i<nfit; i++)
	copy_rvec(xref[ifit[i]],x[i]);
    } else
      WriteXref = eWXR_NO;
  } else {
    /* misuse lambda: -1 for no fit */
    WriteXref = eWXR_NOFIT;
  }

  write_eigenvectors(opt2fn("-v",NFILE,fnm),natoms,mat,TRUE,1,end,
		     WriteXref,x,bDiffMass1,xav,bM);

  thanx(stdout);
  
  return 0;
}
  





