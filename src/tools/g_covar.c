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
    "All structures are fitted to the structure in the generic structure",
    "file. When this is not a run input file, both the fit and the",
    "covariance analysis will be non mass weighted and pbc will not be",
    "taken into account. When the fit and analysis groups are identical",
    "and the analysis is non mass weighted, the fit will also be non mass",
    "weighted, unless [TT]-nomwa[tt] is used.[PAR]",
    "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
    "When the same atoms are used for the fit and the covariance analysis,",
    "the reference structure is written first with t=-1.",
    "The average structure is written with t=0, the eigenvectors",
    "are written as frames with the eigenvector number as timestamp.",
    "The eigenvectors can be analyzed with [TT]g_anaeig[tt]."
  };
  static bool bM=FALSE;
  static int  begin=1,end=-1;
  t_pargs pa[] = {
    { "-mwa",  FALSE, etBOOL, &bM,
      "Mass weighted covariance analysis"},
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
  rvec       *x,*xread,*xref,*xav;
  matrix     box,zerobox;
  real       t,*mat,dev,trace,sum,*rdum1,*rdum2,rdum,avmass,inv_nframes;
  real       xj,*mass,*w_rls,lambda;
  int        ntopatoms,step;
  int        natoms,nat,ndim,count,nframes;
  char       *grpname,*stxfile,*infile,str[STRLEN];
  int        i,j,k,l,d,dj,nfit;
  atom_id    *index,*all_at,*ifit;
  bool       bTop,bDiffMass1,bDiffMass2;
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL, ffREAD }, 
    { efSTX, "-s", "topol.tpr", ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, "eigenval", ffWRITE },
    { efTRN, "-v", "eigenvec", ffWRITE },
    { efSTO, "-av", "average.pdb", ffWRITE }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  clear_mat(zerobox);

  stxfile=ftp2fn(efSTX,NFILE,fnm);

  bTop = (fn2ftp(stxfile)==efTPR) || (fn2ftp(stxfile)==efTPB) || 
    (fn2ftp(stxfile)==efTPA);
  if (bTop) {
    read_tpxheader(stxfile,&header);
    
    ntopatoms=header.natoms;
    snew(xref,ntopatoms);
    
    read_tpx(stxfile,&step,&t,&lambda,&ir,box,
	     &ntopatoms,xref,NULL,NULL,&top);
  }
  else {
    bM = FALSE;
    get_stx_coordnum(stxfile,&ntopatoms);
    init_t_atoms(&top.atoms,ntopatoms,FALSE);
    snew(xref,ntopatoms);
    read_stx_conf(stxfile,str,&top.atoms,xref,NULL,box);
    fprintf(stderr,"Note: will do a non mass weighted fit\n");
  }
  printf("\nChoose a group for the least squares fit\n"); 
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	    &nfit,&ifit,&grpname);
  printf("\nChoose a group for the covariance analysis\n"); 
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
	    &natoms,&index,&grpname);

  if (nfit < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");

  bDiffMass1=FALSE;
  snew(w_rls,ntopatoms);
  for(i=0; (i<nfit); i++)
    if (bTop) {
      w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;
      if (i)
        bDiffMass1 = bDiffMass1 || (w_rls[ifit[i]]!=w_rls[ifit[i-1]]);
    }
    else
      w_rls[ifit[i]]=1.0;
  
  bDiffMass2=FALSE;
  snew(mass,natoms);
  for(i=0; (i<natoms); i++)
    if (bM) {
      mass[i]=top.atoms.atom[index[i]].m;
      if (i)
	bDiffMass2 = bDiffMass2 || (mass[i]!=mass[i-1]);
    }
    else
      mass[i]=1.0;
  
  if (bDiffMass1 && !bDiffMass2 && !opt2parg_bSet("-mwa",asize(pa),pa)) {
    bDiffMass1 = natoms != nfit;
    i=0;
    for (i=0; (i<natoms) && !bDiffMass1; i++)
      bDiffMass1 = index[i] != ifit[i];
    if (!bDiffMass1) {
      fprintf(stderr,"\nNote: the fit and analysis group are identical, while the fit is mass weighted\n"
	               "      and the analysis is not. Making the fit non mass weighted.\n"
	               "      If you don't want this, run again with -nomwa.\n\n");
      for(i=0; (i<nfit); i++)
	w_rls[ifit[i]]=1.0;
    }
  }

  snew(all_at,top.atoms.nr);
  for(i=0; (i<top.atoms.nr); i++)
    all_at[i]=i;

  /* Prepare reference frame */
  if (bTop)
    rm_pbc(&(top.idef),top.atoms.nr,box,xref,xref);
  reset_x(nfit,ifit,top.atoms.nr,all_at,xref,w_rls);

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
    if (bTop)
      rm_pbc(&(top.idef),nat,box,xread,xread);
    reset_x(nfit,ifit,nat,all_at,xread,w_rls);
    do_fit(nat,w_rls,xref,xread);
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
  avmass=0;
  for (i=0; i<natoms; i++) {
    svmul(inv_nframes,xav[i],xav[i]);
    mass[i]=sqrt(mass[i]);
    avmass+=mass[i];
    copy_rvec(xav[i],xread[index[i]]);
  }
  avmass=sqr(natoms/avmass);

  write_sto_conf_indexed(opt2fn("-av",NFILE,fnm),"Average structure",
			   &(top.atoms),xread,NULL,NULL,natoms,index);
  for (j=0; j<natoms; j++) 
    for (dj=0; dj<DIM; dj++) 
      for (i=j; i<natoms; i++) { 
	k=ndim*(DIM*j+dj)+DIM*i;
	for (d=0; d<DIM; d++) {
	  mat[k+d]=(mat[k+d]*inv_nframes-xav[i][d]*xav[j][dj])
	    *mass[i]*mass[j]*avmass;
	}
      }
  for (j=0; j<ndim; j++) 
    for (i=j; i<ndim; i++)
      mat[ndim*i+j]=mat[ndim*j+i];
  
  trace=0;
  for(i=0; i<ndim; i++)
    trace+=mat[i*ndim+i];
  fprintf(stderr,"\nTrace of the covariance matrix: %g (nm^2)\n",trace);

  /* call diagonalization routine. Tested only fortran double precision */

  fprintf(stderr,"\nDiagonalizing...\n");
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

  sum=0;
  for(i=0; i<ndim; i++)
    sum+=rdum1[i];
  fprintf(stderr,"\nSum of the eigenvalues: %g (nm^2)\n",sum);
  if (fabs(trace-sum)>0.01*trace)
    fprintf(stderr,"\nWARNING: eigenvalue sum deviates from the trace of the covariance matrix\n");
  
  fprintf (stderr,"\nWriting eigenvalues to %s\n",ftp2fn(efXVG,NFILE,fnm));

  out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
	       "Eigenvalues of the covariance matrix",
	       "Eigenvector index","(nm\\S2\\N)");  

  for (i=0; (i<ndim); i++) {
    fprintf (out,"%10d %g\n",i+1,rdum1[ndim-1-i]);
  }
  fclose(out);  

  if (end==-1)
    end=ndim;
  fprintf (stderr,
	   "Writing %saverage structure\nand eigenvectors %d to %d to %s\n",
	   (nfit==natoms) ? "reference and " : "",
	   begin,end,opt2fn("-v",NFILE,fnm));

  trjout = open_tpx(opt2fn("-v",NFILE,fnm),"w");
  if (nfit==natoms) {
    for(i=0; i<nfit; i++)
      copy_rvec(xref[ifit[i]],x[i]);
    fwrite_trn(trjout,-1,-1,bDiffMass1 ? 1 : 0,zerobox,natoms,x,NULL,NULL);
  }
  fwrite_trn(trjout,0,0,bDiffMass2 ? 1 : 0,zerobox,natoms,xav,NULL,NULL);
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
  





