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
static char *SRCID_g_rmsdist_c = "$Id$";

#include <math.h>
#include <string.h>
#include "macros.h"
#include "smalloc.h"
#include "typedefs.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "string2.h"
#include "vec.h"
#include "macros.h"
#include "rdgroup.h"
#include "pbc.h"
#include "xvgr.h"
#include "futil.h"
#include "matio.h"
#include "callf77.h"

static void calc_dist(int nind,atom_id index[],
		      rvec x[],matrix box,real **d)
{
  int     i,j;
  real    *xi;
  rvec    dx;

  for(i=0; (i<nind-1); i++) {
    xi=x[index[i]];
    for(j=i+1; (j<nind); j++) {
      pbc_dx(xi,x[index[j]],dx);
      d[i][j]=norm(dx);
    }
  }
}

static real rms_diff(int isize,real **d,real **d_r)
{
  int  i,j;
  real r,r2;
  
  r2=0.0;
  for(i=0; (i<isize-1); i++)
    for(j=i+1; (j<isize); j++) {
      r=d[i][j]-d_r[i][j];
      r2+=r*r;
    }    
  r2/=(isize*(isize-1))/2;
  
  return sqrt(r2);
}

static real **mcpy=NULL;
static int  nnn=0;

int rcomp(const void *a,const void *b)
{
  real **ra,**rb;
  
  ra = (real **)a;
  rb = (real **)b;
  
  if ((*ra)[nnn] < (*rb)[nnn])
    return -1;
  else if ((*ra)[nnn] > (*rb)[nnn])
    return 1;
    
  return 0;
}

static int icomp(const void *a,const void *b)
{
  int ia,ib;
  
  ia = *(int *)a;
  ib = *(int *)b;
  
  if (mcpy[ia][0] < mcpy[ib][0])
    return -1;
  else if (mcpy[ia][0] > mcpy[ib][0])
    return 1;
  else
    return 0;  
}

void sort_matrix(int n1,real **mat)
{
  int i,j,*index;

  snew(index,n1);
  snew(mcpy,n1);
  for(i=0; (i<n1); i++) {
    index[i] = i;
    snew(mcpy[i],n1);
    for(j=0; (j<n1); j++)
      mcpy[i][j] = mat[i][j];
  }
  fprintf(stderr,"Going to sort the RMSD matrix.\n");  
  qsort(index,n1,sizeof(*index),icomp);
  /* Copy the matrix back */
  for(i=0; (i<n1); i++)
    for(j=0; (j<n1); j++)
      mat[i][j] = mcpy[index[i]][index[j]];
  for(i=0; (i<n1); i++)
    sfree(mcpy[i]);
  sfree(mcpy);
  sfree(index);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_rmsdist computes the root mean square deviation of atom distances,",
    "which has the advantage that no fit is needed like in standard RMS",
    "deviation as computed by g_rms.",
    "The reference structure is taken from the structure file.",
    "The rmsd at time t is calculated as the rms",
    "of the differences in distance between atom-pairs in the reference",
    "structure and the structure at time t.[PAR]",
    "g_rmsdist can also produce matrices of the rms distances, rms distances",
    "scaled with the mean distance and the mean distances and matrices with",
    "NMR averaged distances (1/r^3 and 1/r^6 averaging)."
  };
  
  FILE         *fp;
  int          natom1,natom2,i1,i2,i1s,i2s,i,teller,nf;
  real         t,t1,t2;

  matrix       box;
  rvec         *x1,*x2;
  char         *fn1,*fn2;
  real         **rms,maxrms,*resnr;
  real         *rms0,*eigval;
  t_topology   top;
  bool         bSameF;
  t_rgb        rlo,rhi;
  
  int      status1,status2,isize;
  atom_id  *index;
  char     *grpname;
  real     **d1,**d2;
  char     buf[255];
  
  static int  nlevels=40,nframes=100,skip=1;
  static real scalemax=-1.0,rmscut=0.1;
  static bool bEigen=FALSE;
  t_pargs pa[] = {
    { "-nlevels",   FALSE, etINT,  &nlevels,
      "Discretize rms in # levels" },
    { "-rmscut", FALSE, etREAL, &rmscut,
      "RMS cut-off (nm) for two structures to be similar" },
    { "-max",   FALSE, etREAL, &scalemax,
      "Maximum level in matrices" },
    { "-nframes", FALSE, etINT, &nframes,
      "Max number of frames" },
    { "-skip", FALSE, etINT, &skip,
      "Skip N frames after each one used." },
    { "-eigen", FALSE, etBOOL, &bEigen,
      "Diagonalize the rms matrix instead of sorting it" }
  };
  t_filenm fnm[] = {
    { efTRX, "-f1",   NULL,       ffREAD  },
    { efTRX, "-f2",   NULL,       ffOPTRD },
    { efTPS, "-s",    NULL,       ffREAD  },
    { efNDX, NULL,   NULL,        ffOPTRD },
    { efXPM, "-o",   "rmsd",      ffWRITE },
    { efXPM, "-os",  "rmsd-mat",  ffWRITE },
    { efXVG, "-ev",  "rms-eig",   ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  if (skip < 0)
    fatal_error(0,"skip (%d) should be >= 0. 0 means no skip",skip);
  skip++;
  
  /* don't read mass-database as masses (and top) are not used */
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),buf,&top,&x1,NULL,box,FALSE);
  
  get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    1,&isize,&index,&grpname);

  /* Initiate arrays */  
  snew(d1,isize);
  snew(d2,isize);
  for(i=0; (i<isize); i++) {
    snew(d1[i],isize);
    snew(d2[i],isize);
  }
  nf = nframes/skip+1;
  snew(rms,nf);
  if (bEigen)
    snew(rms0,nf*nf); 

  snew(resnr,nf);
  for(i=0; (i<nf); i++) {
    if (bEigen)
      rms[i] = &(rms0[i*nf]);
    else
      snew(rms[i],nf);
    resnr[i] = 0;
  }

  /*set box type*/
  init_pbc(box,FALSE);

  /* Loop over first coordinate file */
  fn1 = opt2fn("-f1",NFILE,fnm);
  fn2 = opt2fn("-f2",NFILE,fnm);
  bSameF = !opt2bSet("-f2",NFILE,fnm) || (strcmp(fn1,fn2) == 0);
  if (bSameF) 
    fprintf(stderr,"Reading the same file twice. "
	    "Will only do half of matrix\n");

  /* Read first frames of both files */  
  natom1 = read_first_x(&status1,fn1,&t1,&x1,box);
  natom2 = read_first_x(&status2,fn2,&t2,&x2,box);
  
  /* Consistency checks */
  if (natom1 != natom2) 
    fprintf(stderr,"Warning: natoms in %s is %d, in %s it is %d\n",
	    fn1,natom1,fn2,natom2);
  if (isize > natom1) 
    fatal_error(0,"Index (from index file or topology)"
		" does not match trajectory %s",fn1);
  if (isize > natom2) 
    fatal_error(0,"Index (from index file or topology)"
		" does not match trajectory %s",fn2);
		
  /* Initiate locals */	
  i1     = 0;
  maxrms = 0;
  nf     = 0;
  do { 
    /* Loop over first file */
    if (i1 >= nframes)
      break;
    else if ((i1 % skip) == 0) {
      /* Print counter */
      if ((i1 % 10) == 0)
	fprintf(stderr,"\ri1 = %d                               \n",i1);
      
      /* Compute distance matrix for first structure */
      calc_dist(isize,index,x1,box,d1);
      
      /* Find the starting frame for the second file */
      if (bSameF)
	i2 = i1+1;
      else
	i2 = 0;
	
      /* Rewind second traj, and skip to first interesting frame */
      rewind_trj(status2);
      for(i=0; (i<=i2); i++)
	if (!read_next_x(status2,&t2,natom2,x2,box))
	  break;
      do {
	/* Loop over second file */
	if (i2 >= nframes)
	  break;
	else if ((i2 % skip) == 0) {
	  /* Indices in RMS matrix */
	  i1s = i1/skip;
	  i2s = i2/skip;
	  
	  /* Compute distance matrix for second structure */  
	  calc_dist(isize,index,x2,box,d2);
	  
	  /* Compute RMS between two distance matrices */
	  rms[i1s][i2s] = rms_diff(isize,d1,d2);
	  
	  if (bSameF)
	    rms[i2s][i1s] = rms[i1s][i2s];
	    
	  /* Store maximum value for plotting */
	  maxrms = max(maxrms,rms[i1s][i2s]);
	}
	i2++;
      } while (read_next_x(status2,&t2,natom2,x2,box));
      nf++;
    }
    i1++;
  } while(read_next_x(status1,&t1,natom1,x1,box));
  fprintf(stderr,"\n");
  
  /* Done with computing, close the files */
  /*close_trx(status1);
  close_trx(status2);
  */

  fprintf(stderr,"The maximum RMS was %g nm, number of frames for matrix %d\n",
	  maxrms,nf);
  
  /* Set colors for plotting: white = zero RMS, black = maximum */
  rlo.r=1.0, rlo.g=1.0, rlo.b=0.0;
  rhi.r=0.0, rhi.g=0.0, rhi.b=1.0;
  
  /* Write out plot file with RMS matrix */
  write_xpm(opt2FILE("-o",NFILE,fnm,"w"),"RMS",
	    "RMS (n)","Confs 1","Confs 2",
	    nf,nf,resnr,resnr,rms,0.0,maxrms,rlo,rhi,&nlevels);

  if (!bEigen) 
    /* Now sort the matrix and write it out again */
    sort_matrix(nf,rms);
  else {
    /* Do a diagonalization */
    snew(eigval,nf);
    ql77(nf,rms0,eigval); 
    fp = xvgropen(opt2fn("-ev",NFILE,fnm),"Eigenvalues","index","value");
    for(i=0; (i<nf); i++)
      fprintf(fp,"%10d  %10g\n",i,eigval[i]);
    fclose(fp);
  }
  write_xpm(opt2FILE("-os",NFILE,fnm,"w"),"RMS Sorted",
	    "RMS (n)","Confs 1","Confs 2",
	    nf,nf,resnr,resnr,rms,0.0,maxrms,rlo,rhi,&nlevels);

	    
  /* Thank the user for her patience */  
  thanx(stdout);
  
  return 0;
}
