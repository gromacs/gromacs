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
#include "assert.h"
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

typedef struct {
  int  i,j;
  real dist;
} t_dist;

typedef struct {
  int  conf,clust;
} t_clustid;

static real **mk_matrix(int nf,bool b1D)
{
  real **rms;
  real *rms0=NULL;
  int  i;
  
  snew(rms,nf);
  if (b1D)
    snew(rms0,nf*nf); 
  
  for(i=0; (i<nf); i++) {
    if (b1D)
      rms[i] = &(rms0[i*nf]);
    else
      snew(rms[i],nf);
  }
  return rms;
}

static void calc_dist_ind(int nind,atom_id index[],
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

static void calc_dist(int nind,rvec x[],matrix box,real **d)
{
  int     i,j;
  real    *xi;
  rvec    dx;

  for(i=0; (i<nind-1); i++) {
    xi=x[i];
    for(j=i+1; (j<nind); j++) {
      pbc_dx(xi,x[j],dx);
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
  
  if (mcpy[ia][nnn] < mcpy[ib][nnn])
    return -1;
  else if (mcpy[ia][nnn] > mcpy[ib][nnn])
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

static int dcomp(const void *a,const void *b)
{
  t_dist *da,*db;
  
  da = (t_dist *)a;
  db = (t_dist *)b;
  
  if (da->dist < db->dist)
    return -1;
  else if (da->dist > db->dist)
    return 1;
  return 0;
}

static int ccomp(const void *a,const void *b)
{
  t_clustid *da,*db;
  
  da = (t_clustid *)a;
  db = (t_clustid *)b;
  
  return da->clust - db->clust;
}

void gather(int n1,real **mat,real cutoff)
{
  t_clustid *c;
  t_dist    *d;
  int    i,j,k,nn,cid;
  
  /* First we sort the entries in the RMS matrix */
  nn = ((n1-1)*n1)/2;
  snew(d,nn);
  for(i=k=0; (i<n1); i++)
    for(j=i+1; (j<n1); j++,k++) {
      d[k].i    = i;
      d[k].j    = j;
      d[k].dist = mat[i][j];
    }
  assert(k == nn);
  qsort(d,nn,sizeof(d[0]),dcomp);
  
  /* Now we make a cluster index for all of the conformations */
  snew(c,n1);
  for(i=0; (i<n1); i++) {
    c[i].conf = i;
    c[i].clust = i;
  }
    
  /* Now we check the closest structures, and equalize their cluster
   * numbers 
   */
  for(k=0; (k<nn) && (d[k].dist < cutoff); k++) {
    i = d[k].i;
    j = d[k].j;
    c[j].clust = c[i].clust;
  }
  printf("Sorting and renumbering clusters\n");
  /* Sort on cluster number */
  qsort(c,n1,sizeof(c[0]),ccomp);

  /* Renumber clusters */
  cid = c[0].clust;
  for(k=1; (k<n1); k++) {
    if (c[k].clust != cid) {
      cid ++;
      c[k].clust = cid;
    }
    else
      c[k].clust = c[k-1].clust;
  }
  for(k=0; (k<n1); k++)
    printf("Cluster index for conformation %d: %d\n",c[k].conf,c[k].clust);

  mcpy = mk_matrix(n1,FALSE);
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++)
      mcpy[c[i].conf][c[j].conf] = mat[i][j];
  }
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++)
      mat[i][j] = mcpy[i][j];
  }
    
  sfree(c);
  sfree(d);
}

rvec **read_whole_trj(char *fn,int isize,atom_id index[],int skip,int *nframe)
{
  rvec   **xx,*x;
  matrix box;
  real   t;
  int    i,i0,j,nf;
  int    natom,nbytes,status;
  
  nf     = *nframe;
  nbytes = nf*isize*sizeof(x[0]);
  printf("Total memory required to read the whole trajectory: %d bytes\n",
	 nbytes);
  snew(xx,nf);
  for(j=0; (j<nf); j++)
    snew(xx[j],isize);
  natom = read_first_x(&status,fn,&t,&x,box);
  i  = 0;
  i0 = 0;
  do {
    if ((i % skip) == 0) {
      /* Store only the interesting atoms */
      for(j=0; (j<isize); j++) 
	copy_rvec(x[index[j]],xx[i0][j]);
      i0 ++;
    }
    i++;
  } while ((read_next_x(status,&t,natom,x,box)) && (i0 < nf));
  fprintf(stderr,"\n");
  printf("Read %d frames from trajectory %s\n",i0,fn);
  *nframe = i0;
  sfree(x);
  
  return xx;
}

void rms_dist(char *fn,int nf,real **mat,real maxrms)
{
  FILE   *fp;
  int    i,j,*histo;
  real   fac;
  
  fac = 100/maxrms;
  snew(histo,101);
  for(i=0; (i<nf); i++) 
    for(j=i+1; (j<nf); j++)
      histo[(int)(fac*mat[i][j])]++;
      
  fp = xvgropen(fn,"RMS Distribution","RMS (nm)","a.u.");
  for(i=0; (i<101); i++)
    fprintf(fp,"%10g  %10d\n",i/fac,histo[i]);
  fclose(fp);
  sfree(histo);
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
  int          natom,i1,i2,i1s,i2s,i,teller,nf,nrms;
  real         t,t1,t2;

  matrix       box;
  rvec         *x1,*x2,**xx;
  char         *fn;
  double       sum;
  real         **rms,maxrms,*resnr;
  real         *rms0,*eigval;
  t_topology   top;
  bool         bSameF;
  t_rgb        rlo,rhi;
  
  int      status1,status2,isize,nbytes;
  atom_id  *index;
  char     *grpname;
  real     **d1,**d2;
  char     buf[255];
  
  static int  nlevels=40,nframes=100,skip=0;
  static real scalemax=-1.0,rmscut=0.1;
  static bool bEigen=FALSE,bMem=TRUE,bLink=FALSE;
  t_pargs pa[] = {
    { "-nlevels",   FALSE, etINT,  &nlevels,
      "Discretize rms in # levels" },
    { "-rmscut", FALSE, etREAL, &rmscut,
      "RMS cut-off (nm) for two structures to be similar" },
    { "-max",   FALSE, etREAL, &scalemax,
      "Maximum level in matrices" },
    { "-nframes", FALSE, etINT, &nframes,
      "Max number of frames" },
    { "-skip",  FALSE, etINT, &skip,
      "Skip N frames after each one used." },
    { "-eigen", FALSE, etBOOL, &bEigen,
      "Diagonalize the rms matrix instead of sorting it" },
    { "-link",  FALSE, etBOOL, &bLink,
      "Use linkage algorithm for clustering" },
    { "-mem",   FALSE, etBOOL, &bMem,
      "Read the whole trajectory in memory (can be large...)" }
  };
  t_filenm fnm[] = {
    { efTRX, "-f",    NULL,       ffREAD  },
    { efTPS, "-s",    NULL,       ffREAD  },
    { efNDX, NULL,    NULL,       ffOPTRD },
    { efXPM, "-o",   "rmsd",      ffWRITE },
    { efXPM, "-os",  "rmsd-mat",  ffWRITE },
    { efXVG, "-ev",  "rms-eig",   ffWRITE },
    { efXVG, "-dist","rms-dist",  ffWRITE }
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
  
  rms = mk_matrix(nf,bEigen);
  if (bEigen)
    rms0 = rms[0];
  else
    rms0 = NULL;
    
  snew(resnr,nf);
  for(i=0; (i<nf); i++) 
    resnr[i] = i;

  /*set box type*/
  init_pbc(box,FALSE);

  /* Loop over first coordinate file */
  fn = opt2fn("-f",NFILE,fnm);
  
  maxrms = 0;
  sum    = 0;
  if (bMem) {
    xx = read_whole_trj(fn,isize,index,skip,&nf);
    fprintf(stderr,"Computing %dX%d RMS matrix\n",nf,nf);
    nrms = (nf*(nf-1))/2;
    for(i1=0; (i1<nf); i1++) {
      calc_dist(isize,xx[i1],box,d1);
      for(i2=i1+1; (i2<nf); i2++) {
	calc_dist(isize,xx[i2],box,d2);
	rms[i1][i2] = rms[i2][i1] = rms_diff(isize,d1,d2);
	
	/* Store maximum value for plotting */
	maxrms = max(maxrms,rms[i1][i2]);
	sum += rms[i1][i2];
      }
      nrms -= (nf-i1-1);
      fprintf(stderr,"\r# RMS calculations left: %d   ",nrms);
    }
    fprintf(stderr,"\nDone\n");
  }
  else {
    /* Read first frames of both files */  
    natom = read_first_x(&status1,fn,&t1,&x1,box);
    natom = read_first_x(&status2,fn,&t1,&x2,box);
    
    /* Consistency checks */
    if (isize > natom) 
      fatal_error(0,"Index (from index file or topology)"
		  " does not match trajectory %s",fn);
    
    /* Initiate locals */	
    i1     = 0;
    i1s    = 0;
    i2s    = 0;
    do { 
      /* Loop over first file */
      if ((i1 % skip) == 0) {
	/* Print counter */
	fprintf(stderr,"\ri1 = %d                               \n",i1);
	
	/* Compute distance matrix for first structure */
	calc_dist_ind(isize,index,x1,box,d1);
	
	/* Find the starting frame for the second file */
	i2 = i1+1;
	
	/* Rewind second traj, and skip to first interesting frame */
	rewind_trj(status2);
	for(i=0; (i<=i2); i++)
	  if (!read_next_x(status2,&t2,natom,x2,box))
	    break;
	
	do {
	  /* Loop over second file */
	  if ((i2 % skip) == 0) {
	    /* Compute distance matrix for second structure */  
	    calc_dist_ind(isize,index,x2,box,d2);
	    
	    /* Compute RMS between two distance matrices */
	    rms[i1s][i2s] = rms[i2s][i1s] = rms_diff(isize,d1,d2);
	    
	    /* Store maximum value for plotting */
	    maxrms = max(maxrms,rms[i1s][i2s]);
	    sum += rms[i1s][i2s];
	    i2s++;
	  }
	  i2++;
	} while ((read_next_x(status2,&t2,natom,x2,box)) && (i2s < nf));
	i1s++;
      }
      i1++;
    } while ((read_next_x(status1,&t1,natom,x1,box)) && (i1s < nf));
    fprintf(stderr,"\n");
    nf = i1;
  }
  fprintf(stderr,"The maximum RMS was %g nm\n"
	  "Average RMS was %g\n"
	  "number of frames for matrix %d\n",
	  maxrms,2*sum/(nf*(nf-1)),nf);
  
  /* Plot the rms distribution */
  rms_dist(opt2fn("-dist",NFILE,fnm),nf,rms,maxrms);
  
  /* Set colors for plotting: white = zero RMS, black = maximum */
  rlo.r=1.0, rlo.g=1.0, rlo.b=0.0;
  rhi.r=0.0, rhi.g=0.0, rhi.b=1.0;
  
  /* Write out plot file with RMS matrix */
  write_xpm(opt2FILE("-o",NFILE,fnm,"w"),"RMS",
	    "RMS (n)","Confs 1","Confs 2",
	    nf,nf,resnr,resnr,rms,0.0,maxrms,rlo,rhi,&nlevels);

  if (bLink) 
    /* Now sort the matrix and write it out again */
    gather(nf,rms,rmscut);
  else {
    /* Do a diagonalization */
    snew(eigval,nf);
    /*for(i1=0; (i1<nf); i1++)
      for(i2=0; (i2<nf); i2++)
      rms[i1][i2] = maxrms - rms[i1][i2];*/
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
