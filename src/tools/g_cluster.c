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
#include "random.h"
#include "pbc.h"
#include "xvgr.h"
#include "futil.h"
#include "matio.h"
#include "ql77.h"
#include "cmat.h"
#include "do_fit.h"
#include "trnio.h"

typedef struct {
  int ncl;
  int *cl;
} t_clusters;

void pr_energy(FILE *fp,real e)
{
  fprintf(fp,"Energy: %8.4f\n",e);  
}

void cp_index(int nn,int from[],int to[])
{
  int i;
  
  for(i=0; (i<nn); i++)
    to[i]=from[i];
}

void mc_optimize(FILE *log,t_mat *m,int maxiter,int *seed,real kT)
{
  real e[2],ei,ej,efac;
  int  *low_index;
  int  cur=0;
#define next (1-cur)
  int  i,isw,jsw,iisw,jjsw,nn;
  
  fprintf(stderr,"Doing Monte Carlo clustering\n");
  nn = m->nn;
  snew(low_index,nn);
  cp_index(nn,m->m_ind,low_index);
  if (getenv("TESTMC")) {
    e[cur] = mat_energy(m);
    pr_energy(log,e[cur]);
    fprintf(log,"Doing 1000 random swaps\n");
    for(i=0; (i<1000); i++) {
      do {
	isw = nn*rando(seed);
	jsw = nn*rando(seed);
      } while ((isw == jsw) || (isw >= nn) || (jsw >= nn));
      iisw = m->m_ind[isw];
      jjsw = m->m_ind[jsw];
      m->m_ind[isw] = jjsw;
      m->m_ind[jsw] = iisw;
    }
  }
  e[cur] = mat_energy(m);
  pr_energy(log,e[cur]);
  for(i=0; (i<maxiter); i++) {
    do {
      isw = nn*rando(seed);
      jsw = nn*rando(seed);
    } while ((isw == jsw) || (isw >= nn) || (jsw >= nn));
    
    iisw = m->m_ind[isw];
    jjsw = m->m_ind[jsw];
    ei   = row_energy(nn,iisw,m->mat[jsw],m->m_ind);
    ej   = row_energy(nn,jjsw,m->mat[isw],m->m_ind);
    
    e[next] = e[cur] + (ei+ej-EROW(m,isw)-EROW(m,jsw))/nn;

    efac = kT ? exp((e[next]-e[cur])/kT) : -1;
    if ((e[next] > e[cur]) || (efac > rando(seed))) {
      
      if (e[next] > e[cur])
	cp_index(nn,m->m_ind,low_index);
      else
	fprintf(log,"Taking uphill step\n");
	
      /* Now swapping rows */
      m->m_ind[isw] = jjsw;
      m->m_ind[jsw] = iisw;
      EROW(m,isw)   = ei;
      EROW(m,jsw)   = ej;
      cur           = next;
      fprintf(log,"Iter: %d Swapped %4d and %4d (now %g)",
	      i,isw,jsw,mat_energy(m));
      pr_energy(log,e[cur]);
    }
  }
  /* Now restore the highest energy index */
  cp_index(nn,low_index,m->m_ind);
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

static real rms_dist(int isize,real **d,real **d_r)
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
  real *ra,*rb;
  
  ra = (real *)a;
  rb = (real *)b;
  
  if ((*ra) < (*rb))
    return -1;
  else if ((*ra) > (*rb))
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
  /* done_mat(n1,mcpy); 
     sfree(mcpy);*/
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

void gather(FILE *log,t_mat *m,real cutoff,t_clusters *clust)
{
  t_clustid *c;
  t_dist    *d;
  int       i,j,k,nn,cid,n1,diff;
  bool      bChange;
  
  /* First we sort the entries in the RMSD matrix */
  n1 = m->nn;
  nn = ((n1-1)*n1)/2;
  snew(d,nn);
  for(i=k=0; (i<n1); i++)
    for(j=i+1; (j<n1); j++,k++) {
      d[k].i    = i;
      d[k].j    = j;
      d[k].dist = m->mat[i][j];
    }
  assert(k == nn);
  qsort(d,nn,sizeof(d[0]),dcomp);
  
  /* Now we make a cluster index for all of the conformations */
  c = new_clustid(n1);
    
  /* Now we check the closest structures, and equalize their cluster
   * numbers 
   */
  fprintf(stderr,"Linking structures ");
  do {
    fprintf(stderr,"*");
    bChange=FALSE;
    for(k=0; (k<nn) && (d[k].dist < cutoff); k++) {
      diff = c[d[k].j].clust - c[d[k].i].clust;
      if (diff) {
	bChange = TRUE;
	if (diff > 0)
	  c[d[k].j].clust = c[d[k].i].clust;
	else
	  c[d[k].i].clust = c[d[k].j].clust;
      }
    }
  } while (bChange);
  fprintf(stderr,"\nSorting and renumbering clusters\n");
  /* Sort on cluster number */
  qsort(c,n1,sizeof(c[0]),ccomp);

  /* Renumber clusters */
  cid = 1;
  for(k=1; k<n1; k++) {
    if (c[k].clust != c[k-1].clust) {
      c[k-1].clust = cid;
      cid ++;
    } else
      c[k-1].clust = cid;
  }
  c[k-1].clust = cid;
  if (debug)
    for(k=0; (k<n1); k++)
      fprintf(debug,"Cluster index for conformation %d: %d\n",
	      c[k].conf,c[k].clust);
  clust->ncl = cid;
  for(k=0; k<n1; k++)
    clust->cl[c[k].conf] = c[k].clust;

  mcpy = mk_matrix(n1,FALSE);
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++)
      mcpy[c[i].conf][c[j].conf] = m->mat[i][j];
  }
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++)
      m->mat[i][j] = mcpy[i][j];
  }
  done_matrix(n1,&mcpy);
  
  sfree(c);
  sfree(d);
}

bool jp_same(int **nnb,int i,int j,int P)
{
  bool bIn;
  int  k,ii,jj,pp;

  bIn = FALSE;
  for(k=0; nnb[i][k]>=0; k++)
    bIn = bIn || (nnb[i][k] == j);
  if (!bIn)
    return FALSE;

  bIn = FALSE;
  for(k=0; nnb[j][k]>=0; k++)
    bIn = bIn || (nnb[j][k] == i);
  if (!bIn)
    return FALSE;

  pp=0;
  for(ii=0; nnb[i][ii]>=0; ii++)
    for(jj=0; nnb[j][jj]>=0; jj++)
      if ((nnb[i][ii] == nnb[j][jj]) && (nnb[i][ii] != -1))
	pp++;

  return (pp >= P);
}

static void jarvis_patrick(FILE *log,int n1,real **mat,int M,int P,
			   real rmsdcut,t_clusters *clust)
{
  t_dist    *tmp;
  t_clustid *c;
  int       **nnb;
  int       i,j,k,cid,diff,max;
  bool      bChange;

  /* First we sort the entries in the RMSD matrix row by row.
   * This gives us the nearest neighbor list.
   */
  snew(nnb,n1);
  snew(tmp,n1);
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++) {
      tmp[j].j    = j;
      tmp[j].dist = mat[i][j];
    }
    qsort(tmp,n1,sizeof(tmp[0]),dcomp);
    if (M>0) {
      /* Put the M nearest neighbors in the list */
      snew(nnb[i],M+1);
      for(j=k=0; (k<M) && (j<n1); j++)
	if (tmp[j].j  != i) {
	  nnb[i][k]  = tmp[j].j;
	  k++;
	}
      nnb[i][M] = -1;
    } else {
      /* Put all neighbors nearer than rmsdcut in the list */
      max=0;
      k=0;
      for(j=0; (j<n1) && (mat[i][tmp[j].j] < rmsdcut); j++)
	if (tmp[j].j != i) {
	  if (k >= max) {
	    max += 10;
	    srenew(nnb[i],max);
	  }
	  nnb[i][k] = tmp[j].j;
	  k++;
	}
      if (k == max)
	srenew(nnb[i],max+1);
      nnb[i][k] = -1;
    }
  }
  sfree(tmp);
  if (log) {
    fprintf(log,"Nearest neighborlist. M = %d, P = %d\n",M,P);
    for(i=0; (i<n1); i++) {
      fprintf(log,"i: %5d nbs:",i);
      for(j=0; nnb[i][j]>=0; j++)
	fprintf(log,"  %5d[%6.3f]",nnb[i][j],mat[i][nnb[i][j]]);
      fprintf(log,"\n");
    }
  }

  c = new_clustid(n1);
  fprintf(stderr,"Linking structures ");
  /* Use mcpy for temporary storage of booleans */
  mcpy = mk_matrix(n1,FALSE);
  for(i=0; i<n1; i++)
    for(j=i+1; j<n1; j++)
      mcpy[i][j] = jp_same(nnb,i,j,P);
  do {
    fprintf(stderr,"*");
    bChange=FALSE;
    for(i=0; i<n1; i++) {
      for(j=i+1; j<n1; j++)
	if (mcpy[i][j]) {
	  diff = c[j].clust - c[i].clust;
	  if (diff) {
	    bChange = TRUE;
	    if (diff > 0)
	      c[j].clust = c[i].clust;
	    else
	      c[i].clust = c[j].clust;
	  }
	}
    }
  } while (bChange);
  
  fprintf(stderr,"\nSorting and renumbering clusters\n");
  /* Sort on cluster number */
  qsort(c,n1,sizeof(c[0]),ccomp);

  /* Renumber clusters */
  cid = 1;
  for(k=1; k<n1; k++) {
    if (c[k].clust != c[k-1].clust) {
      c[k-1].clust = cid;
      cid ++;
    } else
      c[k-1].clust = cid;
  }
  c[k-1].clust = cid;
  clust->ncl = cid;
  for(k=0; k<n1; k++)
    clust->cl[c[k].conf] = c[k].clust;
  if (debug)
    for(k=0; (k<n1); k++)
      fprintf(debug,"Cluster index for conformation %d: %d\n",
	      c[k].conf,c[k].clust);

  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++)
      mcpy[c[i].conf][c[j].conf] = mat[i][j];
  }
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++)
      mat[i][j] = mcpy[i][j];
  }
    
  sfree(c);
  for(i=0; (i<n1); i++)
    sfree(nnb[i]);
  sfree(nnb);
}

rvec **read_whole_trj(char *fn,int isize,atom_id index[],int skip,int *nframe,
		      int *natom, real **time)
{
  rvec   **xx,*x;
  matrix box;
  real   t;
  int    i,i0,j,max_nf;
  int    nbytes,status;
  
  max_nf = 0;
  xx     = NULL;
  *time  = NULL;
  *natom = read_first_x(&status,fn,&t,&x,box);
  i  = 0;
  i0 = 0;
  do {
    if (i0 >= max_nf) {
      max_nf += 10;
      srenew(xx,max_nf);
      srenew(*time,max_nf);
    }
    if ((i % skip) == 0) {
      snew(xx[i0],isize);
      /* Store only the interesting atoms */
      for(j=0; (j<isize); j++) 
	copy_rvec(x[index[j]],xx[i0][j]);
      (*time)[i0] = t;
	i0 ++;
    }
    i++;
  } while (read_next_x(status,&t,*natom,x,box));
  fprintf(stderr,"Allocated %d bytes for frames\n",max_nf*isize*sizeof(**xx));
  fprintf(stderr,"Read %d frames from trajectory %s\n",i0,fn);
  *nframe = i0;
  sfree(x);
  
  return xx;
}

static void plot_clusters(int nf,real **mat,real val,t_clusters *clust)
{
  int i,j;
  
  for(i=0; i<nf; i++)
    for(j=0; j<i; j++)
      if (clust->cl[i] == clust->cl[j])
	mat[i][j] = val;
      else
	mat[i][j] = 0;
}
 
static void analyze_clusters(int nf,t_clusters *clust,real **rmsd,
			     t_atoms *atoms,int natom,
			     int isize,atom_id *index,atom_id *alli,
			     rvec **xx,real *time,real *mass,
			     rvec *xtps,char *trxfn,bool bAv,FILE *log)
{
  FILE *fp;
  char buf[STRLEN],buf1[20],buf2[20];
  int  trxout;
  int  i,i1,i2,cl,nstr,*structure,first,midstr;
  real r,clrmsd,midrmsd;
  rvec *xtpsi,*xav,*xnatom;
  matrix zerobox;

  clear_mat(zerobox);

  sprintf(buf,"\nFound %d clusters\n"
	  "\nWriting %s structure for each cluster to %s\n",
	  clust->ncl,bAv ? "average" : "middle",trxfn);
  fprintf(stderr,buf);
  fprintf(log,buf);
  
  /* Prepare a reference structure for the orientation of the clusters  */
  snew(xtpsi,isize);
  for(i=0; i<isize; i++)
    copy_rvec(xtps[index[i]],xtpsi[i]);
  reset_x(isize,alli,isize,alli,xtpsi,mass);
  trxout = open_trx(trxfn,"w");
  /* Calculate the average structure in each cluster,               *
   * all structures are fitted to the first struture of the cluster */
  snew(xav,isize);
  snew(xnatom,natom);
  snew(structure,nf);
  sprintf(buf,"\n%3s %3s %4s %6s %4s\n","cl.","#st","rmsd","middle","rmsd");
  fprintf(stderr,buf);
  fprintf(log,buf);
  for(cl=1; cl<=clust->ncl; cl++) {
    for(i=0; i<natom;i++)
      clear_rvec(xav[i]);
    nstr=0;
    for(i1=0; i1<nf; i1++)
      if (clust->cl[i1] == cl) {
	structure[nstr] = i1;
	nstr++;
	if (bAv) {
	  reset_x(isize,alli,isize,alli,xx[i1],mass);
	  if (nstr == 1)
	    first = i1;
	  else
	    do_fit(isize,mass,xx[first],xx[i1]);
	  for(i=0; i<isize; i++)
	    rvec_inc(xav[i],xx[i1][i]);
	}
      }
    clrmsd = 0;
    midstr = 0;
    midrmsd = 10000;
    for(i1=0; i1<nstr; i1++) {
      r = 0;
      if (nstr > 1) {
	for(i=0; i<nstr; i++)
	  r += rmsd[structure[i1]][structure[i]];
	r /= (nstr - 1);
      }
      if (r < midrmsd) {
	midstr = structure[i1];
	midrmsd = r;
      }
      clrmsd += r;
    }
    clrmsd /= nstr;
    if (nstr > 1) {
      sprintf(buf1,"%5.3f",clrmsd);
      if (buf1[0] == '0')
	buf1[0] = ' ';
      sprintf(buf2,"%5.3f",midrmsd);
      if (buf2[0] == '0')
	buf2[0] = ' ';
    } else {
      sprintf(buf1,"%5s","");
      sprintf(buf2,"%5s","");
    }
    sprintf(buf,"%3d %3d%s %6.f%s:",
	    cl,nstr,buf1,time[midstr],buf2);
    fprintf(stderr,buf);
    fprintf(log,buf);
    for(i=0; i<nstr; i++) {
      if ((i % 8 == 0) && i)
	sprintf(buf,"\n%25s","");
      else
	buf[0] = '\0';
      i1 = structure[i];
      fprintf(stderr,"%s %6.f",buf,time[i1]);
      fprintf(log,"%s %6.f",buf,time[i1]);
    }
    fprintf(stderr,"\n");
    fprintf(log,"\n");
    
    /* Dump the average structure for this cluster */
    if (bAv) {
      r = 1.0/nstr;
      for(i=0; i<isize; i++)
	svmul(r,xav[i],xav[i]);
    } else {
      for(i=0; i<isize; i++)
	copy_rvec(xx[midstr][i],xav[i]);
      reset_x(isize,alli,isize,alli,xav,mass);
    }
    do_fit(isize,mass,xtpsi,xav);
    for(i=0; i<isize; i++)
      copy_rvec(xav[i],xnatom[index[i]]);
    r = cl;
    write_trx(trxout,isize,index,atoms,0,r,zerobox,xnatom,NULL);
  }
  close_trx(trxout);
  
  sfree(xtpsi);
  sfree(xav);
  sfree(xnatom);
  sfree(structure);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_cluster can clusters structures with several different methods.",
    "RMS deviation after fitting or RMS deviation of atom-pair distances",
    "can be used to define the distance between structures.[PAR]",
    "full linkage: add a structure to a cluster when its distance to any",
    "element of the cluster is less than [TT]cutoff[tt].[PAR]",
    "Jarvis Patrick: add a structure to a cluster when this structure",
    "and a structure in the cluster have each other as neighbors and",
    "they have a least [TT]P[tt] neighbors in common. The neighbors",
    "of a structure are the M closest structures or all structures within",
    "[TT]cutoff[tt].[PAR]",
    "Monte Carlo: reorder the RMSD matrix using Monte Carlo.[PAR]",
    "diagonalization: diagonalize the RMSD matrix.[PAR]"
    "When unique cluster assignments can be determined (full linkage and",
    "Jarvis Patrick), the structure with the smallest average distance to",
    "the others or the average structure for each cluster will be written",
    "to a trajectory file."
  };
  
  FILE         *fp,*log;
  int          natom,i1,i2,i1s,i2s,i,teller,nf,nrms;
  real         t,t1,t2;

  matrix       box;
  rvec         *xtps,*x1,*x2,**xx;
  char         *fn;
  t_clusters   clust;
  t_mat        *rms;
  real         **rmsmat;
  real         *resnr,*eigval;
  t_topology   top;
  bool         bSameF;
  t_rgb        rlo,rhi;
  
  int      status1,status2,isize,nbytes;
  atom_id  *index,*alli;
  char     *grpname;
  real     rmsd,**d1,**d2,*time,*mass;
  char     buf[STRLEN];
  bool     bAnalyze;
  
  static char *method[] = { NULL, "linkage", "jarvis-patrick","monte-carlo",
			    "diagonalization", NULL };
  static int  nlevels=40,skip=1;
  static real scalemax=-1.0,rmsdcut=0.1;
  static bool bRMSdist=FALSE,bBinary=FALSE,bAv=FALSE;
  static int  niter=10000,seed=1993;
  static real kT=1e-3;
  static int  M=10,P=3;
  t_pargs pa[] = {
    { "-dista",FALSE, etBOOL, &bRMSdist,
      "Use RMSD of distances instead of RMS deviation" },
    { "-nlevels",   FALSE, etINT,  &nlevels,
      "Discretize RMSD matrix in # levels" },
    { "-cutoff", FALSE, etREAL, &rmsdcut,
      "RMSD cut-off (nm) for two structures to be similar" },
    { "-max",   FALSE, etREAL, &scalemax,
      "Maximum level in matrices" },
    { "-skip",  FALSE, etINT, &skip,
      "Only analyze every nr-th frame" },
    { "-av",  FALSE, etBOOL, &bAv,
      "Write average iso middle structure for each cluster" },
    { "-method",FALSE, etENUM, method,
      "Method for cluster determination" },
    { "-binary",FALSE, etBOOL, &bBinary,
      "Treat the RMSD matrix as consisting of 0 and 1, where the cut-off is given by the value you put in for rms-cut" },
    { "-M",     FALSE, etINT,  &M,
      "Number of nearest neighbours considered for Jarvis-Patrick algorithm, 0 is use cutoff" },
    { "-P",     FALSE, etINT,  &P,
      "Number of identical nearest neighbors required to form a cluster" },
    { "-seed",  FALSE, etINT,  &seed,
      "Random number seed for Monte Carlo clustering algorithm" },
    { "-niter", FALSE, etINT,  &niter,
      "Number of iterations for MC" },
    { "-kT",    FALSE, etREAL, &kT,
      "Boltzmann weighting factor for Monte Carlo optimization (zero turns off uphill steps)" }
  };
  t_filenm fnm[] = {
    { efTRX, "-f",     NULL,        ffREAD  },
    { efTPS, "-s",     NULL,        ffREAD  },
    { efNDX, NULL,     NULL,        ffOPTRD },
    { efXPM, "-o",    "rmsd-clust", ffWRITE },
    { efLOG, "-g",    "cluster",    ffWRITE },
    { efXVG, "-dist", "rmsd-dist",  ffWRITE },
    { efXVG, "-ev",   "rmsd-eig",   ffOPTWR },
    { efTRX, "-cl",   "clusters.pdb", ffOPTWR }
  };
#define NFILE asize(fnm)

#define m_linkage        (method[0][0] == 'l')
#define m_jarvis_patrick (method[0][0] == 'j')
#define m_monte_carlo    (method[0][0] == 'm')
#define m_diagonalize    (method[0][0] == 'd')

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if ((M<0) || (M == 1))
    fatal_error(0,"M (%d) must be 0 or larger than 1",M);
  if (M < 2)
    fprintf(stderr,"Will use RMSD cutoff (%g) for determining the neighbors\n",
	    rmsdcut);
  else if (P >= M)
    fatal_error(0,"Number of neighbors required (P) must be less than M");
  if (skip < 1)
    fatal_error(0,"skip (%d) should be >= 1",skip);

  bAnalyze = (m_linkage || m_jarvis_patrick);

  /* don't read mass-database as masses (and top) are not used */
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),buf,&top,&xtps,NULL,box,bAnalyze);
  
  fprintf(stderr,"\nSelect group for RMSD calculation:\n");
  get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    1,&isize,&index,&grpname);

  /* Initiate arrays */  
  snew(d1,isize);
  snew(d2,isize);
  for(i=0; (i<isize); i++) {
    snew(d1[i],isize);
    snew(d2[i],isize);
  }
    
  /*set box type*/
  init_pbc(box,FALSE);

  /* Loop over first coordinate file */
  fn = opt2fn("-f",NFILE,fnm);
  
  xx = read_whole_trj(fn,isize,index,skip,&nf,&natom,&time);
  if (!bRMSdist || bAnalyze) {
    /* Center all frames on zero */
    snew(alli,isize);
    for(i=0; i<isize; i++)
      alli[i] = i;
    snew(mass,isize);
    for(i=0; i<isize; i++)
      mass[i] = top.atoms.atom[index[i]].m;
    for(i=0; i<nf; i++)
      reset_x(isize,alli,isize,alli,xx[i],mass);
  }
  rms = init_mat(nf,m_diagonalize);
  nrms = (nf*(nf-1))/2;
  if (!bRMSdist) {
    fprintf(stderr,"Computing %dx%d RMS deviation matrix\n",nf,nf);
    snew(x1,isize);
    for(i1=0; (i1<nf); i1++) {
      for(i2=i1+1; (i2<nf); i2++) {
	for(i=0; i<isize; i++)
	  copy_rvec(xx[i1][i],x1[i]);
	do_fit(isize,mass,xx[i2],x1);
	rmsd = rmsdev(isize,mass,xx[i2],x1);
	set_mat_entry(rms,i1,i2,rmsd);
      }
      nrms -= (nf-i1-1);
      fprintf(stderr,"\r# RMSD calculations left: %d   ",nrms);
    }
  } else {
    fprintf(stderr,"Computing %dx%d RMS distance deviation matrix\n",nf,nf);
    for(i1=0; (i1<nf); i1++) {
      calc_dist(isize,xx[i1],box,d1);
      for(i2=i1+1; (i2<nf); i2++) {
	calc_dist(isize,xx[i2],box,d2);
	set_mat_entry(rms,i1,i2,rms_dist(isize,d1,d2));
      }
      nrms -= (nf-i1-1);
      fprintf(stderr,"\r# RMSD calculations left: %d   ",nrms);
    }
  }
  fprintf(stderr,"\n\n");

  /* Open log file */
  log = ftp2FILE(efLOG,NFILE,fnm,"w");
  
  sprintf(buf,"The maximum RMSD is %g nm\n"
	  "Average RMSD is %g\n"
	  "Number of structures for matrix %d\n"
	  "Energy of the matrix is %g nm\n\n",
	  rms->maxrms,2*rms->sumrms/(nf*(nf-1)),nf,mat_energy(rms));
  fprintf(stderr,buf);
  fprintf(log,buf);
  fprintf(stderr,"Using %s for clustering\n",method[0]);
  fprintf(log,"Using %s for clustering\n",method[0]);

  /* Plot the rmsd distribution */
  rmsd_distribution(opt2fn("-dist",NFILE,fnm),rms);
 
  snew(rmsmat,nf);
  for(i1=0; (i1 < nf); i1++) {
    snew(rmsmat[i1],nf);
    for(i2=0; (i2 < nf); i2++)
      rmsmat[i1][i2] = rms->mat[i1][i2];
  }

  if (bBinary) {
    for(i1=0; (i1 < nf); i1++) 
      for(i2=0; (i2 < nf); i2++)
	if (rms->mat[i1][i2] < rmsdcut)
	  rms->mat[i1][i2] = 0;
	else
	  rms->mat[i1][i2] = 1;
  }

  snew(clust.cl,nf);
  if (m_linkage) 
    /* Now sort the matrix and write it out again */
    gather(log,rms,rmsdcut,&clust);
  else if (m_diagonalize) {
    /* Do a diagonalization */
    snew(eigval,nf);
    /*for(i1=0; (i1<nf); i1++)
      for(i2=0; (i2<nf); i2++)
      rms[i1][i2] = maxrms - rms[i1][i2];*/
    ql77(nf,rms->mat[0],eigval); 
    fp = xvgropen(opt2fn("-ev",NFILE,fnm),"Eigenvalues","index","value");
    for(i=0; (i<nf); i++)
      fprintf(fp,"%10d  %10g\n",i,eigval[i]);
    ffclose(fp);
    xvgr_file(opt2fn("-ev",NFILE,fnm),NULL);
  }
  else if (m_monte_carlo)
    mc_optimize(log,rms,niter,&seed,kT);
  else if (m_jarvis_patrick)
    jarvis_patrick(log,rms->nn,rms->mat,M,P,rmsdcut,&clust);
  else
    fatal_error(0,"unknown method \"%s\"",method[0]);


  if (m_monte_carlo || m_diagonalize)
  fprintf(stderr,"Energy of the matrix after clustering is %g nm\n",
	  mat_energy(rms));
  swap_mat(rms);
  reset_index(rms);

  /* Write the rmsd-matrix to the upper-left part of the matrix */
  for(i1=0; (i1 < nf); i1++)
    for(i2=i1; (i2 < nf); i2++)
      rms->mat[i1][i2] = rmsmat[i1][i2];
  
  if (bAnalyze) {
    plot_clusters(nf,rms->mat,rms->maxrms,&clust);
    analyze_clusters(nf,&clust,rmsmat,&(top.atoms),natom,isize,index,alli,
		     xx,time,mass,xtps,opt2fn("-cl",NFILE,fnm),bAv,log);
  }
  ffclose(log);
  
  if (bBinary && !bAnalyze)
    /* Make the clustering visible */
    for(i2=0; (i2 < nf); i2++)
       for(i1=i2+1; (i1 < nf); i1++)
	 if (rms->mat[i1][i2])
	   rms->mat[i1][i2] = rms->maxrms;

  /* Set colors for plotting: white = zero RMS, black = maximum */
  rlo.r=1.0, rlo.g=1.0, rlo.b=1.0;
  rhi.r=0.0, rhi.g=0.0, rhi.b=0.0;
  
  fp = opt2FILE("-o",NFILE,fnm,"w");
  write_xpm(fp,"RMSD","RMSD (nm)","Time (ps)","Time (ps)",
	    nf,nf,time,time,rms->mat,0.0,rms->maxrms,rlo,rhi,&nlevels);
  ffclose(fp);
  xv_file(opt2fn("-o",NFILE,fnm),NULL);	

  /* Thank the user for her patience */  
  thanx(stdout);
  
  return 0;
}
