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
static char *SRCID_g_cluster_c = "$Id$";

#include <math.h>
#include <string.h>
#include <ctype.h>
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
#include "viewit.h"

typedef struct {
  int ncl;
  int *cl;
} t_clusters;

typedef struct {
  int nr;
  int *nb;
} t_nnb;
  
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
  
  fprintf(stderr,"\nDoing Monte Carlo clustering\n");
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

static void calc_dist(int nind,rvec x[],real **d)
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

static int rms_dist_comp(const void *a,const void *b)
{
  t_dist *da,*db;
  
  da = (t_dist *)a;
  db = (t_dist *)b;
  
  if (da->dist - db->dist < 0)
    return -1;
  else if (da->dist - db->dist > 0)
    return 1;
  return 0;
}

static int clust_id_comp(const void *a,const void *b)
{
  t_clustid *da,*db;
  
  da = (t_clustid *)a;
  db = (t_clustid *)b;
  
  return da->clust - db->clust;
}

static int nrnb_comp(const void *a, const void *b)
{
  t_nnb *da, *db;
  
  da = (t_nnb *)a;
  db = (t_nnb *)b;
  
  /* return the b-a, we want highest first */
  return db->nr - da->nr;
}

void gather(t_mat *m,real cutoff,t_clusters *clust)
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
  qsort(d,nn,sizeof(d[0]),rms_dist_comp);
  
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
  qsort(c,n1,sizeof(c[0]),clust_id_comp);

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

static void jarvis_patrick(int n1,real **mat,int M,int P,
			   real rmsdcut,t_clusters *clust)
{
  t_dist    *row;
  t_clustid *c;
  int       **nnb;
  int       i,j,k,cid,diff,max;
  bool      bChange;

  if (rmsdcut < 0)
    rmsdcut = 10000;

  /* First we sort the entries in the RMSD matrix row by row.
   * This gives us the nearest neighbor list.
   */
  snew(nnb,n1);
  snew(row,n1);
  for(i=0; (i<n1); i++) {
    for(j=0; (j<n1); j++) {
      row[j].j    = j;
      row[j].dist = mat[i][j];
    }
    qsort(row,n1,sizeof(row[0]),rms_dist_comp);
    if (M>0) {
      /* Put the M nearest neighbors in the list */
      snew(nnb[i],M+1);
      for(j=k=0; (k<M) && (j<n1) && (mat[i][row[j].j] < rmsdcut); j++)
	if (row[j].j  != i) {
	  nnb[i][k]  = row[j].j;
	  k++;
	}
      nnb[i][k] = -1;
    } else {
      /* Put all neighbors nearer than rmsdcut in the list */
      max=0;
      k=0;
      for(j=0; (j<n1) && (mat[i][row[j].j] < rmsdcut); j++)
	if (row[j].j != i) {
	  if (k >= max) {
	    max += 10;
	    srenew(nnb[i],max);
	  }
	  nnb[i][k] = row[j].j;
	  k++;
	}
      if (k == max)
	srenew(nnb[i],max+1);
      nnb[i][k] = -1;
    }
  }
  sfree(row);
  if (debug) {
    fprintf(debug,"Nearest neighborlist. M = %d, P = %d\n",M,P);
    for(i=0; (i<n1); i++) {
      fprintf(debug,"i:%5d nbs:",i);
      for(j=0; nnb[i][j]>=0; j++)
	fprintf(debug,"%5d[%5.3f]",nnb[i][j],mat[i][nnb[i][j]]);
      fprintf(debug,"\n");
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
  qsort(c,n1,sizeof(c[0]),clust_id_comp);

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

static void dump_nnb (FILE *fp, char *title, int n1, t_nnb *nnb)
{
  int i,j;
  
  /* dump neighbor list */
  fprintf(fp,title);
  for(i=0; (i<n1); i++) {
    fprintf(fp,"i:%5d #:%5d nbs:",i,nnb[i].nr);
    for(j=0; j<nnb[i].nr; j++)
      fprintf(fp,"%5d",nnb[i].nb[j]);
    fprintf(fp,"\n");
  }
}
  
static void gromos(int n1, real **mat, real rmsdcut, t_clusters *clust)
{
  char buf[STRLEN];
  t_dist *row;
  t_nnb  *nnb;
  int    i,j,k,i1,j1,max;

  /* Put all neighbors nearer than rmsdcut in the list */
  fprintf(stderr,"Making list of neighbors within cutoff ");
  snew(nnb,n1);
  snew(row,n1);
  for(i=0; (i<n1); i++) {
    /* get one row from matrix */
    for(j=0; (j<n1); j++) {
      row[j].j    = j;
      row[j].dist = mat[i][j];
    }
    /* sort it on rms distance, smallest first */
    qsort(row,n1,sizeof(row[0]),rms_dist_comp);
    max=0;
    k=0;
    /* copy neighbors within cut-off to list */
    for(j=0; j<n1 && row[j].dist < rmsdcut; j++) {
      if (k >= max) {
	max += 10;
	srenew(nnb[i].nb,max);
      }
      nnb[i].nb[k] = row[j].j;
      k++;
    }
    /* store nr of neighbors, we'll need that */
    nnb[i].nr = k;
    if (i%(1+n1/100)==0) fprintf(stderr,"%3d%%\b\b\b\b",(i*100+1)/n1);
  }
  fprintf(stderr,"%3d%%\n",100);
  sfree(row);
  
  /* sort neighbor list on number of neighbors, largest first */
  qsort(nnb,n1,sizeof(nnb[0]),nrnb_comp);

  if (debug) dump_nnb(debug, "Nearest neighborlist after sort.\n", n1, nnb);
  
  /* turn first structure with all its neighbors into cluster
     remove them from pool of structures and repeat for all remaining */
  fprintf(stderr,"Finding clusters\n");
  /* cluster id's start at 1 */
  k=1;
  while(nnb[0].nr) {
    /* set cluster id (k) for first item in neighborlist */
    for (j=0; j<nnb[0].nr; j++)
      clust->cl[nnb[0].nb[j]] = k;
    /* mark as done */
    nnb[0].nr=0;
    sfree(nnb[0].nb);
    
    /* adjust number of neighbors for others, taking removals into account: */
    for(i=1; i<n1 && nnb[i].nr; i++) {
      j1=0;
      for(j=0; j<nnb[i].nr; j++)
	/* if this neighbor wasn't removed */
	if ( clust->cl[nnb[i].nb[j]] == 0 ) {
	  /* shift the rest (j1<=j) */
	  nnb[i].nb[j1]=nnb[i].nb[j];
	  /* next */
	  j1++;
	}
      /* now j1 is the new number of neighbors */
      nnb[i].nr=j1;
    }
    /* sort again, because we have new # neighbors: */
    qsort(nnb,n1,sizeof(nnb[0]),nrnb_comp);
    
    /* new cluster id */
    k++;
  }
  if (debug) {
    fprintf(debug,"Clusters (%d):\n", k);
    for(i=0; i<n1; i++)
      fprintf(debug," %3d", clust->cl[i]);
    fprintf(debug,"\n");
  }

  clust->ncl=k-1;
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
  fprintf(stderr,"Allocated %lu bytes for frames\n",
	  (unsigned long) (max_nf*isize*sizeof(**xx)));
  fprintf(stderr,"Read %d frames from trajectory %s\n",i0,fn);
  *nframe = i0;
  sfree(x);
  
  return xx;
}

static void plot_clusters(int nf, real **mat, real val, t_clusters *clust,
			  int nlevels, int keepfree)
{
  int i,j,v;
  
  for(i=0; i<nf; i++) {
    for(j=0; j<i; j++)
      if (clust->cl[i] == clust->cl[j]) {
	/* color different clusters with different colors, as long as
	   we don't run out of colors */
	v = nlevels - clust->cl[i] + 1; /* because cl >= 1 */
	/* don't use some of the colors, otherwise it gets too light */
	if (keepfree>=nlevels)
	  v = nlevels;
	else if ( keepfree<0 )
	  v = max(v, nlevels/(-keepfree) + 1);
	else if ( keepfree>0 )
	  v = max(v, keepfree);
	else 
	  v = 1;
	/* sadly, we have to convert to real now */
	mat[i][j] = val*v/(real)(nlevels-1);
      } else
	mat[i][j] = 0;
  }
}

static void analyze_clusters(int nf,t_clusters *clust,real **rmsd,
			     t_atoms *atoms,int natom,
			     int isize,atom_id *index,atom_id *alli,
			     rvec **xx,real *time,real *mass,
			     rvec *xtps,char *trxfn,char *sizefn,
			     bool bAverage,bool bFirstIsMiddle, bool bWriteAll,
			     FILE *logf)
{
  FILE *fp;
  char buf[STRLEN],buf1[20],buf2[20],*ext;
  int  trxout=0;
  int  i,i1,i2,cl,nstr,*structure,first=0,midstr;
  real r,clrmsd,midrmsd;
  rvec *xtpsi=NULL,*xav=NULL,*xnatom=NULL;
  matrix zerobox;

  clear_mat(zerobox);

  sprintf(buf,"\nFound %d clusters\n\n",clust->ncl);
  fprintf(stderr,buf);
  fprintf(logf,buf);
  if (trxfn) {
    /* do we write all structures? */
    if (bWriteAll) {
      if (strchr(trxfn,'%'))
	fatal_error(0,"will not number filename %s containing '%c'",trxfn,'%');
      /* number of digits needed in numbering */
      i1 = (int)(log(clust->ncl)/log(10)) + 1;
      /* split fn and ext */
      ext = strrchr(trxfn, '.');
      if (!ext)
	fatal_error(0,"cannot separate extension in filename %s",trxfn);
      ext[0] = '\0';
      ext++;
      /* insert e.g. '%03d' between fn and ext */
      sprintf(buf,"%s%%0%dd.%s",trxfn,i1,ext);
      srenew(trxfn,strlen(buf)+1);
      strcpy(trxfn, buf);
    }
    sprintf(buf,"Writing %s structure%s for each cluster to %s\n",
	    bWriteAll ? "all" : bAverage ? "average" : "middle",
	    bWriteAll ? "s" : "", trxfn);
    fprintf(stderr,"%s",buf);
    fprintf(logf,"%s",buf);
  
    /* Prepare a reference structure for the orientation of the clusters  */
    snew(xtpsi,isize);
    for(i=0; i<isize; i++)
      copy_rvec(xtps[index[i]],xtpsi[i]);
    reset_x(isize,alli,isize,alli,xtpsi,mass);
    if (!bWriteAll)
      trxout = open_trx(trxfn,"w");
    /* Calculate the average structure in each cluster,               *
     * all structures are fitted to the first struture of the cluster */
    snew(xav,isize);
    snew(xnatom,natom);
  }
  snew(structure,nf);
  if (sizefn) {
    fp=xvgropen(sizefn,"Cluster Sizes","Cluster #","# Structures");
    fprintf(fp,"@g%d type %s\n",0,"bar");
  }
  fprintf(logf,"\n%3s | %3s %4s | %6s %4s | cluster members\n",
	  "cl.","#st","rmsd","middle","rmsd");
  for(cl=1; cl<=clust->ncl; cl++) {
    if (xav)
      for(i=0; i<isize;i++)
	clear_rvec(xav[i]);
    nstr=0;
    for(i1=0; i1<nf; i1++)
      if (clust->cl[i1] == cl) {
	structure[nstr] = i1;
	nstr++;
	if (trxfn && bAverage) {
	  reset_x(isize,alli,isize,alli,xx[i1],mass);
	  if (nstr == 1)
	    first = i1;
	  else
	    do_fit(isize,mass,xx[first],xx[i1]);
	  if (xav)
	    for(i=0; i<isize; i++)
	      rvec_inc(xav[i],xx[i1][i]);
	}
      }
    if (sizefn)
      fprintf(fp,"%8d %8d\n",cl,nstr);
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
      if ( (  bFirstIsMiddle && (i1==0) ) || 
	   ( !bFirstIsMiddle && (r < midrmsd) ) ) {
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
    fprintf(logf,"%3d | %3d%s | %6g%s |",
	    cl,nstr,buf1,time[midstr],buf2);
    for(i=0; i<nstr; i++) {
      if ((i % 7 == 0) && i)
	sprintf(buf,"\n%3s | %3s %4s | %6s %4s |","","","","","");
      else
	buf[0] = '\0';
      i1 = structure[i];
      fprintf(logf,"%s %6g",buf,time[i1]);
    }
    fprintf(logf,"\n");

    if (trxfn) {
      if (bWriteAll) {
	/* Dump all structures for this cluster */
	/* generate numbered filename (there is a %d in trxfn!) */
	sprintf(buf,trxfn,cl);
	trxout = open_trx(buf,"w");
	for(i=0; i<nstr; i++)
	  write_trx(trxout,isize,index,atoms,0,time[structure[i]],zerobox,
		    xx[structure[i]],NULL);
	close_trx(trxout);
      } else {
	/* Dump the average structure for this cluster */
	if (bAverage) {
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
    }
  }
  if (trxfn) {
    if (!bWriteAll)
      close_trx(trxout);
    sfree(xtpsi);
    sfree(xav);
    sfree(xnatom);
  }
  sfree(structure);
}

static void convert_mat(t_matrix *mat,t_mat *rms)
{
  int i,j;

  rms->n1 = mat->nx;
  matrix2real(mat,rms->mat);
  
  for(i=0; i<mat->nx; i++)
    for(j=i; j<mat->nx; j++) {
      rms->sumrms += rms->mat[i][j];
      if (rms->mat[i][j] > rms->maxrms)
	rms->maxrms = rms->mat[i][j];
    }
  rms->nn = mat->nx;
}  

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_cluster can cluster structures with several different methods.",
    "Distances between structures can be determined from a trajectory",
    "or read from an XPM matrix file with the [TT]-dm[tt] option.",
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
    
    "gromos: use algorithm as described in Daura [IT]et al.[it]",
    "([IT]Angew. Chem. Int. Ed.[it] [BB]1999[b], [IT]38[it], No. 1/2).",
    "Count number of neighbors using cut-off, take structure with",
    "largest number of neighbors with all its neighbors as cluster",
    "and eleminate it from the pool of clusters. Repeat for remaining",
    "structures in pool.[PAR]",
    
    "When the clustering algorithm assigns each structure to exactly one",
    "cluster (full linkage, Jarvis Patrick and gromos) and a trajectory",
    "file is supplied, the structure with",
    "the smallest average distance to the others or the average structure",
    "or all structures for each cluster will be written to a trajectory",
    "file. When writing all structures, separate numbered files are made",
    "for each cluster."
    
  };
  
  FILE         *fp,*log;
  int          natom,i1,i2,i1s,i2s,i,teller,nf,nrms;
  real         t,t1,t2;

  matrix       box;
  rvec         *xtps,*x1,*x2,**xx=NULL;
  char         *fn;
  t_clusters   clust;
  t_mat        *rms;
  real         **rmsmat;
  real         *resnr,*eigval;
  t_topology   top;
  bool         bSameF;
  t_rgb        rlo,rhi;
  t_matrix     *readmat;
  
  int      status1,status2,isize,nbytes;
  atom_id  *index,*alli=NULL;
  char     *grpname;
  real     rmsd,**d1,**d2,*time,*mass=NULL;
  char     buf[STRLEN];
  bool     bAnalyze,bJP_RMSD=FALSE,bReadMat,bReadTraj,bWriteDist;

  int method;  
  static char *methodname[] = { 
    NULL, "linkage", "jarvis-patrick","monte-carlo", 
    "diagonalization", "gromos", NULL
  };
  enum { m_null, m_linkage, m_jarvis_patrick, 
	 m_monte_carlo, m_diagonalize, m_gromos, m_nr };
  static int  nlevels=40,keepfree=-4,skip=1;
  static real scalemax=-1.0,rmsdcut=0.1;
  static bool bRMSdist=FALSE,bBinary=FALSE,bAverage=FALSE,bWriteAll=FALSE;
  static int  niter=10000,seed=1993;
  static real kT=1e-3;
  static int  M=10,P=3;
  t_pargs pa[] = {
    { "-dista", FALSE, etBOOL, {&bRMSdist},
      "Use RMSD of distances instead of RMS deviation" },
    { "-nlevels",FALSE,etINT,  {&nlevels},
      "Discretize RMSD matrix in # levels" },
    { "-keepfree",FALSE,etINT, {&keepfree},
      "if >0 # levels not to use when coloring clusters; "
      "if <0 nlevels/-keepfree+1 levels will not be used"},
    { "-cutoff",FALSE, etREAL, {&rmsdcut},
      "RMSD cut-off (nm) for two structures to be similar" },
    { "-max",   FALSE, etREAL, {&scalemax},
      "Maximum level in RMSD matrix" },
    { "-skip",  FALSE, etINT,  {&skip},
      "Only analyze every nr-th frame" },
    { "-av",    FALSE, etBOOL, {&bAverage},
      "Write average iso middle structure for each cluster" },
    { "-all",   FALSE, etBOOL, {&bWriteAll},
      "Write all structures for each cluster to numbered files" },
    { "-method",FALSE, etENUM, {methodname},
      "Method for cluster determination" },
    { "-binary",FALSE, etBOOL, {&bBinary},
      "Treat the RMSD matrix as consisting of 0 and 1, where the cut-off "
      "is given by -cutoff" },
    { "-M",     FALSE, etINT,  {&M},
      "Number of nearest neighbors considered for Jarvis-Patrick algorithm, "
      "0 is use cutoff" },
    { "-P",     FALSE, etINT,  {&P},
      "Number of identical nearest neighbors required to form a cluster" },
    { "-seed",  FALSE, etINT,  {&seed},
      "Random number seed for Monte Carlo clustering algorithm" },
    { "-niter", FALSE, etINT,  {&niter},
      "Number of iterations for MC" },
    { "-kT",    FALSE, etREAL, {&kT},
      "Boltzmann weighting factor for Monte Carlo optimization "
      "(zero turns off uphill steps)" }
  };
  t_filenm fnm[] = {
    { efTRX, "-f",     NULL,        ffOPTRD },
    { efTPS, "-s",     NULL,        ffOPTRD },
    { efNDX, NULL,     NULL,        ffOPTRD },
    { efXPM, "-dm",   "rmsd",       ffOPTRD },     
    { efXPM, "-o",    "rmsd-clust", ffWRITE },
    { efLOG, "-g",    "cluster",    ffWRITE },
    { efXVG, "-dist", "rmsd-dist",  ffOPTWR },
    { efXVG, "-ev",   "rmsd-eig",   ffOPTWR },
    { efXVG, "-sz",   "clust-size.xvg",ffOPTWR},
    { efTRX, "-cl",   "clusters.pdb", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* parse options */
  bReadMat   = opt2bSet("-dm",NFILE,fnm);
  bReadTraj  = opt2bSet("-f",NFILE,fnm) || !bReadMat;
  bWriteDist = opt2bSet("-dist",NFILE,fnm) || !bReadMat;
  if (opt2bSet("-cl",NFILE,fnm) && !bReadTraj)
    fprintf(stderr,"\n"
 	    "Warning: cannot write cluster structures without "
	    "reading trajectory\n"
	    "         ignoring option -cl %s\n", opt2fn("-cl",NFILE,fnm));

  method=1;
  while ( method < m_nr && strcasecmp(methodname[0], methodname[method])!=0 )
    method++;
  assert(method != m_nr);
  
  bAnalyze = (method == m_linkage || method == m_jarvis_patrick ||
	      method == m_gromos );

  /* Open log file */
  log = ftp2FILE(efLOG,NFILE,fnm,"w");

  fprintf(stderr,"Using %s method for clustering\n",methodname[0]);
  fprintf(log,"Using %s method for clustering\n",methodname[0]);

  /* check input */
  if (method == m_jarvis_patrick) {
    bJP_RMSD = (M == 0) || opt2parg_bSet("-cutoff",asize(pa),pa);
    if ((M<0) || (M == 1))
      fatal_error(0,"M (%d) must be 0 or larger than 1",M);
    if (M < 2)
      sprintf(buf,"Will use P=%d and RMSD cutoff (%g)",P,rmsdcut);
    else {
      if (P >= M)
	fatal_error(0,"Number of neighbors required (P) must be less than M");
      if (bJP_RMSD)
	sprintf(buf,"Will use P=%d, M=%d and RMSD cutoff (%g)",P,M,rmsdcut);
      else
	sprintf(buf,"Will use P=%d, M=%d",P,M);
    }
    fprintf(stderr,"%s for determining the neighbors\n\n",buf);
    fprintf(log,"%s for determining the neighbors\n\n",buf);
  }

  if (skip < 1)
    fatal_error(0,"skip (%d) should be >= 1",skip);

  /* get input */
  if (bReadTraj) {
    /* don't read mass-database as masses (and top) are not used */
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),buf,&top,&xtps,NULL,box,bAnalyze);
    
    fprintf(stderr,"\nSelect group for %soutput:\n",
	    bReadMat?"":"RMSD calculation and ");
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      1,&isize,&index,&grpname);
  }
  /* Initiate arrays */  
  snew(d1,isize);
  snew(d2,isize);
  for(i=0; (i<isize); i++) {
    snew(d1[i],isize);
    snew(d2[i],isize);
  }

  if (bReadTraj) {
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
  }
  if (bReadMat) {
    fprintf(stderr,"Reading rms distance matrix ");
    read_xpm_matrix(opt2fn("-dm",NFILE,fnm),&readmat);
    fprintf(stderr,"\n");
    if (readmat[0].nx != readmat[0].ny)
      fatal_error(0,"Matrix (%dx%d) is not square",
		  readmat[0].nx,readmat[0].ny);
    if (bReadTraj && bAnalyze && (readmat[0].nx != nf))
      fatal_error(0,"Matrix size (%dx%d) does not match the number of frames (%d)",
	      readmat[0].nx,readmat[0].ny,nf);

    nf = readmat[0].nx;
    sfree(time);
    time = readmat[0].axis_x;

    rms = init_mat(readmat[0].nx,method == m_diagonalize);
    convert_mat(&(readmat[0]),rms);
    
    nlevels = readmat[0].nmap;
  } else { /* !bReadMat */
    rms = init_mat(nf,method == m_diagonalize);
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
    } else { /* bRMSdist */
      fprintf(stderr,"Computing %dx%d RMS distance deviation matrix\n",nf,nf);
      for(i1=0; (i1<nf); i1++) {
	calc_dist(isize,xx[i1],d1);
	for(i2=i1+1; (i2<nf); i2++) {
	  calc_dist(isize,xx[i2],d2);
	  set_mat_entry(rms,i1,i2,rms_dist(isize,d1,d2));
      }
	nrms -= (nf-i1-1);
	fprintf(stderr,"\r# RMSD calculations left: %d   ",nrms);
      }
    }
    fprintf(stderr,"\n\n");
  }
  sprintf(buf,"The maximum RMSD is %g nm\n"
	  "Average RMSD is %g\n"
	  "Number of structures for matrix %d\n"
	  "Energy of the matrix is %g nm\n",
	  rms->maxrms,2*rms->sumrms/(nf*(nf-1)),nf,mat_energy(rms));
  fprintf(stderr,buf);
  fprintf(log,buf);

  if (bWriteDist)
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
  switch (method) {
  case m_linkage: 
    /* Now sort the matrix and write it out again */
    gather(rms,rmsdcut,&clust);
    break;
  case m_diagonalize:
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
    break;
  case m_monte_carlo:
    mc_optimize(log,rms,niter,&seed,kT);
    break;
  case m_jarvis_patrick:
    jarvis_patrick(rms->nn,rms->mat,M,P,bJP_RMSD ? rmsdcut : -1,&clust);
    break;
  case m_gromos:
    gromos(rms->nn,rms->mat,rmsdcut,&clust);
    break;
  default:
    fatal_error(0,"unknown method \"%s\"",methodname[0]);
  }

  if (method == m_monte_carlo || method == m_diagonalize)
    fprintf(stderr,"Energy of the matrix after clustering is %g nm\n",
	    mat_energy(rms));
  swap_mat(rms);
  reset_index(rms);
  
  /* Write the rmsd-matrix to the upper-left part of the matrix */
  for(i1=0; (i1 < nf); i1++)
    for(i2=i1; (i2 < nf); i2++)
      rms->mat[i1][i2] = rmsmat[i1][i2];
  
  if (bAnalyze) {
    plot_clusters(nf,rms->mat,rms->maxrms,&clust,nlevels,keepfree);
    analyze_clusters(nf,&clust,rmsmat,&(top.atoms),natom,isize,index,alli,
		     xx,time,mass,xtps,
		     bReadTraj ? opt2fn_null("-cl",NFILE,fnm) : NULL, 
		     opt2fn_null("-sz",NFILE,fnm), 
		     bAverage, method==m_gromos, bWriteAll, log);
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
  fprintf(stderr,"Writing rms distance/clustering matrix ");
  if (bReadMat) {
    write_xpm(fp,readmat[0].title,readmat[0].legend,readmat[0].label_x,
	      readmat[0].label_y,nf,nf,readmat[0].axis_x,readmat[0].axis_y,
	      rms->mat,0.0,rms->maxrms,rlo,rhi,&nlevels);
  } else
    write_xpm(fp,bRMSdist ? "RMS Distance Deviation" : "RMS Deviation",
	      "RMSD (nm)","Time (ps)","Time (ps)",
	      nf,nf,time,time,rms->mat,0.0,rms->maxrms,rlo,rhi,&nlevels);
  fprintf(stderr,"\n");
  ffclose(fp);
  
  /* now show what we've done */
  do_view(opt2fn("-o",NFILE,fnm),NULL);
  do_view(opt2fn_null("-sz",NFILE,fnm),NULL);
  if (method == m_diagonalize)
    do_view(opt2fn_null("-ev",NFILE,fnm),NULL);
  if (bWriteDist)
    do_view(opt2fn("-dist",NFILE,fnm),NULL);

  /* Thank the user for her patience */  
  thanx(stderr);
  
  return 0;
}
