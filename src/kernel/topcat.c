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
 * GROtesk MACabre and Sinister
 */
static char *SRCID_topcat_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "assert.h"
#include "macros.h"
#include "topio.h"
#include "toputil.h"
#include "toppush.h"
#include "topcat.h"

static void bondcat(t_params *dest,t_params *src,int copies,int nrstart,
		    int dstart,bool bEnsemble,char *name,int ftype)
/* append all the bonds from src to dest */
{
  int  i,j,l,m,n0,nrfp,nral;
  int  src_natoms,dest_natoms;
  int  index,max_index;
  real fac,type;
  
  nrfp        = NRFP(ftype);
  nral        = NRAL(ftype);
  src_natoms  = dstart;
  dest_natoms = nrstart;
  
  n0     = nrstart;
  
  /* Add this many entries to array ... */
  pr_alloc(src->nr*copies,dest);
  
  /* First new entry. */
  l=dest->nr;
 
  if (src->nr)
    max_index = src->param[0].c[0];
  for(i=0; i<src->nr; i++)
    max_index = max(src->param[i].c[0],max_index);
      
  /* If we have to do ensemble averaging we interchange the loops! */
  if (bEnsemble) {
    if (src->nr > 0) {
      fac=pow((real)copies,-1.0/6.0);
      fprintf(stderr,"multiplying every type 1 distance bound by %.3f for %d copies of %s\n",fac,copies,name);
      for (i=0; (i<src->nr); i++) {
	for (j=0; (j<copies); j++) {
	  index = src->param[i].c[0];
	  type  = src->param[i].c[1];
	  if ((type == 2) && (j>0)) {
	    max_index++;
	    index = max_index;
	  }
	  dest->param[l].c[0] = index;
	  dest->param[l].c[1] = type;
	  for(m=0; (m<2); m++)
	    dest->param[l].c[m] = src->param[i].c[m];
	  for(m=2; (m<nrfp-1); m++)
	    dest->param[l].c[m] = src->param[i].c[m]*fac;
	  /* do not change the factor for the force-constant */
	  dest->param[l].c[nrfp-1] = src->param[i].c[nrfp-1];
	  for (m=0; (m<nral); m++)
	    dest->param[l].a[m] = src->param[i].a[m]+
	      j*src_natoms+dest_natoms;
	  l++;
	}
      }
    }
  }
  else {
    for (j=0; (j<copies); j++) {
      for (i=0; (i<src->nr); i++) {
	memcpy((char *)&(dest->param[l]),
	       (char *)&(src->param[i]),(size_t)sizeof(src->param[i]));
	for (m=0; (m<nral); m++)
	  dest->param[l].a[m] += n0;
	if (j>0) {
	  max_index++;
	  dest->param[l].c[0] = max_index;
	}
	l++;
      }
      n0 += dstart;
    }
  }
  if (l != dest->nr+src->nr*copies)
    fatal_error(0,"In %s line %d: l = %d, should be %d\n",
		__FILE__,__LINE__,l,dest->nr+src->nr*copies);
  dest->nr = l;
}

static void blockcat(t_block *dest,t_block *src,int copies, 
		     int dnum,int snum)
{
  int i,j,l,size;
  int destnr  = dest->nr;
  int destnra = dest->nra;

  if (src->nr) {
    size=(dest->nr+copies*src->nr+1);
    srenew(dest->index,size);
  }
  if (src->nra) {
    size=(dest->nra+copies*src->nra);
    srenew(dest->a,size);
  }

  for (l=destnr,j=0; (j<copies); j++) {
    for (i=0; (i<src->nr); i++)
      dest->index[l++] = dest->nra+src->index[i];
    dest->nra += src->nra;
  }
  for (l=destnra,j=0; (j<copies); j++) {
    for (i=0; (i<src->nra); i++)
      dest->a[l++] = dnum+src->a[i];
    dnum+=snum;
    dest->nr += src->nr;
  }
  dest->index[dest->nr] = dest->nra;
}

static void atomcat (t_atoms *dest, t_atoms *src, int copies)
{
  int i,j,l,size;
  int srcnr=src->nr;
  int destnr=dest->nr;

  if (srcnr) {
    size=destnr+copies*srcnr;
    srenew(dest->atom,size);
    srenew(dest->atomname,size);
  }
  if (src->nres) {
    size=dest->nres+copies*src->nres;
    srenew(dest->resname,size);
  }

  /* residue information */
  for (l=dest->nres,j=0; (j<copies); j++,l+=src->nres)
    memcpy((char *) &(dest->resname[l]),(char *) &(src->resname[0]),
           (size_t)(src->nres*sizeof(src->resname[0])));

  for (l=destnr,j=0; (j<copies); j++,l+=srcnr) {
    memcpy((char *) &(dest->atomname[l]),(char *) &(src->atomname[0]),
           (size_t)(srcnr*sizeof(src->atomname[0])));
    memcpy((char *) &(dest->atom[l]),(char *) &(src->atom[0]),
           (size_t)(srcnr*sizeof(src->atom[0])));
  }

  /* Increment residue numbers */
  for (l=destnr,j=0; (j<copies); j++)
    for (i=0; (i<srcnr); i++,l++) 
      dest->atom[l].resnr  = dest->nres+j*src->nres+src->atom[i].resnr;

  dest->nres += copies*src->nres;
  dest->nr += copies*src->nr;

  blockcat(&(dest->excl),&(src->excl),copies,destnr,srcnr);
}

static void top1_cat(t_molinfo *dest,t_molinfo *src,
		     int nrcopies,bool bEnsemble)
{
  int srcnr,destnr,i;
  
  if (nrcopies <= 0) 
    return;
  
  /* concat atoms */
  srcnr  = src->atoms.nr;
  destnr = dest->atoms.nr;
  if (debug)
    fprintf(debug,"topcat: srcnr: %d, destnr: %d\n",srcnr,destnr);
  
  atomcat (&(dest->atoms),&(src->atoms),nrcopies);
  
  blockcat(&(dest->cgs),&(src->cgs),nrcopies,destnr,srcnr);
  blockcat(&(dest->mols),&(src->mols),nrcopies,destnr,srcnr);
  
  for (i=0; (i<F_NRE); i++) 
    bondcat(&dest->plist[i],&src->plist[i],
	    nrcopies,destnr,srcnr,(i==F_DISRES) ? bEnsemble : FALSE,
	    *src->name,i);
}
	     
void topcat (t_molinfo *dest,int nsrc,t_molinfo src[],int ntab,
	     int *tab,int Nsim,t_simsystem Sims[],bool bEnsemble)
/* concatenate all copies of src to dest */
{
  int i,n;

  if (ntab > 0) {
    for(i=0; (i<ntab); i++) {
      top1_cat(dest,&(src[tab[i]]),1,FALSE);
    }
  }
  else {
    for(i=0; (i<Nsim); i++) {
      n=Sims[i].whichmol;
      assert((0<=n) && (n<nsrc));
      top1_cat(dest,&(src[n]),Sims[i].nrcopies,bEnsemble);
    }
  }
}

void mi2top(t_topology *dest,t_molinfo *src)
{
  atomcat(&(dest->atoms),&(src->atoms),1);
  blockcat(&(dest->blocks[ebCGS]),&(src->cgs),1,0,src->atoms.nr);
  blockcat(&(dest->blocks[ebMOLS]),&(src->mols),1,0,src->atoms.nr);
  dest->name=src->name;
}

static t_molinfo   *mol_tmp;
static t_simsystem *sim_tmp;

static int sim_comp(const void *a,const void *b)
{
  return (mol_tmp[sim_tmp[*((int *)b)].whichmol].atoms.nr -
	  mol_tmp[sim_tmp[*((int *)a)].whichmol].atoms.nr);
}

static int count_atoms(int *bucket,int nmol,t_molinfo mol[])
{
  int i,n;
  
  n=0;
  for(i=0; (i<nmol); i++)
    n += bucket[i]*mol[i].atoms.nr;
  return n;
}

static int emptiest_bucket(int nbucket,int **bucket,int nmol,t_molinfo mol[])
{
  int i,ca,min_i,min_nat;
  
  min_i   = 0;
  min_nat = count_atoms(bucket[min_i],nmol,mol);
  if (debug)
    fprintf(debug,"BUCKET  %5d",min_nat);
  for(i=1; (i<nbucket); i++) {
    ca = count_atoms(bucket[i],nmol,mol);
    if (debug)
      fprintf(debug,"  %5d",ca);
    if (ca < min_nat) {
      min_nat = ca;
      min_i   = i;
    }
  }
  if (debug)
    fprintf(debug,"  eb=%d\n",min_i);

  return min_i;
}

int *mk_shuffle_tab(int nmol,t_molinfo mol[],int nprocs,int *ntab,
		    int Nsim,t_simsystem Sims[],bool bVerbose)
{
  t_simsystem *sss;
  int  *tab,**bucket,*sim_index;
  int  i,j,k,nm,nmolnz,idum,natom,eb;
  
  nm     = 0;
  for(i=0; (i<Nsim); i++) 
    nm+=Sims[i].nrcopies;
  
  snew(bucket,nprocs);
  for(i=0; (i<nprocs); i++)
    snew(bucket[i],nmol);
    
  /* Sort the simsystems to increasing size of molecules */
  snew(sim_index,Nsim);
  for(i=0; (i<Nsim); i++)
    sim_index[i] = i;
  sim_tmp = Sims;
  mol_tmp = mol;
  qsort(sim_index,Nsim,sizeof(sim_index[0]),sim_comp);
  
  for(i=0; (i<Nsim); i++) {
    sss = &(Sims[sim_index[i]]);
    for(j=0; (j<sss->nrcopies); j++) {
      eb = emptiest_bucket(nprocs,bucket,nmol,mol);
      bucket[eb][sss->whichmol]++;
    }
  }
  
  /* Print a header for the table */
  if (bVerbose) {
    fprintf(stderr,"%7s","Moltype");
    for(i=0; (i<nmol); i++) 
      fprintf(stderr,"  %8s",*(mol[i].name));
    fprintf(stderr,"  %8s\n","#atoms");

    /* Print the table itself */
    for(j=0; (j<nprocs); j++) {
      fprintf(stderr,"CPU%4d",j);
      for(i=0; (i<nmol); i++)
	fprintf(stderr,"  %8d",bucket[j][i]);
      fprintf(stderr,"  %8d\n",count_atoms(bucket[j],nmol,mol));
    }
  }
  snew(tab,nmol*nprocs);
  *ntab = 0;
  for(j=0; (j<nprocs); j++)
    for(k=0; (k<nmol); k++)
      for(i=0; (i<bucket[j][k]); i++)
	tab[(*ntab)++] = k;

  return tab;
}

int *mk_shuffle_tab_old(int nmol,t_molinfo mol[],int nprocs,int *ntab,
			int Nsim,t_simsystem Sims[],bool bVerbose)
{
  int  *tab,*tmol;
  int  i,j,k,nm,ttt,ifrac,nmolnz,idum,natom;
  real frac,rnprocs;
  
  nm     = 0;
  nmolnz = 0;
  for(i=0; (i<Nsim); i++) {
    nm+=Sims[i].nrcopies;
    if (Sims[i].nrcopies)
      nmolnz++;
  }
  fprintf(stderr,"Total number of molecules = %d, Number of Simsystems = %d\n",
	  nm,Nsim);
  
  /* Print a header for the table */
  if (bVerbose) {
    fprintf(stderr,"%10s","Moltype");
    for(i=0; (i<Nsim); i++)
      if (Sims[i].nrcopies > 0)
	fprintf(stderr,"  %8s",*(mol[Sims[i].whichmol].name));
    fprintf(stderr,"  %8s\n","#atoms");
  }
  snew(tab,nm);
  snew(tmol,Nsim);
  rnprocs=nprocs;
  for(i=k=0; (i<nprocs); i++) {
    if (bVerbose) 
      fprintf(stderr,"%-6s%4d","CPU",i);
    natom = 0;
    
    for(j=0; (j<Nsim); j++) {
      frac=((i+1)*Sims[j].nrcopies)/rnprocs;
      if (debug)
	fprintf(debug,"CPU=%3d, MOL=%3d, frac = %g\n",i,j,frac);
      ifrac = frac;
      idum  = tmol[j];
      while(tmol[j] < ifrac) {
	tab[k++] = Sims[j].whichmol;
	tmol[j]++;
	natom   += mol[Sims[j].whichmol].atoms.nr;
      }
      if (bVerbose) {
	fprintf(stderr,"  %8d",tmol[j]-idum);
      }
    }
    if (bVerbose)
      fprintf(stderr,"  %8d\n",natom);
  }
  sfree(tmol);
  
  if (debug)
    for(i=0;(i<nm); i++)
      fprintf(debug,"shuffle_tab[%d] = %d\n",i,tab[i]);
  
  if (bVerbose && 0) {
    for(i=0; (i<nm); ) {
      ttt=tab[i];
      k=0;
      while ((i<nm) && (ttt==tab[i])) {
	i++;
	k++;
      }
      fprintf(stderr,"Mol: %20s  %5d\n",*mol[ttt].name,k);
    }
  }
  *ntab=nm;
  
  return tab;
}
