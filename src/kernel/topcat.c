/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
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
  int  index,max_index=0;
  real fac,type;
  bool bDisres;
  
  nrfp        = NRFP(ftype);
  nral        = NRAL(ftype);
  src_natoms  = dstart;
  dest_natoms = nrstart;
  
  n0     = nrstart;

  /* Add this many entries to array ... */
  pr_alloc(src->nr*copies,dest);
  
  /* First new entry. */
  l=dest->nr;
 
  bDisres = (ftype == F_DISRES);
  if (bDisres) {
    for(i=0; i<src->nr; i++)
      max_index = max(src->param[i].c[0],max_index);
  }

  if (bDisres && bEnsemble) {
    /* If we have to do ensemble averaging we interchange the loops! */
    if (src->nr > 0) {
      fac=pow((real)copies,-1.0/6.0);
      fprintf(stderr,"multiplying every type 1 distance bound by %.3f for %d copies of %s\n",fac,copies,name);
      for (i=0; (i<src->nr); i++) {
	for (j=0; (j<copies); j++) {
	  index = src->param[i].c[0];
	  type  = src->param[i].c[1];
	  if (type == 2) {
	    /* Don't ensemble average type' 2 distance restraints */
	    index += j*(max_index + 1);
	  }
	  dest->param[l].c[0] = index;
	  dest->param[l].c[1] = type;
	  for(m=0; (m<2); m++)
	    dest->param[l].c[m] = src->param[i].c[m];
	  if (type==2)
	    for(m=2; (m<nrfp-1); m++)
	      dest->param[l].c[m] = src->param[i].c[m];
	  else
	    /* Scale the bounds for ensemble averaged distance restraints */
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
	if (bDisres)
	  dest->param[l].c[0] += j*(max_index + 1);
	l++;
      }
      n0 += dstart;
    }
  }
  if (l != dest->nr+src->nr*copies)
    gmx_fatal(FARGS,"In %s line %d: l = %d, should be %d\n",
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
    srenew(dest->atomtype,size);
    srenew(dest->atomtypeB,size);
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
    memcpy((char *) &(dest->atomtype[l]),(char *) &(src->atomtype[0]),
           (size_t)(srcnr*sizeof(src->atomtype[0])));
    memcpy((char *) &(dest->atomtypeB[l]),(char *) &(src->atomtypeB[0]),
           (size_t)(srcnr*sizeof(src->atomtypeB[0])));
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
	    nrcopies,destnr,srcnr,bEnsemble,*src->name,i);
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
      range_check(n,0,nsrc);
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
    fprintf(debug,"  emptiest=%d\n",min_i);

  return min_i;
}

int *mk_shuffle_tab(int nmol,t_molinfo mol[],int nnodes,int *ntab,
		    int Nsim,t_simsystem Sims[],bool bVerbose)
{
  t_simsystem *sss;
  int  *tab,**bucket,*sim_index;
  int  i,j,k,nm,nmolnz,idum,natom,eb;
  
  nm     = 0;
  for(i=0; (i<Nsim); i++) 
    nm+=Sims[i].nrcopies;
  
  snew(bucket,nnodes);
  for(i=0; (i<nnodes); i++)
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
      eb = emptiest_bucket(nnodes,bucket,nmol,mol);
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
    for(j=0; (j<nnodes); j++) {
      fprintf(stderr,"CPU%4d",j);
      for(i=0; (i<nmol); i++)
	fprintf(stderr,"  %8d",bucket[j][i]);
      fprintf(stderr,"  %8d\n",count_atoms(bucket[j],nmol,mol));
    }
  }
  snew(tab,nm);
  *ntab = 0;
  for(j=0; (j<nnodes); j++)
    for(k=0; (k<nmol); k++)
      for(i=0; (i<bucket[j][k]); i++) {
	if (*ntab >= nm)
	  gmx_incons("Shuffle table not large enough");
	tab[(*ntab)++] = k;
      }
  return tab;
}

int *mk_shuffle_tab_old(int nmol,t_molinfo mol[],int nnodes,int *ntab,
			int Nsim,t_simsystem Sims[],bool bVerbose)
{
  int  *tab,*tmol;
  int  i,j,k,nm,ttt,ifrac,nmolnz,idum,natom;
  real frac,rnnodes;
  
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
  rnnodes=nnodes;
  for(i=k=0; (i<nnodes); i++) {
    if (bVerbose) 
      fprintf(stderr,"%-6s%4d","CPU",i);
    natom = 0;
    
    for(j=0; (j<Nsim); j++) {
      frac=((i+1)*Sims[j].nrcopies)/rnnodes;
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
