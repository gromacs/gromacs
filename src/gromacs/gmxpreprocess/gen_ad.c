/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <ctype.h>
#include "sysstuff.h"
#include "macros.h"
#include "smalloc.h"
#include "string2.h"
#include "confio.h"
#include "vec.h"
#include "pbc.h"
#include "toputil.h"
#include "topio.h"
#include "gpp_nextnb.h"
#include "symtab.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "pgutil.h"
#include "resall.h"
#include "gen_ad.h"

typedef gmx_bool (*peq)(t_param *p1, t_param *p2);

static int acomp(const void *a1, const void *a2)
{
  t_param *p1,*p2;
  int     ac;
  
  p1=(t_param *)a1;
  p2=(t_param *)a2;
  if ((ac=(p1->AJ-p2->AJ))!=0)
    return ac;
  else if ((ac=(p1->AI-p2->AI))!=0)
    return ac;
  else 
    return (p1->AK-p2->AK);
}

static int pcomp(const void *a1, const void *a2)
{
  t_param *p1,*p2;
  int     pc;
  
  p1=(t_param *)a1;
  p2=(t_param *)a2;
  if ((pc=(p1->AI-p2->AI))!=0)
    return pc;
  else 
    return (p1->AJ-p2->AJ);
}

static int dcomp(const void *d1, const void *d2)
{
  t_param *p1,*p2;
  int     dc;
  
  p1=(t_param *)d1;
  p2=(t_param *)d2;
  /* First sort by J & K (the two central) atoms */
  if ((dc=(p1->AJ-p2->AJ))!=0)
    return dc;
  else if ((dc=(p1->AK-p2->AK))!=0)
    return dc;
  /* Then make sure to put rtp dihedrals before generated ones */
  else if (p1->c[MAXFORCEPARAM-1]==0 && p2->c[MAXFORCEPARAM-1]==NOTSET)
    return -1;
  else if (p1->c[MAXFORCEPARAM-1]==NOTSET && p2->c[MAXFORCEPARAM-1]==0)
    return 1;
  /* Finally, sort by I and J (two outer) atoms */
  else if ((dc=(p1->AI-p2->AI))!=0)
    return dc;
  else
    return (p1->AL-p2->AL);
}


static gmx_bool aeq(t_param *p1, t_param *p2)
{
  if (p1->AJ!=p2->AJ)
    return FALSE;
  else if (((p1->AI==p2->AI) && (p1->AK==p2->AK)) ||
           ((p1->AI==p2->AK) && (p1->AK==p2->AI)))
    return TRUE;
  else
    return FALSE;
}



static gmx_bool deq(t_param *p1, t_param *p2)
{
  if (((p1->AJ==p2->AJ) && (p1->AK==p2->AK)) ||
      ((p1->AJ==p2->AK) && (p1->AK==p2->AJ)))
    return TRUE;
  else
    return FALSE;
}


static gmx_bool remove_dih(t_param *p, int i, int np)
     /* check if dihedral p[i] should be removed */
{
  gmx_bool bRem;
  int j;

  if (p[i].c[MAXFORCEPARAM-1]==NOTSET) {
    if (i>0)
      bRem = deq(&p[i],&p[i-1]);
    else
      bRem = FALSE;
    /* also remove p[i] if there is a dihedral on the same bond
       which has parameters set */
    j=i+1;
    while (!bRem && (j<np) && deq(&p[i],&p[j])) {
      bRem = (p[j].c[MAXFORCEPARAM-1] != NOTSET);
      j++;
    }
  } else
    bRem = FALSE;

  return bRem;
}

static gmx_bool preq(t_param *p1, t_param *p2)
{
  if ((p1->AI==p2->AI) && (p1->AJ==p2->AJ))
    return TRUE;
  else 
    return FALSE;
}

static void rm2par(t_param p[], int *np, peq eq)
{
  int *index,nind;
  int i,j;

  if ((*np)==0)
    return;

  snew(index,*np);
  nind=0;
    index[nind++]=0;
  for(i=1; (i<(*np)); i++) 
    if (!eq(&p[i],&p[i-1]))
      index[nind++]=i;
  /* Index now holds pointers to all the non-equal params,
   * this only works when p is sorted of course
   */
  for(i=0; (i<nind); i++) {
    for(j=0; (j<MAXATOMLIST); j++)
      p[i].a[j]=p[index[i]].a[j];
    for(j=0; (j<MAXFORCEPARAM); j++)
      p[i].c[j]=p[index[i]].c[j];
    if (p[index[i]].a[0] == p[index[i]].a[1]) {
      if (debug)  
	fprintf(debug,
		"Something VERY strange is going on in rm2par (gen_ad.c)\n"
		"a[0] %u a[1] %u a[2] %u a[3] %u\n",
		p[i].a[0],p[i].a[1],p[i].a[2],p[i].a[3]);
      strcpy(p[i].s,"");
    } else if (index[i] > i) {
      /* Copy the string only if it comes from somewhere else 
       * otherwise we will end up copying a random (newly freed) pointer.
       * Since the index is sorted we only have to test for index[i] > i.
       */ 
      strcpy(p[i].s,p[index[i]].s);
    }
  }
  (*np)=nind;

  sfree(index);
}

static void cppar(t_param p[], int np, t_params plist[], int ftype)
{
  int      i,j,nral,nrfp;
  t_params *ps;

  ps   = &plist[ftype];
  nral = NRAL(ftype);
  nrfp = NRFP(ftype);
  
  /* Keep old stuff */
  pr_alloc(np,ps);
  for(i=0; (i<np); i++) {
    for(j=0; (j<nral); j++)
      ps->param[ps->nr].a[j] = p[i].a[j];
    for(j=0; (j<nrfp); j++)
      ps->param[ps->nr].c[j] = p[i].c[j];
    for(j=0; (j<MAXSLEN); j++)
      ps->param[ps->nr].s[j] = p[i].s[j];
    ps->nr++;
  }
}

static void cpparam(t_param *dest, t_param *src)
{
  int j;

  for(j=0; (j<MAXATOMLIST); j++)
    dest->a[j] = src->a[j];
  for(j=0; (j<MAXFORCEPARAM); j++)
    dest->c[j] = src->c[j];
  for(j=0; (j<MAXSLEN); j++)
    dest->s[j] = src->s[j];
}

static void set_p(t_param *p, atom_id ai[4], real *c, char *s)
{
  int j;

  for(j=0; (j<4); j++)
    p->a[j]=ai[j];
  for(j=0; (j<MAXFORCEPARAM); j++)
    if (c)
      p->c[j]=c[j];
    else
      p->c[j]=NOTSET;

  set_p_string(p,s);
}

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

static int atom_id_comp(const void *a,const void *b)
{
  return (*(atom_id *)a) - (*(atom_id *)b);
}

static int eq_imp(atom_id a1[],atom_id a2[])
{
  int b1[4],b2[4];
  int j;

  for(j=0; (j<4); j++) {
    b1[j]=a1[j];
    b2[j]=a2[j];
  }
  qsort(b1,4,(size_t)sizeof(b1[0]),int_comp);
  qsort(b2,4,(size_t)sizeof(b2[0]),int_comp);

  for(j=0; (j<4); j++)
    if (b1[j] != b2[j])
      return FALSE;

  return TRUE;
}

static gmx_bool ideq(t_param *p1, t_param *p2)
{
  return eq_imp(p1->a,p2->a);
}

static int idcomp(const void *a,const void *b)
{
  t_param *pa,*pb;
  int     d;
  
  pa=(t_param *)a;
  pb=(t_param *)b;
  if ((d=(pa->a[0]-pb->a[0])) != 0)
    return d;
  else if ((d=(pa->a[3]-pb->a[3])) != 0)
    return d;
  else if ((d=(pa->a[1]-pb->a[1])) != 0)
    return d;
  else
    return (int) (pa->a[2]-pb->a[2]);
}

static void sort_id(int nr,t_param ps[])
{
  int i,tmp;
  
  /* First swap order of atoms around if necessary */
  for(i=0; (i<nr); i++) {
    if (ps[i].a[3] < ps[i].a[0]) {
      tmp = ps[i].a[3]; ps[i].a[3] = ps[i].a[0]; ps[i].a[0] = tmp;
      tmp = ps[i].a[2]; ps[i].a[2] = ps[i].a[1]; ps[i].a[1] = tmp;
    }
  }
  /* Now sort it */
  if (nr > 1)
    qsort(ps,nr,(size_t)sizeof(ps[0]),idcomp);
}

static void dump_param(FILE *fp,char *title,int n,t_param ps[])
{
 int i,j;
  
  fprintf(fp,"%s: %d entries\n",title,n);
  for(i=0; (i<n); i++) {
    fprintf(fp,"%3d:  A=[ ",i);
    for(j=0; (j<MAXATOMLIST); j++)
      fprintf(fp," %5d",ps[i].a[j]);
    fprintf(fp,"]  C=[");
    for(j=0; (j<MAXFORCEPARAM); j++)
      fprintf(fp," %10.5e",ps[i].c[j]);
    fprintf(fp,"]\n");  
  }
}

static int n_hydro(atom_id a[],char ***atomname)
{
  int i,nh=0;
  char c0,c1,*aname;

  for(i=0; (i<4); i+=3) {
    aname=*atomname[a[i]];
    c0=toupper(aname[0]);
    if (c0 == 'H')
      nh++;
    else if (((int)strlen(aname) > 1) && (c0 >= '0') && (c0 <= '9')) {
      c1=toupper(aname[1]);
      if (c1 == 'H')
	nh++;
    }
  }
  return nh;
}

static void clean_dih(t_param *dih, int *ndih,t_param idih[],int nidih,
		      t_atoms *atoms,gmx_bool bAlldih, gmx_bool bRemoveDih)
{
  int  i,j,k,l;
  int  *index,nind;
  gmx_bool bIsSet,bKeep;
  int  bestl,nh,minh;
  
  snew(index,*ndih+1);
  if (bAlldih) {
    fprintf(stderr,"Keeping all generated dihedrals\n");
    nind = *ndih;
    for(i=0; i<nind; i++) 
      index[i] = i;
    index[nind] = *ndih;
  } else {
    /* Make an index of all dihedrals over each bond */
    nind = 0;
    for(i=0; i<*ndih; i++) 
      if (!remove_dih(dih,i,*ndih)) 
	index[nind++]=i;
    index[nind] = *ndih;
  }

  /* if we don't want all dihedrals, we need to select the ones with the 
   *  fewest hydrogens
   */
  
  k=0;
  for(i=0; i<nind; i++) {
    bIsSet = (dih[index[i]].c[MAXFORCEPARAM-1] != NOTSET);
    bKeep = TRUE;
    if (!bIsSet && bRemoveDih)
      /* remove the dihedral if there is an improper on the same bond */
      for(j=0; (j<nidih) && bKeep; j++)
	bKeep = !deq(&dih[index[i]],&idih[j]);

    if (bKeep) {
      /* Now select the "fittest" dihedral:
       * the one with as few hydrogens as possible 
       */
      
      /* Best choice to get dihedral from */
      bestl=index[i];
      if (!bAlldih && !bIsSet) {
	/* Minimum number of hydrogens for i and l atoms */
	minh=2;
	for(l=index[i]; (l<index[i+1]) && deq(&dih[index[i]],&dih[l]); l++) {
	  if ((nh=n_hydro(dih[l].a,atoms->atomname)) < minh) {
	    minh=nh;
	    bestl=l;
	  }
	  if (minh == 0)
	    break;
	}
      }
      if (k != bestl)
	cpparam(&(dih[k]),&dih[bestl]);
      k++;
    }
  }

  for (i=k; i<*ndih; i++)
    strcpy(dih[i].s,"");
  *ndih = k;

  sfree(index);
}

static int get_impropers(t_atoms *atoms,t_hackblock hb[],t_param **idih,
			 gmx_bool bAllowMissing)
{
  char      *a0;
  t_rbondeds *idihs;
  t_rbonded  *hbidih;
  int       nidih,i,j,k,r,start,ninc,nalloc;
  atom_id   ai[MAXATOMLIST];
  gmx_bool      bStop;
  
  ninc = 500;
  nalloc = ninc;
  snew(*idih,nalloc);

  /* Add all the impropers from the residue database to the list. */
  nidih = 0;
  start = 0;
  if (hb != NULL) {
    for(i=0; (i<atoms->nres); i++) {
      idihs=&hb[i].rb[ebtsIDIHS];
      for(j=0; (j<idihs->nb); j++) {
	bStop=FALSE;
	for(k=0; (k<4) && !bStop; k++) {
	  ai[k] = search_atom(idihs->b[j].a[k],start,
			      atoms->nr,atoms->atom,atoms->atomname,
			      "improper",bAllowMissing);
	  if (ai[k] == NO_ATID)
	    bStop = TRUE;
	}
	if (!bStop) {
	  if (nidih == nalloc) {
	    nalloc += ninc;
	    srenew(*idih,nalloc);
	  }
	  /* Not broken out */
	  set_p(&((*idih)[nidih]),ai,NULL,idihs->b[j].s);
	  nidih++;
	}
      }
      while ((start<atoms->nr) && (atoms->atom[start].resind == i))
	start++;
    }
  }
  
  return nidih;
}

static int nb_dist(t_nextnb *nnb,int ai,int aj)
{
  int nre,nrx,NRE;
  int *nrexcl;
  int *a;
  
  if (ai == aj)
    return 0;
  
  NRE=-1;
  nrexcl=nnb->nrexcl[ai];
  for(nre=1; (nre < nnb->nrex); nre++) {
    a=nnb->a[ai][nre];
    for(nrx=0; (nrx < nrexcl[nre]); nrx++) {
      if ((aj == a[nrx]) && (NRE == -1))
	NRE=nre;
    }
  }
  return NRE;
}

gmx_bool is_hydro(t_atoms *atoms,int ai)
{
  return ((*(atoms->atomname[ai]))[0] == 'H');
}

static void get_atomnames_min(int n,char **anm,
			      int resind,t_atoms *atoms,atom_id *a)
{
  int m;

  /* Assume ascending residue numbering */
  for(m=0; m<n; m++) {
    if (atoms->atom[a[m]].resind < resind)
      strcpy(anm[m],"-");
    else if (atoms->atom[a[m]].resind > resind)
      strcpy(anm[m],"+");
    else
      strcpy(anm[m],"");
    strcat(anm[m],*(atoms->atomname[a[m]]));
  }
}

static void gen_excls(t_atoms *atoms, t_excls *excls, t_hackblock hb[],
		      gmx_bool bAllowMissing)
{
  int        r;
  atom_id    a,astart,i1,i2,itmp;
  t_rbondeds *hbexcl;
  int        e;
  char       *anm;

  astart = 0;
  for(a=0; a<atoms->nr; a++) {
    r = atoms->atom[a].resind;
    if (a==atoms->nr-1 || atoms->atom[a+1].resind != r) {
      hbexcl = &hb[r].rb[ebtsEXCLS];
      
      for(e=0; e<hbexcl->nb; e++) {
	anm = hbexcl->b[e].a[0];
	i1 = search_atom(anm,astart,atoms->nr,atoms->atom,atoms->atomname,
			 "exclusion",bAllowMissing);
	anm = hbexcl->b[e].a[1];
	i2 = search_atom(anm,astart,atoms->nr,atoms->atom,atoms->atomname,
			 "exclusion",bAllowMissing);
	if (i1!=NO_ATID && i2!=NO_ATID) {
	  if (i1 > i2) {
	    itmp = i1;
	    i1 = i2;
	    i2 = itmp;
	  }
	  srenew(excls[i1].e,excls[i1].nr+1);
	  excls[i1].e[excls[i1].nr] = i2;
	  excls[i1].nr++;
	}
      }
      
      astart = a+1;
    }
  }

  for(a=0; a<atoms->nr; a++)
    if (excls[a].nr > 1)
      qsort(excls[a].e,excls[a].nr,(size_t)sizeof(atom_id),atom_id_comp);
}

static void remove_excl(t_excls *excls, int remove)
{
  int i;

  for(i=remove+1; i<excls->nr; i++)
    excls->e[i-1] = excls->e[i];
  
  excls->nr--;
}

void clean_excls(t_nextnb *nnb, int nrexcl, t_excls excls[])
{
  int i,j,j1,k,k1,l,l1,m,n,e;
  t_excls *excl;

  if (nrexcl >= 1)
    /* extract all i-j-k-l neighbours from nnb struct */
    for(i=0; (i<nnb->nr); i++) {
      /* For all particles */
      excl = &excls[i];
      
      for(j=0; (j<nnb->nrexcl[i][1]); j++) {
	/* For all first neighbours */
	j1=nnb->a[i][1][j];
	
	for(e=0; e<excl->nr; e++)
	  if (excl->e[e] == j1)
	    remove_excl(excl,e);
	
	if (nrexcl >= 2)
	  for(k=0; (k<nnb->nrexcl[j1][1]); k++) {
	    /* For all first neighbours of j1 */
	    k1=nnb->a[j1][1][k];
	  
	    for(e=0; e<excl->nr; e++)
	      if (excl->e[e] == k1)
		remove_excl(excl,e);
	    
	    if (nrexcl >= 3)
	      for(l=0; (l<nnb->nrexcl[k1][1]); l++) {
		/* For all first neighbours of k1 */
		l1=nnb->a[k1][1][l];

		for(e=0; e<excl->nr; e++)
		  if (excl->e[e] == l1)
		    remove_excl(excl,e);
	      }
	  }
      }
    }
}

void generate_excls(t_nextnb *nnb, int nrexcl, t_excls excls[])
{
  int i,j,j1,k,k1,l,l1,m,n,e,N;
  t_excls *excl;

  for(N=1; (N<min(nrexcl,nnb->nrex)); N++) {
    /* extract all i-j-k-l neighbours from nnb struct */
    for(i=0; (i<nnb->nr); i++) {
      /* For all particles */
      excl = &excls[i];
      n = excl->nr;
      excl->nr += nnb->nrexcl[i][N];
      srenew(excl->e,excl->nr);
      for(j=0; (j<nnb->nrexcl[i][N]); j++) 
	/* For all first neighbours */
	if (nnb->a[i][N][j] != i)
	  excl->e[n++] = nnb->a[i][N][j];
    }
  }
}

void gen_pad(t_nextnb *nnb, t_atoms *atoms, int nrexcl, gmx_bool bH14,
	     t_params plist[], t_excls excls[], t_hackblock hb[], 
	     gmx_bool bAlldih, gmx_bool bRemoveDih, gmx_bool bAllowMissing)
{
  t_param *ang,*dih,*pai,*idih;
  t_rbondeds *hbang, *hbdih;
  char    **anm;
  int     res,minres,maxres;
  int     i,j,j1,k,k1,l,l1,m,n,i1,i2;
  int     ninc,maxang,maxdih,maxpai;
  int     nang,ndih,npai,nidih,nbd;
  int     nFound;
  gmx_bool    bFound,bExcl;
  

  /* These are the angles, dihedrals and pairs that we generate
   * from the bonds. The ones that are already there from the rtp file
   * will be retained.
   */
  nang   = 0;
  npai   = 0;
  ndih   = 0;
  ninc   = 500;
  maxang = maxdih = maxpai = ninc;
  snew(ang, maxang);
  snew(dih, maxdih);
  snew(pai, maxpai);

  snew(anm,4);
  for(i=0;i<4;i++)
    snew(anm[i],12);

  if (hb)
    gen_excls(atoms,excls,hb,bAllowMissing);
  
  /* extract all i-j-k-l neighbours from nnb struct */
  for(i=0; (i<nnb->nr); i++) 
    /* For all particles */
    for(j=0; (j<nnb->nrexcl[i][1]); j++) {
      /* For all first neighbours */
      j1=nnb->a[i][1][j];
      for(k=0; (k<nnb->nrexcl[j1][1]); k++) {
	/* For all first neighbours of j1 */
	k1=nnb->a[j1][1][k];
	if (k1 != i) {
	  /* Generate every angle only once */
	  if (i < k1) {
	    if (nang == maxang) {
	      maxang += ninc;
	      srenew(ang,maxang);
	    }
	    ang[nang].AI=i;
	    ang[nang].AJ=j1;
	    ang[nang].AK=k1;
	    ang[nang].C0=NOTSET;
	    ang[nang].C1=NOTSET;
	    set_p_string(&(ang[nang]),"");
	    if (hb) {
	      minres = atoms->atom[ang[nang].a[0]].resind;
	      maxres = minres;
	      for(m=1; m<3; m++) {
		minres = min(minres,atoms->atom[ang[nang].a[m]].resind);
		maxres = max(maxres,atoms->atom[ang[nang].a[m]].resind);
	      }
	      res = 2*minres-maxres;
	      do {
		res += maxres-minres;
		get_atomnames_min(3,anm,res,atoms,ang[nang].a);
		hbang=&hb[res].rb[ebtsANGLES];
		for(l=0; (l<hbang->nb); l++) {
		  if (strcmp(anm[1],hbang->b[l].AJ)==0) {
		    bFound=FALSE;
		    for (m=0; m<3; m+=2)
		      bFound=(bFound ||
			      ((strcmp(anm[m],hbang->b[l].AI)==0) &&
			       (strcmp(anm[2-m],hbang->b[l].AK)==0)));
		    if (bFound) {
		      set_p_string(&(ang[nang]),hbang->b[l].s);
		    }
		  }
		}
	      } while (res < maxres);
	    }
	    nang++;
	  }
	  /* Generate every dihedral, 1-4 exclusion and 1-4 interaction
	     only once */
	  if (j1 < k1) {
	    for(l=0; (l<nnb->nrexcl[k1][1]); l++) {
	      /* For all first neighbours of k1 */
	      l1=nnb->a[k1][1][l];
	      if ((l1 != i) && (l1 != j1)) {
		if (ndih == maxdih) {
		  maxdih += ninc;
		  srenew(dih,maxdih);
		}
		dih[ndih].AI=i;
		dih[ndih].AJ=j1;
		dih[ndih].AK=k1;
		dih[ndih].AL=l1;
		for (m=0; m<MAXFORCEPARAM; m++)
		  dih[ndih].c[m]=NOTSET;
		set_p_string(&(dih[ndih]),"");
		nFound = 0;
		if (hb) {
		  minres = atoms->atom[dih[ndih].a[0]].resind;
		  maxres = minres;
		  for(m=1; m<4; m++) {
		    minres = min(minres,atoms->atom[dih[ndih].a[m]].resind);
		    maxres = max(maxres,atoms->atom[dih[ndih].a[m]].resind);
		  }
		  res = 2*minres-maxres;
		  do {
		    res += maxres-minres;
		    get_atomnames_min(4,anm,res,atoms,dih[ndih].a);
		    hbdih=&hb[res].rb[ebtsPDIHS];
		    for(n=0; (n<hbdih->nb); n++) {
		      bFound=FALSE;
		      for (m=0; m<2; m++)
			bFound=(bFound ||
				((strcmp(anm[3*m],  hbdih->b[n].AI)==0) &&
				 (strcmp(anm[1+m],  hbdih->b[n].AJ)==0) &&
				 (strcmp(anm[2-m],  hbdih->b[n].AK)==0) &&
				 (strcmp(anm[3-3*m],hbdih->b[n].AL)==0)));
		      if (bFound) {
			set_p_string(&dih[ndih],hbdih->b[n].s);
			
			/* Set the last parameter to be able to see
			   if the dihedral was in the rtp list.
			   */
			dih[ndih].c[MAXFORCEPARAM-1] = 0;
			nFound++;
			ndih++;
			/* Set the next direct in case the rtp contains
			   multiple entries for this dihedral.
			   */
			if (ndih == maxdih) {
			  maxdih += ninc;
			  srenew(dih,maxdih);
			}
			dih[ndih].AI=i;
			dih[ndih].AJ=j1;
			dih[ndih].AK=k1;
			dih[ndih].AL=l1;
			for (m=0; m<MAXFORCEPARAM; m++)
			  dih[ndih].c[m]=NOTSET;
		      }
		    }
		  } while (res < maxres);
		}
		if (nFound == 0) {
		  if (ndih == maxdih) {
		    maxdih += ninc;
		    srenew(dih,maxdih);
		  }
		  dih[ndih].AI=i;
		  dih[ndih].AJ=j1;
		  dih[ndih].AK=k1;
		  dih[ndih].AL=l1;
		  for (m=0; m<MAXFORCEPARAM; m++)
		    dih[ndih].c[m]=NOTSET;
		  set_p_string(&(dih[ndih]),"");
		  ndih++;
		}

		nbd=nb_dist(nnb,i,l1);
		if (debug)
		  fprintf(debug,"Distance (%d-%d) = %d\n",i+1,l1+1,nbd);
		if (nbd == 3) {
		  i1 = min(i,l1);
		  i2 = max(i,l1);
		  bExcl = FALSE;
		  for(m=0; m<excls[i1].nr; m++)
		    bExcl = bExcl || excls[i1].e[m]==i2;
		  if (!bExcl) {
		    if (bH14 || !(is_hydro(atoms,i1) && is_hydro(atoms,i2))) {
		      if (npai == maxpai) {
			maxpai += ninc;
			srenew(pai,maxpai);
		      }
		      pai[npai].AI=i1;
		      pai[npai].AJ=i2;
		      pai[npai].C0=NOTSET;
		      pai[npai].C1=NOTSET;
		      set_p_string(&(pai[npai]),"");
		      npai++;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

  /* Sort angles with respect to j-i-k (middle atom first) */
  if (nang > 1)
    qsort(ang,nang,(size_t)sizeof(ang[0]),acomp);
  
  /* Sort dihedrals with respect to j-k-i-l (middle atoms first) */
  if (ndih > 1)
    qsort(dih,ndih,(size_t)sizeof(dih[0]),dcomp);
  
  /* Sort the pairs */
  if (npai > 1)
    qsort(pai,npai,(size_t)sizeof(pai[0]),pcomp);
  if (npai > 0) {
    /* Remove doubles, could occur in 6-rings, such as phenyls,
       maybe one does not want this when fudgeQQ < 1.
       */
    fprintf(stderr,"Before cleaning: %d pairs\n",npai);
    rm2par(pai,&npai,preq);
  }

  /* Get the impropers from the database */
  nidih = get_impropers(atoms,hb,&idih,bAllowMissing);

  /* Sort the impropers */
  sort_id(nidih,idih);
 
  if (ndih > 0) {
    /* Remove dihedrals which are impropers
       and when bAlldih is not set remove multiple dihedrals over one bond.
       */
    fprintf(stderr,"Before cleaning: %d dihedrals\n",ndih);
    clean_dih(dih,&ndih,idih,nidih,atoms,bAlldih,bRemoveDih);
  }

  /* Now we have unique lists of angles and dihedrals 
   * Copy them into the destination struct
   */
  cppar(ang, nang, plist,F_ANGLES);
  cppar(dih, ndih, plist,F_PDIHS);
  cppar(idih,nidih,plist,F_IDIHS);
  cppar(pai, npai, plist,F_LJ14);

  /* Remove all exclusions which are within nrexcl */
  clean_excls(nnb,nrexcl,excls);

  sfree(ang);
  sfree(dih);
  sfree(idih);
  sfree(pai);
}

