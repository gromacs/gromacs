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
 * GROwing Monsters And Cloning Shrimps
 */
static char *SRCID_gen_ad_c = "$Id$";

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
#include "topexcl.h"
#include "symtab.h"
#include "macros.h"
#include "fatal.h"
#include "pgutil.h"
#include "resall.h"
#include "pdb2gmx.h"

typedef bool (*peq)(t_param *p1, t_param *p2);

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
  if ((dc=(p1->AJ-p2->AJ))!=0)
    return dc;
  else if ((dc=(p1->AK-p2->AK))!=0)
    return dc;
  else if ((dc=(p1->AI-p2->AI))!=0)
    return dc;
  else
    return (p1->AL-p2->AL);
}

static bool aeq(t_param *p1, t_param *p2)
{
  if (p1->AJ!=p2->AJ) 
    return FALSE;
  else if (((p1->AI==p2->AI) && (p1->AK==p2->AK)) ||
	   ((p1->AI==p2->AK) && (p1->AK==p2->AI)))
    return TRUE;
  else 
    return FALSE;
}

static bool deq2(t_param *p1, t_param *p2)
{
  /* if bAlldih is true, dihedrals are only equal when 
     ijkl = ijkl or ijkl =lkji*/
  if (((p1->AI==p2->AI) && (p1->AJ==p2->AJ) && 
       (p1->AK==p2->AK) && (p1->AL==p2->AL)) ||
      ((p1->AI==p2->AL) && (p1->AJ==p2->AK) &&
       (p1->AK==p2->AJ) && (p1->AL==p2->AI)))
    return TRUE;
  else 
    return FALSE;
}

static bool deq(t_param *p1, t_param *p2)
{
  if (((p1->AJ==p2->AJ) && (p1->AK==p2->AK)) ||
      ((p1->AJ==p2->AK) && (p1->AK==p2->AJ)))
    return TRUE;
  else 
    return FALSE;
}

static bool preq(t_param *p1, t_param *p2)
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
  }
  (*np)=nind;

  sfree(index);
}

static void cppar(t_param p[], int np, t_params plist[], int ftype)
{
  int      i,j,nral,nrfp;
  t_params *ps;

  ps   = &(plist[ftype]);
  nral = NRAL(ftype);
  nrfp = NRFP(ftype);
  
  /* Keep old stuff */
  srenew(ps->param,np);
  for(i=ps->nr; (i<np); i++) {
    for(j=0; (j<nral); j++)
      ps->param[i].a[j]=p[i].a[j];
    for(j=0; (j<nrfp); j++)
      ps->param[i].c[j]=p[i].c[j];
  }
  ps->nr=np;
}

static void cpparam(t_param *dest,t_param *src)
{
  int j;

  for(j=0; (j<MAXATOMLIST); j++)
    dest->a[j]=src->a[j];
  for(j=0; (j<MAXFORCEPARAM); j++)
    dest->c[j]=src->c[j];
}

static void set_p(t_param *p,atom_id ai[4],real *c)
{
  int j;

  for(j=0; (j<4); j++)
    p->a[j]=ai[j];
  for(j=0; (j<MAXFORCEPARAM); j++)
    p->c[j]=c[j];
}

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

static int eq_imp(atom_id a1[],atom_id a2[])
{
  int b1[MAXATOMLIST],b2[MAXATOMLIST];
  int j;

  for(j=0; (j<MAXATOMLIST); j++) {
    b1[j]=a1[j];
    b2[j]=a2[j];
  }
  qsort(b1,4,sizeof(b1[0]),int_comp);
  qsort(b2,4,sizeof(b2[0]),int_comp);

  for(j=0; (j<4); j++)
    if (b1[j] != b2[j])
      return FALSE;

  return TRUE;
}

static bool ideq(t_param *p1, t_param *p2)
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
  /* int i,j;
     atom_id a[4];
     
     
  for(i=0; (i<nr); i++) {
    if (ps[i].AL < ps[i].AI) {
      for(j=0; (j<4); j++)
	a[j]=ps[i].a[j];
      for(j=0; (j<4); j++)
	ps[i].a[j]=a[3-j];
    }
  }*/
  qsort(ps,nr,sizeof(ps[0]),idcomp);
  
}

static bool is_imp(t_param *p,t_atoms *atoms,int nrdh,t_idihres idih[])
{
  int        j,n,pm,start;
  atom_id    a0[MAXATOMLIST];
  int        aa0;
  char      *atom;
  t_idihres *i0;

  /* Find the max residue number in this dihedral */
  pm=0;
  for(j=0; (j<4); j++)
    pm=max(pm,atoms->atom[p->a[j]].resnr);

  /* Now find the start of this residue */
  for(start=0; (start < atoms->nr); start++)
    if (atoms->atom[start].resnr == pm)
      break;

  /* See if there are any impropers defined for this residue 
   * Impropers are always defined backwards (that is only referring
   * to previous residues, not to the next!)
   */
  i0=search_idih(*(atoms->resname[atoms->atom[start].resnr]),nrdh,idih);

  if (i0 != NULL) {
    for(n=0; (n<i0->nidih); n++) {
      for(j=0; (j<4); j++) {
	atom=i0->idih[n].ai[j];
	aa0=search_atom(atom,start,atoms->nr,atoms->atomname);
	if (aa0 == -1) {
	  if (debug) 
	    fprintf(debug,"Atom %s not found in res %d (pm=%d)\n",atom,
		    atoms->atom[start].resnr,pm);
	  break;
	}
	else 
	  a0[j] = aa0;
      }
      if (j==4) /* Not broken out */
	if (eq_imp(p->a,a0))
	  return TRUE;
    }
  }
  return FALSE;
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

static void pdih2idih(t_param dih[],int *ndih,t_param idih[],int *nidih,
		      t_atoms *atoms,int nrtp,t_restp rtp[],
		      int nrdh,t_idihres idh[],bool bAlldih)
{
  char      *rname,*a0;
  t_idihres *i0;
  t_restp   *r0;
  int       i,j,k,l,start,aa0,nral;
  int       *index,nind;
  atom_id   ai[MAXATOMLIST];
  bool      *bRM,brm;
  
  /* First add all the impropers from the residue database
   * to the list.
   */
  nral = NRAL(F_IDIHS);
  start=0;
  for(i=0; (i<atoms->nres); i++) {
    rname=*(atoms->resname[i]);
    if ((r0=search_rtp(rname,nrtp,rtp)) == NULL) 
      fatal_error(0,"Residue %s not in residue database\n",rname);
    else if ((i0=search_idih(rname,nrdh,idh)) != NULL) {
      for(j=0; (j<i0->nidih); j++) {
	for(k=0; (k<4); k++) {
	  a0=i0->idih[j].ai[k];
	  aa0=search_atom(a0,start,atoms->nr,atoms->atomname);
	  if (aa0 == -1) {
	    if (debug) 
	      fprintf(debug,"Atom %s not found in res %d\n",a0,
		      atoms->atom[start].resnr);
	    break;
	  }
	  else
	    ai[k]=aa0;
	}
	if (k==4) {
	  /* Not broken out */
	  set_p(&(idih[*nidih]),ai,i0->idih[j].c);
	  (*nidih)++;
	}
      }
    }
    while ((start<atoms->nr) && (atoms->atom[start].resnr==i))
      start++;
  }

  if (*ndih == 0)
    return;

  /* Check whether a generated dihedral really is an improper
   * and if so mark it to be removed. */
  snew(bRM,*ndih);
  for(i=0; (i<(*ndih)); i++) 
    if (is_imp(&(dih[i]),atoms,nrdh,idh)) {
      cpparam(&(idih[*nidih]),&(dih[i]));
      (*nidih)++;
      bRM[i]=TRUE;
    }
  /* Now, because this list still contains the double entries,
   * we do the removing of doubles here together with a check on
   * the bRM array, filled just above.
   */
    
  snew(index,*ndih);
  nind=0;
  index[nind++]=0;
  for(i=1; (i<(*ndih)); i++) 
    if (bAlldih)
    {
      fprintf(stderr,"bAlldih = true\n");
      if (!deq2(&dih[i],&dih[i-1]))
	index[nind++]=i;
    } else {
      if (!deq(&dih[i],&dih[i-1]))
	index[nind++]=i;
    }
  /* Index now holds pointers to all the non-equal params,
   * this only works when dih is sorted of course
   */
  for(i=0; (i<nind); i++) {
    brm=FALSE;
    /* Now check if one of the (equal) dihedrals must go */
    for(j=index[i]; (j<index[i+1]); j++)
      brm=brm||bRM[j];
    /* Then all must go */
    for(j=index[i]; (j<index[i+1]); j++)
      bRM[j]=brm;
  }
  
  /* if we don't want all dihedrals, we need to select the ones with the 
   *  fewest hydrogens
   */
  for(i=k=0; (i<nind); i++)
    if (!bRM[index[i]]) {
      /* Now select the "fittest" dihedral:
       * the one with as few hydrogens as possible 
       */
      int bestl,nh,minh;

      /* Best choice to get dihedral from */
      bestl=index[i];
      /* Minimum number of hydrogens for i and l atoms */
      if (!bAlldih)
      {
	minh=2;
	for(l=index[i]; (l<index[i+1]); l++) {
	  if ((nh=n_hydro(dih[l].a,atoms->atomname)) < minh) {
	    minh=nh;
	  bestl=l;
	  }
	  if (minh == 0)
	    break;
	}
      }
      for(j=0; (j<MAXATOMLIST); j++)
	dih[k].a[j]=dih[bestl].a[j];
      for(j=0; (j<MAXFORCEPARAM); j++)
	dih[k].c[j]=dih[bestl].c[j];
      k++;
    }
  (*ndih)=k;
  
  sfree(index);
  sfree(bRM);
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

bool is_hydro(t_atoms *atoms,int ai)
{
  return ((*(atoms->atomname[ai]))[0] == 'H');
}

void gen_pad(t_nextnb *nnb,t_atoms *atoms,bool bH14,t_params plist[],
	     int nrtp,t_restp rtp[],int nra,t_resang ra[],
	     int nid,t_idihres idihs[],bool bAlldih)
{
  t_param *ang,*dih,*pai,*idih;
  t_resang *i0;
  int     i,j,j1,k,k1,l,l1,m;
  int     maxang,maxdih,maxidih,maxpai;
  int     nang,ndih,npai,nidih,nbd;
  bool    bFound;

  nang    = plist[F_ANGLES].nr;
  nidih   = plist[F_IDIHS].nr;
  ndih    = plist[F_PDIHS].nr;
  npai    = plist[F_LJ14].nr;
  maxang  = 6*nnb->nr;
  maxdih  = 24*nnb->nr;
  maxpai  = maxdih;
  maxidih = maxdih;
  snew(ang,maxang);
  snew(dih,maxdih);
  snew(pai,maxpai);
  snew(idih,maxidih);

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
	  if (nang == maxang)
	    fatal_error(0,"Too many angles (%d) generated!\n",nang);
	  ang[nang].AJ=j1;
	  if (i < k1) {
	    ang[nang].AI=i;
	    ang[nang].AK=k1;
	  }
	  else {
	    ang[nang].AI=k1;
	    ang[nang].AK=i;
	  }
	  ang[nang].C0=NOTSET;
	  ang[nang].C1=NOTSET;
	  if (i0=search_rang(*(atoms->resname[atoms->atom[j1].resnr]),
			    nra,ra)) {
	    for(l=0; (l<i0->na); l++) {
	      if (strcmp(*(atoms->atomname[ang[nang].a[1]]),
			 i0->rang[l].aj)==0) {
		bFound=FALSE;
		for (m=0; m<3; m+=2)
		  bFound=(bFound ||
			  ((strcmp(*(atoms->atomname[ang[nang].a[m]]),
				   i0->rang[l].ai)==0) &&
			   (strcmp(*(atoms->atomname[ang[nang].a[2-m]]),
				   i0->rang[l].ak)==0)));
		if (bFound) 
		  for (m=0; m<MAXFORCEPARAM; m++)
		    ang[nang].c[m]=i0->rang[l].c[m];
	      }
	    }
	  }
	  nang++;
	  for(l=0; (l<nnb->nrexcl[k1][1]); l++) {
	    /* For all first neighbours of k1 */
	    l1=nnb->a[k1][1][l];
	    if ((l1 != i) && (l1 != j1)) {
	      if (ndih == maxdih)
		fatal_error(0,"Too many dihedrals (%d) generated!\n",ndih);
	      if (j1 < k1) {
		dih[ndih].AI=i;
		dih[ndih].AJ=j1;
		dih[ndih].AK=k1;
		dih[ndih].AL=l1;
	      }
	      else {
		dih[ndih].AI=l1;
		dih[ndih].AJ=k1;
		dih[ndih].AK=j1;
		dih[ndih].AL=i;
	      }
	      dih[ndih].C0=NOTSET;
	      dih[ndih].C1=NOTSET;
	      dih[ndih].C2=NOTSET;
	      nbd=nb_dist(nnb,i,l1);
	      if (debug)
		fprintf(debug,"Distance (%d-%d) = %d\n",i+1,l1+1,nbd);
	      if (nbd == 3) {
		pai[npai].AI=min(i,l1);
		pai[npai].AJ=max(i,l1);
		pai[npai].C0=NOTSET;
		pai[npai].C1=NOTSET;
		if (bH14 || !(is_hydro(atoms,pai[npai].AI) &&
			      is_hydro(atoms,pai[npai].AJ)))
		  npai++;
	      }
	      ndih++;
	    }
	  }
	}
      }
    }
  
  /* We now have a params list with double entries for each angle,
   * and even more for each dihedral. We will remove these now.
   */
  /* Sort angles with respect to j-i-k (middle atom first) */
  qsort(ang,nang,sizeof(ang[0]),acomp);
  rm2par(ang,&nang,aeq);

  /* Sort dihedrals with respect to j-k-i-l (middle atoms first) */
  fprintf(stderr,"before sorting: %d dihedrals\n",ndih);
  qsort(dih,ndih,sizeof(dih[0]),dcomp);
  pdih2idih(dih,&ndih,idih,&nidih,atoms,nrtp,rtp,nid,idihs,bAlldih);
  fprintf(stderr,"after sorting: %d dihedrals\n",ndih);
  
  /* Now the dihedrals are sorted and doubles removed, this has to be done
   * for impropers too
   */
  sort_id(nidih,idih);
  rm2par(idih,&nidih,ideq);
  
  /* And for the pairs */
  fprintf(stderr,"There are %d pairs before\n",npai);
  qsort(pai,npai,sizeof(pai[0]),pcomp);
  rm2par(pai,&npai,preq);

  /* Now we have unique lists of angles and dihedrals 
   * Copy them into the destination struct
   */
  cppar(ang, nang, plist,F_ANGLES);
  cppar(dih, ndih, plist,F_PDIHS);
  cppar(idih,nidih,plist,F_IDIHS);
  cppar(pai, npai, plist,F_LJ14);

  sfree(ang);
  sfree(dih);
  sfree(idih);
  sfree(pai);
}

