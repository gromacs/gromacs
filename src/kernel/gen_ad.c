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

static bool remove_dih(t_param *p, int i, int np)
     /* check if dihedral p[i] should be removed */
{
  bool bMidEq,bRem;
  int j;

  if (i>0)
    bMidEq = deq(&p[i],&p[i-1]);
  else
    bMidEq = FALSE;

  if (p[i].c[MAXFORCEPARAM-1]==NOTSET) {
    /* also remove p[i] if there is a dihedral on the same bond
       which has parameters set */
    bRem = bMidEq;
    j=i+1;
    while (!bRem && (j<np) && deq(&p[i],&p[j])) {
      bRem = (p[j].c[MAXFORCEPARAM-1] != NOTSET);
      j++;
    }
  } else
    bRem = bMidEq && (((p[i].AI==p[i-1].AI) && (p[i].AL==p[i-1].AL)) ||
		      ((p[i].AI==p[i-1].AL) && (p[i].AL==p[i-1].AI)));

  return bRem;
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
    if (p[index[i]].a[0] == p[index[i]].a[1]) {
      if (debug)  
	fprintf(debug,
		"Something VERY strange is going on in rm2par (gen_ad.c)\n"
		"a[0] %d a[1] %d a[2] %d a[3] %d\n",
		p[i].a[0],p[i].a[1],p[i].a[2],p[i].a[3]);
      p[i].s = strdup(""); 
    } else {
      sfree(p[i].s);
      p[i].s = strdup(p[index[i]].s);
    }
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
    ps->param[i].s=strdup(p[i].s);
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
  dest->s=strdup(src->s);
}

static void set_p(t_param *p,atom_id ai[4],real *c,char *s)
{
  int j;

  for(j=0; (j<4); j++)
    p->a[j]=ai[j];
  for(j=0; (j<MAXFORCEPARAM); j++)
    p->c[j]=c[j];
  p->s=strdup(s);
}

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
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
  qsort(ps,nr,(size_t)sizeof(ps[0]),idcomp);
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
	aa0=search_atom(atom,start,atoms->nr,atoms->atom,atoms->atomname);
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

static void pdih2idih(t_param *alldih,int *nalldih,t_param idih[],int *nidih,
		      t_atoms *atoms,int nrtp,t_restp rtp[],
		      int nrdh,t_idihres idh[],bool bAlldih)
{
  t_param   *dih,tmp_param;
  int       ndih;
  char      *rname,*a0;
  t_idihres *i0;
  int       i,j,k,l,start,aa0;
  int       *index,nind;
  atom_id   ai[MAXATOMLIST];
  bool      bIsSet,bKeep;
  int bestl,nh,minh;

  /* First add all the impropers from the residue database
   * to the list.
   */
  start=0;
  for(i=0; (i<atoms->nres); i++) {
    rname=*(atoms->resname[i]);
    if (search_rtp(rname,nrtp,rtp) == NULL) 
      fatal_error(0,"Residue %s not in residue database\n",rname);
    else if ((i0=search_idih(rname,nrdh,idh)) != NULL) {
      for(j=0; (j<i0->nidih); j++) {
	for(k=0; (k<4); k++) {
	  a0=i0->idih[j].ai[k];
	  aa0=search_atom(a0,start,atoms->nr,atoms->atom,atoms->atomname);
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
	  set_p(&(idih[*nidih]),ai,i0->idih[j].c,i0->idih[j].s);
	  (*nidih)++;
	}
      }
    }
    while ((start<atoms->nr) && (atoms->atom[start].resnr==i))
      start++;
  }

  if (*nalldih == 0)
    return;

  /* Copy the impropers and dihedrals to seperate arrays. */
  snew(dih,*nalldih);
  ndih = 0;
  for(i=0; i<*nalldih; i++) 
    if (is_imp(&(alldih[i]),atoms,nrdh,idh)) {
      cpparam(&(idih[*nidih]),&(alldih[i]));
      (*nidih)++;
    } else {
      cpparam(&(dih[ndih]),&(alldih[i]));
      ndih++;
    }
  
  /* Now, because this list still contains the double entries,
   * keep the dihedral with parameters or the first one.
   */
    
  snew(index,ndih);
  nind=0;
  if (bAlldih) {
    fprintf(stderr,"bAlldih = true\n");
    for(i=0; i<ndih; i++) 
      if ((i==0) || !deq2(&dih[i],&dih[i-1]))
	index[nind++]=i;
  } else {
    for(i=0; i<ndih; i++) 
      if (!remove_dih(dih,i,ndih)) 
	index[nind++]=i;
  }
  index[nind]=ndih;

  /* if we don't want all dihedrals, we need to select the ones with the 
   *  fewest hydrogens
   */
  
  k=0;
  for(i=0; i<nind; i++) {
    bIsSet = (dih[index[i]].c[MAXFORCEPARAM-1] != NOTSET);
    bKeep = TRUE;
    if (!bIsSet)
      /* remove the dihedral if there is an improper on the same bond */
      for(j=0; (j<(*nidih)) && bKeep; j++)
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
      for(j=0; (j<MAXATOMLIST); j++)
	alldih[k].a[j] = dih[bestl].a[j];
      for(j=0; (j<MAXFORCEPARAM); j++)
	alldih[k].c[j] = dih[bestl].c[j];
      sfree(alldih[k].s);
      alldih[k].s = strdup(dih[bestl].s);
      k++;
    }
  }
  for (i=(*nalldih); i<k; i++)
    sfree(alldih[i].s);
  *nalldih = k;

  for(i=0; i<ndih; i++)
    sfree(dih[i].s);
  sfree(dih);
  sfree(index);
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

static void get_atomnames_min(int n,char anm[4][12],
			      int res,t_atoms *atoms,atom_id *a)
{
  int m;

  /* Assume ascending residue numbering */
  for(m=0; m<n; m++) {
    if (atoms->atom[a[m]].resnr < res)
      strcpy(anm[m],"-");
    else if (atoms->atom[a[m]].resnr > res)
      strcpy(anm[m],"+");
    else
      strcpy(anm[m],"");
    strcat(anm[m],*(atoms->atomname[a[m]]));
  }
}

void gen_pad(t_nextnb *nnb,t_atoms *atoms,bool bH14,t_params plist[],
	     int nrtp,t_restp rtp[],
	     int nra,t_resang ra[],int nrd,t_resdih rd[], 
	     int nid,t_idihres idihs[],bool bAlldih)
{
  t_param *ang,*dih,*pai,*idih;
  t_resang *i_ra;
  t_resdih *i_rd;
  char    anm[4][12];
  int     res,minres,maxres;
  int     i,j,j1,k,k1,l,l1,m,n;
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
	  ang[nang].s=strdup("");
	  minres = atoms->atom[ang[nang].a[0]].resnr;
	  maxres = minres;
	  for(m=1; m<3; m++) {
	    minres = min(minres,atoms->atom[ang[nang].a[m]].resnr);
	    maxres = max(maxres,atoms->atom[ang[nang].a[m]].resnr);
	  }
	  res = 2*minres-maxres;
	  do {
	    res += maxres-minres;
	    if ((i_ra=search_rang(*(atoms->resname[res]),nra,ra))) {
	      for(l=0; (l<i_ra->na); l++) {
		get_atomnames_min(3,anm,res,atoms,ang[nang].a); 
		if (strcmp(anm[1],i_ra->rang[l].aj)==0) {
		  bFound=FALSE;
		  for (m=0; m<3; m+=2)
		    bFound=(bFound ||
			    ((strcmp(anm[m],i_ra->rang[l].ai)==0) &&
			     (strcmp(anm[2-m],i_ra->rang[l].ak)==0)));
		  if (bFound) { 
		    for (m=0; m<MAXFORCEPARAM; m++)
		      ang[nang].c[m] = i_ra->rang[l].c[m];
		    sfree(ang[nang].s);
		    ang[nang].s = strdup(i_ra->rang[l].s);
		  }
		}
	      }
	    }
	  } while (res < maxres);
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
	      for (m=0; m<MAXFORCEPARAM; m++)
		dih[ndih].c[m]=NOTSET;
	      dih[ndih].s=strdup("");
	      minres = atoms->atom[dih[ndih].a[0]].resnr;
	      maxres = minres;
	      for(m=1; m<4; m++) {
		minres = min(minres,atoms->atom[dih[ndih].a[m]].resnr);
		maxres = max(maxres,atoms->atom[dih[ndih].a[m]].resnr);
	      }
	      res = 2*minres-maxres;
	      do {
		res += maxres-minres;
		if ((i_rd=search_rdih(*(atoms->resname[res]),nrd,rd))) {
		  for(n=0; (n<i_rd->nd); n++) {
		    get_atomnames_min(4,anm,res,atoms,dih[ndih].a);
		    bFound=FALSE;
		    for (m=0; m<2; m++)
		      bFound=(bFound ||
			      ((strcmp(anm[3*m],  i_rd->rdih[n].ai)==0) &&
			       (strcmp(anm[1+m],  i_rd->rdih[n].aj)==0) &&
			       (strcmp(anm[2-m],  i_rd->rdih[n].ak)==0) &&
			       (strcmp(anm[3-3*m],i_rd->rdih[n].al)==0)));
		    if (bFound) {
		      for (m=0; m<MAXFORCEPARAM-1; m++)
			dih[ndih].c[m] = i_rd->rdih[n].c[m];
		      sfree(dih[ndih].s);
		      dih[ndih].s = strdup(i_rd->rdih[n].s);
		      /* Set the last parameter to be able to see
			 if the dihedral was in the rtp list */
		      dih[ndih].c[MAXFORCEPARAM-1] = 0;
		    }
		  }
		}
	      } while (res < maxres);
	      nbd=nb_dist(nnb,i,l1);
	      if (debug)
		fprintf(debug,"Distance (%d-%d) = %d\n",i+1,l1+1,nbd);
	      if (nbd == 3) {
		pai[npai].AI=min(i,l1);
		pai[npai].AJ=max(i,l1);
		pai[npai].C0=NOTSET;
		pai[npai].C1=NOTSET;
		pai[npai].s =strdup("");
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
  qsort(ang,nang,(size_t)sizeof(ang[0]),acomp);
  rm2par(ang,&nang,aeq);

  /* Sort dihedrals with respect to j-k-i-l (middle atoms first) */
  fprintf(stderr,"before sorting: %d dihedrals\n",ndih);
  qsort(dih,ndih,(size_t)sizeof(dih[0]),dcomp);
  pdih2idih(dih,&ndih,idih,&nidih,atoms,nrtp,rtp,nid,idihs,bAlldih);
  fprintf(stderr,"after sorting: %d dihedrals\n",ndih);
  
  /* Now the dihedrals are sorted and doubles removed, this has to be done
   * for impropers too
   */
  sort_id(nidih,idih);
  rm2par(idih,&nidih,ideq);
  
  /* And for the pairs */
  fprintf(stderr,"There are %d pairs before sorting\n",npai);
  qsort(pai,npai,(size_t)sizeof(pai[0]),pcomp);
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

