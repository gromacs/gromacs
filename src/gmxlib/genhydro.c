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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_genhydro_c = "$Id$";

#include <time.h>
#include <ctype.h>
#include "assert.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "string2.h"
#include "confio.h"
#include "symtab.h"
#include "vec.h"
#include "statutil.h"
#include "futil.h"
#include "fatal.h"
#include "physics.h"
#include "calch.h"
#include "genhydro.h"
#include "h_db.h"

static void copy_atom(t_atoms *atoms1,t_atoms *atoms2,rvec x1[],rvec x2[],
		      int a1,int a2)
{
  atoms2->atom[a2] = atoms1->atom[a1];
  snew(atoms2->atomname[a2],1);
  *atoms2->atomname[a2]=strdup(*atoms1->atomname[a1]);
  copy_rvec(x1[a1],x2[a2]);
}

int pdbasearch_atom(char *name,int resnr,t_atoms *pdba)
{
  int  i;
  
  if (name[0] == '-') {
    name++;
    resnr--;
  }
  for(i=0; (i<pdba->nr) && (pdba->atom[i].resnr != resnr); i++)
    ;
  for(   ; (i<pdba->nr) && (pdba->atom[i].resnr == resnr); i++) {
    if (strcmp(name,*pdba->atomname[i]) == 0)
      return i;
  }
  return -1;
}

void dump_ab(FILE *out,int natom,int nab[], t_hack *ab[], bool bHeader)
{
  int i,j;
  
#define SS(s) (s)?(s):"-"
  /* dump ab */
  if (bHeader)
    fprintf(out,"ADDBLOCK (t_hack)\n"
	    "%4s %2s %-4s %-4s %2s %-4s %-4s %-4s %-4s %1s %s\n",
	    "atom","nr","old","new","tp","ai","aj","ak","al","a","x");
  for(i=0; i<natom; i++)
    for(j=0; j<nab[i]; j++)
      fprintf(out,"%4d %2d %-4s %-4s %2d %-4s %-4s %-4s %-4s %s %g %g %g\n",
	      i+1, ab[i][j].nr, SS(ab[i][j].oname), SS(ab[i][j].nname),
	      ab[i][j].tp, 
	      SS(ab[i][j].AI), SS(ab[i][j].AJ),
	      SS(ab[i][j].AK), SS(ab[i][j].AL),
	      ab[i][j].atom?"+":"", 
	      ab[i][j].newx[XX], ab[i][j].newx[YY], ab[i][j].newx[ZZ]);
#undef SS
}

static t_hackblock *get_hackblocks(t_atoms *pdba, int nah, t_hackblock ah[],
				   t_hackblock *ntdb, t_hackblock *ctdb, 
				   int rN, int rC)
{
  int rnr;
  t_hackblock *hb,*ahptr;

  /* make space */
  snew(hb,pdba->nres);
  /* first the termini */
  if ( rN >= 0 )
    copy_t_hackblock(ntdb, &hb[rN]);
  if ( rC >= 0 )
    copy_t_hackblock(ctdb, &hb[rC]);
  /* then the whole hdb */
  for(rnr=0; rnr < pdba->nres; rnr++) {
    ahptr=search_h_db(nah,ah,*pdba->resname[rnr]);
    if ( ahptr ) {
      if (hb[rnr].name==NULL)
	hb[rnr].name=strdup(ahptr->name);
      merge_hacks(ahptr, &hb[rnr]);
    }
  }
  return hb;
}

static char Hnum[] = "123456";

static void expand_hackblocks_one(t_hackblock *hbr, char *atomname, 
				  int *nabi, t_hack **abi, bool bN, bool bC)
{
  int j, k, l, d;
  bool bIgnore;
  /* recursion depth is recorded for debug purposes only: */
  static int depth=-1;
  
  depth++;
  /* we'll recursively add atoms to atoms */
  if (debug) fprintf(debug,"\n[%d] %s:",depth,atomname);
  for(j=0; j < hbr->nhack; j++) {
    /* first check if we're in the N- or C-terminus, then we should ignore 
       all hacks involving atoms from resp. previous or next residue
       (i.e. which name begins with '-' (N) or '+' (C) */
    bIgnore = FALSE;
    if ( bN ) /* N-terminus: ignore '-' */
      for(k=0; k<4 && hbr->hack[j].a[k] && !bIgnore; k++)
	bIgnore = hbr->hack[j].a[k][0]=='-';
    if ( bC ) /* C-terminus: ignore '+' */
      for(k=0; k<4 && hbr->hack[j].a[k] && !bIgnore; k++)
	bIgnore = hbr->hack[j].a[k][0]=='+';
    /* must be either hdb entry (tp>0) or add from tdb (oname==NULL)
       and first control aton (AI) matches this atom or
       delete/replace from tdb (oname!=NULL) and oname matches this atom */
    if (debug) fprintf(debug," %s",
		       hbr->hack[j].oname?hbr->hack[j].oname:hbr->hack[j].AI);
    if ( !bIgnore && 
	 ( ( ( hbr->hack[j].tp > 0 || hbr->hack[j].oname==NULL ) &&
	     strcmp(atomname, hbr->hack[j].AI) == 0 ) ||
	   ( hbr->hack[j].oname!=NULL && 
	     strcmp(atomname, hbr->hack[j].oname) == 0) ) ) {
      /* now expand all hacks for this atom */
      if (debug) fprintf(debug," +%dh",hbr->hack[j].nr);
      srenew(*abi,*nabi + hbr->hack[j].nr);
      for(k=0; k < hbr->hack[j].nr; k++) {
	copy_t_hack(&hbr->hack[j], &(*abi)[*nabi + k]);
	for(d=0; d<DIM; d++)
	  (*abi)[*nabi + k].newx[d]=NOTSET;
	/* if we're adding (oname==NULL) and don't have a new name (nname) 
	   yet, build it from atomname */
	if ( (*abi)[*nabi + k].nname==NULL ) {
	  if ( (*abi)[*nabi + k].oname==NULL ) {
	    (*abi)[*nabi + k].nname=strdup(atomname);
	    (*abi)[*nabi + k].nname[0]='H';
	  }
	} else {
	  sfree((*abi)[*nabi + k].nname);
	  (*abi)[*nabi + k].nname=strdup(hbr->hack[j].nname);
	}
	/* if adding more than one atom, number them */
	if ( hbr->hack[j].nr > 1 ) {
	  l = strlen((*abi)[*nabi + k].nname);
	  srenew((*abi)[*nabi + k].nname, l+2);
	  (*abi)[*nabi + k].nname[l] = Hnum[k]; /* 1, 2, 3 .... */
	  (*abi)[*nabi + k].nname[l+1] = '\0';
	}
      }
      (*nabi) += hbr->hack[j].nr;
      
      /* add hacks to atoms we've just added */
      if ( hbr->hack[j].tp > 0 || hbr->hack[j].oname==NULL )
	for(k=0; k < hbr->hack[j].nr; k++)
	  expand_hackblocks_one(hbr, (*abi)[*nabi-hbr->hack[j].nr+k].nname, 
				nabi, abi, bN, bC);
    }
  }
  depth--;
}

static void expand_hackblocks(t_atoms *pdba, t_hackblock hb[], 
			      int nab[], t_hack *ab[], int rN, int rC)
{
  int i;
  
  for(i=0; i < pdba->nr; i++)
    /* add hacks to this atom */
    expand_hackblocks_one(&hb[pdba->atom[i].resnr], *pdba->atomname[i], 
			  &nab[i], &ab[i], 
			  pdba->atom[i].resnr==rN, 
			  pdba->atom[i].resnr==rC);
  if (debug) fprintf(debug,"\n");
}

static int check_atoms_present(t_atoms *pdba, int nab[], t_hack *ab[])
{
  int i, j, k, d, rnr, nadd;
  
  nadd=0;
  for(i=0; i < pdba->nr; i++) {
    rnr = pdba->atom[i].resnr;
    for(j=0; j<nab[i]; j++)
      if ( ab[i][j].oname==NULL ) { 
	/* we're adding */
	assert(ab[i][j].nname!=NULL);
	/* check if the atom is already present */
	k=pdbasearch_atom(ab[i][j].nname, rnr, pdba);
	if ( k != -1 ) {
	  /* we found the added atom, so move the hack there: */
	  srenew(ab[k], nab[k]+1);
	  ab[k][nab[k]] = ab[i][j];
	  ab[k][nab[k]].oname = strdup(ab[k][nab[k]].nname);
	  /* reset any possible new coordinates: */
	  for(d=0; d<DIM; d++)
	    ab[k][nab[k]].newx[d]=NOTSET;
	  /* keep count */
	  nab[k]++;
	  /* remove the hack from this atom: */
	  for(k=j+1; k<nab[i]; k++)
	    ab[i][k-1] = ab[i][k];
	  /* keep count */
	  nab[i]--;
	  j--;
	  srenew(ab[i], nab[i]);
	} else
	  /* count how many atoms we'll add */
	  nadd++;
      } else if ( ab[i][j].nname==NULL )
	/* we're deleting */
	nadd--;
  }
  
  return nadd;
}

static void calc_all_pos(t_atoms *pdba, rvec x[], int nab[], t_hack *ab[])
{
  int i, j, m, ia, d, rnr;
  rvec xa[4]; /* control atoms for calc_h_pos */
  rvec xh[3]; /* hydrogen positions from calc_h_pos */
  
  for(i=0; i < pdba->nr; i++) {
    rnr   = pdba->atom[i].resnr;
    for(j=0; j < nab[i]; j+=ab[i][j].nr) {
      /* check if we're adding: */
      if (ab[i][j].oname==NULL && ab[i][j].tp > 0) {
	for(m=0; m<ncontrol[ab[i][j].tp]; m++) {
	  ia = pdbasearch_atom(ab[i][j].a[m], rnr, pdba);
	  if (ia < 0)
	    fatal_error(0,"Atom %s not found in residue %s%d"
			" while adding hydrogens",
			ab[i][j].a[m],*pdba->resname[rnr],rnr+1);
	  copy_rvec(x[ia], xa[m]);
	}
	for(m=0; m<3; m++)
	  for(d=0; d<DIM; d++)
	    if (m<ab[i][j].nr)
	      xh[m][d] = 0;
	    else
	      xh[m][d] = NOTSET;
	calc_h_pos(ab[i][j].tp, xa, xh);
	for(m=0; m<ab[i][j].nr; m++)
	  copy_rvec(xh[m],ab[i][j+m].newx);
      }
    }
  }
}

int add_h(t_atoms **pdbaptr, rvec *xptr[], int nah, t_hackblock ah[],
	  t_hackblock *ntdb, t_hackblock *ctdb, int rN, int rC)
{
  t_atoms     *newpdba,*pdba;
  bool        bSet;
  int         nadd;
  int         i,newi,j,d,natom;
  int         *nab;
  t_hack      **ab;
  t_hackblock *hb;
  rvec        *xn;

  /* set flags for adding hydrogens (accoring to hdb) */
  pdba=*pdbaptr;
  natom=pdba->nr;
  /* first get all the hackblocks for each residue: */
  hb = get_hackblocks(pdba, nah, ah, ntdb, ctdb, rN, rC);
  if (debug) dump_hb(debug, pdba->nres, hb);
  
  /* expand the hackblocks to atom level */
  snew(nab,natom);
  snew(ab,natom);
  expand_hackblocks(pdba, hb, nab, ab, rN, rC);
  free_t_hackblock(pdba->nres, &hb);
  
  if (debug) dump_ab(debug, natom, nab, ab, TRUE);
  
  /* Now calc the positions */
  calc_all_pos(pdba, *xptr, nab, ab);
  
  /* we don't have to add atoms that are already present in pdba,
     so we will remove them from the ab (t_hack) */
  nadd = check_atoms_present(pdba, nab, ab);
  if (debug) {
    fprintf(debug, "removed add hacks that were already in pdba:\n");
    dump_ab(debug, natom, nab, ab, TRUE);
    fprintf(debug, "will be adding %d atoms\n",nadd);
  }
  
  /* Copy old atoms, making space for new ones */
  snew(newpdba,1);
  init_t_atoms(newpdba,natom+nadd,FALSE);
  newpdba->nres    = pdba->nres;   
  sfree(newpdba->resname);
  newpdba->resname = pdba->resname;
  snew(xn,natom+nadd);
  newi=0;
  for(i=0; (i<natom); i++) {
    /* check if this atom wasn't scheduled for deletion */
    if ( nab[i]==0 || ab[i][0].nname!=NULL ) {
      copy_atom(pdba,newpdba,*xptr,xn,i,newi);
      if (debug) fprintf(debug,"(%3d) %3d%4s %4s%3d",
			 i+1,newi+1,*newpdba->atomname[newi],
			 *newpdba->resname[newpdba->atom[newi].resnr],
			 newpdba->atom[newi].resnr+1);
      /* process the hacks for this atom */
      for(j=0; j<nab[i]; j++) {
	if ( ab[i][j].oname==NULL ) { /* add */
	  newi++;
	  newpdba->atom[newi].resnr=pdba->atom[i].resnr;
	  if (debug) fprintf(debug," + %d",newi+1);
	}
	if ( ab[i][j].nname!=NULL ) { /* add or replace */
	  snew(newpdba->atomname[newi],1);
	  *newpdba->atomname[newi]=strdup(ab[i][j].nname);
	  if ( ab[i][j].oname!=NULL && ab[i][j].atom ) { /* replace */
/* 	    newpdba->atom[newi].m    = ab[i][j].atom->m; */
/* 	    newpdba->atom[newi].q    = ab[i][j].atom->q; */
/* 	    newpdba->atom[newi].type = ab[i][j].atom->type; */
	  }
	  bSet=TRUE;
	  for(d=0; d<DIM; d++)
	    bSet = bSet && ab[i][j].newx[d]!=NOTSET;
	  if (bSet)
	    copy_rvec(ab[i][j].newx, xn[newi]);
	  if (debug) 
	    fprintf(debug," %s %g %g",*newpdba->atomname[newi],
		    newpdba->atom[newi].m,newpdba->atom[newi].q);
	}
      }
      if (debug) fprintf(debug,"\n");
      newi++;
    }
  }
  
  /* Clean up */
  for(i=0; i<natom; i++)
    free_t_hack(nab[i], &ab[i]);
  sfree(nab);
  sfree(ab);
  
  for(i=0; i < pdba->nr; i++) {
    sfree(*(pdba->atomname[i]));
    sfree(pdba->atomname[i]);
  }
  sfree(pdba->atomname);
  sfree(pdba->atom);
  sfree(pdba->pdbinfo);
  done_block(&(pdba->excl));
  sfree(pdba);
  
  sfree(*xptr);
  *xptr=xn;
  *pdbaptr=newpdba;
  
  return natom+nadd;
}

void deprotonate(t_atoms *atoms,rvec *x)
{
  int  i,j;
  
  for(i=j=0; (i<atoms->nr); i++) {
    if (((*atoms->atomname[i])[0] != 'H')) {
      atoms->atomname[j]=atoms->atomname[i];
      atoms->atom[j]=atoms->atom[i];
      copy_rvec(x[i],x[j]);
      j++;
    }
  }
  atoms->nr=j;
}

void protonate(t_atoms **atoms,rvec **x)
{
  static t_hackblock *ah=NULL;
  static int      nah;
  static t_symtab tab;
  
  if (ah == NULL) {
    nah=read_h_db("ffgmx2",&ah);
    open_symtab(&tab); 
  }
  deprotonate(*atoms,*x);
  
  add_h(atoms,x,nah,ah,NULL,NULL,0,0);
}
