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

#define MAXHH 3

static int hack_residue(int resnr,int natom,t_pdbatom pdba[],t_add_block *ab[],
			t_terblock *tb,int fac,int nTH[],bool *bDel)
{
  int i,j,dadd=0;
 
  if (tb == NULL)
    return 0;
    
  for(i=0; (i<natom); i++) 
    if (pdba[i].resnr == resnr)
      break;
      
  for(   ; ((i<natom) && (pdba[i].resnr==resnr)); i++) {
    for(j=0; (j<tb->nadd); j++) {
      if (strcmp(pdba[i].atomnm,tb->ab[j].na[0]) == 0) {
	/* Check whether this atom has already been set */
	nTH[i]=fac*(j+1);
	if (ab[i] != NULL)
	  dadd-=ab[i]->nh;
	ab[i]=&(tb->ab[j]);
	dadd+=ab[i]->nh;
	break;
      }
    }
    for(j=0; (j<tb->ndel); j++)
      if (strcmp(pdba[i].atomnm,tb->nm_del[j]) == 0) {
	bDel[i]=TRUE;
	dadd--;
      }
  }
  return dadd;
}

static char *Hnum[] = { "1", "2", "3", "4", "5", "6" };

static void hack_atoms(t_terblock *tdb,int index,
		       int j0,int nadd,int natom,t_pdbatom pdba[])
{
  int j;

  if (tdb == NULL)
    return;
  for(j=j0; (j<j0+nadd); j++) {
    assert(j<natom);
    strcpy(pdba[j].atomnm,tdb->add_nm[index]);
    if (nadd > 1)
      strcat(pdba[j].atomnm,Hnum[j-j0]);
    pdba[j].m    = tdb->adder[index].m;
    pdba[j].q    = tdb->adder[index].q;
    pdba[j].type = tdb->adder[index].type;
  }
}

static void replace_atoms(t_terblock *tdb,int natom,t_pdbatom pdba[],
			  int resnr)
{
  int i,j;
 
  if (tdb == NULL)
    return;
  for(i=0; ((i<natom) && (pdba[i].resnr != resnr)); i++) 
    ;
  for(   ; ((i<natom) && (pdba[i].resnr == resnr)); i++) {
    for(j=0; (j<tdb->nreplace); j++) {
      if (strcmp(tdb->nm_repl[j],pdba[i].atomnm) == 0) {
	strcpy(pdba[i].atomnm,tdb->new_nm[j]);
	pdba[i].m    = tdb->repl_by[j].m;
	pdba[i].q    = tdb->repl_by[j].q;
	pdba[i].type = tdb->repl_by[j].type;
      }
    }
  }
}

int add_h(int natom,t_pdbatom **pdbaptr,int nah,t_addh ah[],rvec **xptr,
	  t_terblock *ntdb,t_terblock *ctdb,
	  int rN,int rC)
{
  /* Global variable from h_db.c */
  extern int ncontrol[];
  
  t_pdbatom   *newpdba,*pdba;
  t_add_block **ab;      /* Array of pointer to add blocks */
  int         *nTH;
  bool        *bDel;
  int         tadd=0;
  int         rnr=-83;
  int         i,j,j0,k,l,m,ncntl;
  int         na[4],nh[3];
  t_addh      *ahptr=NULL;
  rvec        *x;

  snew(ab,natom);
  snew(nTH,natom);
  snew(bDel,natom);
  pdba=*pdbaptr;
  for(i=0; (i<natom); i++) {
    if (rnr != pdba[i].resnr) {
      rnr=pdba[i].resnr;
      ahptr=search_h_db(nah,ah,pdba[i].resnm);
    }
    if (ahptr != NULL) {
      for(j=0; (j<ahptr->n_add); j++) {
	if (strcmp(pdba[i].atomnm,ahptr->ab[j].na[0]) == 0) {
	  /* Check whether this atom has already been set */
	  assert(ab[i] == NULL);
	  ab[i]=&(ahptr->ab[j]);
	  tadd+=ab[i]->nh;
	  break;
	}
      }
    }
  }
  
  /* Modify for termini HERE */
  tadd+=hack_residue(rN,natom,pdba,ab,ntdb,1,nTH,bDel);
  tadd+=hack_residue(rC,natom,pdba,ab,ctdb,-1,nTH,bDel);

  /* Copy old atoms, making space for new ones */
  snew(newpdba,natom+tadd);
  snew(x,natom+tadd);
  j=0;
  for(i=0; (i<natom); i++) {
    if (!bDel[i]) {
      memcpy(&(newpdba[j]),&(pdba[i]),sizeof(pdba[i]));
      j++;
      if (ab[i] != NULL) {
	j0=j;
	for(k=0; (k<ab[i]->nh); k++,j++) {
	  memcpy(&(newpdba[j]),&(pdba[i]),sizeof(pdba[i]));
	  newpdba[j].atomnm[0]='H';
	  if (ab[i]->nh > 1)
	    strcat(newpdba[j].atomnm,Hnum[k]);
	}
	if (nTH[i]!=0) {
	  if (nTH[i] > 0)
	    hack_atoms(ntdb,nTH[i]-1,j0,ab[i]->nh,natom+tadd,newpdba);
	  else if (nTH[i] < 0)
	    hack_atoms(ctdb,-nTH[i]-1,j0,ab[i]->nh,natom+tadd,newpdba);
	  /* check if atoms should be added to added atoms */
	  for(k=j0; (k<j0+ab[i]->nh); k++) {
	    for(l=0; (l<ntdb->nadd); l++) {
	      if (strcmp(newpdba[k].atomnm,ntdb->ab[l].na[0]) == 0) {
		tadd+=ntdb->ab[l].nh;
		srenew(newpdba,natom+tadd);
		srenew(x,natom+tadd);
		for(m=0; (m<ntdb->ab[l].nh); m++)
		  memcpy(&(newpdba[j+m]),&(pdba[i]),sizeof(pdba[i]));
		hack_atoms(ntdb,l,j,ntdb->ab[l].nh,natom+tadd,newpdba);
		j+=ab[i]->nh;
		break;
	      }
	    }
	    for(l=0; (l<ctdb->nadd); l++) {
	      if (strcmp(newpdba[k].atomnm,ctdb->ab[l].na[0]) == 0) {
		tadd+=ctdb->ab[l].nh;
		srenew(newpdba,natom+tadd);
		srenew(x,natom+tadd);
		for(m=0; (m<ctdb->ab[l].nh); m++)
		  memcpy(&(newpdba[j+m]),&(pdba[i]),sizeof(pdba[i]));
		hack_atoms(ctdb,l,j,ctdb->ab[l].nh,natom+tadd,newpdba);
		j+=ab[i]->nh;
		break;
	      }
	    }
	  }
	}
      }
    }
  }
  
  /* Now let me do this... */
  replace_atoms(ntdb,natom+tadd,newpdba,rN);
  replace_atoms(ctdb,natom+tadd,newpdba,rC);
  
  /* Copy the coordinates for calc position */
  for(i=0; (i<natom+tadd); i++)
    copy_rvec(newpdba[i].x,x[i]);
  
  /* Now calc the positions */
  j=0;
  for(i=0; (i<natom); i++,j++) {
    if (ab[i] != NULL) {
      rnr=newpdba[j].resnr;
      ncntl = ncontrol[ab[i]->tp];
      for(m=0; (m<ncntl); m++) {
	na[m]=pdbasearch_atom(ab[i]->na[m],rnr,natom+tadd,newpdba);
	if (na[m] == -1)
	  fatal_error(0,"Atom %s not found in residue %s%s"
		      " while adding hydrogens",
		      ab[i]->na[m],newpdba[j].resnm,newpdba[j].pdbresnr);
      }
      
      for(m=0; (m<ab[i]->nh); m++)
	nh[m]=na[0]+1+m;
      for(   ; (m<3); m++)
	nh[m]=-1;
      calc_h_pos(ab[i]->tp,nh,na,x);
      j+=ab[i]->nh;
    }
  }
  for(i=0; (i<natom+tadd); i++)
    copy_rvec(x[i],newpdba[i].x);
    
  /* Clean up */
  sfree(ab);
  sfree(nTH);
  sfree(bDel);
  sfree(*pdbaptr);
  *xptr=x;
  *pdbaptr=newpdba;
  
  return natom+tadd;
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
  static t_addh   *ah=NULL;
  static int      nah;
  static t_symtab tab;
  
  t_pdbatom     *pdba;
  t_atoms       *newatoms;
  rvec          *newx;
  int           natom;
  
  if (ah == NULL) {
    nah=read_h_db("ffgmx2",&ah);
    open_symtab(&tab); 
  }
  deprotonate(*atoms,*x);
  
  pdba=atoms2pdba(*atoms,*x);
  
  natom=add_h((*atoms)->nr,&pdba,nah,ah,x,NULL,NULL,0,0);

  snew(newatoms,1);
  pdb2atoms(natom,pdba,newatoms,&newx,&tab);
  
  sfree(*x);
  sfree((*atoms)->atomname);
  sfree((*atoms)->resname);
  sfree((*atoms)->atom);
  sfree(*atoms);
  sfree(pdba);
  
  *x=newx;
  *atoms=newatoms;
}

