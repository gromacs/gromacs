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

/* types of hack to be made */
enum {
  ehtNONE, ehtNTER, ehtCTER, ehtDUM
};

char *HackName[] = { "Nothing", "N-terminus", "C-terminus" };

#define MAXHH 3

static void copy_atom(t_atoms *atoms1,t_atoms *atoms2,rvec *x1,rvec *x2,
		      int a1,int a2)
{
  memcpy(&(atoms2->atom[a2]),&(atoms1->atom[a1]),(size_t)sizeof(t_atom));
  snew(atoms2->atomname[a2],1);
  *atoms2->atomname[a2]=strdup(*atoms1->atomname[a1]);
  copy_rvec(x1[a1],x2[a2]);
}

static int hack_residue(int resnr,t_atoms *pdba,t_add_block *ab[],
			t_hackblock *tb,
			int type,int HT[], int nTH[],bool *bDel)
{
  int i,j,dadd=0,natom;
 
  if (debug)
    if (tb)
      fprintf(debug,"Hacking %d atoms of %s into residue %d\n",
	      tb->nadd,HackName[type],resnr);
    else
      fprintf(debug,"Not hacking %s into residue %d\n",HackName[type],resnr);
  
  natom=pdba->nr;
  
  if (tb == NULL)
    return 0;
    
  for(i=0; (i<natom); i++) 
    if (pdba->atom[i].resnr == resnr)
      break;
      
  for(   ; ((i<natom) && (pdba->atom[i].resnr==resnr)); i++) {
    for(j=0; (j<tb->nadd); j++) {
      if (strcmp(*pdba->atomname[i],tb->ab[j].na[0]) == 0) {
	/* Check whether this atom has already been set */
	nTH[i]=j+1;
	HT[i]=type;
	if (ab[i] != NULL)
	  dadd-=ab[i]->nh;
	ab[i]=&(tb->ab[j]);
	dadd+=ab[i]->nh;
	break;
      }
    }
    for(j=0; (j<tb->ndel); j++)
      if (strcmp(*pdba->atomname[i],tb->nm_del[j]) == 0) {
	bDel[i]=TRUE;
	dadd--;
      }
  }
  
  return dadd;
}

static char *Hnum[] = { "1", "2", "3", "4", "5", "6" };

static void hack_atoms(t_hackblock *tdb,int index,
		       int j0,int nadd,t_atoms *pdba)
{
  int j;
  char buf[STRLEN];

  if (debug)
    if (tdb)
      printf("Replacing atoms %d-%d (addition %d)\n",
	     j0,j0+nadd-1,index);
    else
      printf("Not replacing atoms (%d-%d addition %d)\n",
	     j0,j0+nadd-1,index);
  
  if (tdb == NULL)
    return;
  
  for(j=j0; (j<j0+nadd); j++) {
    assert(j<pdba->nr);
    strcpy(buf,tdb->add_nm[index]);
    if (nadd > 1)
      strcat(buf,Hnum[j-j0]);
    sfree(*pdba->atomname[j]);
    *pdba->atomname[j]=strdup(buf);
    pdba->atom[j].m    = tdb->adder[index].m;
    pdba->atom[j].q    = tdb->adder[index].q;
    pdba->atom[j].type = tdb->adder[index].type;
  }
}

static void replace_atoms(t_hackblock *tdb,t_atoms *pdba,int resnr)
{
  int i,j,natom;
  
  natom=pdba->nr;
  if (tdb == NULL)
    return;
  for(i=0; ((i<natom) && (pdba->atom[i].resnr != resnr)); i++) 
    ;
  for(   ; ((i<natom) && (pdba->atom[i].resnr == resnr)); i++) {
    for(j=0; (j<tdb->nreplace); j++) {
      if (strcmp(tdb->nm_repl[j],*pdba->atomname[i]) == 0) {
	sfree(*pdba->atomname[i]);
	*pdba->atomname[i]=strdup(tdb->new_nm[j]);
	pdba->atom[i].m    = tdb->repl_by[j].m;
	pdba->atom[i].q    = tdb->repl_by[j].q;
	pdba->atom[i].type = tdb->repl_by[j].type;
      }
    }
  }
}

int pdbasearch_atom(char *name,int resnr,t_atoms *pdba)
{
  int  i,natom;
  
  natom=pdba->nr;
  if (name[0] == '-') {
    name++;
    resnr--;
  }
  for(i=0; (i<natom) && (pdba->atom[i].resnr != resnr); i++)
    ;
  for( ; (i<natom) && (pdba->atom[i].resnr == resnr); i++) {
    if (strcmp(name,*pdba->atomname[i]) == 0)
      return i;
  }
  return -1;
}

int add_h(t_atoms **pdbaptr,rvec **xptr,int nah,t_addh ah[],
	  t_hackblock *ntdb,t_hackblock *ctdb,
	  int rN,int rC)
{
  /* Global variable from h_db.c */
  extern int ncontrol[];
  
  t_atoms     *newpdba,*pdba;
  t_add_block **ab;      /* Array of pointer to add blocks */
  int         *nTH,*HT;
  bool        *bDel;
  int         tadd=0;
  int         rnr=NOTSET;
  int         i,j,j0,k,l,m,ncntl,natom;
  int         na[4],nh[3];
  t_addh      *ahptr=NULL;
  rvec        *xn;
  char        buf[STRLEN];

  /* set flags for adding hydrogens (accoring to hdb) */
  natom=(*pdbaptr)->nr;
  snew(ab,natom);
  snew(HT,natom);
  snew(nTH,natom);
  snew(bDel,natom);
  pdba=*pdbaptr;
  rnr=-83;
  for(i=0; (i<natom); i++) {
    if (rnr != pdba->atom[i].resnr) {
      rnr=pdba->atom[i].resnr;
      ahptr=search_h_db(nah,ah,*pdba->resname[pdba->atom[i].resnr]);
    }
    if (ahptr != NULL) {
      for(j=0; (j<ahptr->n_add); j++) {
	if (strcmp(*pdba->atomname[i],ahptr->ab[j].na[0]) == 0) {
	  /* Check whether this atom has already been set */
	  assert(ab[i] == NULL);
	  ab[i]=&(ahptr->ab[j]);
	  tadd+=ab[i]->nh;
	  break;
	}
      }
    }
  }
  
  /* set flags for adding atoms to termini (according to tdb) */
  if (rN>=0)
    tadd+=hack_residue(rN,pdba,ab,ntdb,ehtNTER,HT,nTH,bDel);
  if (rC>=0)
    tadd+=hack_residue(rC,pdba,ab,ctdb,ehtCTER,HT,nTH,bDel);
  
  /* Copy old atoms, making space for new ones */
  snew(newpdba,1);
  init_t_atoms(newpdba,natom+tadd,FALSE);
  newpdba->nres   =pdba->nres;   
  newpdba->resname=pdba->resname;
  snew(xn,natom+tadd);
  j=0;
  for(i=0; (i<natom); i++) {
    if (!bDel[i]) {
      copy_atom(pdba,newpdba,*xptr,xn,i,j);
      j++;
      if (ab[i] != NULL) {
	j0=j;
	for(k=0; (k<ab[i]->nh); k++,j++) {
	  strcpy(buf,*pdba->atomname[i]);
	  buf[0]='H';
	  if (ab[i]->nh > 1)
	    strcat(buf,Hnum[k]);
	  snew(newpdba->atomname[j],1);
	  *newpdba->atomname[j]=strdup(buf);
	  newpdba->atom[j].resnr=pdba->atom[i].resnr;
	}
	if (nTH[i]!=0) {
	  switch(HT[i]) {
	  case ehtNTER :
	    hack_atoms(ntdb,nTH[i]-1,j0,ab[i]->nh,newpdba);
	    break;
	  case ehtCTER :
	    hack_atoms(ctdb,nTH[i]-1,j0,ab[i]->nh,newpdba);
	    break;
	  default:
	    fatal_error(0,"Death horror error in genhydro (nTH=%d, HT=%d)\n",
			nTH[i],HT[i]);
	    break;
	  }
	  /* check if atoms should be added to added atoms */
	  for(k=j0; (k<j0+ab[i]->nh); k++) {
	    for(l=0; (l<ntdb->nadd); l++) {
	      if (strcmp(*newpdba->atomname[k],ntdb->ab[l].na[0]) == 0) {
		tadd+=ntdb->ab[l].nh;
		newpdba->nr=natom+tadd;
		srenew(newpdba->atom,natom+tadd);
		srenew(newpdba->atomname,natom+tadd);
		srenew(xn,natom+tadd);
		for(m=0; (m<ntdb->ab[l].nh); m++)
		  copy_atom(pdba,newpdba,*xptr,xn,i,j+m);
		hack_atoms(ntdb,l,j,ntdb->ab[l].nh,newpdba);
		j+=ab[i]->nh;
		break;
	      }
	    }
	    for(l=0; (l<ctdb->nadd); l++) {
	      if (strcmp(*newpdba->atomname[k],ctdb->ab[l].na[0]) == 0) {
		tadd+=ctdb->ab[l].nh;
		newpdba->nr=natom+tadd;
		srenew(newpdba->atom,natom+tadd);
		srenew(newpdba->atomname,natom+tadd);
		srenew(xn,natom+tadd);
		for(m=0; (m<ctdb->ab[l].nh); m++)
		  copy_atom(pdba,newpdba,*xptr,xn,i,j+m);
		hack_atoms(ctdb,l,j,ctdb->ab[l].nh,newpdba);
		j+=ab[i]->nh;
		break;
	      }
	    }
	  }
	}
      }
    }
  }
  
  /* atoms in termini should be replaced (new name/type/mass/charge) */
  if (rN>=0)
    replace_atoms(ntdb,newpdba,rN);
  if (rC>=0)
    replace_atoms(ctdb,newpdba,rC);
  
  /* Now calc the positions */
  j=0;
  for(i=0; (i<pdba->nr); i++,j++) {
    if (ab[i] != NULL) {
      rnr   = /*newpdba->atom[j].resnr;*/
	pdba->atom[i].resnr;
      ncntl = ncontrol[ab[i]->tp];
      for(m=0; (m<ncntl); m++) {
	na[m]=pdbasearch_atom(ab[i]->na[m],rnr,newpdba);
	if (na[m] == -1)
	  fatal_error(0,"Atom %s not found in residue %s%d"
		      " while adding hydrogens",
		      ab[i]->na[m],*newpdba->resname[rnr],rnr+1);
      }
      
      for(m=0; (m<ab[i]->nh); m++)
	nh[m]=na[0]+1+m;
      for(   ; (m<3); m++)
	nh[m]=-1;
      calc_h_pos(ab[i]->tp,nh,na,xn);
      j+=ab[i]->nh;
    }
  }

  /* Clean up */
  sfree(ab);
  sfree(HT);
  sfree(nTH);
  sfree(bDel);
  sfree(*pdbaptr);

  *xptr=xn;
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
  
  if (ah == NULL) {
    nah=read_h_db("ffgmx2",&ah);
    open_symtab(&tab); 
  }
  deprotonate(*atoms,*x);
  
  add_h(atoms,x,nah,ah,NULL,NULL,0,0);
}




