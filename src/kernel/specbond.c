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
static char *SRCID_specbond_c = "$Id$";

#include <ctype.h>
#include <math.h>
#include "typedefs.h"
#include "pdbio.h"
#include "strdb.h"
#include "string2.h"
#include "smalloc.h"
#include "specbond.h"
#include "pdb2top.h"
#include "vec.h"

bool yesno(void)
{
  char c;

  do {
    c=toupper(fgetc(stdin));
  } while ((c != 'Y') && (c != 'N'));
  
  return (c == 'Y');
}

typedef struct {
  char *res1, *res2;
  char *atom1,*atom2;
  char *nres1,*nres2;
  int  nbond1,nbond2;
  real length;
} t_specbond;

static t_specbond *get_specbonds(int *nspecbond)
{
  static char  *sbfile="specbond.dat";
  
  t_specbond *sb=NULL;
  char   r1buf[32],r2buf[32],a1buf[32],a2buf[32],nr1buf[32],nr2buf[32];
  double length;
  int    nb1,nb2;
  char   **lines;
  int    nlines,i,n;
  
  nlines = get_lines(sbfile,&lines);
  if (nlines > 0)
    snew(sb,nlines);
  
  n = 0;
  for(i=0; (i<nlines); i++) {
    if (sscanf(lines[i],"%s%s%d%s%s%d%lf%s%s",
	       r1buf,a1buf,&nb1,r2buf,a2buf,&nb2,&length,nr1buf,nr2buf) != 9) 
      fprintf(stderr,"Invalid line '%s' in %s\n",lines[i],sbfile);
    else {
      sb[n].res1   = strdup(r1buf);
      sb[n].res2   = strdup(r2buf);
      sb[n].nres1  = strdup(nr1buf);
      sb[n].nres2  = strdup(nr2buf);
      sb[n].atom1  = strdup(a1buf);
      sb[n].atom2  = strdup(a2buf);
      sb[n].nbond1 = nb1;
      sb[n].nbond2 = nb2;
      sb[n].length = length;
      n++;
    }
    sfree(lines[i]);
  }
  if (nlines > 0)
    sfree(lines);
  fprintf(stderr,"%d out of %d lines of %s converted succesfully\n",
	  n,nlines,sbfile);
	  
  *nspecbond = n;
  
  return sb;
}

static void done_specbonds(int nsb,t_specbond sb[])
{
  int i;
  
  for(i=0; (i<nsb); i++) {
    sfree(sb[i].res1);
    sfree(sb[i].res2);
    sfree(sb[i].atom1);
    sfree(sb[i].atom2);
    sfree(sb[i].nres1);
    sfree(sb[i].nres2);
  }
}

static bool is_special(int nsb,t_specbond sb[],char *res,char *atom)
{
  int i;
  
  for(i=0; (i<nsb); i++) {
    if (((strcasecmp(sb[i].res1,res) == 0) && 
	 (strcasecmp(sb[i].atom1,atom) == 0)) ||
	((strcasecmp(sb[i].res2,res) == 0) && 
	 (strcasecmp(sb[i].atom2,atom) == 0)))
      return TRUE;
  }
  return FALSE;
}

static bool is_bond(int nsb,t_specbond sb[],t_atoms *pdba,int a1,int a2,
		    real d,int *nb,bool *bSwap)
{
  int i;
  char *at1,*at2,*res1,*res2;
  
  at1=*pdba->atomname[a1];
  at2=*pdba->atomname[a2];
  res1=*pdba->resname[pdba->atom[a1].resnr];
  res2=*pdba->resname[pdba->atom[a2].resnr];
  
  for(i=0; (i<nsb); i++) {
    *nb = i;
    if (((strcasecmp(sb[i].res1,res1) == 0)  && 
	 (strcasecmp(sb[i].atom1,at1) == 0) &&
	 (strcasecmp(sb[i].res2,res2) == 0)  && 
	 (strcasecmp(sb[i].atom2,at2) == 0))) {
      *bSwap = FALSE;
      if ((0.9*sb[i].length < d) && (1.1*sb[i].length > d))
	return TRUE;
    }
    if (((strcasecmp(sb[i].res1,res2) == 0)  && 
	 (strcasecmp(sb[i].atom1,at2) == 0) &&
	 (strcasecmp(sb[i].res2,res1) == 0)  && 
	 (strcasecmp(sb[i].atom2,at1) == 0))) {
      *bSwap = TRUE;
      if ((0.9*sb[i].length < d) && (1.1*sb[i].length > d))
	return TRUE;
    }
  }
  return FALSE;
}

static void rename_1res(t_atoms *pdba,int resnr,char *nres)
{
  sfree(*pdba->resname[resnr]);
  *pdba->resname[resnr]=strdup(nres);
}

int mk_specbonds(t_atoms *pdba,rvec x[],bool bInteractive,
		 t_ssbond **specbonds)
{
  t_specbond *sb=NULL;
  t_ssbond   *bonds=NULL;
  int  nsb,natoms;
  int  nspec,nbonds;
  int  *specp,*sgp;
  bool bDoit,bSwap;
  int  i,j,b,e,e2;
  int  ai,aj,index_sb;
  real **d;
  char buf[10];

  natoms=pdba->nr;
  /* Initiate variables that will be exported */  
  nbonds = 0;
  sb=get_specbonds(&nsb);
  
  if (nsb > 0) {
    snew(specp,natoms);
    snew(sgp,natoms);
    
    nspec = 0;
    for(i=0;(i<natoms);i++) {
      if (is_special(nsb,sb,*pdba->resname[pdba->atom[i].resnr],
		     *pdba->atomname[i])) {
	specp[nspec] = pdba->atom[i].resnr;
	sgp[nspec] = i;
	nspec++;
      }
    }
    /* distance matrix d[nspec][nspec] */
    snew(d,nspec);
    for(i=0; (i<nspec); i++)
      snew(d[i],nspec);
    
    for(i=0; (i<nspec); i++) {
      ai=sgp[i];
      for(j=0; (j<nspec); j++) {
	aj=sgp[j];
	d[i][j]=sqrt(distance2(x[ai],x[aj]));
      }
    }
    if (nspec > 1) {
#define MAXCOL 8
      fprintf(stderr,"Special Atom Distance matrix:\n");
      for(b=0; (b<nspec); b+=MAXCOL) {
	fprintf(stderr,"%8s","");
	e=min(b+MAXCOL, nspec-1);
	for(i=b; (i<e); i++) {
	  sprintf(buf,"%s%d",*pdba->resname[pdba->atom[sgp[i]].resnr],
		  specp[i]+1);
	  fprintf(stderr,"%8s",buf);
	}
	fprintf(stderr,"\n");
	e=min(b+MAXCOL, nspec);
	for(i=b+1; (i<nspec); i++) {
	  sprintf(buf,"%s%d",*pdba->resname[pdba->atom[sgp[i]].resnr],
		  specp[i]+1);
	  fprintf(stderr,"%8s",buf);
	  e2=min(i,e);
	  for(j=b; (j<e2); j++)
	  fprintf(stderr," %7.3f",d[i][j]);
	  fprintf(stderr,"\n");
	}
      }
    }
    i=j=0;
    
    snew(bonds,nspec);
  
    for(i=0; (i<nspec); i++) {
      ai = sgp[i];
      for(j=i+1; (j<nspec); j++) {
	aj = sgp[j];
	if (is_bond(nsb,sb,pdba,ai,aj,d[i][j],
		    &index_sb,&bSwap)) {
	  if (bInteractive) {
	    fprintf(stderr,"Link Res%4d and Res%4d (y/n) ?",specp[i]+1,specp[j]+1);
	    bDoit=yesno();
	  }
	  else {
	    fprintf(stderr,"Linking Res%4d and Res%4d...\n",specp[i]+1,specp[j]+1);
	    bDoit=TRUE;
	  }
	  if (bDoit) {
	    /* Store the residue numbers in the bonds array */
	    bonds[nbonds].res1 = specp[i];
	    bonds[nbonds].res2 = specp[j];
	    bonds[nbonds].a1   = strdup(*pdba->atomname[ai]);
	    bonds[nbonds].a2   = strdup(*pdba->atomname[aj]);
	    if (bSwap) {
	      rename_1res(pdba,specp[i],sb[index_sb].nres2);
	      rename_1res(pdba,specp[j],sb[index_sb].nres1);
	    }
	    else {
	      rename_1res(pdba,specp[i],sb[index_sb].nres1);
	      rename_1res(pdba,specp[j],sb[index_sb].nres2);
	    }
	    nbonds++;
	  }
	}
      }
    }
    
    for(i=0; (i<nspec); i++)
      sfree(d[i]);
    sfree(d);
    sfree(sgp);
    sfree(specp);
 
    done_specbonds(nsb,sb);
    sfree(sb);
  }
  
  *specbonds=bonds;
  
  return nbonds;
}

