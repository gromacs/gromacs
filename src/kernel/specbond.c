/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_specbond_c = "$Id$";

#include "typedefs.h"
#include "pdbio.h"
#include "strdb.h"
#include "string2.h"
#include "smalloc.h"
#include "specbond.h"
#include "pdb2top.h"

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
  
  t_specbond *sb;
  char   r1buf[32],r2buf[32],a1buf[32],a2buf[32],nr1buf[32],nr2buf[32];
  double length;
  int    nb1,nb2;
  char   **lines;
  int    nlines,i,n;
  
  nlines = get_lines(sbfile,&lines);
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

static bool is_bond(int nsb,t_specbond sb[],char *res1,char *at1,
		    char *res2,char *at2,real d,int *nb,bool *bSwap)
{
  int i;
  
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

static void rename_1res(int natom,t_pdbatom pdba[],int resnr,char *nres)
{
  int i;
  
  for(i=0; (i<natom); i++)
    if (pdba[i].resnr == resnr)
      break;
  for( ; (i<natom) && (pdba[i].resnr == resnr); i++)
    strcpy(pdba[i].resnm,nres);
}

int mk_specbonds(int natoms,t_pdbatom pdba[],bool bInteractive,
		 t_ssbond **specbonds)
{
  t_specbond *sb;
  t_ssbond   *bonds;
  int  nsb;
  int  ncys,nbonds;
  int  *cysp,*sgp;
  int  *nBonded;
  bool bDoit,bSwap;
  int  i,j,nres;
  int  ai,aj,index_sb;
  real **d;
  
  sb=get_specbonds(&nsb);
  
  nres=pdba[natoms-1].resnr;
  snew(cysp,nres);
  snew(sgp,nres);
  snew(nBonded,nres);

  ncys = 0;
  for(i=0;(i<natoms);i++) {
    if (is_special(nsb,sb,pdba[i].resnm,pdba[i].atomnm)) {
      cysp[ncys]=pdba[i].resnr;
      sgp[ncys] =i;
      ncys++;
    }
  }
  
  /* distance matrix d[ncys][ncys] */
  snew(d,ncys);
  for(i=0; (i<ncys); i++)
    snew(d[i],ncys);

  for(i=0; (i<ncys); i++) 
    for(j=0; (j<ncys); j++) {
      ai=sgp[i];
      aj=sgp[j];
      d[i][j]=distance(pdba[ai].x,pdba[aj].x);
    }
  if (ncys > 1) {
    fprintf(stderr,"Special Atom Distance matrix\n%8s","");
    for(i=1; (i<ncys); i++) 
      fprintf(stderr," %3s%4d",pdba[sgp[i]].resnm,cysp[i]+1);
    fprintf(stderr,"\n");
    for(i=0; (i<ncys-1); i++) {
      fprintf(stderr," %3s%4d",pdba[sgp[i]].resnm,cysp[i]+1);
      for(j=1; (j<i+1); j++)
	fprintf(stderr,"%8s","");
      for( ; (j<ncys); j++)
	fprintf(stderr,"%6.3f  ",d[i][j]);
      fprintf(stderr,"\n");
    }
  }
  i=j=0;
  
  nbonds = 0;
  snew(bonds,ncys/2);
  
  for(i=0; (i<ncys); i++) {
    ai = sgp[i];
    for(j=i+1; (j<ncys); j++) {
      aj = sgp[j];
      
      if (is_bond(nsb,sb,pdba[ai].resnm,pdba[ai].atomnm,
		  pdba[aj].resnm,pdba[aj].atomnm,d[i][j],
		  &index_sb,&bSwap)) {
	if (bInteractive) {
	  fprintf(stderr,"Link Res%4d and Res%4d (y/n) ?",cysp[i]+1,cysp[j]+1);
	  bDoit=yesno();
	}
	else {
	  fprintf(stderr,"Linking Res%4d and Res%4d...\n",cysp[i]+1,cysp[j]+1);
	  bDoit=TRUE;
	}
	if (bDoit) {
	  /* Store the residue numbers in the bonds array */
	  bonds[nbonds].res1 = cysp[i];
	  bonds[nbonds].res2 = cysp[j];
	  bonds[nbonds].a1   = strdup(pdba[ai].atomnm);
	  bonds[nbonds].a2   = strdup(pdba[aj].atomnm);
	  if (bSwap) {
	    rename_1res(natoms,pdba,cysp[i],sb[index_sb].nres2);
	    rename_1res(natoms,pdba,cysp[j],sb[index_sb].nres1);
	  }
	  else {
	    rename_1res(natoms,pdba,cysp[i],sb[index_sb].nres1);
	    rename_1res(natoms,pdba,cysp[j],sb[index_sb].nres2);
	  }
	  nBonded[i]++;
	  nBonded[j]++;
	  
	  nbonds++;
	}
      }
    }
  }
  
  for(i=0; (i<ncys); i++)
    sfree(d[i]);
  sfree(d);
  sfree(nBonded);
  sfree(sgp);
  sfree(cysp);
 
  done_specbonds(nsb,sb);
  sfree(sb);
  
  *specbonds=bonds;
  
  return nbonds;
}

