/*
 * $Id$
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "strdb.h"
#include "futil.h"
#include "pdbio.h"
#include "string.h"
#include "vec.h"
#include "names.h"
#include "disco.h"
#include "gstat.h"

char *edc_names[edcNR+1] = { "NONBND", "BOND", "DISRE", NULL };

bool newres(int i,t_atom atom[])
{
  return ((i == 0) || (atom[i].resnr != atom[i-1].resnr));
}

void pr_corr(FILE *log,t_correct *c)
{
  fprintf(log,"Parameters for this disco run\n");
  fprintf(log,"maxnit        = %d\n",c->maxnit);
  fprintf(log,"nbcheck       = %d\n",c->nbcheck);
  fprintf(log,"nstprint      = %d\n",c->nstprint);
  fprintf(log,"nstranlist    = %d\n",c->nstranlist);
  fprintf(log,"ngrow         = %d\n",c->ngrow);
  fprintf(log,"bExplicit     = %s\n",BOOL(c->bExplicit));
  fprintf(log,"bChiral       = %s\n",BOOL(c->bChiral));
  fprintf(log,"bPep          = %s\n",BOOL(c->bPep));
  fprintf(log,"bDump         = %s\n",BOOL(c->bDump));
  fprintf(log,"bLowerOnly    = %s\n",BOOL(c->bLowerOnly));
  fprintf(log,"bRanlistFirst = %s\n",BOOL(c->bRanlistFirst));
  fprintf(log,"bCubic        = %s\n",BOOL(c->bCubic));
  fprintf(log,"bBox          = %s\n",BOOL(c->bBox));
  fprintf(log,"bCenter       = %s\n",BOOL(c->bCenter));
  fprintf(log,"ndist         = %d\n",c->ndist);
  fprintf(log,"npep          = %d\n",c->npep);
  fprintf(log,"nimp          = %d\n",c->nimp);
  fprintf(log,"lowdev        = %g\n",c->lodev);
  fflush(log);
}

t_correct *init_corr(int maxnit,int nstprint,int nbcheck,int nstranlist,
		     int ngrow,bool bExplicit,bool bChiral,bool bPep,
		     bool bDump,real lowdev,bool bLowerOnly,
		     bool bRanlistFirst,bool bCubic,bool bBox,bool bCenter)
{
  t_correct *c;
  
  snew(c,1);
  c->maxnit        = maxnit;
  c->nstprint      = nstprint;
  c->nbcheck       = nbcheck;
  c->nstranlist    = nstranlist;
  c->bExplicit     = bExplicit;
  c->bChiral       = bChiral;
  c->bPep          = bPep;
  c->bDump         = bDump;
  c->lodev         = lowdev;
  c->maxdist       = 0;
  c->ndist         = 0;
  c->ngrow         = ngrow;
  c->bLowerOnly    = bLowerOnly;
  c->bRanlistFirst = bRanlistFirst;
  c->bCubic        = bCubic;
  c->bBox          = bBox;
  c->bCenter       = bCenter;

  if (bRanlistFirst) {
    if (nstranlist != 0) {
      fprintf(stderr,
	      "Warning: setting nstranlist to 0 when using ranlistfirst");
      c->nstranlist = 0;
    }
  }
  return c;
}

void init_corr2(t_correct *c,int natom)
{
  int i,n,ai;
  
  if (c->ndist == 0)
    gmx_fatal(FARGS,"Can not make tags without distances. %s, %d",
		__FILE__,__LINE__);
		
  /* Initiate the ip index */
  snew(c->ip,c->ndist);

  /* Init boolean array */
  snew(c->bViol,natom);
  
  /* Initiate omega array */
  snew(c->omega,c->npep);

  /* Initiate improper array */
  srenew(c->idih,c->nimp+1);
  
  /* Now make the tags */
  snew(c->tag,natom+1);
  ai = c->d[0].ai;
  n  = 0;
  for(i=0; (i<c->ndist); i++) {
    c->ip[i] = i;
    if (c->d[i].ai != ai) {
      /* New tag starts here, but skip over possible missing atoms */
      while ((n < natom) && (n<c->d[i].ai)) {
	n++;
	c->tag[n] = i;
	ai = c->d[i].ai;
      }
      /*     else
	     gmx_fatal(FARGS,"Too many atoms, or distances not sorted");*/
    }
  }
  if (n < natom) {
    while (n < natom) {
      n++;
      c->tag[n] = i;
    }
  }
  else
    gmx_fatal(FARGS,"Too many atoms, or distances not sorted");
  if (debug)
    fprintf(debug,"There are %d tags for %d atoms\n",n,natom);
  
  if (debug) 
    for(i=0; (i<=natom); i++)
      fprintf(debug,"tag[%5d] = %d\n",i,c->tag[i]);
}

void center_in_box(int natom,rvec xin[],matrix box,rvec xout[])
{
  int  i,m;
  rvec xcm,dx;
  
  /* Dump the optimization trajectory to an xtc file  */
  /* Put the whole thing in the center of the box 1st */
  clear_rvec(xcm);
  for(i=0; (i<natom); i++) {
    copy_rvec(xin[i],xout[i]);
    rvec_inc(xcm,xin[i]);
  }
  for(m=0; (m<DIM); m++)
    dx[m] = 0.5*box[m][m]-xcm[m]/natom;
  for(i=0; (i<natom); i++) 
    rvec_inc(xout[i],dx);
}

void define_peptide_bonds(FILE *log,t_atoms *atoms,t_correct *c)
{
  int    i,npep,naa,nlist;
  char   **aa;
  t_dlist *dlist;
  
  naa   = get_strings("aminoacids.dat",&aa);
  dlist = mk_dlist(log,atoms,&nlist,TRUE,TRUE,FALSE,FALSE,0,1,naa,aa);
  for(i=0; (i<naa); i++)
    sfree(aa[i]);
  sfree(aa);

  npep  = 0;
  snew(c->pepbond,nlist);
  for(i=0; (i<nlist); i++) {
    if (has_dihedral(edOmega,&dlist[i])) {
      c->pepbond[npep].ai = dlist[i].atm.minC;
      c->pepbond[npep].aj = dlist[i].atm.minO;
      c->pepbond[npep].ak = dlist[i].atm.N;
      c->pepbond[npep].al = dlist[i].atm.H;
      c->pepbond[npep].am = dlist[i].atm.Cn[1];
      npep++;
    }
  }
  c->npep = npep;
  if (debug)
    pr_dlist(debug,nlist,dlist,1.0,0,TRUE,TRUE,FALSE,TRUE,0);
  sfree(dlist);
  fprintf(log,"There are %d peptide bonds\n",npep);
}

static char *aname(t_atoms *atoms,int i)
{
  static char buf[32];
  
  sprintf(buf,"%4s%3d:%4s", 
	  *(atoms->resname[atoms->atom[i].resnr]),
	  atoms->atom[i].resnr+1,
	  *(atoms->atomname[i]));
  
  return buf;
}

int find_atom(t_atoms *atoms,int resnr,int j0,char *nm)
{
  int j,aa=-1;
  
  for(j=j0; (j < atoms->nr) && (atoms->atom[j].resnr == resnr); j++)
    if (strcmp(*atoms->atomname[j],nm) == 0) {
      aa = j;
      break;
    }
  return aa;
}

void define_impropers(FILE *log,t_atoms *atoms,t_correct *c)
{
  typedef struct {
    char *res,*aa[4];
  } t_impdef;
  
  t_impdef id[] = {
    { NULL,   { "CA", "N",  "C",   "CB"  } },
    { NULL,   { "N",  "CA", "H1",  "H3"  } },
    { "LYSH", { "NZ", "CE", "HZ2", "HZ1" } },
    { "LYSH", { "NZ", "CE", "HZ1", "HZ3" } },
    { "LEU",  { "CG", "CD2", "CD1", "CB" } },
    { "VAL",  { "CB", "CG2", "CG1", "CA" } },
    { "ILE",  { "CB", "CG1", "CG2", "CA" } },
    { "THR",  { "CB", "OG1", "CG2", "CA" } }
  };
#define NID asize(id)
  int i,j,k,l,aa[4],nimp;
  
  nimp = 0;
  for(i=0; (i<atoms->nres); i++) {
    /* Skip until residue number */
    for(j=0; (j<atoms->nr) && (atoms->atom[j].resnr != i); j++) 
      ;
    for(k=0; (k<NID); k++) {
      if ((id[k].res == NULL) || 
	  (strcmp(id[k].res,*atoms->resname[i]) == 0)) {
	/* This (i) is the right residue to look for this improper (k) */
	for(l=0; (l<4); l++)
	  aa[l] = find_atom(atoms,i,j,id[k].aa[l]);
	if ((aa[0] != -1) && (aa[1] != -1) && (aa[2] != -1) && (aa[3] != -1)) {
	  srenew(c->imp,nimp+1);
	  c->imp[nimp].ai = aa[0];
	  c->imp[nimp].aj = aa[1];
	  c->imp[nimp].ak = aa[2];
	  c->imp[nimp].al = aa[3];
	  nimp++;
	}
      }
    }
  }
  c->nimp = nimp;

  fprintf(log,"There are %d impropers\n",c->nimp);
  if (debug) {
    fprintf(debug,"Overview of improper dihedrals\n");
    for(i=0; (i<c->nimp); i++) { 
      fprintf(debug,"  %s",aname(atoms,c->imp[i].ai));
      fprintf(debug,"  %s",aname(atoms,c->imp[i].aj));
      fprintf(debug,"  %s",aname(atoms,c->imp[i].ak));
      fprintf(debug,"  %s\n",aname(atoms,c->imp[i].al));
    }
  }
}

void pr_dist(FILE *fp,bool bHeader,t_correct *c,int i)
{
  real ideal=0;
  
  if (bHeader)
    fprintf(fp,"#%4s%5s%10s%10s%10s\n","ai","aj","ideal","lb","ub");
  switch (c->d[i].cons_type) {
  case edcBOND:
    ideal = 5*(c->d[i].lb+c->d[i].ub);
    break;
  case edcNONBOND:
    ideal = 0;
    break;
  case edcDISRE:
    ideal = -1;
    break;
  default:
    gmx_fatal(FARGS,"cons_type for distance %d = %d\n",i,c->d[i].cons_type);
  }
  fprintf(fp,"%5d%5d%10.5f%10.5f%10.5f\n",1+c->d[i].ai,1+c->d[i].aj,
	  ideal,10*c->d[i].lb,10*c->d[i].ub);
}

void pr_distances(FILE *fp,t_correct *c)
{
  int i;
  
  for(i=0; (i<c->ndist); i++) 
    pr_dist(fp,(i==0),c,i);
}

void add_dist(t_correct *c,int ai,int aj,real lb,real ideal,real ub,real w[])
{
  int n = c->ndist;
  
  if ((w[ai] != 0) || (w[aj] != 0)) {
    if (n == c->maxdist) {
      c->maxdist += 100;
      srenew(c->d,c->maxdist);
    }
    c->d[n].ai = ai;
    c->d[n].aj = aj;
    if (ideal > 0)
      c->d[n].cons_type = edcBOND;
    else if (ideal == 0.0)
      c->d[n].cons_type = edcNONBOND;
    else
      c->d[n].cons_type = edcDISRE;
    c->d[n].lb = lb;
    c->d[n].ub = ub;
    c->d[n].wi = w[ai]/(w[ai]+w[aj]);
    c->ndist++;
  }
}

void read_dist(FILE *log,char *fn,int natom,t_correct *c)
{
  FILE   *fp;
  char   buf[256];
  int    ai,aj,i,nline=0;
  double ideal,lb,ub;
  int    nnn[edcNR];
  
  for(i=0; (i<edcNR); i++)
    nnn[i] = 0;
    
  fp = ffopen(fn,"r");
  while (fgets(buf,255,fp)) {
    nline ++;
    if (buf[0] != '#') {
      if (sscanf(buf,"%d%d%lf%lf%lf",&ai,&aj,&ideal,&lb,&ub) != 5)
	fprintf(stderr,"Error in %s, line %d\n",fn,nline);
      else {
	add_dist(c,ai-1,aj-1,lb*0.1,ideal*0.1,ub*0.1,c->weight);
	nnn[c->d[c->ndist-1].cons_type]++;
      }
    }
  }
  ffclose(fp);
  
  fprintf(log,"Read distances from file %s\n",fn);
  for(i=0; (i<edcNR); i++)
    fprintf(log,"There were %d %s distances\n",
	    nnn[i],edc_names[i]);
}

real *read_weights(char *fn,int natom)
{
  t_atoms newatoms;
  int     i;
  rvec    *x;
  matrix  box;
  char    title[256];
  real    *w;

  /* Read the weights from the occupancy field in the pdb file */  
  snew(x,natom);
  init_t_atoms(&newatoms,natom,TRUE);
  read_pdb_conf(fn,title,&newatoms,x,box,FALSE);
  snew(w,newatoms.nr);
  for(i=0; (i<newatoms.nr); i++)
    w[i] = newatoms.pdbinfo[i].occup;
  free_t_atoms(&newatoms);
  sfree(x);
  
  return w;
}

void check_dist(FILE *log,t_correct *c)
{
  int  i;
  real tol=0.001;
  
  fprintf(log,"Checking distances for internal consistency\n");
  for(i=0; (i<c->ndist); i++) {
    if ((c->d[i].ub != 0) && ((c->d[i].ub - c->d[i].lb) < tol)) {
      pr_dist(log,TRUE,c,i);
    }
  }
}

