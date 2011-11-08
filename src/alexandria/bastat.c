/*
 * $Id: bastat.c,v 1.8 2009/04/12 21:24:26 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
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
#include <stdlib.h>
#include <string.h>
#include "typedefs.h"
#include "physics.h"
#include "pdbio.h"
#include "smalloc.h"
#include "atomprop.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "xvgr.h"
#include "futil.h"
#include "gmx_statistics.h"

typedef struct {
  char *a1,*a2;
  int  nhisto;
  int  *histo;
} t_bond;

typedef struct {
  char *a1,*a2,*a3;
  int  nhisto;
  int  *histo;
} t_angle;

typedef struct {
  int nbond;
  t_bond *bond;
  int nangle;
  t_angle *angle;
} t_bonds;

int bondcmp(const void *a,const void *b)
{
  t_bond *ba = (t_bond *)a;
  t_bond *bb = (t_bond *)b;
  int d;
  
  if ((d = strcmp(ba->a1,bb->a1)) == 0)
    return strcmp(ba->a2,bb->a2);
  else
    return d;
}

int anglecmp(const void *a,const void *b)
{
  t_angle *ba = (t_angle *)a;
  t_angle *bb = (t_angle *)b;
  int d;
  
  if ((d = strcmp(ba->a1,bb->a1)) == 0)
    if ((d = strcmp(ba->a2,bb->a2)) == 0)
      return strcmp(ba->a3,bb->a3);
  
  return d;
}

void sort_bonds(t_bonds *b)
{
  qsort(b->bond,b->nbond,sizeof(b->bond[0]),bondcmp);
  qsort(b->angle,b->nangle,sizeof(b->angle[0]),anglecmp);
}

void add_bond(t_bonds *b,char *a1,char *a2,double blen,double spacing)
{
  int i,j,index;
  
  index = gmx_nint(blen/spacing);
  for(i=0; (i<b->nbond); i++) {
    if ((((strcmp(a1,b->bond[i].a1) == 0) && (strcmp(a2,b->bond[i].a2)) == 0)) ||
	(((strcmp(a2,b->bond[i].a1) == 0) && (strcmp(a1,b->bond[i].a2)) == 0))) {
      break;
    }
  }
  if (i == b->nbond) {
    b->nbond++;
    srenew(b->bond,b->nbond);
    if (strcmp(a1,a2) < 0) 
      {
	b->bond[i].a1     = strdup(a1);
	b->bond[i].a2     = strdup(a2);
      }
    else 
      {
	b->bond[i].a1     = strdup(a2);
	b->bond[i].a2     = strdup(a1);
      }
    b->bond[i].nhisto = 2*index+1;
    snew(b->bond[i].histo,b->bond[i].nhisto);
  }
  if (index >= b->bond[i].nhisto) {
    srenew(b->bond[i].histo,index+100);
    for(j=b->bond[i].nhisto; (j<index+100); j++)
      b->bond[i].histo[j] = 0;
    b->bond[i].nhisto = index+100;
  }
  b->bond[i].histo[index]++;
}

void add_angle(t_bonds *b,char *a1,char *a2,char *a3,double angle,double spacing)
{
  int i,j,index;
  
  index = gmx_nint(angle/spacing);
  for(i=0; (i<b->nangle); i++) {
    if ((strcmp(a2,b->angle[i].a2) == 0) && 
	(((strcmp(a1,b->angle[i].a1) == 0) && (strcmp(a3,b->angle[i].a3) == 0)) ||
	 ((strcmp(a3,b->angle[i].a1) == 0) && (strcmp(a1,b->angle[i].a3) == 0)))) {
      break;
    }
  }
  if (i == b->nangle) {
    b->nangle++;
    srenew(b->angle,b->nangle);
    b->angle[i].a2     = strdup(a2);
    if (strcmp(a1,a3) < 0) 
      {
	b->angle[i].a1     = strdup(a1);
	b->angle[i].a3     = strdup(a3);
      }
    else 
      {
	b->angle[i].a1     = strdup(a3);
	b->angle[i].a3     = strdup(a1);
      }
    b->angle[i].nhisto = (int) (180/spacing) + 1;
    snew(b->angle[i].histo,b->angle[i].nhisto);
  }

  b->angle[i].histo[index]++;
}

void lo_dump_histo(FILE *fp,int n,int histo[],double spacing)
{
  int j,j0,j1;
  double sum;
  
  for(j0=0; (j0<n) && (histo[j0] == 0); j0++)
    ;
  j0 = max(j0-1,0);
  for(j1=n-1; (j1>0) && (histo[j1] == 0); j1--)
    ;
  j1 = max(j1+1,n-1);
  sum=0;
  for(j=j0; (j<=j1); j++) 
    sum+=histo[j];
  if (sum == 0) {
    fprintf(stderr,"Empty histogram\n");
    exit(1);
  }
  for(j=j0; (j<=j1); j++) 
    fprintf(fp,"%g  %g\n",spacing*j,histo[j]/sum);
}

void dump_histo(t_bonds *b,double bspacing,double aspacing,output_env_t oenv)
{
  int  i,j,j0,j1;
  FILE *fp;
  char buf[256];
  
  for(i=0; (i<b->nbond); i++) {
    sprintf(buf,"bond-%s-%s.xvg",b->bond[i].a1,b->bond[i].a2);
    fp = xvgropen(buf,buf,"Distance (nm)","P (a.u.)",oenv);
    lo_dump_histo(fp,b->bond[i].nhisto,b->bond[i].histo,bspacing);
    fclose(fp);
  }
  for(i=0; (i<b->nangle); i++) {
    sprintf(buf,"angle-%s-%s-%s.xvg",
	    b->angle[i].a1,b->angle[i].a2,b->angle[i].a3);
    fp = xvgropen(buf,buf,"Angle (deg.)","P (a.u.)",oenv);
    lo_dump_histo(fp,b->angle[i].nhisto,b->angle[i].histo,aspacing);
    fclose(fp);
  }
}

void dump_table(t_bonds *b,double bspacing,double aspacing,
		gmx_atomprop_t aps)
{
  gmx_stats_t lsq1,lsq2;
  int i,j,k;
  real x,y,av,sig;
  FILE *fp;
  
  fp = fopen("smbonds.xml","w");
  fprintf(fp,"  <smbonds length_unit=\"nm\">\n");
  for(i=0; (i<b->nbond); i++) {
    lsq1 = gmx_stats_init();
    lsq2 = gmx_stats_init();
    for(j=0; (j<b->bond[i].nhisto); j++) {
      x = j*bspacing;
      y = b->bond[i].histo[j];
      gmx_stats_add_point(lsq1,y,x,0,0);
      gmx_stats_add_point(lsq2,x,y,0,0);
    }
    gmx_stats_get_average(lsq1,&av);
    av = 0.001*gmx_nint(1000*av);
    fprintf(fp,"    <smbond atom1=\"%s\" atom2=\"%s\" length=\"%g\" bondorder=\"1\"/>\n",
	    b->bond[i].a1,b->bond[i].a2,av);
    gmx_stats_done(lsq1);
    sfree(lsq1);
    gmx_stats_done(lsq2);
    sfree(lsq2);
  }
  fprintf(fp,"  </smbonds>\n");
  fprintf(fp,"  <smangles angle_unit=\"degree\">\n");
  for(i=0; (i<b->nangle); i++) {
    lsq1 = gmx_stats_init();
    for(j=0; (j<b->angle[i].nhisto); j++) {
      x = j*aspacing;
      y = b->angle[i].histo[j];
      gmx_stats_add_point(lsq1,y,x,0,0);
    }
    gmx_stats_get_average(lsq1,&av);
    av = 0.001*gmx_nint(1000*av);
    fprintf(fp,"    <smangle atom1=\"%s\" atom2=\"%s\" atom3=\"%s\" angle=\"%g\"/>\n",
	    b->angle[i].a1,b->angle[i].a2,b->angle[i].a3,av);
    gmx_stats_done(lsq1);
    sfree(lsq1);
  }
  fprintf(fp,"  </smangles>\n");
  fclose(fp);
}

int main(int argc,char *argv[])
{
  FILE *in;
  t_bonds *b;
  int i,j,k,l,model_nr,natoms,ePBC;
  t_atoms atoms;
  char title[2560];
  rvec *x,dx,dx2;
  matrix box;
  gmx_conect conect;
  gmx_atomprop_t aps;
  double ang;
  double bspacing = 0.002; /* nm */
  double aspacing = 0.1; /* degree */
  output_env_t oenv = NULL;
  
  aps = gmx_atomprop_init();
  snew(b,1);
  for(i=1; (i<argc); i++) {
    in = fopen(argv[i],"r");
    if (in) {
      get_pdb_coordnum(in,&natoms);
      frewind(in);
      if (natoms > 0) {
	init_t_atoms(&atoms,natoms,TRUE);
	snew(x,natoms);
	conect = gmx_conect_init();
	read_pdbfile(in,title,&model_nr,&atoms,x,&ePBC,box,FALSE,conect);
	fclose(in);
	get_pdb_atomnumber(&atoms,aps);
	for(j=0; (j<natoms); j++)
	  for(k=j+1; (k<natoms); k++)
	    if (gmx_conect_exist(conect,j,k)) {
	      rvec_sub(x[j],x[k],dx);
	      add_bond(b,*atoms.atomname[j],
		       *atoms.atomname[k],norm(dx),bspacing);
	      for(l=0; (l<natoms); l++) 
		if ((l != j) && gmx_conect_exist(conect,k,l)) {
		  rvec_sub(x[l],x[k],dx2);
		  ang = RAD2DEG*acos(cos_angle(dx,dx2));
		  add_angle(b,*atoms.atomname[j],*atoms.atomname[k],
			    *atoms.atomname[l],ang,aspacing);
		}
	    }
	sfree(x);
	free_t_atoms(&atoms,TRUE);
      }
      else
	fclose(in);
    }
  }
  sort_bonds(b);
  dump_histo(b,bspacing,aspacing,oenv);
  dump_table(b,bspacing,aspacing,aps);
  
  return 0;
}
