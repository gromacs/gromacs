/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_atomprop_c = "$Id$";
#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "fatal.h"
#include "atomprop.h"
#include "macros.h"
#include "index.h"
#include "strdb.h"
#include "copyrite.h"

typedef struct {
  char atomname[10];
  char resname[10];
  real prop;
} t_prop;

static void add_prop(real p, char *resnm, char *atomnm, 
		     int *nprop, t_prop **props, int *maxprop)
{
  if (*nprop >= *maxprop) {
    *maxprop+=10;
    srenew(*props,*maxprop);
  }
  strcpy((*props)[*nprop].atomname,atomnm);
  strcpy((*props)[*nprop].resname, resnm);
  (*props)[*nprop].prop=p;
  (*nprop)++;
}

static int read_props(char *propdata,t_prop **props,int *maxprop) 
{
  FILE   *fp;
  char   line[STRLEN],resname[10],atomname[10];
  double pp;
  int    i;
  
  fp=libopen(propdata);
  *props=NULL;
  *maxprop = 0;
  i=0;
  while(get_a_line(fp,line,STRLEN)) {
    if (sscanf(line,"%s %s %lf",resname,atomname,&pp)==3)
      add_prop(pp,resname,atomname,&i,props,maxprop);
    else
      if (debug)
	fprintf(stderr,"ERROR in file %s at line %d\n",propdata,i);
  }
  fclose(fp);
  
  return i;
}

static void write_props(FILE *fp,int nprop,t_prop props[]) 
{
  int i;
  
  for (i=0; (i<nprop); i++) 
    fprintf(fp,"%10s  %10s  %10g\n",
	    props[i].resname,props[i].atomname,props[i].prop);
}

/* NOTFOUND should be smallest, others larger in increasing priority */
#define NOTFOUND -3
#define WILDCARD NOTFOUND+1
#define PROTEIN  NOTFOUND+2

/* return number of matching characters, 
   or NOTFOUND if not at least all characters in char *database match */
static int dbcmp_len(char *search, char *database)
{
  int i;
  
  i=0;
  while(search[i] && database[i] && (search[i]==database[i]) )
    i++;
  
  if (database[i])
    i=NOTFOUND;
  return i;
}

static bool get_prop(real *p, char *resname, char *atomnm, 
		     int nprop, t_prop *props)
{
  int  i,j=0,alen,rlen,malen,mrlen;
  char *atomname;
  
  malen = NOTFOUND;
  mrlen = NOTFOUND;
  
  if (isdigit(atomnm[0])) {
    /* put digit after atomname */
    snew(atomname,strlen(atomnm)+1);
    for (i=1; (i<strlen(atomnm)); i++)
      atomname[i-1]=atomnm[i];
    atomname[i++]=atomnm[0];
    atomname[i]='\0';
  } else 
    atomname=atomnm;
      
  for(i=0; (i<nprop); i++) {
    if ( (strcmp(props[i].resname,"*")==0) ||
	 (strcmp(props[i].resname,"???")==0) )
      rlen=WILDCARD;
    else if ( is_protein(resname) && (strcmp(props[i].resname,"AAA")==0) )
      rlen=PROTEIN;
    else {
      rlen = dbcmp_len(resname, props[i].resname);
    }
    alen = dbcmp_len(atomname, props[i].atomname);
    if ( alen>=malen && rlen>=mrlen ) {
      malen=alen;
      mrlen=rlen;
      *p = props[i].prop;
      j=i;
    }
  }
  if (debug)
    fprintf(debug,"search: %4s %4s match: %4s %4s %8g\n",
	    resname, atomname, props[j].resname, props[j].atomname, *p);
  
  return (malen!=NOTFOUND && mrlen!=NOTFOUND);
}
#undef NOTFOUND
#undef WILDCARD
#undef PROTEIN

real get_mass(char *resnm, char *atomnm)
{
  real m;
  static t_prop *mass=NULL;
  static int    nmass;
  static int    maxmass;
  
  if (!mass) {
    fprintf(stderr,
	    "WARNING: masses will be determined based on residue and atom names,\n"
	    "         this can deviate from the real mass of the atom type\n");
    nmass = read_props("atommass.dat",&mass,&maxmass);
    if (debug)
      write_props(debug,nmass,mass);
  }
  
  if ( ! get_prop(&m, resnm, atomnm, nmass, mass) ) {
    m=12.0110; /* carbon mass */
    add_prop(m, resnm, atomnm, &nmass, &mass, &maxmass);
    fprintf(stderr,"Mass of atom %s %s set to %g\n",
	    resnm,atomnm,m);
  }
  
  return m;
}

real get_vdw(char *resnm, char *atomnm, real default_r)
{
  real r;
  static t_prop *vdwr=NULL;
  static int    nvdwr;
  static int    maxvdwr;
  
  if (!vdwr) {
    nvdwr = read_props("vdwradii.dat",&vdwr,&maxvdwr);
    if (debug)
      write_props(debug,nvdwr,vdwr);
  }
  
 if ( ! get_prop(&r, resnm, atomnm, nvdwr, vdwr) ) {
    r=default_r;
    add_prop(r, resnm, atomnm, &nvdwr, &vdwr, &maxvdwr);
    fprintf(stderr,"Van der Waals radius of atom %s %s set to %g\n",
	    resnm,atomnm,r);
  }
  return r;
}

real get_dgsolv(char *resnm, char *atomnm, real default_s)
{
  static t_prop *dgs=NULL;
  static int    ndgs;
  static int    maxdgs;
  real   r;
  int    i;
  
  if (!dgs) {
    ndgs = read_props("dgsolv.dat",&dgs,&maxdgs);
    please_cite(stderr,"Wang2001a");
    /* Convert to MD units */
    for(i=0; (i<ndgs); i++)
      dgs[i].prop *= 418.4;
    
    if (debug) {
      fprintf(debug,"Atomic solvation factors in kJ/mol nm^2\n");
      write_props(debug,ndgs,dgs);
    }
  }
  
  if ( ! get_prop(&r, resnm, atomnm, ndgs, dgs) ) {
    r = default_s;
    add_prop(r, resnm, atomnm, &ndgs, &dgs, &maxdgs);
    fprintf(stderr,"Solvation factor of atom %s %s set to %g\n",
	    resnm,atomnm,r);
  }
  return r;
}
