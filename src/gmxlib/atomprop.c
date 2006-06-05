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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "atomprop.h"
#include "macros.h"
#include "index.h"
#include "strdb.h"
#include "copyrite.h"

typedef struct {
  int    nprop,maxprop;
  char   *db;
  double def;
  char   **atomnm;
  char   **resnm;
  bool   *bAvail;
  real   *value;
} t_props;

typedef struct {
  t_props    props[epropNR];
  t_aa_names *aan;
} t_atomprop;

/* NOTFOUND should be smallest, others larger in increasing priority */
enum { NOTFOUND=-4, WILDCARD, WILDPROT, PROTEIN };

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

static int get_prop_index(t_props *ap,t_aa_names *aan,
			  char *resnm,char *atomnm,
			  bool *bExact)
{
  int  i,j=NOTFOUND,alen,rlen,malen,mrlen;
  bool bProtein,bProtWild;
  
  bProtein  = is_protein(aan,resnm);
  bProtWild = (strcmp(resnm,"AAA")==0);
  malen = NOTFOUND;
  mrlen = NOTFOUND;
  for(i=0; (i<ap->nprop); i++) {
    if ( (strcmp(ap->resnm[i],"*")==0) ||
	 (strcmp(ap->resnm[i],"???")==0) )
      rlen=WILDCARD;
    else if (strcmp(ap->resnm[i],"AAA")==0)
      rlen=WILDPROT;
    else {
      rlen = dbcmp_len(resnm, ap->resnm[i]);
    }
    alen = dbcmp_len(atomnm, ap->atomnm[i]);
    if ( (alen > NOTFOUND) && (rlen > NOTFOUND)) {
      if ( (alen>=malen) && (rlen>=mrlen) ) {
	malen = alen;
	mrlen = rlen;
	j     = i;
      }
    }
  }
  
  *bExact = ((malen == strlen(atomnm)) &&
	     ((mrlen == strlen(resnm)) || 
	      ((mrlen == WILDPROT) && bProtWild) ||
	      ((mrlen == WILDCARD) && !bProtein && !bProtWild)));
  
  if (debug) {
    fprintf(debug,"searching residue: %4s atom: %4s\n",resnm,atomnm);
    if (j == NOTFOUND)
      fprintf(debug," not succesful\n");
    else
      fprintf(debug," match: %4s %4s\n",ap->resnm[j],ap->atomnm[j]);
  }
  return j;
}

static void add_prop(t_props *ap,t_aa_names *aan,
		     char *resnm,char *atomnm,
		     real p,int line) 
{
  int  i,j;
  bool bExact;
  
  j = get_prop_index(ap,aan,resnm,atomnm,&bExact);
  
  if (!bExact) {
    if (ap->nprop >= ap->maxprop) {
      ap->maxprop += 10;
      srenew(ap->resnm,ap->maxprop);
      srenew(ap->atomnm,ap->maxprop);
      srenew(ap->value,ap->maxprop);
      srenew(ap->bAvail,ap->maxprop);
      for(i=ap->nprop; (i<ap->maxprop); i++) {
	ap->atomnm[i] = NULL;
	ap->resnm[i]  = NULL;
	ap->value[i]  = 0;
	ap->bAvail[i] = FALSE;
      }
    }
    
    ap->atomnm[ap->nprop] = strdup(atomnm);
    ap->resnm[ap->nprop]  = strdup(resnm);
    j = ap->nprop;
    ap->nprop++;
  }
  if (ap->bAvail[j]) {
    if (ap->value[j] == p)
      fprintf(stderr,"Warning double identical entries for %s %s %g on line %d in file %s\n",
	      resnm,atomnm,p,line,ap->db);
    else {
      fprintf(stderr,"Warning double different entries %s %s %g and %g on line %d in file %s\n"
	      "Using last entry (%g)\n",
	      resnm,atomnm,p,ap->value[j],line,ap->db,p);
      ap->value[j] = p;
    }
  }
  else {
    ap->bAvail[j] = TRUE;
    ap->value[j]  = p;
  }
}

static void read_props(t_atomprop *ap,int eprop,double factor)
{
  FILE   *fp;
  char   line[STRLEN],resnm[32],atomnm[32];
  double pp;
  int    line_no;
  
  fp      = libopen(ap->props[eprop].db);
  line_no = 0;
  while(get_a_line(fp,line,STRLEN)) {
    line_no++;
    if (sscanf(line,"%s %s %lf",resnm,atomnm,&pp) == 3) {
      pp *= factor;
      add_prop(&(ap->props[eprop]),ap->aan,resnm,atomnm,pp,line_no);
    }
    else 
      fprintf(stderr,"WARNING: Error in file %s at line %d ignored\n",
	      ap->props[eprop].db,line_no);
  }
  fclose(fp);
}

void *get_atomprop(void) 
{
  char *fns[epropNR]  = { "atommass.dat", "vdwradii.dat", "dgsolv.dat" };
  double fac[epropNR] = { 1.0,    1.0,  418.4 };
  double def[epropNR] = { 12.011, 0.14, 0.0 };

  t_atomprop *ap;
  int i;
  
  fprintf(stdout,
	  "WARNING: masses will be determined based on residue and atom names,\n"
	  "         this can deviate from the real mass of the atom type\n");
  
  snew(ap,1);

  ap->aan = get_aa_names();
  for(i=0; (i<epropNR); i++) {
    ap->props[i].db  = strdup(fns[i]);
    ap->props[i].def = def[i];
    read_props(ap,i,fac[i]);
  }
  printf("#Entries in");
  for(i=0; (i<epropNR); i++) 
    printf(" %s: %d",ap->props[i].db,ap->props[i].nprop);
  printf("\n");
  
  return (void *)ap;
}

void done_atomprop(void **atomprop)
{
  t_atomprop **ap = (t_atomprop **) atomprop;
}

bool query_atomprop(void *atomprop,int eprop,char *resnm,char *atomnm,
		    real *value)
{
  t_atomprop *ap = (t_atomprop *) atomprop;
  int  i,j;
  char *atomname;
  bool bExact;
  
  if (isdigit(atomnm[0])) {
    /* put digit after atomname */
    snew(atomname,strlen(atomnm)+1);
    for (i=1; (i<strlen(atomnm)); i++)
      atomname[i-1] = atomnm[i];
    atomname[i++] = atomnm[0];
    atomname[i]   = '\0';
  } 
  else 
    atomname = atomnm;

  j = get_prop_index(&(ap->props[eprop]),ap->aan,resnm,atomname,&bExact);
  
  if (j >= 0) {
    *value = ap->props[eprop].value[j];
    return TRUE;
  }
  else {
    *value = ap->props[eprop].def;
    return FALSE;
  }
}

