/*
 * $Id$
 * 
 *                This source code ispart of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Great Red Owns Many ACres of Sand 
 */

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
  bool bAvail;
  real val;
} t_prop;

typedef struct {
  char   *atomnm;
  char   *resnm;
  t_prop p[epropNR];
} t_props;

typedef struct {
  int        nprop,maxprop;
  char       *db[epropNR];
  double     def[epropNR];
  t_props    *props;
  t_aa_names *aan;
} t_atomprop;

/* NOTFOUND should be smallest, others larger in increasing priority */
enum { NOTFOUND=-3, WILDCARD, PROTEIN };

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

static int get_prop_index(t_atomprop *ap,char *resnm,char *atomnm)
{
  int  i,j=NOTFOUND,alen,rlen,malen,mrlen;
  
  malen = NOTFOUND;
  mrlen = NOTFOUND;
  for(i=0; (i<ap->nprop); i++) {
    if ( (strcmp(ap->props[i].resnm,"*")==0) ||
	 (strcmp(ap->props[i].resnm,"???")==0) )
      rlen=WILDCARD;
    else if ( is_protein(ap->aan,resnm) && 
	      (strcmp(ap->props[i].resnm,"AAA")==0) )
      rlen=PROTEIN;
    else {
      rlen = dbcmp_len(resnm, ap->props[i].resnm);
    }
    alen = dbcmp_len(atomnm, ap->props[i].atomnm);
    if ( (alen>=malen) && (rlen>=mrlen) ) {
      malen = alen;
      mrlen = rlen;
      j     = i;
    }
  }
  
  if (debug)
    fprintf(debug,"search: %4s %4s match: %4s %4s\n",
	    resnm,atomnm, ap->props[j].resnm, ap->props[j].atomnm);
  
  return j;
}

static void add_prop(t_atomprop *ap,int eprop,char *resnm,char *atomnm,
		     real p,int line) 
{
  int i,j;
  
  j = get_prop_index(ap,resnm,atomnm);
  
  if (j < 0) {
    if (ap->nprop >= ap->maxprop) {
      ap->maxprop+=10;
      srenew(ap->props,ap->maxprop);
      for(i=ap->nprop; (i<ap->maxprop); i++) {
	ap->props[i].atomnm = NULL;
	ap->props[i].resnm  = NULL;
	memset(ap->props[i].p,0,epropNR*sizeof(ap->props[i].p[0]));
      }
    }
    
    ap->props[ap->nprop].atomnm = strdup(atomnm);
    ap->props[ap->nprop].resnm  = strdup(resnm);
    j = ap->nprop;
    ap->nprop++;
  }
  if (ap->props[j].p[eprop].bAvail) {
    if (ap->props[j].p[eprop].val == p)
      fprintf(stderr,"Warning double identical entries for %s %s %g on line %d in file %s\n",
	      resnm,atomnm,p,line,ap->db[eprop]);
    else {
      fprintf(stderr,"Warning double different entries %s %s %g and %g on line %d in file %s\n"
	      "Using last entry (%g)\n",
	      resnm,atomnm,p,ap->props[j].p[eprop].val,line,ap->db[eprop],p);
      ap->props[j].p[eprop].val    = p;
    }
  }
  else {
    ap->props[j].p[eprop].bAvail = TRUE;
    ap->props[j].p[eprop].val    = p;
  }
}

static void read_props(t_atomprop *ap,int eprop,double factor)
{
  FILE   *fp;
  char   line[STRLEN],resnm[32],atomnm[32];
  double pp;
  int    line_no;
  
  fp      = libopen(ap->db[eprop]);
  line_no = 0;
  while(get_a_line(fp,line,STRLEN)) {
    if (sscanf(line,"%s %s %lf",resnm,atomnm,&pp) == 3) {
      pp *= factor;
      add_prop(ap,eprop,resnm,atomnm,pp,line_no);
    }
    else 
      fprintf(stderr,"WARNING: Error in file %s at line %d ignored\n",
	      ap->db[eprop],line_no);
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
  fprintf(stdout,"In case you use free energy of solvation predictions:\n");
  please_cite(stdout,"Eisenberg86a");
  
  snew(ap,1);

  ap->aan = get_aa_names();
  for(i=0; (i<epropNR); i++) {
    ap->db[i] = strdup(fns[i]);
    ap->def[i] = def[i];
    read_props(ap,i,fac[i]);
  }
    
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

  j = get_prop_index(ap,resnm,atomname);
    
  if ((j >= 0)  && (ap->props[j].p[eprop].bAvail)) {
    *value = ap->props[j].p[eprop].val;
    return TRUE;
  }
  else {
    *value = ap->def[eprop];
    return FALSE;
  }
}

