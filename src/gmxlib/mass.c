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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_mass_c = "$Id$";
#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "mass.h"
#include "macros.h"
#include "index.h"

#define ATOMMASS "atommass.dat"

typedef struct {
  char atomname[10];
  char resname[10];
  real mass;
} t_mass;

static int read_mass(char *massdata,t_mass **mass) 
{
  FILE   *fp;
  char   line[STRLEN];
  t_mass *mdum;
  real   mm;
  int    i;

  fprintf(stderr,
	  "WARNING: masses will be determined based on residue and atom name,\n"
	  "         this can deviate from the real mass of the atom type\n");
  fp=libopen(massdata);
  *mass=NULL;
  i=0;
  while(fgets(line,STRLEN,fp)) {
    srenew(*mass,i+1);
    mdum=*mass;
    if (sscanf(line,"%s %s %f",mdum[i].resname,mdum[i].atomname,&mm)==3) {
      mdum[i].mass=mm;
      i++;
    } else
      if (debug)
	fprintf(stderr,"ERROR in file %s at line %d\n",massdata,i);
  }
  fclose(fp);
  
  return i;
} /*read_mass()*/

real get_mass(char *resnm, char *atomnm)
{
  int i,j,best,mlen,len;
  char atomname[10];
  static t_mass *mass;
  static int    nmass;
  static bool   bRead;
  
  if (!bRead) {
    nmass = read_mass(ATOMMASS,&mass);
    bRead=TRUE;
  }
  
  best=-1;
  mlen=0;
  if (isdigit(atomnm[0])) {
    /* put digit after atomname */
    for (i=1; (i<strlen(atomnm)); i++)
      atomname[i-1]=atomnm[i];
    atomname[i++]=atomnm[0];
    atomname[i]='\0';
  } else 
    strcpy(atomname,atomnm);
      
  for(i=0; (i<nmass); i++)
    if ((strcmp(mass[i].resname,"*")==0) ||
	( (strncmp(mass[i].resname,resnm,3)==0) ) ||
	( (strcmp(mass[i].resname,"AAA")==0) && is_protein(resnm) )) {
      len=min((int) strlen(atomname),(int) strlen(mass[i].atomname));
      for(j=0; (j<len); j++) 
	if (atomname[j] != mass[i].atomname[j])
	  break;
      if (j >= mlen) {
	mlen=j;
	best=i;
      }
    }

  if (best == -1) {
    if (debug)
      fprintf(stderr,"Could not find atom %s %s in %s\n",
	      resnm,atomnm,ATOMMASS);
    return 0.0;
  } else
    return mass[best].mass;
} /* get_mass() */
