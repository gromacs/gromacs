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

#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "mass.h"
#include "macros.h"

int read_mass(char *massdata,t_mass **mass) 
{
  FILE   *fp;
  t_mass *mdum;
  double mm;
  char   line[STRLEN],name[6];
  int    i;

  fprintf(stderr,
	  "WARNING: masses will be determined based on atom name,\n"
	  "         this can deviate from the real mass of the atom type\n");
  fp=libopen(massdata);
  *mass=NULL;
  for (i=0; fgets(line,STRLEN,fp); i++) {
    srenew(*mass,i+1);
    mdum=*mass;
    strncpy(name,line,6);
    name[5]='\0';
    if ((int)strlen(name) < 2) 
      break; 
    sscanf(name,"%s",mdum[i].atomname);
    sscanf(line+5,"%lf",&mm);
    mdum[i].mass=mm;
  }
  fclose(fp);
  
  return i;
} /*read_mass()*/

void write_mass(char *massdata,t_mass mass[],int nmass)
{
  FILE *fp;
  int i;

  fp=ffopen(massdata,"w");
  for(i=0;(i<nmass);i++) {
    fprintf(fp,"%-5s%8.3f\n",mass[i].atomname,mass[i].mass);
  }
  fclose(fp);
}

real get_mass(int nmass,t_mass mass[],char *atom)
{
  int i,j,best,mlen,len;

  best=-1;
  mlen=0;
  for(i=0; (i<nmass); i++) {
    len=min((int) strlen(atom),(int) strlen(mass[i].atomname));
    for(j=0; (j<len); j++) 
      if (atom[j] != mass[i].atomname[j])
	break;
    if (j > mlen) {
      mlen=j;
      best=i;
    }
  }

  if (best == -1) {
    if (debug)
      fprintf(stderr,"Could not find atom %s\n",atom);
    return 0.0;
  }
  else
    return mass[best].mass;
}
