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
 * Grunge ROck MAChoS
 */
static char *SRCID_readcmap_c = "$Id$";

#include "readcmap.h"
#include "smalloc.h"
#include "futil.h"
#include "fatal.h"
#include "string2.h"

int searchcmap(int n,t_mapping map[],char c)
{
  int i;
  
  for(i=0; (i<n); i++) {
    if (map[i].code == c)
      return i;
  }
  return -1;
}

int getcmap(FILE *in,char *fn,t_mapping **map)
{
  int       i,n;
  char      line[STRLEN];
  char      code,desc[STRLEN];
  double    r,g,b;
  t_mapping *m;
  
  if (fgets2(line,STRLEN-1,in) == NULL)
    fatal_error(0,"Not enough lines in colormap file %s"
		"(just wanted to read number of entries)",fn);
  sscanf(line,"%d",&n);
  snew(m,n);
  for(i=0; (i<n); i++) {
    if (fgets2(line,STRLEN-1,in) == NULL)
      fatal_error(0,"Not enough lines in colormap file %s"
		  "(should be %d, found only %d)",fn,n+1,i);
    sscanf(line,"%c%s%lf%lf%lf",&code,desc,&r,&g,&b);
    m[i].code=code;
    m[i].desc=strdup(desc);
    m[i].rgb.r=r;
    m[i].rgb.g=g;
    m[i].rgb.b=b;
  }
  *map=m;
  
  return n;
}

int readcmap(char *fn,t_mapping **map)
{
  FILE      *in;
  int       n;
  t_mapping *m;
  
  in=libopen(fn);
  n=getcmap(in,fn,map);
  fclose(in);
  
  return n;
}

void read2cmap(char *fn,t_mapping **map1,int *n1,t_mapping **map2,int *n2)
{
  FILE      *in;
  
  in=libopen(fn);
  *n1=getcmap(in,fn,map1);
  *n2=getcmap(in,fn,map2);
  fclose(in);
}

void printcmap(FILE *out,int n,t_mapping map[])
{
  int i;
  
  fprintf(out,"%d\n",n);
  for(i=0; (i<n); i++)
    fprintf(out,"%c  %20s  %10g  %10g  %10g\n",map[i].code,map[i].desc,
	    map[i].rgb.r,map[i].rgb.g,map[i].rgb.b);
}

void writecmap(char *fn,int n,t_mapping map[])
{
  FILE *out;
  
  out=ffopen(fn,"w");
  printcmap(out,n,map);
  fclose(out);
}
