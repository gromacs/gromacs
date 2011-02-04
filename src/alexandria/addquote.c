/*
 * $Id$
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include "strdb.h"
#include "copyrite.h"
#include "smalloc.h"

void swap_str(char *str)
{
  int i;
  for(i=0; (i<strlen(str)); i++)
    str[i] = ~str[i];
}

void add_quote(char *q)
{
  FILE *fp;
  int  i,n;
  char **str = NULL;
  char c,*db   = "gurgle.dat";
  
  do {
    fprintf(stderr,"Add quote '%s' (y/n)? ",q);
    c = toupper(fgetc(stdin));
  } while ((c != 'Y') && (c != 'N'));
  if (c == 'Y') {
    n = get_strings(db,&str);
    srenew(str,n+1);
    str[n] = strdup(q);
    swap_str(str[n]);
    n++;
    fp = fopen(db,"w");
    fprintf(fp,"%d\n",n);
    for(i=0; (i<n); i++) 
      fprintf(fp,"%s\n",str[i]);
    fclose(fp);
  }
}

void dump_quotes()
{
  FILE *fp;
  int  i,j,n;
  char **str = NULL;
  char *db   = "gurgle.dat";
  
  n = get_strings(db,&str);
  for(j=0; (j<n); j++) {
    swap_str(str[j]);
    printf("%s\n",str[j]);
  }
}

void undump_quotes(char *fn)
{
  FILE *fp;
  int  i,j,n;
  char **str = NULL;
  char *db   = "gurgle.dat";

  n = get_lines(fn,&str);  
  fp = fopen(db,"w");
  fprintf(fp,"%d\n",n);
  for(j=0; (j<n); j++) {
    swap_str(str[j]);
    fprintf(fp,"%s\n",str[j]);
  }
  fclose(fp);
}

int main(int argc,char *argv[])
{
  int  i;
  char c;
  
  for(i=1; (i<argc); ) {
    if (strcmp(argv[i],"-add") == 0) {
      if (i < argc-1) {
	add_quote(argv[i+1]);
	i+=2;
      }
      else {
	fprintf(stderr,"Expected argument after -add\n");
	exit(1);
      }
    }
    else if (strcmp(argv[i],"-dump") == 0) {
      if (i < argc-1) {
	dump_quotes();
	i += 2;
      }
      else {
	fprintf(stderr,"Expected argument after -dump\n");
	exit(1);
      }
    }
    else if (strcmp(argv[i],"-undump") == 0) {
      if (i < argc-1) {
	undump_quotes(argv[i+1]);
	i += 2;
      }
      else {
	fprintf(stderr,"Expected argument after -undump\n");
	exit(1);
      }
    }
    else {
      fprintf(stderr,"WTF: %s\n",argv[i]);
      i++;
    }
  }
  thanx(stdout);
  
  return 0;
}
