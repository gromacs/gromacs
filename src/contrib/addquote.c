/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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

void add_quote(char *q)
{
  FILE *fp;
  int  i,n;
  char **str = NULL;
  char *db   = "gurgle.dat";
  
  n = get_strings(db,&str);
  srenew(str,n+1);
  snew(str[n],strlen(q)+1);
  for(i=0; (i<strlen(q)); i++)
    str[n][i] = ~q[i];
  str[n][i] = '\0';
  n++;
  fp = fopen(db,"w");
  fprintf(fp,"%d\n",n);
  for(i=0; (i<n); i++) 
    fprintf(fp,"%s\n",str[i]);
  fclose(fp);
}

int main(int argc,char *argv[])
{
  int  i;
  char c;
  
  for(i=1; (i<argc); i++) {
    do {
      fprintf(stderr,"Add quote '%s' (y/n)? ",argv[i]);
      c = toupper(fgetc(stdin));
    } while ((c != 'Y') && (c != 'N'));
    if (c == 'Y') {
      add_quote(argv[i]);
    }
  }
  thanx(stdout);
  
  return 0;
}
