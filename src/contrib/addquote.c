/*
 * This file is part of the GROMACS molecular simulation package,
 * version 4.6.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
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
