/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#include "gmx_lapack.h"


/* LAPACK */


void
F77_FUNC(slamrg,SLAMRG)(int *n1,
                        int *n2,
                        float *a,
                        int *dtrd1,
                        int *dtrd2,
                        int *index)
{
  int n1sv = *n1;
  int n2sv = *n2;
  int i,ind1,ind2;

  if(*dtrd1>0)
    ind1 = 0;
  else
    ind1 = *n1-1;

  if(*dtrd2>0)
    ind2 = *n1;
  else
    ind2 = *n1+*n2-1;

  i = 0;
  
  while(n1sv>0 && n2sv>0) {
    if(a[ind1]<=a[ind2]) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
      n1sv--;
    } else {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
      n2sv--;
    }
  }

  if(n1sv==0) {
    for(n1sv=1;n1sv<=n2sv;n1sv++) {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
    } 
  } else {
    for(n2sv=1;n2sv<=n1sv;n2sv++) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
    } 
  }
  return;
}
