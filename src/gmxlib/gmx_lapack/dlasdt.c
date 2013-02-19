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
#include <math.h>
#include "gmx_lapack.h"

void
F77_FUNC(dlasdt,DLASDT)(int *n,
	int *lvl,
	int *nd,
	int *inode,
	int *ndiml,
	int *ndimr,
	int *msub)
{
  int maxn = (*n > 1) ? *n : 1;
  double temp;
  int i,il,ir,llst,nlvl,ncrnt;

  temp = log( ((double) maxn) / ((double)(*msub+1))) / log(2.0);
  
  *lvl = 1 + (int) temp;

  i = *n / 2;
  inode[0] = i + 1;
  ndiml[0] = i;
  ndimr[0] = *n - i - 1;
  il = -1;
  ir = 0;
  llst = 1;

  for(nlvl=1;nlvl<*lvl;nlvl++) {
    for(i=0;i<llst;i++) {
      il += 2;
      ir += 2;
      ncrnt = llst + i - 1;
      ndiml[il] = ndiml[ncrnt] / 2;
      ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
      inode[il] = inode[ncrnt] - ndimr[il] - 1;
      ndiml[ir] = ndimr[ncrnt] / 2;
      ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
      inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
    }
    llst *= 2;
  }
  *nd = llst*2 - 1;
  return;
}
