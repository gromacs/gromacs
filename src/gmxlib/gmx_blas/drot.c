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
#include "gmx_blas.h"

void
F77_FUNC(drot,DROT)(int *n__,
      double *dx,
      int *incx__,
      double *dy,
      int *incy__,
      double *c__,
      double *s__)
{
  int i,ix,iy;
  double dtemp;

  int n = *n__;
  int incx = *incx__;
  int incy = *incy__;
  double c = *c__;
  double s = *s__;
  
  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*(incx);
    if(incy<0)
      iy = (1-n)*(incy);
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) {
      dtemp  = (c) * dx[ix] + (s) * dy[iy];
      dy[iy] = (c) * dy[iy] - (s) * dx[ix];
      dx[ix] = dtemp;
    }

    return;

  } else {

    /* unit increments */   
    for(i=0;i<n;i++) {
      dtemp = (c) * dx[i] + (s) * dy[i];
      dy[i] = (c) * dy[i] - (s) * dx[i];
      dx[i] = dtemp;      
    }

  }
}
