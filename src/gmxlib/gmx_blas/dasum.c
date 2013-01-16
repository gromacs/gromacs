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
#include "gmx_blas.h"

double
F77_FUNC(dasum,DASUM)(int *n__, 
                      double *dx, 
                      int *incx__)
{
    int i__1, i__2;
    
    int i__, m, mp1;
    double dtemp;
    int nincx;
    
    int n = *n__;
    int incx = *incx__;
    
    --dx;
    
    dtemp = 0.;
    if (n <= 0 || incx <= 0) {
        return 0.0;
    }
    if (incx != 1) {
        nincx = n * incx;
        i__1 = nincx;
        i__2 = incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dtemp += fabs(dx[i__]);
        }
        return dtemp;
    }
    
    m = n % 6;
    if (m != 0) {
        i__2 = m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            dtemp += fabs(dx[i__]);
        }
        if (n < 6) {
            return dtemp;
        }
    }
    mp1 = m + 1;
    i__2 = n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
        dtemp = dtemp + fabs(dx[i__]) + fabs(dx[i__ + 1]) + 
        fabs(dx[i__ + 2]) + fabs(dx[i__+ 3]) + fabs(dx[i__ + 4]) +
        fabs(dx[i__ + 5]);
    }
    return dtemp;
}


