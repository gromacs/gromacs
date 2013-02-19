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


#include <types/simple.h>
#include "gmx_blas.h"

float
F77_FUNC(snrm2,SNRM2)(int  *     n__,
                      float *    x,
                      int    *    incx__)
{
    int ix,max_ix;
    float ssq,scale,absxi,t;
    
    int n = *n__;
    int incx = *incx__;
    
    if(n<1 || incx<1)
        return 0;
    else if (n==1) {
        t = x[0];
        if(t>=0)
            return t;
        else 
            return -t;
    }
    
    scale = 0.0;
    ssq   = 1.0;
    
    max_ix = 1+(n-1)*(incx);
    for(ix=1;ix<=max_ix;ix+=incx) {
        t = x[ix-1];
        if(fabs(t)>GMX_FLOAT_MIN) {
            absxi = (t>=0) ? t : (-t);
            if(scale<absxi) {
                t = scale/absxi;
                t = t*t;
                ssq = ssq*t + 1.0;
                scale = absxi;
            } else {
                t = absxi/scale;
                ssq += t*t;
            }
        }
    }
    return scale*sqrt(ssq);
    
}


 
