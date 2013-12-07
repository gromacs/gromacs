/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "gmx_fatal.h"
#include "gstat.h"

real LegendreP(real x, unsigned long m)

{
    real polynomial = 0, x2, x3;

    switch (m)
    {
        case eacP0:
            polynomial = 1.0;
            break;
        case eacP1:
            polynomial = x;
            break;
        case eacP2:
            x2         = x*x;
            polynomial = 1.5*x2 - 0.5;
            break;
        case eacP3:
            x2         = x*x;
            polynomial = (35*x2*x2 - 30*x2 + 3)/8;
            break;
        case eacP4:
            x2         = x*x;
            x3         = x2*x;
            polynomial = (63*x3*x2 - 70*x3 + 15*x)/8;
            break;
        default:
            gmx_fatal(FARGS, "Legendre polynomials of order %d are not supported, %s %d",
                      m, __FILE__, __LINE__);
    }
    return (polynomial);
}
