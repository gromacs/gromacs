/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
#include <stdio.h>
#include <stdlib.h>

#include "slater_integrals.h"

int main(int argc, char *argv[])
{
    double xi, xj, r, dx, S;
    double S11, S12, S13, S22, S23, S33;
    char   buf[256];
    int    i, nmax = 1000;

    if (argc < 4)
    {
        fprintf(stderr, "Usage: %s xi xj r\n", argv[0]);
        exit(1);
    }
    xi = atof(argv[1]);
    xj = atof(argv[2]);
    dx = 0.001;
    r  = 0;
    for (i = 0; (i <= nmax); i++)
    {
        S11 = Coulomb_SS(r, 1, 1, xi, xj);
        S12 = Coulomb_SS(r, 1, 2, xi, xj);
        S13 = Coulomb_SS(r, 1, 3, xi, xj);
        S22 = Coulomb_SS(r, 2, 2, xi, xj);
        S23 = Coulomb_SS(r, 2, 3, xi, xj);
        S33 = Coulomb_SS(r, 3, 3, xi, xj);
        printf("%8.3f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f\n",
               r, S11, S12, S13, S22, S23, S33);
        r += dx;
    }
    return 0;
}
