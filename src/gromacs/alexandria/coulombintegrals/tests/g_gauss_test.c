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

#include "gaussian_integrals.h"

int main(int argc, char *argv[])
{
    FILE  *fp;
    int    i, n = 1000;
    double xi;
    double r;

    fp = fopen("test_gauss.xvg", "w");
    i  = 0;
    for (xi = 20; (xi <= 100); xi += 20)
    {
        fprintf(fp, "@ s%d legend \"Nuclear xi = %g\"\n", i++, xi);
    }
    for (xi = 20; (xi <= 100); xi += 20)
    {
        fprintf(fp, "@ s%d legend \"Coulomb xi = %g  xj = %g\"\n", i++, xi, xi);
    }
    for (xi = 20; (xi <= 100); xi += 20)
    {
        fprintf(fp, "@type xy\n");
        for (i = 0; (i < n); i++)
        {
            r = i*1.0/n;
            fprintf(fp, "%12e  %12e\n", r, Nuclear_GG(r, xi));
        }
        fprintf(fp, "&\n");
    }
    for (xi = 20; (xi <= 100); xi += 20)
    {
        fprintf(fp, "@type xy\n");
        for (i = 0; (i < n); i++)
        {
            r = i*1.0/n;
            fprintf(fp, "%12e  %12e\n", r, Coulomb_GG(r, xi, xi));
        }
        fprintf(fp, "&\n");
    }
    fclose(fp);

    return 0;
}
