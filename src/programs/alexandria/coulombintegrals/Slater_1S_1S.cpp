/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include "gmxpre.h"
#include "slater_low.h"

cl_R Slater_1S_1S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (5LL*xi)/8LL

            ;
        }
        else
        {
            S = (1LL/r)*((-24LL + 24LL*exp(2LL*rxi) - 33LL*rxi - 18LL*Power(rxi, 2LL) - 4LL*Power(rxi, 3LL))/

                         (24LL*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 2LL) + 3LL*xi*xj + Power(xj, 2LL)))/Power(xi + xj, 3LL)

            ;
        }
        else
        {
            S = (1LL/r)*((exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 3LL) +

                          exp(2LL*rxj)*Power(rxj, 4LL)*

                          (-3LL*Power(rxi, 2LL) - Power(rxi, 3LL) + Power(rxj, 2LL) + rxi*Power(rxj, 2LL)) -

                          exp(2LL*rxi)*Power(rxi, 4LL)*

                          (Power(rxi, 2LL)*(1LL + rxj) - Power(rxj, 2LL)*(3LL + rxj)))/

                         (exp(2LL*(rxi + rxj))*Power(rxi - rxj, 3LL)*Power(rxi + rxj, 3LL))

                         );
        }

    }
    return S;
}
