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

cl_R Slater_1S_3S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (41LL*xi)/128LL

            ;
        }
        else
        {
            S = (1LL/r)*((-120960LL + 120960LL*exp(2LL*rxi) - 203175LL*rxi - 164430LL*Power(rxi, 2LL) -

                          84420LL*Power(rxi, 3LL) - 30240LL*Power(rxi, 4LL) - 7728LL*Power(rxi, 5LL) -

                          1344LL*Power(rxi, 6LL) - 128LL*Power(rxi, 7LL))/(120960LL*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 6LL) + 7LL*Power(xi, 5LL)*xj + 21LL*Power(xi, 4LL)*Power(xj, 2LL) +

                        35LL*Power(xi, 3LL)*Power(xj, 3LL) + 35LL*Power(xi, 2LL)*Power(xj, 4LL) +

                        21LL*xi*Power(xj, 5LL) + 3LL*Power(xj, 6LL)))/(3LL*Power(xi + xj, 7LL))

            ;
        }
        else
        {
            S = (1LL/r)*((45LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 7LL) +

                          15LL*exp(2LL*rxj)*Power(rxj, 8LL)*

                          (-15LL*Power(rxi, 6LL) - 3LL*Power(rxi, 7LL) - 63LL*Power(rxi, 4LL)*Power(rxj, 2LL) -

                           7LL*Power(rxi, 5LL)*Power(rxj, 2LL) - 21LL*Power(rxi, 2LL)*Power(rxj, 4LL) +

                           7LL*Power(rxi, 3LL)*Power(rxj, 4LL) + 3LL*Power(rxj, 6LL) + 3LL*rxi*Power(rxj, 6LL)) +

                          exp(2LL*rxi)*Power(rxi, 4LL)*

                          (-10LL*Power(rxi, 2LL)*Power(rxj, 8LL)*

                           (135LL + 333LL*rxj + 228LL*Power(rxj, 2LL) + 75LL*Power(rxj, 3LL) +

           13LL*Power(rxj, 4LL) + Power(rxj, 5LL)) +

                           2LL*Power(rxj, 10LL)*(945LL + 945LL*rxj + 420LL*Power(rxj, 2LL) +

                                                 105LL*Power(rxj, 3LL) + 15LL*Power(rxj, 4LL) + Power(rxj, 5LL)) -

                           Power(rxi, 10LL)*(45LL + 75LL*rxj + 60LL*Power(rxj, 2LL) + 30LL*Power(rxj, 3LL) +

                                             10LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           5LL*Power(rxi, 8LL)*Power(rxj, 2LL)*

                           (63LL + 105LL*rxj + 84LL*Power(rxj, 2LL) + 42LL*Power(rxj, 3LL) +

           14LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) -

                           5LL*Power(rxi, 6LL)*Power(rxj, 4LL)*

                           (189LL + 315LL*rxj + 252LL*Power(rxj, 2LL) + 132LL*Power(rxj, 3LL) +

           36LL*Power(rxj, 4LL) + 4LL*Power(rxj, 5LL)) +

                           5LL*Power(rxi, 4LL)*Power(rxj, 6LL)*

                           (315LL + 513LL*rxj + 468LL*Power(rxj, 2LL) + 204LL*Power(rxj, 3LL) +

           44LL*Power(rxj, 4LL) + 4LL*Power(rxj, 5LL))))/

                         (45LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 7LL)*Power(rxi + rxj, 7LL))

                         );
        }

    }
    return S;
}


cl_R Slater_3S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_1S_3S(r, xj, xi);
}
