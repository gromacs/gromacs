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

cl_R Slater_1S_2S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (7LL*xi)/16LL

            ;
        }
        else
        {
            S = (1LL/r)*((-240LL + 240LL*exp(2LL*rxi) - 375LL*rxi - 270LL*Power(rxi, 2LL) - 115LL*Power(rxi, 3LL) -

                          30LL*Power(rxi, 4LL) - 4LL*Power(rxi, 5LL))/(240LL*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 4LL) + 5LL*Power(xi, 3LL)*xj + 10LL*Power(xi, 2LL)*Power(xj, 2LL) +

                        10LL*xi*Power(xj, 3LL) + 2LL*Power(xj, 4LL)))/(2LL*Power(xi + xj, 5LL))

            ;
        }
        else
        {
            S = (1LL/r)*((6LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 5LL) +

                          6LL*exp(2LL*rxj)*Power(rxj, 6LL)*

                          (-4LL*Power(rxi, 4LL) - Power(rxi, 5LL) - 5LL*Power(rxi, 2LL)*Power(rxj, 2LL) +

                           Power(rxj, 4LL) + rxi*Power(rxj, 4LL)) -

                          exp(2LL*rxi)*Power(rxi, 4LL)*

                          (Power(rxi, 6LL)*(6LL + 9LL*rxj + 6LL*Power(rxj, 2LL) + 2LL*Power(rxj, 3LL)) -

                           3LL*Power(rxi, 4LL)*Power(rxj, 2LL)*

                           (10LL + 15LL*rxj + 10LL*Power(rxj, 2LL) + 2LL*Power(rxj, 3LL)) +

                           3LL*Power(rxi, 2LL)*Power(rxj, 4LL)*

                           (20LL + 33LL*rxj + 14LL*Power(rxj, 2LL) + 2LL*Power(rxj, 3LL)) -

                           Power(rxj, 6LL)*(84LL + 63LL*rxj + 18LL*Power(rxj, 2LL) + 2LL*Power(rxj, 3LL))))/

                         (6LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 5LL)*Power(rxi + rxj, 5LL))

                         );
        }

    }
    return S;
}


cl_R Slater_2S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_1S_2S(r, xj, xi);
}
