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

cl_R Slater_2S_2S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (93LL*xi)/256LL

            ;
        }
        else
        {
            S = (1LL/r)*((-80640LL + 80640LL*exp(2LL*rxi) - 131985LL*rxi - 102690LL*Power(rxi, 2LL) -

                          49980LL*Power(rxi, 3LL) - 16800LL*Power(rxi, 4LL) - 4032LL*Power(rxi, 5LL) -

                          672LL*Power(rxi, 6LL) - 64LL*Power(rxi, 7LL))/(80640LL*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 6LL) + 7LL*Power(xi, 5LL)*xj + 21LL*Power(xi, 4LL)*Power(xj, 2LL) +

                        35LL*Power(xi, 3LL)*Power(xj, 3LL) + 21LL*Power(xi, 2LL)*Power(xj, 4LL) +

                        7LL*xi*Power(xj, 5LL) + Power(xj, 6LL)))/(2LL*Power(xi + xj, 7LL))

            ;
        }
        else
        {
            S = (1LL/r)*((6LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 7LL) -

                          exp(2LL*rxi)*Power(rxi, 6LL)*

                          (21LL*Power(rxi, 4LL)*Power(rxj, 4LL)*(6LL + 11LL*rxj + 2LL*Power(rxj, 2LL)) -

                           2LL*Power(rxj, 8LL)*(90LL + 54LL*rxj + 12LL*Power(rxj, 2LL) + Power(rxj, 3LL)) +

                           Power(rxi, 8LL)*(6LL + 9LL*rxj + 6LL*Power(rxj, 2LL) + 2LL*Power(rxj, 3LL)) +

                           Power(rxi, 2LL)*Power(rxj, 6LL)*

                           (-390LL - 69LL*rxj + 18LL*Power(rxj, 2LL) + 4LL*Power(rxj, 3LL)) -

                           Power(rxi, 6LL)*Power(rxj, 2LL)*

                           (42LL + 63LL*rxj + 42LL*Power(rxj, 2LL) + 4LL*Power(rxj, 3LL))) +

                          exp(2LL*rxj)*Power(rxj, 6LL)*

                          (-24LL*Power(rxi, 10LL) - 2LL*Power(rxi, 11LL) - 69LL*Power(rxi, 7LL)*Power(rxj, 2LL) +

                           6LL*Power(rxj, 8LL) + 9LL*rxi*Power(rxj, 8LL) +

                           4LL*Power(rxi, 9LL)*(-27LL + Power(rxj, 2LL)) +

                           18LL*Power(rxi, 8LL)*(-10LL + Power(rxj, 2LL)) +

                           6LL*Power(rxi, 2LL)*Power(rxj, 6LL)*(-7LL + Power(rxj, 2LL)) -

                           42LL*Power(rxi, 4LL)*Power(rxj, 4LL)*(-3LL + Power(rxj, 2LL)) +

                           Power(rxi, 3LL)*Power(rxj, 6LL)*(-63LL + 2LL*Power(rxj, 2LL)) +

                           6LL*Power(rxi, 6LL)*Power(rxj, 2LL)*(-65LL + 7LL*Power(rxj, 2LL)) +

                           Power(rxi, 5LL)*(231LL*Power(rxj, 4LL) - 4LL*Power(rxj, 6LL))))/

                         (6LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 7LL)*Power(rxi + rxj, 7LL))

                         );
        }

    }
    return S;
}
