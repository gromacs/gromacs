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

cl_R Slater_2S_3S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (451LL*xi)/1536LL

            ;
        }
        else
        {
            S = (1LL/r)*((-4354560LL + 4354560LL*exp(2LL*rxi) - 7430535LL*rxi - 6151950LL*Power(rxi, 2LL) -

                          3275370LL*Power(rxi, 3LL) - 1251180LL*Power(rxi, 4LL) - 361368LL*Power(rxi, 5LL) -

                          80640LL*Power(rxi, 6LL) - 13824LL*Power(rxi, 7LL) - 1728LL*Power(rxi, 8LL) -

                          128LL*Power(rxi, 9LL))/(4.35456e6*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(2LL*Power(xi, 8LL) + 18LL*Power(xi, 7LL)*xj + 72LL*Power(xi, 6LL)*Power(xj, 2LL) +

                        168LL*Power(xi, 5LL)*Power(xj, 3LL) + 252LL*Power(xi, 4LL)*Power(xj, 4LL) +

                        252LL*Power(xi, 3LL)*Power(xj, 5LL) + 108LL*Power(xi, 2LL)*Power(xj, 6LL) +

                        27LL*xi*Power(xj, 7LL) + 3LL*Power(xj, 8LL)))/(6LL*Power(xi + xj, 9LL))

            ;
        }
        else
        {
            S = (1LL/r)*((90LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 9LL) +

                          5LL*exp(2LL*rxj)*Power(rxj, 8LL)*

                          (-90LL*Power(rxi, 12LL) - 6LL*Power(rxi, 13LL) + 18LL*Power(rxj, 10LL) +

                           27LL*rxi*Power(rxj, 10LL) + 18LL*Power(rxi, 2LL)*Power(rxj, 8LL)*

                           (-9LL + Power(rxj, 2LL)) - 162LL*Power(rxi, 4LL)*Power(rxj, 6LL)*

                           (-4LL + Power(rxj, 2LL)) - 198LL*Power(rxi, 10LL)*(5LL + Power(rxj, 2LL)) -

                           108LL*Power(rxi, 6LL)*Power(rxj, 4LL)*(36LL + Power(rxj, 2LL)) +

                           2LL*Power(rxi, 5LL)*Power(rxj, 6LL)*(675LL + Power(rxj, 2LL)) -

                           18LL*Power(rxi, 7LL)*Power(rxj, 4LL)*(-81LL + 2LL*Power(rxj, 2LL)) +

                           3LL*Power(rxi, 3LL)*Power(rxj, 8LL)*(-81LL + 2LL*Power(rxj, 2LL)) -

                           Power(rxi, 11LL)*(495LL + 2LL*Power(rxj, 2LL)) +

                           9LL*Power(rxi, 9LL)*Power(rxj, 2LL)*(-233LL + 4LL*Power(rxj, 2LL)) +

                           6LL*Power(rxi, 8LL)*Power(rxj, 2LL)*(-1063LL + 90LL*Power(rxj, 2LL))) -

                          2LL*exp(2LL*rxi)*Power(rxi, 6LL)*

                          (-90LL*Power(rxi, 6LL)*Power(rxj, 6LL)*

                           (42LL + 65LL*rxj + 76LL*Power(rxj, 2LL) + 22LL*Power(rxj, 3LL) + 2LL*Power(rxj, 4LL)) -

                           2LL*Power(rxj, 12LL)*(2970LL + 2475LL*rxj + 900LL*Power(rxj, 2LL) +

                                                 180LL*Power(rxj, 3LL) + 20LL*Power(rxj, 4LL) + Power(rxj, 5LL)) +

                           10LL*Power(rxi, 8LL)*Power(rxj, 4LL)*

                           (162LL + 270LL*rxj + 216LL*Power(rxj, 2LL) + 122LL*Power(rxj, 3LL) +

           22LL*Power(rxj, 4LL) + Power(rxj, 5LL)) -

                           5LL*Power(rxi, 4LL)*Power(rxj, 8LL)*

                           (-639LL - 3555LL*rxj - 1452LL*Power(rxj, 2LL) - 174LL*Power(rxj, 3LL) +

           6LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           Power(rxi, 12LL)*(45LL + 75LL*rxj + 60LL*Power(rxj, 2LL) + 30LL*Power(rxj, 3LL) +

                                             10LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) -

                           Power(rxi, 10LL)*Power(rxj, 2LL)*

                           (405LL + 675LL*rxj + 540LL*Power(rxj, 2LL) + 270LL*Power(rxj, 3LL) +

           90LL*Power(rxj, 4LL) + 8LL*Power(rxj, 5LL)) +

                           Power(rxi, 2LL)*Power(rxj, 10LL)*

                           (-21615LL - 9075LL*rxj - 300LL*Power(rxj, 2LL) + 490LL*Power(rxj, 3LL) +

           110LL*Power(rxj, 4LL) + 8LL*Power(rxj, 5LL))))/

                         (90LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 9LL)*Power(rxi + rxj, 9LL))

                         );
        }

    }
    return S;
}


cl_R Slater_3S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_2S_3S(r, xj, xi);
}
