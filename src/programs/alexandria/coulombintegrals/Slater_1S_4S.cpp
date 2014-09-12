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

cl_R Slater_1S_4S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (253LL*xi)/1024LL

            ;
        }
        else
        {
            S = (1LL/r)*((-2903040LL + 2903040LL*exp(2LL*rxi) - 5088825LL*rxi - 4371570LL*Power(rxi, 2LL) -

                          2439990LL*Power(rxi, 3LL) - 986580LL*Power(rxi, 4LL) - 303912LL*Power(rxi, 5LL) -

                          72576LL*Power(rxi, 6LL) - 13248LL*Power(rxi, 7LL) - 1728LL*Power(rxi, 8LL) -

                          128LL*Power(rxi, 9LL))/(2.90304e6*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 8LL) + 9LL*Power(xi, 7LL)*xj + 36LL*Power(xi, 6LL)*Power(xj, 2LL) +

                        84LL*Power(xi, 5LL)*Power(xj, 3LL) + 126LL*Power(xi, 4LL)*Power(xj, 4LL) +

                        126LL*Power(xi, 3LL)*Power(xj, 5LL) + 84LL*Power(xi, 2LL)*Power(xj, 6LL) +

                        36LL*xi*Power(xj, 7LL) + 4LL*Power(xj, 8LL)))/(4LL*Power(xi + xj, 9LL))

            ;
        }
        else
        {
            S = (1LL/r)*((1260LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 9LL) +

                          1260LL*exp(2LL*rxj)*Power(rxj, 10LL)*

                          (-6LL*Power(rxi, 8LL) - Power(rxi, 9LL) - 51LL*Power(rxi, 6LL)*Power(rxj, 2LL) -

                           6LL*Power(rxi, 7LL)*Power(rxj, 2LL) - 63LL*Power(rxi, 4LL)*Power(rxj, 4LL) -

                           9LL*Power(rxi, 2LL)*Power(rxj, 6LL) + 6LL*Power(rxi, 3LL)*Power(rxj, 6LL) +

                           Power(rxj, 8LL) + rxi*Power(rxj, 8LL)) -

                          exp(2LL*rxi)*Power(rxi, 4LL)*

                          (42LL*Power(rxi, 10LL)*Power(rxj, 4LL)*

                           (1080LL + 1890LL*rxj + 1620LL*Power(rxj, 2LL) + 900LL*Power(rxj, 3LL) +

           360LL*Power(rxj, 4LL) + 111LL*Power(rxj, 5LL) + 22LL*Power(rxj, 6LL) +

           2LL*Power(rxj, 7LL)) - 70LL*Power(rxi, 8LL)*Power(rxj, 6LL)*

                           (1512LL + 2646LL*rxj + 2268LL*Power(rxj, 2LL) + 1248LL*Power(rxj, 3LL) +

           528LL*Power(rxj, 4LL) + 153LL*Power(rxj, 5LL) + 26LL*Power(rxj, 6LL) +

           2LL*Power(rxj, 7LL)) + 14LL*Power(rxi, 2LL)*Power(rxj, 12LL)*

                           (2970LL + 16335LL*rxj + 15390LL*Power(rxj, 2LL) + 7110LL*Power(rxj, 3LL) +

           1980LL*Power(rxj, 4LL) + 351LL*Power(rxj, 5LL) + 38LL*Power(rxj, 6LL) +

           2LL*Power(rxj, 7LL)) - 2LL*Power(rxj, 14LL)*

                           (62370LL + 72765LL*rxj + 39690LL*Power(rxj, 2LL) + 13230LL*Power(rxj, 3LL) +

           2940LL*Power(rxj, 4LL) + 441LL*Power(rxj, 5LL) + 42LL*Power(rxj, 6LL) +

           2LL*Power(rxj, 7LL)) + Power(rxi, 14LL)*

                           (1260LL + 2205LL*rxj + 1890LL*Power(rxj, 2LL) + 1050LL*Power(rxj, 3LL) +

           420LL*Power(rxj, 4LL) + 126LL*Power(rxj, 5LL) + 28LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL)) - 7LL*Power(rxi, 12LL)*Power(rxj, 2LL)*

                           (1620LL + 2835LL*rxj + 2430LL*Power(rxj, 2LL) + 1350LL*Power(rxj, 3LL) +

           540LL*Power(rxj, 4LL) + 162LL*Power(rxj, 5LL) + 36LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL)) + 35LL*Power(rxi, 6LL)*Power(rxj, 8LL)*

                           (4536LL + 7983LL*rxj + 6534LL*Power(rxj, 2LL) + 4014LL*Power(rxj, 3LL) +

           1644LL*Power(rxj, 4LL) + 414LL*Power(rxj, 5LL) + 60LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL)) - 21LL*Power(rxi, 4LL)*Power(rxj, 10LL)*

                           (7920LL + 11385LL*rxj + 12330LL*Power(rxj, 2LL) + 7410LL*Power(rxj, 3LL) +

           2580LL*Power(rxj, 4LL) + 546LL*Power(rxj, 5LL) + 68LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL))))/

                         (1260LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 9LL)*Power(rxi + rxj, 9LL))

                         );
        }

    }
    return S;
}


cl_R Slater_4S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_1S_4S(r, xj, xi);
}
