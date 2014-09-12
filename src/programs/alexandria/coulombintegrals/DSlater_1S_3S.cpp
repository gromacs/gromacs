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

cl_R DSlater_1S_3S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = 0LL

            ;
        }
        else
        {
            S = -(-203175LL*xi + 241920LL*exp(2LL*r*xi)*xi - 328860LL*r*Power(xi, 2LL) -

                  253260LL*Power(r, 2LL)*Power(xi, 3LL) - 120960LL*Power(r, 3LL)*Power(xi, 4LL) -

                  38640LL*Power(r, 4LL)*Power(xi, 5LL) - 8064LL*Power(r, 5LL)*Power(xi, 6LL) -

                  896LL*Power(r, 6LL)*Power(xi, 7LL))/(120960LL*exp(2LL*r*xi)*r) +

                (-120960LL + 120960LL*exp(2LL*r*xi) - 203175LL*r*xi -

                 164430LL*Power(r, 2LL)*Power(xi, 2LL) - 84420LL*Power(r, 3LL)*Power(xi, 3LL) -

                 30240LL*Power(r, 4LL)*Power(xi, 4LL) - 7728LL*Power(r, 5LL)*Power(xi, 5LL) -

                 1344LL*Power(r, 6LL)*Power(xi, 6LL) - 128LL*Power(r, 7LL)*Power(xi, 7LL))/

                (120960LL*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-120960LL + 120960LL*exp(2LL*r*xi) - 203175LL*r*xi -

                     164430LL*Power(r, 2LL)*Power(xi, 2LL) - 84420LL*Power(r, 3LL)*Power(xi, 3LL) -

                     30240LL*Power(r, 4LL)*Power(xi, 4LL) - 7728LL*Power(r, 5LL)*Power(xi, 5LL) -

                     1344LL*Power(r, 6LL)*Power(xi, 6LL) - 128LL*Power(r, 7LL)*Power(xi, 7LL)))/

                (60480LL*exp(2LL*r*xi)*r)

            ;
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = 0LL

            ;
        }
        else
        {
            S = (45LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 7LL) +

                 15LL*exp(2LL*r*xj)*Power(xj, 8LL)*

                 (-15LL*Power(xi, 6LL) - 3LL*r*Power(xi, 7LL) - 63LL*Power(xi, 4LL)*Power(xj, 2LL) -

                  7LL*r*Power(xi, 5LL)*Power(xj, 2LL) - 21LL*Power(xi, 2LL)*Power(xj, 4LL) +

                  7LL*r*Power(xi, 3LL)*Power(xj, 4LL) + 3LL*Power(xj, 6LL) + 3LL*r*xi*Power(xj, 6LL)) +

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-10LL*Power(xi, 2LL)*Power(xj, 8LL)*

                  (135LL + 333LL*r*xj + 228LL*Power(r, 2LL)*Power(xj, 2LL) +

            75LL*Power(r, 3LL)*Power(xj, 3LL) + 13LL*Power(r, 4LL)*Power(xj, 4LL) +

            Power(r, 5LL)*Power(xj, 5LL)) +

                  2LL*Power(xj, 10LL)*(945LL + 945LL*r*xj + 420LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       105LL*Power(r, 3LL)*Power(xj, 3LL) + 15LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       Power(r, 5LL)*Power(xj, 5LL)) -

                  Power(xi, 10LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  5LL*Power(xi, 8LL)*Power(xj, 2LL)*

                  (63LL + 105LL*r*xj + 84LL*Power(r, 2LL)*Power(xj, 2LL) +

            42LL*Power(r, 3LL)*Power(xj, 3LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL) +

            2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  5LL*Power(xi, 6LL)*Power(xj, 4LL)*

                  (189LL + 315LL*r*xj + 252LL*Power(r, 2LL)*Power(xj, 2LL) +

            132LL*Power(r, 3LL)*Power(xj, 3LL) + 36LL*Power(r, 4LL)*Power(xj, 4LL) +

            4LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  5LL*Power(xi, 4LL)*Power(xj, 6LL)*

                  (315LL + 513LL*r*xj + 468LL*Power(r, 2LL)*Power(xj, 2LL) +

            204LL*Power(r, 3LL)*Power(xj, 3LL) + 44LL*Power(r, 4LL)*Power(xj, 4LL) +

            4LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (45LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 7LL)*Power(xi + xj, 7LL))

                + (2LL*(45LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 7LL) +

                        15LL*exp(2LL*r*xj)*Power(xj, 8LL)*

                        (-15LL*Power(xi, 6LL) - 3LL*r*Power(xi, 7LL) - 63LL*Power(xi, 4LL)*Power(xj, 2LL) -

                         7LL*r*Power(xi, 5LL)*Power(xj, 2LL) - 21LL*Power(xi, 2LL)*Power(xj, 4LL) +

                         7LL*r*Power(xi, 3LL)*Power(xj, 4LL) + 3LL*Power(xj, 6LL) + 3LL*r*xi*Power(xj, 6LL))

                        + exp(2LL*r*xi)*Power(xi, 4LL)*(-10LL*Power(xi, 2LL)*Power(xj, 8LL)*

                                                        (135LL + 333LL*r*xj + 228LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         75LL*Power(r, 3LL)*Power(xj, 3LL) + 13LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                         Power(r, 5LL)*Power(xj, 5LL)) +

                                                        2LL*Power(xj, 10LL)*(945LL + 945LL*r*xj + 420LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                             105LL*Power(r, 3LL)*Power(xj, 3LL) + 15LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                             Power(r, 5LL)*Power(xj, 5LL)) -

                                                        Power(xi, 10LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                         30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                         2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                                                        5LL*Power(xi, 8LL)*Power(xj, 2LL)*

                                                        (63LL + 105LL*r*xj + 84LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         42LL*Power(r, 3LL)*Power(xj, 3LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                         2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                                                        5LL*Power(xi, 6LL)*Power(xj, 4LL)*

                                                        (189LL + 315LL*r*xj + 252LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         132LL*Power(r, 3LL)*Power(xj, 3LL) + 36LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                         4LL*Power(r, 5LL)*Power(xj, 5LL)) +

                                                        5LL*Power(xi, 4LL)*Power(xj, 6LL)*

                                                        (315LL + 513LL*r*xj + 468LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         204LL*Power(r, 3LL)*Power(xj, 3LL) + 44LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                         4LL*Power(r, 5LL)*Power(xj, 5LL)))))/

                (45LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 7LL)*Power(xi + xj, 6LL)) -

                (90LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi, 2LL) - Power(xj, 2LL), 7LL) +

                 15LL*exp(2LL*r*xj)*Power(xj, 8LL)*

                 (-3LL*Power(xi, 7LL) - 7LL*Power(xi, 5LL)*Power(xj, 2LL) +

        7LL*Power(xi, 3LL)*Power(xj, 4LL) + 3LL*xi*Power(xj, 6LL)) +

                 30LL*exp(2LL*r*xj)*Power(xj, 9LL)*

                 (-15LL*Power(xi, 6LL) - 3LL*r*Power(xi, 7LL) - 63LL*Power(xi, 4LL)*Power(xj, 2LL) -

        7LL*r*Power(xi, 5LL)*Power(xj, 2LL) - 21LL*Power(xi, 2LL)*Power(xj, 4LL) +

        7LL*r*Power(xi, 3LL)*Power(xj, 4LL) + 3LL*Power(xj, 6LL) + 3LL*r*xi*Power(xj, 6LL)) +

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-10LL*Power(xi, 2LL)*Power(xj, 8LL)*

        (333LL*xj + 456LL*r*Power(xj, 2LL) + 225LL*Power(r, 2LL)*Power(xj, 3LL) +

            52LL*Power(r, 3LL)*Power(xj, 4LL) + 5LL*Power(r, 4LL)*Power(xj, 5LL)) +

        2LL*Power(xj, 10LL)*(945LL*xj + 840LL*r*Power(xj, 2LL) +

                             315LL*Power(r, 2LL)*Power(xj, 3LL) + 60LL*Power(r, 3LL)*Power(xj, 4LL) +

                             5LL*Power(r, 4LL)*Power(xj, 5LL)) -

        Power(xi, 10LL)*(75LL*xj + 120LL*r*Power(xj, 2LL) +

                         90LL*Power(r, 2LL)*Power(xj, 3LL) + 40LL*Power(r, 3LL)*Power(xj, 4LL) +

                         10LL*Power(r, 4LL)*Power(xj, 5LL)) +

        5LL*Power(xi, 8LL)*Power(xj, 2LL)*

        (105LL*xj + 168LL*r*Power(xj, 2LL) + 126LL*Power(r, 2LL)*Power(xj, 3LL) +

            56LL*Power(r, 3LL)*Power(xj, 4LL) + 10LL*Power(r, 4LL)*Power(xj, 5LL)) -

        5LL*Power(xi, 6LL)*Power(xj, 4LL)*

        (315LL*xj + 504LL*r*Power(xj, 2LL) + 396LL*Power(r, 2LL)*Power(xj, 3LL) +

            144LL*Power(r, 3LL)*Power(xj, 4LL) + 20LL*Power(r, 4LL)*Power(xj, 5LL)) +

        5LL*Power(xi, 4LL)*Power(xj, 6LL)*

        (513LL*xj + 936LL*r*Power(xj, 2LL) + 612LL*Power(r, 2LL)*Power(xj, 3LL) +

            176LL*Power(r, 3LL)*Power(xj, 4LL) + 20LL*Power(r, 4LL)*Power(xj, 5LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 5LL)*

                 (-10LL*Power(xi, 2LL)*Power(xj, 8LL)*

        (135LL + 333LL*r*xj + 228LL*Power(r, 2LL)*Power(xj, 2LL) +

            75LL*Power(r, 3LL)*Power(xj, 3LL) + 13LL*Power(r, 4LL)*Power(xj, 4LL) +

            Power(r, 5LL)*Power(xj, 5LL)) +

        2LL*Power(xj, 10LL)*(945LL + 945LL*r*xj + 420LL*Power(r, 2LL)*Power(xj, 2LL) +

                             105LL*Power(r, 3LL)*Power(xj, 3LL) + 15LL*Power(r, 4LL)*Power(xj, 4LL) +

                             Power(r, 5LL)*Power(xj, 5LL)) -

        Power(xi, 10LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                         30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                         2LL*Power(r, 5LL)*Power(xj, 5LL)) +

        5LL*Power(xi, 8LL)*Power(xj, 2LL)*

        (63LL + 105LL*r*xj + 84LL*Power(r, 2LL)*Power(xj, 2LL) +

            42LL*Power(r, 3LL)*Power(xj, 3LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL) +

            2LL*Power(r, 5LL)*Power(xj, 5LL)) -

        5LL*Power(xi, 6LL)*Power(xj, 4LL)*

        (189LL + 315LL*r*xj + 252LL*Power(r, 2LL)*Power(xj, 2LL) +

            132LL*Power(r, 3LL)*Power(xj, 3LL) + 36LL*Power(r, 4LL)*Power(xj, 4LL) +

            4LL*Power(r, 5LL)*Power(xj, 5LL)) +

        5LL*Power(xi, 4LL)*Power(xj, 6LL)*

        (315LL + 513LL*r*xj + 468LL*Power(r, 2LL)*Power(xj, 2LL) +

            204LL*Power(r, 3LL)*Power(xj, 3LL) + 44LL*Power(r, 4LL)*Power(xj, 4LL) +

            4LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (45LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 7LL)*Power(xi + xj, 7LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_3S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_1S_3S(r, xj, xi);
}
