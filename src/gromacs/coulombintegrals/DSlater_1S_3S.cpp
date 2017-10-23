/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#if HAVE_LIBCLN
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
            S = -(-203175LL*xi + 241920LL*exp(2LL*r*xi)*xi - 328860LL*r*Pow(xi, 2LL) -

                  253260LL*Pow(r, 2LL)*Pow(xi, 3LL) - 120960LL*Pow(r, 3LL)*Pow(xi, 4LL) -

                  38640LL*Pow(r, 4LL)*Pow(xi, 5LL) - 8064LL*Pow(r, 5LL)*Pow(xi, 6LL) -

                  896LL*Pow(r, 6LL)*Pow(xi, 7LL))/(120960LL*exp(2LL*r*xi)*r) +

                (-120960LL + 120960LL*exp(2LL*r*xi) - 203175LL*r*xi -

                 164430LL*Pow(r, 2LL)*Pow(xi, 2LL) - 84420LL*Pow(r, 3LL)*Pow(xi, 3LL) -

                 30240LL*Pow(r, 4LL)*Pow(xi, 4LL) - 7728LL*Pow(r, 5LL)*Pow(xi, 5LL) -

                 1344LL*Pow(r, 6LL)*Pow(xi, 6LL) - 128LL*Pow(r, 7LL)*Pow(xi, 7LL))/

                (120960LL*exp(2LL*r*xi)*Pow(r, 2LL)) +

                (xi*(-120960LL + 120960LL*exp(2LL*r*xi) - 203175LL*r*xi -

                     164430LL*Pow(r, 2LL)*Pow(xi, 2LL) - 84420LL*Pow(r, 3LL)*Pow(xi, 3LL) -

                     30240LL*Pow(r, 4LL)*Pow(xi, 4LL) - 7728LL*Pow(r, 5LL)*Pow(xi, 5LL) -

                     1344LL*Pow(r, 6LL)*Pow(xi, 6LL) - 128LL*Pow(r, 7LL)*Pow(xi, 7LL)))/

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
            S = (45LL*exp(2LL*r*(xi + xj))*Pow(Pow(xi, 2LL) - Pow(xj, 2LL), 7LL) +

                 15LL*exp(2LL*r*xj)*Pow(xj, 8LL)*

                 (-15LL*Pow(xi, 6LL) - 3LL*r*Pow(xi, 7LL) - 63LL*Pow(xi, 4LL)*Pow(xj, 2LL) -

                  7LL*r*Pow(xi, 5LL)*Pow(xj, 2LL) - 21LL*Pow(xi, 2LL)*Pow(xj, 4LL) +

                  7LL*r*Pow(xi, 3LL)*Pow(xj, 4LL) + 3LL*Pow(xj, 6LL) + 3LL*r*xi*Pow(xj, 6LL)) +

                 exp(2LL*r*xi)*Pow(xi, 4LL)*

                 (-10LL*Pow(xi, 2LL)*Pow(xj, 8LL)*

                  (135LL + 333LL*r*xj + 228LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   75LL*Pow(r, 3LL)*Pow(xj, 3LL) + 13LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   Pow(r, 5LL)*Pow(xj, 5LL)) +

                  2LL*Pow(xj, 10LL)*(945LL + 945LL*r*xj + 420LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                       105LL*Pow(r, 3LL)*Pow(xj, 3LL) + 15LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                       Pow(r, 5LL)*Pow(xj, 5LL)) -

                  Pow(xi, 10LL)*(45LL + 75LL*r*xj + 60LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                   30LL*Pow(r, 3LL)*Pow(xj, 3LL) + 10LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                   2LL*Pow(r, 5LL)*Pow(xj, 5LL)) +

                  5LL*Pow(xi, 8LL)*Pow(xj, 2LL)*

                  (63LL + 105LL*r*xj + 84LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   42LL*Pow(r, 3LL)*Pow(xj, 3LL) + 14LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   2LL*Pow(r, 5LL)*Pow(xj, 5LL)) -

                  5LL*Pow(xi, 6LL)*Pow(xj, 4LL)*

                  (189LL + 315LL*r*xj + 252LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   132LL*Pow(r, 3LL)*Pow(xj, 3LL) + 36LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   4LL*Pow(r, 5LL)*Pow(xj, 5LL)) +

                  5LL*Pow(xi, 4LL)*Pow(xj, 6LL)*

                  (315LL + 513LL*r*xj + 468LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   204LL*Pow(r, 3LL)*Pow(xj, 3LL) + 44LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   4LL*Pow(r, 5LL)*Pow(xj, 5LL))))/

                (45LL*exp(2LL*r*(xi + xj))*Pow(r, 2LL)*Pow(xi - xj, 7LL)*Pow(xi + xj, 7LL))

                + (2LL*(45LL*exp(2LL*r*(xi + xj))*Pow(Pow(xi, 2LL) - Pow(xj, 2LL), 7LL) +

                        15LL*exp(2LL*r*xj)*Pow(xj, 8LL)*

                        (-15LL*Pow(xi, 6LL) - 3LL*r*Pow(xi, 7LL) - 63LL*Pow(xi, 4LL)*Pow(xj, 2LL) -

                         7LL*r*Pow(xi, 5LL)*Pow(xj, 2LL) - 21LL*Pow(xi, 2LL)*Pow(xj, 4LL) +

                         7LL*r*Pow(xi, 3LL)*Pow(xj, 4LL) + 3LL*Pow(xj, 6LL) + 3LL*r*xi*Pow(xj, 6LL))

                        + exp(2LL*r*xi)*Pow(xi, 4LL)*(-10LL*Pow(xi, 2LL)*Pow(xj, 8LL)*

                                                        (135LL + 333LL*r*xj + 228LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                                         75LL*Pow(r, 3LL)*Pow(xj, 3LL) + 13LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                                         Pow(r, 5LL)*Pow(xj, 5LL)) +

                                                        2LL*Pow(xj, 10LL)*(945LL + 945LL*r*xj + 420LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                                                             105LL*Pow(r, 3LL)*Pow(xj, 3LL) + 15LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                                                             Pow(r, 5LL)*Pow(xj, 5LL)) -

                                                        Pow(xi, 10LL)*(45LL + 75LL*r*xj + 60LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                                                         30LL*Pow(r, 3LL)*Pow(xj, 3LL) + 10LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                                                         2LL*Pow(r, 5LL)*Pow(xj, 5LL)) +

                                                        5LL*Pow(xi, 8LL)*Pow(xj, 2LL)*

                                                        (63LL + 105LL*r*xj + 84LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                                         42LL*Pow(r, 3LL)*Pow(xj, 3LL) + 14LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                                         2LL*Pow(r, 5LL)*Pow(xj, 5LL)) -

                                                        5LL*Pow(xi, 6LL)*Pow(xj, 4LL)*

                                                        (189LL + 315LL*r*xj + 252LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                                         132LL*Pow(r, 3LL)*Pow(xj, 3LL) + 36LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                                         4LL*Pow(r, 5LL)*Pow(xj, 5LL)) +

                                                        5LL*Pow(xi, 4LL)*Pow(xj, 6LL)*

                                                        (315LL + 513LL*r*xj + 468LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                                         204LL*Pow(r, 3LL)*Pow(xj, 3LL) + 44LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                                         4LL*Pow(r, 5LL)*Pow(xj, 5LL)))))/

                (45LL*exp(2LL*r*(xi + xj))*r*Pow(xi - xj, 7LL)*Pow(xi + xj, 6LL)) -

                (90LL*exp(2LL*r*(xi + xj))*(xi + xj)*Pow(Pow(xi, 2LL) - Pow(xj, 2LL), 7LL) +

                 15LL*exp(2LL*r*xj)*Pow(xj, 8LL)*

                 (-3LL*Pow(xi, 7LL) - 7LL*Pow(xi, 5LL)*Pow(xj, 2LL) +

                  7LL*Pow(xi, 3LL)*Pow(xj, 4LL) + 3LL*xi*Pow(xj, 6LL)) +

                 30LL*exp(2LL*r*xj)*Pow(xj, 9LL)*

                 (-15LL*Pow(xi, 6LL) - 3LL*r*Pow(xi, 7LL) - 63LL*Pow(xi, 4LL)*Pow(xj, 2LL) -

                  7LL*r*Pow(xi, 5LL)*Pow(xj, 2LL) - 21LL*Pow(xi, 2LL)*Pow(xj, 4LL) +

                  7LL*r*Pow(xi, 3LL)*Pow(xj, 4LL) + 3LL*Pow(xj, 6LL) + 3LL*r*xi*Pow(xj, 6LL)) +

                 exp(2LL*r*xi)*Pow(xi, 4LL)*

                 (-10LL*Pow(xi, 2LL)*Pow(xj, 8LL)*

                  (333LL*xj + 456LL*r*Pow(xj, 2LL) + 225LL*Pow(r, 2LL)*Pow(xj, 3LL) +

                   52LL*Pow(r, 3LL)*Pow(xj, 4LL) + 5LL*Pow(r, 4LL)*Pow(xj, 5LL)) +

                  2LL*Pow(xj, 10LL)*(945LL*xj + 840LL*r*Pow(xj, 2LL) +

                                       315LL*Pow(r, 2LL)*Pow(xj, 3LL) + 60LL*Pow(r, 3LL)*Pow(xj, 4LL) +

                                       5LL*Pow(r, 4LL)*Pow(xj, 5LL)) -

                  Pow(xi, 10LL)*(75LL*xj + 120LL*r*Pow(xj, 2LL) +

                                   90LL*Pow(r, 2LL)*Pow(xj, 3LL) + 40LL*Pow(r, 3LL)*Pow(xj, 4LL) +

                                   10LL*Pow(r, 4LL)*Pow(xj, 5LL)) +

                  5LL*Pow(xi, 8LL)*Pow(xj, 2LL)*

                  (105LL*xj + 168LL*r*Pow(xj, 2LL) + 126LL*Pow(r, 2LL)*Pow(xj, 3LL) +

                   56LL*Pow(r, 3LL)*Pow(xj, 4LL) + 10LL*Pow(r, 4LL)*Pow(xj, 5LL)) -

                  5LL*Pow(xi, 6LL)*Pow(xj, 4LL)*

                  (315LL*xj + 504LL*r*Pow(xj, 2LL) + 396LL*Pow(r, 2LL)*Pow(xj, 3LL) +

                   144LL*Pow(r, 3LL)*Pow(xj, 4LL) + 20LL*Pow(r, 4LL)*Pow(xj, 5LL)) +

                  5LL*Pow(xi, 4LL)*Pow(xj, 6LL)*

                  (513LL*xj + 936LL*r*Pow(xj, 2LL) + 612LL*Pow(r, 2LL)*Pow(xj, 3LL) +

                   176LL*Pow(r, 3LL)*Pow(xj, 4LL) + 20LL*Pow(r, 4LL)*Pow(xj, 5LL))) +

                 2LL*exp(2LL*r*xi)*Pow(xi, 5LL)*

                 (-10LL*Pow(xi, 2LL)*Pow(xj, 8LL)*

                  (135LL + 333LL*r*xj + 228LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   75LL*Pow(r, 3LL)*Pow(xj, 3LL) + 13LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   Pow(r, 5LL)*Pow(xj, 5LL)) +

                  2LL*Pow(xj, 10LL)*(945LL + 945LL*r*xj + 420LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                       105LL*Pow(r, 3LL)*Pow(xj, 3LL) + 15LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                       Pow(r, 5LL)*Pow(xj, 5LL)) -

                  Pow(xi, 10LL)*(45LL + 75LL*r*xj + 60LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                   30LL*Pow(r, 3LL)*Pow(xj, 3LL) + 10LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                                   2LL*Pow(r, 5LL)*Pow(xj, 5LL)) +

                  5LL*Pow(xi, 8LL)*Pow(xj, 2LL)*

                  (63LL + 105LL*r*xj + 84LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   42LL*Pow(r, 3LL)*Pow(xj, 3LL) + 14LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   2LL*Pow(r, 5LL)*Pow(xj, 5LL)) -

                  5LL*Pow(xi, 6LL)*Pow(xj, 4LL)*

                  (189LL + 315LL*r*xj + 252LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   132LL*Pow(r, 3LL)*Pow(xj, 3LL) + 36LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   4LL*Pow(r, 5LL)*Pow(xj, 5LL)) +

                  5LL*Pow(xi, 4LL)*Pow(xj, 6LL)*

                  (315LL + 513LL*r*xj + 468LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   204LL*Pow(r, 3LL)*Pow(xj, 3LL) + 44LL*Pow(r, 4LL)*Pow(xj, 4LL) +

                   4LL*Pow(r, 5LL)*Pow(xj, 5LL))))/

                (45LL*exp(2LL*r*(xi + xj))*r*Pow(xi - xj, 7LL)*Pow(xi + xj, 7LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_3S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_1S_3S(r, xj, xi);
}

#else

double DSlater_1S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = -(-203175*xi + 241920*exp(2*r*xi)*xi - 328860*r*pow(xi, 2) -

                  253260*pow(r, 2)*pow(xi, 3) - 120960*pow(r, 3)*pow(xi, 4) -

                  38640*pow(r, 4)*pow(xi, 5) - 8064*pow(r, 5)*pow(xi, 6) -

                  896*pow(r, 6)*pow(xi, 7))/(120960*exp(2*r*xi)*r) +

                (-120960 + 120960*exp(2*r*xi) - 203175*r*xi -

                 164430*pow(r, 2)*pow(xi, 2) - 84420*pow(r, 3)*pow(xi, 3) -

                 30240*pow(r, 4)*pow(xi, 4) - 7728*pow(r, 5)*pow(xi, 5) -

                 1344*pow(r, 6)*pow(xi, 6) - 128*pow(r, 7)*pow(xi, 7))/

                (120960*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-120960 + 120960*exp(2*r*xi) - 203175*r*xi -

                     164430*pow(r, 2)*pow(xi, 2) - 84420*pow(r, 3)*pow(xi, 3) -

                     30240*pow(r, 4)*pow(xi, 4) - 7728*pow(r, 5)*pow(xi, 5) -

                     1344*pow(r, 6)*pow(xi, 6) - 128*pow(r, 7)*pow(xi, 7)))/

                (60480*exp(2*r*xi)*r)

            ;
        }

    }
    else
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = (45*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) +

                 15*exp(2*r*xj)*pow(xj, 8)*

                 (-15*pow(xi, 6) - 3*r*pow(xi, 7) - 63*pow(xi, 4)*pow(xj, 2) -

                  7*r*pow(xi, 5)*pow(xj, 2) - 21*pow(xi, 2)*pow(xj, 4) +

                  7*r*pow(xi, 3)*pow(xj, 4) + 3*pow(xj, 6) + 3*r*xi*pow(xj, 6)) +

                 exp(2*r*xi)*pow(xi, 4)*

                 (-10*pow(xi, 2)*pow(xj, 8)*

                  (135 + 333*r*xj + 228*pow(r, 2)*pow(xj, 2) +

                   75*pow(r, 3)*pow(xj, 3) + 13*pow(r, 4)*pow(xj, 4) +

                   pow(r, 5)*pow(xj, 5)) +

                  2*pow(xj, 10)*(945 + 945*r*xj + 420*pow(r, 2)*pow(xj, 2) +

                                       105*pow(r, 3)*pow(xj, 3) + 15*pow(r, 4)*pow(xj, 4) +

                                       pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 10)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                   30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                   2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 8)*pow(xj, 2)*

                  (63 + 105*r*xj + 84*pow(r, 2)*pow(xj, 2) +

                   42*pow(r, 3)*pow(xj, 3) + 14*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  5*pow(xi, 6)*pow(xj, 4)*

                  (189 + 315*r*xj + 252*pow(r, 2)*pow(xj, 2) +

                   132*pow(r, 3)*pow(xj, 3) + 36*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 4)*pow(xj, 6)*

                  (315 + 513*r*xj + 468*pow(r, 2)*pow(xj, 2) +

                   204*pow(r, 3)*pow(xj, 3) + 44*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5))))/

                (45*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 7)*pow(xi + xj, 7))

                + (2*(45*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) +

                        15*exp(2*r*xj)*pow(xj, 8)*

                        (-15*pow(xi, 6) - 3*r*pow(xi, 7) - 63*pow(xi, 4)*pow(xj, 2) -

                         7*r*pow(xi, 5)*pow(xj, 2) - 21*pow(xi, 2)*pow(xj, 4) +

                         7*r*pow(xi, 3)*pow(xj, 4) + 3*pow(xj, 6) + 3*r*xi*pow(xj, 6))

                        + exp(2*r*xi)*pow(xi, 4)*(-10*pow(xi, 2)*pow(xj, 8)*

                                                        (135 + 333*r*xj + 228*pow(r, 2)*pow(xj, 2) +

                                                         75*pow(r, 3)*pow(xj, 3) + 13*pow(r, 4)*pow(xj, 4) +

                                                         pow(r, 5)*pow(xj, 5)) +

                                                        2*pow(xj, 10)*(945 + 945*r*xj + 420*pow(r, 2)*pow(xj, 2) +

                                                                             105*pow(r, 3)*pow(xj, 3) + 15*pow(r, 4)*pow(xj, 4) +

                                                                             pow(r, 5)*pow(xj, 5)) -

                                                        pow(xi, 10)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                                                         30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                                                         2*pow(r, 5)*pow(xj, 5)) +

                                                        5*pow(xi, 8)*pow(xj, 2)*

                                                        (63 + 105*r*xj + 84*pow(r, 2)*pow(xj, 2) +

                                                         42*pow(r, 3)*pow(xj, 3) + 14*pow(r, 4)*pow(xj, 4) +

                                                         2*pow(r, 5)*pow(xj, 5)) -

                                                        5*pow(xi, 6)*pow(xj, 4)*

                                                        (189 + 315*r*xj + 252*pow(r, 2)*pow(xj, 2) +

                                                         132*pow(r, 3)*pow(xj, 3) + 36*pow(r, 4)*pow(xj, 4) +

                                                         4*pow(r, 5)*pow(xj, 5)) +

                                                        5*pow(xi, 4)*pow(xj, 6)*

                                                        (315 + 513*r*xj + 468*pow(r, 2)*pow(xj, 2) +

                                                         204*pow(r, 3)*pow(xj, 3) + 44*pow(r, 4)*pow(xj, 4) +

                                                         4*pow(r, 5)*pow(xj, 5)))))/

                (45*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 6)) -

                (90*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 7) +

                 15*exp(2*r*xj)*pow(xj, 8)*

                 (-3*pow(xi, 7) - 7*pow(xi, 5)*pow(xj, 2) +

                  7*pow(xi, 3)*pow(xj, 4) + 3*xi*pow(xj, 6)) +

                 30*exp(2*r*xj)*pow(xj, 9)*

                 (-15*pow(xi, 6) - 3*r*pow(xi, 7) - 63*pow(xi, 4)*pow(xj, 2) -

                  7*r*pow(xi, 5)*pow(xj, 2) - 21*pow(xi, 2)*pow(xj, 4) +

                  7*r*pow(xi, 3)*pow(xj, 4) + 3*pow(xj, 6) + 3*r*xi*pow(xj, 6)) +

                 exp(2*r*xi)*pow(xi, 4)*

                 (-10*pow(xi, 2)*pow(xj, 8)*

                  (333*xj + 456*r*pow(xj, 2) + 225*pow(r, 2)*pow(xj, 3) +

                   52*pow(r, 3)*pow(xj, 4) + 5*pow(r, 4)*pow(xj, 5)) +

                  2*pow(xj, 10)*(945*xj + 840*r*pow(xj, 2) +

                                       315*pow(r, 2)*pow(xj, 3) + 60*pow(r, 3)*pow(xj, 4) +

                                       5*pow(r, 4)*pow(xj, 5)) -

                  pow(xi, 10)*(75*xj + 120*r*pow(xj, 2) +

                                   90*pow(r, 2)*pow(xj, 3) + 40*pow(r, 3)*pow(xj, 4) +

                                   10*pow(r, 4)*pow(xj, 5)) +

                  5*pow(xi, 8)*pow(xj, 2)*

                  (105*xj + 168*r*pow(xj, 2) + 126*pow(r, 2)*pow(xj, 3) +

                   56*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) -

                  5*pow(xi, 6)*pow(xj, 4)*

                  (315*xj + 504*r*pow(xj, 2) + 396*pow(r, 2)*pow(xj, 3) +

                   144*pow(r, 3)*pow(xj, 4) + 20*pow(r, 4)*pow(xj, 5)) +

                  5*pow(xi, 4)*pow(xj, 6)*

                  (513*xj + 936*r*pow(xj, 2) + 612*pow(r, 2)*pow(xj, 3) +

                   176*pow(r, 3)*pow(xj, 4) + 20*pow(r, 4)*pow(xj, 5))) +

                 2*exp(2*r*xi)*pow(xi, 5)*

                 (-10*pow(xi, 2)*pow(xj, 8)*

                  (135 + 333*r*xj + 228*pow(r, 2)*pow(xj, 2) +

                   75*pow(r, 3)*pow(xj, 3) + 13*pow(r, 4)*pow(xj, 4) +

                   pow(r, 5)*pow(xj, 5)) +

                  2*pow(xj, 10)*(945 + 945*r*xj + 420*pow(r, 2)*pow(xj, 2) +

                                       105*pow(r, 3)*pow(xj, 3) + 15*pow(r, 4)*pow(xj, 4) +

                                       pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 10)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                   30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                   2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 8)*pow(xj, 2)*

                  (63 + 105*r*xj + 84*pow(r, 2)*pow(xj, 2) +

                   42*pow(r, 3)*pow(xj, 3) + 14*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  5*pow(xi, 6)*pow(xj, 4)*

                  (189 + 315*r*xj + 252*pow(r, 2)*pow(xj, 2) +

                   132*pow(r, 3)*pow(xj, 3) + 36*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 4)*pow(xj, 6)*

                  (315 + 513*r*xj + 468*pow(r, 2)*pow(xj, 2) +

                   204*pow(r, 3)*pow(xj, 3) + 44*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5))))/

                (45*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 7))

            ;
        }

    }
    return S;
}


double DSlater_3S_1S(double r, double xi, double xj)
{
    return DSlater_1S_3S(r, xj, xi);
}

#endif
