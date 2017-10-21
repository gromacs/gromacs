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

#ifdef HAVE_LIBCLN
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

#else

double Slater_1S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (41*xi)/128

            ;
        }
        else
        {
            S = (1/r)*((-120960 + 120960*exp(2*rxi) - 203175*rxi - 164430*power(rxi, 2) -

                          84420*power(rxi, 3) - 30240*power(rxi, 4) - 7728*power(rxi, 5) -

                          1344*power(rxi, 6) - 128*power(rxi, 7))/(120960*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(power(xi, 6) + 7*power(xi, 5)*xj + 21*power(xi, 4)*power(xj, 2) +

                        35*power(xi, 3)*power(xj, 3) + 35*power(xi, 2)*power(xj, 4) +

                        21*xi*power(xj, 5) + 3*power(xj, 6)))/(3*power(xi + xj, 7))

            ;
        }
        else
        {
            S = (1/r)*((45*exp(2*(rxi + rxj))*power(power(rxi, 2) - power(rxj, 2), 7) +

                          15*exp(2*rxj)*power(rxj, 8)*

                          (-15*power(rxi, 6) - 3*power(rxi, 7) - 63*power(rxi, 4)*power(rxj, 2) -

                           7*power(rxi, 5)*power(rxj, 2) - 21*power(rxi, 2)*power(rxj, 4) +

                           7*power(rxi, 3)*power(rxj, 4) + 3*power(rxj, 6) + 3*rxi*power(rxj, 6)) +

                          exp(2*rxi)*power(rxi, 4)*

                          (-10*power(rxi, 2)*power(rxj, 8)*

                           (135 + 333*rxj + 228*power(rxj, 2) + 75*power(rxj, 3) +

                            13*power(rxj, 4) + power(rxj, 5)) +

                           2*power(rxj, 10)*(945 + 945*rxj + 420*power(rxj, 2) +

                                                 105*power(rxj, 3) + 15*power(rxj, 4) + power(rxj, 5)) -

                           power(rxi, 10)*(45 + 75*rxj + 60*power(rxj, 2) + 30*power(rxj, 3) +

                                             10*power(rxj, 4) + 2*power(rxj, 5)) +

                           5*power(rxi, 8)*power(rxj, 2)*

                           (63 + 105*rxj + 84*power(rxj, 2) + 42*power(rxj, 3) +

                            14*power(rxj, 4) + 2*power(rxj, 5)) -

                           5*power(rxi, 6)*power(rxj, 4)*

                           (189 + 315*rxj + 252*power(rxj, 2) + 132*power(rxj, 3) +

                            36*power(rxj, 4) + 4*power(rxj, 5)) +

                           5*power(rxi, 4)*power(rxj, 6)*

                           (315 + 513*rxj + 468*power(rxj, 2) + 204*power(rxj, 3) +

                            44*power(rxj, 4) + 4*power(rxj, 5))))/

                         (45*exp(2*(rxi + rxj))*power(rxi - rxj, 7)*power(rxi + rxj, 7))

                         );
        }

    }
    return S;
}


double Slater_3S_1S(double r, double xi, double xj)
{
    return Slater_1S_3S(r, xj, xi);
}

#endif
