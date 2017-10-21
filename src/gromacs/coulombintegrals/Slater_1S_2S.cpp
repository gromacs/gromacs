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

#else
double Slater_1S_2S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (7*xi)/16

            ;
        }
        else
        {
            S = (1/r)*((-240 + 240*exp(2*rxi) - 375*rxi - 270*power(rxi, 2) - 115*power(rxi, 3) -

                          30*power(rxi, 4) - 4*power(rxi, 5))/(240*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(power(xi, 4) + 5*power(xi, 3)*xj + 10*power(xi, 2)*power(xj, 2) +

                        10*xi*power(xj, 3) + 2*power(xj, 4)))/(2*power(xi + xj, 5))

            ;
        }
        else
        {
            S = (1/r)*((6*exp(2*(rxi + rxj))*power(power(rxi, 2) - power(rxj, 2), 5) +

                          6*exp(2*rxj)*power(rxj, 6)*

                          (-4*power(rxi, 4) - power(rxi, 5) - 5*power(rxi, 2)*power(rxj, 2) +

                           power(rxj, 4) + rxi*power(rxj, 4)) -

                          exp(2*rxi)*power(rxi, 4)*

                          (power(rxi, 6)*(6 + 9*rxj + 6*power(rxj, 2) + 2*power(rxj, 3)) -

                           3*power(rxi, 4)*power(rxj, 2)*

                           (10 + 15*rxj + 10*power(rxj, 2) + 2*power(rxj, 3)) +

                           3*power(rxi, 2)*power(rxj, 4)*

                           (20 + 33*rxj + 14*power(rxj, 2) + 2*power(rxj, 3)) -

                           power(rxj, 6)*(84 + 63*rxj + 18*power(rxj, 2) + 2*power(rxj, 3))))/

                         (6*exp(2*(rxi + rxj))*power(rxi - rxj, 5)*power(rxi + rxj, 5))

                         );
        }

    }
    return S;
}


double Slater_2S_1S(double r, double xi, double xj)
{
    return Slater_1S_2S(r, xj, xi);
}
#endif
