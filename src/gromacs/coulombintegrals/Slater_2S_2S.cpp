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
            S = (1LL/r)*((-80640LL + 80640LL*exp(2LL*rxi) - 131985LL*rxi - 102690LL*Pow(rxi, 2LL) -

                          49980LL*Pow(rxi, 3LL) - 16800LL*Pow(rxi, 4LL) - 4032LL*Pow(rxi, 5LL) -

                          672LL*Pow(rxi, 6LL) - 64LL*Pow(rxi, 7LL))/(80640LL*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Pow(xi, 6LL) + 7LL*Pow(xi, 5LL)*xj + 21LL*Pow(xi, 4LL)*Pow(xj, 2LL) +

                        35LL*Pow(xi, 3LL)*Pow(xj, 3LL) + 21LL*Pow(xi, 2LL)*Pow(xj, 4LL) +

                        7LL*xi*Pow(xj, 5LL) + Pow(xj, 6LL)))/(2LL*Pow(xi + xj, 7LL))

            ;
        }
        else
        {
            S = (1LL/r)*((6LL*exp(2LL*(rxi + rxj))*Pow(Pow(rxi, 2LL) - Pow(rxj, 2LL), 7LL) -

                          exp(2LL*rxi)*Pow(rxi, 6LL)*

                          (21LL*Pow(rxi, 4LL)*Pow(rxj, 4LL)*(6LL + 11LL*rxj + 2LL*Pow(rxj, 2LL)) -

                           2LL*Pow(rxj, 8LL)*(90LL + 54LL*rxj + 12LL*Pow(rxj, 2LL) + Pow(rxj, 3LL)) +

                           Pow(rxi, 8LL)*(6LL + 9LL*rxj + 6LL*Pow(rxj, 2LL) + 2LL*Pow(rxj, 3LL)) +

                           Pow(rxi, 2LL)*Pow(rxj, 6LL)*

                           (-390LL - 69LL*rxj + 18LL*Pow(rxj, 2LL) + 4LL*Pow(rxj, 3LL)) -

                           Pow(rxi, 6LL)*Pow(rxj, 2LL)*

                           (42LL + 63LL*rxj + 42LL*Pow(rxj, 2LL) + 4LL*Pow(rxj, 3LL))) +

                          exp(2LL*rxj)*Pow(rxj, 6LL)*

                          (-24LL*Pow(rxi, 10LL) - 2LL*Pow(rxi, 11LL) - 69LL*Pow(rxi, 7LL)*Pow(rxj, 2LL) +

                           6LL*Pow(rxj, 8LL) + 9LL*rxi*Pow(rxj, 8LL) +

                           4LL*Pow(rxi, 9LL)*(-27LL + Pow(rxj, 2LL)) +

                           18LL*Pow(rxi, 8LL)*(-10LL + Pow(rxj, 2LL)) +

                           6LL*Pow(rxi, 2LL)*Pow(rxj, 6LL)*(-7LL + Pow(rxj, 2LL)) -

                           42LL*Pow(rxi, 4LL)*Pow(rxj, 4LL)*(-3LL + Pow(rxj, 2LL)) +

                           Pow(rxi, 3LL)*Pow(rxj, 6LL)*(-63LL + 2LL*Pow(rxj, 2LL)) +

                           6LL*Pow(rxi, 6LL)*Pow(rxj, 2LL)*(-65LL + 7LL*Pow(rxj, 2LL)) +

                           Pow(rxi, 5LL)*(231LL*Pow(rxj, 4LL) - 4LL*Pow(rxj, 6LL))))/

                         (6LL*exp(2LL*(rxi + rxj))*Pow(rxi - rxj, 7LL)*Pow(rxi + rxj, 7LL))

                         );
        }

    }
    return S;
}

#else

double Slater_2S_2S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (93*xi)/256

            ;
        }
        else
        {
            S = (1/r)*((-80640 + 80640*exp(2*rxi) - 131985*rxi - 102690*pow(rxi, 2) -

                          49980*pow(rxi, 3) - 16800*pow(rxi, 4) - 4032*pow(rxi, 5) -

                          672*pow(rxi, 6) - 64*pow(rxi, 7))/(80640*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 6) + 7*pow(xi, 5)*xj + 21*pow(xi, 4)*pow(xj, 2) +

                        35*pow(xi, 3)*pow(xj, 3) + 21*pow(xi, 2)*pow(xj, 4) +

                        7*xi*pow(xj, 5) + pow(xj, 6)))/(2*pow(xi + xj, 7))

            ;
        }
        else
        {
            S = (1/r)*((6*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 7) -

                          exp(2*rxi)*pow(rxi, 6)*

                          (21*pow(rxi, 4)*pow(rxj, 4)*(6 + 11*rxj + 2*pow(rxj, 2)) -

                           2*pow(rxj, 8)*(90 + 54*rxj + 12*pow(rxj, 2) + pow(rxj, 3)) +

                           pow(rxi, 8)*(6 + 9*rxj + 6*pow(rxj, 2) + 2*pow(rxj, 3)) +

                           pow(rxi, 2)*pow(rxj, 6)*

                           (-390 - 69*rxj + 18*pow(rxj, 2) + 4*pow(rxj, 3)) -

                           pow(rxi, 6)*pow(rxj, 2)*

                           (42 + 63*rxj + 42*pow(rxj, 2) + 4*pow(rxj, 3))) +

                          exp(2*rxj)*pow(rxj, 6)*

                          (-24*pow(rxi, 10) - 2*pow(rxi, 11) - 69*pow(rxi, 7)*pow(rxj, 2) +

                           6*pow(rxj, 8) + 9*rxi*pow(rxj, 8) +

                           4*pow(rxi, 9)*(-27 + pow(rxj, 2)) +

                           18*pow(rxi, 8)*(-10 + pow(rxj, 2)) +

                           6*pow(rxi, 2)*pow(rxj, 6)*(-7 + pow(rxj, 2)) -

                           42*pow(rxi, 4)*pow(rxj, 4)*(-3 + pow(rxj, 2)) +

                           pow(rxi, 3)*pow(rxj, 6)*(-63 + 2*pow(rxj, 2)) +

                           6*pow(rxi, 6)*pow(rxj, 2)*(-65 + 7*pow(rxj, 2)) +

                           pow(rxi, 5)*(231*pow(rxj, 4) - 4*pow(rxj, 6))))/

                         (6*exp(2*(rxi + rxj))*pow(rxi - rxj, 7)*pow(rxi + rxj, 7))

                         );
        }

    }
    return S;
}

#endif
