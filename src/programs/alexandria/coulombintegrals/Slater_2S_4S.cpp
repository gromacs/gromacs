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

cl_R Slater_2S_4S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (975LL*xi)/4096LL

            ;
        }
        else
        {
            S = (1LL/r)*((-638668800LL + 638668800LL*exp(2LL*rxi) - 1125310725LL*rxi -

                          973283850LL*Power(rxi, 2LL) - 549063900LL*Power(rxi, 3LL) -

                          226195200LL*Power(rxi, 4LL) - 72099720LL*Power(rxi, 5LL) - 18350640LL*Power(rxi, 6LL) -

                          3785760LL*Power(rxi, 7LL) - 633600LL*Power(rxi, 8LL) - 84480LL*Power(rxi, 9LL) -

                          8448LL*Power(rxi, 10LL) - 512LL*Power(rxi, 11LL))/(6.386688e8*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 10LL) + 11LL*Power(xi, 9LL)*xj + 55LL*Power(xi, 8LL)*Power(xj, 2LL) +

                        165LL*Power(xi, 7LL)*Power(xj, 3LL) + 330LL*Power(xi, 6LL)*Power(xj, 4LL) +

                        462LL*Power(xi, 5LL)*Power(xj, 5LL) + 462LL*Power(xi, 4LL)*Power(xj, 6LL) +

                        330LL*Power(xi, 3LL)*Power(xj, 7LL) + 110LL*Power(xi, 2LL)*Power(xj, 8LL) +

                        22LL*xi*Power(xj, 9LL) + 2LL*Power(xj, 10LL)))/(4LL*Power(xi + xj, 11LL))

            ;
        }
        else
        {
            S = (1LL/r)*((1260LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 11LL) +

                          210LL*exp(2LL*rxj)*Power(rxj, 10LL)*

                          (-36LL*Power(rxi, 14LL) - 2LL*Power(rxi, 15LL) -

                           1287LL*Power(rxi, 9LL)*Power(rxj, 4LL) + 6LL*Power(rxj, 12LL) +

                           9LL*rxi*Power(rxj, 12LL) - 22LL*Power(rxi, 7LL)*Power(rxj, 6LL)*

                           (-135LL + Power(rxj, 2LL)) +

                           6LL*Power(rxi, 2LL)*Power(rxj, 10LL)*(-11LL + Power(rxj, 2LL)) -

                           66LL*Power(rxi, 4LL)*Power(rxj, 8LL)*(-5LL + Power(rxj, 2LL)) +

                           8LL*Power(rxi, 5LL)*Power(rxj, 8LL)*(99LL + Power(rxj, 2LL)) +

                           Power(rxi, 3LL)*Power(rxj, 10LL)*(-99LL + 2LL*Power(rxj, 2LL)) -

                           132LL*Power(rxi, 6LL)*Power(rxj, 6LL)*(27LL + 2LL*Power(rxj, 2LL)) -

                           78LL*Power(rxi, 12LL)*(7LL + 3LL*Power(rxj, 2LL)) -

                           2LL*Power(rxi, 13LL)*(117LL + 4LL*Power(rxj, 2LL)) +

                           66LL*Power(rxi, 8LL)*Power(rxj, 4LL)*(-191LL + 6LL*Power(rxj, 2LL)) +

                           Power(rxi, 11LL)*Power(rxj, 2LL)*(-2151LL + 22LL*Power(rxj, 2LL)) +

                           6LL*Power(rxi, 10LL)*Power(rxj, 2LL)*(-1099LL + 33LL*Power(rxj, 2LL))) -

                          exp(2LL*rxi)*Power(rxi, 6LL)*

                          (385LL*Power(rxi, 8LL)*Power(rxj, 8LL)*

                           (1080LL + 1935LL*rxj + 1350LL*Power(rxj, 2LL) + 1170LL*Power(rxj, 3LL) +

           420LL*Power(rxj, 4LL) + 66LL*Power(rxj, 5LL) + 4LL*Power(rxj, 6LL)) -

                           4LL*Power(rxj, 16LL)*(135135LL + 135135LL*rxj + 62370LL*Power(rxj, 2LL) +

                                                 17325LL*Power(rxj, 3LL) + 3150LL*Power(rxj, 4LL) + 378LL*Power(rxj, 5LL) +

                                                 28LL*Power(rxj, 6LL) + Power(rxj, 7LL)) +

                           Power(rxi, 16LL)*(1260LL + 2205LL*rxj + 1890LL*Power(rxj, 2LL) +

                                             1050LL*Power(rxj, 3LL) + 420LL*Power(rxj, 4LL) + 126LL*Power(rxj, 5LL) +

                                             28LL*Power(rxj, 6LL) + 4LL*Power(rxj, 7LL)) +

                           7LL*Power(rxi, 6LL)*Power(rxj, 10LL)*

                           (-99540LL - 58095LL*rxj - 190710LL*Power(rxj, 2LL) - 100950LL*Power(rxj, 3LL) -

           21660LL*Power(rxj, 4LL) - 1938LL*Power(rxj, 5LL) - 4LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) - 7LL*Power(rxi, 4LL)*Power(rxj, 12LL)*

                           (114660LL - 343395LL*rxj - 242910LL*Power(rxj, 2LL) - 61950LL*Power(rxj, 3LL) -

           6060LL*Power(rxj, 4LL) + 282LL*Power(rxj, 5LL) + 116LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) + 7LL*Power(rxi, 12LL)*Power(rxj, 4LL)*

                           (9900LL + 17325LL*rxj + 14850LL*Power(rxj, 2LL) + 8250LL*Power(rxj, 3LL) +

           3300LL*Power(rxj, 4LL) + 1074LL*Power(rxj, 5LL) + 164LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) - 7LL*Power(rxi, 10LL)*Power(rxj, 6LL)*

                           (29700LL + 51975LL*rxj + 44550LL*Power(rxj, 2LL) + 23850LL*Power(rxj, 3LL) +

           11700LL*Power(rxj, 4LL) + 2814LL*Power(rxj, 5LL) + 284LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) - Power(rxi, 14LL)*Power(rxj, 2LL)*

                           (13860LL + 24255LL*rxj + 20790LL*Power(rxj, 2LL) + 11550LL*Power(rxj, 3LL) +

           4620LL*Power(rxj, 4LL) + 1386LL*Power(rxj, 5LL) + 308LL*Power(rxj, 6LL) +

           24LL*Power(rxj, 7LL)) + Power(rxi, 2LL)*Power(rxj, 14LL)*

                           (-3063060LL - 1936935LL*rxj - 408870LL*Power(rxj, 2LL) + 11550LL*Power(rxj, 3LL) +

           23100LL*Power(rxj, 4LL) + 5082LL*Power(rxj, 5LL) + 532LL*Power(rxj, 6LL) +

           24LL*Power(rxj, 7LL))))/

                         (1260LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 11LL)*Power(rxi + rxj, 11LL))

                         );
        }

    }
    return S;
}


cl_R Slater_4S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_2S_4S(r, xj, xi);
}
