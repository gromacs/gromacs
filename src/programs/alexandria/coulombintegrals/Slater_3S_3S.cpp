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

cl_R Slater_3S_3S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (793LL*xi)/3072LL

            ;
        }
        else
        {
            S = (1LL/r)*((-1437004800LL + 1437004800LL*exp(2LL*rxi) - 2503064025LL*rxi -

                          2132118450LL*Power(rxi, 2LL) - 1180664100LL*Power(rxi, 3LL) -

                          476506800LL*Power(rxi, 4LL) - 148856400LL*Power(rxi, 5LL) -

                          37255680LL*Power(rxi, 6LL) - 7603200LL*Power(rxi, 7LL) - 1267200LL*Power(rxi, 8LL) -

                          168960LL*Power(rxi, 9LL) - 16896LL*Power(rxi, 10LL) - 1024LL*Power(rxi, 11LL))/

                         (1.4370048e9*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 10LL) + 11LL*Power(xi, 9LL)*xj + 55LL*Power(xi, 8LL)*Power(xj, 2LL) +

                        165LL*Power(xi, 7LL)*Power(xj, 3LL) + 330LL*Power(xi, 6LL)*Power(xj, 4LL) +

                        462LL*Power(xi, 5LL)*Power(xj, 5LL) + 330LL*Power(xi, 4LL)*Power(xj, 6LL) +

                        165LL*Power(xi, 3LL)*Power(xj, 7LL) + 55LL*Power(xi, 2LL)*Power(xj, 8LL) +

                        11LL*xi*Power(xj, 9LL) + Power(xj, 10LL)))/(3LL*Power(xi + xj, 11LL))

            ;
        }
        else
        {
            S = (1LL/r)*((135LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 11LL) +

                          exp(2LL*rxj)*Power(rxj, 8LL)*

                          (-150LL*Power(rxi, 18LL) - 6LL*Power(rxi, 19LL) + 135LL*Power(rxj, 14LL) +

                           225LL*rxi*Power(rxj, 14LL) + 10LL*Power(rxi, 17LL)*(-165LL + Power(rxj, 2LL)) -

                           30LL*Power(rxi, 16LL)*(330LL + Power(rxj, 2LL)) +

                           45LL*Power(rxi, 3LL)*Power(rxj, 12LL)*(-55LL + 2LL*Power(rxj, 2LL)) +

                           45LL*Power(rxi, 2LL)*Power(rxj, 12LL)*(-33LL + 4LL*Power(rxj, 2LL)) +

                           Power(rxi, 9LL)*Power(rxj, 6LL)*

                           (234135LL - 4950LL*Power(rxj, 2LL) - 34LL*Power(rxj, 4LL)) -

                           5LL*Power(rxi, 7LL)*Power(rxj, 8LL)*

                           (6237LL - 1242LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           3LL*Power(rxi, 5LL)*Power(rxj, 10LL)*

                           (4125LL - 330LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 4LL)*Power(rxj, 10LL)*

                           (495LL - 132LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) -

                           165LL*Power(rxi, 6LL)*Power(rxj, 8LL)*

                           (135LL - 60LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) -

                           5LL*Power(rxi, 13LL)*Power(rxj, 2LL)*

                           (43875LL - 3438LL*Power(rxj, 2LL) + 22LL*Power(rxj, 4LL)) +

                           5LL*Power(rxi, 11LL)*Power(rxj, 4LL)*

                           (7695LL - 2442LL*Power(rxj, 2LL) + 22LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 8LL)*Power(rxj, 6LL)*

                           (-33LL - 3564LL*Power(rxj, 2LL) + 26LL*Power(rxj, 4LL)) +

                           Power(rxi, 15LL)*(-32175LL - 3690LL*Power(rxj, 2LL) + 34LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 10LL)*Power(rxj, 4LL)*

                           (-32277LL + 1364LL*Power(rxj, 2LL) + 66LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 14LL)*(-3003LL - 2932LL*Power(rxj, 2LL) + 94LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 12LL)*Power(rxj, 2LL)*

                           (28119LL - 5252LL*Power(rxj, 2LL) + 154LL*Power(rxj, 4LL))) +

                          exp(2LL*rxi)*Power(rxi, 8LL)*

                          (-5LL*Power(rxi, 2LL)*Power(rxj, 12LL)*

                           (-84357LL - 43875LL*rxj - 8796LL*Power(rxj, 2LL) - 738LL*Power(rxj, 3LL) -

           6LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) -

                           3LL*Power(rxi, 14LL)*(45LL + 75LL*rxj + 60LL*Power(rxj, 2LL) + 30LL*Power(rxj, 3LL) +

                                                 10LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) -

                           55LL*Power(rxi, 8LL)*Power(rxj, 6LL)*

                           (-405LL - 567LL*rxj - 972LL*Power(rxj, 2LL) - 90LL*Power(rxj, 3LL) +

           18LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           55LL*Power(rxi, 6LL)*Power(rxj, 8LL)*

                           (9LL - 4257LL*rxj - 372LL*Power(rxj, 2LL) + 222LL*Power(rxj, 3LL) +

           42LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           3LL*Power(rxj, 14LL)*(15015LL + 10725LL*rxj + 3300LL*Power(rxj, 2LL) +

                                                 550LL*Power(rxj, 3LL) + 50LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           5LL*Power(rxi, 12LL)*Power(rxj, 2LL)*

                           (297LL + 495LL*rxj + 396LL*Power(rxj, 2LL) + 198LL*Power(rxj, 3LL) +

           66LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           Power(rxi, 10LL)*Power(rxj, 4LL)*

                           (-7425LL - 12375LL*rxj - 9900LL*Power(rxj, 2LL) - 6210LL*Power(rxj, 3LL) -

           390LL*Power(rxj, 4LL) + 34LL*Power(rxj, 5LL)) -

                           Power(rxi, 4LL)*Power(rxj, 10LL)*

                           (-484155LL + 38475LL*rxj + 78780LL*Power(rxj, 2LL) + 17190LL*Power(rxj, 3LL) +

           1410LL*Power(rxj, 4LL) + 34LL*Power(rxj, 5LL))))/

                         (135LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 11LL)*Power(rxi + rxj, 11LL))

                         );
        }

    }
    return S;
}
