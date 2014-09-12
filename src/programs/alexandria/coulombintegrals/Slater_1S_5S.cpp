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

cl_R Slater_1S_5S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (2041LL*xi)/10240LL

            ;
        }
        else
        {
            S = (1LL/r)*((-1596672000LL + 1596672000LL*exp(2LL*rxi) - 2875101075LL*rxi -

                          2556858150LL*Power(rxi, 2LL) - 1492929900LL*Power(rxi, 3LL) -

                          641163600LL*Power(rxi, 4LL) - 214719120LL*Power(rxi, 5LL) -

                          57879360LL*Power(rxi, 6LL) - 12735360LL*Power(rxi, 7LL) - 2280960LL*Power(rxi, 8LL) -

                          323840LL*Power(rxi, 9LL) - 33792LL*Power(rxi, 10LL) - 2048LL*Power(rxi, 11LL))/

                         (1.596672e9*exp(2LL*rxi))

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

                        330LL*Power(xi, 3LL)*Power(xj, 7LL) + 165LL*Power(xi, 2LL)*Power(xj, 8LL) +

                        55LL*xi*Power(xj, 9LL) + 5LL*Power(xj, 10LL)))/(5LL*Power(xi + xj, 11LL))

            ;
        }
        else
        {
            S = (1LL/r)*((14175LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 11LL) +

                          2835LL*exp(2LL*rxj)*Power(rxj, 12LL)*

                          (-35LL*Power(rxi, 10LL) - 5LL*Power(rxi, 11LL) - 495LL*Power(rxi, 8LL)*Power(rxj, 2LL) -

                           55LL*Power(rxi, 9LL)*Power(rxj, 2LL) - 1254LL*Power(rxi, 6LL)*Power(rxj, 4LL) -

                           66LL*Power(rxi, 7LL)*Power(rxj, 4LL) - 726LL*Power(rxi, 4LL)*Power(rxj, 6LL) +

                           66LL*Power(rxi, 5LL)*Power(rxj, 6LL) - 55LL*Power(rxi, 2LL)*Power(rxj, 8LL) +

                           55LL*Power(rxi, 3LL)*Power(rxj, 8LL) + 5LL*Power(rxj, 10LL) + 5LL*rxi*Power(rxj, 10LL))

                          - exp(2LL*rxi)*Power(rxi, 4LL)*(Power(rxi, 18LL)*

                                                          (14175LL + 25515LL*rxj + 22680LL*Power(rxj, 2LL) + 13230LL*Power(rxj, 3LL) +

                                                           5670LL*Power(rxj, 4LL) + 1890LL*Power(rxj, 5LL) + 504LL*Power(rxj, 6LL) +

                                                           108LL*Power(rxj, 7LL) + 18LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) -

                                                          9LL*Power(rxi, 16LL)*Power(rxj, 2LL)*

                                                          (17325LL + 31185LL*rxj + 27720LL*Power(rxj, 2LL) + 16170LL*Power(rxj, 3LL) +

                                                           6930LL*Power(rxj, 4LL) + 2310LL*Power(rxj, 5LL) + 616LL*Power(rxj, 6LL) +

                                                           132LL*Power(rxj, 7LL) + 22LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) +

                                                          126LL*Power(rxi, 10LL)*Power(rxj, 8LL)*

                                                          (37125LL + 66825LL*rxj + 59400LL*Power(rxj, 2LL) + 34725LL*Power(rxj, 3LL) +

                                                           14625LL*Power(rxj, 4LL) + 5043LL*Power(rxj, 5LL) + 1396LL*Power(rxj, 6LL) +

                                                           276LL*Power(rxj, 7LL) + 34LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) -

                                                          126LL*Power(rxi, 8LL)*Power(rxj, 10LL)*

                                                          (51975LL + 93420LL*rxj + 84240LL*Power(rxj, 2LL) + 46815LL*Power(rxj, 3LL) +

                                                           20835LL*Power(rxj, 4LL) + 7485LL*Power(rxj, 5LL) + 1964LL*Power(rxj, 6LL) +

                                                           348LL*Power(rxj, 7LL) + 38LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) +

                                                          9LL*Power(rxi, 2LL)*Power(rxj, 16LL)*

                                                          (-135135LL + 405405LL*rxj + 582120LL*Power(rxj, 2LL) + 346500LL*Power(rxj, 3LL) +

                                                           124740LL*Power(rxj, 4LL) + 30492LL*Power(rxj, 5LL) + 5264LL*Power(rxj, 6LL) +

                                                           636LL*Power(rxj, 7LL) + 50LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) -

                                                          Power(rxj, 18LL)*(2837835LL + 3648645LL*rxj + 2245320LL*Power(rxj, 2LL) +

                                                                            873180LL*Power(rxj, 3LL) + 238140LL*Power(rxj, 4LL) + 47628LL*Power(rxj, 5LL) +

                                                                            7056LL*Power(rxj, 6LL) + 756LL*Power(rxj, 7LL) + 54LL*Power(rxj, 8LL) +

                                                                            2LL*Power(rxj, 9LL)) + 9LL*Power(rxi, 14LL)*Power(rxj, 4LL)*

                                                          (86625LL + 155925LL*rxj + 138600LL*Power(rxj, 2LL) + 80850LL*Power(rxj, 3LL) +

                                                           34650LL*Power(rxj, 4LL) + 11550LL*Power(rxj, 5LL) + 3080LL*Power(rxj, 6LL) +

                                                           672LL*Power(rxj, 7LL) + 104LL*Power(rxj, 8LL) + 8LL*Power(rxj, 9LL)) -

                                                          21LL*Power(rxi, 12LL)*Power(rxj, 6LL)*

                                                          (111375LL + 200475LL*rxj + 178200LL*Power(rxj, 2LL) + 103950LL*Power(rxj, 3LL) +

                                                           44550LL*Power(rxj, 4LL) + 14778LL*Power(rxj, 5LL) + 4056LL*Power(rxj, 6LL) +

                                                           864LL*Power(rxj, 7LL) + 120LL*Power(rxj, 8LL) + 8LL*Power(rxj, 9LL)) +

                                                          21LL*Power(rxi, 6LL)*Power(rxj, 12LL)*

                                                          (307125LL + 594945LL*rxj + 456840LL*Power(rxj, 2LL) + 281790LL*Power(rxj, 3LL) +

                                                           137430LL*Power(rxj, 4LL) + 47250LL*Power(rxj, 5LL) + 11064LL*Power(rxj, 6LL) +

                                                           1728LL*Power(rxj, 7LL) + 168LL*Power(rxj, 8LL) + 8LL*Power(rxj, 9LL)) -

                                                          9LL*Power(rxi, 4LL)*Power(rxj, 14LL)*

                                                          (675675LL + 675675LL*rxj + 748440LL*Power(rxj, 2LL) + 561330LL*Power(rxj, 3LL) +

                                                           256410LL*Power(rxj, 4LL) + 76230LL*Power(rxj, 5LL) + 15400LL*Power(rxj, 6LL) +

                                                           2112LL*Power(rxj, 7LL) + 184LL*Power(rxj, 8LL) + 8LL*Power(rxj, 9LL))))/

                         (14175LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 11LL)*Power(rxi + rxj, 11LL))

                         );
        }

    }
    return S;
}


cl_R Slater_5S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_1S_5S(r, xj, xi);
}
