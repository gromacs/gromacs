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

cl_R Slater_2S_5S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (2011LL*xi)/10240LL

            ;
        }
        else
        {
            S = (1LL/r)*((-124540416000LL + 124540416000LL*exp(2LL*rxi) - 224622748350LL*rxi -

                          200164664700LL*Power(rxi, 2LL) - 117249207075LL*Power(rxi, 3LL) -

                          50639138550LL*Power(rxi, 4LL) - 17132415300LL*Power(rxi, 5LL) -

                          4704860160LL*Power(rxi, 6LL) - 1071195840LL*Power(rxi, 7LL) -

                          204478560LL*Power(rxi, 8LL) - 32809920LL*Power(rxi, 9LL) - 4392960LL*Power(rxi, 10LL) -

                          479232LL*Power(rxi, 11LL) - 39936LL*Power(rxi, 12LL) - 2048LL*Power(rxi, 13LL))/

                         (1.24540416e11*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(2LL*Power(xi, 12LL) + 26LL*Power(xi, 11LL)*xj + 156LL*Power(xi, 10LL)*Power(xj, 2LL) +

                        572LL*Power(xi, 9LL)*Power(xj, 3LL) + 1430LL*Power(xi, 8LL)*Power(xj, 4LL) +

                        2574LL*Power(xi, 7LL)*Power(xj, 5LL) + 3432LL*Power(xi, 6LL)*Power(xj, 6LL) +

                        3432LL*Power(xi, 5LL)*Power(xj, 7LL) + 2574LL*Power(xi, 4LL)*Power(xj, 8LL) +

                        1430LL*Power(xi, 3LL)*Power(xj, 9LL) + 390LL*Power(xi, 2LL)*Power(xj, 10LL) +

                        65LL*xi*Power(xj, 11LL) + 5LL*Power(xj, 12LL)))/(10LL*Power(xi + xj, 13LL))

            ;
        }
        else
        {
            S = (1LL/r)*((28350LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 13LL) +

                          945LL*exp(2LL*rxj)*Power(rxj, 12LL)*

                          (-210LL*Power(rxi, 16LL) - 10LL*Power(rxi, 17LL) + 30LL*Power(rxj, 14LL) +

                           45LL*rxi*Power(rxj, 14LL) + 39LL*Power(rxi, 7LL)*Power(rxj, 8LL)*

                           (1309LL - 2LL*Power(rxj, 2LL)) +

                           858LL*Power(rxi, 8LL)*Power(rxj, 6LL)*(-305LL + Power(rxj, 2LL)) +

                           30LL*Power(rxi, 2LL)*Power(rxj, 12LL)*(-13LL + Power(rxj, 2LL)) -

                           390LL*Power(rxi, 4LL)*Power(rxj, 10LL)*(-6LL + Power(rxj, 2LL)) -

                           143LL*Power(rxi, 9LL)*Power(rxj, 6LL)*(-153LL + 2LL*Power(rxj, 2LL)) +

                           5LL*Power(rxi, 3LL)*Power(rxj, 12LL)*(-117LL + 2LL*Power(rxj, 2LL)) -

                           45LL*Power(rxi, 15LL)*(35LL + 2LL*Power(rxj, 2LL)) -

                           138LL*Power(rxi, 12LL)*Power(rxj, 2LL)*(580LL + 13LL*Power(rxj, 2LL)) -

                           150LL*Power(rxi, 14LL)*(28LL + 17LL*Power(rxj, 2LL)) +

                           13LL*Power(rxi, 11LL)*Power(rxj, 4LL)*(-4071LL + 22LL*Power(rxj, 2LL)) +

                           3LL*Power(rxi, 13LL)*Power(rxj, 2LL)*(-8135LL + 26LL*Power(rxj, 2LL)) +

                           3LL*Power(rxi, 5LL)*Power(rxj, 10LL)*(2171LL + 30LL*Power(rxj, 2LL)) +

                           234LL*Power(rxi, 10LL)*Power(rxj, 4LL)*(-1235LL + 33LL*Power(rxj, 2LL)) -

                           78LL*Power(rxi, 6LL)*Power(rxj, 8LL)*(550LL + 47LL*Power(rxj, 2LL))) -

                          2LL*exp(2LL*rxi)*Power(rxi, 6LL)*

                          (-819LL*Power(rxi, 10LL)*Power(rxj, 10LL)*

                           (22275LL + 39780LL*rxj + 38160LL*Power(rxj, 2LL) + 16560LL*Power(rxj, 3LL) +

           9840LL*Power(rxj, 4LL) + 3900LL*Power(rxj, 5LL) + 816LL*Power(rxj, 6LL) +

           88LL*Power(rxj, 7LL) + 4LL*Power(rxj, 8LL)) +

                           Power(rxi, 20LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj, 2LL) +

                                             13230LL*Power(rxj, 3LL) + 5670LL*Power(rxj, 4LL) + 1890LL*Power(rxj, 5LL) +

                                             504LL*Power(rxj, 6LL) + 108LL*Power(rxj, 7LL) + 18LL*Power(rxj, 8LL) +

                                             2LL*Power(rxj, 9LL)) - Power(rxj, 20LL)*

                           (16216200LL + 18243225LL*rxj + 9729720LL*Power(rxj, 2LL) +

           3243240LL*Power(rxj, 3LL) + 748440LL*Power(rxj, 4LL) +

           124740LL*Power(rxj, 5LL) + 15120LL*Power(rxj, 6LL) + 1296LL*Power(rxj, 7LL) +

           72LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) +

                           18LL*Power(rxi, 16LL)*Power(rxj, 4LL)*

                           (61425LL + 110565LL*rxj + 98280LL*Power(rxj, 2LL) + 57330LL*Power(rxj, 3LL) +

           24570LL*Power(rxj, 4LL) + 8190LL*Power(rxj, 5LL) + 2184LL*Power(rxj, 6LL) +

           496LL*Power(rxj, 7LL) + 64LL*Power(rxj, 8LL) + 3LL*Power(rxj, 9LL)) -

                           18LL*Power(rxi, 4LL)*Power(rxj, 16LL)*

                           (6572475LL - 3161340LL*rxj - 4782960LL*Power(rxj, 2LL) -

           1912365LL*Power(rxj, 3LL) - 378105LL*Power(rxj, 4LL) - 34125LL*Power(rxj, 5LL) +

           1092LL*Power(rxj, 6LL) + 650LL*Power(rxj, 7LL) + 71LL*Power(rxj, 8LL) +

           3LL*Power(rxj, 9LL)) - 21LL*Power(rxi, 8LL)*Power(rxj, 12LL)*

                           (-1063800LL - 2775735LL*rxj - 862920LL*Power(rxj, 2LL) -

           1132020LL*Power(rxj, 3LL) - 698580LL*Power(rxj, 4LL) -

           196920LL*Power(rxj, 5LL) - 28992LL*Power(rxj, 6LL) - 2064LL*Power(rxj, 7LL) -

           24LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) +

                           21LL*Power(rxi, 12LL)*Power(rxj, 8LL)*

                           (482625LL + 868725LL*rxj + 772200LL*Power(rxj, 2LL) + 455400LL*Power(rxj, 3LL) +

           178200LL*Power(rxj, 4LL) + 72180LL*Power(rxj, 5LL) + 19920LL*Power(rxj, 6LL) +

           2952LL*Power(rxj, 7LL) + 204LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) +

                           6LL*Power(rxi, 6LL)*Power(rxj, 14LL)*

                           (-10357200LL + 5071815LL*rxj - 6463800LL*Power(rxj, 2LL) -

           7151130LL*Power(rxj, 3LL) - 2572290LL*Power(rxj, 4LL) -

           468720LL*Power(rxj, 5LL) - 42672LL*Power(rxj, 6LL) - 648LL*Power(rxj, 7LL) +

           228LL*Power(rxj, 8LL) + 16LL*Power(rxj, 9LL)) -

                           Power(rxi, 18LL)*Power(rxj, 2LL)*

                           (184275LL + 331695LL*rxj + 294840LL*Power(rxj, 2LL) + 171990LL*Power(rxj, 3LL) +

           73710LL*Power(rxj, 4LL) + 24570LL*Power(rxj, 5LL) + 6552LL*Power(rxj, 6LL) +

           1404LL*Power(rxj, 7LL) + 234LL*Power(rxj, 8LL) + 16LL*Power(rxj, 9LL)) +

                           Power(rxi, 2LL)*Power(rxj, 18LL)*

                           (-133783650LL - 107432325LL*rxj - 35675640LL*Power(rxj, 2LL) -

           5135130LL*Power(rxj, 3LL) + 270270LL*Power(rxj, 4LL) +

           270270LL*Power(rxj, 5LL) + 57960LL*Power(rxj, 6LL) + 6948LL*Power(rxj, 7LL) +

           486LL*Power(rxj, 8LL) + 16LL*Power(rxj, 9LL)) -

                           6LL*Power(rxi, 14LL)*Power(rxj, 6LL)*

                           (675675LL + 1216215LL*rxj + 1081080LL*Power(rxj, 2LL) + 630630LL*Power(rxj, 3LL) +

           270270LL*Power(rxj, 4LL) + 88200LL*Power(rxj, 5LL) + 26544LL*Power(rxj, 6LL) +

           5160LL*Power(rxj, 7LL) + 492LL*Power(rxj, 8LL) + 16LL*Power(rxj, 9LL))))/

                         (28350LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 13LL)*Power(rxi + rxj, 13LL))

                         );
        }

    }
    return S;
}


cl_R Slater_5S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_2S_5S(r, xj, xi);
}
