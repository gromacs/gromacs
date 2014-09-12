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

cl_R DSlater_1S_6S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-137006619750LL*xi + 149448499200LL*exp(2LL*r*xi)*xi -

                  249129480600LL*r*Power(xi, 2LL) - 224263964925LL*Power(r, 2LL)*Power(xi, 3LL) -

                  132956623800LL*Power(r, 3LL)*Power(xi, 4LL) -

                  58224266100LL*Power(r, 4LL)*Power(xi, 5LL) -

                  20004304320LL*Power(r, 5LL)*Power(xi, 6LL) -

                  5582697120LL*Power(r, 6LL)*Power(xi, 7LL) -

                  1289882880LL*Power(r, 7LL)*Power(xi, 8LL) -

                  248339520LL*Power(r, 8LL)*Power(xi, 9LL) - 39536640LL*Power(r, 9LL)*Power(xi, 10LL) -

                  5051904LL*Power(r, 10LL)*Power(xi, 11LL) - 479232LL*Power(r, 11LL)*Power(xi, 12LL) -

                  26624LL*Power(r, 12LL)*Power(xi, 13LL))/(7.47242496e10*exp(2LL*r*xi)*r) +

                (-74724249600LL + 74724249600LL*exp(2LL*r*xi) - 137006619750LL*r*xi -

                 124564740300LL*Power(r, 2LL)*Power(xi, 2LL) -

                 74754654975LL*Power(r, 3LL)*Power(xi, 3LL) -

                 33239155950LL*Power(r, 4LL)*Power(xi, 4LL) -

                 11644853220LL*Power(r, 5LL)*Power(xi, 5LL) -

                 3334050720LL*Power(r, 6LL)*Power(xi, 6LL) - 797528160LL*Power(r, 7LL)*Power(xi, 7LL) -

                 161235360LL*Power(r, 8LL)*Power(xi, 8LL) - 27593280LL*Power(r, 9LL)*Power(xi, 9LL) -

                 3953664LL*Power(r, 10LL)*Power(xi, 10LL) - 459264LL*Power(r, 11LL)*Power(xi, 11LL) -

                 39936LL*Power(r, 12LL)*Power(xi, 12LL) - 2048LL*Power(r, 13LL)*Power(xi, 13LL))/

                (7.47242496e10*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-74724249600LL + 74724249600LL*exp(2LL*r*xi) - 137006619750LL*r*xi -

                     124564740300LL*Power(r, 2LL)*Power(xi, 2LL) -

                     74754654975LL*Power(r, 3LL)*Power(xi, 3LL) -

                     33239155950LL*Power(r, 4LL)*Power(xi, 4LL) -

                     11644853220LL*Power(r, 5LL)*Power(xi, 5LL) -

                     3334050720LL*Power(r, 6LL)*Power(xi, 6LL) -

                     797528160LL*Power(r, 7LL)*Power(xi, 7LL) - 161235360LL*Power(r, 8LL)*Power(xi, 8LL) -

                     27593280LL*Power(r, 9LL)*Power(xi, 9LL) - 3953664LL*Power(r, 10LL)*Power(xi, 10LL) -

                     459264LL*Power(r, 11LL)*Power(xi, 11LL) - 39936LL*Power(r, 12LL)*Power(xi, 12LL) -

                     2048LL*Power(r, 13LL)*Power(xi, 13LL)))/(3.73621248e10*exp(2LL*r*xi)*r)

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
            S = (935550LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                 311850LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-24LL*Power(xi, 12LL) - 3LL*r*Power(xi, 13LL) -

                  507LL*Power(xi, 10LL)*Power(xj, 2LL) - 52LL*r*Power(xi, 11LL)*Power(xj, 2LL) -

                  2145LL*Power(xi, 8LL)*Power(xj, 4LL) - 143LL*r*Power(xi, 9LL)*Power(xj, 4LL) -

                  2574LL*Power(xi, 6LL)*Power(xj, 6LL) - 858LL*Power(xi, 4LL)*Power(xj, 8LL) +

                  143LL*r*Power(xi, 5LL)*Power(xj, 8LL) - 39LL*Power(xi, 2LL)*Power(xj, 10LL) +

                  52LL*r*Power(xi, 3LL)*Power(xj, 10LL) + 3LL*Power(xj, 12LL) +

                  3LL*r*xi*Power(xj, 12LL)) +

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-110LL*Power(xi, 18LL)*Power(xj, 4LL)*

                  (663390LL + 1216215LL*r*xj + 1105650LL*Power(r, 2LL)*Power(xj, 2LL) +

            663390LL*Power(r, 3LL)*Power(xj, 3LL) + 294840LL*Power(r, 4LL)*Power(xj, 4LL) +

            103194LL*Power(r, 5LL)*Power(xj, 5LL) + 29484LL*Power(r, 6LL)*Power(xj, 6LL) +

            7020LL*Power(r, 7LL)*Power(xj, 7LL) + 1404LL*Power(r, 8LL)*Power(xj, 8LL) +

            237LL*Power(r, 9LL)*Power(xj, 9LL) + 30LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  330LL*Power(xi, 16LL)*Power(xj, 6LL)*

                  (810810LL + 1486485LL*r*xj + 1351350LL*Power(r, 2LL)*Power(xj, 2LL) +

            810810LL*Power(r, 3LL)*Power(xj, 3LL) + 360360LL*Power(r, 4LL)*Power(xj, 4LL) +

            126126LL*Power(r, 5LL)*Power(xj, 5LL) + 36036LL*Power(r, 6LL)*Power(xj, 6LL) +

            8556LL*Power(r, 7LL)*Power(xj, 7LL) + 1740LL*Power(r, 8LL)*Power(xj, 8LL) +

            291LL*Power(r, 9LL)*Power(xj, 9LL) + 34LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  330LL*Power(xi, 6LL)*Power(xj, 16LL)*

                  (3169530LL + 7960680LL*r*xj + 5798520LL*Power(r, 2LL)*Power(xj, 2LL) +

            3144960LL*Power(r, 3LL)*Power(xj, 3LL) +

            1572480LL*Power(r, 4LL)*Power(xj, 4LL) +

            638001LL*Power(r, 5LL)*Power(xj, 5LL) + 191646LL*Power(r, 6LL)*Power(xj, 6LL) +

            41886LL*Power(r, 7LL)*Power(xj, 7LL) + 6630LL*Power(r, 8LL)*Power(xj, 8LL) +

            741LL*Power(r, 9LL)*Power(xj, 9LL) + 54LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  110LL*Power(xi, 4LL)*Power(xj, 18LL)*

                  (12162150LL + 8108100LL*r*xj + 6486480LL*Power(r, 2LL)*Power(xj, 2LL) +

            5675670LL*Power(r, 3LL)*Power(xj, 3LL) +

            3243240LL*Power(r, 4LL)*Power(xj, 4LL) +

            1216215LL*Power(r, 5LL)*Power(xj, 5LL) +

            319410LL*Power(r, 6LL)*Power(xj, 6LL) + 61074LL*Power(r, 7LL)*Power(xj, 7LL) +

            8586LL*Power(r, 8LL)*Power(xj, 8LL) + 867LL*Power(r, 9LL)*Power(xj, 9LL) +

            58LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  Power(xi, 22LL)*(935550LL + 1715175LL*r*xj +

                                   1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                   330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                   4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  11LL*Power(xi, 20LL)*Power(xj, 2LL)*

                  (1105650LL + 2027025LL*r*xj + 1842750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1105650LL*Power(r, 3LL)*Power(xj, 3LL) +

            491400LL*Power(r, 4LL)*Power(xj, 4LL) + 171990LL*Power(r, 5LL)*Power(xj, 5LL) +

            49140LL*Power(r, 6LL)*Power(xj, 6LL) + 11700LL*Power(r, 7LL)*Power(xj, 7LL) +

            2340LL*Power(r, 8LL)*Power(xj, 8LL) + 390LL*Power(r, 9LL)*Power(xj, 9LL) +

            52LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  11LL*Power(xi, 2LL)*Power(xj, 20LL)*

                  (-48648600LL + 2027025LL*r*xj + 44594550LL*Power(r, 2LL)*Power(xj, 2LL) +

            36486450LL*Power(r, 3LL)*Power(xj, 3LL) +

            16216200LL*Power(r, 4LL)*Power(xj, 4LL) +

            4864860LL*Power(r, 5LL)*Power(xj, 5LL) +

            1065960LL*Power(r, 6LL)*Power(xj, 6LL) +

            176040LL*Power(r, 7LL)*Power(xj, 7LL) + 21960LL*Power(r, 8LL)*Power(xj, 8LL) +

            2010LL*Power(r, 9LL)*Power(xj, 9LL) + 124LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  Power(xj, 22LL)*(340540200LL + 468242775LL*r*xj +

                                   312161850LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   133783650LL*Power(r, 3LL)*Power(xj, 3LL) +

                                   41164200LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   9604980LL*Power(r, 5LL)*Power(xj, 5LL) +

                                   1746360LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   249480LL*Power(r, 7LL)*Power(xj, 7LL) + 27720LL*Power(r, 8LL)*Power(xj, 8LL) +

                                   2310LL*Power(r, 9LL)*Power(xj, 9LL) + 132LL*Power(r, 10LL)*Power(xj, 10LL) +

                                   4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  165LL*Power(xi, 14LL)*Power(xj, 8LL)*

                  (4054050LL + 7432425LL*r*xj + 6756750LL*Power(r, 2LL)*Power(xj, 2LL) +

            4054050LL*Power(r, 3LL)*Power(xj, 3LL) +

            1801800LL*Power(r, 4LL)*Power(xj, 4LL) +

            631260LL*Power(r, 5LL)*Power(xj, 5LL) + 178920LL*Power(r, 6LL)*Power(xj, 6LL) +

            43176LL*Power(r, 7LL)*Power(xj, 7LL) + 8904LL*Power(r, 8LL)*Power(xj, 8LL) +

            1428LL*Power(r, 9LL)*Power(xj, 9LL) + 152LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  231LL*Power(xi, 12LL)*Power(xj, 10LL)*

                  (5212350LL + 9555975LL*r*xj + 8687250LL*Power(r, 2LL)*Power(xj, 2LL) +

            5209650LL*Power(r, 3LL)*Power(xj, 3LL) +

            2327400LL*Power(r, 4LL)*Power(xj, 4LL) +

            801540LL*Power(r, 5LL)*Power(xj, 5LL) + 230040LL*Power(r, 6LL)*Power(xj, 6LL) +

            57240LL*Power(r, 7LL)*Power(xj, 7LL) + 11640LL*Power(r, 8LL)*Power(xj, 8LL) +

            1740LL*Power(r, 9LL)*Power(xj, 9LL) + 168LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  231LL*Power(xi, 10LL)*Power(xj, 12LL)*

                  (6949800LL + 12746025LL*r*xj + 11535750LL*Power(r, 2LL)*Power(xj, 2LL) +

            7056450LL*Power(r, 3LL)*Power(xj, 3LL) +

            3040200LL*Power(r, 4LL)*Power(xj, 4LL) +

            1051920LL*Power(r, 5LL)*Power(xj, 5LL) +

            316800LL*Power(r, 6LL)*Power(xj, 6LL) + 79680LL*Power(r, 7LL)*Power(xj, 7LL) +

            15360LL*Power(r, 8LL)*Power(xj, 8LL) + 2100LL*Power(r, 9LL)*Power(xj, 9LL) +

            184LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  165LL*Power(xi, 8LL)*Power(xj, 14LL)*

                  (9775080LL + 17424855LL*r*xj + 17019450LL*Power(r, 2LL)*Power(xj, 2LL) +

            9519930LL*Power(r, 3LL)*Power(xj, 3LL) +

            4059720LL*Power(r, 4LL)*Power(xj, 4LL) +

            1519056LL*Power(r, 5LL)*Power(xj, 5LL) + 475776LL*Power(r, 6LL)*Power(xj, 6LL) +

            114720LL*Power(r, 7LL)*Power(xj, 7LL) + 20256LL*Power(r, 8LL)*Power(xj, 8LL) +

            2508LL*Power(r, 9LL)*Power(xj, 9LL) + 200LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (935550LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 13LL)*

                 Power(xi + xj, 13LL)) + (935550LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                                          311850LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                                          (-24LL*Power(xi, 12LL) - 3LL*r*Power(xi, 13LL) -

                                  507LL*Power(xi, 10LL)*Power(xj, 2LL) - 52LL*r*Power(xi, 11LL)*Power(xj, 2LL) -

                                  2145LL*Power(xi, 8LL)*Power(xj, 4LL) - 143LL*r*Power(xi, 9LL)*Power(xj, 4LL) -

                                  2574LL*Power(xi, 6LL)*Power(xj, 6LL) - 858LL*Power(xi, 4LL)*Power(xj, 8LL) +

                                  143LL*r*Power(xi, 5LL)*Power(xj, 8LL) - 39LL*Power(xi, 2LL)*Power(xj, 10LL) +

                                  52LL*r*Power(xi, 3LL)*Power(xj, 10LL) + 3LL*Power(xj, 12LL) +

                                  3LL*r*xi*Power(xj, 12LL)) +

                                          exp(2LL*r*xi)*Power(xi, 4LL)*

                                          (-110LL*Power(xi, 18LL)*Power(xj, 4LL)*

                                  (663390LL + 1216215LL*r*xj + 1105650LL*Power(r, 2LL)*Power(xj, 2LL) +

            663390LL*Power(r, 3LL)*Power(xj, 3LL) + 294840LL*Power(r, 4LL)*Power(xj, 4LL) +

            103194LL*Power(r, 5LL)*Power(xj, 5LL) + 29484LL*Power(r, 6LL)*Power(xj, 6LL) +

            7020LL*Power(r, 7LL)*Power(xj, 7LL) + 1404LL*Power(r, 8LL)*Power(xj, 8LL) +

            237LL*Power(r, 9LL)*Power(xj, 9LL) + 30LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  330LL*Power(xi, 16LL)*Power(xj, 6LL)*

                                  (810810LL + 1486485LL*r*xj + 1351350LL*Power(r, 2LL)*Power(xj, 2LL) +

            810810LL*Power(r, 3LL)*Power(xj, 3LL) + 360360LL*Power(r, 4LL)*Power(xj, 4LL) +

            126126LL*Power(r, 5LL)*Power(xj, 5LL) + 36036LL*Power(r, 6LL)*Power(xj, 6LL) +

            8556LL*Power(r, 7LL)*Power(xj, 7LL) + 1740LL*Power(r, 8LL)*Power(xj, 8LL) +

            291LL*Power(r, 9LL)*Power(xj, 9LL) + 34LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  330LL*Power(xi, 6LL)*Power(xj, 16LL)*

                                  (3169530LL + 7960680LL*r*xj + 5798520LL*Power(r, 2LL)*Power(xj, 2LL) +

            3144960LL*Power(r, 3LL)*Power(xj, 3LL) +

            1572480LL*Power(r, 4LL)*Power(xj, 4LL) +

            638001LL*Power(r, 5LL)*Power(xj, 5LL) + 191646LL*Power(r, 6LL)*Power(xj, 6LL) +

            41886LL*Power(r, 7LL)*Power(xj, 7LL) + 6630LL*Power(r, 8LL)*Power(xj, 8LL) +

            741LL*Power(r, 9LL)*Power(xj, 9LL) + 54LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  110LL*Power(xi, 4LL)*Power(xj, 18LL)*

                                  (12162150LL + 8108100LL*r*xj + 6486480LL*Power(r, 2LL)*Power(xj, 2LL) +

            5675670LL*Power(r, 3LL)*Power(xj, 3LL) +

            3243240LL*Power(r, 4LL)*Power(xj, 4LL) +

            1216215LL*Power(r, 5LL)*Power(xj, 5LL) +

            319410LL*Power(r, 6LL)*Power(xj, 6LL) + 61074LL*Power(r, 7LL)*Power(xj, 7LL) +

            8586LL*Power(r, 8LL)*Power(xj, 8LL) + 867LL*Power(r, 9LL)*Power(xj, 9LL) +

            58LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  Power(xi, 22LL)*(935550LL + 1715175LL*r*xj +

                                                   1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                   330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                   4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  11LL*Power(xi, 20LL)*Power(xj, 2LL)*

                                  (1105650LL + 2027025LL*r*xj + 1842750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1105650LL*Power(r, 3LL)*Power(xj, 3LL) +

            491400LL*Power(r, 4LL)*Power(xj, 4LL) + 171990LL*Power(r, 5LL)*Power(xj, 5LL) +

            49140LL*Power(r, 6LL)*Power(xj, 6LL) + 11700LL*Power(r, 7LL)*Power(xj, 7LL) +

            2340LL*Power(r, 8LL)*Power(xj, 8LL) + 390LL*Power(r, 9LL)*Power(xj, 9LL) +

            52LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  11LL*Power(xi, 2LL)*Power(xj, 20LL)*

                                  (-48648600LL + 2027025LL*r*xj + 44594550LL*Power(r, 2LL)*Power(xj, 2LL) +

            36486450LL*Power(r, 3LL)*Power(xj, 3LL) +

            16216200LL*Power(r, 4LL)*Power(xj, 4LL) +

            4864860LL*Power(r, 5LL)*Power(xj, 5LL) +

            1065960LL*Power(r, 6LL)*Power(xj, 6LL) +

            176040LL*Power(r, 7LL)*Power(xj, 7LL) + 21960LL*Power(r, 8LL)*Power(xj, 8LL) +

            2010LL*Power(r, 9LL)*Power(xj, 9LL) + 124LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  Power(xj, 22LL)*(340540200LL + 468242775LL*r*xj +

                                                   312161850LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   133783650LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                   41164200LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   9604980LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                   1746360LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   249480LL*Power(r, 7LL)*Power(xj, 7LL) + 27720LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                   2310LL*Power(r, 9LL)*Power(xj, 9LL) + 132LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                   4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  165LL*Power(xi, 14LL)*Power(xj, 8LL)*

                                  (4054050LL + 7432425LL*r*xj + 6756750LL*Power(r, 2LL)*Power(xj, 2LL) +

            4054050LL*Power(r, 3LL)*Power(xj, 3LL) +

            1801800LL*Power(r, 4LL)*Power(xj, 4LL) +

            631260LL*Power(r, 5LL)*Power(xj, 5LL) + 178920LL*Power(r, 6LL)*Power(xj, 6LL) +

            43176LL*Power(r, 7LL)*Power(xj, 7LL) + 8904LL*Power(r, 8LL)*Power(xj, 8LL) +

            1428LL*Power(r, 9LL)*Power(xj, 9LL) + 152LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  231LL*Power(xi, 12LL)*Power(xj, 10LL)*

                                  (5212350LL + 9555975LL*r*xj + 8687250LL*Power(r, 2LL)*Power(xj, 2LL) +

            5209650LL*Power(r, 3LL)*Power(xj, 3LL) +

            2327400LL*Power(r, 4LL)*Power(xj, 4LL) +

            801540LL*Power(r, 5LL)*Power(xj, 5LL) + 230040LL*Power(r, 6LL)*Power(xj, 6LL) +

            57240LL*Power(r, 7LL)*Power(xj, 7LL) + 11640LL*Power(r, 8LL)*Power(xj, 8LL) +

            1740LL*Power(r, 9LL)*Power(xj, 9LL) + 168LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  231LL*Power(xi, 10LL)*Power(xj, 12LL)*

                                  (6949800LL + 12746025LL*r*xj + 11535750LL*Power(r, 2LL)*Power(xj, 2LL) +

            7056450LL*Power(r, 3LL)*Power(xj, 3LL) +

            3040200LL*Power(r, 4LL)*Power(xj, 4LL) +

            1051920LL*Power(r, 5LL)*Power(xj, 5LL) +

            316800LL*Power(r, 6LL)*Power(xj, 6LL) + 79680LL*Power(r, 7LL)*Power(xj, 7LL) +

            15360LL*Power(r, 8LL)*Power(xj, 8LL) + 2100LL*Power(r, 9LL)*Power(xj, 9LL) +

            184LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  165LL*Power(xi, 8LL)*Power(xj, 14LL)*

                                  (9775080LL + 17424855LL*r*xj + 17019450LL*Power(r, 2LL)*Power(xj, 2LL) +

            9519930LL*Power(r, 3LL)*Power(xj, 3LL) +

            4059720LL*Power(r, 4LL)*Power(xj, 4LL) +

            1519056LL*Power(r, 5LL)*Power(xj, 5LL) + 475776LL*Power(r, 6LL)*Power(xj, 6LL) +

            114720LL*Power(r, 7LL)*Power(xj, 7LL) + 20256LL*Power(r, 8LL)*Power(xj, 8LL) +

            2508LL*Power(r, 9LL)*Power(xj, 9LL) + 200LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (467775LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 13LL)*Power(xi + xj, 12LL)) -

                (1871100LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                 311850LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-3LL*Power(xi, 13LL) - 52LL*Power(xi, 11LL)*Power(xj, 2LL) -

        143LL*Power(xi, 9LL)*Power(xj, 4LL) + 143LL*Power(xi, 5LL)*Power(xj, 8LL) +

        52LL*Power(xi, 3LL)*Power(xj, 10LL) + 3LL*xi*Power(xj, 12LL)) +

                 623700LL*exp(2LL*r*xj)*Power(xj, 15LL)*

                 (-24LL*Power(xi, 12LL) - 3LL*r*Power(xi, 13LL) - 507LL*Power(xi, 10LL)*Power(xj, 2LL) -

        52LL*r*Power(xi, 11LL)*Power(xj, 2LL) - 2145LL*Power(xi, 8LL)*Power(xj, 4LL) -

        143LL*r*Power(xi, 9LL)*Power(xj, 4LL) - 2574LL*Power(xi, 6LL)*Power(xj, 6LL) -

        858LL*Power(xi, 4LL)*Power(xj, 8LL) + 143LL*r*Power(xi, 5LL)*Power(xj, 8LL) -

        39LL*Power(xi, 2LL)*Power(xj, 10LL) + 52LL*r*Power(xi, 3LL)*Power(xj, 10LL) +

        3LL*Power(xj, 12LL) + 3LL*r*xi*Power(xj, 12LL)) +

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-110LL*Power(xi, 18LL)*Power(xj, 4LL)*

        (1216215LL*xj + 2211300LL*r*Power(xj, 2LL) +

            1990170LL*Power(r, 2LL)*Power(xj, 3LL) +

            1179360LL*Power(r, 3LL)*Power(xj, 4LL) +

            515970LL*Power(r, 4LL)*Power(xj, 5LL) + 176904LL*Power(r, 5LL)*Power(xj, 6LL) +

            49140LL*Power(r, 6LL)*Power(xj, 7LL) + 11232LL*Power(r, 7LL)*Power(xj, 8LL) +

            2133LL*Power(r, 8LL)*Power(xj, 9LL) + 300LL*Power(r, 9LL)*Power(xj, 10LL) +

            22LL*Power(r, 10LL)*Power(xj, 11LL)) +

        330LL*Power(xi, 16LL)*Power(xj, 6LL)*

        (1486485LL*xj + 2702700LL*r*Power(xj, 2LL) +

            2432430LL*Power(r, 2LL)*Power(xj, 3LL) +

            1441440LL*Power(r, 3LL)*Power(xj, 4LL) +

            630630LL*Power(r, 4LL)*Power(xj, 5LL) + 216216LL*Power(r, 5LL)*Power(xj, 6LL) +

            59892LL*Power(r, 6LL)*Power(xj, 7LL) + 13920LL*Power(r, 7LL)*Power(xj, 8LL) +

            2619LL*Power(r, 8LL)*Power(xj, 9LL) + 340LL*Power(r, 9LL)*Power(xj, 10LL) +

            22LL*Power(r, 10LL)*Power(xj, 11LL)) -

        330LL*Power(xi, 6LL)*Power(xj, 16LL)*

        (7960680LL*xj + 11597040LL*r*Power(xj, 2LL) +

            9434880LL*Power(r, 2LL)*Power(xj, 3LL) +

            6289920LL*Power(r, 3LL)*Power(xj, 4LL) +

            3190005LL*Power(r, 4LL)*Power(xj, 5LL) +

            1149876LL*Power(r, 5LL)*Power(xj, 6LL) +

            293202LL*Power(r, 6LL)*Power(xj, 7LL) + 53040LL*Power(r, 7LL)*Power(xj, 8LL) +

            6669LL*Power(r, 8LL)*Power(xj, 9LL) + 540LL*Power(r, 9LL)*Power(xj, 10LL) +

            22LL*Power(r, 10LL)*Power(xj, 11LL)) +

        110LL*Power(xi, 4LL)*Power(xj, 18LL)*

        (8108100LL*xj + 12972960LL*r*Power(xj, 2LL) +

            17027010LL*Power(r, 2LL)*Power(xj, 3LL) +

            12972960LL*Power(r, 3LL)*Power(xj, 4LL) +

            6081075LL*Power(r, 4LL)*Power(xj, 5LL) +

            1916460LL*Power(r, 5LL)*Power(xj, 6LL) +

            427518LL*Power(r, 6LL)*Power(xj, 7LL) + 68688LL*Power(r, 7LL)*Power(xj, 8LL) +

            7803LL*Power(r, 8LL)*Power(xj, 9LL) + 580LL*Power(r, 9LL)*Power(xj, 10LL) +

            22LL*Power(r, 10LL)*Power(xj, 11LL)) -

        Power(xi, 22LL)*(1715175LL*xj + 3118500LL*r*Power(xj, 2LL) +

                         2806650LL*Power(r, 2LL)*Power(xj, 3LL) +

                         1663200LL*Power(r, 3LL)*Power(xj, 4LL) +

                         727650LL*Power(r, 4LL)*Power(xj, 5LL) + 249480LL*Power(r, 5LL)*Power(xj, 6LL) +

                         69300LL*Power(r, 6LL)*Power(xj, 7LL) + 15840LL*Power(r, 7LL)*Power(xj, 8LL) +

                         2970LL*Power(r, 8LL)*Power(xj, 9LL) + 440LL*Power(r, 9LL)*Power(xj, 10LL) +

                         44LL*Power(r, 10LL)*Power(xj, 11LL)) +

        11LL*Power(xi, 20LL)*Power(xj, 2LL)*

        (2027025LL*xj + 3685500LL*r*Power(xj, 2LL) +

            3316950LL*Power(r, 2LL)*Power(xj, 3LL) +

            1965600LL*Power(r, 3LL)*Power(xj, 4LL) +

            859950LL*Power(r, 4LL)*Power(xj, 5LL) + 294840LL*Power(r, 5LL)*Power(xj, 6LL) +

            81900LL*Power(r, 6LL)*Power(xj, 7LL) + 18720LL*Power(r, 7LL)*Power(xj, 8LL) +

            3510LL*Power(r, 8LL)*Power(xj, 9LL) + 520LL*Power(r, 9LL)*Power(xj, 10LL) +

            44LL*Power(r, 10LL)*Power(xj, 11LL)) -

        11LL*Power(xi, 2LL)*Power(xj, 20LL)*

        (2027025LL*xj + 89189100LL*r*Power(xj, 2LL) +

            109459350LL*Power(r, 2LL)*Power(xj, 3LL) +

            64864800LL*Power(r, 3LL)*Power(xj, 4LL) +

            24324300LL*Power(r, 4LL)*Power(xj, 5LL) +

            6395760LL*Power(r, 5LL)*Power(xj, 6LL) +

            1232280LL*Power(r, 6LL)*Power(xj, 7LL) +

            175680LL*Power(r, 7LL)*Power(xj, 8LL) + 18090LL*Power(r, 8LL)*Power(xj, 9LL) +

            1240LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) +

        Power(xj, 22LL)*(468242775LL*xj + 624323700LL*r*Power(xj, 2LL) +

                         401350950LL*Power(r, 2LL)*Power(xj, 3LL) +

                         164656800LL*Power(r, 3LL)*Power(xj, 4LL) +

                         48024900LL*Power(r, 4LL)*Power(xj, 5LL) +

                         10478160LL*Power(r, 5LL)*Power(xj, 6LL) +

                         1746360LL*Power(r, 6LL)*Power(xj, 7LL) +

                         221760LL*Power(r, 7LL)*Power(xj, 8LL) + 20790LL*Power(r, 8LL)*Power(xj, 9LL) +

                         1320LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) -

        165LL*Power(xi, 14LL)*Power(xj, 8LL)*

        (7432425LL*xj + 13513500LL*r*Power(xj, 2LL) +

            12162150LL*Power(r, 2LL)*Power(xj, 3LL) +

            7207200LL*Power(r, 3LL)*Power(xj, 4LL) +

            3156300LL*Power(r, 4LL)*Power(xj, 5LL) +

            1073520LL*Power(r, 5LL)*Power(xj, 6LL) +

            302232LL*Power(r, 6LL)*Power(xj, 7LL) + 71232LL*Power(r, 7LL)*Power(xj, 8LL) +

            12852LL*Power(r, 8LL)*Power(xj, 9LL) + 1520LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) +

        231LL*Power(xi, 12LL)*Power(xj, 10LL)*

        (9555975LL*xj + 17374500LL*r*Power(xj, 2LL) +

            15628950LL*Power(r, 2LL)*Power(xj, 3LL) +

            9309600LL*Power(r, 3LL)*Power(xj, 4LL) +

            4007700LL*Power(r, 4LL)*Power(xj, 5LL) +

            1380240LL*Power(r, 5LL)*Power(xj, 6LL) +

            400680LL*Power(r, 6LL)*Power(xj, 7LL) + 93120LL*Power(r, 7LL)*Power(xj, 8LL) +

            15660LL*Power(r, 8LL)*Power(xj, 9LL) + 1680LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) -

        231LL*Power(xi, 10LL)*Power(xj, 12LL)*

        (12746025LL*xj + 23071500LL*r*Power(xj, 2LL) +

            21169350LL*Power(r, 2LL)*Power(xj, 3LL) +

            12160800LL*Power(r, 3LL)*Power(xj, 4LL) +

            5259600LL*Power(r, 4LL)*Power(xj, 5LL) +

            1900800LL*Power(r, 5LL)*Power(xj, 6LL) +

            557760LL*Power(r, 6LL)*Power(xj, 7LL) + 122880LL*Power(r, 7LL)*Power(xj, 8LL) +

            18900LL*Power(r, 8LL)*Power(xj, 9LL) + 1840LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) +

        165LL*Power(xi, 8LL)*Power(xj, 14LL)*

        (17424855LL*xj + 34038900LL*r*Power(xj, 2LL) +

            28559790LL*Power(r, 2LL)*Power(xj, 3LL) +

            16238880LL*Power(r, 3LL)*Power(xj, 4LL) +

            7595280LL*Power(r, 4LL)*Power(xj, 5LL) +

            2854656LL*Power(r, 5LL)*Power(xj, 6LL) + 803040LL*Power(r, 6LL)*Power(xj, 7LL) +

            162048LL*Power(r, 7LL)*Power(xj, 8LL) + 22572LL*Power(r, 8LL)*Power(xj, 9LL) +

            2000LL*Power(r, 9LL)*Power(xj, 10LL) + 88LL*Power(r, 10LL)*Power(xj, 11LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 5LL)*

                 (-110LL*Power(xi, 18LL)*Power(xj, 4LL)*

        (663390LL + 1216215LL*r*xj + 1105650LL*Power(r, 2LL)*Power(xj, 2LL) +

            663390LL*Power(r, 3LL)*Power(xj, 3LL) + 294840LL*Power(r, 4LL)*Power(xj, 4LL) +

            103194LL*Power(r, 5LL)*Power(xj, 5LL) + 29484LL*Power(r, 6LL)*Power(xj, 6LL) +

            7020LL*Power(r, 7LL)*Power(xj, 7LL) + 1404LL*Power(r, 8LL)*Power(xj, 8LL) +

            237LL*Power(r, 9LL)*Power(xj, 9LL) + 30LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

        330LL*Power(xi, 16LL)*Power(xj, 6LL)*

        (810810LL + 1486485LL*r*xj + 1351350LL*Power(r, 2LL)*Power(xj, 2LL) +

            810810LL*Power(r, 3LL)*Power(xj, 3LL) + 360360LL*Power(r, 4LL)*Power(xj, 4LL) +

            126126LL*Power(r, 5LL)*Power(xj, 5LL) + 36036LL*Power(r, 6LL)*Power(xj, 6LL) +

            8556LL*Power(r, 7LL)*Power(xj, 7LL) + 1740LL*Power(r, 8LL)*Power(xj, 8LL) +

            291LL*Power(r, 9LL)*Power(xj, 9LL) + 34LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) -

        330LL*Power(xi, 6LL)*Power(xj, 16LL)*

        (3169530LL + 7960680LL*r*xj + 5798520LL*Power(r, 2LL)*Power(xj, 2LL) +

            3144960LL*Power(r, 3LL)*Power(xj, 3LL) +

            1572480LL*Power(r, 4LL)*Power(xj, 4LL) + 638001LL*Power(r, 5LL)*Power(xj, 5LL) +

            191646LL*Power(r, 6LL)*Power(xj, 6LL) + 41886LL*Power(r, 7LL)*Power(xj, 7LL) +

            6630LL*Power(r, 8LL)*Power(xj, 8LL) + 741LL*Power(r, 9LL)*Power(xj, 9LL) +

            54LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

        110LL*Power(xi, 4LL)*Power(xj, 18LL)*

        (12162150LL + 8108100LL*r*xj + 6486480LL*Power(r, 2LL)*Power(xj, 2LL) +

            5675670LL*Power(r, 3LL)*Power(xj, 3LL) +

            3243240LL*Power(r, 4LL)*Power(xj, 4LL) +

            1216215LL*Power(r, 5LL)*Power(xj, 5LL) + 319410LL*Power(r, 6LL)*Power(xj, 6LL) +

            61074LL*Power(r, 7LL)*Power(xj, 7LL) + 8586LL*Power(r, 8LL)*Power(xj, 8LL) +

            867LL*Power(r, 9LL)*Power(xj, 9LL) + 58LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) -

        Power(xi, 22LL)*(935550LL + 1715175LL*r*xj +

                         1559250LL*Power(r, 2LL)*Power(xj, 2LL) + 935550LL*Power(r, 3LL)*Power(xj, 3LL) +

                         415800LL*Power(r, 4LL)*Power(xj, 4LL) + 145530LL*Power(r, 5LL)*Power(xj, 5LL) +

                         41580LL*Power(r, 6LL)*Power(xj, 6LL) + 9900LL*Power(r, 7LL)*Power(xj, 7LL) +

                         1980LL*Power(r, 8LL)*Power(xj, 8LL) + 330LL*Power(r, 9LL)*Power(xj, 9LL) +

                         44LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

        11LL*Power(xi, 20LL)*Power(xj, 2LL)*

        (1105650LL + 2027025LL*r*xj + 1842750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1105650LL*Power(r, 3LL)*Power(xj, 3LL) + 491400LL*Power(r, 4LL)*Power(xj, 4LL) +

            171990LL*Power(r, 5LL)*Power(xj, 5LL) + 49140LL*Power(r, 6LL)*Power(xj, 6LL) +

            11700LL*Power(r, 7LL)*Power(xj, 7LL) + 2340LL*Power(r, 8LL)*Power(xj, 8LL) +

            390LL*Power(r, 9LL)*Power(xj, 9LL) + 52LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) -

        11LL*Power(xi, 2LL)*Power(xj, 20LL)*

        (-48648600LL + 2027025LL*r*xj + 44594550LL*Power(r, 2LL)*Power(xj, 2LL) +

            36486450LL*Power(r, 3LL)*Power(xj, 3LL) +

            16216200LL*Power(r, 4LL)*Power(xj, 4LL) +

            4864860LL*Power(r, 5LL)*Power(xj, 5LL) +

            1065960LL*Power(r, 6LL)*Power(xj, 6LL) + 176040LL*Power(r, 7LL)*Power(xj, 7LL) +

            21960LL*Power(r, 8LL)*Power(xj, 8LL) + 2010LL*Power(r, 9LL)*Power(xj, 9LL) +

            124LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

        Power(xj, 22LL)*(340540200LL + 468242775LL*r*xj +

                         312161850LL*Power(r, 2LL)*Power(xj, 2LL) +

                         133783650LL*Power(r, 3LL)*Power(xj, 3LL) +

                         41164200LL*Power(r, 4LL)*Power(xj, 4LL) +

                         9604980LL*Power(r, 5LL)*Power(xj, 5LL) +

                         1746360LL*Power(r, 6LL)*Power(xj, 6LL) + 249480LL*Power(r, 7LL)*Power(xj, 7LL) +

                         27720LL*Power(r, 8LL)*Power(xj, 8LL) + 2310LL*Power(r, 9LL)*Power(xj, 9LL) +

                         132LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

        165LL*Power(xi, 14LL)*Power(xj, 8LL)*

        (4054050LL + 7432425LL*r*xj + 6756750LL*Power(r, 2LL)*Power(xj, 2LL) +

            4054050LL*Power(r, 3LL)*Power(xj, 3LL) +

            1801800LL*Power(r, 4LL)*Power(xj, 4LL) + 631260LL*Power(r, 5LL)*Power(xj, 5LL) +

            178920LL*Power(r, 6LL)*Power(xj, 6LL) + 43176LL*Power(r, 7LL)*Power(xj, 7LL) +

            8904LL*Power(r, 8LL)*Power(xj, 8LL) + 1428LL*Power(r, 9LL)*Power(xj, 9LL) +

            152LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

        231LL*Power(xi, 12LL)*Power(xj, 10LL)*

        (5212350LL + 9555975LL*r*xj + 8687250LL*Power(r, 2LL)*Power(xj, 2LL) +

            5209650LL*Power(r, 3LL)*Power(xj, 3LL) +

            2327400LL*Power(r, 4LL)*Power(xj, 4LL) + 801540LL*Power(r, 5LL)*Power(xj, 5LL) +

            230040LL*Power(r, 6LL)*Power(xj, 6LL) + 57240LL*Power(r, 7LL)*Power(xj, 7LL) +

            11640LL*Power(r, 8LL)*Power(xj, 8LL) + 1740LL*Power(r, 9LL)*Power(xj, 9LL) +

            168LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) -

        231LL*Power(xi, 10LL)*Power(xj, 12LL)*

        (6949800LL + 12746025LL*r*xj + 11535750LL*Power(r, 2LL)*Power(xj, 2LL) +

            7056450LL*Power(r, 3LL)*Power(xj, 3LL) +

            3040200LL*Power(r, 4LL)*Power(xj, 4LL) +

            1051920LL*Power(r, 5LL)*Power(xj, 5LL) + 316800LL*Power(r, 6LL)*Power(xj, 6LL) +

            79680LL*Power(r, 7LL)*Power(xj, 7LL) + 15360LL*Power(r, 8LL)*Power(xj, 8LL) +

            2100LL*Power(r, 9LL)*Power(xj, 9LL) + 184LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

        165LL*Power(xi, 8LL)*Power(xj, 14LL)*

        (9775080LL + 17424855LL*r*xj + 17019450LL*Power(r, 2LL)*Power(xj, 2LL) +

            9519930LL*Power(r, 3LL)*Power(xj, 3LL) + 4059720LL*Power(r, 4LL)*Power(xj, 4LL) +

            1519056LL*Power(r, 5LL)*Power(xj, 5LL) + 475776LL*Power(r, 6LL)*Power(xj, 6LL) +

            114720LL*Power(r, 7LL)*Power(xj, 7LL) + 20256LL*Power(r, 8LL)*Power(xj, 8LL) +

            2508LL*Power(r, 9LL)*Power(xj, 9LL) + 200LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (935550LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 13LL)*Power(xi + xj, 13LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_6S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_1S_6S(r, xj, xi);
}
