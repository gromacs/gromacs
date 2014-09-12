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

cl_R DSlater_5S_5S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-(1e9*cl_float(223247882, precision)+3524011562LL)*xi + (1e9*cl_float(243290200, precision)+8176640000LL)*exp(2LL*r*xi)*xi -

                  (1e9*cl_float(406411127, precision)+7742766250LL)*r*Power(xi, 2LL) -

                  (1e9*cl_float(366777722, precision)+5928565750LL)*Power(r, 2LL)*Power(xi, 3LL) -

                  (1e9*cl_float(218697853, precision)+9397652000LL)*Power(r, 3LL)*Power(xi, 4LL) -

                  (1e9*cl_float(968809971, precision)+432303000LL)*Power(r, 4LL)*Power(xi, 5LL) -

                  (1e9*cl_float(339917275, precision)+693195200LL)*Power(r, 5LL)*Power(xi, 6LL) -

                  983239817883523200LL*Power(r, 6LL)*Power(xi, 7LL) -

                  240924879420825600LL*Power(r, 7LL)*Power(xi, 8LL) -

                  50973581199340800LL*Power(r, 8LL)*Power(xi, 9LL) -

                  9439831425024000LL*Power(r, 9LL)*Power(xi, 10LL) -

                  1544699687731200LL*Power(r, 10LL)*Power(xi, 11LL) -

                  224683590942720LL*Power(r, 11LL)*Power(xi, 12LL) -

                  29125650677760LL*Power(r, 12LL)*Power(xi, 13LL) -

                  3360652001280LL*Power(r, 13LL)*Power(xi, 14LL) -

                  342923673600LL*Power(r, 14LL)*Power(xi, 15LL) -

                  30482104320LL*Power(r, 15LL)*Power(xi, 16LL) -

                  2286157824LL*Power(r, 16LL)*Power(xi, 17LL) -

                  134479872LL*Power(r, 17LL)*Power(xi, 18LL) - 4980736LL*Power(r, 18LL)*Power(xi, 19LL))/

                (1.21645100408832e19*exp(2LL*r*xi)*r) +

                (-(1e9*cl_float(121645100, precision)+4088320000LL) + (1e9*cl_float(121645100, precision)+4088320000LL)*exp(2LL*r*xi) -

                 (1e9*cl_float(223247882, precision)+3524011562LL)*r*xi -

                 (1e9*cl_float(203205563, precision)+8871383125LL)*Power(r, 2LL)*Power(xi, 2LL) -

                 (1e9*cl_float(122259240, precision)+8642855250LL)*Power(r, 3LL)*Power(xi, 3LL) -

                 (1e9*cl_float(546744634, precision)+849413000LL)*Power(r, 4LL)*Power(xi, 4LL) -

                 (1e9*cl_float(193761994, precision)+286460600LL)*Power(r, 5LL)*Power(xi, 5LL) -

                 566528792821992000LL*Power(r, 6LL)*Power(xi, 6LL) -

                 140462831126217600LL*Power(r, 7LL)*Power(xi, 7LL) -

                 30115609927603200LL*Power(r, 8LL)*Power(xi, 8LL) -

                 5663731244371200LL*Power(r, 9LL)*Power(xi, 9LL) -

                 943983142502400LL*Power(r, 10LL)*Power(xi, 10LL) -

                 140427244339200LL*Power(r, 11LL)*Power(xi, 11LL) -

                 18723632578560LL*Power(r, 12LL)*Power(xi, 12LL) -

                 2240434667520LL*Power(r, 13LL)*Power(xi, 13LL) -

                 240046571520LL*Power(r, 14LL)*Power(xi, 14LL) -

                 22861578240LL*Power(r, 15LL)*Power(xi, 15LL) -

                 1905131520LL*Power(r, 16LL)*Power(xi, 16LL) -

                 134479872LL*Power(r, 17LL)*Power(xi, 17LL) -

                 7471104LL*Power(r, 18LL)*Power(xi, 18LL) - 262144LL*Power(r, 19LL)*Power(xi, 19LL))/

                (1.21645100408832e19*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-(1e9*cl_float(121645100, precision)+4088320000LL) + (1e9*cl_float(121645100, precision)+4088320000LL)*exp(2LL*r*xi) -

                     (1e9*cl_float(223247882, precision)+3524011562LL)*r*xi -

                     (1e9*cl_float(203205563, precision)+8871383125LL)*Power(r, 2LL)*Power(xi, 2LL) -

                     (1e9*cl_float(122259240, precision)+8642855250LL)*Power(r, 3LL)*Power(xi, 3LL) -

                     (1e9*cl_float(546744634, precision)+849413000LL)*Power(r, 4LL)*Power(xi, 4LL) -

                     (1e9*cl_float(193761994, precision)+286460600LL)*Power(r, 5LL)*Power(xi, 5LL) -

                     566528792821992000LL*Power(r, 6LL)*Power(xi, 6LL) -

                     140462831126217600LL*Power(r, 7LL)*Power(xi, 7LL) -

                     30115609927603200LL*Power(r, 8LL)*Power(xi, 8LL) -

                     5663731244371200LL*Power(r, 9LL)*Power(xi, 9LL) -

                     943983142502400LL*Power(r, 10LL)*Power(xi, 10LL) -

                     140427244339200LL*Power(r, 11LL)*Power(xi, 11LL) -

                     18723632578560LL*Power(r, 12LL)*Power(xi, 12LL) -

                     2240434667520LL*Power(r, 13LL)*Power(xi, 13LL) -

                     240046571520LL*Power(r, 14LL)*Power(xi, 14LL) -

                     22861578240LL*Power(r, 15LL)*Power(xi, 15LL) -

                     1905131520LL*Power(r, 16LL)*Power(xi, 16LL) -

                     134479872LL*Power(r, 17LL)*Power(xi, 17LL) -

                     7471104LL*Power(r, 18LL)*Power(xi, 18LL) - 262144LL*Power(r, 19LL)*Power(xi, 19LL)))/

                (6.0822550204416e18*exp(2LL*r*xi)*r)

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
            S = (70875LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 19LL) +

                 exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-630LL*Power(r, 8LL)*Power(xi, 34LL) - 10LL*Power(r, 9LL)*Power(xi, 35LL) +

                  70875LL*Power(xj, 26LL) + 127575LL*r*xi*Power(xj, 26LL) -

                  30LL*Power(r, 7LL)*Power(xi, 33LL)*(630LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  14175LL*Power(xi, 2LL)*Power(xj, 24LL)*(-95LL + 8LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  4725LL*r*Power(xi, 3LL)*Power(xj, 24LL)*

                  (-513LL + 14LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  90LL*Power(r, 6LL)*Power(xi, 32LL)*(3920LL + 43LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  4725LL*r*Power(xi, 5LL)*Power(xj, 22LL)*

                  (4617LL - 266LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  14175LL*Power(xi, 4LL)*Power(xj, 22LL)*

                  (855LL - 152LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  36LL*Power(r, 5LL)*Power(xi, 31LL)*

                  (-124950LL - 4985LL*Power(r, 2LL)*Power(xj, 2LL) +

            13LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  36LL*Power(r, 4LL)*Power(xi, 30LL)*

                  (-1124550LL - 127960LL*Power(r, 2LL)*Power(xj, 2LL) +

            863LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  135LL*r*Power(xi, 7LL)*Power(xj, 20LL)*

                  (-915705LL + 83790LL*Power(r, 2LL)*Power(xj, 2LL) -

            1330LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  315LL*Power(xi, 6LL)*Power(xj, 20LL)*

                  (-218025LL + 61560LL*Power(r, 2LL)*Power(xj, 2LL) -

            1710LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  36LL*Power(r, 3LL)*Power(xi, 29LL)*

                  (7122150LL + 2102730LL*Power(r, 2LL)*Power(xj, 2LL) -

            23294LL*Power(r, 4LL)*Power(xj, 4LL) + 37LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  36LL*Power(r, 2LL)*Power(xi, 28LL)*

                  (30523500LL + 23401350LL*Power(r, 2LL)*Power(xj, 2LL) -

            299250LL*Power(r, 4LL)*Power(xj, 4LL) + 1297LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  r*Power(xi, 17LL)*Power(xj, 10LL)*

                  (1073961177975LL - 21753487980LL*Power(r, 2LL)*Power(xj, 2LL) -

            745994340LL*Power(r, 4LL)*Power(xj, 4LL) +

            5307156LL*Power(r, 6LL)*Power(xj, 6LL) - 818LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  10LL*r*Power(xi, 9LL)*Power(xj, 18LL)*

                  (49448070LL - 6409935LL*Power(r, 2LL)*Power(xj, 2LL) +

            161595LL*Power(r, 4LL)*Power(xj, 4LL) - 1026LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) +

                  90LL*Power(xi, 8LL)*Power(xj, 18LL)*

                  (3052350LL - 1220940LL*Power(r, 2LL)*Power(xj, 2LL) +

            53865LL*Power(r, 4LL)*Power(xj, 4LL) - 532LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) -

                  1710LL*Power(xi, 10LL)*Power(xj, 16LL)*

                  (481950LL - 257040LL*Power(r, 2LL)*Power(xj, 2LL) +

            16065LL*Power(r, 4LL)*Power(xj, 4LL) - 252LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) +

                  6LL*r*Power(xi, 11LL)*Power(xj, 16LL)*

                  (-207559800LL + 50390550LL*Power(r, 2LL)*Power(xj, 2LL) -

            1165815LL*Power(r, 4LL)*Power(xj, 4LL) +

            21396LL*Power(r, 6LL)*Power(xj, 6LL) + 5LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  18LL*r*Power(xi, 13LL)*Power(xj, 14LL)*

                  (-1703720025LL - 155669850LL*Power(r, 2LL)*Power(xj, 2LL) -

            7410270LL*Power(r, 4LL)*Power(xj, 4LL) - 1532LL*Power(r, 6LL)*Power(xj, 6LL) +

            26LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  18LL*r*Power(xi, 15LL)*Power(xj, 12LL)*

                  (19380896325LL + 1329128850LL*Power(r, 2LL)*Power(xj, 2LL) -

            7608930LL*Power(r, 4LL)*Power(xj, 4LL) -

            116238LL*Power(r, 6LL)*Power(xj, 6LL) + 74LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  18LL*Power(xi, 12LL)*Power(xj, 14LL)*

                  (89026875LL + 179071200LL*Power(r, 2LL)*Power(xj, 2LL) +

            1552950LL*Power(r, 4LL)*Power(xj, 4LL) +

            295820LL*Power(r, 6LL)*Power(xj, 6LL) + 146LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  18LL*r*Power(xi, 25LL)*Power(xj, 2LL)*

                  (-5449970925LL - 1137574935LL*Power(r, 2LL)*Power(xj, 2LL) +

            37834755LL*Power(r, 4LL)*Power(xj, 4LL) -

            273062LL*Power(r, 6LL)*Power(xj, 6LL) + 171LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  9LL*r*Power(xi, 19LL)*Power(xj, 8LL)*

                  (-37914907275LL + 7613889570LL*Power(r, 2LL)*Power(xj, 2LL) -

            170524620LL*Power(r, 4LL)*Power(xj, 4LL) +

            397936LL*Power(r, 6LL)*Power(xj, 6LL) + 342LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  3LL*r*Power(xi, 23LL)*Power(xj, 4LL)*

                  (219130630425LL - 11118046590LL*Power(r, 2LL)*Power(xj, 2LL) +

            327611970LL*Power(r, 4LL)*Power(xj, 4LL) -

            2920908LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  3LL*r*Power(xi, 21LL)*Power(xj, 6LL)*

                  (-345162539925LL + 19030764690LL*Power(r, 2LL)*Power(xj, 2LL) -

            141976170LL*Power(r, 4LL)*Power(xj, 4LL) -

            1441872LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  63LL*Power(xi, 20LL)*Power(xj, 6LL)*

                  (-50980542525LL + 6240202920LL*Power(r, 2LL)*Power(xj, 2LL) -

            201314310LL*Power(r, 4LL)*Power(xj, 4LL) +

            956080LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  18LL*Power(xi, 14LL)*Power(xj, 12LL)*

                  (-7803332775LL - 2519206200LL*Power(r, 2LL)*Power(xj, 2LL) -

            119719950LL*Power(r, 4LL)*Power(xj, 4LL) +

            182280LL*Power(r, 6LL)*Power(xj, 6LL) + 2734LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  18LL*Power(xi, 26LL)*(195859125LL + 1794781800LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        67337235LL*Power(r, 4LL)*Power(xj, 4LL) -

                                        1659700LL*Power(r, 6LL)*Power(xj, 6LL) + 4089LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  9LL*Power(xi, 18LL)*Power(xj, 8LL)*

                  (-357591274425LL + 8328390840LL*Power(r, 2LL)*Power(xj, 2LL) +

            912042180LL*Power(r, 4LL)*Power(xj, 4LL) -

            12842480LL*Power(r, 6LL)*Power(xj, 6LL) + 10678LL*Power(r, 8LL)*Power(xj, 8LL))

                  - 9LL*Power(xi, 16LL)*Power(xj, 10LL)*

                  (128599724925LL + 21298077360LL*Power(r, 2LL)*Power(xj, 2LL) -

            267928500LL*Power(r, 4LL)*Power(xj, 4LL) -

            5458320LL*Power(r, 6LL)*Power(xj, 6LL) + 14722LL*Power(r, 8LL)*Power(xj, 8LL))

                  + 18LL*Power(xi, 24LL)*Power(xj, 2LL)*

                  (-7604930025LL - 8866107180LL*Power(r, 2LL)*Power(xj, 2LL) +

            399272265LL*Power(r, 4LL)*Power(xj, 4LL) -

            5925780LL*Power(r, 6LL)*Power(xj, 6LL) + 17651LL*Power(r, 8LL)*Power(xj, 8LL))

                  - 9LL*Power(xi, 22LL)*Power(xj, 4LL)*(129194933175LL +

                                                        3909863160LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        91420770LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                        8762040LL*Power(r, 6LL)*Power(xj, 6LL) + 43928LL*Power(r, 8LL)*Power(xj, 8LL))

                  + Power(xi, 27LL)*(-2884470750LL*r - 6409935000LL*Power(r, 3LL)*Power(xj, 2LL) +

                                     28332990LL*Power(r, 5LL)*Power(xj, 4LL) +

                                     58104LL*Power(r, 7LL)*Power(xj, 6LL) + 818LL*Power(r, 9LL)*Power(xj, 8LL))) +

                 exp(2LL*r*xi)*Power(xi, 12LL)*

                 (Power(xi, 8LL)*Power(xj, 18LL)*

                  (3218321469825LL - 341234165475LL*r*xj -

            393132783960LL*Power(r, 2LL)*Power(xj, 2LL) -

            57092294070LL*Power(r, 3LL)*Power(xj, 3LL) +

            822786930LL*Power(r, 4LL)*Power(xj, 4LL) +

            982835910LL*Power(r, 5LL)*Power(xj, 5LL) +

            106664040LL*Power(r, 6LL)*Power(xj, 6LL) +

            4915116LL*Power(r, 7LL)*Power(xj, 7LL) + 73602LL*Power(r, 8LL)*Power(xj, 8LL) -

            818LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  10LL*Power(xj, 26LL)*(352546425LL + 288447075LL*r*xj +

                                        109884600LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        25639740LL*Power(r, 3LL)*Power(xj, 3LL) +

                                        4048380LL*Power(r, 4LL)*Power(xj, 4LL) +

                                        449820LL*Power(r, 5LL)*Power(xj, 5LL) + 35280LL*Power(r, 6LL)*Power(xj, 6LL) +

                                        1890LL*Power(r, 7LL)*Power(xj, 7LL) + 63LL*Power(r, 8LL)*Power(xj, 8LL) +

                                        Power(r, 9LL)*Power(xj, 9LL)) +

                  30LL*Power(xi, 2LL)*Power(xj, 24LL)*

                  (4562958015LL + 3269982555LL*r*xj +

            1076869080LL*Power(r, 2LL)*Power(xj, 2LL) +

            213664500LL*Power(r, 3LL)*Power(xj, 3LL) +

            28081620LL*Power(r, 4LL)*Power(xj, 4LL) +

            2523276LL*Power(r, 5LL)*Power(xj, 5LL) +

            153552LL*Power(r, 6LL)*Power(xj, 6LL) + 5982LL*Power(r, 7LL)*Power(xj, 7LL) +

            129LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) -

                  15LL*Power(xi, 24LL)*Power(xj, 2LL)*

                  (-89775LL - 161595LL*r*xj - 143640LL*Power(r, 2LL)*Power(xj, 2LL) -

            83790LL*Power(r, 3LL)*Power(xj, 3LL) - 35910LL*Power(r, 4LL)*Power(xj, 4LL) -

            11970LL*Power(r, 5LL)*Power(xj, 5LL) - 3192LL*Power(r, 6LL)*Power(xj, 6LL) -

            684LL*Power(r, 7LL)*Power(xj, 7LL) - 114LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  5LL*Power(xi, 26LL)*(14175LL + 25515LL*r*xj +

                                       22680LL*Power(r, 2LL)*Power(xj, 2LL) + 13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       5670LL*Power(r, 4LL)*Power(xj, 4LL) + 1890LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       504LL*Power(r, 6LL)*Power(xj, 6LL) + 108LL*Power(r, 7LL)*Power(xj, 7LL) +

                                       18LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  1938LL*Power(xi, 14LL)*Power(xj, 12LL)*

                  (-826875LL + 15824025LL*r*xj - 23398200LL*Power(r, 2LL)*Power(xj, 2LL) +

            12344850LL*Power(r, 3LL)*Power(xj, 3LL) +

            1244250LL*Power(r, 4LL)*Power(xj, 4LL) -

            384930LL*Power(r, 5LL)*Power(xj, 5LL) - 59640LL*Power(r, 6LL)*Power(xj, 6LL) -

            1848LL*Power(r, 7LL)*Power(xj, 7LL) + 84LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  1938LL*Power(xi, 12LL)*Power(xj, 14LL)*

                  (72476775LL - 180008325LL*r*xj + 98907480LL*Power(r, 2LL)*Power(xj, 2LL) +

            11224710LL*Power(r, 3LL)*Power(xj, 3LL) -

            4235490LL*Power(r, 4LL)*Power(xj, 4LL) -

            791910LL*Power(r, 5LL)*Power(xj, 5LL) - 31080LL*Power(r, 6LL)*Power(xj, 6LL) +

            2232LL*Power(r, 7LL)*Power(xj, 7LL) + 204LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  342LL*Power(xi, 16LL)*Power(xj, 10LL)*

                  (2409750LL + 3641400LL*r*xj + 9424800LL*Power(r, 2LL)*Power(xj, 2LL) -

            8193150LL*Power(r, 3LL)*Power(xj, 3LL) +

            6301050LL*Power(r, 4LL)*Power(xj, 4LL) +

            400470LL*Power(r, 5LL)*Power(xj, 5LL) - 143640LL*Power(r, 6LL)*Power(xj, 6LL) -

            15518LL*Power(r, 7LL)*Power(xj, 7LL) - 281LL*Power(r, 8LL)*Power(xj, 8LL) +

            9LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  171LL*Power(xi, 10LL)*Power(xj, 16LL)*

                  (-6768406575LL + 6280474725LL*r*xj +

            438336360LL*Power(r, 2LL)*Power(xj, 2LL) -

            400731030LL*Power(r, 3LL)*Power(xj, 3LL) -

            74168430LL*Power(r, 4LL)*Power(xj, 4LL) -

            2490810LL*Power(r, 5LL)*Power(xj, 5LL) +

            461160LL*Power(r, 6LL)*Power(xj, 6LL) + 51244LL*Power(r, 7LL)*Power(xj, 7LL) +

            1858LL*Power(r, 8LL)*Power(xj, 8LL) + 18LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  9LL*Power(xi, 22LL)*Power(xj, 4LL)*

                  (-1346625LL - 2423925LL*r*xj - 2154600LL*Power(r, 2LL)*Power(xj, 2LL) -

            1256850LL*Power(r, 3LL)*Power(xj, 3LL) -

            538650LL*Power(r, 4LL)*Power(xj, 4LL) - 179550LL*Power(r, 5LL)*Power(xj, 5LL) -

            47880LL*Power(r, 6LL)*Power(xj, 6LL) - 14264LL*Power(r, 7LL)*Power(xj, 7LL) +

            292LL*Power(r, 8LL)*Power(xj, 8LL) + 52LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  9LL*Power(xi, 4LL)*Power(xj, 22LL)*

                  (-129194933175LL - 73043543475LL*r*xj -

            17732214360LL*Power(r, 2LL)*Power(xj, 2LL) -

            2275149870LL*Power(r, 3LL)*Power(xj, 3LL) -

            134674470LL*Power(r, 4LL)*Power(xj, 4LL) +

            3148110LL*Power(r, 5LL)*Power(xj, 5LL) +

            1197000LL*Power(r, 6LL)*Power(xj, 6LL) + 93176LL*Power(r, 7LL)*Power(xj, 7LL) +

            3452LL*Power(r, 8LL)*Power(xj, 8LL) + 52LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  9LL*Power(xi, 6LL)*Power(xj, 20LL)*

                  (356863797675LL + 115054179975LL*r*xj +

            3909863160LL*Power(r, 2LL)*Power(xj, 2LL) -

            3706015530LL*Power(r, 3LL)*Power(xj, 3LL) -

            798544530LL*Power(r, 4LL)*Power(xj, 4LL) -

            75669510LL*Power(r, 5LL)*Power(xj, 5LL) -

            3319400LL*Power(r, 6LL)*Power(xj, 6LL) - 6456LL*Power(r, 7LL)*Power(xj, 7LL) +

            5188LL*Power(r, 8LL)*Power(xj, 8LL) + 148LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  9LL*Power(xi, 20LL)*Power(xj, 6LL)*

                  (-7630875LL - 13735575LL*r*xj - 12209400LL*Power(r, 2LL)*Power(xj, 2LL) -

            7122150LL*Power(r, 3LL)*Power(xj, 3LL) -

            3052350LL*Power(r, 4LL)*Power(xj, 4LL) -

            777210LL*Power(r, 5LL)*Power(xj, 5LL) - 591640LL*Power(r, 6LL)*Power(xj, 6LL) +

            3064LL*Power(r, 7LL)*Power(xj, 7LL) + 5468LL*Power(r, 8LL)*Power(xj, 8LL) +

            148LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  2LL*Power(xi, 18LL)*Power(xj, 8LL)*

                  (-137355750LL - 247240350LL*r*xj - 219769200LL*Power(r, 2LL)*Power(xj, 2LL) -

            151171650LL*Power(r, 3LL)*Power(xj, 3LL) +

            13976550LL*Power(r, 4LL)*Power(xj, 4LL) -

            66692430LL*Power(r, 5LL)*Power(xj, 5LL) -

            1640520LL*Power(r, 6LL)*Power(xj, 6LL) +

            1046142LL*Power(r, 7LL)*Power(xj, 7LL) + 66249LL*Power(r, 8LL)*Power(xj, 8LL) +

            409LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (70875LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 19LL)*

                 Power(xi + xj, 19LL)) + (2LL*(70875LL*exp(2LL*r*(xi + xj))*

                                               Power(Power(xi, 2LL) - Power(xj, 2LL), 19LL) +

                                               exp(2LL*r*xj)*Power(xj, 12LL)*

                                               (-630LL*Power(r, 8LL)*Power(xi, 34LL) - 10LL*Power(r, 9LL)*Power(xi, 35LL) +

                                       70875LL*Power(xj, 26LL) + 127575LL*r*xi*Power(xj, 26LL) -

                                       30LL*Power(r, 7LL)*Power(xi, 33LL)*(630LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                       14175LL*Power(xi, 2LL)*Power(xj, 24LL)*

                                       (-95LL + 8LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                       4725LL*r*Power(xi, 3LL)*Power(xj, 24LL)*

                                       (-513LL + 14LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                       90LL*Power(r, 6LL)*Power(xi, 32LL)*(3920LL + 43LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                       4725LL*r*Power(xi, 5LL)*Power(xj, 22LL)*

                                       (4617LL - 266LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

                                       + 14175LL*Power(xi, 4LL)*Power(xj, 22LL)*

                                       (855LL - 152LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       36LL*Power(r, 5LL)*Power(xi, 31LL)*

                                       (-124950LL - 4985LL*Power(r, 2LL)*Power(xj, 2LL) +

              13LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       36LL*Power(r, 4LL)*Power(xi, 30LL)*

                                       (-1124550LL - 127960LL*Power(r, 2LL)*Power(xj, 2LL) +

              863LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       135LL*r*Power(xi, 7LL)*Power(xj, 20LL)*

                                       (-915705LL + 83790LL*Power(r, 2LL)*Power(xj, 2LL) -

              1330LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                       315LL*Power(xi, 6LL)*Power(xj, 20LL)*

                                       (-218025LL + 61560LL*Power(r, 2LL)*Power(xj, 2LL) -

              1710LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                       36LL*Power(r, 3LL)*Power(xi, 29LL)*

                                       (7122150LL + 2102730LL*Power(r, 2LL)*Power(xj, 2LL) -

              23294LL*Power(r, 4LL)*Power(xj, 4LL) + 37LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                       36LL*Power(r, 2LL)*Power(xi, 28LL)*

                                       (30523500LL + 23401350LL*Power(r, 2LL)*Power(xj, 2LL) -

              299250LL*Power(r, 4LL)*Power(xj, 4LL) + 1297LL*Power(r, 6LL)*Power(xj, 6LL))

                                       + r*Power(xi, 17LL)*Power(xj, 10LL)*

                                       (1073961177975LL - 21753487980LL*Power(r, 2LL)*Power(xj, 2LL) -

              745994340LL*Power(r, 4LL)*Power(xj, 4LL) +

              5307156LL*Power(r, 6LL)*Power(xj, 6LL) - 818LL*Power(r, 8LL)*Power(xj, 8LL))

                                       + 10LL*r*Power(xi, 9LL)*Power(xj, 18LL)*

                                       (49448070LL - 6409935LL*Power(r, 2LL)*Power(xj, 2LL) +

              161595LL*Power(r, 4LL)*Power(xj, 4LL) -

              1026LL*Power(r, 6LL)*Power(xj, 6LL) + Power(r, 8LL)*Power(xj, 8LL)) +

                                       90LL*Power(xi, 8LL)*Power(xj, 18LL)*

                                       (3052350LL - 1220940LL*Power(r, 2LL)*Power(xj, 2LL) +

              53865LL*Power(r, 4LL)*Power(xj, 4LL) - 532LL*Power(r, 6LL)*Power(xj, 6LL) +

              Power(r, 8LL)*Power(xj, 8LL)) -

                                       1710LL*Power(xi, 10LL)*Power(xj, 16LL)*

                                       (481950LL - 257040LL*Power(r, 2LL)*Power(xj, 2LL) +

              16065LL*Power(r, 4LL)*Power(xj, 4LL) - 252LL*Power(r, 6LL)*Power(xj, 6LL) +

              Power(r, 8LL)*Power(xj, 8LL)) +

                                       6LL*r*Power(xi, 11LL)*Power(xj, 16LL)*

                                       (-207559800LL + 50390550LL*Power(r, 2LL)*Power(xj, 2LL) -

              1165815LL*Power(r, 4LL)*Power(xj, 4LL) +

              21396LL*Power(r, 6LL)*Power(xj, 6LL) + 5LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                       18LL*r*Power(xi, 13LL)*Power(xj, 14LL)*

                                       (-1703720025LL - 155669850LL*Power(r, 2LL)*Power(xj, 2LL) -

              7410270LL*Power(r, 4LL)*Power(xj, 4LL) -

              1532LL*Power(r, 6LL)*Power(xj, 6LL) + 26LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                       18LL*r*Power(xi, 15LL)*Power(xj, 12LL)*

                                       (19380896325LL + 1329128850LL*Power(r, 2LL)*Power(xj, 2LL) -

              7608930LL*Power(r, 4LL)*Power(xj, 4LL) -

              116238LL*Power(r, 6LL)*Power(xj, 6LL) + 74LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                       18LL*Power(xi, 12LL)*Power(xj, 14LL)*

                                       (89026875LL + 179071200LL*Power(r, 2LL)*Power(xj, 2LL) +

              1552950LL*Power(r, 4LL)*Power(xj, 4LL) +

              295820LL*Power(r, 6LL)*Power(xj, 6LL) + 146LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                       18LL*r*Power(xi, 25LL)*Power(xj, 2LL)*

                                       (-5449970925LL - 1137574935LL*Power(r, 2LL)*Power(xj, 2LL) +

              37834755LL*Power(r, 4LL)*Power(xj, 4LL) -

              273062LL*Power(r, 6LL)*Power(xj, 6LL) + 171LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                       9LL*r*Power(xi, 19LL)*Power(xj, 8LL)*

                                       (-37914907275LL + 7613889570LL*Power(r, 2LL)*Power(xj, 2LL) -

              170524620LL*Power(r, 4LL)*Power(xj, 4LL) +

              397936LL*Power(r, 6LL)*Power(xj, 6LL) + 342LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                       3LL*r*Power(xi, 23LL)*Power(xj, 4LL)*

                                       (219130630425LL - 11118046590LL*Power(r, 2LL)*Power(xj, 2LL) +

              327611970LL*Power(r, 4LL)*Power(xj, 4LL) -

              2920908LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL))

                                       + 3LL*r*Power(xi, 21LL)*Power(xj, 6LL)*

                                       (-345162539925LL + 19030764690LL*Power(r, 2LL)*Power(xj, 2LL) -

              141976170LL*Power(r, 4LL)*Power(xj, 4LL) -

              1441872LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL))

                                       + 63LL*Power(xi, 20LL)*Power(xj, 6LL)*

                                       (-50980542525LL + 6240202920LL*Power(r, 2LL)*Power(xj, 2LL) -

              201314310LL*Power(r, 4LL)*Power(xj, 4LL) +

              956080LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL))

                                       + 18LL*Power(xi, 14LL)*Power(xj, 12LL)*

                                       (-7803332775LL - 2519206200LL*Power(r, 2LL)*Power(xj, 2LL) -

              119719950LL*Power(r, 4LL)*Power(xj, 4LL) +

              182280LL*Power(r, 6LL)*Power(xj, 6LL) + 2734LL*Power(r, 8LL)*Power(xj, 8LL))

                                       - 18LL*Power(xi, 26LL)*(195859125LL + 1794781800LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                               67337235LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                               1659700LL*Power(r, 6LL)*Power(xj, 6LL) + 4089LL*Power(r, 8LL)*Power(xj, 8LL))

                                       + 9LL*Power(xi, 18LL)*Power(xj, 8LL)*(-357591274425LL +

                                                                             8328390840LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                             912042180LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                             12842480LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                             10678LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                       9LL*Power(xi, 16LL)*Power(xj, 10LL)*

                                       (128599724925LL + 21298077360LL*Power(r, 2LL)*Power(xj, 2LL) -

              267928500LL*Power(r, 4LL)*Power(xj, 4LL) -

              5458320LL*Power(r, 6LL)*Power(xj, 6LL) + 14722LL*Power(r, 8LL)*Power(xj, 8LL)

                                       ) + 18LL*Power(xi, 24LL)*Power(xj, 2LL)*

                                       (-7604930025LL - 8866107180LL*Power(r, 2LL)*Power(xj, 2LL) +

              399272265LL*Power(r, 4LL)*Power(xj, 4LL) -

              5925780LL*Power(r, 6LL)*Power(xj, 6LL) + 17651LL*Power(r, 8LL)*Power(xj, 8LL)

                                       ) - 9LL*Power(xi, 22LL)*Power(xj, 4LL)*

                                       (129194933175LL + 3909863160LL*Power(r, 2LL)*Power(xj, 2LL) +

              91420770LL*Power(r, 4LL)*Power(xj, 4LL) -

              8762040LL*Power(r, 6LL)*Power(xj, 6LL) + 43928LL*Power(r, 8LL)*Power(xj, 8LL)

                                       ) + Power(xi, 27LL)*(-2884470750LL*r - 6409935000LL*Power(r, 3LL)*Power(xj, 2LL) +

                                                            28332990LL*Power(r, 5LL)*Power(xj, 4LL) +

                                                            58104LL*Power(r, 7LL)*Power(xj, 6LL) + 818LL*Power(r, 9LL)*Power(xj, 8LL))) +

                                               exp(2LL*r*xi)*Power(xi, 12LL)*

                                               (Power(xi, 8LL)*Power(xj, 18LL)*

                                       (3218321469825LL - 341234165475LL*r*xj -

              393132783960LL*Power(r, 2LL)*Power(xj, 2LL) -

              57092294070LL*Power(r, 3LL)*Power(xj, 3LL) +

              822786930LL*Power(r, 4LL)*Power(xj, 4LL) +

              982835910LL*Power(r, 5LL)*Power(xj, 5LL) +

              106664040LL*Power(r, 6LL)*Power(xj, 6LL) +

              4915116LL*Power(r, 7LL)*Power(xj, 7LL) +

              73602LL*Power(r, 8LL)*Power(xj, 8LL) - 818LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       10LL*Power(xj, 26LL)*(352546425LL + 288447075LL*r*xj +

                                                             109884600LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                             25639740LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                             4048380LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                             449820LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                             35280LL*Power(r, 6LL)*Power(xj, 6LL) + 1890LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                             63LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) +

                                       30LL*Power(xi, 2LL)*Power(xj, 24LL)*

                                       (4562958015LL + 3269982555LL*r*xj +

              1076869080LL*Power(r, 2LL)*Power(xj, 2LL) +

              213664500LL*Power(r, 3LL)*Power(xj, 3LL) +

              28081620LL*Power(r, 4LL)*Power(xj, 4LL) +

              2523276LL*Power(r, 5LL)*Power(xj, 5LL) +

              153552LL*Power(r, 6LL)*Power(xj, 6LL) + 5982LL*Power(r, 7LL)*Power(xj, 7LL) +

              129LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) -

                                       15LL*Power(xi, 24LL)*Power(xj, 2LL)*

                                       (-89775LL - 161595LL*r*xj - 143640LL*Power(r, 2LL)*Power(xj, 2LL) -

              83790LL*Power(r, 3LL)*Power(xj, 3LL) - 35910LL*Power(r, 4LL)*Power(xj, 4LL) -

              11970LL*Power(r, 5LL)*Power(xj, 5LL) - 3192LL*Power(r, 6LL)*Power(xj, 6LL) -

              684LL*Power(r, 7LL)*Power(xj, 7LL) - 114LL*Power(r, 8LL)*Power(xj, 8LL) +

              2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       5LL*Power(xi, 26LL)*(14175LL + 25515LL*r*xj +

                                                            22680LL*Power(r, 2LL)*Power(xj, 2LL) + 13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                            5670LL*Power(r, 4LL)*Power(xj, 4LL) + 1890LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                            504LL*Power(r, 6LL)*Power(xj, 6LL) + 108LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                            18LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       1938LL*Power(xi, 14LL)*Power(xj, 12LL)*

                                       (-826875LL + 15824025LL*r*xj - 23398200LL*Power(r, 2LL)*Power(xj, 2LL) +

              12344850LL*Power(r, 3LL)*Power(xj, 3LL) +

              1244250LL*Power(r, 4LL)*Power(xj, 4LL) -

              384930LL*Power(r, 5LL)*Power(xj, 5LL) -

              59640LL*Power(r, 6LL)*Power(xj, 6LL) - 1848LL*Power(r, 7LL)*Power(xj, 7LL) +

              84LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       1938LL*Power(xi, 12LL)*Power(xj, 14LL)*

                                       (72476775LL - 180008325LL*r*xj + 98907480LL*Power(r, 2LL)*Power(xj, 2LL) +

              11224710LL*Power(r, 3LL)*Power(xj, 3LL) -

              4235490LL*Power(r, 4LL)*Power(xj, 4LL) -

              791910LL*Power(r, 5LL)*Power(xj, 5LL) -

              31080LL*Power(r, 6LL)*Power(xj, 6LL) + 2232LL*Power(r, 7LL)*Power(xj, 7LL) +

              204LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       342LL*Power(xi, 16LL)*Power(xj, 10LL)*

                                       (2409750LL + 3641400LL*r*xj + 9424800LL*Power(r, 2LL)*Power(xj, 2LL) -

              8193150LL*Power(r, 3LL)*Power(xj, 3LL) +

              6301050LL*Power(r, 4LL)*Power(xj, 4LL) +

              400470LL*Power(r, 5LL)*Power(xj, 5LL) -

              143640LL*Power(r, 6LL)*Power(xj, 6LL) -

              15518LL*Power(r, 7LL)*Power(xj, 7LL) - 281LL*Power(r, 8LL)*Power(xj, 8LL) +

              9LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       171LL*Power(xi, 10LL)*Power(xj, 16LL)*

                                       (-6768406575LL + 6280474725LL*r*xj +

              438336360LL*Power(r, 2LL)*Power(xj, 2LL) -

              400731030LL*Power(r, 3LL)*Power(xj, 3LL) -

              74168430LL*Power(r, 4LL)*Power(xj, 4LL) -

              2490810LL*Power(r, 5LL)*Power(xj, 5LL) +

              461160LL*Power(r, 6LL)*Power(xj, 6LL) +

              51244LL*Power(r, 7LL)*Power(xj, 7LL) + 1858LL*Power(r, 8LL)*Power(xj, 8LL) +

              18LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       9LL*Power(xi, 22LL)*Power(xj, 4LL)*

                                       (-1346625LL - 2423925LL*r*xj - 2154600LL*Power(r, 2LL)*Power(xj, 2LL) -

              1256850LL*Power(r, 3LL)*Power(xj, 3LL) -

              538650LL*Power(r, 4LL)*Power(xj, 4LL) -

              179550LL*Power(r, 5LL)*Power(xj, 5LL) -

              47880LL*Power(r, 6LL)*Power(xj, 6LL) - 14264LL*Power(r, 7LL)*Power(xj, 7LL) +

              292LL*Power(r, 8LL)*Power(xj, 8LL) + 52LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       9LL*Power(xi, 4LL)*Power(xj, 22LL)*

                                       (-129194933175LL - 73043543475LL*r*xj -

              17732214360LL*Power(r, 2LL)*Power(xj, 2LL) -

              2275149870LL*Power(r, 3LL)*Power(xj, 3LL) -

              134674470LL*Power(r, 4LL)*Power(xj, 4LL) +

              3148110LL*Power(r, 5LL)*Power(xj, 5LL) +

              1197000LL*Power(r, 6LL)*Power(xj, 6LL) +

              93176LL*Power(r, 7LL)*Power(xj, 7LL) + 3452LL*Power(r, 8LL)*Power(xj, 8LL) +

              52LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       9LL*Power(xi, 6LL)*Power(xj, 20LL)*

                                       (356863797675LL + 115054179975LL*r*xj +

              3909863160LL*Power(r, 2LL)*Power(xj, 2LL) -

              3706015530LL*Power(r, 3LL)*Power(xj, 3LL) -

              798544530LL*Power(r, 4LL)*Power(xj, 4LL) -

              75669510LL*Power(r, 5LL)*Power(xj, 5LL) -

              3319400LL*Power(r, 6LL)*Power(xj, 6LL) -

              6456LL*Power(r, 7LL)*Power(xj, 7LL) + 5188LL*Power(r, 8LL)*Power(xj, 8LL) +

              148LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       9LL*Power(xi, 20LL)*Power(xj, 6LL)*

                                       (-7630875LL - 13735575LL*r*xj - 12209400LL*Power(r, 2LL)*Power(xj, 2LL) -

              7122150LL*Power(r, 3LL)*Power(xj, 3LL) -

              3052350LL*Power(r, 4LL)*Power(xj, 4LL) -

              777210LL*Power(r, 5LL)*Power(xj, 5LL) -

              591640LL*Power(r, 6LL)*Power(xj, 6LL) + 3064LL*Power(r, 7LL)*Power(xj, 7LL) +

              5468LL*Power(r, 8LL)*Power(xj, 8LL) + 148LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       2LL*Power(xi, 18LL)*Power(xj, 8LL)*

                                       (-137355750LL - 247240350LL*r*xj -

              219769200LL*Power(r, 2LL)*Power(xj, 2LL) -

              151171650LL*Power(r, 3LL)*Power(xj, 3LL) +

              13976550LL*Power(r, 4LL)*Power(xj, 4LL) -

              66692430LL*Power(r, 5LL)*Power(xj, 5LL) -

              1640520LL*Power(r, 6LL)*Power(xj, 6LL) +

              1046142LL*Power(r, 7LL)*Power(xj, 7LL) +

              66249LL*Power(r, 8LL)*Power(xj, 8LL) + 409LL*Power(r, 9LL)*Power(xj, 9LL)))))/

                (70875LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 19LL)*Power(xi + xj, 18LL)) -

                (141750LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 19LL) +

                 exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-5040LL*Power(r, 7LL)*Power(xi, 34LL) - 90LL*Power(r, 8LL)*Power(xi, 35LL) -

        7740LL*Power(r, 7LL)*Power(xi, 32LL)*Power(xj, 2LL) -

        60LL*Power(r, 8LL)*Power(xi, 33LL)*Power(xj, 2LL) + 127575LL*xi*Power(xj, 26LL) +

        226800LL*r*Power(xi, 2LL)*Power(xj, 26LL) +

        132300LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 26LL) -

        210LL*Power(r, 6LL)*Power(xi, 33LL)*(630LL + Power(r, 2LL)*Power(xj, 2LL)) +

        4725LL*Power(xi, 3LL)*Power(xj, 24LL)*(-513LL + 14LL*Power(r, 2LL)*Power(xj, 2LL)) -

        540LL*Power(r, 5LL)*Power(xi, 32LL)*(3920LL + 43LL*Power(r, 2LL)*Power(xj, 2LL)) +

        4725LL*r*Power(xi, 5LL)*Power(xj, 22LL)*

        (-532LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        14175LL*Power(xi, 4LL)*Power(xj, 22LL)*

        (-304LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        36LL*Power(r, 5LL)*Power(xi, 31LL)*

        (-9970LL*r*Power(xj, 2LL) + 52LL*Power(r, 3LL)*Power(xj, 4LL)) +

        36LL*Power(r, 4LL)*Power(xi, 30LL)*

        (-255920LL*r*Power(xj, 2LL) + 3452LL*Power(r, 3LL)*Power(xj, 4LL)) +

        4725LL*Power(xi, 5LL)*Power(xj, 22LL)*

        (4617LL - 266LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        180LL*Power(r, 4LL)*Power(xi, 31LL)*

        (-124950LL - 4985LL*Power(r, 2LL)*Power(xj, 2LL) +

            13LL*Power(r, 4LL)*Power(xj, 4LL)) +

        144LL*Power(r, 3LL)*Power(xi, 30LL)*

        (-1124550LL - 127960LL*Power(r, 2LL)*Power(xj, 2LL) +

            863LL*Power(r, 4LL)*Power(xj, 4LL)) +

        135LL*r*Power(xi, 7LL)*Power(xj, 20LL)*

        (167580LL*r*Power(xj, 2LL) - 5320LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) +

        315LL*Power(xi, 6LL)*Power(xj, 20LL)*

        (123120LL*r*Power(xj, 2LL) - 6840LL*Power(r, 3LL)*Power(xj, 4LL) +

            48LL*Power(r, 5LL)*Power(xj, 6LL)) -

        36LL*Power(r, 3LL)*Power(xi, 29LL)*

        (4205460LL*r*Power(xj, 2LL) - 93176LL*Power(r, 3LL)*Power(xj, 4LL) +

            222LL*Power(r, 5LL)*Power(xj, 6LL)) -

        36LL*Power(r, 2LL)*Power(xi, 28LL)*

        (46802700LL*r*Power(xj, 2LL) - 1197000LL*Power(r, 3LL)*Power(xj, 4LL) +

            7782LL*Power(r, 5LL)*Power(xj, 6LL)) +

        135LL*Power(xi, 7LL)*Power(xj, 20LL)*

        (-915705LL + 83790LL*Power(r, 2LL)*Power(xj, 2LL) -

            1330LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

        108LL*Power(r, 2LL)*Power(xi, 29LL)*

        (7122150LL + 2102730LL*Power(r, 2LL)*Power(xj, 2LL) -

            23294LL*Power(r, 4LL)*Power(xj, 4LL) + 37LL*Power(r, 6LL)*Power(xj, 6LL)) -

        72LL*r*Power(xi, 28LL)*(30523500LL + 23401350LL*Power(r, 2LL)*Power(xj, 2LL) -

                                299250LL*Power(r, 4LL)*Power(xj, 4LL) + 1297LL*Power(r, 6LL)*Power(xj, 6LL)) +

        r*Power(xi, 17LL)*Power(xj, 10LL)*

        (-43506975960LL*r*Power(xj, 2LL) - 2983977360LL*Power(r, 3LL)*Power(xj, 4LL) +

            31842936LL*Power(r, 5LL)*Power(xj, 6LL) - 6544LL*Power(r, 7LL)*Power(xj, 8LL)) +

        10LL*r*Power(xi, 9LL)*Power(xj, 18LL)*

        (-12819870LL*r*Power(xj, 2LL) + 646380LL*Power(r, 3LL)*Power(xj, 4LL) -

            6156LL*Power(r, 5LL)*Power(xj, 6LL) + 8LL*Power(r, 7LL)*Power(xj, 8LL)) +

        90LL*Power(xi, 8LL)*Power(xj, 18LL)*

        (-2441880LL*r*Power(xj, 2LL) + 215460LL*Power(r, 3LL)*Power(xj, 4LL) -

            3192LL*Power(r, 5LL)*Power(xj, 6LL) + 8LL*Power(r, 7LL)*Power(xj, 8LL)) -

        1710LL*Power(xi, 10LL)*Power(xj, 16LL)*

        (-514080LL*r*Power(xj, 2LL) + 64260LL*Power(r, 3LL)*Power(xj, 4LL) -

            1512LL*Power(r, 5LL)*Power(xj, 6LL) + 8LL*Power(r, 7LL)*Power(xj, 8LL)) +

        6LL*r*Power(xi, 11LL)*Power(xj, 16LL)*

        (100781100LL*r*Power(xj, 2LL) - 4663260LL*Power(r, 3LL)*Power(xj, 4LL) +

            128376LL*Power(r, 5LL)*Power(xj, 6LL) + 40LL*Power(r, 7LL)*Power(xj, 8LL)) -

        18LL*r*Power(xi, 13LL)*Power(xj, 14LL)*

        (-311339700LL*r*Power(xj, 2LL) - 29641080LL*Power(r, 3LL)*Power(xj, 4LL) -

            9192LL*Power(r, 5LL)*Power(xj, 6LL) + 208LL*Power(r, 7LL)*Power(xj, 8LL)) +

        18LL*r*Power(xi, 15LL)*Power(xj, 12LL)*

        (2658257700LL*r*Power(xj, 2LL) - 30435720LL*Power(r, 3LL)*Power(xj, 4LL) -

            697428LL*Power(r, 5LL)*Power(xj, 6LL) + 592LL*Power(r, 7LL)*Power(xj, 8LL)) -

        18LL*Power(xi, 12LL)*Power(xj, 14LL)*

        (358142400LL*r*Power(xj, 2LL) + 6211800LL*Power(r, 3LL)*Power(xj, 4LL) +

            1774920LL*Power(r, 5LL)*Power(xj, 6LL) + 1168LL*Power(r, 7LL)*Power(xj, 8LL)) +

        18LL*r*Power(xi, 25LL)*Power(xj, 2LL)*

        (-2275149870LL*r*Power(xj, 2LL) + 151339020LL*Power(r, 3LL)*Power(xj, 4LL) -

            1638372LL*Power(r, 5LL)*Power(xj, 6LL) + 1368LL*Power(r, 7LL)*Power(xj, 8LL)) -

        9LL*r*Power(xi, 19LL)*Power(xj, 8LL)*

        (15227779140LL*r*Power(xj, 2LL) - 682098480LL*Power(r, 3LL)*Power(xj, 4LL) +

            2387616LL*Power(r, 5LL)*Power(xj, 6LL) + 2736LL*Power(r, 7LL)*Power(xj, 8LL)) -

        3LL*r*Power(xi, 23LL)*Power(xj, 4LL)*

        (-22236093180LL*r*Power(xj, 2LL) + 1310447880LL*Power(r, 3LL)*Power(xj, 4LL) -

            17525448LL*Power(r, 5LL)*Power(xj, 6LL) + 20672LL*Power(r, 7LL)*Power(xj, 8LL))

        + 3LL*r*Power(xi, 21LL)*Power(xj, 6LL)*

        (38061529380LL*r*Power(xj, 2LL) - 567904680LL*Power(r, 3LL)*Power(xj, 4LL) -

            8651232LL*Power(r, 5LL)*Power(xj, 6LL) + 20672LL*Power(r, 7LL)*Power(xj, 8LL)) +

        63LL*Power(xi, 20LL)*Power(xj, 6LL)*

        (12480405840LL*r*Power(xj, 2LL) - 805257240LL*Power(r, 3LL)*Power(xj, 4LL) +

            5736480LL*Power(r, 5LL)*Power(xj, 6LL) + 20672LL*Power(r, 7LL)*Power(xj, 8LL)) +

        18LL*Power(xi, 14LL)*Power(xj, 12LL)*

        (-5038412400LL*r*Power(xj, 2LL) - 478879800LL*Power(r, 3LL)*Power(xj, 4LL) +

            1093680LL*Power(r, 5LL)*Power(xj, 6LL) + 21872LL*Power(r, 7LL)*Power(xj, 8LL)) -

        18LL*Power(xi, 26LL)*(3589563600LL*r*Power(xj, 2LL) +

                              269348940LL*Power(r, 3LL)*Power(xj, 4LL) -

                              9958200LL*Power(r, 5LL)*Power(xj, 6LL) + 32712LL*Power(r, 7LL)*Power(xj, 8LL)) +

        9LL*Power(xi, 18LL)*Power(xj, 8LL)*

        (16656781680LL*r*Power(xj, 2LL) + 3648168720LL*Power(r, 3LL)*Power(xj, 4LL) -

            77054880LL*Power(r, 5LL)*Power(xj, 6LL) + 85424LL*Power(r, 7LL)*Power(xj, 8LL))

        - 9LL*Power(xi, 16LL)*Power(xj, 10LL)*(42596154720LL*r*Power(xj, 2LL) -

                                               1071714000LL*Power(r, 3LL)*Power(xj, 4LL) -

                                               32749920LL*Power(r, 5LL)*Power(xj, 6LL) + 117776LL*Power(r, 7LL)*Power(xj, 8LL))

        + 18LL*Power(xi, 24LL)*Power(xj, 2LL)*(-17732214360LL*r*Power(xj, 2LL) +

                                               1597089060LL*Power(r, 3LL)*Power(xj, 4LL) -

                                               35554680LL*Power(r, 5LL)*Power(xj, 6LL) + 141208LL*Power(r, 7LL)*Power(xj, 8LL))

        - 9LL*Power(xi, 22LL)*Power(xj, 4LL)*(7819726320LL*r*Power(xj, 2LL) +

                                              365683080LL*Power(r, 3LL)*Power(xj, 4LL) -

                                              52572240LL*Power(r, 5LL)*Power(xj, 6LL) + 351424LL*Power(r, 7LL)*Power(xj, 8LL))

        + Power(xi, 17LL)*Power(xj, 10LL)*(1073961177975LL -

                                           21753487980LL*Power(r, 2LL)*Power(xj, 2LL) -

                                           745994340LL*Power(r, 4LL)*Power(xj, 4LL) +

                                           5307156LL*Power(r, 6LL)*Power(xj, 6LL) - 818LL*Power(r, 8LL)*Power(xj, 8LL)) +

        10LL*Power(xi, 9LL)*Power(xj, 18LL)*

        (49448070LL - 6409935LL*Power(r, 2LL)*Power(xj, 2LL) +

            161595LL*Power(r, 4LL)*Power(xj, 4LL) - 1026LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) +

        6LL*Power(xi, 11LL)*Power(xj, 16LL)*

        (-207559800LL + 50390550LL*Power(r, 2LL)*Power(xj, 2LL) -

            1165815LL*Power(r, 4LL)*Power(xj, 4LL) + 21396LL*Power(r, 6LL)*Power(xj, 6LL) +

            5LL*Power(r, 8LL)*Power(xj, 8LL)) -

        18LL*Power(xi, 13LL)*Power(xj, 14LL)*

        (-1703720025LL - 155669850LL*Power(r, 2LL)*Power(xj, 2LL) -

            7410270LL*Power(r, 4LL)*Power(xj, 4LL) - 1532LL*Power(r, 6LL)*Power(xj, 6LL) +

            26LL*Power(r, 8LL)*Power(xj, 8LL)) +

        18LL*Power(xi, 15LL)*Power(xj, 12LL)*

        (19380896325LL + 1329128850LL*Power(r, 2LL)*Power(xj, 2LL) -

            7608930LL*Power(r, 4LL)*Power(xj, 4LL) -

            116238LL*Power(r, 6LL)*Power(xj, 6LL) + 74LL*Power(r, 8LL)*Power(xj, 8LL)) +

        18LL*Power(xi, 25LL)*Power(xj, 2LL)*

        (-5449970925LL - 1137574935LL*Power(r, 2LL)*Power(xj, 2LL) +

            37834755LL*Power(r, 4LL)*Power(xj, 4LL) -

            273062LL*Power(r, 6LL)*Power(xj, 6LL) + 171LL*Power(r, 8LL)*Power(xj, 8LL)) -

        9LL*Power(xi, 19LL)*Power(xj, 8LL)*

        (-37914907275LL + 7613889570LL*Power(r, 2LL)*Power(xj, 2LL) -

            170524620LL*Power(r, 4LL)*Power(xj, 4LL) +

            397936LL*Power(r, 6LL)*Power(xj, 6LL) + 342LL*Power(r, 8LL)*Power(xj, 8LL)) -

        3LL*Power(xi, 23LL)*Power(xj, 4LL)*

        (219130630425LL - 11118046590LL*Power(r, 2LL)*Power(xj, 2LL) +

            327611970LL*Power(r, 4LL)*Power(xj, 4LL) -

            2920908LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

        3LL*Power(xi, 21LL)*Power(xj, 6LL)*

        (-345162539925LL + 19030764690LL*Power(r, 2LL)*Power(xj, 2LL) -

            141976170LL*Power(r, 4LL)*Power(xj, 4LL) -

            1441872LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

        Power(xi, 27LL)*(-2884470750LL - 19229805000LL*Power(r, 2LL)*Power(xj, 2LL) +

                         141664950LL*Power(r, 4LL)*Power(xj, 4LL) +

                         406728LL*Power(r, 6LL)*Power(xj, 6LL) + 7362LL*Power(r, 8LL)*Power(xj, 8LL))) +

                 2LL*exp(2LL*r*xj)*Power(xj, 13LL)*

                 (-630LL*Power(r, 8LL)*Power(xi, 34LL) - 10LL*Power(r, 9LL)*Power(xi, 35LL) +

        70875LL*Power(xj, 26LL) + 127575LL*r*xi*Power(xj, 26LL) -

        30LL*Power(r, 7LL)*Power(xi, 33LL)*(630LL + Power(r, 2LL)*Power(xj, 2LL)) +

        14175LL*Power(xi, 2LL)*Power(xj, 24LL)*(-95LL + 8LL*Power(r, 2LL)*Power(xj, 2LL)) +

        4725LL*r*Power(xi, 3LL)*Power(xj, 24LL)*

        (-513LL + 14LL*Power(r, 2LL)*Power(xj, 2LL)) -

        90LL*Power(r, 6LL)*Power(xi, 32LL)*(3920LL + 43LL*Power(r, 2LL)*Power(xj, 2LL)) +

        4725LL*r*Power(xi, 5LL)*Power(xj, 22LL)*

        (4617LL - 266LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        14175LL*Power(xi, 4LL)*Power(xj, 22LL)*

        (855LL - 152LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        36LL*Power(r, 5LL)*Power(xi, 31LL)*

        (-124950LL - 4985LL*Power(r, 2LL)*Power(xj, 2LL) +

            13LL*Power(r, 4LL)*Power(xj, 4LL)) +

        36LL*Power(r, 4LL)*Power(xi, 30LL)*

        (-1124550LL - 127960LL*Power(r, 2LL)*Power(xj, 2LL) +

            863LL*Power(r, 4LL)*Power(xj, 4LL)) +

        135LL*r*Power(xi, 7LL)*Power(xj, 20LL)*

        (-915705LL + 83790LL*Power(r, 2LL)*Power(xj, 2LL) -

            1330LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

        315LL*Power(xi, 6LL)*Power(xj, 20LL)*

        (-218025LL + 61560LL*Power(r, 2LL)*Power(xj, 2LL) -

            1710LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) -

        36LL*Power(r, 3LL)*Power(xi, 29LL)*

        (7122150LL + 2102730LL*Power(r, 2LL)*Power(xj, 2LL) -

            23294LL*Power(r, 4LL)*Power(xj, 4LL) + 37LL*Power(r, 6LL)*Power(xj, 6LL)) -

        36LL*Power(r, 2LL)*Power(xi, 28LL)*

        (30523500LL + 23401350LL*Power(r, 2LL)*Power(xj, 2LL) -

            299250LL*Power(r, 4LL)*Power(xj, 4LL) + 1297LL*Power(r, 6LL)*Power(xj, 6LL)) +

        r*Power(xi, 17LL)*Power(xj, 10LL)*

        (1073961177975LL - 21753487980LL*Power(r, 2LL)*Power(xj, 2LL) -

            745994340LL*Power(r, 4LL)*Power(xj, 4LL) +

            5307156LL*Power(r, 6LL)*Power(xj, 6LL) - 818LL*Power(r, 8LL)*Power(xj, 8LL)) +

        10LL*r*Power(xi, 9LL)*Power(xj, 18LL)*

        (49448070LL - 6409935LL*Power(r, 2LL)*Power(xj, 2LL) +

            161595LL*Power(r, 4LL)*Power(xj, 4LL) - 1026LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) +

        90LL*Power(xi, 8LL)*Power(xj, 18LL)*

        (3052350LL - 1220940LL*Power(r, 2LL)*Power(xj, 2LL) +

            53865LL*Power(r, 4LL)*Power(xj, 4LL) - 532LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) -

        1710LL*Power(xi, 10LL)*Power(xj, 16LL)*

        (481950LL - 257040LL*Power(r, 2LL)*Power(xj, 2LL) +

            16065LL*Power(r, 4LL)*Power(xj, 4LL) - 252LL*Power(r, 6LL)*Power(xj, 6LL) +

            Power(r, 8LL)*Power(xj, 8LL)) +

        6LL*r*Power(xi, 11LL)*Power(xj, 16LL)*

        (-207559800LL + 50390550LL*Power(r, 2LL)*Power(xj, 2LL) -

            1165815LL*Power(r, 4LL)*Power(xj, 4LL) + 21396LL*Power(r, 6LL)*Power(xj, 6LL) +

            5LL*Power(r, 8LL)*Power(xj, 8LL)) -

        18LL*r*Power(xi, 13LL)*Power(xj, 14LL)*

        (-1703720025LL - 155669850LL*Power(r, 2LL)*Power(xj, 2LL) -

            7410270LL*Power(r, 4LL)*Power(xj, 4LL) - 1532LL*Power(r, 6LL)*Power(xj, 6LL) +

            26LL*Power(r, 8LL)*Power(xj, 8LL)) +

        18LL*r*Power(xi, 15LL)*Power(xj, 12LL)*

        (19380896325LL + 1329128850LL*Power(r, 2LL)*Power(xj, 2LL) -

            7608930LL*Power(r, 4LL)*Power(xj, 4LL) -

            116238LL*Power(r, 6LL)*Power(xj, 6LL) + 74LL*Power(r, 8LL)*Power(xj, 8LL)) -

        18LL*Power(xi, 12LL)*Power(xj, 14LL)*

        (89026875LL + 179071200LL*Power(r, 2LL)*Power(xj, 2LL) +

            1552950LL*Power(r, 4LL)*Power(xj, 4LL) +

            295820LL*Power(r, 6LL)*Power(xj, 6LL) + 146LL*Power(r, 8LL)*Power(xj, 8LL)) +

        18LL*r*Power(xi, 25LL)*Power(xj, 2LL)*

        (-5449970925LL - 1137574935LL*Power(r, 2LL)*Power(xj, 2LL) +

            37834755LL*Power(r, 4LL)*Power(xj, 4LL) -

            273062LL*Power(r, 6LL)*Power(xj, 6LL) + 171LL*Power(r, 8LL)*Power(xj, 8LL)) -

        9LL*r*Power(xi, 19LL)*Power(xj, 8LL)*

        (-37914907275LL + 7613889570LL*Power(r, 2LL)*Power(xj, 2LL) -

            170524620LL*Power(r, 4LL)*Power(xj, 4LL) +

            397936LL*Power(r, 6LL)*Power(xj, 6LL) + 342LL*Power(r, 8LL)*Power(xj, 8LL)) -

        3LL*r*Power(xi, 23LL)*Power(xj, 4LL)*

        (219130630425LL - 11118046590LL*Power(r, 2LL)*Power(xj, 2LL) +

            327611970LL*Power(r, 4LL)*Power(xj, 4LL) -

            2920908LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

        3LL*r*Power(xi, 21LL)*Power(xj, 6LL)*

        (-345162539925LL + 19030764690LL*Power(r, 2LL)*Power(xj, 2LL) -

            141976170LL*Power(r, 4LL)*Power(xj, 4LL) -

            1441872LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

        63LL*Power(xi, 20LL)*Power(xj, 6LL)*

        (-50980542525LL + 6240202920LL*Power(r, 2LL)*Power(xj, 2LL) -

            201314310LL*Power(r, 4LL)*Power(xj, 4LL) +

            956080LL*Power(r, 6LL)*Power(xj, 6LL) + 2584LL*Power(r, 8LL)*Power(xj, 8LL)) +

        18LL*Power(xi, 14LL)*Power(xj, 12LL)*

        (-7803332775LL - 2519206200LL*Power(r, 2LL)*Power(xj, 2LL) -

            119719950LL*Power(r, 4LL)*Power(xj, 4LL) +

            182280LL*Power(r, 6LL)*Power(xj, 6LL) + 2734LL*Power(r, 8LL)*Power(xj, 8LL)) -

        18LL*Power(xi, 26LL)*(195859125LL + 1794781800LL*Power(r, 2LL)*Power(xj, 2LL) +

                              67337235LL*Power(r, 4LL)*Power(xj, 4LL) -

                              1659700LL*Power(r, 6LL)*Power(xj, 6LL) + 4089LL*Power(r, 8LL)*Power(xj, 8LL)) +

        9LL*Power(xi, 18LL)*Power(xj, 8LL)*

        (-357591274425LL + 8328390840LL*Power(r, 2LL)*Power(xj, 2LL) +

            912042180LL*Power(r, 4LL)*Power(xj, 4LL) -

            12842480LL*Power(r, 6LL)*Power(xj, 6LL) + 10678LL*Power(r, 8LL)*Power(xj, 8LL))

        - 9LL*Power(xi, 16LL)*Power(xj, 10LL)*(128599724925LL +

                                               21298077360LL*Power(r, 2LL)*Power(xj, 2LL) -

                                               267928500LL*Power(r, 4LL)*Power(xj, 4LL) -

                                               5458320LL*Power(r, 6LL)*Power(xj, 6LL) + 14722LL*Power(r, 8LL)*Power(xj, 8LL)) +

        18LL*Power(xi, 24LL)*Power(xj, 2LL)*

        (-7604930025LL - 8866107180LL*Power(r, 2LL)*Power(xj, 2LL) +

            399272265LL*Power(r, 4LL)*Power(xj, 4LL) -

            5925780LL*Power(r, 6LL)*Power(xj, 6LL) + 17651LL*Power(r, 8LL)*Power(xj, 8LL)) -

        9LL*Power(xi, 22LL)*Power(xj, 4LL)*

        (129194933175LL + 3909863160LL*Power(r, 2LL)*Power(xj, 2LL) +

            91420770LL*Power(r, 4LL)*Power(xj, 4LL) -

            8762040LL*Power(r, 6LL)*Power(xj, 6LL) + 43928LL*Power(r, 8LL)*Power(xj, 8LL)) +

        Power(xi, 27LL)*(-2884470750LL*r - 6409935000LL*Power(r, 3LL)*Power(xj, 2LL) +

                         28332990LL*Power(r, 5LL)*Power(xj, 4LL) + 58104LL*Power(r, 7LL)*Power(xj, 6LL) +

                         818LL*Power(r, 9LL)*Power(xj, 8LL))) +

                 exp(2LL*r*xi)*Power(xi, 12LL)*

                 (Power(xi, 8LL)*Power(xj, 18LL)*

        (-341234165475LL*xj - 786265567920LL*r*Power(xj, 2LL) -

            171276882210LL*Power(r, 2LL)*Power(xj, 3LL) +

            3291147720LL*Power(r, 3LL)*Power(xj, 4LL) +

            4914179550LL*Power(r, 4LL)*Power(xj, 5LL) +

            639984240LL*Power(r, 5LL)*Power(xj, 6LL) +

            34405812LL*Power(r, 6LL)*Power(xj, 7LL) +

            588816LL*Power(r, 7LL)*Power(xj, 8LL) - 7362LL*Power(r, 8LL)*Power(xj, 9LL)) +

        10LL*Power(xj, 26LL)*(288447075LL*xj + 219769200LL*r*Power(xj, 2LL) +

                              76919220LL*Power(r, 2LL)*Power(xj, 3LL) +

                              16193520LL*Power(r, 3LL)*Power(xj, 4LL) +

                              2249100LL*Power(r, 4LL)*Power(xj, 5LL) +

                              211680LL*Power(r, 5LL)*Power(xj, 6LL) + 13230LL*Power(r, 6LL)*Power(xj, 7LL) +

                              504LL*Power(r, 7LL)*Power(xj, 8LL) + 9LL*Power(r, 8LL)*Power(xj, 9LL)) +

        30LL*Power(xi, 2LL)*Power(xj, 24LL)*

        (3269982555LL*xj + 2153738160LL*r*Power(xj, 2LL) +

            640993500LL*Power(r, 2LL)*Power(xj, 3LL) +

            112326480LL*Power(r, 3LL)*Power(xj, 4LL) +

            12616380LL*Power(r, 4LL)*Power(xj, 5LL) +

            921312LL*Power(r, 5LL)*Power(xj, 6LL) + 41874LL*Power(r, 6LL)*Power(xj, 7LL) +

            1032LL*Power(r, 7LL)*Power(xj, 8LL) + 9LL*Power(r, 8LL)*Power(xj, 9LL)) -

        15LL*Power(xi, 24LL)*Power(xj, 2LL)*

        (-161595LL*xj - 287280LL*r*Power(xj, 2LL) -

            251370LL*Power(r, 2LL)*Power(xj, 3LL) - 143640LL*Power(r, 3LL)*Power(xj, 4LL) -

            59850LL*Power(r, 4LL)*Power(xj, 5LL) - 19152LL*Power(r, 5LL)*Power(xj, 6LL) -

            4788LL*Power(r, 6LL)*Power(xj, 7LL) - 912LL*Power(r, 7LL)*Power(xj, 8LL) +

            18LL*Power(r, 8LL)*Power(xj, 9LL)) -

        5LL*Power(xi, 26LL)*(25515LL*xj + 45360LL*r*Power(xj, 2LL) +

                             39690LL*Power(r, 2LL)*Power(xj, 3LL) + 22680LL*Power(r, 3LL)*Power(xj, 4LL) +

                             9450LL*Power(r, 4LL)*Power(xj, 5LL) + 3024LL*Power(r, 5LL)*Power(xj, 6LL) +

                             756LL*Power(r, 6LL)*Power(xj, 7LL) + 144LL*Power(r, 7LL)*Power(xj, 8LL) +

                             18LL*Power(r, 8LL)*Power(xj, 9LL)) -

        1938LL*Power(xi, 14LL)*Power(xj, 12LL)*

        (15824025LL*xj - 46796400LL*r*Power(xj, 2LL) +

            37034550LL*Power(r, 2LL)*Power(xj, 3LL) +

            4977000LL*Power(r, 3LL)*Power(xj, 4LL) -

            1924650LL*Power(r, 4LL)*Power(xj, 5LL) -

            357840LL*Power(r, 5LL)*Power(xj, 6LL) - 12936LL*Power(r, 6LL)*Power(xj, 7LL) +

            672LL*Power(r, 7LL)*Power(xj, 8LL) + 36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        1938LL*Power(xi, 12LL)*Power(xj, 14LL)*

        (-180008325LL*xj + 197814960LL*r*Power(xj, 2LL) +

            33674130LL*Power(r, 2LL)*Power(xj, 3LL) -

            16941960LL*Power(r, 3LL)*Power(xj, 4LL) -

            3959550LL*Power(r, 4LL)*Power(xj, 5LL) -

            186480LL*Power(r, 5LL)*Power(xj, 6LL) + 15624LL*Power(r, 6LL)*Power(xj, 7LL) +

            1632LL*Power(r, 7LL)*Power(xj, 8LL) + 36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        342LL*Power(xi, 16LL)*Power(xj, 10LL)*

        (3641400LL*xj + 18849600LL*r*Power(xj, 2LL) -

            24579450LL*Power(r, 2LL)*Power(xj, 3LL) +

            25204200LL*Power(r, 3LL)*Power(xj, 4LL) +

            2002350LL*Power(r, 4LL)*Power(xj, 5LL) -

            861840LL*Power(r, 5LL)*Power(xj, 6LL) - 108626LL*Power(r, 6LL)*Power(xj, 7LL) -

            2248LL*Power(r, 7LL)*Power(xj, 8LL) + 81LL*Power(r, 8LL)*Power(xj, 9LL)) -

        171LL*Power(xi, 10LL)*Power(xj, 16LL)*

        (6280474725LL*xj + 876672720LL*r*Power(xj, 2LL) -

            1202193090LL*Power(r, 2LL)*Power(xj, 3LL) -

            296673720LL*Power(r, 3LL)*Power(xj, 4LL) -

            12454050LL*Power(r, 4LL)*Power(xj, 5LL) +

            2766960LL*Power(r, 5LL)*Power(xj, 6LL) +

            358708LL*Power(r, 6LL)*Power(xj, 7LL) + 14864LL*Power(r, 7LL)*Power(xj, 8LL) +

            162LL*Power(r, 8LL)*Power(xj, 9LL)) +

        9LL*Power(xi, 22LL)*Power(xj, 4LL)*

        (-2423925LL*xj - 4309200LL*r*Power(xj, 2LL) -

            3770550LL*Power(r, 2LL)*Power(xj, 3LL) -

            2154600LL*Power(r, 3LL)*Power(xj, 4LL) -

            897750LL*Power(r, 4LL)*Power(xj, 5LL) - 287280LL*Power(r, 5LL)*Power(xj, 6LL) -

            99848LL*Power(r, 6LL)*Power(xj, 7LL) + 2336LL*Power(r, 7LL)*Power(xj, 8LL) +

            468LL*Power(r, 8LL)*Power(xj, 9LL)) -

        9LL*Power(xi, 4LL)*Power(xj, 22LL)*

        (-73043543475LL*xj - 35464428720LL*r*Power(xj, 2LL) -

            6825449610LL*Power(r, 2LL)*Power(xj, 3LL) -

            538697880LL*Power(r, 3LL)*Power(xj, 4LL) +

            15740550LL*Power(r, 4LL)*Power(xj, 5LL) +

            7182000LL*Power(r, 5LL)*Power(xj, 6LL) +

            652232LL*Power(r, 6LL)*Power(xj, 7LL) + 27616LL*Power(r, 7LL)*Power(xj, 8LL) +

            468LL*Power(r, 8LL)*Power(xj, 9LL)) +

        9LL*Power(xi, 6LL)*Power(xj, 20LL)*

        (115054179975LL*xj + 7819726320LL*r*Power(xj, 2LL) -

            11118046590LL*Power(r, 2LL)*Power(xj, 3LL) -

            3194178120LL*Power(r, 3LL)*Power(xj, 4LL) -

            378347550LL*Power(r, 4LL)*Power(xj, 5LL) -

            19916400LL*Power(r, 5LL)*Power(xj, 6LL) -

            45192LL*Power(r, 6LL)*Power(xj, 7LL) + 41504LL*Power(r, 7LL)*Power(xj, 8LL) +

            1332LL*Power(r, 8LL)*Power(xj, 9LL)) -

        9LL*Power(xi, 20LL)*Power(xj, 6LL)*

        (-13735575LL*xj - 24418800LL*r*Power(xj, 2LL) -

            21366450LL*Power(r, 2LL)*Power(xj, 3LL) -

            12209400LL*Power(r, 3LL)*Power(xj, 4LL) -

            3886050LL*Power(r, 4LL)*Power(xj, 5LL) -

            3549840LL*Power(r, 5LL)*Power(xj, 6LL) + 21448LL*Power(r, 6LL)*Power(xj, 7LL) +

            43744LL*Power(r, 7LL)*Power(xj, 8LL) + 1332LL*Power(r, 8LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 18LL)*Power(xj, 8LL)*

        (-247240350LL*xj - 439538400LL*r*Power(xj, 2LL) -

            453514950LL*Power(r, 2LL)*Power(xj, 3LL) +

            55906200LL*Power(r, 3LL)*Power(xj, 4LL) -

            333462150LL*Power(r, 4LL)*Power(xj, 5LL) -

            9843120LL*Power(r, 5LL)*Power(xj, 6LL) +

            7322994LL*Power(r, 6LL)*Power(xj, 7LL) + 529992LL*Power(r, 7LL)*Power(xj, 8LL) +

            3681LL*Power(r, 8LL)*Power(xj, 9LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 13LL)*

                 (Power(xi, 8LL)*Power(xj, 18LL)*

        (3218321469825LL - 341234165475LL*r*xj -

            393132783960LL*Power(r, 2LL)*Power(xj, 2LL) -

            57092294070LL*Power(r, 3LL)*Power(xj, 3LL) +

            822786930LL*Power(r, 4LL)*Power(xj, 4LL) +

            982835910LL*Power(r, 5LL)*Power(xj, 5LL) +

            106664040LL*Power(r, 6LL)*Power(xj, 6LL) +

            4915116LL*Power(r, 7LL)*Power(xj, 7LL) + 73602LL*Power(r, 8LL)*Power(xj, 8LL) -

            818LL*Power(r, 9LL)*Power(xj, 9LL)) +

        10LL*Power(xj, 26LL)*(352546425LL + 288447075LL*r*xj +

                              109884600LL*Power(r, 2LL)*Power(xj, 2LL) +

                              25639740LL*Power(r, 3LL)*Power(xj, 3LL) +

                              4048380LL*Power(r, 4LL)*Power(xj, 4LL) + 449820LL*Power(r, 5LL)*Power(xj, 5LL) +

                              35280LL*Power(r, 6LL)*Power(xj, 6LL) + 1890LL*Power(r, 7LL)*Power(xj, 7LL) +

                              63LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) +

        30LL*Power(xi, 2LL)*Power(xj, 24LL)*

        (4562958015LL + 3269982555LL*r*xj +

            1076869080LL*Power(r, 2LL)*Power(xj, 2LL) +

            213664500LL*Power(r, 3LL)*Power(xj, 3LL) +

            28081620LL*Power(r, 4LL)*Power(xj, 4LL) +

            2523276LL*Power(r, 5LL)*Power(xj, 5LL) + 153552LL*Power(r, 6LL)*Power(xj, 6LL) +

            5982LL*Power(r, 7LL)*Power(xj, 7LL) + 129LL*Power(r, 8LL)*Power(xj, 8LL) +

            Power(r, 9LL)*Power(xj, 9LL)) -

        15LL*Power(xi, 24LL)*Power(xj, 2LL)*

        (-89775LL - 161595LL*r*xj - 143640LL*Power(r, 2LL)*Power(xj, 2LL) -

            83790LL*Power(r, 3LL)*Power(xj, 3LL) - 35910LL*Power(r, 4LL)*Power(xj, 4LL) -

            11970LL*Power(r, 5LL)*Power(xj, 5LL) - 3192LL*Power(r, 6LL)*Power(xj, 6LL) -

            684LL*Power(r, 7LL)*Power(xj, 7LL) - 114LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        5LL*Power(xi, 26LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                             13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                             1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                             108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                             2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        1938LL*Power(xi, 14LL)*Power(xj, 12LL)*

        (-826875LL + 15824025LL*r*xj - 23398200LL*Power(r, 2LL)*Power(xj, 2LL) +

            12344850LL*Power(r, 3LL)*Power(xj, 3LL) +

            1244250LL*Power(r, 4LL)*Power(xj, 4LL) - 384930LL*Power(r, 5LL)*Power(xj, 5LL) -

            59640LL*Power(r, 6LL)*Power(xj, 6LL) - 1848LL*Power(r, 7LL)*Power(xj, 7LL) +

            84LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        1938LL*Power(xi, 12LL)*Power(xj, 14LL)*

        (72476775LL - 180008325LL*r*xj + 98907480LL*Power(r, 2LL)*Power(xj, 2LL) +

            11224710LL*Power(r, 3LL)*Power(xj, 3LL) -

            4235490LL*Power(r, 4LL)*Power(xj, 4LL) - 791910LL*Power(r, 5LL)*Power(xj, 5LL) -

            31080LL*Power(r, 6LL)*Power(xj, 6LL) + 2232LL*Power(r, 7LL)*Power(xj, 7LL) +

            204LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        342LL*Power(xi, 16LL)*Power(xj, 10LL)*

        (2409750LL + 3641400LL*r*xj + 9424800LL*Power(r, 2LL)*Power(xj, 2LL) -

            8193150LL*Power(r, 3LL)*Power(xj, 3LL) +

            6301050LL*Power(r, 4LL)*Power(xj, 4LL) + 400470LL*Power(r, 5LL)*Power(xj, 5LL) -

            143640LL*Power(r, 6LL)*Power(xj, 6LL) - 15518LL*Power(r, 7LL)*Power(xj, 7LL) -

            281LL*Power(r, 8LL)*Power(xj, 8LL) + 9LL*Power(r, 9LL)*Power(xj, 9LL)) -

        171LL*Power(xi, 10LL)*Power(xj, 16LL)*

        (-6768406575LL + 6280474725LL*r*xj +

            438336360LL*Power(r, 2LL)*Power(xj, 2LL) -

            400731030LL*Power(r, 3LL)*Power(xj, 3LL) -

            74168430LL*Power(r, 4LL)*Power(xj, 4LL) -

            2490810LL*Power(r, 5LL)*Power(xj, 5LL) + 461160LL*Power(r, 6LL)*Power(xj, 6LL) +

            51244LL*Power(r, 7LL)*Power(xj, 7LL) + 1858LL*Power(r, 8LL)*Power(xj, 8LL) +

            18LL*Power(r, 9LL)*Power(xj, 9LL)) +

        9LL*Power(xi, 22LL)*Power(xj, 4LL)*

        (-1346625LL - 2423925LL*r*xj - 2154600LL*Power(r, 2LL)*Power(xj, 2LL) -

            1256850LL*Power(r, 3LL)*Power(xj, 3LL) - 538650LL*Power(r, 4LL)*Power(xj, 4LL) -

            179550LL*Power(r, 5LL)*Power(xj, 5LL) - 47880LL*Power(r, 6LL)*Power(xj, 6LL) -

            14264LL*Power(r, 7LL)*Power(xj, 7LL) + 292LL*Power(r, 8LL)*Power(xj, 8LL) +

            52LL*Power(r, 9LL)*Power(xj, 9LL)) -

        9LL*Power(xi, 4LL)*Power(xj, 22LL)*

        (-129194933175LL - 73043543475LL*r*xj -

            17732214360LL*Power(r, 2LL)*Power(xj, 2LL) -

            2275149870LL*Power(r, 3LL)*Power(xj, 3LL) -

            134674470LL*Power(r, 4LL)*Power(xj, 4LL) +

            3148110LL*Power(r, 5LL)*Power(xj, 5LL) +

            1197000LL*Power(r, 6LL)*Power(xj, 6LL) + 93176LL*Power(r, 7LL)*Power(xj, 7LL) +

            3452LL*Power(r, 8LL)*Power(xj, 8LL) + 52LL*Power(r, 9LL)*Power(xj, 9LL)) +

        9LL*Power(xi, 6LL)*Power(xj, 20LL)*

        (356863797675LL + 115054179975LL*r*xj +

            3909863160LL*Power(r, 2LL)*Power(xj, 2LL) -

            3706015530LL*Power(r, 3LL)*Power(xj, 3LL) -

            798544530LL*Power(r, 4LL)*Power(xj, 4LL) -

            75669510LL*Power(r, 5LL)*Power(xj, 5LL) -

            3319400LL*Power(r, 6LL)*Power(xj, 6LL) - 6456LL*Power(r, 7LL)*Power(xj, 7LL) +

            5188LL*Power(r, 8LL)*Power(xj, 8LL) + 148LL*Power(r, 9LL)*Power(xj, 9LL)) -

        9LL*Power(xi, 20LL)*Power(xj, 6LL)*

        (-7630875LL - 13735575LL*r*xj - 12209400LL*Power(r, 2LL)*Power(xj, 2LL) -

            7122150LL*Power(r, 3LL)*Power(xj, 3LL) -

            3052350LL*Power(r, 4LL)*Power(xj, 4LL) - 777210LL*Power(r, 5LL)*Power(xj, 5LL) -

            591640LL*Power(r, 6LL)*Power(xj, 6LL) + 3064LL*Power(r, 7LL)*Power(xj, 7LL) +

            5468LL*Power(r, 8LL)*Power(xj, 8LL) + 148LL*Power(r, 9LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 18LL)*Power(xj, 8LL)*

        (-137355750LL - 247240350LL*r*xj - 219769200LL*Power(r, 2LL)*Power(xj, 2LL) -

            151171650LL*Power(r, 3LL)*Power(xj, 3LL) +

            13976550LL*Power(r, 4LL)*Power(xj, 4LL) -

            66692430LL*Power(r, 5LL)*Power(xj, 5LL) -

            1640520LL*Power(r, 6LL)*Power(xj, 6LL) + 1046142LL*Power(r, 7LL)*Power(xj, 7LL) +

            66249LL*Power(r, 8LL)*Power(xj, 8LL) + 409LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (70875LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 19LL)*Power(xi + xj, 19LL))

            ;
        }

    }
    return S;
}
