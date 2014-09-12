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

cl_R DSlater_3S_5S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-568188982486125LL*xi + 627683696640000LL*exp(2LL*r*xi)*xi -

                  1017388536664500LL*r*Power(xi, 2LL) -

                  899677411132500LL*Power(r, 2LL)*Power(xi, 3LL) -

                  523015260768000LL*Power(r, 3LL)*Power(xi, 4LL) -

                  224405775594000LL*Power(r, 4LL)*Power(xi, 5LL) -

                  75610821686400LL*Power(r, 5LL)*Power(xi, 6LL) -

                  20775676521600LL*Power(r, 6LL)*Power(xi, 7LL) -

                  4769897932800LL*Power(r, 7LL)*Power(xi, 8LL) -

                  929382854400LL*Power(r, 8LL)*Power(xi, 9LL) -

                  154983628800LL*Power(r, 9LL)*Power(xi, 10LL) -

                  22140518400LL*Power(r, 10LL)*Power(xi, 11LL) -

                  2683699200LL*Power(r, 11LL)*Power(xi, 12LL) -

                  268369920LL*Power(r, 12LL)*Power(xi, 13LL) -

                  20643840LL*Power(r, 13LL)*Power(xi, 14LL) - 983040LL*Power(r, 14LL)*Power(xi, 15LL))/

                (3.1384184832e14*exp(2LL*r*xi)*r) +

                (-313841848320000LL + 313841848320000LL*exp(2LL*r*xi) -

                 568188982486125LL*r*xi - 508694268332250LL*Power(r, 2LL)*Power(xi, 2LL) -

                 299892470377500LL*Power(r, 3LL)*Power(xi, 3LL) -

                 130753815192000LL*Power(r, 4LL)*Power(xi, 4LL) -

                 44881155118800LL*Power(r, 5LL)*Power(xi, 5LL) -

                 12601803614400LL*Power(r, 6LL)*Power(xi, 6LL) -

                 2967953788800LL*Power(r, 7LL)*Power(xi, 7LL) -

                 596237241600LL*Power(r, 8LL)*Power(xi, 8LL) -

                 103264761600LL*Power(r, 9LL)*Power(xi, 9LL) -

                 15498362880LL*Power(r, 10LL)*Power(xi, 10LL) -

                 2012774400LL*Power(r, 11LL)*Power(xi, 11LL) -

                 223641600LL*Power(r, 12LL)*Power(xi, 12LL) -

                 20643840LL*Power(r, 13LL)*Power(xi, 13LL) - 1474560LL*Power(r, 14LL)*Power(xi, 14LL) -

                 65536LL*Power(r, 15LL)*Power(xi, 15LL))/

                (3.1384184832e14*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-313841848320000LL + 313841848320000LL*exp(2LL*r*xi) -

                     568188982486125LL*r*xi - 508694268332250LL*Power(r, 2LL)*Power(xi, 2LL) -

                     299892470377500LL*Power(r, 3LL)*Power(xi, 3LL) -

                     130753815192000LL*Power(r, 4LL)*Power(xi, 4LL) -

                     44881155118800LL*Power(r, 5LL)*Power(xi, 5LL) -

                     12601803614400LL*Power(r, 6LL)*Power(xi, 6LL) -

                     2967953788800LL*Power(r, 7LL)*Power(xi, 7LL) -

                     596237241600LL*Power(r, 8LL)*Power(xi, 8LL) -

                     103264761600LL*Power(r, 9LL)*Power(xi, 9LL) -

                     15498362880LL*Power(r, 10LL)*Power(xi, 10LL) -

                     2012774400LL*Power(r, 11LL)*Power(xi, 11LL) -

                     223641600LL*Power(r, 12LL)*Power(xi, 12LL) -

                     20643840LL*Power(r, 13LL)*Power(xi, 13LL) -

                     1474560LL*Power(r, 14LL)*Power(xi, 14LL) - 65536LL*Power(r, 15LL)*Power(xi, 15LL)))/

                (1.5692092416e14*exp(2LL*r*xi)*r)

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
            S = (42525LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                 189LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-350LL*Power(r, 4LL)*Power(xi, 22LL) - 10LL*Power(r, 5LL)*Power(xi, 23LL) +

                  225LL*Power(xj, 18LL) + 375LL*r*xi*Power(xj, 18LL) -

                  70LL*Power(r, 3LL)*Power(xi, 21LL)*(75LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  75LL*r*Power(xi, 3LL)*Power(xj, 16LL)*(-75LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  75LL*Power(xi, 2LL)*Power(xj, 16LL)*(-45LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  50LL*Power(r, 2LL)*Power(xi, 20LL)*(840LL + 71LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 9LL)*Power(xj, 10LL)*

                  (4694625LL + 124800LL*Power(r, 2LL)*Power(xj, 2LL) -

            248LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  20LL*r*Power(xi, 17LL)*Power(xj, 2LL)*

                  (-185895LL - 948LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 5LL*r*Power(xi, 5LL)*Power(xj, 14LL)*

                  (7875LL - 450LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  25LL*Power(xi, 4LL)*Power(xj, 14LL)*

                  (945LL - 180LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  375LL*Power(xi, 6LL)*Power(xj, 12LL)*

                  (273LL - 84LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  5LL*r*Power(xi, 11LL)*Power(xj, 8LL)*

                  (-2803125LL + 49140LL*Power(r, 2LL)*Power(xj, 2LL) +

            8LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  5LL*r*Power(xi, 7LL)*Power(xj, 12LL)*

                  (-16965LL + 5152LL*Power(r, 2LL)*Power(xj, 2LL) +

            14LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  325LL*Power(xi, 10LL)*Power(xj, 8LL)*

                  (-60117LL - 5340LL*Power(r, 2LL)*Power(xj, 2LL) +

            40LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*r*Power(xi, 15LL)*Power(xj, 4LL)*

                  (845085LL - 22960LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*r*Power(xi, 13LL)*Power(xj, 6LL)*

                  (-139125LL - 10140LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  75LL*Power(xi, 12LL)*Power(xj, 6LL)*

                  (-729687LL + 25532LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  60LL*Power(xi, 18LL)*(-5355LL - 11940LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        86LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  2LL*r*Power(xi, 19LL)*(-89250LL - 35425LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         124LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  100LL*Power(xi, 16LL)*Power(xj, 2LL)*

                  (-79713LL - 13311LL*Power(r, 2LL)*Power(xj, 2LL) +

            146LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  5LL*Power(xi, 8LL)*Power(xj, 10LL)*

                  (157365LL + 95940LL*Power(r, 2LL)*Power(xj, 2LL) +

            952LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*Power(xi, 14LL)*Power(xj, 4LL)*

                  (2638467LL - 157500LL*Power(r, 2LL)*Power(xj, 2LL) +

            1820LL*Power(r, 4LL)*Power(xj, 4LL))) +

                 exp(2LL*r*xi)*Power(xi, 8LL)*

                 (2LL*Power(xi, 2LL)*Power(xj, 20LL)*

                  (1782492075LL + 1449175455LL*r*xj +

            533365560LL*Power(r, 2LL)*Power(xj, 2LL) +

            114631335LL*Power(r, 3LL)*Power(xj, 3LL) +

            15221115LL*Power(r, 4LL)*Power(xj, 4LL) +

            1142505LL*Power(r, 5LL)*Power(xj, 5LL) + 18396LL*Power(r, 6LL)*Power(xj, 6LL) -

            5238LL*Power(r, 7LL)*Power(xj, 7LL) - 513LL*Power(r, 8LL)*Power(xj, 8LL) -

            17LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  42LL*Power(xi, 4LL)*Power(xj, 18LL)*

                  (251336925LL + 104824125LL*r*xj + 340200LL*Power(r, 2LL)*Power(xj, 2LL) -

            9122085LL*Power(r, 3LL)*Power(xj, 3LL) -

            2798145LL*Power(r, 4LL)*Power(xj, 4LL) -

            433755LL*Power(r, 5LL)*Power(xj, 5LL) - 39060LL*Power(r, 6LL)*Power(xj, 6LL) -

            1890LL*Power(r, 7LL)*Power(xj, 7LL) - 27LL*Power(r, 8LL)*Power(xj, 8LL) +

            Power(r, 9LL)*Power(xj, 9LL)) +

                  6LL*Power(xj, 22LL)*(34459425LL + 34459425LL*r*xj +

                                       16216200LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       4729725LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       945945LL*Power(r, 4LL)*Power(xj, 4LL) + 135135LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       13860LL*Power(r, 6LL)*Power(xj, 6LL) + 990LL*Power(r, 7LL)*Power(xj, 7LL) +

                                       45LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) -

                  3LL*Power(xi, 22LL)*(14175LL + 25515LL*r*xj +

                                       22680LL*Power(r, 2LL)*Power(xj, 2LL) + 13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       5670LL*Power(r, 4LL)*Power(xj, 4LL) + 1890LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       504LL*Power(r, 6LL)*Power(xj, 6LL) + 108LL*Power(r, 7LL)*Power(xj, 7LL) +

                                       18LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  21LL*Power(xi, 18LL)*Power(xj, 4LL)*

                  (212625LL + 382725LL*r*xj + 340200LL*Power(r, 2LL)*Power(xj, 2LL) +

            198450LL*Power(r, 3LL)*Power(xj, 3LL) + 85050LL*Power(r, 4LL)*Power(xj, 4LL) +

            28350LL*Power(r, 5LL)*Power(xj, 5LL) + 7560LL*Power(r, 6LL)*Power(xj, 6LL) +

            1836LL*Power(r, 7LL)*Power(xj, 7LL) + 162LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  54LL*Power(xi, 6LL)*Power(xj, 16LL)*

                  (133451955LL - 73700865LL*r*xj - 54096840LL*Power(r, 2LL)*Power(xj, 2LL) -

            8306235LL*Power(r, 3LL)*Power(xj, 3LL) +

            966945LL*Power(r, 4LL)*Power(xj, 4LL) + 516747LL*Power(r, 5LL)*Power(xj, 5LL) +

            80724LL*Power(r, 6LL)*Power(xj, 6LL) + 6434LL*Power(r, 7LL)*Power(xj, 7LL) +

            251LL*Power(r, 8LL)*Power(xj, 8LL) + 3LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  315LL*Power(xi, 12LL)*Power(xj, 10LL)*

                  (-405405LL - 710073LL*r*xj - 805896LL*Power(r, 2LL)*Power(xj, 2LL) -

            101556LL*Power(r, 3LL)*Power(xj, 3LL) - 258804LL*Power(r, 4LL)*Power(xj, 4LL) -

            90972LL*Power(r, 5LL)*Power(xj, 5LL) - 9744LL*Power(r, 6LL)*Power(xj, 6LL) +

            120LL*Power(r, 7LL)*Power(xj, 7LL) + 84LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  315LL*Power(xi, 10LL)*Power(xj, 12LL)*

                  (-482895LL - 2656395LL*r*xj + 1186920LL*Power(r, 2LL)*Power(xj, 2LL) -

            1155420LL*Power(r, 3LL)*Power(xj, 3LL) -

            643356LL*Power(r, 4LL)*Power(xj, 4LL) - 93492LL*Power(r, 5LL)*Power(xj, 5LL) +

            336LL*Power(r, 6LL)*Power(xj, 6LL) + 1368LL*Power(r, 7LL)*Power(xj, 7LL) +

            132LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  27LL*Power(xi, 16LL)*Power(xj, 6LL)*

                  (-716625LL - 1289925LL*r*xj - 1146600LL*Power(r, 2LL)*Power(xj, 2LL) -

            668850LL*Power(r, 3LL)*Power(xj, 3LL) - 286650LL*Power(r, 4LL)*Power(xj, 4LL) -

            90006LL*Power(r, 5LL)*Power(xj, 5LL) - 32872LL*Power(r, 6LL)*Power(xj, 6LL) -

            4812LL*Power(r, 7LL)*Power(xj, 7LL) - 178LL*Power(r, 8LL)*Power(xj, 8LL) +

            6LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  Power(xi, 20LL)*Power(xj, 2LL)*

                  (637875LL + 1148175LL*r*xj + 1020600LL*Power(r, 2LL)*Power(xj, 2LL) +

            595350LL*Power(r, 3LL)*Power(xj, 3LL) + 255150LL*Power(r, 4LL)*Power(xj, 4LL) +

            85050LL*Power(r, 5LL)*Power(xj, 5LL) + 22680LL*Power(r, 6LL)*Power(xj, 6LL) +

            4860LL*Power(r, 7LL)*Power(xj, 7LL) + 810LL*Power(r, 8LL)*Power(xj, 8LL) +

            34LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  3LL*Power(xi, 14LL)*Power(xj, 8LL)*

                  (-19348875LL - 34827975LL*r*xj - 30958200LL*Power(r, 2LL)*Power(xj, 2LL) -

            18689580LL*Power(r, 3LL)*Power(xj, 3LL) -

            5847660LL*Power(r, 4LL)*Power(xj, 4LL) -

            3723300LL*Power(r, 5LL)*Power(xj, 5LL) -

            845040LL*Power(r, 6LL)*Power(xj, 6LL) - 58680LL*Power(r, 7LL)*Power(xj, 7LL) +

            1548LL*Power(r, 8LL)*Power(xj, 8LL) + 236LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  3LL*Power(xi, 8LL)*Power(xj, 14LL)*

                  (-593408025LL + 946053675LL*r*xj - 394427880LL*Power(r, 2LL)*Power(xj, 2LL) -

            315870660LL*Power(r, 3LL)*Power(xj, 3LL) -

            53891460LL*Power(r, 4LL)*Power(xj, 4LL) +

            910980LL*Power(r, 5LL)*Power(xj, 5LL) + 1409520LL*Power(r, 6LL)*Power(xj, 6LL) +

            192168LL*Power(r, 7LL)*Power(xj, 7LL) + 11196LL*Power(r, 8LL)*Power(xj, 8LL) +

            236LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (42525LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 15LL)*

                 Power(xi + xj, 15LL)) + (2LL*(42525LL*exp(2LL*r*(xi + xj))*

                                               Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                                               189LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                                               (-350LL*Power(r, 4LL)*Power(xi, 22LL) - 10LL*Power(r, 5LL)*Power(xi, 23LL) +

                                       225LL*Power(xj, 18LL) + 375LL*r*xi*Power(xj, 18LL) -

                                       70LL*Power(r, 3LL)*Power(xi, 21LL)*(75LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                       75LL*r*Power(xi, 3LL)*Power(xj, 16LL)*(-75LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                       75LL*Power(xi, 2LL)*Power(xj, 16LL)*(-45LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                       50LL*Power(r, 2LL)*Power(xi, 20LL)*(840LL + 71LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                       r*Power(xi, 9LL)*Power(xj, 10LL)*

                                       (4694625LL + 124800LL*Power(r, 2LL)*Power(xj, 2LL) -

              248LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       20LL*r*Power(xi, 17LL)*Power(xj, 2LL)*

                                       (-185895LL - 948LL*Power(r, 2LL)*Power(xj, 2LL) +

              2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       5LL*r*Power(xi, 5LL)*Power(xj, 14LL)*

                                       (7875LL - 450LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

                                       + 25LL*Power(xi, 4LL)*Power(xj, 14LL)*

                                       (945LL - 180LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                       375LL*Power(xi, 6LL)*Power(xj, 12LL)*

                                       (273LL - 84LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                       5LL*r*Power(xi, 11LL)*Power(xj, 8LL)*

                                       (-2803125LL + 49140LL*Power(r, 2LL)*Power(xj, 2LL) +

              8LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       5LL*r*Power(xi, 7LL)*Power(xj, 12LL)*

                                       (-16965LL + 5152LL*Power(r, 2LL)*Power(xj, 2LL) +

              14LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       325LL*Power(xi, 10LL)*Power(xj, 8LL)*

                                       (-60117LL - 5340LL*Power(r, 2LL)*Power(xj, 2LL) +

              40LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                       15LL*r*Power(xi, 15LL)*Power(xj, 4LL)*

                                       (845085LL - 22960LL*Power(r, 2LL)*Power(xj, 2LL) +

              52LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       15LL*r*Power(xi, 13LL)*Power(xj, 6LL)*

                                       (-139125LL - 10140LL*Power(r, 2LL)*Power(xj, 2LL) +

              52LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       75LL*Power(xi, 12LL)*Power(xj, 6LL)*

                                       (-729687LL + 25532LL*Power(r, 2LL)*Power(xj, 2LL) +

              52LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       60LL*Power(xi, 18LL)*(-5355LL - 11940LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                             86LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       2LL*r*Power(xi, 19LL)*(-89250LL - 35425LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                              124LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                       100LL*Power(xi, 16LL)*Power(xj, 2LL)*

                                       (-79713LL - 13311LL*Power(r, 2LL)*Power(xj, 2LL) +

              146LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                       5LL*Power(xi, 8LL)*Power(xj, 10LL)*

                                       (157365LL + 95940LL*Power(r, 2LL)*Power(xj, 2LL) +

              952LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                       15LL*Power(xi, 14LL)*Power(xj, 4LL)*

                                       (2638467LL - 157500LL*Power(r, 2LL)*Power(xj, 2LL) +

              1820LL*Power(r, 4LL)*Power(xj, 4LL))) +

                                               exp(2LL*r*xi)*Power(xi, 8LL)*

                                               (2LL*Power(xi, 2LL)*Power(xj, 20LL)*

                                       (1782492075LL + 1449175455LL*r*xj +

              533365560LL*Power(r, 2LL)*Power(xj, 2LL) +

              114631335LL*Power(r, 3LL)*Power(xj, 3LL) +

              15221115LL*Power(r, 4LL)*Power(xj, 4LL) +

              1142505LL*Power(r, 5LL)*Power(xj, 5LL) +

              18396LL*Power(r, 6LL)*Power(xj, 6LL) - 5238LL*Power(r, 7LL)*Power(xj, 7LL) -

              513LL*Power(r, 8LL)*Power(xj, 8LL) - 17LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       42LL*Power(xi, 4LL)*Power(xj, 18LL)*

                                       (251336925LL + 104824125LL*r*xj + 340200LL*Power(r, 2LL)*Power(xj, 2LL) -

              9122085LL*Power(r, 3LL)*Power(xj, 3LL) -

              2798145LL*Power(r, 4LL)*Power(xj, 4LL) -

              433755LL*Power(r, 5LL)*Power(xj, 5LL) -

              39060LL*Power(r, 6LL)*Power(xj, 6LL) - 1890LL*Power(r, 7LL)*Power(xj, 7LL) -

              27LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) +

                                       6LL*Power(xj, 22LL)*(34459425LL + 34459425LL*r*xj +

                                                            16216200LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                            4729725LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                            945945LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                            135135LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                            13860LL*Power(r, 6LL)*Power(xj, 6LL) + 990LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                            45LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) -

                                       3LL*Power(xi, 22LL)*(14175LL + 25515LL*r*xj +

                                                            22680LL*Power(r, 2LL)*Power(xj, 2LL) + 13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                            5670LL*Power(r, 4LL)*Power(xj, 4LL) + 1890LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                            504LL*Power(r, 6LL)*Power(xj, 6LL) + 108LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                            18LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       21LL*Power(xi, 18LL)*Power(xj, 4LL)*

                                       (212625LL + 382725LL*r*xj + 340200LL*Power(r, 2LL)*Power(xj, 2LL) +

              198450LL*Power(r, 3LL)*Power(xj, 3LL) +

              85050LL*Power(r, 4LL)*Power(xj, 4LL) + 28350LL*Power(r, 5LL)*Power(xj, 5LL) +

              7560LL*Power(r, 6LL)*Power(xj, 6LL) + 1836LL*Power(r, 7LL)*Power(xj, 7LL) +

              162LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       54LL*Power(xi, 6LL)*Power(xj, 16LL)*

                                       (133451955LL - 73700865LL*r*xj - 54096840LL*Power(r, 2LL)*Power(xj, 2LL) -

              8306235LL*Power(r, 3LL)*Power(xj, 3LL) +

              966945LL*Power(r, 4LL)*Power(xj, 4LL) +

              516747LL*Power(r, 5LL)*Power(xj, 5LL) +

              80724LL*Power(r, 6LL)*Power(xj, 6LL) + 6434LL*Power(r, 7LL)*Power(xj, 7LL) +

              251LL*Power(r, 8LL)*Power(xj, 8LL) + 3LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       315LL*Power(xi, 12LL)*Power(xj, 10LL)*

                                       (-405405LL - 710073LL*r*xj - 805896LL*Power(r, 2LL)*Power(xj, 2LL) -

              101556LL*Power(r, 3LL)*Power(xj, 3LL) -

              258804LL*Power(r, 4LL)*Power(xj, 4LL) -

              90972LL*Power(r, 5LL)*Power(xj, 5LL) - 9744LL*Power(r, 6LL)*Power(xj, 6LL) +

              120LL*Power(r, 7LL)*Power(xj, 7LL) + 84LL*Power(r, 8LL)*Power(xj, 8LL) +

              4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       315LL*Power(xi, 10LL)*Power(xj, 12LL)*

                                       (-482895LL - 2656395LL*r*xj + 1186920LL*Power(r, 2LL)*Power(xj, 2LL) -

              1155420LL*Power(r, 3LL)*Power(xj, 3LL) -

              643356LL*Power(r, 4LL)*Power(xj, 4LL) -

              93492LL*Power(r, 5LL)*Power(xj, 5LL) + 336LL*Power(r, 6LL)*Power(xj, 6LL) +

              1368LL*Power(r, 7LL)*Power(xj, 7LL) + 132LL*Power(r, 8LL)*Power(xj, 8LL) +

              4LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       27LL*Power(xi, 16LL)*Power(xj, 6LL)*

                                       (-716625LL - 1289925LL*r*xj - 1146600LL*Power(r, 2LL)*Power(xj, 2LL) -

              668850LL*Power(r, 3LL)*Power(xj, 3LL) -

              286650LL*Power(r, 4LL)*Power(xj, 4LL) -

              90006LL*Power(r, 5LL)*Power(xj, 5LL) - 32872LL*Power(r, 6LL)*Power(xj, 6LL) -

              4812LL*Power(r, 7LL)*Power(xj, 7LL) - 178LL*Power(r, 8LL)*Power(xj, 8LL) +

              6LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       Power(xi, 20LL)*Power(xj, 2LL)*

                                       (637875LL + 1148175LL*r*xj + 1020600LL*Power(r, 2LL)*Power(xj, 2LL) +

              595350LL*Power(r, 3LL)*Power(xj, 3LL) +

              255150LL*Power(r, 4LL)*Power(xj, 4LL) +

              85050LL*Power(r, 5LL)*Power(xj, 5LL) + 22680LL*Power(r, 6LL)*Power(xj, 6LL) +

              4860LL*Power(r, 7LL)*Power(xj, 7LL) + 810LL*Power(r, 8LL)*Power(xj, 8LL) +

              34LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       3LL*Power(xi, 14LL)*Power(xj, 8LL)*

                                       (-19348875LL - 34827975LL*r*xj - 30958200LL*Power(r, 2LL)*Power(xj, 2LL) -

              18689580LL*Power(r, 3LL)*Power(xj, 3LL) -

              5847660LL*Power(r, 4LL)*Power(xj, 4LL) -

              3723300LL*Power(r, 5LL)*Power(xj, 5LL) -

              845040LL*Power(r, 6LL)*Power(xj, 6LL) -

              58680LL*Power(r, 7LL)*Power(xj, 7LL) + 1548LL*Power(r, 8LL)*Power(xj, 8LL) +

              236LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       3LL*Power(xi, 8LL)*Power(xj, 14LL)*

                                       (-593408025LL + 946053675LL*r*xj -

              394427880LL*Power(r, 2LL)*Power(xj, 2LL) -

              315870660LL*Power(r, 3LL)*Power(xj, 3LL) -

              53891460LL*Power(r, 4LL)*Power(xj, 4LL) +

              910980LL*Power(r, 5LL)*Power(xj, 5LL) +

              1409520LL*Power(r, 6LL)*Power(xj, 6LL) +

              192168LL*Power(r, 7LL)*Power(xj, 7LL) + 11196LL*Power(r, 8LL)*Power(xj, 8LL) +

              236LL*Power(r, 9LL)*Power(xj, 9LL)))))/

                (42525LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 15LL)*Power(xi + xj, 14LL)) -

                (85050LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                 189LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-1400LL*Power(r, 3LL)*Power(xi, 22LL) - 50LL*Power(r, 4LL)*Power(xi, 23LL) -

        7100LL*Power(r, 3LL)*Power(xi, 20LL)*Power(xj, 2LL) -

        140LL*Power(r, 4LL)*Power(xi, 21LL)*Power(xj, 2LL) + 375LL*xi*Power(xj, 18LL) +

        600LL*r*Power(xi, 2LL)*Power(xj, 18LL) +

        300LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 18LL) -

        210LL*Power(r, 2LL)*Power(xi, 21LL)*(75LL + Power(r, 2LL)*Power(xj, 2LL)) +

        75LL*Power(xi, 3LL)*Power(xj, 16LL)*(-75LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        100LL*r*Power(xi, 20LL)*(840LL + 71LL*Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 9LL)*Power(xj, 10LL)*

        (249600LL*r*Power(xj, 2LL) - 992LL*Power(r, 3LL)*Power(xj, 4LL)) +

        20LL*r*Power(xi, 17LL)*Power(xj, 2LL)*

        (-1896LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        5LL*r*Power(xi, 5LL)*Power(xj, 14LL)*

        (-900LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        25LL*Power(xi, 4LL)*Power(xj, 14LL)*

        (-360LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) -

        375LL*Power(xi, 6LL)*Power(xj, 12LL)*

        (-168LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) -

        5LL*r*Power(xi, 11LL)*Power(xj, 8LL)*

        (98280LL*r*Power(xj, 2LL) + 32LL*Power(r, 3LL)*Power(xj, 4LL)) +

        5LL*r*Power(xi, 7LL)*Power(xj, 12LL)*

        (10304LL*r*Power(xj, 2LL) + 56LL*Power(r, 3LL)*Power(xj, 4LL)) +

        325LL*Power(xi, 10LL)*Power(xj, 8LL)*

        (-10680LL*r*Power(xj, 2LL) + 160LL*Power(r, 3LL)*Power(xj, 4LL)) -

        15LL*r*Power(xi, 15LL)*Power(xj, 4LL)*

        (-45920LL*r*Power(xj, 2LL) + 208LL*Power(r, 3LL)*Power(xj, 4LL)) +

        15LL*r*Power(xi, 13LL)*Power(xj, 6LL)*

        (-20280LL*r*Power(xj, 2LL) + 208LL*Power(r, 3LL)*Power(xj, 4LL)) +

        75LL*Power(xi, 12LL)*Power(xj, 6LL)*

        (51064LL*r*Power(xj, 2LL) + 208LL*Power(r, 3LL)*Power(xj, 4LL)) +

        60LL*Power(xi, 18LL)*(-23880LL*r*Power(xj, 2LL) +

                              344LL*Power(r, 3LL)*Power(xj, 4LL)) +

        2LL*r*Power(xi, 19LL)*(-70850LL*r*Power(xj, 2LL) +

                               496LL*Power(r, 3LL)*Power(xj, 4LL)) +

        100LL*Power(xi, 16LL)*Power(xj, 2LL)*

        (-26622LL*r*Power(xj, 2LL) + 584LL*Power(r, 3LL)*Power(xj, 4LL)) -

        5LL*Power(xi, 8LL)*Power(xj, 10LL)*

        (191880LL*r*Power(xj, 2LL) + 3808LL*Power(r, 3LL)*Power(xj, 4LL)) -

        15LL*Power(xi, 14LL)*Power(xj, 4LL)*

        (-315000LL*r*Power(xj, 2LL) + 7280LL*Power(r, 3LL)*Power(xj, 4LL)) +

        Power(xi, 9LL)*Power(xj, 10LL)*

        (4694625LL + 124800LL*Power(r, 2LL)*Power(xj, 2LL) -

            248LL*Power(r, 4LL)*Power(xj, 4LL)) +

        20LL*Power(xi, 17LL)*Power(xj, 2LL)*

        (-185895LL - 948LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

        + 5LL*Power(xi, 5LL)*Power(xj, 14LL)*(7875LL - 450LL*Power(r, 2LL)*Power(xj, 2LL) +

                                              2LL*Power(r, 4LL)*Power(xj, 4LL)) -

        5LL*Power(xi, 11LL)*Power(xj, 8LL)*

        (-2803125LL + 49140LL*Power(r, 2LL)*Power(xj, 2LL) +

            8LL*Power(r, 4LL)*Power(xj, 4LL)) +

        5LL*Power(xi, 7LL)*Power(xj, 12LL)*

        (-16965LL + 5152LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

        - 15LL*Power(xi, 15LL)*Power(xj, 4LL)*(845085LL - 22960LL*Power(r, 2LL)*Power(xj, 2LL) +

                                               52LL*Power(r, 4LL)*Power(xj, 4LL)) +

        15LL*Power(xi, 13LL)*Power(xj, 6LL)*

        (-139125LL - 10140LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

        2LL*Power(xi, 19LL)*(-89250LL - 35425LL*Power(r, 2LL)*Power(xj, 2LL) +

                             124LL*Power(r, 4LL)*Power(xj, 4LL))) +

                 378LL*exp(2LL*r*xj)*Power(xj, 13LL)*

                 (-350LL*Power(r, 4LL)*Power(xi, 22LL) - 10LL*Power(r, 5LL)*Power(xi, 23LL) +

        225LL*Power(xj, 18LL) + 375LL*r*xi*Power(xj, 18LL) -

        70LL*Power(r, 3LL)*Power(xi, 21LL)*(75LL + Power(r, 2LL)*Power(xj, 2LL)) +

        75LL*r*Power(xi, 3LL)*Power(xj, 16LL)*(-75LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        75LL*Power(xi, 2LL)*Power(xj, 16LL)*(-45LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

        50LL*Power(r, 2LL)*Power(xi, 20LL)*(840LL + 71LL*Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 9LL)*Power(xj, 10LL)*

        (4694625LL + 124800LL*Power(r, 2LL)*Power(xj, 2LL) -

            248LL*Power(r, 4LL)*Power(xj, 4LL)) +

        20LL*r*Power(xi, 17LL)*Power(xj, 2LL)*

        (-185895LL - 948LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

        + 5LL*r*Power(xi, 5LL)*Power(xj, 14LL)*

        (7875LL - 450LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        25LL*Power(xi, 4LL)*Power(xj, 14LL)*

        (945LL - 180LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

        375LL*Power(xi, 6LL)*Power(xj, 12LL)*

        (273LL - 84LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

        5LL*r*Power(xi, 11LL)*Power(xj, 8LL)*

        (-2803125LL + 49140LL*Power(r, 2LL)*Power(xj, 2LL) +

            8LL*Power(r, 4LL)*Power(xj, 4LL)) +

        5LL*r*Power(xi, 7LL)*Power(xj, 12LL)*

        (-16965LL + 5152LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

        + 325LL*Power(xi, 10LL)*Power(xj, 8LL)*

        (-60117LL - 5340LL*Power(r, 2LL)*Power(xj, 2LL) + 40LL*Power(r, 4LL)*Power(xj, 4LL))

        - 15LL*r*Power(xi, 15LL)*Power(xj, 4LL)*

        (845085LL - 22960LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

        15LL*r*Power(xi, 13LL)*Power(xj, 6LL)*

        (-139125LL - 10140LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

        75LL*Power(xi, 12LL)*Power(xj, 6LL)*

        (-729687LL + 25532LL*Power(r, 2LL)*Power(xj, 2LL) +

            52LL*Power(r, 4LL)*Power(xj, 4LL)) +

        60LL*Power(xi, 18LL)*(-5355LL - 11940LL*Power(r, 2LL)*Power(xj, 2LL) +

                              86LL*Power(r, 4LL)*Power(xj, 4LL)) +

        2LL*r*Power(xi, 19LL)*(-89250LL - 35425LL*Power(r, 2LL)*Power(xj, 2LL) +

                               124LL*Power(r, 4LL)*Power(xj, 4LL)) +

        100LL*Power(xi, 16LL)*Power(xj, 2LL)*

        (-79713LL - 13311LL*Power(r, 2LL)*Power(xj, 2LL) +

            146LL*Power(r, 4LL)*Power(xj, 4LL)) -

        5LL*Power(xi, 8LL)*Power(xj, 10LL)*

        (157365LL + 95940LL*Power(r, 2LL)*Power(xj, 2LL) +

            952LL*Power(r, 4LL)*Power(xj, 4LL)) -

        15LL*Power(xi, 14LL)*Power(xj, 4LL)*

        (2638467LL - 157500LL*Power(r, 2LL)*Power(xj, 2LL) +

            1820LL*Power(r, 4LL)*Power(xj, 4LL))) +

                 exp(2LL*r*xi)*Power(xi, 8LL)*

                 (2LL*Power(xi, 2LL)*Power(xj, 20LL)*

        (1449175455LL*xj + 1066731120LL*r*Power(xj, 2LL) +

            343894005LL*Power(r, 2LL)*Power(xj, 3LL) +

            60884460LL*Power(r, 3LL)*Power(xj, 4LL) +

            5712525LL*Power(r, 4LL)*Power(xj, 5LL) +

            110376LL*Power(r, 5LL)*Power(xj, 6LL) - 36666LL*Power(r, 6LL)*Power(xj, 7LL) -

            4104LL*Power(r, 7LL)*Power(xj, 8LL) - 153LL*Power(r, 8LL)*Power(xj, 9LL)) +

        42LL*Power(xi, 4LL)*Power(xj, 18LL)*

        (104824125LL*xj + 680400LL*r*Power(xj, 2LL) -

            27366255LL*Power(r, 2LL)*Power(xj, 3LL) -

            11192580LL*Power(r, 3LL)*Power(xj, 4LL) -

            2168775LL*Power(r, 4LL)*Power(xj, 5LL) -

            234360LL*Power(r, 5LL)*Power(xj, 6LL) - 13230LL*Power(r, 6LL)*Power(xj, 7LL) -

            216LL*Power(r, 7LL)*Power(xj, 8LL) + 9LL*Power(r, 8LL)*Power(xj, 9LL)) +

        6LL*Power(xj, 22LL)*(34459425LL*xj + 32432400LL*r*Power(xj, 2LL) +

                             14189175LL*Power(r, 2LL)*Power(xj, 3LL) +

                             3783780LL*Power(r, 3LL)*Power(xj, 4LL) +

                             675675LL*Power(r, 4LL)*Power(xj, 5LL) + 83160LL*Power(r, 5LL)*Power(xj, 6LL) +

                             6930LL*Power(r, 6LL)*Power(xj, 7LL) + 360LL*Power(r, 7LL)*Power(xj, 8LL) +

                             9LL*Power(r, 8LL)*Power(xj, 9LL)) -

        3LL*Power(xi, 22LL)*(25515LL*xj + 45360LL*r*Power(xj, 2LL) +

                             39690LL*Power(r, 2LL)*Power(xj, 3LL) + 22680LL*Power(r, 3LL)*Power(xj, 4LL) +

                             9450LL*Power(r, 4LL)*Power(xj, 5LL) + 3024LL*Power(r, 5LL)*Power(xj, 6LL) +

                             756LL*Power(r, 6LL)*Power(xj, 7LL) + 144LL*Power(r, 7LL)*Power(xj, 8LL) +

                             18LL*Power(r, 8LL)*Power(xj, 9LL)) -

        21LL*Power(xi, 18LL)*Power(xj, 4LL)*

        (382725LL*xj + 680400LL*r*Power(xj, 2LL) +

            595350LL*Power(r, 2LL)*Power(xj, 3LL) + 340200LL*Power(r, 3LL)*Power(xj, 4LL) +

            141750LL*Power(r, 4LL)*Power(xj, 5LL) + 45360LL*Power(r, 5LL)*Power(xj, 6LL) +

            12852LL*Power(r, 6LL)*Power(xj, 7LL) + 1296LL*Power(r, 7LL)*Power(xj, 8LL) +

            18LL*Power(r, 8LL)*Power(xj, 9LL)) +

        54LL*Power(xi, 6LL)*Power(xj, 16LL)*

        (-73700865LL*xj - 108193680LL*r*Power(xj, 2LL) -

            24918705LL*Power(r, 2LL)*Power(xj, 3LL) +

            3867780LL*Power(r, 3LL)*Power(xj, 4LL) +

            2583735LL*Power(r, 4LL)*Power(xj, 5LL) +

            484344LL*Power(r, 5LL)*Power(xj, 6LL) + 45038LL*Power(r, 6LL)*Power(xj, 7LL) +

            2008LL*Power(r, 7LL)*Power(xj, 8LL) + 27LL*Power(r, 8LL)*Power(xj, 9LL)) -

        315LL*Power(xi, 12LL)*Power(xj, 10LL)*

        (-710073LL*xj - 1611792LL*r*Power(xj, 2LL) -

            304668LL*Power(r, 2LL)*Power(xj, 3LL) -

            1035216LL*Power(r, 3LL)*Power(xj, 4LL) -

            454860LL*Power(r, 4LL)*Power(xj, 5LL) - 58464LL*Power(r, 5LL)*Power(xj, 6LL) +

            840LL*Power(r, 6LL)*Power(xj, 7LL) + 672LL*Power(r, 7LL)*Power(xj, 8LL) +

            36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        315LL*Power(xi, 10LL)*Power(xj, 12LL)*

        (-2656395LL*xj + 2373840LL*r*Power(xj, 2LL) -

            3466260LL*Power(r, 2LL)*Power(xj, 3LL) -

            2573424LL*Power(r, 3LL)*Power(xj, 4LL) -

            467460LL*Power(r, 4LL)*Power(xj, 5LL) + 2016LL*Power(r, 5LL)*Power(xj, 6LL) +

            9576LL*Power(r, 6LL)*Power(xj, 7LL) + 1056LL*Power(r, 7LL)*Power(xj, 8LL) +

            36LL*Power(r, 8LL)*Power(xj, 9LL)) -

        27LL*Power(xi, 16LL)*Power(xj, 6LL)*

        (-1289925LL*xj - 2293200LL*r*Power(xj, 2LL) -

            2006550LL*Power(r, 2LL)*Power(xj, 3LL) -

            1146600LL*Power(r, 3LL)*Power(xj, 4LL) -

            450030LL*Power(r, 4LL)*Power(xj, 5LL) - 197232LL*Power(r, 5LL)*Power(xj, 6LL) -

            33684LL*Power(r, 6LL)*Power(xj, 7LL) - 1424LL*Power(r, 7LL)*Power(xj, 8LL) +

            54LL*Power(r, 8LL)*Power(xj, 9LL)) +

        Power(xi, 20LL)*Power(xj, 2LL)*

        (1148175LL*xj + 2041200LL*r*Power(xj, 2LL) +

            1786050LL*Power(r, 2LL)*Power(xj, 3LL) +

            1020600LL*Power(r, 3LL)*Power(xj, 4LL) +

            425250LL*Power(r, 4LL)*Power(xj, 5LL) + 136080LL*Power(r, 5LL)*Power(xj, 6LL) +

            34020LL*Power(r, 6LL)*Power(xj, 7LL) + 6480LL*Power(r, 7LL)*Power(xj, 8LL) +

            306LL*Power(r, 8LL)*Power(xj, 9LL)) +

        3LL*Power(xi, 14LL)*Power(xj, 8LL)*

        (-34827975LL*xj - 61916400LL*r*Power(xj, 2LL) -

            56068740LL*Power(r, 2LL)*Power(xj, 3LL) -

            23390640LL*Power(r, 3LL)*Power(xj, 4LL) -

            18616500LL*Power(r, 4LL)*Power(xj, 5LL) -

            5070240LL*Power(r, 5LL)*Power(xj, 6LL) -

            410760LL*Power(r, 6LL)*Power(xj, 7LL) + 12384LL*Power(r, 7LL)*Power(xj, 8LL) +

            2124LL*Power(r, 8LL)*Power(xj, 9LL)) -

        3LL*Power(xi, 8LL)*Power(xj, 14LL)*

        (946053675LL*xj - 788855760LL*r*Power(xj, 2LL) -

            947611980LL*Power(r, 2LL)*Power(xj, 3LL) -

            215565840LL*Power(r, 3LL)*Power(xj, 4LL) +

            4554900LL*Power(r, 4LL)*Power(xj, 5LL) +

            8457120LL*Power(r, 5LL)*Power(xj, 6LL) +

            1345176LL*Power(r, 6LL)*Power(xj, 7LL) + 89568LL*Power(r, 7LL)*Power(xj, 8LL) +

            2124LL*Power(r, 8LL)*Power(xj, 9LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 9LL)*

                 (2LL*Power(xi, 2LL)*Power(xj, 20LL)*

        (1782492075LL + 1449175455LL*r*xj + 533365560LL*Power(r, 2LL)*Power(xj, 2LL) +

            114631335LL*Power(r, 3LL)*Power(xj, 3LL) +

            15221115LL*Power(r, 4LL)*Power(xj, 4LL) +

            1142505LL*Power(r, 5LL)*Power(xj, 5LL) + 18396LL*Power(r, 6LL)*Power(xj, 6LL) -

            5238LL*Power(r, 7LL)*Power(xj, 7LL) - 513LL*Power(r, 8LL)*Power(xj, 8LL) -

            17LL*Power(r, 9LL)*Power(xj, 9LL)) +

        42LL*Power(xi, 4LL)*Power(xj, 18LL)*

        (251336925LL + 104824125LL*r*xj + 340200LL*Power(r, 2LL)*Power(xj, 2LL) -

            9122085LL*Power(r, 3LL)*Power(xj, 3LL) -

            2798145LL*Power(r, 4LL)*Power(xj, 4LL) - 433755LL*Power(r, 5LL)*Power(xj, 5LL) -

            39060LL*Power(r, 6LL)*Power(xj, 6LL) - 1890LL*Power(r, 7LL)*Power(xj, 7LL) -

            27LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) +

        6LL*Power(xj, 22LL)*(34459425LL + 34459425LL*r*xj +

                             16216200LL*Power(r, 2LL)*Power(xj, 2LL) +

                             4729725LL*Power(r, 3LL)*Power(xj, 3LL) + 945945LL*Power(r, 4LL)*Power(xj, 4LL) +

                             135135LL*Power(r, 5LL)*Power(xj, 5LL) + 13860LL*Power(r, 6LL)*Power(xj, 6LL) +

                             990LL*Power(r, 7LL)*Power(xj, 7LL) + 45LL*Power(r, 8LL)*Power(xj, 8LL) +

                             Power(r, 9LL)*Power(xj, 9LL)) -

        3LL*Power(xi, 22LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                             13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                             1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                             108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                             2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        21LL*Power(xi, 18LL)*Power(xj, 4LL)*

        (212625LL + 382725LL*r*xj + 340200LL*Power(r, 2LL)*Power(xj, 2LL) +

            198450LL*Power(r, 3LL)*Power(xj, 3LL) + 85050LL*Power(r, 4LL)*Power(xj, 4LL) +

            28350LL*Power(r, 5LL)*Power(xj, 5LL) + 7560LL*Power(r, 6LL)*Power(xj, 6LL) +

            1836LL*Power(r, 7LL)*Power(xj, 7LL) + 162LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

        54LL*Power(xi, 6LL)*Power(xj, 16LL)*

        (133451955LL - 73700865LL*r*xj - 54096840LL*Power(r, 2LL)*Power(xj, 2LL) -

            8306235LL*Power(r, 3LL)*Power(xj, 3LL) + 966945LL*Power(r, 4LL)*Power(xj, 4LL) +

            516747LL*Power(r, 5LL)*Power(xj, 5LL) + 80724LL*Power(r, 6LL)*Power(xj, 6LL) +

            6434LL*Power(r, 7LL)*Power(xj, 7LL) + 251LL*Power(r, 8LL)*Power(xj, 8LL) +

            3LL*Power(r, 9LL)*Power(xj, 9LL)) -

        315LL*Power(xi, 12LL)*Power(xj, 10LL)*

        (-405405LL - 710073LL*r*xj - 805896LL*Power(r, 2LL)*Power(xj, 2LL) -

            101556LL*Power(r, 3LL)*Power(xj, 3LL) - 258804LL*Power(r, 4LL)*Power(xj, 4LL) -

            90972LL*Power(r, 5LL)*Power(xj, 5LL) - 9744LL*Power(r, 6LL)*Power(xj, 6LL) +

            120LL*Power(r, 7LL)*Power(xj, 7LL) + 84LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        315LL*Power(xi, 10LL)*Power(xj, 12LL)*

        (-482895LL - 2656395LL*r*xj + 1186920LL*Power(r, 2LL)*Power(xj, 2LL) -

            1155420LL*Power(r, 3LL)*Power(xj, 3LL) - 643356LL*Power(r, 4LL)*Power(xj, 4LL) -

            93492LL*Power(r, 5LL)*Power(xj, 5LL) + 336LL*Power(r, 6LL)*Power(xj, 6LL) +

            1368LL*Power(r, 7LL)*Power(xj, 7LL) + 132LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) -

        27LL*Power(xi, 16LL)*Power(xj, 6LL)*

        (-716625LL - 1289925LL*r*xj - 1146600LL*Power(r, 2LL)*Power(xj, 2LL) -

            668850LL*Power(r, 3LL)*Power(xj, 3LL) - 286650LL*Power(r, 4LL)*Power(xj, 4LL) -

            90006LL*Power(r, 5LL)*Power(xj, 5LL) - 32872LL*Power(r, 6LL)*Power(xj, 6LL) -

            4812LL*Power(r, 7LL)*Power(xj, 7LL) - 178LL*Power(r, 8LL)*Power(xj, 8LL) +

            6LL*Power(r, 9LL)*Power(xj, 9LL)) +

        Power(xi, 20LL)*Power(xj, 2LL)*

        (637875LL + 1148175LL*r*xj + 1020600LL*Power(r, 2LL)*Power(xj, 2LL) +

            595350LL*Power(r, 3LL)*Power(xj, 3LL) + 255150LL*Power(r, 4LL)*Power(xj, 4LL) +

            85050LL*Power(r, 5LL)*Power(xj, 5LL) + 22680LL*Power(r, 6LL)*Power(xj, 6LL) +

            4860LL*Power(r, 7LL)*Power(xj, 7LL) + 810LL*Power(r, 8LL)*Power(xj, 8LL) +

            34LL*Power(r, 9LL)*Power(xj, 9LL)) +

        3LL*Power(xi, 14LL)*Power(xj, 8LL)*

        (-19348875LL - 34827975LL*r*xj - 30958200LL*Power(r, 2LL)*Power(xj, 2LL) -

            18689580LL*Power(r, 3LL)*Power(xj, 3LL) -

            5847660LL*Power(r, 4LL)*Power(xj, 4LL) -

            3723300LL*Power(r, 5LL)*Power(xj, 5LL) - 845040LL*Power(r, 6LL)*Power(xj, 6LL) -

            58680LL*Power(r, 7LL)*Power(xj, 7LL) + 1548LL*Power(r, 8LL)*Power(xj, 8LL) +

            236LL*Power(r, 9LL)*Power(xj, 9LL)) -

        3LL*Power(xi, 8LL)*Power(xj, 14LL)*

        (-593408025LL + 946053675LL*r*xj - 394427880LL*Power(r, 2LL)*Power(xj, 2LL) -

            315870660LL*Power(r, 3LL)*Power(xj, 3LL) -

            53891460LL*Power(r, 4LL)*Power(xj, 4LL) + 910980LL*Power(r, 5LL)*Power(xj, 5LL) +

            1409520LL*Power(r, 6LL)*Power(xj, 6LL) + 192168LL*Power(r, 7LL)*Power(xj, 7LL) +

            11196LL*Power(r, 8LL)*Power(xj, 8LL) + 236LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (42525LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 15LL)*Power(xi + xj, 15LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_5S_3S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_3S_5S(r, xj, xi);
}
