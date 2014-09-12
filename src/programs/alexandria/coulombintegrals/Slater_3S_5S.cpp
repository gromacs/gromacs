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

cl_R Slater_3S_5S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (31059LL*xi)/163840LL

            ;
        }
        else
        {
            S = (1LL/r)*((-313841848320000LL + 313841848320000LL*exp(2LL*rxi) - 568188982486125LL*rxi -

                          508694268332250LL*Power(rxi, 2LL) - 299892470377500LL*Power(rxi, 3LL) -

                          130753815192000LL*Power(rxi, 4LL) - 44881155118800LL*Power(rxi, 5LL) -

                          12601803614400LL*Power(rxi, 6LL) - 2967953788800LL*Power(rxi, 7LL) -

                          596237241600LL*Power(rxi, 8LL) - 103264761600LL*Power(rxi, 9LL) -

                          15498362880LL*Power(rxi, 10LL) - 2012774400LL*Power(rxi, 11LL) -

                          223641600LL*Power(rxi, 12LL) - 20643840LL*Power(rxi, 13LL) -

                          1474560LL*Power(rxi, 14LL) - 65536LL*Power(rxi, 15LL))/

                         (3.1384184832e14*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(3LL*Power(xi, 14LL) + 45LL*Power(xi, 13LL)*xj + 315LL*Power(xi, 12LL)*Power(xj, 2LL) +

                        1365LL*Power(xi, 11LL)*Power(xj, 3LL) + 4095LL*Power(xi, 10LL)*Power(xj, 4LL) +

                        9009LL*Power(xi, 9LL)*Power(xj, 5LL) + 15015LL*Power(xi, 8LL)*Power(xj, 6LL) +

                        19305LL*Power(xi, 7LL)*Power(xj, 7LL) + 19305LL*Power(xi, 6LL)*Power(xj, 8LL) +

                        15015LL*Power(xi, 5LL)*Power(xj, 9LL) + 6825LL*Power(xi, 4LL)*Power(xj, 10LL) +

                        2275LL*Power(xi, 3LL)*Power(xj, 11LL) + 525LL*Power(xi, 2LL)*Power(xj, 12LL) +

                        75LL*xi*Power(xj, 13LL) + 5LL*Power(xj, 14LL)))/(15LL*Power(xi + xj, 15LL))

            ;
        }
        else
        {
            S = (1LL/r)*((42525LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 15LL) +

                          189LL*exp(2LL*rxj)*Power(rxj, 12LL)*

                          (-350LL*Power(rxi, 22LL) - 10LL*Power(rxi, 23LL) + 225LL*Power(rxj, 18LL) +

                           375LL*rxi*Power(rxj, 18LL) - 70LL*Power(rxi, 21LL)*(75LL + Power(rxj, 2LL)) +

                           75LL*Power(rxi, 3LL)*Power(rxj, 16LL)*(-75LL + 2LL*Power(rxj, 2LL)) +

                           75LL*Power(rxi, 2LL)*Power(rxj, 16LL)*(-45LL + 4LL*Power(rxj, 2LL)) -

                           50LL*Power(rxi, 20LL)*(840LL + 71LL*Power(rxj, 2LL)) +

                           Power(rxi, 9LL)*Power(rxj, 10LL)*

                           (4694625LL + 124800LL*Power(rxj, 2LL) - 248LL*Power(rxj, 4LL)) +

                           20LL*Power(rxi, 17LL)*Power(rxj, 2LL)*

                           (-185895LL - 948LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           5LL*Power(rxi, 5LL)*Power(rxj, 14LL)*

                           (7875LL - 450LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           25LL*Power(rxi, 4LL)*Power(rxj, 14LL)*

                           (945LL - 180LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) -

                           375LL*Power(rxi, 6LL)*Power(rxj, 12LL)*

                           (273LL - 84LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) -

                           5LL*Power(rxi, 11LL)*Power(rxj, 8LL)*

                           (-2803125LL + 49140LL*Power(rxj, 2LL) + 8LL*Power(rxj, 4LL)) +

                           5LL*Power(rxi, 7LL)*Power(rxj, 12LL)*

                           (-16965LL + 5152LL*Power(rxj, 2LL) + 14LL*Power(rxj, 4LL)) +

                           325LL*Power(rxi, 10LL)*Power(rxj, 8LL)*

                           (-60117LL - 5340LL*Power(rxj, 2LL) + 40LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 15LL)*Power(rxj, 4LL)*

                           (845085LL - 22960LL*Power(rxj, 2LL) + 52LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 13LL)*Power(rxj, 6LL)*

                           (-139125LL - 10140LL*Power(rxj, 2LL) + 52LL*Power(rxj, 4LL)) +

                           75LL*Power(rxi, 12LL)*Power(rxj, 6LL)*

                           (-729687LL + 25532LL*Power(rxj, 2LL) + 52LL*Power(rxj, 4LL)) +

                           60LL*Power(rxi, 18LL)*(-5355LL - 11940LL*Power(rxj, 2LL) + 86LL*Power(rxj, 4LL)) +

                           2LL*Power(rxi, 19LL)*(-89250LL - 35425LL*Power(rxj, 2LL) + 124LL*Power(rxj, 4LL)) +

                           100LL*Power(rxi, 16LL)*Power(rxj, 2LL)*

                           (-79713LL - 13311LL*Power(rxj, 2LL) + 146LL*Power(rxj, 4LL)) -

                           5LL*Power(rxi, 8LL)*Power(rxj, 10LL)*

                           (157365LL + 95940LL*Power(rxj, 2LL) + 952LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 14LL)*Power(rxj, 4LL)*

                           (2638467LL - 157500LL*Power(rxj, 2LL) + 1820LL*Power(rxj, 4LL))) -

                          exp(2LL*rxi)*Power(rxi, 8LL)*

                          (3LL*Power(rxi, 14LL)*Power(rxj, 8LL)*

                           (19348875LL + 34827975LL*rxj + 30958200LL*Power(rxj, 2LL) +

           18689580LL*Power(rxj, 3LL) + 5847660LL*Power(rxj, 4LL) +

           3723300LL*Power(rxj, 5LL) + 845040LL*Power(rxj, 6LL) + 58680LL*Power(rxj, 7LL) -

           1548LL*Power(rxj, 8LL) - 236LL*Power(rxj, 9LL)) -

                           42LL*Power(rxi, 4LL)*Power(rxj, 18LL)*

                           (251336925LL + 104824125LL*rxj + 340200LL*Power(rxj, 2LL) -

           9122085LL*Power(rxj, 3LL) - 2798145LL*Power(rxj, 4LL) -

           433755LL*Power(rxj, 5LL) - 39060LL*Power(rxj, 6LL) - 1890LL*Power(rxj, 7LL) -

           27LL*Power(rxj, 8LL) + Power(rxj, 9LL)) -

                           6LL*Power(rxj, 22LL)*(34459425LL + 34459425LL*rxj + 16216200LL*Power(rxj, 2LL) +

                                                 4729725LL*Power(rxj, 3LL) + 945945LL*Power(rxj, 4LL) +

                                                 135135LL*Power(rxj, 5LL) + 13860LL*Power(rxj, 6LL) + 990LL*Power(rxj, 7LL) +

                                                 45LL*Power(rxj, 8LL) + Power(rxj, 9LL)) +

                           3LL*Power(rxi, 22LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj, 2LL) +

                                                 13230LL*Power(rxj, 3LL) + 5670LL*Power(rxj, 4LL) + 1890LL*Power(rxj, 5LL) +

                                                 504LL*Power(rxj, 6LL) + 108LL*Power(rxj, 7LL) + 18LL*Power(rxj, 8LL) +

                                                 2LL*Power(rxj, 9LL)) + 21LL*Power(rxi, 18LL)*Power(rxj, 4LL)*

                           (212625LL + 382725LL*rxj + 340200LL*Power(rxj, 2LL) + 198450LL*Power(rxj, 3LL) +

           85050LL*Power(rxj, 4LL) + 28350LL*Power(rxj, 5LL) + 7560LL*Power(rxj, 6LL) +

           1836LL*Power(rxj, 7LL) + 162LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) -

                           54LL*Power(rxi, 6LL)*Power(rxj, 16LL)*

                           (133451955LL - 73700865LL*rxj - 54096840LL*Power(rxj, 2LL) -

           8306235LL*Power(rxj, 3LL) + 966945LL*Power(rxj, 4LL) +

           516747LL*Power(rxj, 5LL) + 80724LL*Power(rxj, 6LL) + 6434LL*Power(rxj, 7LL) +

           251LL*Power(rxj, 8LL) + 3LL*Power(rxj, 9LL)) +

                           315LL*Power(rxi, 12LL)*Power(rxj, 10LL)*

                           (-405405LL - 710073LL*rxj - 805896LL*Power(rxj, 2LL) - 101556LL*Power(rxj, 3LL) -

           258804LL*Power(rxj, 4LL) - 90972LL*Power(rxj, 5LL) - 9744LL*Power(rxj, 6LL) +

           120LL*Power(rxj, 7LL) + 84LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) -

                           315LL*Power(rxi, 10LL)*Power(rxj, 12LL)*

                           (-482895LL - 2656395LL*rxj + 1186920LL*Power(rxj, 2LL) -

           1155420LL*Power(rxj, 3LL) - 643356LL*Power(rxj, 4LL) - 93492LL*Power(rxj, 5LL) +

           336LL*Power(rxj, 6LL) + 1368LL*Power(rxj, 7LL) + 132LL*Power(rxj, 8LL) +

           4LL*Power(rxj, 9LL)) + 27LL*Power(rxi, 16LL)*Power(rxj, 6LL)*

                           (-716625LL - 1289925LL*rxj - 1146600LL*Power(rxj, 2LL) -

           668850LL*Power(rxj, 3LL) - 286650LL*Power(rxj, 4LL) - 90006LL*Power(rxj, 5LL) -

           32872LL*Power(rxj, 6LL) - 4812LL*Power(rxj, 7LL) - 178LL*Power(rxj, 8LL) +

           6LL*Power(rxj, 9LL)) + 2LL*Power(rxi, 2LL)*Power(rxj, 20LL)*

                           (-1782492075LL - 1449175455LL*rxj - 533365560LL*Power(rxj, 2LL) -

           114631335LL*Power(rxj, 3LL) - 15221115LL*Power(rxj, 4LL) -

           1142505LL*Power(rxj, 5LL) - 18396LL*Power(rxj, 6LL) + 5238LL*Power(rxj, 7LL) +

           513LL*Power(rxj, 8LL) + 17LL*Power(rxj, 9LL)) -

                           Power(rxi, 20LL)*Power(rxj, 2LL)*

                           (637875LL + 1148175LL*rxj + 1020600LL*Power(rxj, 2LL) +

           595350LL*Power(rxj, 3LL) + 255150LL*Power(rxj, 4LL) + 85050LL*Power(rxj, 5LL) +

           22680LL*Power(rxj, 6LL) + 4860LL*Power(rxj, 7LL) + 810LL*Power(rxj, 8LL) +

           34LL*Power(rxj, 9LL)) + 3LL*Power(rxi, 8LL)*Power(rxj, 14LL)*

                           (-593408025LL + 946053675LL*rxj - 394427880LL*Power(rxj, 2LL) -

           315870660LL*Power(rxj, 3LL) - 53891460LL*Power(rxj, 4LL) +

           910980LL*Power(rxj, 5LL) + 1409520LL*Power(rxj, 6LL) + 192168LL*Power(rxj, 7LL) +

           11196LL*Power(rxj, 8LL) + 236LL*Power(rxj, 9LL))))/

                         (42525LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 15LL)*Power(rxi + rxj, 15LL))

                         );
        }

    }
    return S;
}


cl_R Slater_5S_3S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_3S_5S(r, xj, xi);
}
