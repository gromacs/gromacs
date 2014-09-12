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

cl_R DSlater_1S_5S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-2875101075LL*xi + 3193344000LL*exp(2LL*r*xi)*xi -

                  5113716300LL*r*Power(xi, 2LL) - 4478789700LL*Power(r, 2LL)*Power(xi, 3LL) -

                  2564654400LL*Power(r, 3LL)*Power(xi, 4LL) -

                  1073595600LL*Power(r, 4LL)*Power(xi, 5LL) -

                  347276160LL*Power(r, 5LL)*Power(xi, 6LL) - 89147520LL*Power(r, 6LL)*Power(xi, 7LL) -

                  18247680LL*Power(r, 7LL)*Power(xi, 8LL) - 2914560LL*Power(r, 8LL)*Power(xi, 9LL) -

                  337920LL*Power(r, 9LL)*Power(xi, 10LL) - 22528LL*Power(r, 10LL)*Power(xi, 11LL))/

                (1.596672e9*exp(2LL*r*xi)*r) +

                (-1596672000LL + 1596672000LL*exp(2LL*r*xi) - 2875101075LL*r*xi -

                 2556858150LL*Power(r, 2LL)*Power(xi, 2LL) -

                 1492929900LL*Power(r, 3LL)*Power(xi, 3LL) - 641163600LL*Power(r, 4LL)*Power(xi, 4LL) -

                 214719120LL*Power(r, 5LL)*Power(xi, 5LL) - 57879360LL*Power(r, 6LL)*Power(xi, 6LL) -

                 12735360LL*Power(r, 7LL)*Power(xi, 7LL) - 2280960LL*Power(r, 8LL)*Power(xi, 8LL) -

                 323840LL*Power(r, 9LL)*Power(xi, 9LL) - 33792LL*Power(r, 10LL)*Power(xi, 10LL) -

                 2048LL*Power(r, 11LL)*Power(xi, 11LL))/(1.596672e9*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-1596672000LL + 1596672000LL*exp(2LL*r*xi) - 2875101075LL*r*xi -

                     2556858150LL*Power(r, 2LL)*Power(xi, 2LL) -

                     1492929900LL*Power(r, 3LL)*Power(xi, 3LL) -

                     641163600LL*Power(r, 4LL)*Power(xi, 4LL) - 214719120LL*Power(r, 5LL)*Power(xi, 5LL) -

                     57879360LL*Power(r, 6LL)*Power(xi, 6LL) - 12735360LL*Power(r, 7LL)*Power(xi, 7LL) -

                     2280960LL*Power(r, 8LL)*Power(xi, 8LL) - 323840LL*Power(r, 9LL)*Power(xi, 9LL) -

                     33792LL*Power(r, 10LL)*Power(xi, 10LL) - 2048LL*Power(r, 11LL)*Power(xi, 11LL)))/

                (7.98336e8*exp(2LL*r*xi)*r)

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
            S = (14175LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                 2835LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-35LL*Power(xi, 10LL) - 5LL*r*Power(xi, 11LL) - 495LL*Power(xi, 8LL)*Power(xj, 2LL) -

                  55LL*r*Power(xi, 9LL)*Power(xj, 2LL) - 1254LL*Power(xi, 6LL)*Power(xj, 4LL) -

                  66LL*r*Power(xi, 7LL)*Power(xj, 4LL) - 726LL*Power(xi, 4LL)*Power(xj, 6LL) +

                  66LL*r*Power(xi, 5LL)*Power(xj, 6LL) - 55LL*Power(xi, 2LL)*Power(xj, 8LL) +

                  55LL*r*Power(xi, 3LL)*Power(xj, 8LL) + 5LL*Power(xj, 10LL) + 5LL*r*xi*Power(xj, 10LL)

                 ) + exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-(Power(xi, 18LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                     13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                                     1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                                     108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                                     2LL*Power(r, 9LL)*Power(xj, 9LL))) +

                  9LL*Power(xi, 16LL)*Power(xj, 2LL)*

                  (17325LL + 31185LL*r*xj + 27720LL*Power(r, 2LL)*Power(xj, 2LL) +

            16170LL*Power(r, 3LL)*Power(xj, 3LL) + 6930LL*Power(r, 4LL)*Power(xj, 4LL) +

            2310LL*Power(r, 5LL)*Power(xj, 5LL) + 616LL*Power(r, 6LL)*Power(xj, 6LL) +

            132LL*Power(r, 7LL)*Power(xj, 7LL) + 22LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  126LL*Power(xi, 10LL)*Power(xj, 8LL)*

                  (37125LL + 66825LL*r*xj + 59400LL*Power(r, 2LL)*Power(xj, 2LL) +

            34725LL*Power(r, 3LL)*Power(xj, 3LL) + 14625LL*Power(r, 4LL)*Power(xj, 4LL) +

            5043LL*Power(r, 5LL)*Power(xj, 5LL) + 1396LL*Power(r, 6LL)*Power(xj, 6LL) +

            276LL*Power(r, 7LL)*Power(xj, 7LL) + 34LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  126LL*Power(xi, 8LL)*Power(xj, 10LL)*

                  (51975LL + 93420LL*r*xj + 84240LL*Power(r, 2LL)*Power(xj, 2LL) +

            46815LL*Power(r, 3LL)*Power(xj, 3LL) + 20835LL*Power(r, 4LL)*Power(xj, 4LL) +

            7485LL*Power(r, 5LL)*Power(xj, 5LL) + 1964LL*Power(r, 6LL)*Power(xj, 6LL) +

            348LL*Power(r, 7LL)*Power(xj, 7LL) + 38LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  9LL*Power(xi, 2LL)*Power(xj, 16LL)*

                  (-135135LL + 405405LL*r*xj + 582120LL*Power(r, 2LL)*Power(xj, 2LL) +

            346500LL*Power(r, 3LL)*Power(xj, 3LL) + 124740LL*Power(r, 4LL)*Power(xj, 4LL) +

            30492LL*Power(r, 5LL)*Power(xj, 5LL) + 5264LL*Power(r, 6LL)*Power(xj, 6LL) +

            636LL*Power(r, 7LL)*Power(xj, 7LL) + 50LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  Power(xj, 18LL)*(2837835LL + 3648645LL*r*xj +

                                   2245320LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   873180LL*Power(r, 3LL)*Power(xj, 3LL) + 238140LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   47628LL*Power(r, 5LL)*Power(xj, 5LL) + 7056LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   756LL*Power(r, 7LL)*Power(xj, 7LL) + 54LL*Power(r, 8LL)*Power(xj, 8LL) +

                                   2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  9LL*Power(xi, 14LL)*Power(xj, 4LL)*

                  (86625LL + 155925LL*r*xj + 138600LL*Power(r, 2LL)*Power(xj, 2LL) +

            80850LL*Power(r, 3LL)*Power(xj, 3LL) + 34650LL*Power(r, 4LL)*Power(xj, 4LL) +

            11550LL*Power(r, 5LL)*Power(xj, 5LL) + 3080LL*Power(r, 6LL)*Power(xj, 6LL) +

            672LL*Power(r, 7LL)*Power(xj, 7LL) + 104LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  21LL*Power(xi, 12LL)*Power(xj, 6LL)*

                  (111375LL + 200475LL*r*xj + 178200LL*Power(r, 2LL)*Power(xj, 2LL) +

            103950LL*Power(r, 3LL)*Power(xj, 3LL) + 44550LL*Power(r, 4LL)*Power(xj, 4LL) +

            14778LL*Power(r, 5LL)*Power(xj, 5LL) + 4056LL*Power(r, 6LL)*Power(xj, 6LL) +

            864LL*Power(r, 7LL)*Power(xj, 7LL) + 120LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  21LL*Power(xi, 6LL)*Power(xj, 12LL)*

                  (307125LL + 594945LL*r*xj + 456840LL*Power(r, 2LL)*Power(xj, 2LL) +

            281790LL*Power(r, 3LL)*Power(xj, 3LL) + 137430LL*Power(r, 4LL)*Power(xj, 4LL) +

            47250LL*Power(r, 5LL)*Power(xj, 5LL) + 11064LL*Power(r, 6LL)*Power(xj, 6LL) +

            1728LL*Power(r, 7LL)*Power(xj, 7LL) + 168LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  9LL*Power(xi, 4LL)*Power(xj, 14LL)*

                  (675675LL + 675675LL*r*xj + 748440LL*Power(r, 2LL)*Power(xj, 2LL) +

            561330LL*Power(r, 3LL)*Power(xj, 3LL) + 256410LL*Power(r, 4LL)*Power(xj, 4LL) +

            76230LL*Power(r, 5LL)*Power(xj, 5LL) + 15400LL*Power(r, 6LL)*Power(xj, 6LL) +

            2112LL*Power(r, 7LL)*Power(xj, 7LL) + 184LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (14175LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 11LL)*

                 Power(xi + xj, 11LL)) + (2LL*(14175LL*exp(2LL*r*(xi + xj))*

                                               Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                                               2835LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                                               (-35LL*Power(xi, 10LL) - 5LL*r*Power(xi, 11LL) -

                                       495LL*Power(xi, 8LL)*Power(xj, 2LL) - 55LL*r*Power(xi, 9LL)*Power(xj, 2LL) -

                                       1254LL*Power(xi, 6LL)*Power(xj, 4LL) - 66LL*r*Power(xi, 7LL)*Power(xj, 4LL) -

                                       726LL*Power(xi, 4LL)*Power(xj, 6LL) + 66LL*r*Power(xi, 5LL)*Power(xj, 6LL) -

                                       55LL*Power(xi, 2LL)*Power(xj, 8LL) + 55LL*r*Power(xi, 3LL)*Power(xj, 8LL) +

                                       5LL*Power(xj, 10LL) + 5LL*r*xi*Power(xj, 10LL)) +

                                               exp(2LL*r*xi)*Power(xi, 4LL)*

                                               (-(Power(xi, 18LL)*(14175LL + 25515LL*r*xj +

                                                                   22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                   13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                   5670LL*Power(r, 4LL)*Power(xj, 4LL) + 1890LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                   504LL*Power(r, 6LL)*Power(xj, 6LL) + 108LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                   18LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL))) +

                                       9LL*Power(xi, 16LL)*Power(xj, 2LL)*

                                       (17325LL + 31185LL*r*xj + 27720LL*Power(r, 2LL)*Power(xj, 2LL) +

              16170LL*Power(r, 3LL)*Power(xj, 3LL) + 6930LL*Power(r, 4LL)*Power(xj, 4LL) +

              2310LL*Power(r, 5LL)*Power(xj, 5LL) + 616LL*Power(r, 6LL)*Power(xj, 6LL) +

              132LL*Power(r, 7LL)*Power(xj, 7LL) + 22LL*Power(r, 8LL)*Power(xj, 8LL) +

              2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       126LL*Power(xi, 10LL)*Power(xj, 8LL)*

                                       (37125LL + 66825LL*r*xj + 59400LL*Power(r, 2LL)*Power(xj, 2LL) +

              34725LL*Power(r, 3LL)*Power(xj, 3LL) + 14625LL*Power(r, 4LL)*Power(xj, 4LL) +

              5043LL*Power(r, 5LL)*Power(xj, 5LL) + 1396LL*Power(r, 6LL)*Power(xj, 6LL) +

              276LL*Power(r, 7LL)*Power(xj, 7LL) + 34LL*Power(r, 8LL)*Power(xj, 8LL) +

              2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       126LL*Power(xi, 8LL)*Power(xj, 10LL)*

                                       (51975LL + 93420LL*r*xj + 84240LL*Power(r, 2LL)*Power(xj, 2LL) +

              46815LL*Power(r, 3LL)*Power(xj, 3LL) + 20835LL*Power(r, 4LL)*Power(xj, 4LL) +

              7485LL*Power(r, 5LL)*Power(xj, 5LL) + 1964LL*Power(r, 6LL)*Power(xj, 6LL) +

              348LL*Power(r, 7LL)*Power(xj, 7LL) + 38LL*Power(r, 8LL)*Power(xj, 8LL) +

              2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       9LL*Power(xi, 2LL)*Power(xj, 16LL)*

                                       (-135135LL + 405405LL*r*xj + 582120LL*Power(r, 2LL)*Power(xj, 2LL) +

              346500LL*Power(r, 3LL)*Power(xj, 3LL) +

              124740LL*Power(r, 4LL)*Power(xj, 4LL) +

              30492LL*Power(r, 5LL)*Power(xj, 5LL) + 5264LL*Power(r, 6LL)*Power(xj, 6LL) +

              636LL*Power(r, 7LL)*Power(xj, 7LL) + 50LL*Power(r, 8LL)*Power(xj, 8LL) +

              2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       Power(xj, 18LL)*(2837835LL + 3648645LL*r*xj +

                                                        2245320LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        873180LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                        238140LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                        47628LL*Power(r, 5LL)*Power(xj, 5LL) + 7056LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                        756LL*Power(r, 7LL)*Power(xj, 7LL) + 54LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                        2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       9LL*Power(xi, 14LL)*Power(xj, 4LL)*

                                       (86625LL + 155925LL*r*xj + 138600LL*Power(r, 2LL)*Power(xj, 2LL) +

              80850LL*Power(r, 3LL)*Power(xj, 3LL) + 34650LL*Power(r, 4LL)*Power(xj, 4LL) +

              11550LL*Power(r, 5LL)*Power(xj, 5LL) + 3080LL*Power(r, 6LL)*Power(xj, 6LL) +

              672LL*Power(r, 7LL)*Power(xj, 7LL) + 104LL*Power(r, 8LL)*Power(xj, 8LL) +

              8LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       21LL*Power(xi, 12LL)*Power(xj, 6LL)*

                                       (111375LL + 200475LL*r*xj + 178200LL*Power(r, 2LL)*Power(xj, 2LL) +

              103950LL*Power(r, 3LL)*Power(xj, 3LL) +

              44550LL*Power(r, 4LL)*Power(xj, 4LL) + 14778LL*Power(r, 5LL)*Power(xj, 5LL) +

              4056LL*Power(r, 6LL)*Power(xj, 6LL) + 864LL*Power(r, 7LL)*Power(xj, 7LL) +

              120LL*Power(r, 8LL)*Power(xj, 8LL) + 8LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                       21LL*Power(xi, 6LL)*Power(xj, 12LL)*

                                       (307125LL + 594945LL*r*xj + 456840LL*Power(r, 2LL)*Power(xj, 2LL) +

              281790LL*Power(r, 3LL)*Power(xj, 3LL) +

              137430LL*Power(r, 4LL)*Power(xj, 4LL) +

              47250LL*Power(r, 5LL)*Power(xj, 5LL) + 11064LL*Power(r, 6LL)*Power(xj, 6LL) +

              1728LL*Power(r, 7LL)*Power(xj, 7LL) + 168LL*Power(r, 8LL)*Power(xj, 8LL) +

              8LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                       9LL*Power(xi, 4LL)*Power(xj, 14LL)*

                                       (675675LL + 675675LL*r*xj + 748440LL*Power(r, 2LL)*Power(xj, 2LL) +

              561330LL*Power(r, 3LL)*Power(xj, 3LL) +

              256410LL*Power(r, 4LL)*Power(xj, 4LL) + 76230LL*Power(r, 5LL)*Power(xj, 5LL) +

              15400LL*Power(r, 6LL)*Power(xj, 6LL) + 2112LL*Power(r, 7LL)*Power(xj, 7LL) +

              184LL*Power(r, 8LL)*Power(xj, 8LL) + 8LL*Power(r, 9LL)*Power(xj, 9LL)))))/

                (14175LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 11LL)*Power(xi + xj, 10LL)) -

                (28350LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                 2835LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-5LL*Power(xi, 11LL) - 55LL*Power(xi, 9LL)*Power(xj, 2LL) -

        66LL*Power(xi, 7LL)*Power(xj, 4LL) + 66LL*Power(xi, 5LL)*Power(xj, 6LL) +

        55LL*Power(xi, 3LL)*Power(xj, 8LL) + 5LL*xi*Power(xj, 10LL)) +

                 5670LL*exp(2LL*r*xj)*Power(xj, 13LL)*

                 (-35LL*Power(xi, 10LL) - 5LL*r*Power(xi, 11LL) - 495LL*Power(xi, 8LL)*Power(xj, 2LL) -

        55LL*r*Power(xi, 9LL)*Power(xj, 2LL) - 1254LL*Power(xi, 6LL)*Power(xj, 4LL) -

        66LL*r*Power(xi, 7LL)*Power(xj, 4LL) - 726LL*Power(xi, 4LL)*Power(xj, 6LL) +

        66LL*r*Power(xi, 5LL)*Power(xj, 6LL) - 55LL*Power(xi, 2LL)*Power(xj, 8LL) +

        55LL*r*Power(xi, 3LL)*Power(xj, 8LL) + 5LL*Power(xj, 10LL) + 5LL*r*xi*Power(xj, 10LL))

                 + exp(2LL*r*xi)*Power(xi, 4LL)*(-(Power(xi, 18LL)*

                                                   (25515LL*xj + 45360LL*r*Power(xj, 2LL) +

                                          39690LL*Power(r, 2LL)*Power(xj, 3LL) + 22680LL*Power(r, 3LL)*Power(xj, 4LL) +

                                          9450LL*Power(r, 4LL)*Power(xj, 5LL) + 3024LL*Power(r, 5LL)*Power(xj, 6LL) +

                                          756LL*Power(r, 6LL)*Power(xj, 7LL) + 144LL*Power(r, 7LL)*Power(xj, 8LL) +

                                          18LL*Power(r, 8LL)*Power(xj, 9LL))) +

                                                 9LL*Power(xi, 16LL)*Power(xj, 2LL)*

                                                 (31185LL*xj + 55440LL*r*Power(xj, 2LL) + 48510LL*Power(r, 2LL)*Power(xj, 3LL) +

                                        27720LL*Power(r, 3LL)*Power(xj, 4LL) + 11550LL*Power(r, 4LL)*Power(xj, 5LL) +

                                        3696LL*Power(r, 5LL)*Power(xj, 6LL) + 924LL*Power(r, 6LL)*Power(xj, 7LL) +

                                        176LL*Power(r, 7LL)*Power(xj, 8LL) + 18LL*Power(r, 8LL)*Power(xj, 9LL)) -

                                                 126LL*Power(xi, 10LL)*Power(xj, 8LL)*

                                                 (66825LL*xj + 118800LL*r*Power(xj, 2LL) +

                                        104175LL*Power(r, 2LL)*Power(xj, 3LL) + 58500LL*Power(r, 3LL)*Power(xj, 4LL) +

                                        25215LL*Power(r, 4LL)*Power(xj, 5LL) + 8376LL*Power(r, 5LL)*Power(xj, 6LL) +

                                        1932LL*Power(r, 6LL)*Power(xj, 7LL) + 272LL*Power(r, 7LL)*Power(xj, 8LL) +

                                        18LL*Power(r, 8LL)*Power(xj, 9LL)) +

                                                 126LL*Power(xi, 8LL)*Power(xj, 10LL)*

                                                 (93420LL*xj + 168480LL*r*Power(xj, 2LL) +

                                        140445LL*Power(r, 2LL)*Power(xj, 3LL) + 83340LL*Power(r, 3LL)*Power(xj, 4LL) +

                                        37425LL*Power(r, 4LL)*Power(xj, 5LL) + 11784LL*Power(r, 5LL)*Power(xj, 6LL) +

                                        2436LL*Power(r, 6LL)*Power(xj, 7LL) + 304LL*Power(r, 7LL)*Power(xj, 8LL) +

                                        18LL*Power(r, 8LL)*Power(xj, 9LL)) -

                                                 9LL*Power(xi, 2LL)*Power(xj, 16LL)*

                                                 (405405LL*xj + 1164240LL*r*Power(xj, 2LL) +

                                        1039500LL*Power(r, 2LL)*Power(xj, 3LL) +

                                        498960LL*Power(r, 3LL)*Power(xj, 4LL) + 152460LL*Power(r, 4LL)*Power(xj, 5LL) +

                                        31584LL*Power(r, 5LL)*Power(xj, 6LL) + 4452LL*Power(r, 6LL)*Power(xj, 7LL) +

                                        400LL*Power(r, 7LL)*Power(xj, 8LL) + 18LL*Power(r, 8LL)*Power(xj, 9LL)) +

                                                 Power(xj, 18LL)*(3648645LL*xj + 4490640LL*r*Power(xj, 2LL) +

                                                                  2619540LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                                  952560LL*Power(r, 3LL)*Power(xj, 4LL) + 238140LL*Power(r, 4LL)*Power(xj, 5LL) +

                                                                  42336LL*Power(r, 5LL)*Power(xj, 6LL) + 5292LL*Power(r, 6LL)*Power(xj, 7LL) +

                                                                  432LL*Power(r, 7LL)*Power(xj, 8LL) + 18LL*Power(r, 8LL)*Power(xj, 9LL)) -

                                                 9LL*Power(xi, 14LL)*Power(xj, 4LL)*

                                                 (155925LL*xj + 277200LL*r*Power(xj, 2LL) +

                                        242550LL*Power(r, 2LL)*Power(xj, 3LL) + 138600LL*Power(r, 3LL)*Power(xj, 4LL) +

                                        57750LL*Power(r, 4LL)*Power(xj, 5LL) + 18480LL*Power(r, 5LL)*Power(xj, 6LL) +

                                        4704LL*Power(r, 6LL)*Power(xj, 7LL) + 832LL*Power(r, 7LL)*Power(xj, 8LL) +

                                        72LL*Power(r, 8LL)*Power(xj, 9LL)) +

                                                 21LL*Power(xi, 12LL)*Power(xj, 6LL)*

                                                 (200475LL*xj + 356400LL*r*Power(xj, 2LL) +

                                        311850LL*Power(r, 2LL)*Power(xj, 3LL) + 178200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                        73890LL*Power(r, 4LL)*Power(xj, 5LL) + 24336LL*Power(r, 5LL)*Power(xj, 6LL) +

                                        6048LL*Power(r, 6LL)*Power(xj, 7LL) + 960LL*Power(r, 7LL)*Power(xj, 8LL) +

                                        72LL*Power(r, 8LL)*Power(xj, 9LL)) -

                                                 21LL*Power(xi, 6LL)*Power(xj, 12LL)*

                                                 (594945LL*xj + 913680LL*r*Power(xj, 2LL) +

                                        845370LL*Power(r, 2LL)*Power(xj, 3LL) + 549720LL*Power(r, 3LL)*Power(xj, 4LL) +

                                        236250LL*Power(r, 4LL)*Power(xj, 5LL) + 66384LL*Power(r, 5LL)*Power(xj, 6LL) +

                                        12096LL*Power(r, 6LL)*Power(xj, 7LL) + 1344LL*Power(r, 7LL)*Power(xj, 8LL) +

                                        72LL*Power(r, 8LL)*Power(xj, 9LL)) +

                                                 9LL*Power(xi, 4LL)*Power(xj, 14LL)*

                                                 (675675LL*xj + 1496880LL*r*Power(xj, 2LL) +

                                        1683990LL*Power(r, 2LL)*Power(xj, 3LL) +

                                        1025640LL*Power(r, 3LL)*Power(xj, 4LL) + 381150LL*Power(r, 4LL)*Power(xj, 5LL) +

                                        92400LL*Power(r, 5LL)*Power(xj, 6LL) + 14784LL*Power(r, 6LL)*Power(xj, 7LL) +

                                        1472LL*Power(r, 7LL)*Power(xj, 8LL) + 72LL*Power(r, 8LL)*Power(xj, 9LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 5LL)*

                 (-(Power(xi, 18LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                     13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                                     1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                                     108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                                     2LL*Power(r, 9LL)*Power(xj, 9LL))) +

        9LL*Power(xi, 16LL)*Power(xj, 2LL)*

        (17325LL + 31185LL*r*xj + 27720LL*Power(r, 2LL)*Power(xj, 2LL) +

            16170LL*Power(r, 3LL)*Power(xj, 3LL) + 6930LL*Power(r, 4LL)*Power(xj, 4LL) +

            2310LL*Power(r, 5LL)*Power(xj, 5LL) + 616LL*Power(r, 6LL)*Power(xj, 6LL) +

            132LL*Power(r, 7LL)*Power(xj, 7LL) + 22LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        126LL*Power(xi, 10LL)*Power(xj, 8LL)*

        (37125LL + 66825LL*r*xj + 59400LL*Power(r, 2LL)*Power(xj, 2LL) +

            34725LL*Power(r, 3LL)*Power(xj, 3LL) + 14625LL*Power(r, 4LL)*Power(xj, 4LL) +

            5043LL*Power(r, 5LL)*Power(xj, 5LL) + 1396LL*Power(r, 6LL)*Power(xj, 6LL) +

            276LL*Power(r, 7LL)*Power(xj, 7LL) + 34LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

        126LL*Power(xi, 8LL)*Power(xj, 10LL)*

        (51975LL + 93420LL*r*xj + 84240LL*Power(r, 2LL)*Power(xj, 2LL) +

            46815LL*Power(r, 3LL)*Power(xj, 3LL) + 20835LL*Power(r, 4LL)*Power(xj, 4LL) +

            7485LL*Power(r, 5LL)*Power(xj, 5LL) + 1964LL*Power(r, 6LL)*Power(xj, 6LL) +

            348LL*Power(r, 7LL)*Power(xj, 7LL) + 38LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        9LL*Power(xi, 2LL)*Power(xj, 16LL)*

        (-135135LL + 405405LL*r*xj + 582120LL*Power(r, 2LL)*Power(xj, 2LL) +

            346500LL*Power(r, 3LL)*Power(xj, 3LL) + 124740LL*Power(r, 4LL)*Power(xj, 4LL) +

            30492LL*Power(r, 5LL)*Power(xj, 5LL) + 5264LL*Power(r, 6LL)*Power(xj, 6LL) +

            636LL*Power(r, 7LL)*Power(xj, 7LL) + 50LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

        Power(xj, 18LL)*(2837835LL + 3648645LL*r*xj +

                         2245320LL*Power(r, 2LL)*Power(xj, 2LL) + 873180LL*Power(r, 3LL)*Power(xj, 3LL) +

                         238140LL*Power(r, 4LL)*Power(xj, 4LL) + 47628LL*Power(r, 5LL)*Power(xj, 5LL) +

                         7056LL*Power(r, 6LL)*Power(xj, 6LL) + 756LL*Power(r, 7LL)*Power(xj, 7LL) +

                         54LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        9LL*Power(xi, 14LL)*Power(xj, 4LL)*

        (86625LL + 155925LL*r*xj + 138600LL*Power(r, 2LL)*Power(xj, 2LL) +

            80850LL*Power(r, 3LL)*Power(xj, 3LL) + 34650LL*Power(r, 4LL)*Power(xj, 4LL) +

            11550LL*Power(r, 5LL)*Power(xj, 5LL) + 3080LL*Power(r, 6LL)*Power(xj, 6LL) +

            672LL*Power(r, 7LL)*Power(xj, 7LL) + 104LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL)) +

        21LL*Power(xi, 12LL)*Power(xj, 6LL)*

        (111375LL + 200475LL*r*xj + 178200LL*Power(r, 2LL)*Power(xj, 2LL) +

            103950LL*Power(r, 3LL)*Power(xj, 3LL) + 44550LL*Power(r, 4LL)*Power(xj, 4LL) +

            14778LL*Power(r, 5LL)*Power(xj, 5LL) + 4056LL*Power(r, 6LL)*Power(xj, 6LL) +

            864LL*Power(r, 7LL)*Power(xj, 7LL) + 120LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL)) -

        21LL*Power(xi, 6LL)*Power(xj, 12LL)*

        (307125LL + 594945LL*r*xj + 456840LL*Power(r, 2LL)*Power(xj, 2LL) +

            281790LL*Power(r, 3LL)*Power(xj, 3LL) + 137430LL*Power(r, 4LL)*Power(xj, 4LL) +

            47250LL*Power(r, 5LL)*Power(xj, 5LL) + 11064LL*Power(r, 6LL)*Power(xj, 6LL) +

            1728LL*Power(r, 7LL)*Power(xj, 7LL) + 168LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL)) +

        9LL*Power(xi, 4LL)*Power(xj, 14LL)*

        (675675LL + 675675LL*r*xj + 748440LL*Power(r, 2LL)*Power(xj, 2LL) +

            561330LL*Power(r, 3LL)*Power(xj, 3LL) + 256410LL*Power(r, 4LL)*Power(xj, 4LL) +

            76230LL*Power(r, 5LL)*Power(xj, 5LL) + 15400LL*Power(r, 6LL)*Power(xj, 6LL) +

            2112LL*Power(r, 7LL)*Power(xj, 7LL) + 184LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (14175LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 11LL)*Power(xi + xj, 11LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_5S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_1S_5S(r, xj, xi);
}
