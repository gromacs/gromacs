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

cl_R DSlater_2S_3S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-7430535LL*xi + 8709120LL*exp(2LL*r*xi)*xi - 12303900LL*r*Power(xi, 2LL) -

                  9826110LL*Power(r, 2LL)*Power(xi, 3LL) - 5004720LL*Power(r, 3LL)*Power(xi, 4LL) -

                  1806840LL*Power(r, 4LL)*Power(xi, 5LL) - 483840LL*Power(r, 5LL)*Power(xi, 6LL) -

                  96768LL*Power(r, 6LL)*Power(xi, 7LL) - 13824LL*Power(r, 7LL)*Power(xi, 8LL) -

                  1152LL*Power(r, 8LL)*Power(xi, 9LL))/(4.35456e6*exp(2LL*r*xi)*r) +

                (-4354560LL + 4354560LL*exp(2LL*r*xi) - 7430535LL*r*xi -

                 6151950LL*Power(r, 2LL)*Power(xi, 2LL) - 3275370LL*Power(r, 3LL)*Power(xi, 3LL) -

                 1251180LL*Power(r, 4LL)*Power(xi, 4LL) - 361368LL*Power(r, 5LL)*Power(xi, 5LL) -

                 80640LL*Power(r, 6LL)*Power(xi, 6LL) - 13824LL*Power(r, 7LL)*Power(xi, 7LL) -

                 1728LL*Power(r, 8LL)*Power(xi, 8LL) - 128LL*Power(r, 9LL)*Power(xi, 9LL))/

                (4.35456e6*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-4354560LL + 4354560LL*exp(2LL*r*xi) - 7430535LL*r*xi -

                     6151950LL*Power(r, 2LL)*Power(xi, 2LL) - 3275370LL*Power(r, 3LL)*Power(xi, 3LL) -

                     1251180LL*Power(r, 4LL)*Power(xi, 4LL) - 361368LL*Power(r, 5LL)*Power(xi, 5LL) -

                     80640LL*Power(r, 6LL)*Power(xi, 6LL) - 13824LL*Power(r, 7LL)*Power(xi, 7LL) -

                     1728LL*Power(r, 8LL)*Power(xi, 8LL) - 128LL*Power(r, 9LL)*Power(xi, 9LL)))/

                (2.17728e6*exp(2LL*r*xi)*r)

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
            S = (90LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 9LL) +

                 5LL*exp(2LL*r*xj)*Power(xj, 8LL)*

                 (-90LL*Power(r, 2LL)*Power(xi, 12LL) - 6LL*Power(r, 3LL)*Power(xi, 13LL) +

                  18LL*Power(xj, 10LL) + 27LL*r*xi*Power(xj, 10LL) +

                  18LL*Power(xi, 2LL)*Power(xj, 8LL)*(-9LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  162LL*Power(xi, 4LL)*Power(xj, 6LL)*(-4LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  198LL*Power(xi, 10LL)*(5LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  108LL*Power(xi, 6LL)*Power(xj, 4LL)*(36LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  2LL*r*Power(xi, 5LL)*Power(xj, 6LL)*(675LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  18LL*r*Power(xi, 7LL)*Power(xj, 4LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  3LL*r*Power(xi, 3LL)*Power(xj, 8LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  r*Power(xi, 11LL)*(495LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  9LL*r*Power(xi, 9LL)*Power(xj, 2LL)*(-233LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  6LL*Power(xi, 8LL)*Power(xj, 2LL)*(-1063LL + 90LL*Power(r, 2LL)*Power(xj, 2LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-90LL*Power(xi, 6LL)*Power(xj, 6LL)*

                  (42LL + 65LL*r*xj + 76LL*Power(r, 2LL)*Power(xj, 2LL) +

            22LL*Power(r, 3LL)*Power(xj, 3LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  2LL*Power(xj, 12LL)*(2970LL + 2475LL*r*xj + 900LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       180LL*Power(r, 3LL)*Power(xj, 3LL) + 20LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       Power(r, 5LL)*Power(xj, 5LL)) +

                  10LL*Power(xi, 8LL)*Power(xj, 4LL)*

                  (162LL + 270LL*r*xj + 216LL*Power(r, 2LL)*Power(xj, 2LL) +

            122LL*Power(r, 3LL)*Power(xj, 3LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL) +

            Power(r, 5LL)*Power(xj, 5LL)) -

                  5LL*Power(xi, 4LL)*Power(xj, 8LL)*

                  (-639LL - 3555LL*r*xj - 1452LL*Power(r, 2LL)*Power(xj, 2LL) -

            174LL*Power(r, 3LL)*Power(xj, 3LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL) +

            2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  Power(xi, 12LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  Power(xi, 10LL)*Power(xj, 2LL)*

                  (405LL + 675LL*r*xj + 540LL*Power(r, 2LL)*Power(xj, 2LL) +

            270LL*Power(r, 3LL)*Power(xj, 3LL) + 90LL*Power(r, 4LL)*Power(xj, 4LL) +

            8LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  Power(xi, 2LL)*Power(xj, 10LL)*

                  (-21615LL - 9075LL*r*xj - 300LL*Power(r, 2LL)*Power(xj, 2LL) +

            490LL*Power(r, 3LL)*Power(xj, 3LL) + 110LL*Power(r, 4LL)*Power(xj, 4LL) +

            8LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (90LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 9LL)*Power(xi + xj, 9LL))

                + (90LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 9LL) +

                   5LL*exp(2LL*r*xj)*Power(xj, 8LL)*

                   (-90LL*Power(r, 2LL)*Power(xi, 12LL) - 6LL*Power(r, 3LL)*Power(xi, 13LL) +

                    18LL*Power(xj, 10LL) + 27LL*r*xi*Power(xj, 10LL) +

                    18LL*Power(xi, 2LL)*Power(xj, 8LL)*(-9LL + Power(r, 2LL)*Power(xj, 2LL)) -

                    162LL*Power(xi, 4LL)*Power(xj, 6LL)*(-4LL + Power(r, 2LL)*Power(xj, 2LL)) -

                    198LL*Power(xi, 10LL)*(5LL + Power(r, 2LL)*Power(xj, 2LL)) -

                    108LL*Power(xi, 6LL)*Power(xj, 4LL)*(36LL + Power(r, 2LL)*Power(xj, 2LL)) +

                    2LL*r*Power(xi, 5LL)*Power(xj, 6LL)*(675LL + Power(r, 2LL)*Power(xj, 2LL)) -

                    18LL*r*Power(xi, 7LL)*Power(xj, 4LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    3LL*r*Power(xi, 3LL)*Power(xj, 8LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    r*Power(xi, 11LL)*(495LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    9LL*r*Power(xi, 9LL)*Power(xj, 2LL)*(-233LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    6LL*Power(xi, 8LL)*Power(xj, 2LL)*(-1063LL + 90LL*Power(r, 2LL)*Power(xj, 2LL))) -

                   2LL*exp(2LL*r*xi)*Power(xi, 6LL)*

                   (-90LL*Power(xi, 6LL)*Power(xj, 6LL)*

                    (42LL + 65LL*r*xj + 76LL*Power(r, 2LL)*Power(xj, 2LL) +

            22LL*Power(r, 3LL)*Power(xj, 3LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    2LL*Power(xj, 12LL)*(2970LL + 2475LL*r*xj + 900LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         180LL*Power(r, 3LL)*Power(xj, 3LL) + 20LL*Power(r, 4LL)*Power(xj, 4LL) +

                                         Power(r, 5LL)*Power(xj, 5LL)) +

                    10LL*Power(xi, 8LL)*Power(xj, 4LL)*

                    (162LL + 270LL*r*xj + 216LL*Power(r, 2LL)*Power(xj, 2LL) +

            122LL*Power(r, 3LL)*Power(xj, 3LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL) +

            Power(r, 5LL)*Power(xj, 5LL)) -

                    5LL*Power(xi, 4LL)*Power(xj, 8LL)*

                    (-639LL - 3555LL*r*xj - 1452LL*Power(r, 2LL)*Power(xj, 2LL) -

            174LL*Power(r, 3LL)*Power(xj, 3LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL) +

            2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                    Power(xi, 12LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                     30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                     2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                    Power(xi, 10LL)*Power(xj, 2LL)*

                    (405LL + 675LL*r*xj + 540LL*Power(r, 2LL)*Power(xj, 2LL) +

            270LL*Power(r, 3LL)*Power(xj, 3LL) + 90LL*Power(r, 4LL)*Power(xj, 4LL) +

            8LL*Power(r, 5LL)*Power(xj, 5LL)) +

                    Power(xi, 2LL)*Power(xj, 10LL)*

                    (-21615LL - 9075LL*r*xj - 300LL*Power(r, 2LL)*Power(xj, 2LL) +

            490LL*Power(r, 3LL)*Power(xj, 3LL) + 110LL*Power(r, 4LL)*Power(xj, 4LL) +

            8LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (45LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 9LL)*Power(xi + xj, 8LL)) -

                (180LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi, 2LL) - Power(xj, 2LL), 9LL) +

                 5LL*exp(2LL*r*xj)*Power(xj, 8LL)*

                 (-180LL*r*Power(xi, 12LL) - 18LL*Power(r, 2LL)*Power(xi, 13LL) -

        396LL*r*Power(xi, 10LL)*Power(xj, 2LL) -

        4LL*Power(r, 2LL)*Power(xi, 11LL)*Power(xj, 2LL) +

        1080LL*r*Power(xi, 8LL)*Power(xj, 4LL) +

        72LL*Power(r, 2LL)*Power(xi, 9LL)*Power(xj, 4LL) -

        216LL*r*Power(xi, 6LL)*Power(xj, 6LL) -

        72LL*Power(r, 2LL)*Power(xi, 7LL)*Power(xj, 6LL) -

        324LL*r*Power(xi, 4LL)*Power(xj, 8LL) +

        4LL*Power(r, 2LL)*Power(xi, 5LL)*Power(xj, 8LL) + 27LL*xi*Power(xj, 10LL) +

        36LL*r*Power(xi, 2LL)*Power(xj, 10LL) +

        12LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 10LL) +

        2LL*Power(xi, 5LL)*Power(xj, 6LL)*(675LL + Power(r, 2LL)*Power(xj, 2LL)) -

        18LL*Power(xi, 7LL)*Power(xj, 4LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*Power(xi, 3LL)*Power(xj, 8LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        Power(xi, 11LL)*(495LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        9LL*Power(xi, 9LL)*Power(xj, 2LL)*(-233LL + 4LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 10LL*exp(2LL*r*xj)*Power(xj, 9LL)*

                 (-90LL*Power(r, 2LL)*Power(xi, 12LL) - 6LL*Power(r, 3LL)*Power(xi, 13LL) +

        18LL*Power(xj, 10LL) + 27LL*r*xi*Power(xj, 10LL) +

        18LL*Power(xi, 2LL)*Power(xj, 8LL)*(-9LL + Power(r, 2LL)*Power(xj, 2LL)) -

        162LL*Power(xi, 4LL)*Power(xj, 6LL)*(-4LL + Power(r, 2LL)*Power(xj, 2LL)) -

        198LL*Power(xi, 10LL)*(5LL + Power(r, 2LL)*Power(xj, 2LL)) -

        108LL*Power(xi, 6LL)*Power(xj, 4LL)*(36LL + Power(r, 2LL)*Power(xj, 2LL)) +

        2LL*r*Power(xi, 5LL)*Power(xj, 6LL)*(675LL + Power(r, 2LL)*Power(xj, 2LL)) -

        18LL*r*Power(xi, 7LL)*Power(xj, 4LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*r*Power(xi, 3LL)*Power(xj, 8LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        r*Power(xi, 11LL)*(495LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        9LL*r*Power(xi, 9LL)*Power(xj, 2LL)*(-233LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 8LL)*Power(xj, 2LL)*(-1063LL + 90LL*Power(r, 2LL)*Power(xj, 2LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-90LL*Power(xi, 6LL)*Power(xj, 6LL)*

        (65LL*xj + 152LL*r*Power(xj, 2LL) + 66LL*Power(r, 2LL)*Power(xj, 3LL) +

            8LL*Power(r, 3LL)*Power(xj, 4LL)) -

        2LL*Power(xj, 12LL)*(2475LL*xj + 1800LL*r*Power(xj, 2LL) +

                             540LL*Power(r, 2LL)*Power(xj, 3LL) + 80LL*Power(r, 3LL)*Power(xj, 4LL) +

                             5LL*Power(r, 4LL)*Power(xj, 5LL)) +

        10LL*Power(xi, 8LL)*Power(xj, 4LL)*

        (270LL*xj + 432LL*r*Power(xj, 2LL) + 366LL*Power(r, 2LL)*Power(xj, 3LL) +

            88LL*Power(r, 3LL)*Power(xj, 4LL) + 5LL*Power(r, 4LL)*Power(xj, 5LL)) -

        5LL*Power(xi, 4LL)*Power(xj, 8LL)*

        (-3555LL*xj - 2904LL*r*Power(xj, 2LL) - 522LL*Power(r, 2LL)*Power(xj, 3LL) +

            24LL*Power(r, 3LL)*Power(xj, 4LL) + 10LL*Power(r, 4LL)*Power(xj, 5LL)) +

        Power(xi, 12LL)*(75LL*xj + 120LL*r*Power(xj, 2LL) +

                         90LL*Power(r, 2LL)*Power(xj, 3LL) + 40LL*Power(r, 3LL)*Power(xj, 4LL) +

                         10LL*Power(r, 4LL)*Power(xj, 5LL)) -

        Power(xi, 10LL)*Power(xj, 2LL)*

        (675LL*xj + 1080LL*r*Power(xj, 2LL) + 810LL*Power(r, 2LL)*Power(xj, 3LL) +

            360LL*Power(r, 3LL)*Power(xj, 4LL) + 40LL*Power(r, 4LL)*Power(xj, 5LL)) +

        Power(xi, 2LL)*Power(xj, 10LL)*

        (-9075LL*xj - 600LL*r*Power(xj, 2LL) + 1470LL*Power(r, 2LL)*Power(xj, 3LL) +

            440LL*Power(r, 3LL)*Power(xj, 4LL) + 40LL*Power(r, 4LL)*Power(xj, 5LL))) -

                 4LL*exp(2LL*r*xi)*Power(xi, 7LL)*

                 (-90LL*Power(xi, 6LL)*Power(xj, 6LL)*

        (42LL + 65LL*r*xj + 76LL*Power(r, 2LL)*Power(xj, 2LL) +

            22LL*Power(r, 3LL)*Power(xj, 3LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

        2LL*Power(xj, 12LL)*(2970LL + 2475LL*r*xj + 900LL*Power(r, 2LL)*Power(xj, 2LL) +

                             180LL*Power(r, 3LL)*Power(xj, 3LL) + 20LL*Power(r, 4LL)*Power(xj, 4LL) +

                             Power(r, 5LL)*Power(xj, 5LL)) +

        10LL*Power(xi, 8LL)*Power(xj, 4LL)*

        (162LL + 270LL*r*xj + 216LL*Power(r, 2LL)*Power(xj, 2LL) +

            122LL*Power(r, 3LL)*Power(xj, 3LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL) +

            Power(r, 5LL)*Power(xj, 5LL)) -

        5LL*Power(xi, 4LL)*Power(xj, 8LL)*

        (-639LL - 3555LL*r*xj - 1452LL*Power(r, 2LL)*Power(xj, 2LL) -

            174LL*Power(r, 3LL)*Power(xj, 3LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL) +

            2LL*Power(r, 5LL)*Power(xj, 5LL)) +

        Power(xi, 12LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                         30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                         2LL*Power(r, 5LL)*Power(xj, 5LL)) -

        Power(xi, 10LL)*Power(xj, 2LL)*

        (405LL + 675LL*r*xj + 540LL*Power(r, 2LL)*Power(xj, 2LL) +

            270LL*Power(r, 3LL)*Power(xj, 3LL) + 90LL*Power(r, 4LL)*Power(xj, 4LL) +

            8LL*Power(r, 5LL)*Power(xj, 5LL)) +

        Power(xi, 2LL)*Power(xj, 10LL)*

        (-21615LL - 9075LL*r*xj - 300LL*Power(r, 2LL)*Power(xj, 2LL) +

            490LL*Power(r, 3LL)*Power(xj, 3LL) + 110LL*Power(r, 4LL)*Power(xj, 4LL) +

            8LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (90LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 9LL)*Power(xi + xj, 9LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_3S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_2S_3S(r, xj, xi);
}
