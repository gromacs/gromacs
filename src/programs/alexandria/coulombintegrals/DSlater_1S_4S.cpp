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

cl_R DSlater_1S_4S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-5088825LL*xi + 5806080LL*exp(2LL*r*xi)*xi - 8743140LL*r*Power(xi, 2LL) -

                  7319970LL*Power(r, 2LL)*Power(xi, 3LL) - 3946320LL*Power(r, 3LL)*Power(xi, 4LL) -

                  1519560LL*Power(r, 4LL)*Power(xi, 5LL) - 435456LL*Power(r, 5LL)*Power(xi, 6LL) -

                  92736LL*Power(r, 6LL)*Power(xi, 7LL) - 13824LL*Power(r, 7LL)*Power(xi, 8LL) -

                  1152LL*Power(r, 8LL)*Power(xi, 9LL))/(2.90304e6*exp(2LL*r*xi)*r) +

                (-2903040LL + 2903040LL*exp(2LL*r*xi) - 5088825LL*r*xi -

                 4371570LL*Power(r, 2LL)*Power(xi, 2LL) - 2439990LL*Power(r, 3LL)*Power(xi, 3LL) -

                 986580LL*Power(r, 4LL)*Power(xi, 4LL) - 303912LL*Power(r, 5LL)*Power(xi, 5LL) -

                 72576LL*Power(r, 6LL)*Power(xi, 6LL) - 13248LL*Power(r, 7LL)*Power(xi, 7LL) -

                 1728LL*Power(r, 8LL)*Power(xi, 8LL) - 128LL*Power(r, 9LL)*Power(xi, 9LL))/

                (2.90304e6*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-2903040LL + 2903040LL*exp(2LL*r*xi) - 5088825LL*r*xi -

                     4371570LL*Power(r, 2LL)*Power(xi, 2LL) - 2439990LL*Power(r, 3LL)*Power(xi, 3LL) -

                     986580LL*Power(r, 4LL)*Power(xi, 4LL) - 303912LL*Power(r, 5LL)*Power(xi, 5LL) -

                     72576LL*Power(r, 6LL)*Power(xi, 6LL) - 13248LL*Power(r, 7LL)*Power(xi, 7LL) -

                     1728LL*Power(r, 8LL)*Power(xi, 8LL) - 128LL*Power(r, 9LL)*Power(xi, 9LL)))/

                (1.45152e6*exp(2LL*r*xi)*r)

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
            S = (1260LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 9LL) +

                 1260LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-6LL*Power(xi, 8LL) - r*Power(xi, 9LL) - 51LL*Power(xi, 6LL)*Power(xj, 2LL) -

                  6LL*r*Power(xi, 7LL)*Power(xj, 2LL) - 63LL*Power(xi, 4LL)*Power(xj, 4LL) -

                  9LL*Power(xi, 2LL)*Power(xj, 6LL) + 6LL*r*Power(xi, 3LL)*Power(xj, 6LL) +

                  Power(xj, 8LL) + r*xi*Power(xj, 8LL)) +

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-42LL*Power(xi, 10LL)*Power(xj, 4LL)*

                  (1080LL + 1890LL*r*xj + 1620LL*Power(r, 2LL)*Power(xj, 2LL) +

            900LL*Power(r, 3LL)*Power(xj, 3LL) + 360LL*Power(r, 4LL)*Power(xj, 4LL) +

            111LL*Power(r, 5LL)*Power(xj, 5LL) + 22LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  70LL*Power(xi, 8LL)*Power(xj, 6LL)*

                  (1512LL + 2646LL*r*xj + 2268LL*Power(r, 2LL)*Power(xj, 2LL) +

            1248LL*Power(r, 3LL)*Power(xj, 3LL) + 528LL*Power(r, 4LL)*Power(xj, 4LL) +

            153LL*Power(r, 5LL)*Power(xj, 5LL) + 26LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  14LL*Power(xi, 2LL)*Power(xj, 12LL)*

                  (2970LL + 16335LL*r*xj + 15390LL*Power(r, 2LL)*Power(xj, 2LL) +

            7110LL*Power(r, 3LL)*Power(xj, 3LL) + 1980LL*Power(r, 4LL)*Power(xj, 4LL) +

            351LL*Power(r, 5LL)*Power(xj, 5LL) + 38LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  2LL*Power(xj, 14LL)*(62370LL + 72765LL*r*xj +

                                       39690LL*Power(r, 2LL)*Power(xj, 2LL) + 13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       2940LL*Power(r, 4LL)*Power(xj, 4LL) + 441LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       42LL*Power(r, 6LL)*Power(xj, 6LL) + 2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  Power(xi, 14LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  7LL*Power(xi, 12LL)*Power(xj, 2LL)*

                  (1620LL + 2835LL*r*xj + 2430LL*Power(r, 2LL)*Power(xj, 2LL) +

            1350LL*Power(r, 3LL)*Power(xj, 3LL) + 540LL*Power(r, 4LL)*Power(xj, 4LL) +

            162LL*Power(r, 5LL)*Power(xj, 5LL) + 36LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  35LL*Power(xi, 6LL)*Power(xj, 8LL)*

                  (4536LL + 7983LL*r*xj + 6534LL*Power(r, 2LL)*Power(xj, 2LL) +

            4014LL*Power(r, 3LL)*Power(xj, 3LL) + 1644LL*Power(r, 4LL)*Power(xj, 4LL) +

            414LL*Power(r, 5LL)*Power(xj, 5LL) + 60LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  21LL*Power(xi, 4LL)*Power(xj, 10LL)*

                  (7920LL + 11385LL*r*xj + 12330LL*Power(r, 2LL)*Power(xj, 2LL) +

            7410LL*Power(r, 3LL)*Power(xj, 3LL) + 2580LL*Power(r, 4LL)*Power(xj, 4LL) +

            546LL*Power(r, 5LL)*Power(xj, 5LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1260LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 9LL)*

                 Power(xi + xj, 9LL)) + (1260LL*exp(2LL*r*(xi + xj))*

                                         Power(Power(xi, 2LL) - Power(xj, 2LL), 9LL) +

                                         1260LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                                         (-6LL*Power(xi, 8LL) - r*Power(xi, 9LL) - 51LL*Power(xi, 6LL)*Power(xj, 2LL) -

                                 6LL*r*Power(xi, 7LL)*Power(xj, 2LL) - 63LL*Power(xi, 4LL)*Power(xj, 4LL) -

                                 9LL*Power(xi, 2LL)*Power(xj, 6LL) + 6LL*r*Power(xi, 3LL)*Power(xj, 6LL) +

                                 Power(xj, 8LL) + r*xi*Power(xj, 8LL)) +

                                         exp(2LL*r*xi)*Power(xi, 4LL)*

                                         (-42LL*Power(xi, 10LL)*Power(xj, 4LL)*

                                 (1080LL + 1890LL*r*xj + 1620LL*Power(r, 2LL)*Power(xj, 2LL) +

            900LL*Power(r, 3LL)*Power(xj, 3LL) + 360LL*Power(r, 4LL)*Power(xj, 4LL) +

            111LL*Power(r, 5LL)*Power(xj, 5LL) + 22LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                 70LL*Power(xi, 8LL)*Power(xj, 6LL)*

                                 (1512LL + 2646LL*r*xj + 2268LL*Power(r, 2LL)*Power(xj, 2LL) +

            1248LL*Power(r, 3LL)*Power(xj, 3LL) + 528LL*Power(r, 4LL)*Power(xj, 4LL) +

            153LL*Power(r, 5LL)*Power(xj, 5LL) + 26LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                 14LL*Power(xi, 2LL)*Power(xj, 12LL)*

                                 (2970LL + 16335LL*r*xj + 15390LL*Power(r, 2LL)*Power(xj, 2LL) +

            7110LL*Power(r, 3LL)*Power(xj, 3LL) + 1980LL*Power(r, 4LL)*Power(xj, 4LL) +

            351LL*Power(r, 5LL)*Power(xj, 5LL) + 38LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                 2LL*Power(xj, 14LL)*(62370LL + 72765LL*r*xj +

                                                      39690LL*Power(r, 2LL)*Power(xj, 2LL) + 13230LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                      2940LL*Power(r, 4LL)*Power(xj, 4LL) + 441LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                      42LL*Power(r, 6LL)*Power(xj, 6LL) + 2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                 Power(xi, 14LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                  4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                 7LL*Power(xi, 12LL)*Power(xj, 2LL)*

                                 (1620LL + 2835LL*r*xj + 2430LL*Power(r, 2LL)*Power(xj, 2LL) +

            1350LL*Power(r, 3LL)*Power(xj, 3LL) + 540LL*Power(r, 4LL)*Power(xj, 4LL) +

            162LL*Power(r, 5LL)*Power(xj, 5LL) + 36LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                 35LL*Power(xi, 6LL)*Power(xj, 8LL)*

                                 (4536LL + 7983LL*r*xj + 6534LL*Power(r, 2LL)*Power(xj, 2LL) +

            4014LL*Power(r, 3LL)*Power(xj, 3LL) + 1644LL*Power(r, 4LL)*Power(xj, 4LL) +

            414LL*Power(r, 5LL)*Power(xj, 5LL) + 60LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                 21LL*Power(xi, 4LL)*Power(xj, 10LL)*

                                 (7920LL + 11385LL*r*xj + 12330LL*Power(r, 2LL)*Power(xj, 2LL) +

            7410LL*Power(r, 3LL)*Power(xj, 3LL) + 2580LL*Power(r, 4LL)*Power(xj, 4LL) +

            546LL*Power(r, 5LL)*Power(xj, 5LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (630LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 9LL)*Power(xi + xj, 8LL)) -

                (2520LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 9LL) +

                 1260LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-Power(xi, 9LL) - 6LL*Power(xi, 7LL)*Power(xj, 2LL) +

        6LL*Power(xi, 3LL)*Power(xj, 6LL) + xi*Power(xj, 8LL)) +

                 2520LL*exp(2LL*r*xj)*Power(xj, 11LL)*

                 (-6LL*Power(xi, 8LL) - r*Power(xi, 9LL) - 51LL*Power(xi, 6LL)*Power(xj, 2LL) -

        6LL*r*Power(xi, 7LL)*Power(xj, 2LL) - 63LL*Power(xi, 4LL)*Power(xj, 4LL) -

        9LL*Power(xi, 2LL)*Power(xj, 6LL) + 6LL*r*Power(xi, 3LL)*Power(xj, 6LL) +

        Power(xj, 8LL) + r*xi*Power(xj, 8LL)) +

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (-42LL*Power(xi, 10LL)*Power(xj, 4LL)*

        (1890LL*xj + 3240LL*r*Power(xj, 2LL) + 2700LL*Power(r, 2LL)*Power(xj, 3LL) +

            1440LL*Power(r, 3LL)*Power(xj, 4LL) + 555LL*Power(r, 4LL)*Power(xj, 5LL) +

            132LL*Power(r, 5LL)*Power(xj, 6LL) + 14LL*Power(r, 6LL)*Power(xj, 7LL)) +

        70LL*Power(xi, 8LL)*Power(xj, 6LL)*

        (2646LL*xj + 4536LL*r*Power(xj, 2LL) + 3744LL*Power(r, 2LL)*Power(xj, 3LL) +

            2112LL*Power(r, 3LL)*Power(xj, 4LL) + 765LL*Power(r, 4LL)*Power(xj, 5LL) +

            156LL*Power(r, 5LL)*Power(xj, 6LL) + 14LL*Power(r, 6LL)*Power(xj, 7LL)) -

        14LL*Power(xi, 2LL)*Power(xj, 12LL)*

        (16335LL*xj + 30780LL*r*Power(xj, 2LL) + 21330LL*Power(r, 2LL)*Power(xj, 3LL) +

            7920LL*Power(r, 3LL)*Power(xj, 4LL) + 1755LL*Power(r, 4LL)*Power(xj, 5LL) +

            228LL*Power(r, 5LL)*Power(xj, 6LL) + 14LL*Power(r, 6LL)*Power(xj, 7LL)) +

        2LL*Power(xj, 14LL)*(72765LL*xj + 79380LL*r*Power(xj, 2LL) +

                             39690LL*Power(r, 2LL)*Power(xj, 3LL) + 11760LL*Power(r, 3LL)*Power(xj, 4LL) +

                             2205LL*Power(r, 4LL)*Power(xj, 5LL) + 252LL*Power(r, 5LL)*Power(xj, 6LL) +

                             14LL*Power(r, 6LL)*Power(xj, 7LL)) -

        Power(xi, 14LL)*(2205LL*xj + 3780LL*r*Power(xj, 2LL) +

                         3150LL*Power(r, 2LL)*Power(xj, 3LL) + 1680LL*Power(r, 3LL)*Power(xj, 4LL) +

                         630LL*Power(r, 4LL)*Power(xj, 5LL) + 168LL*Power(r, 5LL)*Power(xj, 6LL) +

                         28LL*Power(r, 6LL)*Power(xj, 7LL)) +

        7LL*Power(xi, 12LL)*Power(xj, 2LL)*

        (2835LL*xj + 4860LL*r*Power(xj, 2LL) + 4050LL*Power(r, 2LL)*Power(xj, 3LL) +

            2160LL*Power(r, 3LL)*Power(xj, 4LL) + 810LL*Power(r, 4LL)*Power(xj, 5LL) +

            216LL*Power(r, 5LL)*Power(xj, 6LL) + 28LL*Power(r, 6LL)*Power(xj, 7LL)) -

        35LL*Power(xi, 6LL)*Power(xj, 8LL)*

        (7983LL*xj + 13068LL*r*Power(xj, 2LL) + 12042LL*Power(r, 2LL)*Power(xj, 3LL) +

            6576LL*Power(r, 3LL)*Power(xj, 4LL) + 2070LL*Power(r, 4LL)*Power(xj, 5LL) +

            360LL*Power(r, 5LL)*Power(xj, 6LL) + 28LL*Power(r, 6LL)*Power(xj, 7LL)) +

        21LL*Power(xi, 4LL)*Power(xj, 10LL)*

        (11385LL*xj + 24660LL*r*Power(xj, 2LL) + 22230LL*Power(r, 2LL)*Power(xj, 3LL) +

            10320LL*Power(r, 3LL)*Power(xj, 4LL) + 2730LL*Power(r, 4LL)*Power(xj, 5LL) +

            408LL*Power(r, 5LL)*Power(xj, 6LL) + 28LL*Power(r, 6LL)*Power(xj, 7LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 5LL)*

                 (-42LL*Power(xi, 10LL)*Power(xj, 4LL)*

        (1080LL + 1890LL*r*xj + 1620LL*Power(r, 2LL)*Power(xj, 2LL) +

            900LL*Power(r, 3LL)*Power(xj, 3LL) + 360LL*Power(r, 4LL)*Power(xj, 4LL) +

            111LL*Power(r, 5LL)*Power(xj, 5LL) + 22LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

        70LL*Power(xi, 8LL)*Power(xj, 6LL)*

        (1512LL + 2646LL*r*xj + 2268LL*Power(r, 2LL)*Power(xj, 2LL) +

            1248LL*Power(r, 3LL)*Power(xj, 3LL) + 528LL*Power(r, 4LL)*Power(xj, 4LL) +

            153LL*Power(r, 5LL)*Power(xj, 5LL) + 26LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

        14LL*Power(xi, 2LL)*Power(xj, 12LL)*

        (2970LL + 16335LL*r*xj + 15390LL*Power(r, 2LL)*Power(xj, 2LL) +

            7110LL*Power(r, 3LL)*Power(xj, 3LL) + 1980LL*Power(r, 4LL)*Power(xj, 4LL) +

            351LL*Power(r, 5LL)*Power(xj, 5LL) + 38LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

        2LL*Power(xj, 14LL)*(62370LL + 72765LL*r*xj + 39690LL*Power(r, 2LL)*Power(xj, 2LL) +

                             13230LL*Power(r, 3LL)*Power(xj, 3LL) + 2940LL*Power(r, 4LL)*Power(xj, 4LL) +

                             441LL*Power(r, 5LL)*Power(xj, 5LL) + 42LL*Power(r, 6LL)*Power(xj, 6LL) +

                             2LL*Power(r, 7LL)*Power(xj, 7LL)) -

        Power(xi, 14LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                         1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                         126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                         4LL*Power(r, 7LL)*Power(xj, 7LL)) +

        7LL*Power(xi, 12LL)*Power(xj, 2LL)*

        (1620LL + 2835LL*r*xj + 2430LL*Power(r, 2LL)*Power(xj, 2LL) +

            1350LL*Power(r, 3LL)*Power(xj, 3LL) + 540LL*Power(r, 4LL)*Power(xj, 4LL) +

            162LL*Power(r, 5LL)*Power(xj, 5LL) + 36LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) -

        35LL*Power(xi, 6LL)*Power(xj, 8LL)*

        (4536LL + 7983LL*r*xj + 6534LL*Power(r, 2LL)*Power(xj, 2LL) +

            4014LL*Power(r, 3LL)*Power(xj, 3LL) + 1644LL*Power(r, 4LL)*Power(xj, 4LL) +

            414LL*Power(r, 5LL)*Power(xj, 5LL) + 60LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) +

        21LL*Power(xi, 4LL)*Power(xj, 10LL)*

        (7920LL + 11385LL*r*xj + 12330LL*Power(r, 2LL)*Power(xj, 2LL) +

            7410LL*Power(r, 3LL)*Power(xj, 3LL) + 2580LL*Power(r, 4LL)*Power(xj, 4LL) +

            546LL*Power(r, 5LL)*Power(xj, 5LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1260LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 9LL)*Power(xi + xj, 9LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_4S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_1S_4S(r, xj, xi);
}
