/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#include "gromacs/coulombintegrals/slater_low.h"

#if HAVE_LIBCLN
cl_R DSlater_3S_3S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-2503064025LL*xi + 2874009600LL*exp(2LL*r*xi)*xi -

                  4264236900LL*r*Power(xi, 2LL) - 3541992300LL*Power(r, 2LL)*Power(xi, 3LL) -

                  1906027200LL*Power(r, 3LL)*Power(xi, 4LL) -

                  744282000LL*Power(r, 4LL)*Power(xi, 5LL) - 223534080LL*Power(r, 5LL)*Power(xi, 6LL) -

                  53222400LL*Power(r, 6LL)*Power(xi, 7LL) - 10137600LL*Power(r, 7LL)*Power(xi, 8LL) -

                  1520640LL*Power(r, 8LL)*Power(xi, 9LL) - 168960LL*Power(r, 9LL)*Power(xi, 10LL) -

                  11264LL*Power(r, 10LL)*Power(xi, 11LL))/(1.4370048e9*exp(2LL*r*xi)*r) +

                (-1437004800LL + 1437004800LL*exp(2LL*r*xi) - 2503064025LL*r*xi -

                 2132118450LL*Power(r, 2LL)*Power(xi, 2LL) -

                 1180664100LL*Power(r, 3LL)*Power(xi, 3LL) - 476506800LL*Power(r, 4LL)*Power(xi, 4LL) -

                 148856400LL*Power(r, 5LL)*Power(xi, 5LL) - 37255680LL*Power(r, 6LL)*Power(xi, 6LL) -

                 7603200LL*Power(r, 7LL)*Power(xi, 7LL) - 1267200LL*Power(r, 8LL)*Power(xi, 8LL) -

                 168960LL*Power(r, 9LL)*Power(xi, 9LL) - 16896LL*Power(r, 10LL)*Power(xi, 10LL) -

                 1024LL*Power(r, 11LL)*Power(xi, 11LL))/(1.4370048e9*exp(2LL*r*xi)*Power(r, 2LL))

                + (xi*(-1437004800LL + 1437004800LL*exp(2LL*r*xi) - 2503064025LL*r*xi -

                       2132118450LL*Power(r, 2LL)*Power(xi, 2LL) -

                       1180664100LL*Power(r, 3LL)*Power(xi, 3LL) -

                       476506800LL*Power(r, 4LL)*Power(xi, 4LL) - 148856400LL*Power(r, 5LL)*Power(xi, 5LL) -

                       37255680LL*Power(r, 6LL)*Power(xi, 6LL) - 7603200LL*Power(r, 7LL)*Power(xi, 7LL) -

                       1267200LL*Power(r, 8LL)*Power(xi, 8LL) - 168960LL*Power(r, 9LL)*Power(xi, 9LL) -

                       16896LL*Power(r, 10LL)*Power(xi, 10LL) - 1024LL*Power(r, 11LL)*Power(xi, 11LL)))/

                (7.185024e8*exp(2LL*r*xi)*r)

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
            S = (135LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                 exp(2LL*r*xj)*Power(xj, 8LL)*

                 (-150LL*Power(r, 4LL)*Power(xi, 18LL) - 6LL*Power(r, 5LL)*Power(xi, 19LL) +

                  135LL*Power(xj, 14LL) + 225LL*r*xi*Power(xj, 14LL) +

                  10LL*Power(r, 3LL)*Power(xi, 17LL)*(-165LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  30LL*Power(r, 2LL)*Power(xi, 16LL)*(330LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*r*Power(xi, 3LL)*Power(xj, 12LL)*(-55LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*Power(xi, 2LL)*Power(xj, 12LL)*(-33LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 9LL)*Power(xj, 6LL)*

                  (234135LL - 4950LL*Power(r, 2LL)*Power(xj, 2LL) -

                   34LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  5LL*r*Power(xi, 7LL)*Power(xj, 8LL)*

                  (6237LL - 1242LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*

                  (4125LL - 330LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 4LL)*Power(xj, 10LL)*

                  (495LL - 132LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  165LL*Power(xi, 6LL)*Power(xj, 8LL)*

                  (135LL - 60LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  5LL*r*Power(xi, 13LL)*Power(xj, 2LL)*

                  (43875LL - 3438LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 5LL*r*Power(xi, 11LL)*Power(xj, 4LL)*

                  (7695LL - 2442LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 15LL*Power(xi, 8LL)*Power(xj, 6LL)*(-33LL - 3564LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        26LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  r*Power(xi, 15LL)*(-32175LL - 3690LL*Power(r, 2LL)*Power(xj, 2LL) +

                                     34LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 10LL)*Power(xj, 4LL)*

                  (-32277LL + 1364LL*Power(r, 2LL)*Power(xj, 2LL) +

                   66LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 14LL)*(-3003LL - 2932LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        94LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*Power(xi, 12LL)*Power(xj, 2LL)*

                  (28119LL - 5252LL*Power(r, 2LL)*Power(xj, 2LL) + 154LL*Power(r, 4LL)*Power(xj, 4LL))

                 ) + exp(2LL*r*xi)*Power(xi, 8LL)*

                 (-5LL*Power(xi, 2LL)*Power(xj, 12LL)*

                  (-84357LL - 43875LL*r*xj - 8796LL*Power(r, 2LL)*Power(xj, 2LL) -

                   738LL*Power(r, 3LL)*Power(xj, 3LL) - 6LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  3LL*Power(xi, 14LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  55LL*Power(xi, 8LL)*Power(xj, 6LL)*

                  (-405LL - 567LL*r*xj - 972LL*Power(r, 2LL)*Power(xj, 2LL) -

                   90LL*Power(r, 3LL)*Power(xj, 3LL) + 18LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  55LL*Power(xi, 6LL)*Power(xj, 8LL)*

                  (9LL - 4257LL*r*xj - 372LL*Power(r, 2LL)*Power(xj, 2LL) +

                   222LL*Power(r, 3LL)*Power(xj, 3LL) + 42LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  3LL*Power(xj, 14LL)*(15015LL + 10725LL*r*xj + 3300LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       550LL*Power(r, 3LL)*Power(xj, 3LL) + 50LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  5LL*Power(xi, 12LL)*Power(xj, 2LL)*

                  (297LL + 495LL*r*xj + 396LL*Power(r, 2LL)*Power(xj, 2LL) +

                   198LL*Power(r, 3LL)*Power(xj, 3LL) + 66LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  Power(xi, 10LL)*Power(xj, 4LL)*

                  (-7425LL - 12375LL*r*xj - 9900LL*Power(r, 2LL)*Power(xj, 2LL) -

                   6210LL*Power(r, 3LL)*Power(xj, 3LL) - 390LL*Power(r, 4LL)*Power(xj, 4LL) +

                   34LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  Power(xi, 4LL)*Power(xj, 10LL)*

                  (-484155LL + 38475LL*r*xj + 78780LL*Power(r, 2LL)*Power(xj, 2LL) +

                   17190LL*Power(r, 3LL)*Power(xj, 3LL) + 1410LL*Power(r, 4LL)*Power(xj, 4LL) +

                   34LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (135LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 11LL)*

                 Power(xi + xj, 11LL)) + (2LL*(135LL*exp(2LL*r*(xi + xj))*

                                               Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                                               exp(2LL*r*xj)*Power(xj, 8LL)*

                                               (-150LL*Power(r, 4LL)*Power(xi, 18LL) - 6LL*Power(r, 5LL)*Power(xi, 19LL) +

                                                135LL*Power(xj, 14LL) + 225LL*r*xi*Power(xj, 14LL) +

                                                10LL*Power(r, 3LL)*Power(xi, 17LL)*(-165LL + Power(r, 2LL)*Power(xj, 2LL)) -

                                                30LL*Power(r, 2LL)*Power(xi, 16LL)*(330LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                                45LL*r*Power(xi, 3LL)*Power(xj, 12LL)*(-55LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                                45LL*Power(xi, 2LL)*Power(xj, 12LL)*(-33LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                                r*Power(xi, 9LL)*Power(xj, 6LL)*

                                                (234135LL - 4950LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                 34LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                                5LL*r*Power(xi, 7LL)*Power(xj, 8LL)*

                                                (6237LL - 1242LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

                                                + 3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*

                                                (4125LL - 330LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL))

                                                + 15LL*Power(xi, 4LL)*Power(xj, 10LL)*

                                                (495LL - 132LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                                165LL*Power(xi, 6LL)*Power(xj, 8LL)*

                                                (135LL - 60LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                                5LL*r*Power(xi, 13LL)*Power(xj, 2LL)*

                                                (43875LL - 3438LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 22LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                5LL*r*Power(xi, 11LL)*Power(xj, 4LL)*

                                                (7695LL - 2442LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 22LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                15LL*Power(xi, 8LL)*Power(xj, 6LL)*

                                                (-33LL - 3564LL*Power(r, 2LL)*Power(xj, 2LL) + 26LL*Power(r, 4LL)*Power(xj, 4LL))

                                                + r*Power(xi, 15LL)*(-32175LL - 3690LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                     34LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                15LL*Power(xi, 10LL)*Power(xj, 4LL)*

                                                (-32277LL + 1364LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 66LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                15LL*Power(xi, 14LL)*(-3003LL - 2932LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                      94LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                                15LL*Power(xi, 12LL)*Power(xj, 2LL)*

                                                (28119LL - 5252LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 154LL*Power(r, 4LL)*Power(xj, 4LL))) +

                                               exp(2LL*r*xi)*Power(xi, 8LL)*

                                               (-5LL*Power(xi, 2LL)*Power(xj, 12LL)*

                                                (-84357LL - 43875LL*r*xj - 8796LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                 738LL*Power(r, 3LL)*Power(xj, 3LL) - 6LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                 2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                                                3LL*Power(xi, 14LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                     30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                     2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                                                55LL*Power(xi, 8LL)*Power(xj, 6LL)*

                                                (-405LL - 567LL*r*xj - 972LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                 90LL*Power(r, 3LL)*Power(xj, 3LL) + 18LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                 2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                                                55LL*Power(xi, 6LL)*Power(xj, 8LL)*

                                                (9LL - 4257LL*r*xj - 372LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 222LL*Power(r, 3LL)*Power(xj, 3LL) + 42LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                 2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                                                3LL*Power(xj, 14LL)*(15015LL + 10725LL*r*xj +

                                                                     3300LL*Power(r, 2LL)*Power(xj, 2LL) + 550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                     50LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                                                5LL*Power(xi, 12LL)*Power(xj, 2LL)*

                                                (297LL + 495LL*r*xj + 396LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 198LL*Power(r, 3LL)*Power(xj, 3LL) + 66LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                 2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                                                Power(xi, 10LL)*Power(xj, 4LL)*

                                                (-7425LL - 12375LL*r*xj - 9900LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                 6210LL*Power(r, 3LL)*Power(xj, 3LL) - 390LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                 34LL*Power(r, 5LL)*Power(xj, 5LL)) -

                                                Power(xi, 4LL)*Power(xj, 10LL)*

                                                (-484155LL + 38475LL*r*xj + 78780LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                 17190LL*Power(r, 3LL)*Power(xj, 3LL) + 1410LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                 34LL*Power(r, 5LL)*Power(xj, 5LL)))))/

                (135LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 11LL)*Power(xi + xj, 10LL)) -

                (270LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                 exp(2LL*r*xj)*Power(xj, 8LL)*

                 (-600LL*Power(r, 3LL)*Power(xi, 18LL) - 30LL*Power(r, 4LL)*Power(xi, 19LL) -

                  60LL*Power(r, 3LL)*Power(xi, 16LL)*Power(xj, 2LL) +

                  20LL*Power(r, 4LL)*Power(xi, 17LL)*Power(xj, 2LL) + 225LL*xi*Power(xj, 14LL) +

                  360LL*r*Power(xi, 2LL)*Power(xj, 14LL) +

                  180LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 14LL) +

                  30LL*Power(r, 2LL)*Power(xi, 17LL)*(-165LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  60LL*r*Power(xi, 16LL)*(330LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*Power(xi, 3LL)*Power(xj, 12LL)*(-55LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 9LL)*Power(xj, 6LL)*

                  (-9900LL*r*Power(xj, 2LL) - 136LL*Power(r, 3LL)*Power(xj, 4LL)) -

                  5LL*r*Power(xi, 7LL)*Power(xj, 8LL)*

                  (-2484LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*

                  (-660LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 4LL)*Power(xj, 10LL)*

                  (-264LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) -

                  165LL*Power(xi, 6LL)*Power(xj, 8LL)*

                  (-120LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) -

                  5LL*r*Power(xi, 13LL)*Power(xj, 2LL)*

                  (-6876LL*r*Power(xj, 2LL) + 88LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  5LL*r*Power(xi, 11LL)*Power(xj, 4LL)*

                  (-4884LL*r*Power(xj, 2LL) + 88LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 8LL)*Power(xj, 6LL)*

                  (-7128LL*r*Power(xj, 2LL) + 104LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  r*Power(xi, 15LL)*(-7380LL*r*Power(xj, 2LL) + 136LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 10LL)*Power(xj, 4LL)*

                  (2728LL*r*Power(xj, 2LL) + 264LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 14LL)*(-5864LL*r*Power(xj, 2LL) +

                                        376LL*Power(r, 3LL)*Power(xj, 4LL)) -

                  15LL*Power(xi, 12LL)*Power(xj, 2LL)*

                  (-10504LL*r*Power(xj, 2LL) + 616LL*Power(r, 3LL)*Power(xj, 4LL)) +

                  Power(xi, 9LL)*Power(xj, 6LL)*

                  (234135LL - 4950LL*Power(r, 2LL)*Power(xj, 2LL) - 34LL*Power(r, 4LL)*Power(xj, 4LL))

                  - 5LL*Power(xi, 7LL)*Power(xj, 8LL)*(6237LL - 1242LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  3LL*Power(xi, 5LL)*Power(xj, 10LL)*

                  (4125LL - 330LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  5LL*Power(xi, 13LL)*Power(xj, 2LL)*

                  (43875LL - 3438LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 5LL*Power(xi, 11LL)*Power(xj, 4LL)*(7695LL - 2442LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        22LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  Power(xi, 15LL)*(-32175LL - 3690LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   34LL*Power(r, 4LL)*Power(xj, 4LL))) +

                 2LL*exp(2LL*r*xj)*Power(xj, 9LL)*

                 (-150LL*Power(r, 4LL)*Power(xi, 18LL) - 6LL*Power(r, 5LL)*Power(xi, 19LL) +

                  135LL*Power(xj, 14LL) + 225LL*r*xi*Power(xj, 14LL) +

                  10LL*Power(r, 3LL)*Power(xi, 17LL)*(-165LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  30LL*Power(r, 2LL)*Power(xi, 16LL)*(330LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*r*Power(xi, 3LL)*Power(xj, 12LL)*(-55LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*Power(xi, 2LL)*Power(xj, 12LL)*(-33LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 9LL)*Power(xj, 6LL)*

                  (234135LL - 4950LL*Power(r, 2LL)*Power(xj, 2LL) - 34LL*Power(r, 4LL)*Power(xj, 4LL))

                  - 5LL*r*Power(xi, 7LL)*Power(xj, 8LL)*(6237LL - 1242LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*

                  (4125LL - 330LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 4LL)*Power(xj, 10LL)*

                  (495LL - 132LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  165LL*Power(xi, 6LL)*Power(xj, 8LL)*

                  (135LL - 60LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  5LL*r*Power(xi, 13LL)*Power(xj, 2LL)*

                  (43875LL - 3438LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 5LL*r*Power(xi, 11LL)*Power(xj, 4LL)*

                  (7695LL - 2442LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 8LL)*Power(xj, 6LL)*

                  (-33LL - 3564LL*Power(r, 2LL)*Power(xj, 2LL) + 26LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  r*Power(xi, 15LL)*(-32175LL - 3690LL*Power(r, 2LL)*Power(xj, 2LL) +

                                     34LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  15LL*Power(xi, 10LL)*Power(xj, 4LL)*

                  (-32277LL + 1364LL*Power(r, 2LL)*Power(xj, 2LL) + 66LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 15LL*Power(xi, 14LL)*(-3003LL - 2932LL*Power(r, 2LL)*Power(xj, 2LL) +

                                          94LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*Power(xi, 12LL)*Power(xj, 2LL)*

                  (28119LL - 5252LL*Power(r, 2LL)*Power(xj, 2LL) + 154LL*Power(r, 4LL)*Power(xj, 4LL)))

                 + exp(2LL*r*xi)*Power(xi, 8LL)*(-5LL*Power(xi, 2LL)*Power(xj, 12LL)*

                                                 (-43875LL*xj - 17592LL*r*Power(xj, 2LL) - 2214LL*Power(r, 2LL)*Power(xj, 3LL) -

                                                  24LL*Power(r, 3LL)*Power(xj, 4LL) + 10LL*Power(r, 4LL)*Power(xj, 5LL)) -

                                                 3LL*Power(xi, 14LL)*(75LL*xj + 120LL*r*Power(xj, 2LL) +

                                                                      90LL*Power(r, 2LL)*Power(xj, 3LL) + 40LL*Power(r, 3LL)*Power(xj, 4LL) +

                                                                      10LL*Power(r, 4LL)*Power(xj, 5LL)) -

                                                 55LL*Power(xi, 8LL)*Power(xj, 6LL)*

                                                 (-567LL*xj - 1944LL*r*Power(xj, 2LL) - 270LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                  72LL*Power(r, 3LL)*Power(xj, 4LL) + 10LL*Power(r, 4LL)*Power(xj, 5LL)) +

                                                 55LL*Power(xi, 6LL)*Power(xj, 8LL)*

                                                 (-4257LL*xj - 744LL*r*Power(xj, 2LL) + 666LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                  168LL*Power(r, 3LL)*Power(xj, 4LL) + 10LL*Power(r, 4LL)*Power(xj, 5LL)) +

                                                 3LL*Power(xj, 14LL)*(10725LL*xj + 6600LL*r*Power(xj, 2LL) +

                                                                      1650LL*Power(r, 2LL)*Power(xj, 3LL) + 200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                                                      10LL*Power(r, 4LL)*Power(xj, 5LL)) +

                                                 5LL*Power(xi, 12LL)*Power(xj, 2LL)*

                                                 (495LL*xj + 792LL*r*Power(xj, 2LL) + 594LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                  264LL*Power(r, 3LL)*Power(xj, 4LL) + 10LL*Power(r, 4LL)*Power(xj, 5LL)) +

                                                 Power(xi, 10LL)*Power(xj, 4LL)*

                                                 (-12375LL*xj - 19800LL*r*Power(xj, 2LL) - 18630LL*Power(r, 2LL)*Power(xj, 3LL) -

                                                  1560LL*Power(r, 3LL)*Power(xj, 4LL) + 170LL*Power(r, 4LL)*Power(xj, 5LL)) -

                                                 Power(xi, 4LL)*Power(xj, 10LL)*

                                                 (38475LL*xj + 157560LL*r*Power(xj, 2LL) + 51570LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                  5640LL*Power(r, 3LL)*Power(xj, 4LL) + 170LL*Power(r, 4LL)*Power(xj, 5LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 9LL)*

                 (-5LL*Power(xi, 2LL)*Power(xj, 12LL)*

                  (-84357LL - 43875LL*r*xj - 8796LL*Power(r, 2LL)*Power(xj, 2LL) -

                   738LL*Power(r, 3LL)*Power(xj, 3LL) - 6LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  3LL*Power(xi, 14LL)*(45LL + 75LL*r*xj + 60LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       30LL*Power(r, 3LL)*Power(xj, 3LL) + 10LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       2LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  55LL*Power(xi, 8LL)*Power(xj, 6LL)*

                  (-405LL - 567LL*r*xj - 972LL*Power(r, 2LL)*Power(xj, 2LL) -

                   90LL*Power(r, 3LL)*Power(xj, 3LL) + 18LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  55LL*Power(xi, 6LL)*Power(xj, 8LL)*

                  (9LL - 4257LL*r*xj - 372LL*Power(r, 2LL)*Power(xj, 2LL) +

                   222LL*Power(r, 3LL)*Power(xj, 3LL) + 42LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  3LL*Power(xj, 14LL)*(15015LL + 10725LL*r*xj + 3300LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       550LL*Power(r, 3LL)*Power(xj, 3LL) + 50LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  5LL*Power(xi, 12LL)*Power(xj, 2LL)*

                  (297LL + 495LL*r*xj + 396LL*Power(r, 2LL)*Power(xj, 2LL) +

                   198LL*Power(r, 3LL)*Power(xj, 3LL) + 66LL*Power(r, 4LL)*Power(xj, 4LL) +

                   2LL*Power(r, 5LL)*Power(xj, 5LL)) +

                  Power(xi, 10LL)*Power(xj, 4LL)*

                  (-7425LL - 12375LL*r*xj - 9900LL*Power(r, 2LL)*Power(xj, 2LL) -

                   6210LL*Power(r, 3LL)*Power(xj, 3LL) - 390LL*Power(r, 4LL)*Power(xj, 4LL) +

                   34LL*Power(r, 5LL)*Power(xj, 5LL)) -

                  Power(xi, 4LL)*Power(xj, 10LL)*

                  (-484155LL + 38475LL*r*xj + 78780LL*Power(r, 2LL)*Power(xj, 2LL) +

                   17190LL*Power(r, 3LL)*Power(xj, 3LL) + 1410LL*Power(r, 4LL)*Power(xj, 4LL) +

                   34LL*Power(r, 5LL)*Power(xj, 5LL))))/

                (135LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 11LL)*Power(xi + xj, 11LL))

            ;
        }

    }
    return S;
}

#else

double DSlater_3S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = -(-2503064025*xi + 2874009600*exp(2*r*xi)*xi -

                  4264236900*r*pow(xi, 2) - 3541992300*pow(r, 2)*pow(xi, 3) -

                  1906027200*pow(r, 3)*pow(xi, 4) -

                  744282000*pow(r, 4)*pow(xi, 5) - 223534080*pow(r, 5)*pow(xi, 6) -

                  53222400*pow(r, 6)*pow(xi, 7) - 10137600*pow(r, 7)*pow(xi, 8) -

                  1520640*pow(r, 8)*pow(xi, 9) - 168960*pow(r, 9)*pow(xi, 10) -

                  11264*pow(r, 10)*pow(xi, 11))/(1.4370048e9*exp(2*r*xi)*r) +

                (-1437004800 + 1437004800*exp(2*r*xi) - 2503064025*r*xi -

                 2132118450*pow(r, 2)*pow(xi, 2) -

                 1180664100*pow(r, 3)*pow(xi, 3) - 476506800*pow(r, 4)*pow(xi, 4) -

                 148856400*pow(r, 5)*pow(xi, 5) - 37255680*pow(r, 6)*pow(xi, 6) -

                 7603200*pow(r, 7)*pow(xi, 7) - 1267200*pow(r, 8)*pow(xi, 8) -

                 168960*pow(r, 9)*pow(xi, 9) - 16896*pow(r, 10)*pow(xi, 10) -

                 1024*pow(r, 11)*pow(xi, 11))/(1.4370048e9*exp(2*r*xi)*pow(r, 2))

                + (xi*(-1437004800 + 1437004800*exp(2*r*xi) - 2503064025*r*xi -

                       2132118450*pow(r, 2)*pow(xi, 2) -

                       1180664100*pow(r, 3)*pow(xi, 3) -

                       476506800*pow(r, 4)*pow(xi, 4) - 148856400*pow(r, 5)*pow(xi, 5) -

                       37255680*pow(r, 6)*pow(xi, 6) - 7603200*pow(r, 7)*pow(xi, 7) -

                       1267200*pow(r, 8)*pow(xi, 8) - 168960*pow(r, 9)*pow(xi, 9) -

                       16896*pow(r, 10)*pow(xi, 10) - 1024*pow(r, 11)*pow(xi, 11)))/

                (7.185024e8*exp(2*r*xi)*r)

            ;
        }

    }
    else
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = (135*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 11) +

                 exp(2*r*xj)*pow(xj, 8)*

                 (-150*pow(r, 4)*pow(xi, 18) - 6*pow(r, 5)*pow(xi, 19) +

                  135*pow(xj, 14) + 225*r*xi*pow(xj, 14) +

                  10*pow(r, 3)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                  30*pow(r, 2)*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                  45*r*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                  45*pow(xi, 2)*pow(xj, 12)*(-33 + 4*pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 9)*pow(xj, 6)*

                  (234135 - 4950*pow(r, 2)*pow(xj, 2) -

                   34*pow(r, 4)*pow(xj, 4)) -

                  5*r*pow(xi, 7)*pow(xj, 8)*

                  (6237 - 1242*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) +

                  3*r*pow(xi, 5)*pow(xj, 10)*

                  (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 4)*pow(xj, 10)*

                  (495 - 132*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  165*pow(xi, 6)*pow(xj, 8)*

                  (135 - 60*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  5*r*pow(xi, 13)*pow(xj, 2)*

                  (43875 - 3438*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 5*r*pow(xi, 11)*pow(xj, 4)*

                  (7695 - 2442*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 15*pow(xi, 8)*pow(xj, 6)*(-33 - 3564*pow(r, 2)*pow(xj, 2) +

                                                        26*pow(r, 4)*pow(xj, 4)) +

                  r*pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                     34*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 10)*pow(xj, 4)*

                  (-32277 + 1364*pow(r, 2)*pow(xj, 2) +

                   66*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 14)*(-3003 - 2932*pow(r, 2)*pow(xj, 2) +

                                        94*pow(r, 4)*pow(xj, 4)) -

                  15*pow(xi, 12)*pow(xj, 2)*

                  (28119 - 5252*pow(r, 2)*pow(xj, 2) + 154*pow(r, 4)*pow(xj, 4))

                 ) + exp(2*r*xi)*pow(xi, 8)*

                 (-5*pow(xi, 2)*pow(xj, 12)*

                  (-84357 - 43875*r*xj - 8796*pow(r, 2)*pow(xj, 2) -

                   738*pow(r, 3)*pow(xj, 3) - 6*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  3*pow(xi, 14)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                       30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) -

                  55*pow(xi, 8)*pow(xj, 6)*

                  (-405 - 567*r*xj - 972*pow(r, 2)*pow(xj, 2) -

                   90*pow(r, 3)*pow(xj, 3) + 18*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  55*pow(xi, 6)*pow(xj, 8)*

                  (9 - 4257*r*xj - 372*pow(r, 2)*pow(xj, 2) +

                   222*pow(r, 3)*pow(xj, 3) + 42*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  3*pow(xj, 14)*(15015 + 10725*r*xj + 3300*pow(r, 2)*pow(xj, 2) +

                                       550*pow(r, 3)*pow(xj, 3) + 50*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 12)*pow(xj, 2)*

                  (297 + 495*r*xj + 396*pow(r, 2)*pow(xj, 2) +

                   198*pow(r, 3)*pow(xj, 3) + 66*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 10)*pow(xj, 4)*

                  (-7425 - 12375*r*xj - 9900*pow(r, 2)*pow(xj, 2) -

                   6210*pow(r, 3)*pow(xj, 3) - 390*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 4)*pow(xj, 10)*

                  (-484155 + 38475*r*xj + 78780*pow(r, 2)*pow(xj, 2) +

                   17190*pow(r, 3)*pow(xj, 3) + 1410*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5))))/

                (135*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 11)*

                 pow(xi + xj, 11)) + (2*(135*exp(2*r*(xi + xj))*

                                               pow(pow(xi, 2) - pow(xj, 2), 11) +

                                               exp(2*r*xj)*pow(xj, 8)*

                                               (-150*pow(r, 4)*pow(xi, 18) - 6*pow(r, 5)*pow(xi, 19) +

                                                135*pow(xj, 14) + 225*r*xi*pow(xj, 14) +

                                                10*pow(r, 3)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                                                30*pow(r, 2)*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                                                45*r*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                                                45*pow(xi, 2)*pow(xj, 12)*(-33 + 4*pow(r, 2)*pow(xj, 2)) +

                                                r*pow(xi, 9)*pow(xj, 6)*

                                                (234135 - 4950*pow(r, 2)*pow(xj, 2) -

                                                 34*pow(r, 4)*pow(xj, 4)) -

                                                5*r*pow(xi, 7)*pow(xj, 8)*

                                                (6237 - 1242*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4))

                                                + 3*r*pow(xi, 5)*pow(xj, 10)*

                                                (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4))

                                                + 15*pow(xi, 4)*pow(xj, 10)*

                                                (495 - 132*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                                                165*pow(xi, 6)*pow(xj, 8)*

                                                (135 - 60*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                                                5*r*pow(xi, 13)*pow(xj, 2)*

                                                (43875 - 3438*pow(r, 2)*pow(xj, 2) +

                                                 22*pow(r, 4)*pow(xj, 4)) +

                                                5*r*pow(xi, 11)*pow(xj, 4)*

                                                (7695 - 2442*pow(r, 2)*pow(xj, 2) +

                                                 22*pow(r, 4)*pow(xj, 4)) +

                                                15*pow(xi, 8)*pow(xj, 6)*

                                                (-33 - 3564*pow(r, 2)*pow(xj, 2) + 26*pow(r, 4)*pow(xj, 4))

                                                + r*pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                                                     34*pow(r, 4)*pow(xj, 4)) +

                                                15*pow(xi, 10)*pow(xj, 4)*

                                                (-32277 + 1364*pow(r, 2)*pow(xj, 2) +

                                                 66*pow(r, 4)*pow(xj, 4)) +

                                                15*pow(xi, 14)*(-3003 - 2932*pow(r, 2)*pow(xj, 2) +

                                                                      94*pow(r, 4)*pow(xj, 4)) -

                                                15*pow(xi, 12)*pow(xj, 2)*

                                                (28119 - 5252*pow(r, 2)*pow(xj, 2) +

                                                 154*pow(r, 4)*pow(xj, 4))) +

                                               exp(2*r*xi)*pow(xi, 8)*

                                               (-5*pow(xi, 2)*pow(xj, 12)*

                                                (-84357 - 43875*r*xj - 8796*pow(r, 2)*pow(xj, 2) -

                                                 738*pow(r, 3)*pow(xj, 3) - 6*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) -

                                                3*pow(xi, 14)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                                                     30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                                                     2*pow(r, 5)*pow(xj, 5)) -

                                                55*pow(xi, 8)*pow(xj, 6)*

                                                (-405 - 567*r*xj - 972*pow(r, 2)*pow(xj, 2) -

                                                 90*pow(r, 3)*pow(xj, 3) + 18*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) +

                                                55*pow(xi, 6)*pow(xj, 8)*

                                                (9 - 4257*r*xj - 372*pow(r, 2)*pow(xj, 2) +

                                                 222*pow(r, 3)*pow(xj, 3) + 42*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) +

                                                3*pow(xj, 14)*(15015 + 10725*r*xj +

                                                                     3300*pow(r, 2)*pow(xj, 2) + 550*pow(r, 3)*pow(xj, 3) +

                                                                     50*pow(r, 4)*pow(xj, 4) + 2*pow(r, 5)*pow(xj, 5)) +

                                                5*pow(xi, 12)*pow(xj, 2)*

                                                (297 + 495*r*xj + 396*pow(r, 2)*pow(xj, 2) +

                                                 198*pow(r, 3)*pow(xj, 3) + 66*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) +

                                                pow(xi, 10)*pow(xj, 4)*

                                                (-7425 - 12375*r*xj - 9900*pow(r, 2)*pow(xj, 2) -

                                                 6210*pow(r, 3)*pow(xj, 3) - 390*pow(r, 4)*pow(xj, 4) +

                                                 34*pow(r, 5)*pow(xj, 5)) -

                                                pow(xi, 4)*pow(xj, 10)*

                                                (-484155 + 38475*r*xj + 78780*pow(r, 2)*pow(xj, 2) +

                                                 17190*pow(r, 3)*pow(xj, 3) + 1410*pow(r, 4)*pow(xj, 4) +

                                                 34*pow(r, 5)*pow(xj, 5)))))/

                (135*exp(2*r*(xi + xj))*r*pow(xi - xj, 11)*pow(xi + xj, 10)) -

                (270*exp(2*r*(xi + xj))*(xi + xj)*

                 pow(pow(xi, 2) - pow(xj, 2), 11) +

                 exp(2*r*xj)*pow(xj, 8)*

                 (-600*pow(r, 3)*pow(xi, 18) - 30*pow(r, 4)*pow(xi, 19) -

                  60*pow(r, 3)*pow(xi, 16)*pow(xj, 2) +

                  20*pow(r, 4)*pow(xi, 17)*pow(xj, 2) + 225*xi*pow(xj, 14) +

                  360*r*pow(xi, 2)*pow(xj, 14) +

                  180*pow(r, 2)*pow(xi, 3)*pow(xj, 14) +

                  30*pow(r, 2)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                  60*r*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                  45*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 9)*pow(xj, 6)*

                  (-9900*r*pow(xj, 2) - 136*pow(r, 3)*pow(xj, 4)) -

                  5*r*pow(xi, 7)*pow(xj, 8)*

                  (-2484*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) +

                  3*r*pow(xi, 5)*pow(xj, 10)*

                  (-660*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 4)*pow(xj, 10)*

                  (-264*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) -

                  165*pow(xi, 6)*pow(xj, 8)*

                  (-120*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) -

                  5*r*pow(xi, 13)*pow(xj, 2)*

                  (-6876*r*pow(xj, 2) + 88*pow(r, 3)*pow(xj, 4)) +

                  5*r*pow(xi, 11)*pow(xj, 4)*

                  (-4884*r*pow(xj, 2) + 88*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 8)*pow(xj, 6)*

                  (-7128*r*pow(xj, 2) + 104*pow(r, 3)*pow(xj, 4)) +

                  r*pow(xi, 15)*(-7380*r*pow(xj, 2) + 136*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 10)*pow(xj, 4)*

                  (2728*r*pow(xj, 2) + 264*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 14)*(-5864*r*pow(xj, 2) +

                                        376*pow(r, 3)*pow(xj, 4)) -

                  15*pow(xi, 12)*pow(xj, 2)*

                  (-10504*r*pow(xj, 2) + 616*pow(r, 3)*pow(xj, 4)) +

                  pow(xi, 9)*pow(xj, 6)*

                  (234135 - 4950*pow(r, 2)*pow(xj, 2) - 34*pow(r, 4)*pow(xj, 4))

                  - 5*pow(xi, 7)*pow(xj, 8)*(6237 - 1242*pow(r, 2)*pow(xj, 2) +

                                                       2*pow(r, 4)*pow(xj, 4)) +

                  3*pow(xi, 5)*pow(xj, 10)*

                  (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  5*pow(xi, 13)*pow(xj, 2)*

                  (43875 - 3438*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 5*pow(xi, 11)*pow(xj, 4)*(7695 - 2442*pow(r, 2)*pow(xj, 2) +

                                                        22*pow(r, 4)*pow(xj, 4)) +

                  pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                   34*pow(r, 4)*pow(xj, 4))) +

                 2*exp(2*r*xj)*pow(xj, 9)*

                 (-150*pow(r, 4)*pow(xi, 18) - 6*pow(r, 5)*pow(xi, 19) +

                  135*pow(xj, 14) + 225*r*xi*pow(xj, 14) +

                  10*pow(r, 3)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                  30*pow(r, 2)*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                  45*r*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                  45*pow(xi, 2)*pow(xj, 12)*(-33 + 4*pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 9)*pow(xj, 6)*

                  (234135 - 4950*pow(r, 2)*pow(xj, 2) - 34*pow(r, 4)*pow(xj, 4))

                  - 5*r*pow(xi, 7)*pow(xj, 8)*(6237 - 1242*pow(r, 2)*pow(xj, 2) +

                                                         2*pow(r, 4)*pow(xj, 4)) +

                  3*r*pow(xi, 5)*pow(xj, 10)*

                  (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 4)*pow(xj, 10)*

                  (495 - 132*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  165*pow(xi, 6)*pow(xj, 8)*

                  (135 - 60*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  5*r*pow(xi, 13)*pow(xj, 2)*

                  (43875 - 3438*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 5*r*pow(xi, 11)*pow(xj, 4)*

                  (7695 - 2442*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 8)*pow(xj, 6)*

                  (-33 - 3564*pow(r, 2)*pow(xj, 2) + 26*pow(r, 4)*pow(xj, 4)) +

                  r*pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                     34*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 10)*pow(xj, 4)*

                  (-32277 + 1364*pow(r, 2)*pow(xj, 2) + 66*pow(r, 4)*pow(xj, 4))

                  + 15*pow(xi, 14)*(-3003 - 2932*pow(r, 2)*pow(xj, 2) +

                                          94*pow(r, 4)*pow(xj, 4)) -

                  15*pow(xi, 12)*pow(xj, 2)*

                  (28119 - 5252*pow(r, 2)*pow(xj, 2) + 154*pow(r, 4)*pow(xj, 4)))

                 + exp(2*r*xi)*pow(xi, 8)*(-5*pow(xi, 2)*pow(xj, 12)*

                                                 (-43875*xj - 17592*r*pow(xj, 2) - 2214*pow(r, 2)*pow(xj, 3) -

                                                  24*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) -

                                                 3*pow(xi, 14)*(75*xj + 120*r*pow(xj, 2) +

                                                                      90*pow(r, 2)*pow(xj, 3) + 40*pow(r, 3)*pow(xj, 4) +

                                                                      10*pow(r, 4)*pow(xj, 5)) -

                                                 55*pow(xi, 8)*pow(xj, 6)*

                                                 (-567*xj - 1944*r*pow(xj, 2) - 270*pow(r, 2)*pow(xj, 3) +

                                                  72*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                                                 55*pow(xi, 6)*pow(xj, 8)*

                                                 (-4257*xj - 744*r*pow(xj, 2) + 666*pow(r, 2)*pow(xj, 3) +

                                                  168*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                                                 3*pow(xj, 14)*(10725*xj + 6600*r*pow(xj, 2) +

                                                                      1650*pow(r, 2)*pow(xj, 3) + 200*pow(r, 3)*pow(xj, 4) +

                                                                      10*pow(r, 4)*pow(xj, 5)) +

                                                 5*pow(xi, 12)*pow(xj, 2)*

                                                 (495*xj + 792*r*pow(xj, 2) + 594*pow(r, 2)*pow(xj, 3) +

                                                  264*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                                                 pow(xi, 10)*pow(xj, 4)*

                                                 (-12375*xj - 19800*r*pow(xj, 2) - 18630*pow(r, 2)*pow(xj, 3) -

                                                  1560*pow(r, 3)*pow(xj, 4) + 170*pow(r, 4)*pow(xj, 5)) -

                                                 pow(xi, 4)*pow(xj, 10)*

                                                 (38475*xj + 157560*r*pow(xj, 2) + 51570*pow(r, 2)*pow(xj, 3) +

                                                  5640*pow(r, 3)*pow(xj, 4) + 170*pow(r, 4)*pow(xj, 5))) +

                 2*exp(2*r*xi)*pow(xi, 9)*

                 (-5*pow(xi, 2)*pow(xj, 12)*

                  (-84357 - 43875*r*xj - 8796*pow(r, 2)*pow(xj, 2) -

                   738*pow(r, 3)*pow(xj, 3) - 6*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  3*pow(xi, 14)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                       30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) -

                  55*pow(xi, 8)*pow(xj, 6)*

                  (-405 - 567*r*xj - 972*pow(r, 2)*pow(xj, 2) -

                   90*pow(r, 3)*pow(xj, 3) + 18*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  55*pow(xi, 6)*pow(xj, 8)*

                  (9 - 4257*r*xj - 372*pow(r, 2)*pow(xj, 2) +

                   222*pow(r, 3)*pow(xj, 3) + 42*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  3*pow(xj, 14)*(15015 + 10725*r*xj + 3300*pow(r, 2)*pow(xj, 2) +

                                       550*pow(r, 3)*pow(xj, 3) + 50*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 12)*pow(xj, 2)*

                  (297 + 495*r*xj + 396*pow(r, 2)*pow(xj, 2) +

                   198*pow(r, 3)*pow(xj, 3) + 66*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 10)*pow(xj, 4)*

                  (-7425 - 12375*r*xj - 9900*pow(r, 2)*pow(xj, 2) -

                   6210*pow(r, 3)*pow(xj, 3) - 390*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 4)*pow(xj, 10)*

                  (-484155 + 38475*r*xj + 78780*pow(r, 2)*pow(xj, 2) +

                   17190*pow(r, 3)*pow(xj, 3) + 1410*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5))))/

                (135*exp(2*r*(xi + xj))*r*pow(xi - xj, 11)*pow(xi + xj, 11))

            ;
        }

    }
    return S;
}

#endif
