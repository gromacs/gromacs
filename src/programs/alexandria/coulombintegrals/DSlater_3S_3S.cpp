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
