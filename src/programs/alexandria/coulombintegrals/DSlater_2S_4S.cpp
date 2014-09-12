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

cl_R DSlater_2S_4S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-1125310725LL*xi + 1277337600LL*exp(2LL*r*xi)*xi -

                  1946567700LL*r*Power(xi, 2LL) - 1647191700LL*Power(r, 2LL)*Power(xi, 3LL) -

                  904780800LL*Power(r, 3LL)*Power(xi, 4LL) - 360498600LL*Power(r, 4LL)*Power(xi, 5LL) -

                  110103840LL*Power(r, 5LL)*Power(xi, 6LL) - 26500320LL*Power(r, 6LL)*Power(xi, 7LL) -

                  5068800LL*Power(r, 7LL)*Power(xi, 8LL) - 760320LL*Power(r, 8LL)*Power(xi, 9LL) -

                  84480LL*Power(r, 9LL)*Power(xi, 10LL) - 5632LL*Power(r, 10LL)*Power(xi, 11LL))/

                (6.386688e8*exp(2LL*r*xi)*r) +

                (-638668800LL + 638668800LL*exp(2LL*r*xi) - 1125310725LL*r*xi -

                 973283850LL*Power(r, 2LL)*Power(xi, 2LL) - 549063900LL*Power(r, 3LL)*Power(xi, 3LL) -

                 226195200LL*Power(r, 4LL)*Power(xi, 4LL) - 72099720LL*Power(r, 5LL)*Power(xi, 5LL) -

                 18350640LL*Power(r, 6LL)*Power(xi, 6LL) - 3785760LL*Power(r, 7LL)*Power(xi, 7LL) -

                 633600LL*Power(r, 8LL)*Power(xi, 8LL) - 84480LL*Power(r, 9LL)*Power(xi, 9LL) -

                 8448LL*Power(r, 10LL)*Power(xi, 10LL) - 512LL*Power(r, 11LL)*Power(xi, 11LL))/

                (6.386688e8*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-638668800LL + 638668800LL*exp(2LL*r*xi) - 1125310725LL*r*xi -

                     973283850LL*Power(r, 2LL)*Power(xi, 2LL) - 549063900LL*Power(r, 3LL)*Power(xi, 3LL) -

                     226195200LL*Power(r, 4LL)*Power(xi, 4LL) - 72099720LL*Power(r, 5LL)*Power(xi, 5LL) -

                     18350640LL*Power(r, 6LL)*Power(xi, 6LL) - 3785760LL*Power(r, 7LL)*Power(xi, 7LL) -

                     633600LL*Power(r, 8LL)*Power(xi, 8LL) - 84480LL*Power(r, 9LL)*Power(xi, 9LL) -

                     8448LL*Power(r, 10LL)*Power(xi, 10LL) - 512LL*Power(r, 11LL)*Power(xi, 11LL)))/

                (3.193344e8*exp(2LL*r*xi)*r)

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
            S = (1260LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                 210LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-36LL*Power(r, 2LL)*Power(xi, 14LL) - 2LL*Power(r, 3LL)*Power(xi, 15LL) -

                  1287LL*r*Power(xi, 9LL)*Power(xj, 4LL) + 6LL*Power(xj, 12LL) +

                  9LL*r*xi*Power(xj, 12LL) - 22LL*r*Power(xi, 7LL)*Power(xj, 6LL)*

                  (-135LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  6LL*Power(xi, 2LL)*Power(xj, 10LL)*(-11LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  66LL*Power(xi, 4LL)*Power(xj, 8LL)*(-5LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  8LL*r*Power(xi, 5LL)*Power(xj, 8LL)*(99LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 3LL)*Power(xj, 10LL)*(-99LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  132LL*Power(xi, 6LL)*Power(xj, 6LL)*(27LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  78LL*Power(xi, 12LL)*(7LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  2LL*r*Power(xi, 13LL)*(117LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  66LL*Power(xi, 8LL)*Power(xj, 4LL)*(-191LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 11LL)*Power(xj, 2LL)*(-2151LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  6LL*Power(xi, 10LL)*Power(xj, 2LL)*(-1099LL + 33LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-385LL*Power(xi, 8LL)*Power(xj, 8LL)*

                  (1080LL + 1935LL*r*xj + 1350LL*Power(r, 2LL)*Power(xj, 2LL) +

            1170LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

            66LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  7LL*Power(xi, 6LL)*Power(xj, 10LL)*

                  (99540LL + 58095LL*r*xj + 190710LL*Power(r, 2LL)*Power(xj, 2LL) +

            100950LL*Power(r, 3LL)*Power(xj, 3LL) + 21660LL*Power(r, 4LL)*Power(xj, 4LL) +

            1938LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL) -

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  4LL*Power(xj, 16LL)*(135135LL + 135135LL*r*xj +

                                       62370LL*Power(r, 2LL)*Power(xj, 2LL) + 17325LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       3150LL*Power(r, 4LL)*Power(xj, 4LL) + 378LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       28LL*Power(r, 6LL)*Power(xj, 6LL) + Power(r, 7LL)*Power(xj, 7LL)) -

                  Power(xi, 16LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  7LL*Power(xi, 4LL)*Power(xj, 12LL)*

                  (114660LL - 343395LL*r*xj - 242910LL*Power(r, 2LL)*Power(xj, 2LL) -

            61950LL*Power(r, 3LL)*Power(xj, 3LL) - 6060LL*Power(r, 4LL)*Power(xj, 4LL) +

            282LL*Power(r, 5LL)*Power(xj, 5LL) + 116LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  7LL*Power(xi, 12LL)*Power(xj, 4LL)*

                  (9900LL + 17325LL*r*xj + 14850LL*Power(r, 2LL)*Power(xj, 2LL) +

            8250LL*Power(r, 3LL)*Power(xj, 3LL) + 3300LL*Power(r, 4LL)*Power(xj, 4LL) +

            1074LL*Power(r, 5LL)*Power(xj, 5LL) + 164LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  7LL*Power(xi, 10LL)*Power(xj, 6LL)*

                  (29700LL + 51975LL*r*xj + 44550LL*Power(r, 2LL)*Power(xj, 2LL) +

            23850LL*Power(r, 3LL)*Power(xj, 3LL) + 11700LL*Power(r, 4LL)*Power(xj, 4LL) +

            2814LL*Power(r, 5LL)*Power(xj, 5LL) + 284LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  Power(xi, 14LL)*Power(xj, 2LL)*

                  (13860LL + 24255LL*r*xj + 20790LL*Power(r, 2LL)*Power(xj, 2LL) +

            11550LL*Power(r, 3LL)*Power(xj, 3LL) + 4620LL*Power(r, 4LL)*Power(xj, 4LL) +

            1386LL*Power(r, 5LL)*Power(xj, 5LL) + 308LL*Power(r, 6LL)*Power(xj, 6LL) +

            24LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  Power(xi, 2LL)*Power(xj, 14LL)*

                  (-3063060LL - 1936935LL*r*xj - 408870LL*Power(r, 2LL)*Power(xj, 2LL) +

            11550LL*Power(r, 3LL)*Power(xj, 3LL) + 23100LL*Power(r, 4LL)*Power(xj, 4LL) +

            5082LL*Power(r, 5LL)*Power(xj, 5LL) + 532LL*Power(r, 6LL)*Power(xj, 6LL) +

            24LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1260LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 11LL)*

                 Power(xi + xj, 11LL)) + (1260LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                                          210LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                                          (-36LL*Power(r, 2LL)*Power(xi, 14LL) - 2LL*Power(r, 3LL)*Power(xi, 15LL) -

                                  1287LL*r*Power(xi, 9LL)*Power(xj, 4LL) + 6LL*Power(xj, 12LL) +

                                  9LL*r*xi*Power(xj, 12LL) - 22LL*r*Power(xi, 7LL)*Power(xj, 6LL)*

                                  (-135LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  6LL*Power(xi, 2LL)*Power(xj, 10LL)*(-11LL + Power(r, 2LL)*Power(xj, 2LL)) -

                                  66LL*Power(xi, 4LL)*Power(xj, 8LL)*(-5LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  8LL*r*Power(xi, 5LL)*Power(xj, 8LL)*(99LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  r*Power(xi, 3LL)*Power(xj, 10LL)*(-99LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  132LL*Power(xi, 6LL)*Power(xj, 6LL)*(27LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  78LL*Power(xi, 12LL)*(7LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  2LL*r*Power(xi, 13LL)*(117LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  66LL*Power(xi, 8LL)*Power(xj, 4LL)*(-191LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  r*Power(xi, 11LL)*Power(xj, 2LL)*(-2151LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  6LL*Power(xi, 10LL)*Power(xj, 2LL)*(-1099LL + 33LL*Power(r, 2LL)*Power(xj, 2LL))) +

                                          exp(2LL*r*xi)*Power(xi, 6LL)*

                                          (-385LL*Power(xi, 8LL)*Power(xj, 8LL)*

                                  (1080LL + 1935LL*r*xj + 1350LL*Power(r, 2LL)*Power(xj, 2LL) +

            1170LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

            66LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  7LL*Power(xi, 6LL)*Power(xj, 10LL)*

                                  (99540LL + 58095LL*r*xj + 190710LL*Power(r, 2LL)*Power(xj, 2LL) +

            100950LL*Power(r, 3LL)*Power(xj, 3LL) + 21660LL*Power(r, 4LL)*Power(xj, 4LL) +

            1938LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL) -

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  4LL*Power(xj, 16LL)*(135135LL + 135135LL*r*xj +

                                                       62370LL*Power(r, 2LL)*Power(xj, 2LL) + 17325LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                       3150LL*Power(r, 4LL)*Power(xj, 4LL) + 378LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                       28LL*Power(r, 6LL)*Power(xj, 6LL) + Power(r, 7LL)*Power(xj, 7LL)) -

                                  Power(xi, 16LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  7LL*Power(xi, 4LL)*Power(xj, 12LL)*

                                  (114660LL - 343395LL*r*xj - 242910LL*Power(r, 2LL)*Power(xj, 2LL) -

            61950LL*Power(r, 3LL)*Power(xj, 3LL) - 6060LL*Power(r, 4LL)*Power(xj, 4LL) +

            282LL*Power(r, 5LL)*Power(xj, 5LL) + 116LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  7LL*Power(xi, 12LL)*Power(xj, 4LL)*

                                  (9900LL + 17325LL*r*xj + 14850LL*Power(r, 2LL)*Power(xj, 2LL) +

            8250LL*Power(r, 3LL)*Power(xj, 3LL) + 3300LL*Power(r, 4LL)*Power(xj, 4LL) +

            1074LL*Power(r, 5LL)*Power(xj, 5LL) + 164LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  7LL*Power(xi, 10LL)*Power(xj, 6LL)*

                                  (29700LL + 51975LL*r*xj + 44550LL*Power(r, 2LL)*Power(xj, 2LL) +

            23850LL*Power(r, 3LL)*Power(xj, 3LL) + 11700LL*Power(r, 4LL)*Power(xj, 4LL) +

            2814LL*Power(r, 5LL)*Power(xj, 5LL) + 284LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  Power(xi, 14LL)*Power(xj, 2LL)*

                                  (13860LL + 24255LL*r*xj + 20790LL*Power(r, 2LL)*Power(xj, 2LL) +

            11550LL*Power(r, 3LL)*Power(xj, 3LL) + 4620LL*Power(r, 4LL)*Power(xj, 4LL) +

            1386LL*Power(r, 5LL)*Power(xj, 5LL) + 308LL*Power(r, 6LL)*Power(xj, 6LL) +

            24LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  Power(xi, 2LL)*Power(xj, 14LL)*

                                  (-3063060LL - 1936935LL*r*xj - 408870LL*Power(r, 2LL)*Power(xj, 2LL) +

            11550LL*Power(r, 3LL)*Power(xj, 3LL) + 23100LL*Power(r, 4LL)*Power(xj, 4LL) +

            5082LL*Power(r, 5LL)*Power(xj, 5LL) + 532LL*Power(r, 6LL)*Power(xj, 6LL) +

            24LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (630LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 11LL)*Power(xi + xj, 10LL)) -

                (2520LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 11LL) +

                 210LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-72LL*r*Power(xi, 14LL) - 6LL*Power(r, 2LL)*Power(xi, 15LL) -

        468LL*r*Power(xi, 12LL)*Power(xj, 2LL) -

        16LL*Power(r, 2LL)*Power(xi, 13LL)*Power(xj, 2LL) -

        1287LL*Power(xi, 9LL)*Power(xj, 4LL) + 396LL*r*Power(xi, 10LL)*Power(xj, 4LL) +

        44LL*Power(r, 2LL)*Power(xi, 11LL)*Power(xj, 4LL) +

        792LL*r*Power(xi, 8LL)*Power(xj, 6LL) - 528LL*r*Power(xi, 6LL)*Power(xj, 8LL) -

        44LL*Power(r, 2LL)*Power(xi, 7LL)*Power(xj, 8LL) -

        132LL*r*Power(xi, 4LL)*Power(xj, 10LL) +

        16LL*Power(r, 2LL)*Power(xi, 5LL)*Power(xj, 10LL) + 9LL*xi*Power(xj, 12LL) +

        12LL*r*Power(xi, 2LL)*Power(xj, 12LL) +

        4LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 12LL) -

        22LL*Power(xi, 7LL)*Power(xj, 6LL)*(-135LL + Power(r, 2LL)*Power(xj, 2LL)) +

        8LL*Power(xi, 5LL)*Power(xj, 8LL)*(99LL + Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 3LL)*Power(xj, 10LL)*(-99LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        2LL*Power(xi, 13LL)*(117LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 11LL)*Power(xj, 2LL)*(-2151LL + 22LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 420LL*exp(2LL*r*xj)*Power(xj, 11LL)*

                 (-36LL*Power(r, 2LL)*Power(xi, 14LL) - 2LL*Power(r, 3LL)*Power(xi, 15LL) -

        1287LL*r*Power(xi, 9LL)*Power(xj, 4LL) + 6LL*Power(xj, 12LL) +

        9LL*r*xi*Power(xj, 12LL) - 22LL*r*Power(xi, 7LL)*Power(xj, 6LL)*

        (-135LL + Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 2LL)*Power(xj, 10LL)*(-11LL + Power(r, 2LL)*Power(xj, 2LL)) -

        66LL*Power(xi, 4LL)*Power(xj, 8LL)*(-5LL + Power(r, 2LL)*Power(xj, 2LL)) +

        8LL*r*Power(xi, 5LL)*Power(xj, 8LL)*(99LL + Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 3LL)*Power(xj, 10LL)*(-99LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        132LL*Power(xi, 6LL)*Power(xj, 6LL)*(27LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        78LL*Power(xi, 12LL)*(7LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) -

        2LL*r*Power(xi, 13LL)*(117LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) +

        66LL*Power(xi, 8LL)*Power(xj, 4LL)*(-191LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 11LL)*Power(xj, 2LL)*(-2151LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 10LL)*Power(xj, 2LL)*(-1099LL + 33LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-385LL*Power(xi, 8LL)*Power(xj, 8LL)*

        (1935LL*xj + 2700LL*r*Power(xj, 2LL) + 3510LL*Power(r, 2LL)*Power(xj, 3LL) +

            1680LL*Power(r, 3LL)*Power(xj, 4LL) + 330LL*Power(r, 4LL)*Power(xj, 5LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) +

        7LL*Power(xi, 6LL)*Power(xj, 10LL)*

        (58095LL*xj + 381420LL*r*Power(xj, 2LL) +

            302850LL*Power(r, 2LL)*Power(xj, 3LL) + 86640LL*Power(r, 3LL)*Power(xj, 4LL) +

            9690LL*Power(r, 4LL)*Power(xj, 5LL) + 24LL*Power(r, 5LL)*Power(xj, 6LL) -

            56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        4LL*Power(xj, 16LL)*(135135LL*xj + 124740LL*r*Power(xj, 2LL) +

                             51975LL*Power(r, 2LL)*Power(xj, 3LL) + 12600LL*Power(r, 3LL)*Power(xj, 4LL) +

                             1890LL*Power(r, 4LL)*Power(xj, 5LL) + 168LL*Power(r, 5LL)*Power(xj, 6LL) +

                             7LL*Power(r, 6LL)*Power(xj, 7LL)) -

        Power(xi, 16LL)*(2205LL*xj + 3780LL*r*Power(xj, 2LL) +

                         3150LL*Power(r, 2LL)*Power(xj, 3LL) + 1680LL*Power(r, 3LL)*Power(xj, 4LL) +

                         630LL*Power(r, 4LL)*Power(xj, 5LL) + 168LL*Power(r, 5LL)*Power(xj, 6LL) +

                         28LL*Power(r, 6LL)*Power(xj, 7LL)) +

        7LL*Power(xi, 4LL)*Power(xj, 12LL)*

        (-343395LL*xj - 485820LL*r*Power(xj, 2LL) -

            185850LL*Power(r, 2LL)*Power(xj, 3LL) - 24240LL*Power(r, 3LL)*Power(xj, 4LL) +

            1410LL*Power(r, 4LL)*Power(xj, 5LL) + 696LL*Power(r, 5LL)*Power(xj, 6LL) +

            56LL*Power(r, 6LL)*Power(xj, 7LL)) -

        7LL*Power(xi, 12LL)*Power(xj, 4LL)*

        (17325LL*xj + 29700LL*r*Power(xj, 2LL) + 24750LL*Power(r, 2LL)*Power(xj, 3LL) +

            13200LL*Power(r, 3LL)*Power(xj, 4LL) + 5370LL*Power(r, 4LL)*Power(xj, 5LL) +

            984LL*Power(r, 5LL)*Power(xj, 6LL) + 56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        7LL*Power(xi, 10LL)*Power(xj, 6LL)*

        (51975LL*xj + 89100LL*r*Power(xj, 2LL) + 71550LL*Power(r, 2LL)*Power(xj, 3LL) +

            46800LL*Power(r, 3LL)*Power(xj, 4LL) + 14070LL*Power(r, 4LL)*Power(xj, 5LL) +

            1704LL*Power(r, 5LL)*Power(xj, 6LL) + 56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        Power(xi, 14LL)*Power(xj, 2LL)*

        (24255LL*xj + 41580LL*r*Power(xj, 2LL) + 34650LL*Power(r, 2LL)*Power(xj, 3LL) +

            18480LL*Power(r, 3LL)*Power(xj, 4LL) + 6930LL*Power(r, 4LL)*Power(xj, 5LL) +

            1848LL*Power(r, 5LL)*Power(xj, 6LL) + 168LL*Power(r, 6LL)*Power(xj, 7LL)) -

        Power(xi, 2LL)*Power(xj, 14LL)*

        (-1936935LL*xj - 817740LL*r*Power(xj, 2LL) +

            34650LL*Power(r, 2LL)*Power(xj, 3LL) + 92400LL*Power(r, 3LL)*Power(xj, 4LL) +

            25410LL*Power(r, 4LL)*Power(xj, 5LL) + 3192LL*Power(r, 5LL)*Power(xj, 6LL) +

            168LL*Power(r, 6LL)*Power(xj, 7LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 7LL)*

                 (-385LL*Power(xi, 8LL)*Power(xj, 8LL)*

        (1080LL + 1935LL*r*xj + 1350LL*Power(r, 2LL)*Power(xj, 2LL) +

            1170LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

            66LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

        7LL*Power(xi, 6LL)*Power(xj, 10LL)*

        (99540LL + 58095LL*r*xj + 190710LL*Power(r, 2LL)*Power(xj, 2LL) +

            100950LL*Power(r, 3LL)*Power(xj, 3LL) + 21660LL*Power(r, 4LL)*Power(xj, 4LL) +

            1938LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL) -

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        4LL*Power(xj, 16LL)*(135135LL + 135135LL*r*xj +

                             62370LL*Power(r, 2LL)*Power(xj, 2LL) + 17325LL*Power(r, 3LL)*Power(xj, 3LL) +

                             3150LL*Power(r, 4LL)*Power(xj, 4LL) + 378LL*Power(r, 5LL)*Power(xj, 5LL) +

                             28LL*Power(r, 6LL)*Power(xj, 6LL) + Power(r, 7LL)*Power(xj, 7LL)) -

        Power(xi, 16LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                         1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                         126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                         4LL*Power(r, 7LL)*Power(xj, 7LL)) +

        7LL*Power(xi, 4LL)*Power(xj, 12LL)*

        (114660LL - 343395LL*r*xj - 242910LL*Power(r, 2LL)*Power(xj, 2LL) -

            61950LL*Power(r, 3LL)*Power(xj, 3LL) - 6060LL*Power(r, 4LL)*Power(xj, 4LL) +

            282LL*Power(r, 5LL)*Power(xj, 5LL) + 116LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

        7LL*Power(xi, 12LL)*Power(xj, 4LL)*

        (9900LL + 17325LL*r*xj + 14850LL*Power(r, 2LL)*Power(xj, 2LL) +

            8250LL*Power(r, 3LL)*Power(xj, 3LL) + 3300LL*Power(r, 4LL)*Power(xj, 4LL) +

            1074LL*Power(r, 5LL)*Power(xj, 5LL) + 164LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        7LL*Power(xi, 10LL)*Power(xj, 6LL)*

        (29700LL + 51975LL*r*xj + 44550LL*Power(r, 2LL)*Power(xj, 2LL) +

            23850LL*Power(r, 3LL)*Power(xj, 3LL) + 11700LL*Power(r, 4LL)*Power(xj, 4LL) +

            2814LL*Power(r, 5LL)*Power(xj, 5LL) + 284LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        Power(xi, 14LL)*Power(xj, 2LL)*

        (13860LL + 24255LL*r*xj + 20790LL*Power(r, 2LL)*Power(xj, 2LL) +

            11550LL*Power(r, 3LL)*Power(xj, 3LL) + 4620LL*Power(r, 4LL)*Power(xj, 4LL) +

            1386LL*Power(r, 5LL)*Power(xj, 5LL) + 308LL*Power(r, 6LL)*Power(xj, 6LL) +

            24LL*Power(r, 7LL)*Power(xj, 7LL)) -

        Power(xi, 2LL)*Power(xj, 14LL)*

        (-3063060LL - 1936935LL*r*xj - 408870LL*Power(r, 2LL)*Power(xj, 2LL) +

            11550LL*Power(r, 3LL)*Power(xj, 3LL) + 23100LL*Power(r, 4LL)*Power(xj, 4LL) +

            5082LL*Power(r, 5LL)*Power(xj, 5LL) + 532LL*Power(r, 6LL)*Power(xj, 6LL) +

            24LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1260LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 11LL)*Power(xi + xj, 11LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_4S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_2S_4S(r, xj, xi);
}
