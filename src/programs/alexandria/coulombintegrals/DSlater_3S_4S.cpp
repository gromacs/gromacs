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

cl_R DSlater_3S_4S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-132871488750LL*xi + 149448499200LL*exp(2LL*r*xi)*xi -

                  232588956600LL*r*Power(xi, 2LL) - 200036962125LL*Power(r, 2LL)*Power(xi, 3LL) -

                  112459347000LL*Power(r, 3LL)*Power(xi, 4LL) -

                  46370223900LL*Power(r, 4LL)*Power(xi, 5LL) -

                  14905931040LL*Power(r, 5LL)*Power(xi, 6LL) -

                  3872428560LL*Power(r, 6LL)*Power(xi, 7LL) -

                  830269440LL*Power(r, 7LL)*Power(xi, 8LL) - 148262400LL*Power(r, 8LL)*Power(xi, 9LL) -

                  21964800LL*Power(r, 9LL)*Power(xi, 10LL) - 2635776LL*Power(r, 10LL)*Power(xi, 11LL) -

                  239616LL*Power(r, 11LL)*Power(xi, 12LL) - 13312LL*Power(r, 12LL)*Power(xi, 13LL))/

                (7.47242496e10*exp(2LL*r*xi)*r) +

                (-74724249600LL + 74724249600LL*exp(2LL*r*xi) - 132871488750LL*r*xi -

                 116294478300LL*Power(r, 2LL)*Power(xi, 2LL) -

                 66678987375LL*Power(r, 3LL)*Power(xi, 3LL) -

                 28114836750LL*Power(r, 4LL)*Power(xi, 4LL) -

                 9274044780LL*Power(r, 5LL)*Power(xi, 5LL) -

                 2484321840LL*Power(r, 6LL)*Power(xi, 6LL) - 553204080LL*Power(r, 7LL)*Power(xi, 7LL) -

                 103783680LL*Power(r, 8LL)*Power(xi, 8LL) - 16473600LL*Power(r, 9LL)*Power(xi, 9LL) -

                 2196480LL*Power(r, 10LL)*Power(xi, 10LL) - 239616LL*Power(r, 11LL)*Power(xi, 11LL) -

                 19968LL*Power(r, 12LL)*Power(xi, 12LL) - 1024LL*Power(r, 13LL)*Power(xi, 13LL))/

                (7.47242496e10*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-74724249600LL + 74724249600LL*exp(2LL*r*xi) - 132871488750LL*r*xi -

                     116294478300LL*Power(r, 2LL)*Power(xi, 2LL) -

                     66678987375LL*Power(r, 3LL)*Power(xi, 3LL) -

                     28114836750LL*Power(r, 4LL)*Power(xi, 4LL) -

                     9274044780LL*Power(r, 5LL)*Power(xi, 5LL) -

                     2484321840LL*Power(r, 6LL)*Power(xi, 6LL) -

                     553204080LL*Power(r, 7LL)*Power(xi, 7LL) - 103783680LL*Power(r, 8LL)*Power(xi, 8LL) -

                     16473600LL*Power(r, 9LL)*Power(xi, 9LL) - 2196480LL*Power(r, 10LL)*Power(xi, 10LL) -

                     239616LL*Power(r, 11LL)*Power(xi, 11LL) - 19968LL*Power(r, 12LL)*Power(xi, 12LL) -

                     1024LL*Power(r, 13LL)*Power(xi, 13LL)))/(3.73621248e10*exp(2LL*r*xi)*r)

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
            S = (3780LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                 84LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-60LL*Power(r, 4LL)*Power(xi, 20LL) - 2LL*Power(r, 5LL)*Power(xi, 21LL) +

                  45LL*Power(xj, 16LL) + 75LL*r*xi*Power(xj, 16LL) -

                  4LL*Power(r, 3LL)*Power(xi, 19LL)*(195LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  15LL*r*Power(xi, 3LL)*Power(xj, 14LL)*(-65LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  15LL*Power(xi, 2LL)*Power(xj, 14LL)*(-39LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  30LL*Power(r, 2LL)*Power(xi, 18LL)*(182LL + 9LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  30LL*r*Power(xi, 13LL)*Power(xj, 4LL)*

                  (-13047LL + 377LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  2LL*r*Power(xi, 5LL)*Power(xj, 12LL)*

                  (2925LL - 195LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  10LL*Power(xi, 4LL)*Power(xj, 12LL)*

                  (351LL - 78LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) -

                  130LL*Power(xi, 6LL)*Power(xj, 10LL)*

                  (99LL - 36LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  13LL*r*Power(xi, 11LL)*Power(xj, 6LL)*

                  (30735LL - 1650LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL))

                  + r*Power(xi, 7LL)*Power(xj, 10LL)*(-15015LL + 3330LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                      4LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  210LL*Power(xi, 16LL)*(-156LL - 262LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         5LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  6LL*r*Power(xi, 9LL)*Power(xj, 8LL)*

                  (-48620LL - 715LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 3LL*r*Power(xi, 17LL)*(-6825LL - 1870LL*Power(r, 2LL)*Power(xj, 2LL) +

                                           12LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  30LL*Power(xi, 14LL)*Power(xj, 2LL)*

                  (17934LL - 12LL*Power(r, 2LL)*Power(xj, 2LL) + 13LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*Power(xi, 8LL)*Power(xj, 8LL)*

                  (2145LL + 2860LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 65LL*Power(xi, 10LL)*Power(xj, 6LL)*

                  (-13725LL - 792LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

                  - 10LL*Power(xi, 12LL)*Power(xj, 4LL)*

                  (153630LL - 15054LL*Power(r, 2LL)*Power(xj, 2LL) +

            143LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  Power(xi, 15LL)*(-269325LL*r*Power(xj, 2LL) +

                                   9270LL*Power(r, 3LL)*Power(xj, 4LL) - 52LL*Power(r, 5LL)*Power(xj, 6LL))) +

                 exp(2LL*r*xi)*Power(xi, 8LL)*

                 (Power(xi, 2LL)*Power(xj, 16LL)*

                  (70073640LL + 47669895LL*r*xj + 13931190LL*Power(r, 2LL)*Power(xj, 2LL) +

            2170350LL*Power(r, 3LL)*Power(xj, 3LL) +

            169260LL*Power(r, 4LL)*Power(xj, 4LL) + 1638LL*Power(r, 5LL)*Power(xj, 5LL) -

            756LL*Power(r, 6LL)*Power(xj, 6LL) - 44LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  364LL*Power(xi, 10LL)*Power(xj, 8LL)*

                  (-7425LL - 13860LL*r*xj - 5940LL*Power(r, 2LL)*Power(xj, 2LL) -

            11880LL*Power(r, 3LL)*Power(xj, 3LL) - 2640LL*Power(r, 4LL)*Power(xj, 4LL) -

            45LL*Power(r, 5LL)*Power(xj, 5LL) + 30LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  364LL*Power(xi, 8LL)*Power(xj, 10LL)*

                  (-20925LL + 18270LL*r*xj - 58320LL*Power(r, 2LL)*Power(xj, 2LL) -

            17730LL*Power(r, 3LL)*Power(xj, 3LL) - 300LL*Power(r, 4LL)*Power(xj, 4LL) +

            423LL*Power(r, 5LL)*Power(xj, 5LL) + 54LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  3LL*Power(xi, 18LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  3LL*Power(xj, 18LL)*(1801800LL + 1576575LL*r*xj +

                                       630630LL*Power(r, 2LL)*Power(xj, 2LL) + 150150LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       23100LL*Power(r, 4LL)*Power(xj, 4LL) + 2310LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       140LL*Power(r, 6LL)*Power(xj, 6LL) + 4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  2LL*Power(xi, 14LL)*Power(xj, 4LL)*

                  (-147420LL - 257985LL*r*xj - 221130LL*Power(r, 2LL)*Power(xj, 2LL) -

            122850LL*Power(r, 3LL)*Power(xj, 3LL) - 49140LL*Power(r, 4LL)*Power(xj, 4LL) -

            17388LL*Power(r, 5LL)*Power(xj, 5LL) - 1512LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  42LL*Power(xi, 12LL)*Power(xj, 6LL)*

                  (-25740LL - 45045LL*r*xj - 38610LL*Power(r, 2LL)*Power(xj, 2LL) -

            19470LL*Power(r, 3LL)*Power(xj, 3LL) - 12540LL*Power(r, 4LL)*Power(xj, 4LL) -

            1836LL*Power(r, 5LL)*Power(xj, 5LL) - 8LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  42LL*Power(xi, 6LL)*Power(xj, 12LL)*

                  (921600LL - 1640835LL*r*xj - 546030LL*Power(r, 2LL)*Power(xj, 2LL) +

            20730LL*Power(r, 3LL)*Power(xj, 3LL) + 30180LL*Power(r, 4LL)*Power(xj, 4LL) +

            5028LL*Power(r, 5LL)*Power(xj, 5LL) + 344LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  2LL*Power(xi, 4LL)*Power(xj, 14LL)*

                  (-67767840LL - 13377735LL*r*xj + 6601770LL*Power(r, 2LL)*Power(xj, 2LL) +

            3115350LL*Power(r, 3LL)*Power(xj, 3LL) +

            548940LL*Power(r, 4LL)*Power(xj, 4LL) + 48132LL*Power(r, 5LL)*Power(xj, 5LL) +

            1848LL*Power(r, 6LL)*Power(xj, 6LL) + 8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  Power(xi, 16LL)*Power(xj, 2LL)*

                  (49140LL + 85995LL*r*xj + 73710LL*Power(r, 2LL)*Power(xj, 2LL) +

            40950LL*Power(r, 3LL)*Power(xj, 3LL) + 16380LL*Power(r, 4LL)*Power(xj, 4LL) +

            4914LL*Power(r, 5LL)*Power(xj, 5LL) + 1092LL*Power(r, 6LL)*Power(xj, 6LL) +

            44LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (3780LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 13LL)*

                 Power(xi + xj, 13LL)) + (3780LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                                          84LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                                          (-60LL*Power(r, 4LL)*Power(xi, 20LL) - 2LL*Power(r, 5LL)*Power(xi, 21LL) +

                                  45LL*Power(xj, 16LL) + 75LL*r*xi*Power(xj, 16LL) -

                                  4LL*Power(r, 3LL)*Power(xi, 19LL)*(195LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  15LL*r*Power(xi, 3LL)*Power(xj, 14LL)*(-65LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  15LL*Power(xi, 2LL)*Power(xj, 14LL)*(-39LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  30LL*Power(r, 2LL)*Power(xi, 18LL)*(182LL + 9LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  30LL*r*Power(xi, 13LL)*Power(xj, 4LL)*

                                  (-13047LL + 377LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  2LL*r*Power(xi, 5LL)*Power(xj, 12LL)*

                                  (2925LL - 195LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  10LL*Power(xi, 4LL)*Power(xj, 12LL)*

                                  (351LL - 78LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) -

                                  130LL*Power(xi, 6LL)*Power(xj, 10LL)*

                                  (99LL - 36LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  13LL*r*Power(xi, 11LL)*Power(xj, 6LL)*

                                  (30735LL - 1650LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL))

                                  + r*Power(xi, 7LL)*Power(xj, 10LL)*(-15015LL + 3330LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                      4LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  210LL*Power(xi, 16LL)*(-156LL - 262LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         5LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  6LL*r*Power(xi, 9LL)*Power(xj, 8LL)*

                                  (-48620LL - 715LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL))

                                  + 3LL*r*Power(xi, 17LL)*(-6825LL - 1870LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                           12LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  30LL*Power(xi, 14LL)*Power(xj, 2LL)*

                                  (17934LL - 12LL*Power(r, 2LL)*Power(xj, 2LL) + 13LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  15LL*Power(xi, 8LL)*Power(xj, 8LL)*

                                  (2145LL + 2860LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

                                  + 65LL*Power(xi, 10LL)*Power(xj, 6LL)*

                                  (-13725LL - 792LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

                                  - 10LL*Power(xi, 12LL)*Power(xj, 4LL)*

                                  (153630LL - 15054LL*Power(r, 2LL)*Power(xj, 2LL) +

            143LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  Power(xi, 15LL)*(-269325LL*r*Power(xj, 2LL) +

                                                   9270LL*Power(r, 3LL)*Power(xj, 4LL) - 52LL*Power(r, 5LL)*Power(xj, 6LL))) +

                                          exp(2LL*r*xi)*Power(xi, 8LL)*

                                          (Power(xi, 2LL)*Power(xj, 16LL)*

                                  (70073640LL + 47669895LL*r*xj + 13931190LL*Power(r, 2LL)*Power(xj, 2LL) +

            2170350LL*Power(r, 3LL)*Power(xj, 3LL) +

            169260LL*Power(r, 4LL)*Power(xj, 4LL) + 1638LL*Power(r, 5LL)*Power(xj, 5LL) -

            756LL*Power(r, 6LL)*Power(xj, 6LL) - 44LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  364LL*Power(xi, 10LL)*Power(xj, 8LL)*

                                  (-7425LL - 13860LL*r*xj - 5940LL*Power(r, 2LL)*Power(xj, 2LL) -

            11880LL*Power(r, 3LL)*Power(xj, 3LL) - 2640LL*Power(r, 4LL)*Power(xj, 4LL) -

            45LL*Power(r, 5LL)*Power(xj, 5LL) + 30LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  364LL*Power(xi, 8LL)*Power(xj, 10LL)*

                                  (-20925LL + 18270LL*r*xj - 58320LL*Power(r, 2LL)*Power(xj, 2LL) -

            17730LL*Power(r, 3LL)*Power(xj, 3LL) - 300LL*Power(r, 4LL)*Power(xj, 4LL) +

            423LL*Power(r, 5LL)*Power(xj, 5LL) + 54LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  3LL*Power(xi, 18LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                       126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                       4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  3LL*Power(xj, 18LL)*(1801800LL + 1576575LL*r*xj +

                                                       630630LL*Power(r, 2LL)*Power(xj, 2LL) + 150150LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                       23100LL*Power(r, 4LL)*Power(xj, 4LL) + 2310LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                       140LL*Power(r, 6LL)*Power(xj, 6LL) + 4LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  2LL*Power(xi, 14LL)*Power(xj, 4LL)*

                                  (-147420LL - 257985LL*r*xj - 221130LL*Power(r, 2LL)*Power(xj, 2LL) -

            122850LL*Power(r, 3LL)*Power(xj, 3LL) - 49140LL*Power(r, 4LL)*Power(xj, 4LL) -

            17388LL*Power(r, 5LL)*Power(xj, 5LL) - 1512LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  42LL*Power(xi, 12LL)*Power(xj, 6LL)*

                                  (-25740LL - 45045LL*r*xj - 38610LL*Power(r, 2LL)*Power(xj, 2LL) -

            19470LL*Power(r, 3LL)*Power(xj, 3LL) - 12540LL*Power(r, 4LL)*Power(xj, 4LL) -

            1836LL*Power(r, 5LL)*Power(xj, 5LL) - 8LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  42LL*Power(xi, 6LL)*Power(xj, 12LL)*

                                  (921600LL - 1640835LL*r*xj - 546030LL*Power(r, 2LL)*Power(xj, 2LL) +

            20730LL*Power(r, 3LL)*Power(xj, 3LL) + 30180LL*Power(r, 4LL)*Power(xj, 4LL) +

            5028LL*Power(r, 5LL)*Power(xj, 5LL) + 344LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  2LL*Power(xi, 4LL)*Power(xj, 14LL)*

                                  (-67767840LL - 13377735LL*r*xj + 6601770LL*Power(r, 2LL)*Power(xj, 2LL) +

            3115350LL*Power(r, 3LL)*Power(xj, 3LL) +

            548940LL*Power(r, 4LL)*Power(xj, 4LL) + 48132LL*Power(r, 5LL)*Power(xj, 5LL) +

            1848LL*Power(r, 6LL)*Power(xj, 6LL) + 8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  Power(xi, 16LL)*Power(xj, 2LL)*

                                  (49140LL + 85995LL*r*xj + 73710LL*Power(r, 2LL)*Power(xj, 2LL) +

            40950LL*Power(r, 3LL)*Power(xj, 3LL) + 16380LL*Power(r, 4LL)*Power(xj, 4LL) +

            4914LL*Power(r, 5LL)*Power(xj, 5LL) + 1092LL*Power(r, 6LL)*Power(xj, 6LL) +

            44LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1890LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 13LL)*Power(xi + xj, 12LL)) -

                (7560LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                 84LL*exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-240LL*Power(r, 3LL)*Power(xi, 20LL) - 10LL*Power(r, 4LL)*Power(xi, 21LL) -

        540LL*Power(r, 3LL)*Power(xi, 18LL)*Power(xj, 2LL) -

        8LL*Power(r, 4LL)*Power(xi, 19LL)*Power(xj, 2LL) +

        22620LL*Power(r, 2LL)*Power(xi, 13LL)*Power(xj, 6LL) + 75LL*xi*Power(xj, 16LL) +

        120LL*r*Power(xi, 2LL)*Power(xj, 16LL) +

        60LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 16LL) -

        12LL*Power(r, 2LL)*Power(xi, 19LL)*(195LL + Power(r, 2LL)*Power(xj, 2LL)) +

        15LL*Power(xi, 3LL)*Power(xj, 14LL)*(-65LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        60LL*r*Power(xi, 18LL)*(182LL + 9LL*Power(r, 2LL)*Power(xj, 2LL)) +

        30LL*Power(xi, 13LL)*Power(xj, 4LL)*(-13047LL + 377LL*Power(r, 2LL)*Power(xj, 2LL)) +

        2LL*r*Power(xi, 5LL)*Power(xj, 12LL)*

        (-390LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

        10LL*Power(xi, 4LL)*Power(xj, 12LL)*

        (-156LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) -

        130LL*Power(xi, 6LL)*Power(xj, 10LL)*

        (-72LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

        13LL*r*Power(xi, 11LL)*Power(xj, 6LL)*

        (-3300LL*r*Power(xj, 2LL) + 16LL*Power(r, 3LL)*Power(xj, 4LL)) +

        r*Power(xi, 7LL)*Power(xj, 10LL)*

        (6660LL*r*Power(xj, 2LL) + 16LL*Power(r, 3LL)*Power(xj, 4LL)) +

        210LL*Power(xi, 16LL)*(-524LL*r*Power(xj, 2LL) + 20LL*Power(r, 3LL)*Power(xj, 4LL)) -

        6LL*r*Power(xi, 9LL)*Power(xj, 8LL)*

        (-1430LL*r*Power(xj, 2LL) + 24LL*Power(r, 3LL)*Power(xj, 4LL)) +

        3LL*r*Power(xi, 17LL)*(-3740LL*r*Power(xj, 2LL) +

                               48LL*Power(r, 3LL)*Power(xj, 4LL)) -

        30LL*Power(xi, 14LL)*Power(xj, 2LL)*

        (-24LL*r*Power(xj, 2LL) + 52LL*Power(r, 3LL)*Power(xj, 4LL)) -

        15LL*Power(xi, 8LL)*Power(xj, 8LL)*

        (5720LL*r*Power(xj, 2LL) + 56LL*Power(r, 3LL)*Power(xj, 4LL)) +

        65LL*Power(xi, 10LL)*Power(xj, 6LL)*

        (-1584LL*r*Power(xj, 2LL) + 88LL*Power(r, 3LL)*Power(xj, 4LL)) -

        10LL*Power(xi, 12LL)*Power(xj, 4LL)*

        (-30108LL*r*Power(xj, 2LL) + 572LL*Power(r, 3LL)*Power(xj, 4LL)) +

        2LL*Power(xi, 5LL)*Power(xj, 12LL)*

        (2925LL - 195LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

        13LL*Power(xi, 11LL)*Power(xj, 6LL)*

        (30735LL - 1650LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) +

        Power(xi, 7LL)*Power(xj, 10LL)*

        (-15015LL + 3330LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL))

        - 6LL*Power(xi, 9LL)*Power(xj, 8LL)*(-48620LL - 715LL*Power(r, 2LL)*Power(xj, 2LL) +

                                             6LL*Power(r, 4LL)*Power(xj, 4LL)) +

        3LL*Power(xi, 17LL)*(-6825LL - 1870LL*Power(r, 2LL)*Power(xj, 2LL) +

                             12LL*Power(r, 4LL)*Power(xj, 4LL)) +

        Power(xi, 15LL)*(-269325LL*Power(xj, 2LL) + 27810LL*Power(r, 2LL)*Power(xj, 4LL) -

                         260LL*Power(r, 4LL)*Power(xj, 6LL))) +

                 168LL*exp(2LL*r*xj)*Power(xj, 11LL)*

                 (-60LL*Power(r, 4LL)*Power(xi, 20LL) - 2LL*Power(r, 5LL)*Power(xi, 21LL) +

        45LL*Power(xj, 16LL) + 75LL*r*xi*Power(xj, 16LL) -

        4LL*Power(r, 3LL)*Power(xi, 19LL)*(195LL + Power(r, 2LL)*Power(xj, 2LL)) +

        15LL*r*Power(xi, 3LL)*Power(xj, 14LL)*(-65LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        15LL*Power(xi, 2LL)*Power(xj, 14LL)*(-39LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

        30LL*Power(r, 2LL)*Power(xi, 18LL)*(182LL + 9LL*Power(r, 2LL)*Power(xj, 2LL)) +

        30LL*r*Power(xi, 13LL)*Power(xj, 4LL)*

        (-13047LL + 377LL*Power(r, 2LL)*Power(xj, 2LL)) +

        2LL*r*Power(xi, 5LL)*Power(xj, 12LL)*

        (2925LL - 195LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

        10LL*Power(xi, 4LL)*Power(xj, 12LL)*

        (351LL - 78LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) -

        130LL*Power(xi, 6LL)*Power(xj, 10LL)*

        (99LL - 36LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

        13LL*r*Power(xi, 11LL)*Power(xj, 6LL)*

        (30735LL - 1650LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) +

        r*Power(xi, 7LL)*Power(xj, 10LL)*

        (-15015LL + 3330LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL))

        + 210LL*Power(xi, 16LL)*(-156LL - 262LL*Power(r, 2LL)*Power(xj, 2LL) +

                                 5LL*Power(r, 4LL)*Power(xj, 4LL)) -

        6LL*r*Power(xi, 9LL)*Power(xj, 8LL)*

        (-48620LL - 715LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL)) +

        3LL*r*Power(xi, 17LL)*(-6825LL - 1870LL*Power(r, 2LL)*Power(xj, 2LL) +

                               12LL*Power(r, 4LL)*Power(xj, 4LL)) -

        30LL*Power(xi, 14LL)*Power(xj, 2LL)*

        (17934LL - 12LL*Power(r, 2LL)*Power(xj, 2LL) + 13LL*Power(r, 4LL)*Power(xj, 4LL)) -

        15LL*Power(xi, 8LL)*Power(xj, 8LL)*

        (2145LL + 2860LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL)) +

        65LL*Power(xi, 10LL)*Power(xj, 6LL)*

        (-13725LL - 792LL*Power(r, 2LL)*Power(xj, 2LL) + 22LL*Power(r, 4LL)*Power(xj, 4LL))

        - 10LL*Power(xi, 12LL)*Power(xj, 4LL)*(153630LL - 15054LL*Power(r, 2LL)*Power(xj, 2LL) +

                                               143LL*Power(r, 4LL)*Power(xj, 4LL)) +

        Power(xi, 15LL)*(-269325LL*r*Power(xj, 2LL) + 9270LL*Power(r, 3LL)*Power(xj, 4LL) -

                         52LL*Power(r, 5LL)*Power(xj, 6LL))) +

                 exp(2LL*r*xi)*Power(xi, 8LL)*

                 (Power(xi, 2LL)*Power(xj, 16LL)*

        (47669895LL*xj + 27862380LL*r*Power(xj, 2LL) +

            6511050LL*Power(r, 2LL)*Power(xj, 3LL) +

            677040LL*Power(r, 3LL)*Power(xj, 4LL) + 8190LL*Power(r, 4LL)*Power(xj, 5LL) -

            4536LL*Power(r, 5LL)*Power(xj, 6LL) - 308LL*Power(r, 6LL)*Power(xj, 7LL)) +

        364LL*Power(xi, 10LL)*Power(xj, 8LL)*

        (-13860LL*xj - 11880LL*r*Power(xj, 2LL) - 35640LL*Power(r, 2LL)*Power(xj, 3LL) -

            10560LL*Power(r, 3LL)*Power(xj, 4LL) - 225LL*Power(r, 4LL)*Power(xj, 5LL) +

            180LL*Power(r, 5LL)*Power(xj, 6LL) + 14LL*Power(r, 6LL)*Power(xj, 7LL)) -

        364LL*Power(xi, 8LL)*Power(xj, 10LL)*

        (18270LL*xj - 116640LL*r*Power(xj, 2LL) - 53190LL*Power(r, 2LL)*Power(xj, 3LL) -

            1200LL*Power(r, 3LL)*Power(xj, 4LL) + 2115LL*Power(r, 4LL)*Power(xj, 5LL) +

            324LL*Power(r, 5LL)*Power(xj, 6LL) + 14LL*Power(r, 6LL)*Power(xj, 7LL)) -

        3LL*Power(xi, 18LL)*(2205LL*xj + 3780LL*r*Power(xj, 2LL) +

                             3150LL*Power(r, 2LL)*Power(xj, 3LL) + 1680LL*Power(r, 3LL)*Power(xj, 4LL) +

                             630LL*Power(r, 4LL)*Power(xj, 5LL) + 168LL*Power(r, 5LL)*Power(xj, 6LL) +

                             28LL*Power(r, 6LL)*Power(xj, 7LL)) +

        3LL*Power(xj, 18LL)*(1576575LL*xj + 1261260LL*r*Power(xj, 2LL) +

                             450450LL*Power(r, 2LL)*Power(xj, 3LL) + 92400LL*Power(r, 3LL)*Power(xj, 4LL) +

                             11550LL*Power(r, 4LL)*Power(xj, 5LL) + 840LL*Power(r, 5LL)*Power(xj, 6LL) +

                             28LL*Power(r, 6LL)*Power(xj, 7LL)) +

        2LL*Power(xi, 14LL)*Power(xj, 4LL)*

        (-257985LL*xj - 442260LL*r*Power(xj, 2LL) -

            368550LL*Power(r, 2LL)*Power(xj, 3LL) - 196560LL*Power(r, 3LL)*Power(xj, 4LL) -

            86940LL*Power(r, 4LL)*Power(xj, 5LL) - 9072LL*Power(r, 5LL)*Power(xj, 6LL) +

            56LL*Power(r, 6LL)*Power(xj, 7LL)) -

        42LL*Power(xi, 12LL)*Power(xj, 6LL)*

        (-45045LL*xj - 77220LL*r*Power(xj, 2LL) - 58410LL*Power(r, 2LL)*Power(xj, 3LL) -

            50160LL*Power(r, 3LL)*Power(xj, 4LL) - 9180LL*Power(r, 4LL)*Power(xj, 5LL) -

            48LL*Power(r, 5LL)*Power(xj, 6LL) + 56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        42LL*Power(xi, 6LL)*Power(xj, 12LL)*

        (-1640835LL*xj - 1092060LL*r*Power(xj, 2LL) +

            62190LL*Power(r, 2LL)*Power(xj, 3LL) + 120720LL*Power(r, 3LL)*Power(xj, 4LL) +

            25140LL*Power(r, 4LL)*Power(xj, 5LL) + 2064LL*Power(r, 5LL)*Power(xj, 6LL) +

            56LL*Power(r, 6LL)*Power(xj, 7LL)) -

        2LL*Power(xi, 4LL)*Power(xj, 14LL)*

        (-13377735LL*xj + 13203540LL*r*Power(xj, 2LL) +

            9346050LL*Power(r, 2LL)*Power(xj, 3LL) +

            2195760LL*Power(r, 3LL)*Power(xj, 4LL) +

            240660LL*Power(r, 4LL)*Power(xj, 5LL) + 11088LL*Power(r, 5LL)*Power(xj, 6LL) +

            56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        Power(xi, 16LL)*Power(xj, 2LL)*

        (85995LL*xj + 147420LL*r*Power(xj, 2LL) + 122850LL*Power(r, 2LL)*Power(xj, 3LL) +

            65520LL*Power(r, 3LL)*Power(xj, 4LL) + 24570LL*Power(r, 4LL)*Power(xj, 5LL) +

            6552LL*Power(r, 5LL)*Power(xj, 6LL) + 308LL*Power(r, 6LL)*Power(xj, 7LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 9LL)*

                 (Power(xi, 2LL)*Power(xj, 16LL)*

        (70073640LL + 47669895LL*r*xj + 13931190LL*Power(r, 2LL)*Power(xj, 2LL) +

            2170350LL*Power(r, 3LL)*Power(xj, 3LL) + 169260LL*Power(r, 4LL)*Power(xj, 4LL) +

            1638LL*Power(r, 5LL)*Power(xj, 5LL) - 756LL*Power(r, 6LL)*Power(xj, 6LL) -

            44LL*Power(r, 7LL)*Power(xj, 7LL)) +

        364LL*Power(xi, 10LL)*Power(xj, 8LL)*

        (-7425LL - 13860LL*r*xj - 5940LL*Power(r, 2LL)*Power(xj, 2LL) -

            11880LL*Power(r, 3LL)*Power(xj, 3LL) - 2640LL*Power(r, 4LL)*Power(xj, 4LL) -

            45LL*Power(r, 5LL)*Power(xj, 5LL) + 30LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

        364LL*Power(xi, 8LL)*Power(xj, 10LL)*

        (-20925LL + 18270LL*r*xj - 58320LL*Power(r, 2LL)*Power(xj, 2LL) -

            17730LL*Power(r, 3LL)*Power(xj, 3LL) - 300LL*Power(r, 4LL)*Power(xj, 4LL) +

            423LL*Power(r, 5LL)*Power(xj, 5LL) + 54LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) -

        3LL*Power(xi, 18LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                             1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                             126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                             4LL*Power(r, 7LL)*Power(xj, 7LL)) +

        3LL*Power(xj, 18LL)*(1801800LL + 1576575LL*r*xj +

                             630630LL*Power(r, 2LL)*Power(xj, 2LL) + 150150LL*Power(r, 3LL)*Power(xj, 3LL) +

                             23100LL*Power(r, 4LL)*Power(xj, 4LL) + 2310LL*Power(r, 5LL)*Power(xj, 5LL) +

                             140LL*Power(r, 6LL)*Power(xj, 6LL) + 4LL*Power(r, 7LL)*Power(xj, 7LL)) +

        2LL*Power(xi, 14LL)*Power(xj, 4LL)*

        (-147420LL - 257985LL*r*xj - 221130LL*Power(r, 2LL)*Power(xj, 2LL) -

            122850LL*Power(r, 3LL)*Power(xj, 3LL) - 49140LL*Power(r, 4LL)*Power(xj, 4LL) -

            17388LL*Power(r, 5LL)*Power(xj, 5LL) - 1512LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

        42LL*Power(xi, 12LL)*Power(xj, 6LL)*

        (-25740LL - 45045LL*r*xj - 38610LL*Power(r, 2LL)*Power(xj, 2LL) -

            19470LL*Power(r, 3LL)*Power(xj, 3LL) - 12540LL*Power(r, 4LL)*Power(xj, 4LL) -

            1836LL*Power(r, 5LL)*Power(xj, 5LL) - 8LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        42LL*Power(xi, 6LL)*Power(xj, 12LL)*

        (921600LL - 1640835LL*r*xj - 546030LL*Power(r, 2LL)*Power(xj, 2LL) +

            20730LL*Power(r, 3LL)*Power(xj, 3LL) + 30180LL*Power(r, 4LL)*Power(xj, 4LL) +

            5028LL*Power(r, 5LL)*Power(xj, 5LL) + 344LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) -

        2LL*Power(xi, 4LL)*Power(xj, 14LL)*

        (-67767840LL - 13377735LL*r*xj + 6601770LL*Power(r, 2LL)*Power(xj, 2LL) +

            3115350LL*Power(r, 3LL)*Power(xj, 3LL) + 548940LL*Power(r, 4LL)*Power(xj, 4LL) +

            48132LL*Power(r, 5LL)*Power(xj, 5LL) + 1848LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        Power(xi, 16LL)*Power(xj, 2LL)*

        (49140LL + 85995LL*r*xj + 73710LL*Power(r, 2LL)*Power(xj, 2LL) +

            40950LL*Power(r, 3LL)*Power(xj, 3LL) + 16380LL*Power(r, 4LL)*Power(xj, 4LL) +

            4914LL*Power(r, 5LL)*Power(xj, 5LL) + 1092LL*Power(r, 6LL)*Power(xj, 6LL) +

            44LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (3780LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 13LL)*Power(xi + xj, 13LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_4S_3S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_3S_4S(r, xj, xi);
}
