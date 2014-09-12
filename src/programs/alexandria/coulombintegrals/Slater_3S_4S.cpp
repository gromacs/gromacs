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

cl_R Slater_3S_4S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (1363LL*xi)/6144LL

            ;
        }
        else
        {
            S = (1LL/r)*((-74724249600LL + 74724249600LL*exp(2LL*rxi) - 132871488750LL*rxi -

                          116294478300LL*Power(rxi, 2LL) - 66678987375LL*Power(rxi, 3LL) -

                          28114836750LL*Power(rxi, 4LL) - 9274044780LL*Power(rxi, 5LL) -

                          2484321840LL*Power(rxi, 6LL) - 553204080LL*Power(rxi, 7LL) -

                          103783680LL*Power(rxi, 8LL) - 16473600LL*Power(rxi, 9LL) - 2196480LL*Power(rxi, 10LL) -

                          239616LL*Power(rxi, 11LL) - 19968LL*Power(rxi, 12LL) - 1024LL*Power(rxi, 13LL))/

                         (7.47242496e10*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(3LL*Power(xi, 12LL) + 39LL*Power(xi, 11LL)*xj + 234LL*Power(xi, 10LL)*Power(xj, 2LL) +

                        858LL*Power(xi, 9LL)*Power(xj, 3LL) + 2145LL*Power(xi, 8LL)*Power(xj, 4LL) +

                        3861LL*Power(xi, 7LL)*Power(xj, 5LL) + 5148LL*Power(xi, 6LL)*Power(xj, 6LL) +

                        5148LL*Power(xi, 5LL)*Power(xj, 7LL) + 2860LL*Power(xi, 4LL)*Power(xj, 8LL) +

                        1144LL*Power(xi, 3LL)*Power(xj, 9LL) + 312LL*Power(xi, 2LL)*Power(xj, 10LL) +

                        52LL*xi*Power(xj, 11LL) + 4LL*Power(xj, 12LL)))/(12LL*Power(xi + xj, 13LL))

            ;
        }
        else
        {
            S = (1LL/r)*((3780LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 13LL) +

                          84LL*exp(2LL*rxj)*Power(rxj, 10LL)*

                          (-60LL*Power(rxi, 20LL) - 2LL*Power(rxi, 21LL) + 45LL*Power(rxj, 16LL) +

                           75LL*rxi*Power(rxj, 16LL) - 4LL*Power(rxi, 19LL)*(195LL + Power(rxj, 2LL)) +

                           15LL*Power(rxi, 3LL)*Power(rxj, 14LL)*(-65LL + 2LL*Power(rxj, 2LL)) +

                           15LL*Power(rxi, 2LL)*Power(rxj, 14LL)*(-39LL + 4LL*Power(rxj, 2LL)) -

                           30LL*Power(rxi, 18LL)*(182LL + 9LL*Power(rxj, 2LL)) +

                           30LL*Power(rxi, 13LL)*Power(rxj, 4LL)*(-13047LL + 377LL*Power(rxj, 2LL)) +

                           2LL*Power(rxi, 5LL)*Power(rxj, 12LL)*

                           (2925LL - 195LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           10LL*Power(rxi, 4LL)*Power(rxj, 12LL)*

                           (351LL - 78LL*Power(rxj, 2LL) + Power(rxj, 4LL)) -

                           130LL*Power(rxi, 6LL)*Power(rxj, 10LL)*

                           (99LL - 36LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           13LL*Power(rxi, 11LL)*Power(rxj, 6LL)*

                           (30735LL - 1650LL*Power(rxj, 2LL) + 4LL*Power(rxj, 4LL)) +

                           Power(rxi, 7LL)*Power(rxj, 10LL)*

                           (-15015LL + 3330LL*Power(rxj, 2LL) + 4LL*Power(rxj, 4LL)) +

                           210LL*Power(rxi, 16LL)*(-156LL - 262LL*Power(rxj, 2LL) + 5LL*Power(rxj, 4LL)) -

                           6LL*Power(rxi, 9LL)*Power(rxj, 8LL)*

                           (-48620LL - 715LL*Power(rxj, 2LL) + 6LL*Power(rxj, 4LL)) +

                           3LL*Power(rxi, 17LL)*(-6825LL - 1870LL*Power(rxj, 2LL) + 12LL*Power(rxj, 4LL)) -

                           30LL*Power(rxi, 14LL)*Power(rxj, 2LL)*

                           (17934LL - 12LL*Power(rxj, 2LL) + 13LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 8LL)*Power(rxj, 8LL)*

                           (2145LL + 2860LL*Power(rxj, 2LL) + 14LL*Power(rxj, 4LL)) +

                           65LL*Power(rxi, 10LL)*Power(rxj, 6LL)*

                           (-13725LL - 792LL*Power(rxj, 2LL) + 22LL*Power(rxj, 4LL)) -

                           10LL*Power(rxi, 12LL)*Power(rxj, 4LL)*

                           (153630LL - 15054LL*Power(rxj, 2LL) + 143LL*Power(rxj, 4LL)) +

                           Power(rxi, 15LL)*(-269325LL*Power(rxj, 2LL) + 9270LL*Power(rxj, 4LL) -

                                             52LL*Power(rxj, 6LL))) + exp(2LL*rxi)*Power(rxi, 8LL)*

                          (Power(rxi, 2LL)*Power(rxj, 16LL)*

                           (70073640LL + 47669895LL*rxj + 13931190LL*Power(rxj, 2LL) +

           2170350LL*Power(rxj, 3LL) + 169260LL*Power(rxj, 4LL) + 1638LL*Power(rxj, 5LL) -

           756LL*Power(rxj, 6LL) - 44LL*Power(rxj, 7LL)) +

                           364LL*Power(rxi, 10LL)*Power(rxj, 8LL)*

                           (-7425LL - 13860LL*rxj - 5940LL*Power(rxj, 2LL) - 11880LL*Power(rxj, 3LL) -

           2640LL*Power(rxj, 4LL) - 45LL*Power(rxj, 5LL) + 30LL*Power(rxj, 6LL) +

           2LL*Power(rxj, 7LL)) - 364LL*Power(rxi, 8LL)*Power(rxj, 10LL)*

                           (-20925LL + 18270LL*rxj - 58320LL*Power(rxj, 2LL) - 17730LL*Power(rxj, 3LL) -

           300LL*Power(rxj, 4LL) + 423LL*Power(rxj, 5LL) + 54LL*Power(rxj, 6LL) +

           2LL*Power(rxj, 7LL)) - 3LL*Power(rxi, 18LL)*

                           (1260LL + 2205LL*rxj + 1890LL*Power(rxj, 2LL) + 1050LL*Power(rxj, 3LL) +

           420LL*Power(rxj, 4LL) + 126LL*Power(rxj, 5LL) + 28LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL)) + 3LL*Power(rxj, 18LL)*

                           (1801800LL + 1576575LL*rxj + 630630LL*Power(rxj, 2LL) +

           150150LL*Power(rxj, 3LL) + 23100LL*Power(rxj, 4LL) + 2310LL*Power(rxj, 5LL) +

           140LL*Power(rxj, 6LL) + 4LL*Power(rxj, 7LL)) +

                           2LL*Power(rxi, 14LL)*Power(rxj, 4LL)*

                           (-147420LL - 257985LL*rxj - 221130LL*Power(rxj, 2LL) - 122850LL*Power(rxj, 3LL) -

           49140LL*Power(rxj, 4LL) - 17388LL*Power(rxj, 5LL) - 1512LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) - 42LL*Power(rxi, 12LL)*Power(rxj, 6LL)*

                           (-25740LL - 45045LL*rxj - 38610LL*Power(rxj, 2LL) - 19470LL*Power(rxj, 3LL) -

           12540LL*Power(rxj, 4LL) - 1836LL*Power(rxj, 5LL) - 8LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) + 42LL*Power(rxi, 6LL)*Power(rxj, 12LL)*

                           (921600LL - 1640835LL*rxj - 546030LL*Power(rxj, 2LL) + 20730LL*Power(rxj, 3LL) +

           30180LL*Power(rxj, 4LL) + 5028LL*Power(rxj, 5LL) + 344LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) - 2LL*Power(rxi, 4LL)*Power(rxj, 14LL)*

                           (-67767840LL - 13377735LL*rxj + 6601770LL*Power(rxj, 2LL) +

           3115350LL*Power(rxj, 3LL) + 548940LL*Power(rxj, 4LL) + 48132LL*Power(rxj, 5LL) +

           1848LL*Power(rxj, 6LL) + 8LL*Power(rxj, 7LL)) +

                           Power(rxi, 16LL)*Power(rxj, 2LL)*

                           (49140LL + 85995LL*rxj + 73710LL*Power(rxj, 2LL) + 40950LL*Power(rxj, 3LL) +

           16380LL*Power(rxj, 4LL) + 4914LL*Power(rxj, 5LL) + 1092LL*Power(rxj, 6LL) +

           44LL*Power(rxj, 7LL))))/

                         (3780LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 13LL)*Power(rxi + rxj, 13LL))

                         );
        }

    }
    return S;
}


cl_R Slater_4S_3S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_3S_4S(r, xj, xi);
}
