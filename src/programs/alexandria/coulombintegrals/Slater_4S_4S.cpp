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

cl_R Slater_4S_4S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (26333LL*xi)/131072LL

            ;
        }
        else
        {
            S = (1LL/r)*((-83691159552000LL + 83691159552000LL*exp(2LL*rxi) - 150568359566625LL*rxi -

                          133754400029250LL*Power(rxi, 2LL) - 78142908343500LL*Power(rxi, 3LL) -

                          33740723016000LL*Power(rxi, 4LL) - 11470756096800LL*Power(rxi, 5LL) -

                          3193358968800LL*Power(rxi, 6LL) - 747112766400LL*Power(rxi, 7LL) -

                          149448499200LL*Power(rxi, 8LL) - 25830604800LL*Power(rxi, 9LL) -

                          3874590720LL*Power(rxi, 10LL) - 503193600LL*Power(rxi, 11LL) -

                          55910400LL*Power(rxi, 12LL) - 5160960LL*Power(rxi, 13LL) - 368640LL*Power(rxi, 14LL) -

                          16384LL*Power(rxi, 15LL))/(8.3691159552e13*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 14LL) + 15LL*Power(xi, 13LL)*xj + 105LL*Power(xi, 12LL)*Power(xj, 2LL) +

                        455LL*Power(xi, 11LL)*Power(xj, 3LL) + 1365LL*Power(xi, 10LL)*Power(xj, 4LL) +

                        3003LL*Power(xi, 9LL)*Power(xj, 5LL) + 5005LL*Power(xi, 8LL)*Power(xj, 6LL) +

                        6435LL*Power(xi, 7LL)*Power(xj, 7LL) + 5005LL*Power(xi, 6LL)*Power(xj, 8LL) +

                        3003LL*Power(xi, 5LL)*Power(xj, 9LL) + 1365LL*Power(xi, 4LL)*Power(xj, 10LL) +

                        455LL*Power(xi, 3LL)*Power(xj, 11LL) + 105LL*Power(xi, 2LL)*Power(xj, 12LL) +

                        15LL*xi*Power(xj, 13LL) + Power(xj, 14LL)))/(4LL*Power(xi + xj, 15LL))

            ;
        }
        else
        {
            S = (1LL/r)*((1260LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 15LL) +

                          exp(2LL*rxj)*Power(rxj, 10LL)*

                          (-3276LL*Power(rxi, 25LL) - 168LL*Power(rxi, 26LL) - 4LL*Power(rxi, 27LL) +

                           1260LL*Power(rxj, 20LL) + 2205LL*rxi*Power(rxj, 20LL) +

                           1890LL*Power(rxi, 2LL)*Power(rxj, 18LL)*(-10LL + Power(rxj, 2LL)) -

                           420LL*Power(rxi, 24LL)*(91LL + Power(rxj, 2LL)) +

                           525LL*Power(rxi, 3LL)*Power(rxj, 18LL)*(-63LL + 2LL*Power(rxj, 2LL)) +

                           42LL*Power(rxi, 23LL)*(-6825LL - 405LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           63LL*Power(rxi, 5LL)*Power(rxj, 16LL)*

                           (3675LL - 250LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           210LL*Power(rxi, 4LL)*Power(rxj, 16LL)*

                           (630LL - 135LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           252LL*Power(rxi, 22LL)*(-5460LL - 1225LL*Power(rxj, 2LL) + 17LL*Power(rxj, 4LL)) -

                           1260LL*Power(rxi, 17LL)*Power(rxj, 4LL)*

                           (141729LL - 10145LL*Power(rxj, 2LL) + 116LL*Power(rxj, 4LL)) +

                           21LL*Power(rxi, 9LL)*Power(rxj, 12LL)*

                           (164775LL - 18460LL*Power(rxj, 2LL) + 828LL*Power(rxj, 4LL)) +

                           14LL*Power(rxi, 6LL)*Power(rxj, 14LL)*

                           (-40950LL + 14175LL*Power(rxj, 2LL) - 450LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL)) -

                           210LL*Power(rxi, 8LL)*Power(rxj, 12LL)*

                           (-8190LL + 4095LL*Power(rxj, 2LL) - 210LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL)) +

                           42LL*Power(rxi, 10LL)*Power(rxj, 10LL)*

                           (-209430LL - 2925LL*Power(rxj, 2LL) - 8840LL*Power(rxj, 4LL) + 4LL*Power(rxj, 6LL))

                           + Power(rxi, 7LL)*Power(rxj, 14LL)*(-1003275LL + 110250LL*Power(rxj, 2LL) -

                                                               1890LL*Power(rxj, 4LL) + 4LL*Power(rxj, 6LL)) -

                           21LL*Power(rxi, 11LL)*Power(rxj, 10LL)*

                           (-1033695LL - 218400LL*Power(rxj, 2LL) + 552LL*Power(rxj, 4LL) +

           4LL*Power(rxj, 6LL)) + 280LL*Power(rxi, 18LL)*Power(rxj, 2LL)*

                           (-385560LL - 73953LL*Power(rxj, 2LL) + 2370LL*Power(rxj, 4LL) + 4LL*Power(rxj, 6LL))

                           - 35LL*Power(rxi, 15LL)*Power(rxj, 6LL)*

                           (-1565613LL + 359520LL*Power(rxj, 2LL) - 7020LL*Power(rxj, 4LL) +

           8LL*Power(rxj, 6LL)) + 14LL*Power(rxi, 19LL)*Power(rxj, 2LL)*

                           (-4980150LL + 126765LL*Power(rxj, 2LL) - 3852LL*Power(rxj, 4LL) +

           20LL*Power(rxj, 6LL)) - 630LL*Power(rxi, 14LL)*Power(rxj, 6LL)*

                           (708714LL - 14385LL*Power(rxj, 2LL) - 2340LL*Power(rxj, 4LL) + 20LL*Power(rxj, 6LL))

                           + 210LL*Power(rxi, 16LL)*Power(rxj, 4LL)*

                           (-2087532LL + 328491LL*Power(rxj, 2LL) - 11740LL*Power(rxj, 4LL) +

           52LL*Power(rxj, 6LL)) - 84LL*Power(rxi, 20LL)*

                           (59670LL + 236250LL*Power(rxj, 2LL) - 8745LL*Power(rxj, 4LL) + 92LL*Power(rxj, 6LL))

                           - 2LL*Power(rxi, 21LL)*(1949220LL + 1598625LL*Power(rxj, 2LL) - 41391LL*Power(rxj, 4LL) +

                                                   128LL*Power(rxj, 6LL)) + Power(rxi, 13LL)*Power(rxj, 8LL)*

                           (173037375LL - 2784600LL*Power(rxj, 2LL) - 112140LL*Power(rxj, 4LL) +

           256LL*Power(rxj, 6LL)) + 14LL*Power(rxi, 12LL)*Power(rxj, 8LL)*

                           (-7260750LL - 2521935LL*Power(rxj, 2LL) + 19500LL*Power(rxj, 4LL) +

           344LL*Power(rxj, 6LL))) +

                          exp(2LL*rxi)*Power(rxi, 10LL)*

                          (210LL*Power(rxi, 2LL)*Power(rxj, 18LL)*

                           (514080LL + 332010LL*rxj + 94500LL*Power(rxj, 2LL) + 15225LL*Power(rxj, 3LL) +

           1470LL*Power(rxj, 4LL) + 81LL*Power(rxj, 5LL) + 2LL*Power(rxj, 6LL)) +

                           105LL*Power(rxi, 18LL)*Power(rxj, 2LL)*

                           (180LL + 315LL*rxj + 270LL*Power(rxj, 2LL) + 150LL*Power(rxj, 3LL) +

           60LL*Power(rxj, 4LL) + 18LL*Power(rxj, 5LL) + 4LL*Power(rxj, 6LL)) -

                           1365LL*Power(rxi, 10LL)*Power(rxj, 10LL)*

                           (-6444LL + 15903LL*rxj - 25866LL*Power(rxj, 2LL) - 2040LL*Power(rxj, 3LL) +

           1080LL*Power(rxj, 4LL) + 180LL*Power(rxj, 5LL) + 8LL*Power(rxj, 6LL)) +

                           Power(rxi, 14LL)*Power(rxj, 6LL)*

                           (573300LL + 1003275LL*rxj + 859950LL*Power(rxj, 2LL) + 387660LL*Power(rxj, 3LL) +

           371280LL*Power(rxj, 4LL) + 11592LL*Power(rxj, 5LL) - 4816LL*Power(rxj, 6LL) -

           256LL*Power(rxj, 7LL)) + 2LL*Power(rxj, 20LL)*

                           (2506140LL + 1949220LL*rxj + 687960LL*Power(rxj, 2LL) +

           143325LL*Power(rxj, 3LL) + 19110LL*Power(rxj, 4LL) + 1638LL*Power(rxj, 5LL) +

           84LL*Power(rxj, 6LL) + 2LL*Power(rxj, 7LL)) -

                           42LL*Power(rxi, 4LL)*Power(rxj, 16LL)*

                           (-10437660LL - 4251870LL*rxj - 493020LL*Power(rxj, 2LL) +

           42255LL*Power(rxj, 3LL) + 17490LL*Power(rxj, 4LL) + 1971LL*Power(rxj, 5LL) +

           102LL*Power(rxj, 6LL) + 2LL*Power(rxj, 7LL)) +

                           21LL*Power(rxi, 16LL)*Power(rxj, 4LL)*

                           (-6300LL - 11025LL*rxj - 9450LL*Power(rxj, 2LL) - 5250LL*Power(rxj, 3LL) -

           2100LL*Power(rxj, 4LL) - 828LL*Power(rxj, 5LL) - 8LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL)) - Power(rxi, 20LL)*

                           (1260LL + 2205LL*rxj + 1890LL*Power(rxj, 2LL) + 1050LL*Power(rxj, 3LL) +

           420LL*Power(rxj, 4LL) + 126LL*Power(rxj, 5LL) + 28LL*Power(rxj, 6LL) +

           4LL*Power(rxj, 7LL)) - 35LL*Power(rxi, 8LL)*Power(rxj, 12LL)*

                           (-2904300LL + 4943925LL*rxj + 258930LL*Power(rxj, 2LL) -

           359520LL*Power(rxj, 3LL) - 70440LL*Power(rxj, 4LL) - 4176LL*Power(rxj, 5LL) +

           32LL*Power(rxj, 6LL) + 8LL*Power(rxj, 7LL)) +

                           35LL*Power(rxi, 12LL)*Power(rxj, 8LL)*

                           (-49140LL - 98865LL*rxj + 3510LL*Power(rxj, 2LL) - 131040LL*Power(rxj, 3LL) -

           7800LL*Power(rxj, 4LL) + 3204LL*Power(rxj, 5LL) + 360LL*Power(rxj, 6LL) +

           8LL*Power(rxj, 7LL)) + Power(rxi, 6LL)*Power(rxj, 14LL)*

                           (446489820LL - 54796455LL*rxj - 68983110LL*Power(rxj, 2LL) -

           12782700LL*Power(rxj, 3LL) - 663600LL*Power(rxj, 4LL) + 53928LL*Power(rxj, 5LL) +

           7728LL*Power(rxj, 6LL) + 256LL*Power(rxj, 7LL))))/

                         (1260LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 15LL)*Power(rxi + rxj, 15LL))

                         );
        }

    }
    return S;
}
