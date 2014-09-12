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

cl_R DSlater_4S_4S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-150568359566625LL*xi + 167382319104000LL*exp(2LL*r*xi)*xi -

                  267508800058500LL*r*Power(xi, 2LL) -

                  234428725030500LL*Power(r, 2LL)*Power(xi, 3LL) -

                  134962892064000LL*Power(r, 3LL)*Power(xi, 4LL) -

                  57353780484000LL*Power(r, 4LL)*Power(xi, 5LL) -

                  19160153812800LL*Power(r, 5LL)*Power(xi, 6LL) -

                  5229789364800LL*Power(r, 6LL)*Power(xi, 7LL) -

                  1195587993600LL*Power(r, 7LL)*Power(xi, 8LL) -

                  232475443200LL*Power(r, 8LL)*Power(xi, 9LL) -

                  38745907200LL*Power(r, 9LL)*Power(xi, 10LL) -

                  5535129600LL*Power(r, 10LL)*Power(xi, 11LL) -

                  670924800LL*Power(r, 11LL)*Power(xi, 12LL) -

                  67092480LL*Power(r, 12LL)*Power(xi, 13LL) -

                  5160960LL*Power(r, 13LL)*Power(xi, 14LL) - 245760LL*Power(r, 14LL)*Power(xi, 15LL))/

                (8.3691159552e13*exp(2LL*r*xi)*r) +

                (-83691159552000LL + 83691159552000LL*exp(2LL*r*xi) -

                 150568359566625LL*r*xi - 133754400029250LL*Power(r, 2LL)*Power(xi, 2LL) -

                 78142908343500LL*Power(r, 3LL)*Power(xi, 3LL) -

                 33740723016000LL*Power(r, 4LL)*Power(xi, 4LL) -

                 11470756096800LL*Power(r, 5LL)*Power(xi, 5LL) -

                 3193358968800LL*Power(r, 6LL)*Power(xi, 6LL) -

                 747112766400LL*Power(r, 7LL)*Power(xi, 7LL) -

                 149448499200LL*Power(r, 8LL)*Power(xi, 8LL) -

                 25830604800LL*Power(r, 9LL)*Power(xi, 9LL) -

                 3874590720LL*Power(r, 10LL)*Power(xi, 10LL) -

                 503193600LL*Power(r, 11LL)*Power(xi, 11LL) -

                 55910400LL*Power(r, 12LL)*Power(xi, 12LL) - 5160960LL*Power(r, 13LL)*Power(xi, 13LL) -

                 368640LL*Power(r, 14LL)*Power(xi, 14LL) - 16384LL*Power(r, 15LL)*Power(xi, 15LL))/

                (8.3691159552e13*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-83691159552000LL + 83691159552000LL*exp(2LL*r*xi) -

                     150568359566625LL*r*xi - 133754400029250LL*Power(r, 2LL)*Power(xi, 2LL) -

                     78142908343500LL*Power(r, 3LL)*Power(xi, 3LL) -

                     33740723016000LL*Power(r, 4LL)*Power(xi, 4LL) -

                     11470756096800LL*Power(r, 5LL)*Power(xi, 5LL) -

                     3193358968800LL*Power(r, 6LL)*Power(xi, 6LL) -

                     747112766400LL*Power(r, 7LL)*Power(xi, 7LL) -

                     149448499200LL*Power(r, 8LL)*Power(xi, 8LL) -

                     25830604800LL*Power(r, 9LL)*Power(xi, 9LL) -

                     3874590720LL*Power(r, 10LL)*Power(xi, 10LL) -

                     503193600LL*Power(r, 11LL)*Power(xi, 11LL) -

                     55910400LL*Power(r, 12LL)*Power(xi, 12LL) -

                     5160960LL*Power(r, 13LL)*Power(xi, 13LL) - 368640LL*Power(r, 14LL)*Power(xi, 14LL) -

                     16384LL*Power(r, 15LL)*Power(xi, 15LL)))/(4.1845579776e13*exp(2LL*r*xi)*r)

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
            S = (1260LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                 exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-3276LL*Power(r, 5LL)*Power(xi, 25LL) - 168LL*Power(r, 6LL)*Power(xi, 26LL) -

                  4LL*Power(r, 7LL)*Power(xi, 27LL) + 1260LL*Power(xj, 20LL) +

                  2205LL*r*xi*Power(xj, 20LL) +

                  1890LL*Power(xi, 2LL)*Power(xj, 18LL)*(-10LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  420LL*Power(r, 4LL)*Power(xi, 24LL)*(91LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  525LL*r*Power(xi, 3LL)*Power(xj, 18LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  42LL*Power(r, 3LL)*Power(xi, 23LL)*

                  (-6825LL - 405LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  63LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

                  (3675LL - 250LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  210LL*Power(xi, 4LL)*Power(xj, 16LL)*

                  (630LL - 135LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  252LL*Power(r, 2LL)*Power(xi, 22LL)*

                  (-5460LL - 1225LL*Power(r, 2LL)*Power(xj, 2LL) + 17LL*Power(r, 4LL)*Power(xj, 4LL))

                  - 1260LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

                  (141729LL - 10145LL*Power(r, 2LL)*Power(xj, 2LL) +

            116LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  21LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

                  (164775LL - 18460LL*Power(r, 2LL)*Power(xj, 2LL) +

            828LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  14LL*Power(xi, 6LL)*Power(xj, 14LL)*

                  (-40950LL + 14175LL*Power(r, 2LL)*Power(xj, 2LL) -

            450LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  210LL*Power(xi, 8LL)*Power(xj, 12LL)*

                  (-8190LL + 4095LL*Power(r, 2LL)*Power(xj, 2LL) -

            210LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  42LL*Power(xi, 10LL)*Power(xj, 10LL)*

                  (-209430LL - 2925LL*Power(r, 2LL)*Power(xj, 2LL) -

            8840LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  r*Power(xi, 7LL)*Power(xj, 14LL)*

                  (-1003275LL + 110250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1890LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  21LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

                  (-1033695LL - 218400LL*Power(r, 2LL)*Power(xj, 2LL) +

            552LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  280LL*Power(xi, 18LL)*Power(xj, 2LL)*

                  (-385560LL - 73953LL*Power(r, 2LL)*Power(xj, 2LL) +

            2370LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  35LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

                  (-1565613LL + 359520LL*Power(r, 2LL)*Power(xj, 2LL) -

            7020LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  14LL*r*Power(xi, 19LL)*Power(xj, 2LL)*

                  (-4980150LL + 126765LL*Power(r, 2LL)*Power(xj, 2LL) -

            3852LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  630LL*Power(xi, 14LL)*Power(xj, 6LL)*

                  (708714LL - 14385LL*Power(r, 2LL)*Power(xj, 2LL) -

            2340LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  210LL*Power(xi, 16LL)*Power(xj, 4LL)*

                  (-2087532LL + 328491LL*Power(r, 2LL)*Power(xj, 2LL) -

            11740LL*Power(r, 4LL)*Power(xj, 4LL) + 52LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  84LL*Power(xi, 20LL)*(59670LL + 236250LL*Power(r, 2LL)*Power(xj, 2LL) -

                                        8745LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  2LL*r*Power(xi, 21LL)*(1949220LL + 1598625LL*Power(r, 2LL)*Power(xj, 2LL) -

                                         41391LL*Power(r, 4LL)*Power(xj, 4LL) + 128LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  r*Power(xi, 13LL)*Power(xj, 8LL)*

                  (173037375LL - 2784600LL*Power(r, 2LL)*Power(xj, 2LL) -

            112140LL*Power(r, 4LL)*Power(xj, 4LL) + 256LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  14LL*Power(xi, 12LL)*Power(xj, 8LL)*

                  (-7260750LL - 2521935LL*Power(r, 2LL)*Power(xj, 2LL) +

            19500LL*Power(r, 4LL)*Power(xj, 4LL) + 344LL*Power(r, 6LL)*Power(xj, 6LL))) +

                 exp(2LL*r*xi)*Power(xi, 10LL)*

                 (210LL*Power(xi, 2LL)*Power(xj, 18LL)*

                  (514080LL + 332010LL*r*xj + 94500LL*Power(r, 2LL)*Power(xj, 2LL) +

            15225LL*Power(r, 3LL)*Power(xj, 3LL) + 1470LL*Power(r, 4LL)*Power(xj, 4LL) +

            81LL*Power(r, 5LL)*Power(xj, 5LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  105LL*Power(xi, 18LL)*Power(xj, 2LL)*

                  (180LL + 315LL*r*xj + 270LL*Power(r, 2LL)*Power(xj, 2LL) +

            150LL*Power(r, 3LL)*Power(xj, 3LL) + 60LL*Power(r, 4LL)*Power(xj, 4LL) +

            18LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  1365LL*Power(xi, 10LL)*Power(xj, 10LL)*

                  (-6444LL + 15903LL*r*xj - 25866LL*Power(r, 2LL)*Power(xj, 2LL) -

            2040LL*Power(r, 3LL)*Power(xj, 3LL) + 1080LL*Power(r, 4LL)*Power(xj, 4LL) +

            180LL*Power(r, 5LL)*Power(xj, 5LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  Power(xi, 14LL)*Power(xj, 6LL)*

                  (573300LL + 1003275LL*r*xj + 859950LL*Power(r, 2LL)*Power(xj, 2LL) +

            387660LL*Power(r, 3LL)*Power(xj, 3LL) + 371280LL*Power(r, 4LL)*Power(xj, 4LL) +

            11592LL*Power(r, 5LL)*Power(xj, 5LL) - 4816LL*Power(r, 6LL)*Power(xj, 6LL) -

            256LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  2LL*Power(xj, 20LL)*(2506140LL + 1949220LL*r*xj +

                                       687960LL*Power(r, 2LL)*Power(xj, 2LL) + 143325LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       19110LL*Power(r, 4LL)*Power(xj, 4LL) + 1638LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       84LL*Power(r, 6LL)*Power(xj, 6LL) + 2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  42LL*Power(xi, 4LL)*Power(xj, 16LL)*

                  (-10437660LL - 4251870LL*r*xj - 493020LL*Power(r, 2LL)*Power(xj, 2LL) +

            42255LL*Power(r, 3LL)*Power(xj, 3LL) + 17490LL*Power(r, 4LL)*Power(xj, 4LL) +

            1971LL*Power(r, 5LL)*Power(xj, 5LL) + 102LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  21LL*Power(xi, 16LL)*Power(xj, 4LL)*

                  (-6300LL - 11025LL*r*xj - 9450LL*Power(r, 2LL)*Power(xj, 2LL) -

            5250LL*Power(r, 3LL)*Power(xj, 3LL) - 2100LL*Power(r, 4LL)*Power(xj, 4LL) -

            828LL*Power(r, 5LL)*Power(xj, 5LL) - 8LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  Power(xi, 20LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   4LL*Power(r, 7LL)*Power(xj, 7LL)) -

                  35LL*Power(xi, 8LL)*Power(xj, 12LL)*

                  (-2904300LL + 4943925LL*r*xj + 258930LL*Power(r, 2LL)*Power(xj, 2LL) -

            359520LL*Power(r, 3LL)*Power(xj, 3LL) - 70440LL*Power(r, 4LL)*Power(xj, 4LL) -

            4176LL*Power(r, 5LL)*Power(xj, 5LL) + 32LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  35LL*Power(xi, 12LL)*Power(xj, 8LL)*

                  (-49140LL - 98865LL*r*xj + 3510LL*Power(r, 2LL)*Power(xj, 2LL) -

            131040LL*Power(r, 3LL)*Power(xj, 3LL) - 7800LL*Power(r, 4LL)*Power(xj, 4LL) +

            3204LL*Power(r, 5LL)*Power(xj, 5LL) + 360LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                  Power(xi, 6LL)*Power(xj, 14LL)*

                  (446489820LL - 54796455LL*r*xj - 68983110LL*Power(r, 2LL)*Power(xj, 2LL) -

            12782700LL*Power(r, 3LL)*Power(xj, 3LL) -

            663600LL*Power(r, 4LL)*Power(xj, 4LL) + 53928LL*Power(r, 5LL)*Power(xj, 5LL) +

            7728LL*Power(r, 6LL)*Power(xj, 6LL) + 256LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1260LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 15LL)*

                 Power(xi + xj, 15LL)) + (1260LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                                          exp(2LL*r*xj)*Power(xj, 10LL)*

                                          (-3276LL*Power(r, 5LL)*Power(xi, 25LL) - 168LL*Power(r, 6LL)*Power(xi, 26LL) -

                                  4LL*Power(r, 7LL)*Power(xi, 27LL) + 1260LL*Power(xj, 20LL) +

                                  2205LL*r*xi*Power(xj, 20LL) +

                                  1890LL*Power(xi, 2LL)*Power(xj, 18LL)*(-10LL + Power(r, 2LL)*Power(xj, 2LL)) -

                                  420LL*Power(r, 4LL)*Power(xi, 24LL)*(91LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  525LL*r*Power(xi, 3LL)*Power(xj, 18LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  42LL*Power(r, 3LL)*Power(xi, 23LL)*

                                  (-6825LL - 405LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  63LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

                                  (3675LL - 250LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  210LL*Power(xi, 4LL)*Power(xj, 16LL)*

                                  (630LL - 135LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  252LL*Power(r, 2LL)*Power(xi, 22LL)*

                                  (-5460LL - 1225LL*Power(r, 2LL)*Power(xj, 2LL) + 17LL*Power(r, 4LL)*Power(xj, 4LL))

                                  - 1260LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

                                  (141729LL - 10145LL*Power(r, 2LL)*Power(xj, 2LL) +

            116LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  21LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

                                  (164775LL - 18460LL*Power(r, 2LL)*Power(xj, 2LL) +

            828LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  14LL*Power(xi, 6LL)*Power(xj, 14LL)*

                                  (-40950LL + 14175LL*Power(r, 2LL)*Power(xj, 2LL) -

            450LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  210LL*Power(xi, 8LL)*Power(xj, 12LL)*

                                  (-8190LL + 4095LL*Power(r, 2LL)*Power(xj, 2LL) -

            210LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  42LL*Power(xi, 10LL)*Power(xj, 10LL)*

                                  (-209430LL - 2925LL*Power(r, 2LL)*Power(xj, 2LL) -

            8840LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  r*Power(xi, 7LL)*Power(xj, 14LL)*

                                  (-1003275LL + 110250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1890LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  21LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

                                  (-1033695LL - 218400LL*Power(r, 2LL)*Power(xj, 2LL) +

            552LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  280LL*Power(xi, 18LL)*Power(xj, 2LL)*

                                  (-385560LL - 73953LL*Power(r, 2LL)*Power(xj, 2LL) +

            2370LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  35LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

                                  (-1565613LL + 359520LL*Power(r, 2LL)*Power(xj, 2LL) -

            7020LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  14LL*r*Power(xi, 19LL)*Power(xj, 2LL)*

                                  (-4980150LL + 126765LL*Power(r, 2LL)*Power(xj, 2LL) -

            3852LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  630LL*Power(xi, 14LL)*Power(xj, 6LL)*

                                  (708714LL - 14385LL*Power(r, 2LL)*Power(xj, 2LL) -

            2340LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  210LL*Power(xi, 16LL)*Power(xj, 4LL)*

                                  (-2087532LL + 328491LL*Power(r, 2LL)*Power(xj, 2LL) -

            11740LL*Power(r, 4LL)*Power(xj, 4LL) + 52LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  84LL*Power(xi, 20LL)*(59670LL + 236250LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                        8745LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  2LL*r*Power(xi, 21LL)*(1949220LL + 1598625LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                         41391LL*Power(r, 4LL)*Power(xj, 4LL) + 128LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  r*Power(xi, 13LL)*Power(xj, 8LL)*

                                  (173037375LL - 2784600LL*Power(r, 2LL)*Power(xj, 2LL) -

            112140LL*Power(r, 4LL)*Power(xj, 4LL) + 256LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  14LL*Power(xi, 12LL)*Power(xj, 8LL)*

                                  (-7260750LL - 2521935LL*Power(r, 2LL)*Power(xj, 2LL) +

            19500LL*Power(r, 4LL)*Power(xj, 4LL) + 344LL*Power(r, 6LL)*Power(xj, 6LL))) +

                                          exp(2LL*r*xi)*Power(xi, 10LL)*

                                          (210LL*Power(xi, 2LL)*Power(xj, 18LL)*

                                  (514080LL + 332010LL*r*xj + 94500LL*Power(r, 2LL)*Power(xj, 2LL) +

            15225LL*Power(r, 3LL)*Power(xj, 3LL) + 1470LL*Power(r, 4LL)*Power(xj, 4LL) +

            81LL*Power(r, 5LL)*Power(xj, 5LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  105LL*Power(xi, 18LL)*Power(xj, 2LL)*

                                  (180LL + 315LL*r*xj + 270LL*Power(r, 2LL)*Power(xj, 2LL) +

            150LL*Power(r, 3LL)*Power(xj, 3LL) + 60LL*Power(r, 4LL)*Power(xj, 4LL) +

            18LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  1365LL*Power(xi, 10LL)*Power(xj, 10LL)*

                                  (-6444LL + 15903LL*r*xj - 25866LL*Power(r, 2LL)*Power(xj, 2LL) -

            2040LL*Power(r, 3LL)*Power(xj, 3LL) + 1080LL*Power(r, 4LL)*Power(xj, 4LL) +

            180LL*Power(r, 5LL)*Power(xj, 5LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  Power(xi, 14LL)*Power(xj, 6LL)*

                                  (573300LL + 1003275LL*r*xj + 859950LL*Power(r, 2LL)*Power(xj, 2LL) +

            387660LL*Power(r, 3LL)*Power(xj, 3LL) + 371280LL*Power(r, 4LL)*Power(xj, 4LL) +

            11592LL*Power(r, 5LL)*Power(xj, 5LL) - 4816LL*Power(r, 6LL)*Power(xj, 6LL) -

            256LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  2LL*Power(xj, 20LL)*(2506140LL + 1949220LL*r*xj +

                                                       687960LL*Power(r, 2LL)*Power(xj, 2LL) + 143325LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                       19110LL*Power(r, 4LL)*Power(xj, 4LL) + 1638LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                       84LL*Power(r, 6LL)*Power(xj, 6LL) + 2LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  42LL*Power(xi, 4LL)*Power(xj, 16LL)*

                                  (-10437660LL - 4251870LL*r*xj - 493020LL*Power(r, 2LL)*Power(xj, 2LL) +

            42255LL*Power(r, 3LL)*Power(xj, 3LL) + 17490LL*Power(r, 4LL)*Power(xj, 4LL) +

            1971LL*Power(r, 5LL)*Power(xj, 5LL) + 102LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  21LL*Power(xi, 16LL)*Power(xj, 4LL)*

                                  (-6300LL - 11025LL*r*xj - 9450LL*Power(r, 2LL)*Power(xj, 2LL) -

            5250LL*Power(r, 3LL)*Power(xj, 3LL) - 2100LL*Power(r, 4LL)*Power(xj, 4LL) -

            828LL*Power(r, 5LL)*Power(xj, 5LL) - 8LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  Power(xi, 20LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   4LL*Power(r, 7LL)*Power(xj, 7LL)) -

                                  35LL*Power(xi, 8LL)*Power(xj, 12LL)*

                                  (-2904300LL + 4943925LL*r*xj + 258930LL*Power(r, 2LL)*Power(xj, 2LL) -

            359520LL*Power(r, 3LL)*Power(xj, 3LL) - 70440LL*Power(r, 4LL)*Power(xj, 4LL) -

            4176LL*Power(r, 5LL)*Power(xj, 5LL) + 32LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  35LL*Power(xi, 12LL)*Power(xj, 8LL)*

                                  (-49140LL - 98865LL*r*xj + 3510LL*Power(r, 2LL)*Power(xj, 2LL) -

            131040LL*Power(r, 3LL)*Power(xj, 3LL) - 7800LL*Power(r, 4LL)*Power(xj, 4LL) +

            3204LL*Power(r, 5LL)*Power(xj, 5LL) + 360LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

                                  Power(xi, 6LL)*Power(xj, 14LL)*

                                  (446489820LL - 54796455LL*r*xj - 68983110LL*Power(r, 2LL)*Power(xj, 2LL) -

            12782700LL*Power(r, 3LL)*Power(xj, 3LL) -

            663600LL*Power(r, 4LL)*Power(xj, 4LL) + 53928LL*Power(r, 5LL)*Power(xj, 5LL) +

            7728LL*Power(r, 6LL)*Power(xj, 6LL) + 256LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (630LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 15LL)*Power(xi + xj, 14LL)) -

                (2520LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                 exp(2LL*r*xj)*Power(xj, 10LL)*

                 (-16380LL*Power(r, 4LL)*Power(xi, 25LL) - 1008LL*Power(r, 5LL)*Power(xi, 26LL) -

        28LL*Power(r, 6LL)*Power(xi, 27LL) -

        840LL*Power(r, 5LL)*Power(xi, 24LL)*Power(xj, 2LL) + 2205LL*xi*Power(xj, 20LL) +

        3780LL*r*Power(xi, 2LL)*Power(xj, 20LL) +

        2100LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 20LL) -

        1680LL*Power(r, 3LL)*Power(xi, 24LL)*(91LL + Power(r, 2LL)*Power(xj, 2LL)) +

        525LL*Power(xi, 3LL)*Power(xj, 18LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        42LL*Power(r, 3LL)*Power(xi, 23LL)*

        (-810LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        63LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

        (-500LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        210LL*Power(xi, 4LL)*Power(xj, 16LL)*

        (-270LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        252LL*Power(r, 2LL)*Power(xi, 22LL)*

        (-2450LL*r*Power(xj, 2LL) + 68LL*Power(r, 3LL)*Power(xj, 4LL)) -

        1260LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

        (-20290LL*r*Power(xj, 2LL) + 464LL*Power(r, 3LL)*Power(xj, 4LL)) +

        21LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

        (-36920LL*r*Power(xj, 2LL) + 3312LL*Power(r, 3LL)*Power(xj, 4LL)) +

        126LL*Power(r, 2LL)*Power(xi, 23LL)*

        (-6825LL - 405LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        63LL*Power(xi, 5LL)*Power(xj, 16LL)*

        (3675LL - 250LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        504LL*r*Power(xi, 22LL)*(-5460LL - 1225LL*Power(r, 2LL)*Power(xj, 2LL) +

                                 17LL*Power(r, 4LL)*Power(xj, 4LL)) -

        1260LL*Power(xi, 17LL)*Power(xj, 4LL)*

        (141729LL - 10145LL*Power(r, 2LL)*Power(xj, 2LL) +

            116LL*Power(r, 4LL)*Power(xj, 4LL)) +

        21LL*Power(xi, 9LL)*Power(xj, 12LL)*

        (164775LL - 18460LL*Power(r, 2LL)*Power(xj, 2LL) +

            828LL*Power(r, 4LL)*Power(xj, 4LL)) +

        14LL*Power(xi, 6LL)*Power(xj, 14LL)*

        (28350LL*r*Power(xj, 2LL) - 1800LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) -

        210LL*Power(xi, 8LL)*Power(xj, 12LL)*

        (8190LL*r*Power(xj, 2LL) - 840LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 10LL)*Power(xj, 10LL)*

        (-5850LL*r*Power(xj, 2LL) - 35360LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) +

        r*Power(xi, 7LL)*Power(xj, 14LL)*

        (220500LL*r*Power(xj, 2LL) - 7560LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) -

        21LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

        (-436800LL*r*Power(xj, 2LL) + 2208LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) +

        280LL*Power(xi, 18LL)*Power(xj, 2LL)*

        (-147906LL*r*Power(xj, 2LL) + 9480LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) -

        35LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

        (719040LL*r*Power(xj, 2LL) - 28080LL*Power(r, 3LL)*Power(xj, 4LL) +

            48LL*Power(r, 5LL)*Power(xj, 6LL)) +

        14LL*r*Power(xi, 19LL)*Power(xj, 2LL)*

        (253530LL*r*Power(xj, 2LL) - 15408LL*Power(r, 3LL)*Power(xj, 4LL) +

            120LL*Power(r, 5LL)*Power(xj, 6LL)) -

        630LL*Power(xi, 14LL)*Power(xj, 6LL)*

        (-28770LL*r*Power(xj, 2LL) - 9360LL*Power(r, 3LL)*Power(xj, 4LL) +

            120LL*Power(r, 5LL)*Power(xj, 6LL)) +

        210LL*Power(xi, 16LL)*Power(xj, 4LL)*

        (656982LL*r*Power(xj, 2LL) - 46960LL*Power(r, 3LL)*Power(xj, 4LL) +

            312LL*Power(r, 5LL)*Power(xj, 6LL)) -

        84LL*Power(xi, 20LL)*(472500LL*r*Power(xj, 2LL) -

                              34980LL*Power(r, 3LL)*Power(xj, 4LL) + 552LL*Power(r, 5LL)*Power(xj, 6LL)) -

        2LL*r*Power(xi, 21LL)*(3197250LL*r*Power(xj, 2LL) -

                               165564LL*Power(r, 3LL)*Power(xj, 4LL) + 768LL*Power(r, 5LL)*Power(xj, 6LL)) +

        r*Power(xi, 13LL)*Power(xj, 8LL)*

        (-5569200LL*r*Power(xj, 2LL) - 448560LL*Power(r, 3LL)*Power(xj, 4LL) +

            1536LL*Power(r, 5LL)*Power(xj, 6LL)) +

        14LL*Power(xi, 12LL)*Power(xj, 8LL)*

        (-5043870LL*r*Power(xj, 2LL) + 78000LL*Power(r, 3LL)*Power(xj, 4LL) +

            2064LL*Power(r, 5LL)*Power(xj, 6LL)) +

        Power(xi, 7LL)*Power(xj, 14LL)*

        (-1003275LL + 110250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1890LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

        21LL*Power(xi, 11LL)*Power(xj, 10LL)*

        (-1033695LL - 218400LL*Power(r, 2LL)*Power(xj, 2LL) +

            552LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

        35LL*Power(xi, 15LL)*Power(xj, 6LL)*

        (-1565613LL + 359520LL*Power(r, 2LL)*Power(xj, 2LL) -

            7020LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

        14LL*Power(xi, 19LL)*Power(xj, 2LL)*

        (-4980150LL + 126765LL*Power(r, 2LL)*Power(xj, 2LL) -

            3852LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

        2LL*Power(xi, 21LL)*(1949220LL + 1598625LL*Power(r, 2LL)*Power(xj, 2LL) -

                             41391LL*Power(r, 4LL)*Power(xj, 4LL) + 128LL*Power(r, 6LL)*Power(xj, 6LL)) +

        Power(xi, 13LL)*Power(xj, 8LL)*

        (173037375LL - 2784600LL*Power(r, 2LL)*Power(xj, 2LL) -

            112140LL*Power(r, 4LL)*Power(xj, 4LL) + 256LL*Power(r, 6LL)*Power(xj, 6LL))) +

                 2LL*exp(2LL*r*xj)*Power(xj, 11LL)*

                 (-3276LL*Power(r, 5LL)*Power(xi, 25LL) - 168LL*Power(r, 6LL)*Power(xi, 26LL) -

        4LL*Power(r, 7LL)*Power(xi, 27LL) + 1260LL*Power(xj, 20LL) +

        2205LL*r*xi*Power(xj, 20LL) +

        1890LL*Power(xi, 2LL)*Power(xj, 18LL)*(-10LL + Power(r, 2LL)*Power(xj, 2LL)) -

        420LL*Power(r, 4LL)*Power(xi, 24LL)*(91LL + Power(r, 2LL)*Power(xj, 2LL)) +

        525LL*r*Power(xi, 3LL)*Power(xj, 18LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        42LL*Power(r, 3LL)*Power(xi, 23LL)*

        (-6825LL - 405LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        63LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

        (3675LL - 250LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        210LL*Power(xi, 4LL)*Power(xj, 16LL)*

        (630LL - 135LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        252LL*Power(r, 2LL)*Power(xi, 22LL)*

        (-5460LL - 1225LL*Power(r, 2LL)*Power(xj, 2LL) + 17LL*Power(r, 4LL)*Power(xj, 4LL))

        - 1260LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

        (141729LL - 10145LL*Power(r, 2LL)*Power(xj, 2LL) +

            116LL*Power(r, 4LL)*Power(xj, 4LL)) +

        21LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

        (164775LL - 18460LL*Power(r, 2LL)*Power(xj, 2LL) +

            828LL*Power(r, 4LL)*Power(xj, 4LL)) +

        14LL*Power(xi, 6LL)*Power(xj, 14LL)*

        (-40950LL + 14175LL*Power(r, 2LL)*Power(xj, 2LL) -

            450LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

        210LL*Power(xi, 8LL)*Power(xj, 12LL)*

        (-8190LL + 4095LL*Power(r, 2LL)*Power(xj, 2LL) -

            210LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 10LL)*Power(xj, 10LL)*

        (-209430LL - 2925LL*Power(r, 2LL)*Power(xj, 2LL) -

            8840LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

        r*Power(xi, 7LL)*Power(xj, 14LL)*

        (-1003275LL + 110250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1890LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

        21LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

        (-1033695LL - 218400LL*Power(r, 2LL)*Power(xj, 2LL) +

            552LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

        280LL*Power(xi, 18LL)*Power(xj, 2LL)*

        (-385560LL - 73953LL*Power(r, 2LL)*Power(xj, 2LL) +

            2370LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

        35LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

        (-1565613LL + 359520LL*Power(r, 2LL)*Power(xj, 2LL) -

            7020LL*Power(r, 4LL)*Power(xj, 4LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

        14LL*r*Power(xi, 19LL)*Power(xj, 2LL)*

        (-4980150LL + 126765LL*Power(r, 2LL)*Power(xj, 2LL) -

            3852LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

        630LL*Power(xi, 14LL)*Power(xj, 6LL)*

        (708714LL - 14385LL*Power(r, 2LL)*Power(xj, 2LL) -

            2340LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

        210LL*Power(xi, 16LL)*Power(xj, 4LL)*

        (-2087532LL + 328491LL*Power(r, 2LL)*Power(xj, 2LL) -

            11740LL*Power(r, 4LL)*Power(xj, 4LL) + 52LL*Power(r, 6LL)*Power(xj, 6LL)) -

        84LL*Power(xi, 20LL)*(59670LL + 236250LL*Power(r, 2LL)*Power(xj, 2LL) -

                              8745LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

        2LL*r*Power(xi, 21LL)*(1949220LL + 1598625LL*Power(r, 2LL)*Power(xj, 2LL) -

                               41391LL*Power(r, 4LL)*Power(xj, 4LL) + 128LL*Power(r, 6LL)*Power(xj, 6LL)) +

        r*Power(xi, 13LL)*Power(xj, 8LL)*

        (173037375LL - 2784600LL*Power(r, 2LL)*Power(xj, 2LL) -

            112140LL*Power(r, 4LL)*Power(xj, 4LL) + 256LL*Power(r, 6LL)*Power(xj, 6LL)) +

        14LL*Power(xi, 12LL)*Power(xj, 8LL)*

        (-7260750LL - 2521935LL*Power(r, 2LL)*Power(xj, 2LL) +

            19500LL*Power(r, 4LL)*Power(xj, 4LL) + 344LL*Power(r, 6LL)*Power(xj, 6LL))) +

                 exp(2LL*r*xi)*Power(xi, 10LL)*

                 (210LL*Power(xi, 2LL)*Power(xj, 18LL)*

        (332010LL*xj + 189000LL*r*Power(xj, 2LL) +

            45675LL*Power(r, 2LL)*Power(xj, 3LL) + 5880LL*Power(r, 3LL)*Power(xj, 4LL) +

            405LL*Power(r, 4LL)*Power(xj, 5LL) + 12LL*Power(r, 5LL)*Power(xj, 6LL)) +

        105LL*Power(xi, 18LL)*Power(xj, 2LL)*

        (315LL*xj + 540LL*r*Power(xj, 2LL) + 450LL*Power(r, 2LL)*Power(xj, 3LL) +

            240LL*Power(r, 3LL)*Power(xj, 4LL) + 90LL*Power(r, 4LL)*Power(xj, 5LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) -

        1365LL*Power(xi, 10LL)*Power(xj, 10LL)*

        (15903LL*xj - 51732LL*r*Power(xj, 2LL) - 6120LL*Power(r, 2LL)*Power(xj, 3LL) +

            4320LL*Power(r, 3LL)*Power(xj, 4LL) + 900LL*Power(r, 4LL)*Power(xj, 5LL) +

            48LL*Power(r, 5LL)*Power(xj, 6LL)) +

        Power(xi, 14LL)*Power(xj, 6LL)*

        (1003275LL*xj + 1719900LL*r*Power(xj, 2LL) +

            1162980LL*Power(r, 2LL)*Power(xj, 3LL) +

            1485120LL*Power(r, 3LL)*Power(xj, 4LL) + 57960LL*Power(r, 4LL)*Power(xj, 5LL) -

            28896LL*Power(r, 5LL)*Power(xj, 6LL) - 1792LL*Power(r, 6LL)*Power(xj, 7LL)) +

        2LL*Power(xj, 20LL)*(1949220LL*xj + 1375920LL*r*Power(xj, 2LL) +

                             429975LL*Power(r, 2LL)*Power(xj, 3LL) + 76440LL*Power(r, 3LL)*Power(xj, 4LL) +

                             8190LL*Power(r, 4LL)*Power(xj, 5LL) + 504LL*Power(r, 5LL)*Power(xj, 6LL) +

                             14LL*Power(r, 6LL)*Power(xj, 7LL)) -

        42LL*Power(xi, 4LL)*Power(xj, 16LL)*

        (-4251870LL*xj - 986040LL*r*Power(xj, 2LL) +

            126765LL*Power(r, 2LL)*Power(xj, 3LL) + 69960LL*Power(r, 3LL)*Power(xj, 4LL) +

            9855LL*Power(r, 4LL)*Power(xj, 5LL) + 612LL*Power(r, 5LL)*Power(xj, 6LL) +

            14LL*Power(r, 6LL)*Power(xj, 7LL)) +

        21LL*Power(xi, 16LL)*Power(xj, 4LL)*

        (-11025LL*xj - 18900LL*r*Power(xj, 2LL) - 15750LL*Power(r, 2LL)*Power(xj, 3LL) -

            8400LL*Power(r, 3LL)*Power(xj, 4LL) - 4140LL*Power(r, 4LL)*Power(xj, 5LL) -

            48LL*Power(r, 5LL)*Power(xj, 6LL) + 28LL*Power(r, 6LL)*Power(xj, 7LL)) -

        Power(xi, 20LL)*(2205LL*xj + 3780LL*r*Power(xj, 2LL) +

                         3150LL*Power(r, 2LL)*Power(xj, 3LL) + 1680LL*Power(r, 3LL)*Power(xj, 4LL) +

                         630LL*Power(r, 4LL)*Power(xj, 5LL) + 168LL*Power(r, 5LL)*Power(xj, 6LL) +

                         28LL*Power(r, 6LL)*Power(xj, 7LL)) -

        35LL*Power(xi, 8LL)*Power(xj, 12LL)*

        (4943925LL*xj + 517860LL*r*Power(xj, 2LL) -

            1078560LL*Power(r, 2LL)*Power(xj, 3LL) -

            281760LL*Power(r, 3LL)*Power(xj, 4LL) - 20880LL*Power(r, 4LL)*Power(xj, 5LL) +

            192LL*Power(r, 5LL)*Power(xj, 6LL) + 56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        35LL*Power(xi, 12LL)*Power(xj, 8LL)*

        (-98865LL*xj + 7020LL*r*Power(xj, 2LL) - 393120LL*Power(r, 2LL)*Power(xj, 3LL) -

            31200LL*Power(r, 3LL)*Power(xj, 4LL) + 16020LL*Power(r, 4LL)*Power(xj, 5LL) +

            2160LL*Power(r, 5LL)*Power(xj, 6LL) + 56LL*Power(r, 6LL)*Power(xj, 7LL)) +

        Power(xi, 6LL)*Power(xj, 14LL)*

        (-54796455LL*xj - 137966220LL*r*Power(xj, 2LL) -

            38348100LL*Power(r, 2LL)*Power(xj, 3LL) -

            2654400LL*Power(r, 3LL)*Power(xj, 4LL) + 269640LL*Power(r, 4LL)*Power(xj, 5LL) +

            46368LL*Power(r, 5LL)*Power(xj, 6LL) + 1792LL*Power(r, 6LL)*Power(xj, 7LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 11LL)*

                 (210LL*Power(xi, 2LL)*Power(xj, 18LL)*

        (514080LL + 332010LL*r*xj + 94500LL*Power(r, 2LL)*Power(xj, 2LL) +

            15225LL*Power(r, 3LL)*Power(xj, 3LL) + 1470LL*Power(r, 4LL)*Power(xj, 4LL) +

            81LL*Power(r, 5LL)*Power(xj, 5LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

        105LL*Power(xi, 18LL)*Power(xj, 2LL)*

        (180LL + 315LL*r*xj + 270LL*Power(r, 2LL)*Power(xj, 2LL) +

            150LL*Power(r, 3LL)*Power(xj, 3LL) + 60LL*Power(r, 4LL)*Power(xj, 4LL) +

            18LL*Power(r, 5LL)*Power(xj, 5LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

        1365LL*Power(xi, 10LL)*Power(xj, 10LL)*

        (-6444LL + 15903LL*r*xj - 25866LL*Power(r, 2LL)*Power(xj, 2LL) -

            2040LL*Power(r, 3LL)*Power(xj, 3LL) + 1080LL*Power(r, 4LL)*Power(xj, 4LL) +

            180LL*Power(r, 5LL)*Power(xj, 5LL) + 8LL*Power(r, 6LL)*Power(xj, 6LL)) +

        Power(xi, 14LL)*Power(xj, 6LL)*

        (573300LL + 1003275LL*r*xj + 859950LL*Power(r, 2LL)*Power(xj, 2LL) +

            387660LL*Power(r, 3LL)*Power(xj, 3LL) + 371280LL*Power(r, 4LL)*Power(xj, 4LL) +

            11592LL*Power(r, 5LL)*Power(xj, 5LL) - 4816LL*Power(r, 6LL)*Power(xj, 6LL) -

            256LL*Power(r, 7LL)*Power(xj, 7LL)) +

        2LL*Power(xj, 20LL)*(2506140LL + 1949220LL*r*xj +

                             687960LL*Power(r, 2LL)*Power(xj, 2LL) + 143325LL*Power(r, 3LL)*Power(xj, 3LL) +

                             19110LL*Power(r, 4LL)*Power(xj, 4LL) + 1638LL*Power(r, 5LL)*Power(xj, 5LL) +

                             84LL*Power(r, 6LL)*Power(xj, 6LL) + 2LL*Power(r, 7LL)*Power(xj, 7LL)) -

        42LL*Power(xi, 4LL)*Power(xj, 16LL)*

        (-10437660LL - 4251870LL*r*xj - 493020LL*Power(r, 2LL)*Power(xj, 2LL) +

            42255LL*Power(r, 3LL)*Power(xj, 3LL) + 17490LL*Power(r, 4LL)*Power(xj, 4LL) +

            1971LL*Power(r, 5LL)*Power(xj, 5LL) + 102LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 7LL)*Power(xj, 7LL)) +

        21LL*Power(xi, 16LL)*Power(xj, 4LL)*

        (-6300LL - 11025LL*r*xj - 9450LL*Power(r, 2LL)*Power(xj, 2LL) -

            5250LL*Power(r, 3LL)*Power(xj, 3LL) - 2100LL*Power(r, 4LL)*Power(xj, 4LL) -

            828LL*Power(r, 5LL)*Power(xj, 5LL) - 8LL*Power(r, 6LL)*Power(xj, 6LL) +

            4LL*Power(r, 7LL)*Power(xj, 7LL)) -

        Power(xi, 20LL)*(1260LL + 2205LL*r*xj + 1890LL*Power(r, 2LL)*Power(xj, 2LL) +

                         1050LL*Power(r, 3LL)*Power(xj, 3LL) + 420LL*Power(r, 4LL)*Power(xj, 4LL) +

                         126LL*Power(r, 5LL)*Power(xj, 5LL) + 28LL*Power(r, 6LL)*Power(xj, 6LL) +

                         4LL*Power(r, 7LL)*Power(xj, 7LL)) -

        35LL*Power(xi, 8LL)*Power(xj, 12LL)*

        (-2904300LL + 4943925LL*r*xj + 258930LL*Power(r, 2LL)*Power(xj, 2LL) -

            359520LL*Power(r, 3LL)*Power(xj, 3LL) - 70440LL*Power(r, 4LL)*Power(xj, 4LL) -

            4176LL*Power(r, 5LL)*Power(xj, 5LL) + 32LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        35LL*Power(xi, 12LL)*Power(xj, 8LL)*

        (-49140LL - 98865LL*r*xj + 3510LL*Power(r, 2LL)*Power(xj, 2LL) -

            131040LL*Power(r, 3LL)*Power(xj, 3LL) - 7800LL*Power(r, 4LL)*Power(xj, 4LL) +

            3204LL*Power(r, 5LL)*Power(xj, 5LL) + 360LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 7LL)*Power(xj, 7LL)) +

        Power(xi, 6LL)*Power(xj, 14LL)*

        (446489820LL - 54796455LL*r*xj - 68983110LL*Power(r, 2LL)*Power(xj, 2LL) -

            12782700LL*Power(r, 3LL)*Power(xj, 3LL) - 663600LL*Power(r, 4LL)*Power(xj, 4LL) +

            53928LL*Power(r, 5LL)*Power(xj, 5LL) + 7728LL*Power(r, 6LL)*Power(xj, 6LL) +

            256LL*Power(r, 7LL)*Power(xj, 7LL))))/

                (1260LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 15LL)*Power(xi + xj, 15LL))

            ;
        }

    }
    return S;
}
