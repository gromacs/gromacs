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

cl_R DSlater_4S_5S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-25913502934444125LL*xi + 28454994247680000LL*exp(2LL*r*xi)*xi -

                  46744023242416500LL*r*Power(xi, 2LL) -

                  41723129607909750LL*Power(r, 2LL)*Power(xi, 3LL) -

                  24550942638222000LL*Power(r, 3LL)*Power(xi, 4LL) -

                  10704286944351000LL*Power(r, 4LL)*Power(xi, 5LL) -

                  3684699450432000LL*Power(r, 5LL)*Power(xi, 6LL) -

                  1041667066440000LL*Power(r, 6LL)*Power(xi, 7LL) -

                  248293113868800LL*Power(r, 7LL)*Power(xi, 8LL) -

                  50808078921600LL*Power(r, 8LL)*Power(xi, 9LL) -

                  9033331507200LL*Power(r, 9LL)*Power(xi, 10LL) -

                  1405184901120LL*Power(r, 10LL)*Power(xi, 11LL) -

                  191616122880LL*Power(r, 11LL)*Power(xi, 12LL) -

                  22811443200LL*Power(r, 12LL)*Power(xi, 13LL) -

                  2339635200LL*Power(r, 13LL)*Power(xi, 14LL) -

                  200540160LL*Power(r, 14LL)*Power(xi, 15LL) -

                  13369344LL*Power(r, 15LL)*Power(xi, 16LL) - 557056LL*Power(r, 16LL)*Power(xi, 17LL))/

                (1.422749712384e16*exp(2LL*r*xi)*r) +

                (-14227497123840000LL + 14227497123840000LL*exp(2LL*r*xi) -

                 25913502934444125LL*r*xi - 23372011621208250LL*Power(r, 2LL)*Power(xi, 2LL) -

                 13907709869303250LL*Power(r, 3LL)*Power(xi, 3LL) -

                 6137735659555500LL*Power(r, 4LL)*Power(xi, 4LL) -

                 2140857388870200LL*Power(r, 5LL)*Power(xi, 5LL) -

                 614116575072000LL*Power(r, 6LL)*Power(xi, 6LL) -

                 148809580920000LL*Power(r, 7LL)*Power(xi, 7LL) -

                 31036639233600LL*Power(r, 8LL)*Power(xi, 8LL) -

                 5645342102400LL*Power(r, 9LL)*Power(xi, 9LL) -

                 903333150720LL*Power(r, 10LL)*Power(xi, 10LL) -

                 127744081920LL*Power(r, 11LL)*Power(xi, 11LL) -

                 15968010240LL*Power(r, 12LL)*Power(xi, 12LL) -

                 1754726400LL*Power(r, 13LL)*Power(xi, 13LL) -

                 167116800LL*Power(r, 14LL)*Power(xi, 14LL) -

                 13369344LL*Power(r, 15LL)*Power(xi, 15LL) - 835584LL*Power(r, 16LL)*Power(xi, 16LL) -

                 32768LL*Power(r, 17LL)*Power(xi, 17LL))/

                (1.422749712384e16*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-14227497123840000LL + 14227497123840000LL*exp(2LL*r*xi) -

                     25913502934444125LL*r*xi - 23372011621208250LL*Power(r, 2LL)*Power(xi, 2LL) -

                     13907709869303250LL*Power(r, 3LL)*Power(xi, 3LL) -

                     6137735659555500LL*Power(r, 4LL)*Power(xi, 4LL) -

                     2140857388870200LL*Power(r, 5LL)*Power(xi, 5LL) -

                     614116575072000LL*Power(r, 6LL)*Power(xi, 6LL) -

                     148809580920000LL*Power(r, 7LL)*Power(xi, 7LL) -

                     31036639233600LL*Power(r, 8LL)*Power(xi, 8LL) -

                     5645342102400LL*Power(r, 9LL)*Power(xi, 9LL) -

                     903333150720LL*Power(r, 10LL)*Power(xi, 10LL) -

                     127744081920LL*Power(r, 11LL)*Power(xi, 11LL) -

                     15968010240LL*Power(r, 12LL)*Power(xi, 12LL) -

                     1754726400LL*Power(r, 13LL)*Power(xi, 13LL) -

                     167116800LL*Power(r, 14LL)*Power(xi, 14LL) -

                     13369344LL*Power(r, 15LL)*Power(xi, 15LL) - 835584LL*Power(r, 16LL)*Power(xi, 16LL) -

                     32768LL*Power(r, 17LL)*Power(xi, 17LL)))/(7.11374856192e15*exp(2LL*r*xi)*r)

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
            S = (56700LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 17LL) +

                 9LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-980LL*Power(r, 6LL)*Power(xi, 28LL) - 20LL*Power(r, 7LL)*Power(xi, 29LL) +

                  6300LL*Power(xj, 22LL) + 11025LL*r*xi*Power(xj, 22LL) -

                  50LL*Power(r, 5LL)*Power(xi, 27LL)*(441LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  3150LL*Power(xi, 2LL)*Power(xj, 20LL)*(-34LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  525LL*r*Power(xi, 3LL)*Power(xj, 20LL)*

                  (-357LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  420LL*Power(r, 4LL)*Power(xi, 26LL)*(700LL + 19LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  1050LL*Power(xi, 4LL)*Power(xj, 18LL)*

                  (816LL - 153LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  210LL*r*Power(xi, 5LL)*Power(xj, 18LL)*

                  (7140LL - 425LL*Power(r, 2LL)*Power(xj, 2LL) + 3LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  42LL*Power(r, 3LL)*Power(xi, 25LL)*

                  (-59500LL - 6035LL*Power(r, 2LL)*Power(xj, 2LL) +

            18LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  84LL*Power(r, 2LL)*Power(xi, 24LL)*

                  (-160650LL - 52700LL*Power(r, 2LL)*Power(xj, 2LL) +

            397LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  28LL*Power(xi, 12LL)*Power(xj, 10LL)*

                  (100849950LL + 27100125LL*Power(r, 2LL)*Power(xj, 2LL) +

            186150LL*Power(r, 4LL)*Power(xj, 4LL) - 2177LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  140LL*Power(xi, 6LL)*Power(xj, 16LL)*

                  (-30600LL + 9180LL*Power(r, 2LL)*Power(xj, 2LL) -

            255LL*Power(r, 4LL)*Power(xj, 4LL) + Power(r, 6LL)*Power(xj, 6LL)) -

                  2380LL*Power(xi, 8LL)*Power(xj, 14LL)*

                  (-6300LL + 2700LL*Power(r, 2LL)*Power(xj, 2LL) -

            120LL*Power(r, 4LL)*Power(xj, 4LL) + Power(r, 6LL)*Power(xj, 6LL)) +

                  10LL*r*Power(xi, 7LL)*Power(xj, 16LL)*

                  (-749700LL + 71400LL*Power(r, 2LL)*Power(xj, 2LL) -

            1071LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  204LL*r*Power(xi, 15LL)*Power(xj, 8LL)*

                  (28962255LL - 1744750LL*Power(r, 2LL)*Power(xj, 2LL) +

            9555LL*Power(r, 4LL)*Power(xj, 4LL) + 6LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  42LL*r*Power(xi, 11LL)*Power(xj, 12LL)*

                  (-12911925LL - 1634550LL*Power(r, 2LL)*Power(xj, 2LL) -

            7103LL*Power(r, 4LL)*Power(xj, 4LL) + 18LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  2LL*r*Power(xi, 9LL)*Power(xj, 14LL)*

                  (16948575LL - 1184400LL*Power(r, 2LL)*Power(xj, 2LL) +

            63861LL*Power(r, 4LL)*Power(xj, 4LL) + 50LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  28LL*Power(xi, 22LL)*(-2180250LL - 10993050LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        14925LL*Power(r, 4LL)*Power(xj, 4LL) + 73LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  952LL*Power(xi, 14LL)*Power(xj, 8LL)*

                  (16966215LL + 725175LL*Power(r, 2LL)*Power(xj, 2LL) -

            36075LL*Power(r, 4LL)*Power(xj, 4LL) + 79LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  84LL*Power(xi, 10LL)*Power(xj, 12LL)*

                  (1723800LL + 279225LL*Power(r, 2LL)*Power(xj, 2LL) +

            45600LL*Power(r, 4LL)*Power(xj, 4LL) + 107LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  35LL*r*Power(xi, 17LL)*Power(xj, 6LL)*

                  (132637869LL - 2205240LL*Power(r, 2LL)*Power(xj, 2LL) -

            48348LL*Power(r, 4LL)*Power(xj, 4LL) + 136LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  6LL*r*Power(xi, 21LL)*Power(xj, 2LL)*

                  (192298050LL + 12644275LL*Power(r, 2LL)*Power(xj, 2LL) -

            218029LL*Power(r, 4LL)*Power(xj, 4LL) + 204LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  4LL*r*Power(xi, 13LL)*Power(xj, 10LL)*

                  (1259522775LL + 15895425LL*Power(r, 2LL)*Power(xj, 2LL) -

            493017LL*Power(r, 4LL)*Power(xj, 4LL) + 263LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  140LL*Power(xi, 16LL)*Power(xj, 6LL)*

                  (180826281LL - 15101406LL*Power(r, 2LL)*Power(xj, 2LL) +

            160140LL*Power(r, 4LL)*Power(xj, 4LL) + 442LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  2LL*r*Power(xi, 23LL)*(21366450LL + 23526300LL*Power(r, 2LL)*Power(xj, 2LL) -

                                         246729LL*Power(r, 4LL)*Power(xj, 4LL) + 526LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  7LL*r*Power(xi, 19LL)*Power(xj, 4LL)*

                  (-811081215LL + 39095550LL*Power(r, 2LL)*Power(xj, 2LL) -

            515916LL*Power(r, 4LL)*Power(xj, 4LL) + 680LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  70LL*Power(xi, 18LL)*Power(xj, 4LL)*

                  (-180554454LL + 9873711LL*Power(r, 2LL)*Power(xj, 2LL) -

            414120LL*Power(r, 4LL)*Power(xj, 4LL) + 2924LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  14LL*Power(xi, 20LL)*Power(xj, 2LL)*

                  (136919700LL + 71867115LL*Power(r, 2LL)*Power(xj, 2LL) -

            2154150LL*Power(r, 4LL)*Power(xj, 4LL) + 10268LL*Power(r, 6LL)*Power(xj, 6LL)))

                 - 4LL*exp(2LL*r*xi)*Power(xi, 10LL)*

                 (-10710LL*Power(xi, 12LL)*Power(xj, 12LL)*

                  (-3555LL - 127008LL*r*xj + 138384LL*Power(r, 2LL)*Power(xj, 2LL) -

            74556LL*Power(r, 3LL)*Power(xj, 3LL) - 22284LL*Power(r, 4LL)*Power(xj, 4LL) +

            408LL*Power(r, 5LL)*Power(xj, 5LL) + 576LL*Power(r, 6LL)*Power(xj, 6LL) +

            60LL*Power(r, 7LL)*Power(xj, 7LL) + 2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  2LL*Power(xi, 20LL)*Power(xj, 4LL)*

                  (963900LL + 1735020LL*r*xj + 1542240LL*Power(r, 2LL)*Power(xj, 2LL) +

            899640LL*Power(r, 3LL)*Power(xj, 3LL) + 385560LL*Power(r, 4LL)*Power(xj, 4LL) +

            128520LL*Power(r, 5LL)*Power(xj, 5LL) + 34272LL*Power(r, 6LL)*Power(xj, 6LL) +

            9126LL*Power(r, 7LL)*Power(xj, 7LL) + 333LL*Power(r, 8LL)*Power(xj, 8LL) -

            20LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  2LL*Power(xj, 24LL)*(119041650LL + 107137485LL*r*xj +

                                       45110520LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       11695320LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       2063880LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       257985LL*Power(r, 5LL)*Power(xj, 5LL) + 22932LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       1404LL*Power(r, 7LL)*Power(xj, 7LL) + 54LL*Power(r, 8LL)*Power(xj, 8LL) +

                                       Power(r, 9LL)*Power(xj, 9LL)) +

                  2LL*Power(xi, 2LL)*Power(xj, 22LL)*

                  (-3264488325LL - 2505368880LL*r*xj -

            881390160LL*Power(r, 2LL)*Power(xj, 2LL) -

            185775660LL*Power(r, 3LL)*Power(xj, 3LL) -

            25639740LL*Power(r, 4LL)*Power(xj, 4LL) -

            2361555LL*Power(r, 5LL)*Power(xj, 5LL) -

            139356LL*Power(r, 6LL)*Power(xj, 6LL) - 4482LL*Power(r, 7LL)*Power(xj, 7LL) -

            27LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  Power(xi, 24LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                                   2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  102LL*Power(xi, 10LL)*Power(xj, 14LL)*

                  (44986725LL - 97433280LL*r*xj + 44467920LL*Power(r, 2LL)*Power(xj, 2LL) +

            15857100LL*Power(r, 3LL)*Power(xj, 3LL) -

            457380LL*Power(r, 4LL)*Power(xj, 4LL) - 620550LL*Power(r, 5LL)*Power(xj, 5LL) -

            83160LL*Power(r, 6LL)*Power(xj, 6LL) - 4068LL*Power(r, 7LL)*Power(xj, 7LL) -

            6LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  102LL*Power(xi, 14LL)*Power(xj, 10LL)*

                  (-859950LL - 1437345LL*r*xj - 2260440LL*Power(r, 2LL)*Power(xj, 2LL) +

            810810LL*Power(r, 3LL)*Power(xj, 3LL) -

            1056510LL*Power(r, 4LL)*Power(xj, 4LL) -

            217854LL*Power(r, 5LL)*Power(xj, 5LL) + 6552LL*Power(r, 6LL)*Power(xj, 6LL) +

            3852LL*Power(r, 7LL)*Power(xj, 7LL) + 258LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  Power(xi, 22LL)*Power(xj, 2LL)*

                  (240975LL + 433755LL*r*xj + 385560LL*Power(r, 2LL)*Power(xj, 2LL) +

            224910LL*Power(r, 3LL)*Power(xj, 3LL) + 96390LL*Power(r, 4LL)*Power(xj, 4LL) +

            32130LL*Power(r, 5LL)*Power(xj, 5LL) + 8568LL*Power(r, 6LL)*Power(xj, 6LL) +

            1836LL*Power(r, 7LL)*Power(xj, 7LL) + 306LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  2LL*Power(xi, 4LL)*Power(xj, 20LL)*

                  (-18032978565LL - 9823683240LL*r*xj -

            2047323600LL*Power(r, 2LL)*Power(xj, 2LL) -

            129098340LL*Power(r, 3LL)*Power(xj, 3LL) +

            26410860LL*Power(r, 4LL)*Power(xj, 4LL) +

            7094304LL*Power(r, 5LL)*Power(xj, 5LL) +

            788256LL*Power(r, 6LL)*Power(xj, 6LL) + 48654LL*Power(r, 7LL)*Power(xj, 7LL) +

            1593LL*Power(r, 8LL)*Power(xj, 8LL) + 20LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  6LL*Power(xi, 16LL)*Power(xj, 8LL)*

                  (-5622750LL - 10120950LL*r*xj - 8996400LL*Power(r, 2LL)*Power(xj, 2LL) -

            5698350LL*Power(r, 3LL)*Power(xj, 3LL) -

            897750LL*Power(r, 4LL)*Power(xj, 4LL) -

            1641591LL*Power(r, 5LL)*Power(xj, 5LL) -

            211932LL*Power(r, 6LL)*Power(xj, 6LL) + 10224LL*Power(r, 7LL)*Power(xj, 7LL) +

            2364LL*Power(r, 8LL)*Power(xj, 8LL) + 73LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  2LL*Power(xi, 18LL)*Power(xj, 6LL)*

                  (-4819500LL - 8675100LL*r*xj - 7711200LL*Power(r, 2LL)*Power(xj, 2LL) -

            4498200LL*Power(r, 3LL)*Power(xj, 3LL) -

            1927800LL*Power(r, 4LL)*Power(xj, 4LL) -

            561519LL*Power(r, 5LL)*Power(xj, 5LL) - 279468LL*Power(r, 6LL)*Power(xj, 6LL) -

            20682LL*Power(r, 7LL)*Power(xj, 7LL) + 1305LL*Power(r, 8LL)*Power(xj, 8LL) +

            106LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  3LL*Power(xi, 8LL)*Power(xj, 16LL)*

                  (-9364244085LL + 6940428705LL*r*xj +

            2117684520LL*Power(r, 2LL)*Power(xj, 2LL) -

            230268150LL*Power(r, 3LL)*Power(xj, 3LL) -

            149610510LL*Power(r, 4LL)*Power(xj, 4LL) -

            21824334LL*Power(r, 5LL)*Power(xj, 5LL) -

            1223208LL*Power(r, 6LL)*Power(xj, 6LL) + 12708LL*Power(r, 7LL)*Power(xj, 7LL) +

            4470LL*Power(r, 8LL)*Power(xj, 8LL) + 146LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  Power(xi, 6LL)*Power(xj, 18LL)*

                  (57304872765LL + 7147185255LL*r*xj -

            5801702760LL*Power(r, 2LL)*Power(xj, 2LL) -

            2053388610LL*Power(r, 3LL)*Power(xj, 3LL) -

            271655370LL*Power(r, 4LL)*Power(xj, 4LL) -

            10864854LL*Power(r, 5LL)*Power(xj, 5LL) +

            1337112LL*Power(r, 6LL)*Power(xj, 6LL) + 202716LL*Power(r, 7LL)*Power(xj, 7LL) +

            10746LL*Power(r, 8LL)*Power(xj, 8LL) + 212LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (56700LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 17LL)*

                 Power(xi + xj, 17LL)) + (56700LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 17LL) +

                                          9LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                                          (-980LL*Power(r, 6LL)*Power(xi, 28LL) - 20LL*Power(r, 7LL)*Power(xi, 29LL) +

                                  6300LL*Power(xj, 22LL) + 11025LL*r*xi*Power(xj, 22LL) -

                                  50LL*Power(r, 5LL)*Power(xi, 27LL)*(441LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  3150LL*Power(xi, 2LL)*Power(xj, 20LL)*(-34LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  525LL*r*Power(xi, 3LL)*Power(xj, 20LL)*

                                  (-357LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  420LL*Power(r, 4LL)*Power(xi, 26LL)*(700LL + 19LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  1050LL*Power(xi, 4LL)*Power(xj, 18LL)*

                                  (816LL - 153LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  210LL*r*Power(xi, 5LL)*Power(xj, 18LL)*

                                  (7140LL - 425LL*Power(r, 2LL)*Power(xj, 2LL) + 3LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  42LL*Power(r, 3LL)*Power(xi, 25LL)*

                                  (-59500LL - 6035LL*Power(r, 2LL)*Power(xj, 2LL) +

            18LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  84LL*Power(r, 2LL)*Power(xi, 24LL)*

                                  (-160650LL - 52700LL*Power(r, 2LL)*Power(xj, 2LL) +

            397LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  28LL*Power(xi, 12LL)*Power(xj, 10LL)*

                                  (100849950LL + 27100125LL*Power(r, 2LL)*Power(xj, 2LL) +

            186150LL*Power(r, 4LL)*Power(xj, 4LL) - 2177LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  140LL*Power(xi, 6LL)*Power(xj, 16LL)*

                                  (-30600LL + 9180LL*Power(r, 2LL)*Power(xj, 2LL) -

            255LL*Power(r, 4LL)*Power(xj, 4LL) + Power(r, 6LL)*Power(xj, 6LL)) -

                                  2380LL*Power(xi, 8LL)*Power(xj, 14LL)*

                                  (-6300LL + 2700LL*Power(r, 2LL)*Power(xj, 2LL) -

            120LL*Power(r, 4LL)*Power(xj, 4LL) + Power(r, 6LL)*Power(xj, 6LL)) +

                                  10LL*r*Power(xi, 7LL)*Power(xj, 16LL)*

                                  (-749700LL + 71400LL*Power(r, 2LL)*Power(xj, 2LL) -

            1071LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  204LL*r*Power(xi, 15LL)*Power(xj, 8LL)*

                                  (28962255LL - 1744750LL*Power(r, 2LL)*Power(xj, 2LL) +

            9555LL*Power(r, 4LL)*Power(xj, 4LL) + 6LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  42LL*r*Power(xi, 11LL)*Power(xj, 12LL)*

                                  (-12911925LL - 1634550LL*Power(r, 2LL)*Power(xj, 2LL) -

            7103LL*Power(r, 4LL)*Power(xj, 4LL) + 18LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  2LL*r*Power(xi, 9LL)*Power(xj, 14LL)*

                                  (16948575LL - 1184400LL*Power(r, 2LL)*Power(xj, 2LL) +

            63861LL*Power(r, 4LL)*Power(xj, 4LL) + 50LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  28LL*Power(xi, 22LL)*(-2180250LL - 10993050LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        14925LL*Power(r, 4LL)*Power(xj, 4LL) + 73LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  952LL*Power(xi, 14LL)*Power(xj, 8LL)*

                                  (16966215LL + 725175LL*Power(r, 2LL)*Power(xj, 2LL) -

            36075LL*Power(r, 4LL)*Power(xj, 4LL) + 79LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  84LL*Power(xi, 10LL)*Power(xj, 12LL)*

                                  (1723800LL + 279225LL*Power(r, 2LL)*Power(xj, 2LL) +

            45600LL*Power(r, 4LL)*Power(xj, 4LL) + 107LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  35LL*r*Power(xi, 17LL)*Power(xj, 6LL)*

                                  (132637869LL - 2205240LL*Power(r, 2LL)*Power(xj, 2LL) -

            48348LL*Power(r, 4LL)*Power(xj, 4LL) + 136LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  6LL*r*Power(xi, 21LL)*Power(xj, 2LL)*

                                  (192298050LL + 12644275LL*Power(r, 2LL)*Power(xj, 2LL) -

            218029LL*Power(r, 4LL)*Power(xj, 4LL) + 204LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  4LL*r*Power(xi, 13LL)*Power(xj, 10LL)*

                                  (1259522775LL + 15895425LL*Power(r, 2LL)*Power(xj, 2LL) -

            493017LL*Power(r, 4LL)*Power(xj, 4LL) + 263LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  140LL*Power(xi, 16LL)*Power(xj, 6LL)*

                                  (180826281LL - 15101406LL*Power(r, 2LL)*Power(xj, 2LL) +

            160140LL*Power(r, 4LL)*Power(xj, 4LL) + 442LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  2LL*r*Power(xi, 23LL)*(21366450LL + 23526300LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                         246729LL*Power(r, 4LL)*Power(xj, 4LL) + 526LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  7LL*r*Power(xi, 19LL)*Power(xj, 4LL)*

                                  (-811081215LL + 39095550LL*Power(r, 2LL)*Power(xj, 2LL) -

            515916LL*Power(r, 4LL)*Power(xj, 4LL) + 680LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  70LL*Power(xi, 18LL)*Power(xj, 4LL)*

                                  (-180554454LL + 9873711LL*Power(r, 2LL)*Power(xj, 2LL) -

            414120LL*Power(r, 4LL)*Power(xj, 4LL) + 2924LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  14LL*Power(xi, 20LL)*Power(xj, 2LL)*

                                  (136919700LL + 71867115LL*Power(r, 2LL)*Power(xj, 2LL) -

            2154150LL*Power(r, 4LL)*Power(xj, 4LL) + 10268LL*Power(r, 6LL)*Power(xj, 6LL)))

                                          - 4LL*exp(2LL*r*xi)*Power(xi, 10LL)*

                                          (-10710LL*Power(xi, 12LL)*Power(xj, 12LL)*

                                  (-3555LL - 127008LL*r*xj + 138384LL*Power(r, 2LL)*Power(xj, 2LL) -

            74556LL*Power(r, 3LL)*Power(xj, 3LL) - 22284LL*Power(r, 4LL)*Power(xj, 4LL) +

            408LL*Power(r, 5LL)*Power(xj, 5LL) + 576LL*Power(r, 6LL)*Power(xj, 6LL) +

            60LL*Power(r, 7LL)*Power(xj, 7LL) + 2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  2LL*Power(xi, 20LL)*Power(xj, 4LL)*

                                  (963900LL + 1735020LL*r*xj + 1542240LL*Power(r, 2LL)*Power(xj, 2LL) +

            899640LL*Power(r, 3LL)*Power(xj, 3LL) + 385560LL*Power(r, 4LL)*Power(xj, 4LL) +

            128520LL*Power(r, 5LL)*Power(xj, 5LL) + 34272LL*Power(r, 6LL)*Power(xj, 6LL) +

            9126LL*Power(r, 7LL)*Power(xj, 7LL) + 333LL*Power(r, 8LL)*Power(xj, 8LL) -

            20LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  2LL*Power(xj, 24LL)*(119041650LL + 107137485LL*r*xj +

                                                       45110520LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       11695320LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                       2063880LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                       257985LL*Power(r, 5LL)*Power(xj, 5LL) + 22932LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                       1404LL*Power(r, 7LL)*Power(xj, 7LL) + 54LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                       Power(r, 9LL)*Power(xj, 9LL)) +

                                  2LL*Power(xi, 2LL)*Power(xj, 22LL)*

                                  (-3264488325LL - 2505368880LL*r*xj -

            881390160LL*Power(r, 2LL)*Power(xj, 2LL) -

            185775660LL*Power(r, 3LL)*Power(xj, 3LL) -

            25639740LL*Power(r, 4LL)*Power(xj, 4LL) -

            2361555LL*Power(r, 5LL)*Power(xj, 5LL) -

            139356LL*Power(r, 6LL)*Power(xj, 6LL) - 4482LL*Power(r, 7LL)*Power(xj, 7LL) -

            27LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  Power(xi, 24LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                   2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  102LL*Power(xi, 10LL)*Power(xj, 14LL)*

                                  (44986725LL - 97433280LL*r*xj + 44467920LL*Power(r, 2LL)*Power(xj, 2LL) +

            15857100LL*Power(r, 3LL)*Power(xj, 3LL) -

            457380LL*Power(r, 4LL)*Power(xj, 4LL) - 620550LL*Power(r, 5LL)*Power(xj, 5LL) -

            83160LL*Power(r, 6LL)*Power(xj, 6LL) - 4068LL*Power(r, 7LL)*Power(xj, 7LL) -

            6LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  102LL*Power(xi, 14LL)*Power(xj, 10LL)*

                                  (-859950LL - 1437345LL*r*xj - 2260440LL*Power(r, 2LL)*Power(xj, 2LL) +

            810810LL*Power(r, 3LL)*Power(xj, 3LL) -

            1056510LL*Power(r, 4LL)*Power(xj, 4LL) -

            217854LL*Power(r, 5LL)*Power(xj, 5LL) + 6552LL*Power(r, 6LL)*Power(xj, 6LL) +

            3852LL*Power(r, 7LL)*Power(xj, 7LL) + 258LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  Power(xi, 22LL)*Power(xj, 2LL)*

                                  (240975LL + 433755LL*r*xj + 385560LL*Power(r, 2LL)*Power(xj, 2LL) +

            224910LL*Power(r, 3LL)*Power(xj, 3LL) + 96390LL*Power(r, 4LL)*Power(xj, 4LL) +

            32130LL*Power(r, 5LL)*Power(xj, 5LL) + 8568LL*Power(r, 6LL)*Power(xj, 6LL) +

            1836LL*Power(r, 7LL)*Power(xj, 7LL) + 306LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  2LL*Power(xi, 4LL)*Power(xj, 20LL)*

                                  (-18032978565LL - 9823683240LL*r*xj -

            2047323600LL*Power(r, 2LL)*Power(xj, 2LL) -

            129098340LL*Power(r, 3LL)*Power(xj, 3LL) +

            26410860LL*Power(r, 4LL)*Power(xj, 4LL) +

            7094304LL*Power(r, 5LL)*Power(xj, 5LL) +

            788256LL*Power(r, 6LL)*Power(xj, 6LL) + 48654LL*Power(r, 7LL)*Power(xj, 7LL) +

            1593LL*Power(r, 8LL)*Power(xj, 8LL) + 20LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  6LL*Power(xi, 16LL)*Power(xj, 8LL)*

                                  (-5622750LL - 10120950LL*r*xj - 8996400LL*Power(r, 2LL)*Power(xj, 2LL) -

            5698350LL*Power(r, 3LL)*Power(xj, 3LL) -

            897750LL*Power(r, 4LL)*Power(xj, 4LL) -

            1641591LL*Power(r, 5LL)*Power(xj, 5LL) -

            211932LL*Power(r, 6LL)*Power(xj, 6LL) + 10224LL*Power(r, 7LL)*Power(xj, 7LL) +

            2364LL*Power(r, 8LL)*Power(xj, 8LL) + 73LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  2LL*Power(xi, 18LL)*Power(xj, 6LL)*

                                  (-4819500LL - 8675100LL*r*xj - 7711200LL*Power(r, 2LL)*Power(xj, 2LL) -

            4498200LL*Power(r, 3LL)*Power(xj, 3LL) -

            1927800LL*Power(r, 4LL)*Power(xj, 4LL) -

            561519LL*Power(r, 5LL)*Power(xj, 5LL) - 279468LL*Power(r, 6LL)*Power(xj, 6LL) -

            20682LL*Power(r, 7LL)*Power(xj, 7LL) + 1305LL*Power(r, 8LL)*Power(xj, 8LL) +

            106LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  3LL*Power(xi, 8LL)*Power(xj, 16LL)*

                                  (-9364244085LL + 6940428705LL*r*xj +

            2117684520LL*Power(r, 2LL)*Power(xj, 2LL) -

            230268150LL*Power(r, 3LL)*Power(xj, 3LL) -

            149610510LL*Power(r, 4LL)*Power(xj, 4LL) -

            21824334LL*Power(r, 5LL)*Power(xj, 5LL) -

            1223208LL*Power(r, 6LL)*Power(xj, 6LL) + 12708LL*Power(r, 7LL)*Power(xj, 7LL) +

            4470LL*Power(r, 8LL)*Power(xj, 8LL) + 146LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  Power(xi, 6LL)*Power(xj, 18LL)*

                                  (57304872765LL + 7147185255LL*r*xj -

            5801702760LL*Power(r, 2LL)*Power(xj, 2LL) -

            2053388610LL*Power(r, 3LL)*Power(xj, 3LL) -

            271655370LL*Power(r, 4LL)*Power(xj, 4LL) -

            10864854LL*Power(r, 5LL)*Power(xj, 5LL) +

            1337112LL*Power(r, 6LL)*Power(xj, 6LL) + 202716LL*Power(r, 7LL)*Power(xj, 7LL) +

            10746LL*Power(r, 8LL)*Power(xj, 8LL) + 212LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (28350LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 17LL)*Power(xi + xj, 16LL)) -

                (113400LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 17LL) +

                 9LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-5880LL*Power(r, 5LL)*Power(xi, 28LL) - 140LL*Power(r, 6LL)*Power(xi, 29LL) -

        15960LL*Power(r, 5LL)*Power(xi, 26LL)*Power(xj, 2LL) -

        200LL*Power(r, 6LL)*Power(xi, 27LL)*Power(xj, 2LL) + 11025LL*xi*Power(xj, 22LL) +

        18900LL*r*Power(xi, 2LL)*Power(xj, 22LL) +

        10500LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 22LL) -

        250LL*Power(r, 4LL)*Power(xi, 27LL)*(441LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        525LL*Power(xi, 3LL)*Power(xj, 20LL)*(-357LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

        1680LL*Power(r, 3LL)*Power(xi, 26LL)*(700LL + 19LL*Power(r, 2LL)*Power(xj, 2LL)) +

        1050LL*Power(xi, 4LL)*Power(xj, 18LL)*

        (-306LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        210LL*r*Power(xi, 5LL)*Power(xj, 18LL)*

        (-850LL*r*Power(xj, 2LL) + 12LL*Power(r, 3LL)*Power(xj, 4LL)) +

        42LL*Power(r, 3LL)*Power(xi, 25LL)*

        (-12070LL*r*Power(xj, 2LL) + 72LL*Power(r, 3LL)*Power(xj, 4LL)) +

        84LL*Power(r, 2LL)*Power(xi, 24LL)*

        (-105400LL*r*Power(xj, 2LL) + 1588LL*Power(r, 3LL)*Power(xj, 4LL)) +

        210LL*Power(xi, 5LL)*Power(xj, 18LL)*

        (7140LL - 425LL*Power(r, 2LL)*Power(xj, 2LL) + 3LL*Power(r, 4LL)*Power(xj, 4LL)) +

        126LL*Power(r, 2LL)*Power(xi, 25LL)*

        (-59500LL - 6035LL*Power(r, 2LL)*Power(xj, 2LL) + 18LL*Power(r, 4LL)*Power(xj, 4LL))

        + 168LL*r*Power(xi, 24LL)*(-160650LL - 52700LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   397LL*Power(r, 4LL)*Power(xj, 4LL)) -

        28LL*Power(xi, 12LL)*Power(xj, 10LL)*

        (54200250LL*r*Power(xj, 2LL) + 744600LL*Power(r, 3LL)*Power(xj, 4LL) -

            13062LL*Power(r, 5LL)*Power(xj, 6LL)) +

        140LL*Power(xi, 6LL)*Power(xj, 16LL)*

        (18360LL*r*Power(xj, 2LL) - 1020LL*Power(r, 3LL)*Power(xj, 4LL) +

            6LL*Power(r, 5LL)*Power(xj, 6LL)) -

        2380LL*Power(xi, 8LL)*Power(xj, 14LL)*

        (5400LL*r*Power(xj, 2LL) - 480LL*Power(r, 3LL)*Power(xj, 4LL) +

            6LL*Power(r, 5LL)*Power(xj, 6LL)) +

        10LL*r*Power(xi, 7LL)*Power(xj, 16LL)*

        (142800LL*r*Power(xj, 2LL) - 4284LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) +

        204LL*r*Power(xi, 15LL)*Power(xj, 8LL)*

        (-3489500LL*r*Power(xj, 2LL) + 38220LL*Power(r, 3LL)*Power(xj, 4LL) +

            36LL*Power(r, 5LL)*Power(xj, 6LL)) -

        42LL*r*Power(xi, 11LL)*Power(xj, 12LL)*

        (-3269100LL*r*Power(xj, 2LL) - 28412LL*Power(r, 3LL)*Power(xj, 4LL) +

            108LL*Power(r, 5LL)*Power(xj, 6LL)) +

        2LL*r*Power(xi, 9LL)*Power(xj, 14LL)*

        (-2368800LL*r*Power(xj, 2LL) + 255444LL*Power(r, 3LL)*Power(xj, 4LL) +

            300LL*Power(r, 5LL)*Power(xj, 6LL)) +

        28LL*Power(xi, 22LL)*(-21986100LL*r*Power(xj, 2LL) +

                              59700LL*Power(r, 3LL)*Power(xj, 4LL) + 438LL*Power(r, 5LL)*Power(xj, 6LL)) -

        952LL*Power(xi, 14LL)*Power(xj, 8LL)*

        (1450350LL*r*Power(xj, 2LL) - 144300LL*Power(r, 3LL)*Power(xj, 4LL) +

            474LL*Power(r, 5LL)*Power(xj, 6LL)) -

        84LL*Power(xi, 10LL)*Power(xj, 12LL)*

        (558450LL*r*Power(xj, 2LL) + 182400LL*Power(r, 3LL)*Power(xj, 4LL) +

            642LL*Power(r, 5LL)*Power(xj, 6LL)) -

        35LL*r*Power(xi, 17LL)*Power(xj, 6LL)*

        (-4410480LL*r*Power(xj, 2LL) - 193392LL*Power(r, 3LL)*Power(xj, 4LL) +

            816LL*Power(r, 5LL)*Power(xj, 6LL)) -

        6LL*r*Power(xi, 21LL)*Power(xj, 2LL)*

        (25288550LL*r*Power(xj, 2LL) - 872116LL*Power(r, 3LL)*Power(xj, 4LL) +

            1224LL*Power(r, 5LL)*Power(xj, 6LL)) +

        4LL*r*Power(xi, 13LL)*Power(xj, 10LL)*

        (31790850LL*r*Power(xj, 2LL) - 1972068LL*Power(r, 3LL)*Power(xj, 4LL) +

            1578LL*Power(r, 5LL)*Power(xj, 6LL)) -

        140LL*Power(xi, 16LL)*Power(xj, 6LL)*

        (-30202812LL*r*Power(xj, 2LL) + 640560LL*Power(r, 3LL)*Power(xj, 4LL) +

            2652LL*Power(r, 5LL)*Power(xj, 6LL)) -

        2LL*r*Power(xi, 23LL)*(47052600LL*r*Power(xj, 2LL) -

                               986916LL*Power(r, 3LL)*Power(xj, 4LL) + 3156LL*Power(r, 5LL)*Power(xj, 6LL)) +

        7LL*r*Power(xi, 19LL)*Power(xj, 4LL)*

        (78191100LL*r*Power(xj, 2LL) - 2063664LL*Power(r, 3LL)*Power(xj, 4LL) +

            4080LL*Power(r, 5LL)*Power(xj, 6LL)) +

        70LL*Power(xi, 18LL)*Power(xj, 4LL)*

        (19747422LL*r*Power(xj, 2LL) - 1656480LL*Power(r, 3LL)*Power(xj, 4LL) +

            17544LL*Power(r, 5LL)*Power(xj, 6LL)) -

        14LL*Power(xi, 20LL)*Power(xj, 2LL)*

        (143734230LL*r*Power(xj, 2LL) - 8616600LL*Power(r, 3LL)*Power(xj, 4LL) +

            61608LL*Power(r, 5LL)*Power(xj, 6LL)) +

        10LL*Power(xi, 7LL)*Power(xj, 16LL)*

        (-749700LL + 71400LL*Power(r, 2LL)*Power(xj, 2LL) -

            1071LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

        204LL*Power(xi, 15LL)*Power(xj, 8LL)*

        (28962255LL - 1744750LL*Power(r, 2LL)*Power(xj, 2LL) +

            9555LL*Power(r, 4LL)*Power(xj, 4LL) + 6LL*Power(r, 6LL)*Power(xj, 6LL)) -

        42LL*Power(xi, 11LL)*Power(xj, 12LL)*

        (-12911925LL - 1634550LL*Power(r, 2LL)*Power(xj, 2LL) -

            7103LL*Power(r, 4LL)*Power(xj, 4LL) + 18LL*Power(r, 6LL)*Power(xj, 6LL)) +

        2LL*Power(xi, 9LL)*Power(xj, 14LL)*

        (16948575LL - 1184400LL*Power(r, 2LL)*Power(xj, 2LL) +

            63861LL*Power(r, 4LL)*Power(xj, 4LL) + 50LL*Power(r, 6LL)*Power(xj, 6LL)) -

        35LL*Power(xi, 17LL)*Power(xj, 6LL)*

        (132637869LL - 2205240LL*Power(r, 2LL)*Power(xj, 2LL) -

            48348LL*Power(r, 4LL)*Power(xj, 4LL) + 136LL*Power(r, 6LL)*Power(xj, 6LL)) -

        6LL*Power(xi, 21LL)*Power(xj, 2LL)*

        (192298050LL + 12644275LL*Power(r, 2LL)*Power(xj, 2LL) -

            218029LL*Power(r, 4LL)*Power(xj, 4LL) + 204LL*Power(r, 6LL)*Power(xj, 6LL)) +

        4LL*Power(xi, 13LL)*Power(xj, 10LL)*

        (1259522775LL + 15895425LL*Power(r, 2LL)*Power(xj, 2LL) -

            493017LL*Power(r, 4LL)*Power(xj, 4LL) + 263LL*Power(r, 6LL)*Power(xj, 6LL)) -

        2LL*Power(xi, 23LL)*(21366450LL + 23526300LL*Power(r, 2LL)*Power(xj, 2LL) -

                             246729LL*Power(r, 4LL)*Power(xj, 4LL) + 526LL*Power(r, 6LL)*Power(xj, 6LL)) +

        7LL*Power(xi, 19LL)*Power(xj, 4LL)*

        (-811081215LL + 39095550LL*Power(r, 2LL)*Power(xj, 2LL) -

            515916LL*Power(r, 4LL)*Power(xj, 4LL) + 680LL*Power(r, 6LL)*Power(xj, 6LL))) +

                 18LL*exp(2LL*r*xj)*Power(xj, 13LL)*

                 (-980LL*Power(r, 6LL)*Power(xi, 28LL) - 20LL*Power(r, 7LL)*Power(xi, 29LL) +

        6300LL*Power(xj, 22LL) + 11025LL*r*xi*Power(xj, 22LL) -

        50LL*Power(r, 5LL)*Power(xi, 27LL)*(441LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3150LL*Power(xi, 2LL)*Power(xj, 20LL)*(-34LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) +

        525LL*r*Power(xi, 3LL)*Power(xj, 20LL)*(-357LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

        420LL*Power(r, 4LL)*Power(xi, 26LL)*(700LL + 19LL*Power(r, 2LL)*Power(xj, 2LL)) +

        1050LL*Power(xi, 4LL)*Power(xj, 18LL)*

        (816LL - 153LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        210LL*r*Power(xi, 5LL)*Power(xj, 18LL)*

        (7140LL - 425LL*Power(r, 2LL)*Power(xj, 2LL) + 3LL*Power(r, 4LL)*Power(xj, 4LL)) +

        42LL*Power(r, 3LL)*Power(xi, 25LL)*

        (-59500LL - 6035LL*Power(r, 2LL)*Power(xj, 2LL) + 18LL*Power(r, 4LL)*Power(xj, 4LL))

        + 84LL*Power(r, 2LL)*Power(xi, 24LL)*(-160650LL - 52700LL*Power(r, 2LL)*Power(xj, 2LL) +

                                              397LL*Power(r, 4LL)*Power(xj, 4LL)) -

        28LL*Power(xi, 12LL)*Power(xj, 10LL)*

        (100849950LL + 27100125LL*Power(r, 2LL)*Power(xj, 2LL) +

            186150LL*Power(r, 4LL)*Power(xj, 4LL) - 2177LL*Power(r, 6LL)*Power(xj, 6LL)) +

        140LL*Power(xi, 6LL)*Power(xj, 16LL)*

        (-30600LL + 9180LL*Power(r, 2LL)*Power(xj, 2LL) -

            255LL*Power(r, 4LL)*Power(xj, 4LL) + Power(r, 6LL)*Power(xj, 6LL)) -

        2380LL*Power(xi, 8LL)*Power(xj, 14LL)*

        (-6300LL + 2700LL*Power(r, 2LL)*Power(xj, 2LL) -

            120LL*Power(r, 4LL)*Power(xj, 4LL) + Power(r, 6LL)*Power(xj, 6LL)) +

        10LL*r*Power(xi, 7LL)*Power(xj, 16LL)*

        (-749700LL + 71400LL*Power(r, 2LL)*Power(xj, 2LL) -

            1071LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

        204LL*r*Power(xi, 15LL)*Power(xj, 8LL)*

        (28962255LL - 1744750LL*Power(r, 2LL)*Power(xj, 2LL) +

            9555LL*Power(r, 4LL)*Power(xj, 4LL) + 6LL*Power(r, 6LL)*Power(xj, 6LL)) -

        42LL*r*Power(xi, 11LL)*Power(xj, 12LL)*

        (-12911925LL - 1634550LL*Power(r, 2LL)*Power(xj, 2LL) -

            7103LL*Power(r, 4LL)*Power(xj, 4LL) + 18LL*Power(r, 6LL)*Power(xj, 6LL)) +

        2LL*r*Power(xi, 9LL)*Power(xj, 14LL)*

        (16948575LL - 1184400LL*Power(r, 2LL)*Power(xj, 2LL) +

            63861LL*Power(r, 4LL)*Power(xj, 4LL) + 50LL*Power(r, 6LL)*Power(xj, 6LL)) +

        28LL*Power(xi, 22LL)*(-2180250LL - 10993050LL*Power(r, 2LL)*Power(xj, 2LL) +

                              14925LL*Power(r, 4LL)*Power(xj, 4LL) + 73LL*Power(r, 6LL)*Power(xj, 6LL)) -

        952LL*Power(xi, 14LL)*Power(xj, 8LL)*

        (16966215LL + 725175LL*Power(r, 2LL)*Power(xj, 2LL) -

            36075LL*Power(r, 4LL)*Power(xj, 4LL) + 79LL*Power(r, 6LL)*Power(xj, 6LL)) -

        84LL*Power(xi, 10LL)*Power(xj, 12LL)*

        (1723800LL + 279225LL*Power(r, 2LL)*Power(xj, 2LL) +

            45600LL*Power(r, 4LL)*Power(xj, 4LL) + 107LL*Power(r, 6LL)*Power(xj, 6LL)) -

        35LL*r*Power(xi, 17LL)*Power(xj, 6LL)*

        (132637869LL - 2205240LL*Power(r, 2LL)*Power(xj, 2LL) -

            48348LL*Power(r, 4LL)*Power(xj, 4LL) + 136LL*Power(r, 6LL)*Power(xj, 6LL)) -

        6LL*r*Power(xi, 21LL)*Power(xj, 2LL)*

        (192298050LL + 12644275LL*Power(r, 2LL)*Power(xj, 2LL) -

            218029LL*Power(r, 4LL)*Power(xj, 4LL) + 204LL*Power(r, 6LL)*Power(xj, 6LL)) +

        4LL*r*Power(xi, 13LL)*Power(xj, 10LL)*

        (1259522775LL + 15895425LL*Power(r, 2LL)*Power(xj, 2LL) -

            493017LL*Power(r, 4LL)*Power(xj, 4LL) + 263LL*Power(r, 6LL)*Power(xj, 6LL)) -

        140LL*Power(xi, 16LL)*Power(xj, 6LL)*

        (180826281LL - 15101406LL*Power(r, 2LL)*Power(xj, 2LL) +

            160140LL*Power(r, 4LL)*Power(xj, 4LL) + 442LL*Power(r, 6LL)*Power(xj, 6LL)) -

        2LL*r*Power(xi, 23LL)*(21366450LL + 23526300LL*Power(r, 2LL)*Power(xj, 2LL) -

                               246729LL*Power(r, 4LL)*Power(xj, 4LL) + 526LL*Power(r, 6LL)*Power(xj, 6LL)) +

        7LL*r*Power(xi, 19LL)*Power(xj, 4LL)*

        (-811081215LL + 39095550LL*Power(r, 2LL)*Power(xj, 2LL) -

            515916LL*Power(r, 4LL)*Power(xj, 4LL) + 680LL*Power(r, 6LL)*Power(xj, 6LL)) +

        70LL*Power(xi, 18LL)*Power(xj, 4LL)*

        (-180554454LL + 9873711LL*Power(r, 2LL)*Power(xj, 2LL) -

            414120LL*Power(r, 4LL)*Power(xj, 4LL) + 2924LL*Power(r, 6LL)*Power(xj, 6LL)) -

        14LL*Power(xi, 20LL)*Power(xj, 2LL)*

        (136919700LL + 71867115LL*Power(r, 2LL)*Power(xj, 2LL) -

            2154150LL*Power(r, 4LL)*Power(xj, 4LL) + 10268LL*Power(r, 6LL)*Power(xj, 6LL))) -

                 4LL*exp(2LL*r*xi)*Power(xi, 10LL)*

                 (-10710LL*Power(xi, 12LL)*Power(xj, 12LL)*

        (-127008LL*xj + 276768LL*r*Power(xj, 2LL) -

            223668LL*Power(r, 2LL)*Power(xj, 3LL) - 89136LL*Power(r, 3LL)*Power(xj, 4LL) +

            2040LL*Power(r, 4LL)*Power(xj, 5LL) + 3456LL*Power(r, 5LL)*Power(xj, 6LL) +

            420LL*Power(r, 6LL)*Power(xj, 7LL) + 16LL*Power(r, 7LL)*Power(xj, 8LL)) +

        2LL*Power(xi, 20LL)*Power(xj, 4LL)*

        (1735020LL*xj + 3084480LL*r*Power(xj, 2LL) +

            2698920LL*Power(r, 2LL)*Power(xj, 3LL) +

            1542240LL*Power(r, 3LL)*Power(xj, 4LL) +

            642600LL*Power(r, 4LL)*Power(xj, 5LL) + 205632LL*Power(r, 5LL)*Power(xj, 6LL) +

            63882LL*Power(r, 6LL)*Power(xj, 7LL) + 2664LL*Power(r, 7LL)*Power(xj, 8LL) -

            180LL*Power(r, 8LL)*Power(xj, 9LL)) -

        2LL*Power(xj, 24LL)*(107137485LL*xj + 90221040LL*r*Power(xj, 2LL) +

                             35085960LL*Power(r, 2LL)*Power(xj, 3LL) +

                             8255520LL*Power(r, 3LL)*Power(xj, 4LL) +

                             1289925LL*Power(r, 4LL)*Power(xj, 5LL) +

                             137592LL*Power(r, 5LL)*Power(xj, 6LL) + 9828LL*Power(r, 6LL)*Power(xj, 7LL) +

                             432LL*Power(r, 7LL)*Power(xj, 8LL) + 9LL*Power(r, 8LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 2LL)*Power(xj, 22LL)*

        (-2505368880LL*xj - 1762780320LL*r*Power(xj, 2LL) -

            557326980LL*Power(r, 2LL)*Power(xj, 3LL) -

            102558960LL*Power(r, 3LL)*Power(xj, 4LL) -

            11807775LL*Power(r, 4LL)*Power(xj, 5LL) -

            836136LL*Power(r, 5LL)*Power(xj, 6LL) - 31374LL*Power(r, 6LL)*Power(xj, 7LL) -

            216LL*Power(r, 7LL)*Power(xj, 8LL) + 18LL*Power(r, 8LL)*Power(xj, 9LL)) +

        Power(xi, 24LL)*(25515LL*xj + 45360LL*r*Power(xj, 2LL) +

                         39690LL*Power(r, 2LL)*Power(xj, 3LL) + 22680LL*Power(r, 3LL)*Power(xj, 4LL) +

                         9450LL*Power(r, 4LL)*Power(xj, 5LL) + 3024LL*Power(r, 5LL)*Power(xj, 6LL) +

                         756LL*Power(r, 6LL)*Power(xj, 7LL) + 144LL*Power(r, 7LL)*Power(xj, 8LL) +

                         18LL*Power(r, 8LL)*Power(xj, 9LL)) -

        102LL*Power(xi, 10LL)*Power(xj, 14LL)*

        (-97433280LL*xj + 88935840LL*r*Power(xj, 2LL) +

            47571300LL*Power(r, 2LL)*Power(xj, 3LL) -

            1829520LL*Power(r, 3LL)*Power(xj, 4LL) -

            3102750LL*Power(r, 4LL)*Power(xj, 5LL) -

            498960LL*Power(r, 5LL)*Power(xj, 6LL) - 28476LL*Power(r, 6LL)*Power(xj, 7LL) -

            48LL*Power(r, 7LL)*Power(xj, 8LL) + 36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        102LL*Power(xi, 14LL)*Power(xj, 10LL)*

        (-1437345LL*xj - 4520880LL*r*Power(xj, 2LL) +

            2432430LL*Power(r, 2LL)*Power(xj, 3LL) -

            4226040LL*Power(r, 3LL)*Power(xj, 4LL) -

            1089270LL*Power(r, 4LL)*Power(xj, 5LL) + 39312LL*Power(r, 5LL)*Power(xj, 6LL) +

            26964LL*Power(r, 6LL)*Power(xj, 7LL) + 2064LL*Power(r, 7LL)*Power(xj, 8LL) +

            36LL*Power(r, 8LL)*Power(xj, 9LL)) -

        Power(xi, 22LL)*Power(xj, 2LL)*

        (433755LL*xj + 771120LL*r*Power(xj, 2LL) +

            674730LL*Power(r, 2LL)*Power(xj, 3LL) + 385560LL*Power(r, 3LL)*Power(xj, 4LL) +

            160650LL*Power(r, 4LL)*Power(xj, 5LL) + 51408LL*Power(r, 5LL)*Power(xj, 6LL) +

            12852LL*Power(r, 6LL)*Power(xj, 7LL) + 2448LL*Power(r, 7LL)*Power(xj, 8LL) +

            36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 4LL)*Power(xj, 20LL)*

        (-9823683240LL*xj - 4094647200LL*r*Power(xj, 2LL) -

            387295020LL*Power(r, 2LL)*Power(xj, 3LL) +

            105643440LL*Power(r, 3LL)*Power(xj, 4LL) +

            35471520LL*Power(r, 4LL)*Power(xj, 5LL) +

            4729536LL*Power(r, 5LL)*Power(xj, 6LL) +

            340578LL*Power(r, 6LL)*Power(xj, 7LL) + 12744LL*Power(r, 7LL)*Power(xj, 8LL) +

            180LL*Power(r, 8LL)*Power(xj, 9LL)) -

        6LL*Power(xi, 16LL)*Power(xj, 8LL)*

        (-10120950LL*xj - 17992800LL*r*Power(xj, 2LL) -

            17095050LL*Power(r, 2LL)*Power(xj, 3LL) -

            3591000LL*Power(r, 3LL)*Power(xj, 4LL) -

            8207955LL*Power(r, 4LL)*Power(xj, 5LL) -

            1271592LL*Power(r, 5LL)*Power(xj, 6LL) + 71568LL*Power(r, 6LL)*Power(xj, 7LL) +

            18912LL*Power(r, 7LL)*Power(xj, 8LL) + 657LL*Power(r, 8LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 18LL)*Power(xj, 6LL)*

        (-8675100LL*xj - 15422400LL*r*Power(xj, 2LL) -

            13494600LL*Power(r, 2LL)*Power(xj, 3LL) -

            7711200LL*Power(r, 3LL)*Power(xj, 4LL) -

            2807595LL*Power(r, 4LL)*Power(xj, 5LL) -

            1676808LL*Power(r, 5LL)*Power(xj, 6LL) -

            144774LL*Power(r, 6LL)*Power(xj, 7LL) + 10440LL*Power(r, 7LL)*Power(xj, 8LL) +

            954LL*Power(r, 8LL)*Power(xj, 9LL)) +

        3LL*Power(xi, 8LL)*Power(xj, 16LL)*

        (6940428705LL*xj + 4235369040LL*r*Power(xj, 2LL) -

            690804450LL*Power(r, 2LL)*Power(xj, 3LL) -

            598442040LL*Power(r, 3LL)*Power(xj, 4LL) -

            109121670LL*Power(r, 4LL)*Power(xj, 5LL) -

            7339248LL*Power(r, 5LL)*Power(xj, 6LL) + 88956LL*Power(r, 6LL)*Power(xj, 7LL) +

            35760LL*Power(r, 7LL)*Power(xj, 8LL) + 1314LL*Power(r, 8LL)*Power(xj, 9LL)) -

        Power(xi, 6LL)*Power(xj, 18LL)*

        (7147185255LL*xj - 11603405520LL*r*Power(xj, 2LL) -

            6160165830LL*Power(r, 2LL)*Power(xj, 3LL) -

            1086621480LL*Power(r, 3LL)*Power(xj, 4LL) -

            54324270LL*Power(r, 4LL)*Power(xj, 5LL) +

            8022672LL*Power(r, 5LL)*Power(xj, 6LL) +

            1419012LL*Power(r, 6LL)*Power(xj, 7LL) + 85968LL*Power(r, 7LL)*Power(xj, 8LL) +

            1908LL*Power(r, 8LL)*Power(xj, 9LL))) -

                 8LL*exp(2LL*r*xi)*Power(xi, 11LL)*

                 (-10710LL*Power(xi, 12LL)*Power(xj, 12LL)*

        (-3555LL - 127008LL*r*xj + 138384LL*Power(r, 2LL)*Power(xj, 2LL) -

            74556LL*Power(r, 3LL)*Power(xj, 3LL) - 22284LL*Power(r, 4LL)*Power(xj, 4LL) +

            408LL*Power(r, 5LL)*Power(xj, 5LL) + 576LL*Power(r, 6LL)*Power(xj, 6LL) +

            60LL*Power(r, 7LL)*Power(xj, 7LL) + 2LL*Power(r, 8LL)*Power(xj, 8LL)) +

        2LL*Power(xi, 20LL)*Power(xj, 4LL)*

        (963900LL + 1735020LL*r*xj + 1542240LL*Power(r, 2LL)*Power(xj, 2LL) +

            899640LL*Power(r, 3LL)*Power(xj, 3LL) + 385560LL*Power(r, 4LL)*Power(xj, 4LL) +

            128520LL*Power(r, 5LL)*Power(xj, 5LL) + 34272LL*Power(r, 6LL)*Power(xj, 6LL) +

            9126LL*Power(r, 7LL)*Power(xj, 7LL) + 333LL*Power(r, 8LL)*Power(xj, 8LL) -

            20LL*Power(r, 9LL)*Power(xj, 9LL)) -

        2LL*Power(xj, 24LL)*(119041650LL + 107137485LL*r*xj +

                             45110520LL*Power(r, 2LL)*Power(xj, 2LL) +

                             11695320LL*Power(r, 3LL)*Power(xj, 3LL) +

                             2063880LL*Power(r, 4LL)*Power(xj, 4LL) + 257985LL*Power(r, 5LL)*Power(xj, 5LL) +

                             22932LL*Power(r, 6LL)*Power(xj, 6LL) + 1404LL*Power(r, 7LL)*Power(xj, 7LL) +

                             54LL*Power(r, 8LL)*Power(xj, 8LL) + Power(r, 9LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 2LL)*Power(xj, 22LL)*

        (-3264488325LL - 2505368880LL*r*xj -

            881390160LL*Power(r, 2LL)*Power(xj, 2LL) -

            185775660LL*Power(r, 3LL)*Power(xj, 3LL) -

            25639740LL*Power(r, 4LL)*Power(xj, 4LL) -

            2361555LL*Power(r, 5LL)*Power(xj, 5LL) - 139356LL*Power(r, 6LL)*Power(xj, 6LL) -

            4482LL*Power(r, 7LL)*Power(xj, 7LL) - 27LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 9LL)*Power(xj, 9LL)) +

        Power(xi, 24LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                         13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                         1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                         108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                         2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        102LL*Power(xi, 10LL)*Power(xj, 14LL)*

        (44986725LL - 97433280LL*r*xj + 44467920LL*Power(r, 2LL)*Power(xj, 2LL) +

            15857100LL*Power(r, 3LL)*Power(xj, 3LL) -

            457380LL*Power(r, 4LL)*Power(xj, 4LL) - 620550LL*Power(r, 5LL)*Power(xj, 5LL) -

            83160LL*Power(r, 6LL)*Power(xj, 6LL) - 4068LL*Power(r, 7LL)*Power(xj, 7LL) -

            6LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        102LL*Power(xi, 14LL)*Power(xj, 10LL)*

        (-859950LL - 1437345LL*r*xj - 2260440LL*Power(r, 2LL)*Power(xj, 2LL) +

            810810LL*Power(r, 3LL)*Power(xj, 3LL) - 1056510LL*Power(r, 4LL)*Power(xj, 4LL) -

            217854LL*Power(r, 5LL)*Power(xj, 5LL) + 6552LL*Power(r, 6LL)*Power(xj, 6LL) +

            3852LL*Power(r, 7LL)*Power(xj, 7LL) + 258LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) -

        Power(xi, 22LL)*Power(xj, 2LL)*

        (240975LL + 433755LL*r*xj + 385560LL*Power(r, 2LL)*Power(xj, 2LL) +

            224910LL*Power(r, 3LL)*Power(xj, 3LL) + 96390LL*Power(r, 4LL)*Power(xj, 4LL) +

            32130LL*Power(r, 5LL)*Power(xj, 5LL) + 8568LL*Power(r, 6LL)*Power(xj, 6LL) +

            1836LL*Power(r, 7LL)*Power(xj, 7LL) + 306LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 4LL)*Power(xj, 20LL)*

        (-18032978565LL - 9823683240LL*r*xj -

            2047323600LL*Power(r, 2LL)*Power(xj, 2LL) -

            129098340LL*Power(r, 3LL)*Power(xj, 3LL) +

            26410860LL*Power(r, 4LL)*Power(xj, 4LL) +

            7094304LL*Power(r, 5LL)*Power(xj, 5LL) + 788256LL*Power(r, 6LL)*Power(xj, 6LL) +

            48654LL*Power(r, 7LL)*Power(xj, 7LL) + 1593LL*Power(r, 8LL)*Power(xj, 8LL) +

            20LL*Power(r, 9LL)*Power(xj, 9LL)) -

        6LL*Power(xi, 16LL)*Power(xj, 8LL)*

        (-5622750LL - 10120950LL*r*xj - 8996400LL*Power(r, 2LL)*Power(xj, 2LL) -

            5698350LL*Power(r, 3LL)*Power(xj, 3LL) - 897750LL*Power(r, 4LL)*Power(xj, 4LL) -

            1641591LL*Power(r, 5LL)*Power(xj, 5LL) - 211932LL*Power(r, 6LL)*Power(xj, 6LL) +

            10224LL*Power(r, 7LL)*Power(xj, 7LL) + 2364LL*Power(r, 8LL)*Power(xj, 8LL) +

            73LL*Power(r, 9LL)*Power(xj, 9LL)) +

        2LL*Power(xi, 18LL)*Power(xj, 6LL)*

        (-4819500LL - 8675100LL*r*xj - 7711200LL*Power(r, 2LL)*Power(xj, 2LL) -

            4498200LL*Power(r, 3LL)*Power(xj, 3LL) -

            1927800LL*Power(r, 4LL)*Power(xj, 4LL) - 561519LL*Power(r, 5LL)*Power(xj, 5LL) -

            279468LL*Power(r, 6LL)*Power(xj, 6LL) - 20682LL*Power(r, 7LL)*Power(xj, 7LL) +

            1305LL*Power(r, 8LL)*Power(xj, 8LL) + 106LL*Power(r, 9LL)*Power(xj, 9LL)) +

        3LL*Power(xi, 8LL)*Power(xj, 16LL)*

        (-9364244085LL + 6940428705LL*r*xj +

            2117684520LL*Power(r, 2LL)*Power(xj, 2LL) -

            230268150LL*Power(r, 3LL)*Power(xj, 3LL) -

            149610510LL*Power(r, 4LL)*Power(xj, 4LL) -

            21824334LL*Power(r, 5LL)*Power(xj, 5LL) -

            1223208LL*Power(r, 6LL)*Power(xj, 6LL) + 12708LL*Power(r, 7LL)*Power(xj, 7LL) +

            4470LL*Power(r, 8LL)*Power(xj, 8LL) + 146LL*Power(r, 9LL)*Power(xj, 9LL)) -

        Power(xi, 6LL)*Power(xj, 18LL)*

        (57304872765LL + 7147185255LL*r*xj -

            5801702760LL*Power(r, 2LL)*Power(xj, 2LL) -

            2053388610LL*Power(r, 3LL)*Power(xj, 3LL) -

            271655370LL*Power(r, 4LL)*Power(xj, 4LL) -

            10864854LL*Power(r, 5LL)*Power(xj, 5LL) +

            1337112LL*Power(r, 6LL)*Power(xj, 6LL) + 202716LL*Power(r, 7LL)*Power(xj, 7LL) +

            10746LL*Power(r, 8LL)*Power(xj, 8LL) + 212LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (56700LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 17LL)*Power(xi + xj, 17LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_5S_4S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_4S_5S(r, xj, xi);
}
