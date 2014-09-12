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

cl_R DSlater_2S_6S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-230286692010375LL*xi + 251073478656000LL*exp(2LL*r*xi)*xi -

                  418999810729500LL*r*Power(xi, 2LL) -

                  377542446781500LL*Power(r, 2LL)*Power(xi, 3LL) -

                  224211667680000LL*Power(r, 3LL)*Power(xi, 4LL) -

                  98491036644000LL*Power(r, 4LL)*Power(xi, 5LL) -

                  34029501105600LL*Power(r, 5LL)*Power(xi, 6LL) -

                  9594152568000LL*Power(r, 6LL)*Power(xi, 7LL) -

                  2258332876800LL*Power(r, 7LL)*Power(xi, 8LL) -

                  449902252800LL*Power(r, 8LL)*Power(xi, 9LL) -

                  76338662400LL*Power(r, 9LL)*Power(xi, 10LL) -

                  11024133120LL*Power(r, 10LL)*Power(xi, 11LL) -

                  1341849600LL*Power(r, 11LL)*Power(xi, 12LL) -

                  134184960LL*Power(r, 12LL)*Power(xi, 13LL) -

                  10321920LL*Power(r, 13LL)*Power(xi, 14LL) - 491520LL*Power(r, 14LL)*Power(xi, 15LL))/

                (1.25536739328e14*exp(2LL*r*xi)*r) +

                (-125536739328000LL + 125536739328000LL*exp(2LL*r*xi) -

                 230286692010375LL*r*xi - 209499905364750LL*Power(r, 2LL)*Power(xi, 2LL) -

                 125847482260500LL*Power(r, 3LL)*Power(xi, 3LL) -

                 56052916920000LL*Power(r, 4LL)*Power(xi, 4LL) -

                 19698207328800LL*Power(r, 5LL)*Power(xi, 5LL) -

                 5671583517600LL*Power(r, 6LL)*Power(xi, 6LL) -

                 1370593224000LL*Power(r, 7LL)*Power(xi, 7LL) -

                 282291609600LL*Power(r, 8LL)*Power(xi, 8LL) -

                 49989139200LL*Power(r, 9LL)*Power(xi, 9LL) -

                 7633866240LL*Power(r, 10LL)*Power(xi, 10LL) -

                 1002193920LL*Power(r, 11LL)*Power(xi, 11LL) -

                 111820800LL*Power(r, 12LL)*Power(xi, 12LL) -

                 10321920LL*Power(r, 13LL)*Power(xi, 13LL) - 737280LL*Power(r, 14LL)*Power(xi, 14LL) -

                 32768LL*Power(r, 15LL)*Power(xi, 15LL))/

                (1.25536739328e14*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-125536739328000LL + 125536739328000LL*exp(2LL*r*xi) -

                     230286692010375LL*r*xi - 209499905364750LL*Power(r, 2LL)*Power(xi, 2LL) -

                     125847482260500LL*Power(r, 3LL)*Power(xi, 3LL) -

                     56052916920000LL*Power(r, 4LL)*Power(xi, 4LL) -

                     19698207328800LL*Power(r, 5LL)*Power(xi, 5LL) -

                     5671583517600LL*Power(r, 6LL)*Power(xi, 6LL) -

                     1370593224000LL*Power(r, 7LL)*Power(xi, 7LL) -

                     282291609600LL*Power(r, 8LL)*Power(xi, 8LL) -

                     49989139200LL*Power(r, 9LL)*Power(xi, 9LL) -

                     7633866240LL*Power(r, 10LL)*Power(xi, 10LL) -

                     1002193920LL*Power(r, 11LL)*Power(xi, 11LL) -

                     111820800LL*Power(r, 12LL)*Power(xi, 12LL) -

                     10321920LL*Power(r, 13LL)*Power(xi, 13LL) - 737280LL*Power(r, 14LL)*Power(xi, 14LL) -

                     32768LL*Power(r, 15LL)*Power(xi, 15LL)))/(6.2768369664e13*exp(2LL*r*xi)*r)

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
            S = (935550LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                 51975LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-144LL*Power(r, 2LL)*Power(xi, 18LL) - 6LL*Power(r, 3LL)*Power(xi, 19LL) -

                  63999LL*r*Power(xi, 11LL)*Power(xj, 6LL) + 18LL*Power(xj, 16LL) +

                  27LL*r*xi*Power(xj, 16LL) +

                  18LL*Power(xi, 2LL)*Power(xj, 14LL)*(-15LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  270LL*Power(xi, 4LL)*Power(xj, 12LL)*(-7LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  3LL*r*Power(xi, 3LL)*Power(xj, 14LL)*(-135LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  918LL*Power(xi, 16LL)*(4LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  117LL*r*Power(xi, 9LL)*Power(xj, 8LL)*(-1045LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  4LL*r*Power(xi, 17LL)*(306LL + 23LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  3LL*r*Power(xi, 15LL)*Power(xj, 2LL)*(9441LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  3LL*r*Power(xi, 7LL)*Power(xj, 10LL)*(27261LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  9LL*r*Power(xi, 13LL)*Power(xj, 4LL)*

                  (-12915LL + 52LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  234LL*Power(xi, 10LL)*Power(xj, 6LL)*(-4209LL + 55LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  78LL*Power(xi, 8LL)*Power(xj, 8LL)*(6655LL + 69LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  90LL*Power(xi, 14LL)*Power(xj, 2LL)*(1117LL + 77LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 5LL)*Power(xj, 12LL)*(6111LL + 92LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  18LL*Power(xi, 6LL)*Power(xj, 10LL)*(3107LL + 259LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  18LL*Power(xi, 12LL)*Power(xj, 4LL)*(-31885LL + 403LL*Power(r, 2LL)*Power(xj, 2LL)))

                 + exp(2LL*r*xi)*Power(xi, 6LL)*(-3465LL*Power(xi, 12LL)*Power(xj, 12LL)*

                                                 (1351350LL + 2483775LL*r*xj + 2189250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  1499400LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                  512400LL*Power(r, 4LL)*Power(xj, 4LL) + 191940LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                  73080LL*Power(r, 6LL)*Power(xj, 6LL) + 18200LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                  2680LL*Power(r, 8LL)*Power(xj, 8LL) + 220LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                  8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                 330LL*Power(xi, 8LL)*Power(xj, 16LL)*

                                                 (-2409750LL - 79762725LL*r*xj - 9440550LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                  6036975LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                  10098900LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                  4800285LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                  1163190LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                  164670LL*Power(r, 7LL)*Power(xj, 7LL) - 13110LL*Power(r, 8LL)*Power(xj, 8LL) -

                                                  365LL*Power(r, 9LL)*Power(xj, 9LL) + 26LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                  2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                 2LL*Power(xj, 24LL)*(1240539300LL + 1516214700LL*r*xj +

                                                                      891891000LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                      334459125LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                      89189100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                      17837820LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                      2744280LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                      326700LL*Power(r, 7LL)*Power(xj, 7LL) + 29700LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                      1980LL*Power(r, 9LL)*Power(xj, 9LL) + 88LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                      2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                 Power(xi, 24LL)*(935550LL + 1715175LL*r*xj +

                                                                  1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                  9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                  330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                  4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                 110LL*Power(xi, 6LL)*Power(xj, 18LL)*

                                                 (-313749450LL + 140006475LL*r*xj + 40682250LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                  63603225LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                  41107500LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                  11688705LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                  1918350LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                  179550LL*Power(r, 7LL)*Power(xj, 7LL) - 5670LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                  735LL*Power(r, 9LL)*Power(xj, 9LL) + 98LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                  4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                 10LL*Power(xi, 2LL)*Power(xj, 22LL)*

                                                 (-2825672850LL - 2653375725LL*r*xj -

                                                  1114863750LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                  260134875LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                  29729700LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  1486485LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                  1295910LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                  272250LL*Power(r, 7LL)*Power(xj, 7LL) + 34650LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                  2915LL*Power(r, 9LL)*Power(xj, 9LL) + 154LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                  4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                 165LL*Power(xi, 16LL)*Power(xj, 8LL)*

                                                 (7739550LL + 14189175LL*r*xj + 12899250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  7739550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                  3439800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  1210860LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                  330120LL*Power(r, 6LL)*Power(xj, 6LL) + 86400LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                  18480LL*Power(r, 8LL)*Power(xj, 8LL) + 2460LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                  168LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                 5LL*Power(xi, 22LL)*Power(xj, 2LL)*

                                                 (2806650LL + 5145525LL*r*xj + 4677750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  2806650LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                  1247400LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  436590LL*Power(r, 5LL)*Power(xj, 5LL) + 124740LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                  29700LL*Power(r, 7LL)*Power(xj, 7LL) + 5940LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                  990LL*Power(r, 9LL)*Power(xj, 9LL) + 132LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                  8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                 55LL*Power(xi, 18LL)*Power(xj, 6LL)*

                                                 (7739550LL + 14189175LL*r*xj + 12899250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  7739550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                  3439800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  1203930LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                  343980LL*Power(r, 6LL)*Power(xj, 6LL) + 80820LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                  17460LL*Power(r, 8LL)*Power(xj, 8LL) + 2790LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                  244LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                 22LL*Power(xi, 4LL)*Power(xj, 20LL)*

                                                 (2199137850LL + 366522975LL*r*xj - 665232750LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                  422542575LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                  123095700LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                  20724795LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                  1838970LL*Power(r, 6LL)*Power(xj, 6LL) + 12150LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                  26910LL*Power(r, 8LL)*Power(xj, 8LL) + 3735LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                  258LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                 33LL*Power(xi, 10LL)*Power(xj, 14LL)*

                                                 (-188215650LL - 280764225LL*r*xj - 416886750LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                  131922000LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                  59043600LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                  34671420LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                  11740680LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                  2266200LL*Power(r, 7LL)*Power(xj, 7LL) -

                                                  255000LL*Power(r, 8LL)*Power(xj, 8LL) - 15060LL*Power(r, 9LL)*Power(xj, 9LL) -

                                                  216LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                 11LL*Power(xi, 20LL)*Power(xj, 4LL)*

                                                 (8930250LL + 16372125LL*r*xj + 14883750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  8930250LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                  3969000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  1389150LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                  396900LL*Power(r, 6LL)*Power(xj, 6LL) + 94500LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                  18900LL*Power(r, 8LL)*Power(xj, 8LL) + 3290LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                  364LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                 33LL*Power(xi, 14LL)*Power(xj, 10LL)*

                                                 (85135050LL + 156080925LL*r*xj + 141891750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                  84848400LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                  38984400LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                  12157740LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                  3814440LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                  1072200LL*Power(r, 7LL)*Power(xj, 7LL) + 198120LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                  21020LL*Power(r, 9LL)*Power(xj, 9LL) + 1096LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                  16LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (935550LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 15LL)*

                 Power(xi + xj, 15LL)) + (935550LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                                          51975LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                                          (-144LL*Power(r, 2LL)*Power(xi, 18LL) - 6LL*Power(r, 3LL)*Power(xi, 19LL) -

                                  63999LL*r*Power(xi, 11LL)*Power(xj, 6LL) + 18LL*Power(xj, 16LL) +

                                  27LL*r*xi*Power(xj, 16LL) +

                                  18LL*Power(xi, 2LL)*Power(xj, 14LL)*(-15LL + Power(r, 2LL)*Power(xj, 2LL)) -

                                  270LL*Power(xi, 4LL)*Power(xj, 12LL)*(-7LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  3LL*r*Power(xi, 3LL)*Power(xj, 14LL)*(-135LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  918LL*Power(xi, 16LL)*(4LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  117LL*r*Power(xi, 9LL)*Power(xj, 8LL)*(-1045LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  4LL*r*Power(xi, 17LL)*(306LL + 23LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  3LL*r*Power(xi, 15LL)*Power(xj, 2LL)*(9441LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  3LL*r*Power(xi, 7LL)*Power(xj, 10LL)*(27261LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  9LL*r*Power(xi, 13LL)*Power(xj, 4LL)*

                                  (-12915LL + 52LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  234LL*Power(xi, 10LL)*Power(xj, 6LL)*(-4209LL + 55LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  78LL*Power(xi, 8LL)*Power(xj, 8LL)*(6655LL + 69LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  90LL*Power(xi, 14LL)*Power(xj, 2LL)*(1117LL + 77LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  r*Power(xi, 5LL)*Power(xj, 12LL)*(6111LL + 92LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  18LL*Power(xi, 6LL)*Power(xj, 10LL)*(3107LL + 259LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  18LL*Power(xi, 12LL)*Power(xj, 4LL)*(-31885LL + 403LL*Power(r, 2LL)*Power(xj, 2LL)))

                                          + exp(2LL*r*xi)*Power(xi, 6LL)*(-3465LL*Power(xi, 12LL)*Power(xj, 12LL)*

                                                                          (1351350LL + 2483775LL*r*xj + 2189250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  1499400LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                  512400LL*Power(r, 4LL)*Power(xj, 4LL) + 191940LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                  73080LL*Power(r, 6LL)*Power(xj, 6LL) + 18200LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                  2680LL*Power(r, 8LL)*Power(xj, 8LL) + 220LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                  8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                          330LL*Power(xi, 8LL)*Power(xj, 16LL)*

                                                                          (-2409750LL - 79762725LL*r*xj - 9440550LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                  6036975LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                                  10098900LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                  4800285LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                                  1163190LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                  164670LL*Power(r, 7LL)*Power(xj, 7LL) - 13110LL*Power(r, 8LL)*Power(xj, 8LL) -

                                                                  365LL*Power(r, 9LL)*Power(xj, 9LL) + 26LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                  2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                          2LL*Power(xj, 24LL)*(1240539300LL + 1516214700LL*r*xj +

                                                                                               891891000LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                                               334459125LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                                               89189100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                                               17837820LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                                               2744280LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                                               326700LL*Power(r, 7LL)*Power(xj, 7LL) + 29700LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                                               1980LL*Power(r, 9LL)*Power(xj, 9LL) + 88LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                                               2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                          Power(xi, 24LL)*(935550LL + 1715175LL*r*xj +

                                                                                           1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                                           935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                                           145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                                           9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                                           330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                                           4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                          110LL*Power(xi, 6LL)*Power(xj, 18LL)*

                                                                          (-313749450LL + 140006475LL*r*xj + 40682250LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                  63603225LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                                  41107500LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                  11688705LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                                  1918350LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                  179550LL*Power(r, 7LL)*Power(xj, 7LL) - 5670LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                  735LL*Power(r, 9LL)*Power(xj, 9LL) + 98LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                  4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                          10LL*Power(xi, 2LL)*Power(xj, 22LL)*

                                                                          (-2825672850LL - 2653375725LL*r*xj -

                                                                  1114863750LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                  260134875LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                                  29729700LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  1486485LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                  1295910LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                  272250LL*Power(r, 7LL)*Power(xj, 7LL) + 34650LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                  2915LL*Power(r, 9LL)*Power(xj, 9LL) + 154LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                  4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                          165LL*Power(xi, 16LL)*Power(xj, 8LL)*

                                                                          (7739550LL + 14189175LL*r*xj + 12899250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  7739550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                  3439800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  1210860LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                  330120LL*Power(r, 6LL)*Power(xj, 6LL) + 86400LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                  18480LL*Power(r, 8LL)*Power(xj, 8LL) + 2460LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                  168LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                          5LL*Power(xi, 22LL)*Power(xj, 2LL)*

                                                                          (2806650LL + 5145525LL*r*xj + 4677750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  2806650LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                  1247400LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  436590LL*Power(r, 5LL)*Power(xj, 5LL) + 124740LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                  29700LL*Power(r, 7LL)*Power(xj, 7LL) + 5940LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                  990LL*Power(r, 9LL)*Power(xj, 9LL) + 132LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                  8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                          55LL*Power(xi, 18LL)*Power(xj, 6LL)*

                                                                          (7739550LL + 14189175LL*r*xj + 12899250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  7739550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                  3439800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  1203930LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                  343980LL*Power(r, 6LL)*Power(xj, 6LL) + 80820LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                  17460LL*Power(r, 8LL)*Power(xj, 8LL) + 2790LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                  244LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                          22LL*Power(xi, 4LL)*Power(xj, 20LL)*

                                                                          (2199137850LL + 366522975LL*r*xj - 665232750LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                  422542575LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                                  123095700LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                  20724795LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                                  1838970LL*Power(r, 6LL)*Power(xj, 6LL) + 12150LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                  26910LL*Power(r, 8LL)*Power(xj, 8LL) + 3735LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                  258LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                          33LL*Power(xi, 10LL)*Power(xj, 14LL)*

                                                                          (-188215650LL - 280764225LL*r*xj - 416886750LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                  131922000LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                                  59043600LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                  34671420LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                                  11740680LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                  2266200LL*Power(r, 7LL)*Power(xj, 7LL) -

                                                                  255000LL*Power(r, 8LL)*Power(xj, 8LL) - 15060LL*Power(r, 9LL)*Power(xj, 9LL) -

                                                                  216LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                          11LL*Power(xi, 20LL)*Power(xj, 4LL)*

                                                                          (8930250LL + 16372125LL*r*xj + 14883750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  8930250LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                  3969000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  1389150LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                  396900LL*Power(r, 6LL)*Power(xj, 6LL) + 94500LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                  18900LL*Power(r, 8LL)*Power(xj, 8LL) + 3290LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                  364LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                          33LL*Power(xi, 14LL)*Power(xj, 10LL)*

                                                                          (85135050LL + 156080925LL*r*xj + 141891750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                  84848400LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                  38984400LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                  12157740LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                  3814440LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                  1072200LL*Power(r, 7LL)*Power(xj, 7LL) + 198120LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                  21020LL*Power(r, 9LL)*Power(xj, 9LL) + 1096LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                  16LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (467775LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 15LL)*Power(xi + xj, 14LL)) -

                (1871100LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 15LL) +

                 51975LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-288LL*r*Power(xi, 18LL) - 18LL*Power(r, 2LL)*Power(xi, 19LL) -

        5508LL*r*Power(xi, 16LL)*Power(xj, 2LL) -

        184LL*Power(r, 2LL)*Power(xi, 17LL)*Power(xj, 2LL) -

        13860LL*r*Power(xi, 14LL)*Power(xj, 4LL) -

        168LL*Power(r, 2LL)*Power(xi, 15LL)*Power(xj, 4LL) -

        63999LL*Power(xi, 11LL)*Power(xj, 6LL) + 14508LL*r*Power(xi, 12LL)*Power(xj, 6LL) +

        936LL*Power(r, 2LL)*Power(xi, 13LL)*Power(xj, 6LL) +

        25740LL*r*Power(xi, 10LL)*Power(xj, 8LL) -

        10764LL*r*Power(xi, 8LL)*Power(xj, 10LL) -

        936LL*Power(r, 2LL)*Power(xi, 9LL)*Power(xj, 10LL) -

        9324LL*r*Power(xi, 6LL)*Power(xj, 12LL) +

        168LL*Power(r, 2LL)*Power(xi, 7LL)*Power(xj, 12LL) -

        540LL*r*Power(xi, 4LL)*Power(xj, 14LL) +

        184LL*Power(r, 2LL)*Power(xi, 5LL)*Power(xj, 14LL) + 27LL*xi*Power(xj, 16LL) +

        36LL*r*Power(xi, 2LL)*Power(xj, 16LL) +

        12LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 16LL) +

        3LL*Power(xi, 3LL)*Power(xj, 14LL)*(-135LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        117LL*Power(xi, 9LL)*Power(xj, 8LL)*(-1045LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

        4LL*Power(xi, 17LL)*(306LL + 23LL*Power(r, 2LL)*Power(xj, 2LL)) -

        3LL*Power(xi, 15LL)*Power(xj, 2LL)*(9441LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*Power(xi, 7LL)*Power(xj, 10LL)*(27261LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

        9LL*Power(xi, 13LL)*Power(xj, 4LL)*(-12915LL + 52LL*Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 5LL)*Power(xj, 12LL)*(6111LL + 92LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 103950LL*exp(2LL*r*xj)*Power(xj, 15LL)*

                 (-144LL*Power(r, 2LL)*Power(xi, 18LL) - 6LL*Power(r, 3LL)*Power(xi, 19LL) -

        63999LL*r*Power(xi, 11LL)*Power(xj, 6LL) + 18LL*Power(xj, 16LL) +

        27LL*r*xi*Power(xj, 16LL) +

        18LL*Power(xi, 2LL)*Power(xj, 14LL)*(-15LL + Power(r, 2LL)*Power(xj, 2LL)) -

        270LL*Power(xi, 4LL)*Power(xj, 12LL)*(-7LL + Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*r*Power(xi, 3LL)*Power(xj, 14LL)*(-135LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        918LL*Power(xi, 16LL)*(4LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) -

        117LL*r*Power(xi, 9LL)*Power(xj, 8LL)*(-1045LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

        4LL*r*Power(xi, 17LL)*(306LL + 23LL*Power(r, 2LL)*Power(xj, 2LL)) -

        3LL*r*Power(xi, 15LL)*Power(xj, 2LL)*(9441LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*r*Power(xi, 7LL)*Power(xj, 10LL)*(27261LL + 28LL*Power(r, 2LL)*Power(xj, 2LL)) +

        9LL*r*Power(xi, 13LL)*Power(xj, 4LL)*(-12915LL + 52LL*Power(r, 2LL)*Power(xj, 2LL)) +

        234LL*Power(xi, 10LL)*Power(xj, 6LL)*(-4209LL + 55LL*Power(r, 2LL)*Power(xj, 2LL)) -

        78LL*Power(xi, 8LL)*Power(xj, 8LL)*(6655LL + 69LL*Power(r, 2LL)*Power(xj, 2LL)) -

        90LL*Power(xi, 14LL)*Power(xj, 2LL)*(1117LL + 77LL*Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 5LL)*Power(xj, 12LL)*(6111LL + 92LL*Power(r, 2LL)*Power(xj, 2LL)) -

        18LL*Power(xi, 6LL)*Power(xj, 10LL)*(3107LL + 259LL*Power(r, 2LL)*Power(xj, 2LL)) +

        18LL*Power(xi, 12LL)*Power(xj, 4LL)*(-31885LL + 403LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-3465LL*Power(xi, 12LL)*Power(xj, 12LL)*

        (2483775LL*xj + 4378500LL*r*Power(xj, 2LL) +

            4498200LL*Power(r, 2LL)*Power(xj, 3LL) +

            2049600LL*Power(r, 3LL)*Power(xj, 4LL) +

            959700LL*Power(r, 4LL)*Power(xj, 5LL) + 438480LL*Power(r, 5LL)*Power(xj, 6LL) +

            127400LL*Power(r, 6LL)*Power(xj, 7LL) + 21440LL*Power(r, 7LL)*Power(xj, 8LL) +

            1980LL*Power(r, 8LL)*Power(xj, 9LL) + 80LL*Power(r, 9LL)*Power(xj, 10LL)) +

        330LL*Power(xi, 8LL)*Power(xj, 16LL)*

        (-79762725LL*xj - 18881100LL*r*Power(xj, 2LL) -

            18110925LL*Power(r, 2LL)*Power(xj, 3LL) -

            40395600LL*Power(r, 3LL)*Power(xj, 4LL) -

            24001425LL*Power(r, 4LL)*Power(xj, 5LL) -

            6979140LL*Power(r, 5LL)*Power(xj, 6LL) -

            1152690LL*Power(r, 6LL)*Power(xj, 7LL) -

            104880LL*Power(r, 7LL)*Power(xj, 8LL) - 3285LL*Power(r, 8LL)*Power(xj, 9LL) +

            260LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

        2LL*Power(xj, 24LL)*(1516214700LL*xj + 1783782000LL*r*Power(xj, 2LL) +

                             1003377375LL*Power(r, 2LL)*Power(xj, 3LL) +

                             356756400LL*Power(r, 3LL)*Power(xj, 4LL) +

                             89189100LL*Power(r, 4LL)*Power(xj, 5LL) +

                             16465680LL*Power(r, 5LL)*Power(xj, 6LL) +

                             2286900LL*Power(r, 6LL)*Power(xj, 7LL) +

                             237600LL*Power(r, 7LL)*Power(xj, 8LL) + 17820LL*Power(r, 8LL)*Power(xj, 9LL) +

                             880LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) -

        Power(xi, 24LL)*(1715175LL*xj + 3118500LL*r*Power(xj, 2LL) +

                         2806650LL*Power(r, 2LL)*Power(xj, 3LL) +

                         1663200LL*Power(r, 3LL)*Power(xj, 4LL) +

                         727650LL*Power(r, 4LL)*Power(xj, 5LL) + 249480LL*Power(r, 5LL)*Power(xj, 6LL) +

                         69300LL*Power(r, 6LL)*Power(xj, 7LL) + 15840LL*Power(r, 7LL)*Power(xj, 8LL) +

                         2970LL*Power(r, 8LL)*Power(xj, 9LL) + 440LL*Power(r, 9LL)*Power(xj, 10LL) +

                         44LL*Power(r, 10LL)*Power(xj, 11LL)) -

        110LL*Power(xi, 6LL)*Power(xj, 18LL)*

        (140006475LL*xj + 81364500LL*r*Power(xj, 2LL) -

            190809675LL*Power(r, 2LL)*Power(xj, 3LL) -

            164430000LL*Power(r, 3LL)*Power(xj, 4LL) -

            58443525LL*Power(r, 4LL)*Power(xj, 5LL) -

            11510100LL*Power(r, 5LL)*Power(xj, 6LL) -

            1256850LL*Power(r, 6LL)*Power(xj, 7LL) - 45360LL*Power(r, 7LL)*Power(xj, 8LL) +

            6615LL*Power(r, 8LL)*Power(xj, 9LL) + 980LL*Power(r, 9LL)*Power(xj, 10LL) +

            44LL*Power(r, 10LL)*Power(xj, 11LL)) -

        10LL*Power(xi, 2LL)*Power(xj, 22LL)*

        (-2653375725LL*xj - 2229727500LL*r*Power(xj, 2LL) -

            780404625LL*Power(r, 2LL)*Power(xj, 3LL) -

            118918800LL*Power(r, 3LL)*Power(xj, 4LL) +

            7432425LL*Power(r, 4LL)*Power(xj, 5LL) +

            7775460LL*Power(r, 5LL)*Power(xj, 6LL) +

            1905750LL*Power(r, 6LL)*Power(xj, 7LL) +

            277200LL*Power(r, 7LL)*Power(xj, 8LL) + 26235LL*Power(r, 8LL)*Power(xj, 9LL) +

            1540LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) -

        165LL*Power(xi, 16LL)*Power(xj, 8LL)*

        (14189175LL*xj + 25798500LL*r*Power(xj, 2LL) +

            23218650LL*Power(r, 2LL)*Power(xj, 3LL) +

            13759200LL*Power(r, 3LL)*Power(xj, 4LL) +

            6054300LL*Power(r, 4LL)*Power(xj, 5LL) +

            1980720LL*Power(r, 5LL)*Power(xj, 6LL) +

            604800LL*Power(r, 6LL)*Power(xj, 7LL) + 147840LL*Power(r, 7LL)*Power(xj, 8LL) +

            22140LL*Power(r, 8LL)*Power(xj, 9LL) + 1680LL*Power(r, 9LL)*Power(xj, 10LL) +

            44LL*Power(r, 10LL)*Power(xj, 11LL)) +

        5LL*Power(xi, 22LL)*Power(xj, 2LL)*

        (5145525LL*xj + 9355500LL*r*Power(xj, 2LL) +

            8419950LL*Power(r, 2LL)*Power(xj, 3LL) +

            4989600LL*Power(r, 3LL)*Power(xj, 4LL) +

            2182950LL*Power(r, 4LL)*Power(xj, 5LL) +

            748440LL*Power(r, 5LL)*Power(xj, 6LL) + 207900LL*Power(r, 6LL)*Power(xj, 7LL) +

            47520LL*Power(r, 7LL)*Power(xj, 8LL) + 8910LL*Power(r, 8LL)*Power(xj, 9LL) +

            1320LL*Power(r, 9LL)*Power(xj, 10LL) + 88LL*Power(r, 10LL)*Power(xj, 11LL)) +

        55LL*Power(xi, 18LL)*Power(xj, 6LL)*

        (14189175LL*xj + 25798500LL*r*Power(xj, 2LL) +

            23218650LL*Power(r, 2LL)*Power(xj, 3LL) +

            13759200LL*Power(r, 3LL)*Power(xj, 4LL) +

            6019650LL*Power(r, 4LL)*Power(xj, 5LL) +

            2063880LL*Power(r, 5LL)*Power(xj, 6LL) +

            565740LL*Power(r, 6LL)*Power(xj, 7LL) + 139680LL*Power(r, 7LL)*Power(xj, 8LL) +

            25110LL*Power(r, 8LL)*Power(xj, 9LL) + 2440LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) +

        22LL*Power(xi, 4LL)*Power(xj, 20LL)*

        (366522975LL*xj - 1330465500LL*r*Power(xj, 2LL) -

            1267627725LL*Power(r, 2LL)*Power(xj, 3LL) -

            492382800LL*Power(r, 3LL)*Power(xj, 4LL) -

            103623975LL*Power(r, 4LL)*Power(xj, 5LL) -

            11033820LL*Power(r, 5LL)*Power(xj, 6LL) +

            85050LL*Power(r, 6LL)*Power(xj, 7LL) + 215280LL*Power(r, 7LL)*Power(xj, 8LL) +

            33615LL*Power(r, 8LL)*Power(xj, 9LL) + 2580LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) -

        33LL*Power(xi, 10LL)*Power(xj, 14LL)*

        (-280764225LL*xj - 833773500LL*r*Power(xj, 2LL) -

            395766000LL*Power(r, 2LL)*Power(xj, 3LL) -

            236174400LL*Power(r, 3LL)*Power(xj, 4LL) -

            173357100LL*Power(r, 4LL)*Power(xj, 5LL) -

            70444080LL*Power(r, 5LL)*Power(xj, 6LL) -

            15863400LL*Power(r, 6LL)*Power(xj, 7LL) -

            2040000LL*Power(r, 7LL)*Power(xj, 8LL) -

            135540LL*Power(r, 8LL)*Power(xj, 9LL) - 2160LL*Power(r, 9LL)*Power(xj, 10LL) +

            176LL*Power(r, 10LL)*Power(xj, 11LL)) -

        11LL*Power(xi, 20LL)*Power(xj, 4LL)*

        (16372125LL*xj + 29767500LL*r*Power(xj, 2LL) +

            26790750LL*Power(r, 2LL)*Power(xj, 3LL) +

            15876000LL*Power(r, 3LL)*Power(xj, 4LL) +

            6945750LL*Power(r, 4LL)*Power(xj, 5LL) +

            2381400LL*Power(r, 5LL)*Power(xj, 6LL) +

            661500LL*Power(r, 6LL)*Power(xj, 7LL) + 151200LL*Power(r, 7LL)*Power(xj, 8LL) +

            29610LL*Power(r, 8LL)*Power(xj, 9LL) + 3640LL*Power(r, 9LL)*Power(xj, 10LL) +

            176LL*Power(r, 10LL)*Power(xj, 11LL)) +

        33LL*Power(xi, 14LL)*Power(xj, 10LL)*

        (156080925LL*xj + 283783500LL*r*Power(xj, 2LL) +

            254545200LL*Power(r, 2LL)*Power(xj, 3LL) +

            155937600LL*Power(r, 3LL)*Power(xj, 4LL) +

            60788700LL*Power(r, 4LL)*Power(xj, 5LL) +

            22886640LL*Power(r, 5LL)*Power(xj, 6LL) +

            7505400LL*Power(r, 6LL)*Power(xj, 7LL) +

            1584960LL*Power(r, 7LL)*Power(xj, 8LL) + 189180LL*Power(r, 8LL)*Power(xj, 9LL) +

            10960LL*Power(r, 9LL)*Power(xj, 10LL) + 176LL*Power(r, 10LL)*Power(xj, 11LL))) +

                 2LL*exp(2LL*r*xi)*Power(xi, 7LL)*

                 (-3465LL*Power(xi, 12LL)*Power(xj, 12LL)*

        (1351350LL + 2483775LL*r*xj + 2189250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1499400LL*Power(r, 3LL)*Power(xj, 3LL) + 512400LL*Power(r, 4LL)*Power(xj, 4LL) +

            191940LL*Power(r, 5LL)*Power(xj, 5LL) + 73080LL*Power(r, 6LL)*Power(xj, 6LL) +

            18200LL*Power(r, 7LL)*Power(xj, 7LL) + 2680LL*Power(r, 8LL)*Power(xj, 8LL) +

            220LL*Power(r, 9LL)*Power(xj, 9LL) + 8LL*Power(r, 10LL)*Power(xj, 10LL)) +

        330LL*Power(xi, 8LL)*Power(xj, 16LL)*

        (-2409750LL - 79762725LL*r*xj - 9440550LL*Power(r, 2LL)*Power(xj, 2LL) -

            6036975LL*Power(r, 3LL)*Power(xj, 3LL) -

            10098900LL*Power(r, 4LL)*Power(xj, 4LL) -

            4800285LL*Power(r, 5LL)*Power(xj, 5LL) -

            1163190LL*Power(r, 6LL)*Power(xj, 6LL) - 164670LL*Power(r, 7LL)*Power(xj, 7LL) -

            13110LL*Power(r, 8LL)*Power(xj, 8LL) - 365LL*Power(r, 9LL)*Power(xj, 9LL) +

            26LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

        2LL*Power(xj, 24LL)*(1240539300LL + 1516214700LL*r*xj +

                             891891000LL*Power(r, 2LL)*Power(xj, 2LL) +

                             334459125LL*Power(r, 3LL)*Power(xj, 3LL) +

                             89189100LL*Power(r, 4LL)*Power(xj, 4LL) +

                             17837820LL*Power(r, 5LL)*Power(xj, 5LL) +

                             2744280LL*Power(r, 6LL)*Power(xj, 6LL) + 326700LL*Power(r, 7LL)*Power(xj, 7LL) +

                             29700LL*Power(r, 8LL)*Power(xj, 8LL) + 1980LL*Power(r, 9LL)*Power(xj, 9LL) +

                             88LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

        Power(xi, 24LL)*(935550LL + 1715175LL*r*xj +

                         1559250LL*Power(r, 2LL)*Power(xj, 2LL) + 935550LL*Power(r, 3LL)*Power(xj, 3LL) +

                         415800LL*Power(r, 4LL)*Power(xj, 4LL) + 145530LL*Power(r, 5LL)*Power(xj, 5LL) +

                         41580LL*Power(r, 6LL)*Power(xj, 6LL) + 9900LL*Power(r, 7LL)*Power(xj, 7LL) +

                         1980LL*Power(r, 8LL)*Power(xj, 8LL) + 330LL*Power(r, 9LL)*Power(xj, 9LL) +

                         44LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

        110LL*Power(xi, 6LL)*Power(xj, 18LL)*

        (-313749450LL + 140006475LL*r*xj + 40682250LL*Power(r, 2LL)*Power(xj, 2LL) -

            63603225LL*Power(r, 3LL)*Power(xj, 3LL) -

            41107500LL*Power(r, 4LL)*Power(xj, 4LL) -

            11688705LL*Power(r, 5LL)*Power(xj, 5LL) -

            1918350LL*Power(r, 6LL)*Power(xj, 6LL) - 179550LL*Power(r, 7LL)*Power(xj, 7LL) -

            5670LL*Power(r, 8LL)*Power(xj, 8LL) + 735LL*Power(r, 9LL)*Power(xj, 9LL) +

            98LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

        10LL*Power(xi, 2LL)*Power(xj, 22LL)*

        (-2825672850LL - 2653375725LL*r*xj -

            1114863750LL*Power(r, 2LL)*Power(xj, 2LL) -

            260134875LL*Power(r, 3LL)*Power(xj, 3LL) -

            29729700LL*Power(r, 4LL)*Power(xj, 4LL) +

            1486485LL*Power(r, 5LL)*Power(xj, 5LL) +

            1295910LL*Power(r, 6LL)*Power(xj, 6LL) + 272250LL*Power(r, 7LL)*Power(xj, 7LL) +

            34650LL*Power(r, 8LL)*Power(xj, 8LL) + 2915LL*Power(r, 9LL)*Power(xj, 9LL) +

            154LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

        165LL*Power(xi, 16LL)*Power(xj, 8LL)*

        (7739550LL + 14189175LL*r*xj + 12899250LL*Power(r, 2LL)*Power(xj, 2LL) +

            7739550LL*Power(r, 3LL)*Power(xj, 3LL) +

            3439800LL*Power(r, 4LL)*Power(xj, 4LL) +

            1210860LL*Power(r, 5LL)*Power(xj, 5LL) + 330120LL*Power(r, 6LL)*Power(xj, 6LL) +

            86400LL*Power(r, 7LL)*Power(xj, 7LL) + 18480LL*Power(r, 8LL)*Power(xj, 8LL) +

            2460LL*Power(r, 9LL)*Power(xj, 9LL) + 168LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

        5LL*Power(xi, 22LL)*Power(xj, 2LL)*

        (2806650LL + 5145525LL*r*xj + 4677750LL*Power(r, 2LL)*Power(xj, 2LL) +

            2806650LL*Power(r, 3LL)*Power(xj, 3LL) +

            1247400LL*Power(r, 4LL)*Power(xj, 4LL) + 436590LL*Power(r, 5LL)*Power(xj, 5LL) +

            124740LL*Power(r, 6LL)*Power(xj, 6LL) + 29700LL*Power(r, 7LL)*Power(xj, 7LL) +

            5940LL*Power(r, 8LL)*Power(xj, 8LL) + 990LL*Power(r, 9LL)*Power(xj, 9LL) +

            132LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

        55LL*Power(xi, 18LL)*Power(xj, 6LL)*

        (7739550LL + 14189175LL*r*xj + 12899250LL*Power(r, 2LL)*Power(xj, 2LL) +

            7739550LL*Power(r, 3LL)*Power(xj, 3LL) +

            3439800LL*Power(r, 4LL)*Power(xj, 4LL) +

            1203930LL*Power(r, 5LL)*Power(xj, 5LL) + 343980LL*Power(r, 6LL)*Power(xj, 6LL) +

            80820LL*Power(r, 7LL)*Power(xj, 7LL) + 17460LL*Power(r, 8LL)*Power(xj, 8LL) +

            2790LL*Power(r, 9LL)*Power(xj, 9LL) + 244LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

        22LL*Power(xi, 4LL)*Power(xj, 20LL)*

        (2199137850LL + 366522975LL*r*xj - 665232750LL*Power(r, 2LL)*Power(xj, 2LL) -

            422542575LL*Power(r, 3LL)*Power(xj, 3LL) -

            123095700LL*Power(r, 4LL)*Power(xj, 4LL) -

            20724795LL*Power(r, 5LL)*Power(xj, 5LL) -

            1838970LL*Power(r, 6LL)*Power(xj, 6LL) + 12150LL*Power(r, 7LL)*Power(xj, 7LL) +

            26910LL*Power(r, 8LL)*Power(xj, 8LL) + 3735LL*Power(r, 9LL)*Power(xj, 9LL) +

            258LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) -

        33LL*Power(xi, 10LL)*Power(xj, 14LL)*

        (-188215650LL - 280764225LL*r*xj - 416886750LL*Power(r, 2LL)*Power(xj, 2LL) -

            131922000LL*Power(r, 3LL)*Power(xj, 3LL) -

            59043600LL*Power(r, 4LL)*Power(xj, 4LL) -

            34671420LL*Power(r, 5LL)*Power(xj, 5LL) -

            11740680LL*Power(r, 6LL)*Power(xj, 6LL) -

            2266200LL*Power(r, 7LL)*Power(xj, 7LL) - 255000LL*Power(r, 8LL)*Power(xj, 8LL) -

            15060LL*Power(r, 9LL)*Power(xj, 9LL) - 216LL*Power(r, 10LL)*Power(xj, 10LL) +

            16LL*Power(r, 11LL)*Power(xj, 11LL)) -

        11LL*Power(xi, 20LL)*Power(xj, 4LL)*

        (8930250LL + 16372125LL*r*xj + 14883750LL*Power(r, 2LL)*Power(xj, 2LL) +

            8930250LL*Power(r, 3LL)*Power(xj, 3LL) +

            3969000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1389150LL*Power(r, 5LL)*Power(xj, 5LL) + 396900LL*Power(r, 6LL)*Power(xj, 6LL) +

            94500LL*Power(r, 7LL)*Power(xj, 7LL) + 18900LL*Power(r, 8LL)*Power(xj, 8LL) +

            3290LL*Power(r, 9LL)*Power(xj, 9LL) + 364LL*Power(r, 10LL)*Power(xj, 10LL) +

            16LL*Power(r, 11LL)*Power(xj, 11LL)) +

        33LL*Power(xi, 14LL)*Power(xj, 10LL)*

        (85135050LL + 156080925LL*r*xj + 141891750LL*Power(r, 2LL)*Power(xj, 2LL) +

            84848400LL*Power(r, 3LL)*Power(xj, 3LL) +

            38984400LL*Power(r, 4LL)*Power(xj, 4LL) +

            12157740LL*Power(r, 5LL)*Power(xj, 5LL) +

            3814440LL*Power(r, 6LL)*Power(xj, 6LL) + 1072200LL*Power(r, 7LL)*Power(xj, 7LL) +

            198120LL*Power(r, 8LL)*Power(xj, 8LL) + 21020LL*Power(r, 9LL)*Power(xj, 9LL) +

            1096LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (935550LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 15LL)*Power(xi + xj, 15LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_6S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_2S_6S(r, xj, xi);
}
