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

cl_R DSlater_3S_6S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-23523793155237375LL*xi + 25609494822912000LL*exp(2LL*r*xi)*xi -

                  42876182975125500LL*r*Power(xi, 2LL) -

                  38728486345799250LL*Power(r, 2LL)*Power(xi, 3LL) -

                  23085468752346000LL*Power(r, 3LL)*Power(xi, 4LL) -

                  10200338529381000LL*Power(r, 4LL)*Power(xi, 5LL) -

                  3556874280960000LL*Power(r, 5LL)*Power(xi, 6LL) -

                  1017321620515200LL*Power(r, 6LL)*Power(xi, 7LL) -

                  244835041651200LL*Power(r, 7LL)*Power(xi, 8LL) -

                  50455214409600LL*Power(r, 8LL)*Power(xi, 9LL) -

                  9009807206400LL*Power(r, 9LL)*Power(xi, 10LL) -

                  1404400757760LL*Power(r, 10LL)*Power(xi, 11LL) -

                  191616122880LL*Power(r, 11LL)*Power(xi, 12LL) -

                  22811443200LL*Power(r, 12LL)*Power(xi, 13LL) -

                  2339635200LL*Power(r, 13LL)*Power(xi, 14LL) -

                  200540160LL*Power(r, 14LL)*Power(xi, 15LL) -

                  13369344LL*Power(r, 15LL)*Power(xi, 16LL) - 557056LL*Power(r, 16LL)*Power(xi, 17LL))/

                (1.2804747411456e16*exp(2LL*r*xi)*r) +

                (-12804747411456000LL + 12804747411456000LL*exp(2LL*r*xi) -

                 23523793155237375LL*r*xi - 21438091487562750LL*Power(r, 2LL)*Power(xi, 2LL) -

                 12909495448599750LL*Power(r, 3LL)*Power(xi, 3LL) -

                 5771367188086500LL*Power(r, 4LL)*Power(xi, 4LL) -

                 2040067705876200LL*Power(r, 5LL)*Power(xi, 5LL) -

                 592812380160000LL*Power(r, 6LL)*Power(xi, 6LL) -

                 145331660073600LL*Power(r, 7LL)*Power(xi, 7LL) -

                 30604380206400LL*Power(r, 8LL)*Power(xi, 8LL) -

                 5606134934400LL*Power(r, 9LL)*Power(xi, 9LL) -

                 900980720640LL*Power(r, 10LL)*Power(xi, 10LL) -

                 127672796160LL*Power(r, 11LL)*Power(xi, 11LL) -

                 15968010240LL*Power(r, 12LL)*Power(xi, 12LL) -

                 1754726400LL*Power(r, 13LL)*Power(xi, 13LL) -

                 167116800LL*Power(r, 14LL)*Power(xi, 14LL) -

                 13369344LL*Power(r, 15LL)*Power(xi, 15LL) - 835584LL*Power(r, 16LL)*Power(xi, 16LL) -

                 32768LL*Power(r, 17LL)*Power(xi, 17LL))/

                (1.2804747411456e16*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-12804747411456000LL + 12804747411456000LL*exp(2LL*r*xi) -

                     23523793155237375LL*r*xi - 21438091487562750LL*Power(r, 2LL)*Power(xi, 2LL) -

                     12909495448599750LL*Power(r, 3LL)*Power(xi, 3LL) -

                     5771367188086500LL*Power(r, 4LL)*Power(xi, 4LL) -

                     2040067705876200LL*Power(r, 5LL)*Power(xi, 5LL) -

                     592812380160000LL*Power(r, 6LL)*Power(xi, 6LL) -

                     145331660073600LL*Power(r, 7LL)*Power(xi, 7LL) -

                     30604380206400LL*Power(r, 8LL)*Power(xi, 8LL) -

                     5606134934400LL*Power(r, 9LL)*Power(xi, 9LL) -

                     900980720640LL*Power(r, 10LL)*Power(xi, 10LL) -

                     127672796160LL*Power(r, 11LL)*Power(xi, 11LL) -

                     15968010240LL*Power(r, 12LL)*Power(xi, 12LL) -

                     1754726400LL*Power(r, 13LL)*Power(xi, 13LL) -

                     167116800LL*Power(r, 14LL)*Power(xi, 14LL) -

                     13369344LL*Power(r, 15LL)*Power(xi, 15LL) - 835584LL*Power(r, 16LL)*Power(xi, 16LL) -

                     32768LL*Power(r, 17LL)*Power(xi, 17LL)))/(6.402373705728e15*exp(2LL*r*xi)*r)

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
            S = (2806650LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 17LL) +

                 20790LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-240LL*Power(r, 4LL)*Power(xi, 24LL) - 6LL*Power(r, 5LL)*Power(xi, 25LL) +

                  135LL*Power(xj, 20LL) + 225LL*r*xi*Power(xj, 20LL) -

                  80LL*Power(r, 3LL)*Power(xi, 23LL)*(51LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*r*Power(xi, 3LL)*Power(xj, 18LL)*(-85LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  45LL*Power(xi, 2LL)*Power(xj, 18LL)*(-51LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  30LL*Power(r, 2LL)*Power(xi, 22LL)*(1224LL + 137LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  3060LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

                  (-11875LL + 146LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  2LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

                  (3977235LL + 115260LL*Power(r, 2LL)*Power(xj, 2LL) -

            47LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  1020LL*r*Power(xi, 13LL)*Power(xj, 8LL)*

                  (23775LL - 741LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  6LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

                  (5100LL - 255LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  30LL*Power(xi, 4LL)*Power(xj, 16LL)*

                  (612LL - 102LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) -

                  510LL*Power(xi, 6LL)*Power(xj, 14LL)*

                  (180LL - 48LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  20LL*r*Power(xi, 7LL)*Power(xj, 14LL)*

                  (-1683LL + 1158LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 510LL*Power(xi, 10LL)*Power(xj, 10LL)*

                  (-83889LL - 7948LL*Power(r, 2LL)*Power(xj, 2LL) +

            12LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  34LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

                  (-1158885LL + 3450LL*Power(r, 2LL)*Power(xj, 2LL) +

            16LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  90LL*Power(xi, 20LL)*(3876LL + 10354LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        29LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  1020LL*Power(xi, 12LL)*Power(xj, 8LL)*

                  (-172098LL - 26LL*Power(r, 2LL)*Power(xj, 2LL) + 31LL*Power(r, 4LL)*Power(xj, 4LL))

                  - 1020LL*Power(xi, 14LL)*Power(xj, 6LL)*

                  (210168LL - 8596LL*Power(r, 2LL)*Power(xj, 2LL) +

            39LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  2LL*r*Power(xi, 21LL)*(-87210LL - 43125LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         47LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

                  (1992273LL - 31144LL*Power(r, 2LL)*Power(xj, 2LL) +

            68LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  90LL*Power(xi, 8LL)*Power(xj, 12LL)*

                  (17425LL + 6664LL*Power(r, 2LL)*Power(xj, 2LL) + 76LL*Power(r, 4LL)*Power(xj, 4LL))

                  + r*Power(xi, 19LL)*Power(xj, 2LL)*(-5204385LL - 202710LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                      544LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  45LL*Power(xi, 18LL)*Power(xj, 2LL)*

                  (-267615LL - 83676LL*Power(r, 2LL)*Power(xj, 2LL) +

            680LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  15LL*Power(xi, 16LL)*Power(xj, 4LL)*

                  (6000651LL - 41616LL*Power(r, 2LL)*Power(xj, 2LL) +

            952LL*Power(r, 4LL)*Power(xj, 4LL))) +

                 exp(2LL*r*xi)*Power(xi, 8LL)*

                 (2LL*Power(xi, 2LL)*Power(xj, 24LL)*

                  (436049563950LL + 402658381125LL*r*xj +

            173330907750LL*Power(r, 2LL)*Power(xj, 2LL) +

            45555359850LL*Power(r, 3LL)*Power(xj, 3LL) +

            7994586600LL*Power(r, 4LL)*Power(xj, 4LL) +

            948782835LL*Power(r, 5LL)*Power(xj, 5LL) +

            69999930LL*Power(r, 6LL)*Power(xj, 6LL) +

            1737450LL*Power(r, 7LL)*Power(xj, 7LL) -

            254430LL*Power(r, 8LL)*Power(xj, 8LL) - 34155LL*Power(r, 9LL)*Power(xj, 9LL) -

            1914LL*Power(r, 10LL)*Power(xj, 10LL) - 46LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  44LL*Power(xi, 20LL)*Power(xj, 6LL)*

                  (-43375500LL - 79521750LL*r*xj - 72292500LL*Power(r, 2LL)*Power(xj, 2LL) -

            43375500LL*Power(r, 3LL)*Power(xj, 3LL) -

            19278000LL*Power(r, 4LL)*Power(xj, 4LL) -

            6747300LL*Power(r, 5LL)*Power(xj, 5LL) -

            1927800LL*Power(r, 6LL)*Power(xj, 6LL) -

            441180LL*Power(r, 7LL)*Power(xj, 7LL) - 109620LL*Power(r, 8LL)*Power(xj, 8LL) -

            14715LL*Power(r, 9LL)*Power(xj, 9LL) - 690LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  6LL*Power(xj, 26LL)*(6547290750LL + 7202019825LL*r*xj +

                                       3790536750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       1263512250LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       297297000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       52026975LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       6936930LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       707850LL*Power(r, 7LL)*Power(xj, 7LL) + 54450LL*Power(r, 8LL)*Power(xj, 8LL) +

                                       3025LL*Power(r, 9LL)*Power(xj, 9LL) + 110LL*Power(r, 10LL)*Power(xj, 10LL) +

                                       2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  44LL*Power(xi, 6LL)*Power(xj, 20LL)*

                  (100049928300LL - 5205782925LL*r*xj -

            25852279950LL*Power(r, 2LL)*Power(xj, 2LL) -

            8238935250LL*Power(r, 3LL)*Power(xj, 3LL) -

            784614600LL*Power(r, 4LL)*Power(xj, 4LL) +

            136745280LL*Power(r, 5LL)*Power(xj, 5LL) +

            52950240LL*Power(r, 6LL)*Power(xj, 6LL) +

            7931520LL*Power(r, 7LL)*Power(xj, 7LL) +

            685440LL*Power(r, 8LL)*Power(xj, 8LL) + 34425LL*Power(r, 9LL)*Power(xj, 9LL) +

            822LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  3LL*Power(xi, 26LL)*(935550LL + 1715175LL*r*xj +

                                       1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                       330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                       4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  2244LL*Power(xi, 14LL)*Power(xj, 12LL)*

                  (-15479100LL - 28676025LL*r*xj - 22821750LL*Power(r, 2LL)*Power(xj, 2LL) -

            22689450LL*Power(r, 3LL)*Power(xj, 3LL) -

            1852200LL*Power(r, 4LL)*Power(xj, 4LL) -

            2372580LL*Power(r, 5LL)*Power(xj, 5LL) -

            1252440LL*Power(r, 6LL)*Power(xj, 6LL) -

            228600LL*Power(r, 7LL)*Power(xj, 7LL) - 15000LL*Power(r, 8LL)*Power(xj, 8LL) +

            450LL*Power(r, 9LL)*Power(xj, 9LL) + 108LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  2244LL*Power(xi, 12LL)*Power(xj, 14LL)*

                  (-27556200LL - 14104125LL*r*xj - 108438750LL*Power(r, 2LL)*Power(xj, 2LL) +

            15375150LL*Power(r, 3LL)*Power(xj, 3LL) -

            5632200LL*Power(r, 4LL)*Power(xj, 4LL) -

            8370180LL*Power(r, 5LL)*Power(xj, 5LL) -

            2119320LL*Power(r, 6LL)*Power(xj, 6LL) -

            198000LL*Power(r, 7LL)*Power(xj, 7LL) + 2400LL*Power(r, 8LL)*Power(xj, 8LL) +

            2010LL*Power(r, 9LL)*Power(xj, 9LL) + 156LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  330LL*Power(xi, 18LL)*Power(xj, 8LL)*

                  (-20241900LL - 37110150LL*r*xj - 33736500LL*Power(r, 2LL)*Power(xj, 2LL) -

            20241900LL*Power(r, 3LL)*Power(xj, 3LL) -

            8996400LL*Power(r, 4LL)*Power(xj, 4LL) -

            3211803LL*Power(r, 5LL)*Power(xj, 5LL) -

            773514LL*Power(r, 6LL)*Power(xj, 6LL) - 263898LL*Power(r, 7LL)*Power(xj, 7LL) -

            53202LL*Power(r, 8LL)*Power(xj, 8LL) - 4393LL*Power(r, 9LL)*Power(xj, 9LL) -

            62LL*Power(r, 10LL)*Power(xj, 10LL) + 6LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  165LL*Power(xi, 8LL)*Power(xj, 18LL)*

                  (-11754743490LL + 11330341155LL*r*xj +

            1384290810LL*Power(r, 2LL)*Power(xj, 2LL) -

            2116476810LL*Power(r, 3LL)*Power(xj, 3LL) -

            782225640LL*Power(r, 4LL)*Power(xj, 4LL) -

            97437186LL*Power(r, 5LL)*Power(xj, 5LL) +

            2679012LL*Power(r, 6LL)*Power(xj, 6LL) +

            2436804LL*Power(r, 7LL)*Power(xj, 7LL) +

            347316LL*Power(r, 8LL)*Power(xj, 8LL) + 25014LL*Power(r, 9LL)*Power(xj, 9LL) +

            916LL*Power(r, 10LL)*Power(xj, 10LL) + 12LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  4LL*Power(xi, 4LL)*Power(xj, 22LL)*

                  (921052717200LL + 543777678675LL*r*xj +

            99905825250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10883876850LL*Power(r, 3LL)*Power(xj, 3LL) -

            9266934600LL*Power(r, 4LL)*Power(xj, 4LL) -

            2236505040LL*Power(r, 5LL)*Power(xj, 5LL) -

            316673280LL*Power(r, 6LL)*Power(xj, 6LL) -

            28779300LL*Power(r, 7LL)*Power(xj, 7LL) -

            1601820LL*Power(r, 8LL)*Power(xj, 8LL) - 40095LL*Power(r, 9LL)*Power(xj, 9LL) +

            726LL*Power(r, 10LL)*Power(xj, 10LL) + 58LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  4LL*Power(xi, 22LL)*Power(xj, 4LL)*

                  (95426100LL + 174947850LL*r*xj + 159043500LL*Power(r, 2LL)*Power(xj, 2LL) +

            95426100LL*Power(r, 3LL)*Power(xj, 3LL) +

            42411600LL*Power(r, 4LL)*Power(xj, 4LL) +

            14844060LL*Power(r, 5LL)*Power(xj, 5LL) +

            4241160LL*Power(r, 6LL)*Power(xj, 6LL) +

            1009800LL*Power(r, 7LL)*Power(xj, 7LL) +

            201960LL*Power(r, 8LL)*Power(xj, 8LL) + 37125LL*Power(r, 9LL)*Power(xj, 9LL) +

            3102LL*Power(r, 10LL)*Power(xj, 10LL) + 58LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  66LL*Power(xi, 16LL)*Power(xj, 10LL)*

                  (-263144700LL - 482431950LL*r*xj - 438574500LL*Power(r, 2LL)*Power(xj, 2LL) -

            259704900LL*Power(r, 3LL)*Power(xj, 3LL) -

            130712400LL*Power(r, 4LL)*Power(xj, 4LL) -

            27031095LL*Power(r, 5LL)*Power(xj, 5LL) -

            13816530LL*Power(r, 6LL)*Power(xj, 6LL) -

            4240170LL*Power(r, 7LL)*Power(xj, 7LL) -

            537330LL*Power(r, 8LL)*Power(xj, 8LL) - 20565LL*Power(r, 9LL)*Power(xj, 9LL) +

            1146LL*Power(r, 10LL)*Power(xj, 10LL) + 86LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  Power(xi, 24LL)*Power(xj, 2LL)*

                  (47713050LL + 87473925LL*r*xj + 79521750LL*Power(r, 2LL)*Power(xj, 2LL) +

            47713050LL*Power(r, 3LL)*Power(xj, 3LL) +

            21205800LL*Power(r, 4LL)*Power(xj, 4LL) +

            7422030LL*Power(r, 5LL)*Power(xj, 5LL) +

            2120580LL*Power(r, 6LL)*Power(xj, 6LL) +

            504900LL*Power(r, 7LL)*Power(xj, 7LL) + 100980LL*Power(r, 8LL)*Power(xj, 8LL) +

            16830LL*Power(r, 9LL)*Power(xj, 9LL) + 2244LL*Power(r, 10LL)*Power(xj, 10LL) +

            92LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  33LL*Power(xi, 10LL)*Power(xj, 16LL)*

                  (5519319750LL - 27722883825LL*r*xj +

            11646151650LL*Power(r, 2LL)*Power(xj, 2LL) +

            955234350LL*Power(r, 3LL)*Power(xj, 3LL) -

            2729953800LL*Power(r, 4LL)*Power(xj, 4LL) -

            902572650LL*Power(r, 5LL)*Power(xj, 5LL) -

            105286860LL*Power(r, 6LL)*Power(xj, 6LL) +

            622260LL*Power(r, 7LL)*Power(xj, 7LL) + 1538340LL*Power(r, 8LL)*Power(xj, 8LL) +

            178830LL*Power(r, 9LL)*Power(xj, 9LL) + 9060LL*Power(r, 10LL)*Power(xj, 10LL) +

            172LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (2.80665e6*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 17LL)*

                 Power(xi + xj, 17LL)) + (2806650LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 17LL) +

                                          20790LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                                          (-240LL*Power(r, 4LL)*Power(xi, 24LL) - 6LL*Power(r, 5LL)*Power(xi, 25LL) +

                                  135LL*Power(xj, 20LL) + 225LL*r*xi*Power(xj, 20LL) -

                                  80LL*Power(r, 3LL)*Power(xi, 23LL)*(51LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  45LL*r*Power(xi, 3LL)*Power(xj, 18LL)*(-85LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  45LL*Power(xi, 2LL)*Power(xj, 18LL)*(-51LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  30LL*Power(r, 2LL)*Power(xi, 22LL)*(1224LL + 137LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  3060LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

                                  (-11875LL + 146LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  2LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

                                  (3977235LL + 115260LL*Power(r, 2LL)*Power(xj, 2LL) -

            47LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  1020LL*r*Power(xi, 13LL)*Power(xj, 8LL)*

                                  (23775LL - 741LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  6LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

                                  (5100LL - 255LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  30LL*Power(xi, 4LL)*Power(xj, 16LL)*

                                  (612LL - 102LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) -

                                  510LL*Power(xi, 6LL)*Power(xj, 14LL)*

                                  (180LL - 48LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  20LL*r*Power(xi, 7LL)*Power(xj, 14LL)*

                                  (-1683LL + 1158LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL))

                                  + 510LL*Power(xi, 10LL)*Power(xj, 10LL)*

                                  (-83889LL - 7948LL*Power(r, 2LL)*Power(xj, 2LL) +

            12LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  34LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

                                  (-1158885LL + 3450LL*Power(r, 2LL)*Power(xj, 2LL) +

            16LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  90LL*Power(xi, 20LL)*(3876LL + 10354LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        29LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  1020LL*Power(xi, 12LL)*Power(xj, 8LL)*

                                  (-172098LL - 26LL*Power(r, 2LL)*Power(xj, 2LL) + 31LL*Power(r, 4LL)*Power(xj, 4LL))

                                  - 1020LL*Power(xi, 14LL)*Power(xj, 6LL)*

                                  (210168LL - 8596LL*Power(r, 2LL)*Power(xj, 2LL) +

            39LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  2LL*r*Power(xi, 21LL)*(-87210LL - 43125LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         47LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  15LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

                                  (1992273LL - 31144LL*Power(r, 2LL)*Power(xj, 2LL) +

            68LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  90LL*Power(xi, 8LL)*Power(xj, 12LL)*

                                  (17425LL + 6664LL*Power(r, 2LL)*Power(xj, 2LL) + 76LL*Power(r, 4LL)*Power(xj, 4LL))

                                  + r*Power(xi, 19LL)*Power(xj, 2LL)*(-5204385LL - 202710LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                      544LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  45LL*Power(xi, 18LL)*Power(xj, 2LL)*

                                  (-267615LL - 83676LL*Power(r, 2LL)*Power(xj, 2LL) +

            680LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  15LL*Power(xi, 16LL)*Power(xj, 4LL)*

                                  (6000651LL - 41616LL*Power(r, 2LL)*Power(xj, 2LL) +

            952LL*Power(r, 4LL)*Power(xj, 4LL))) +

                                          exp(2LL*r*xi)*Power(xi, 8LL)*

                                          (2LL*Power(xi, 2LL)*Power(xj, 24LL)*

                                  (436049563950LL + 402658381125LL*r*xj +

            173330907750LL*Power(r, 2LL)*Power(xj, 2LL) +

            45555359850LL*Power(r, 3LL)*Power(xj, 3LL) +

            7994586600LL*Power(r, 4LL)*Power(xj, 4LL) +

            948782835LL*Power(r, 5LL)*Power(xj, 5LL) +

            69999930LL*Power(r, 6LL)*Power(xj, 6LL) +

            1737450LL*Power(r, 7LL)*Power(xj, 7LL) -

            254430LL*Power(r, 8LL)*Power(xj, 8LL) - 34155LL*Power(r, 9LL)*Power(xj, 9LL) -

            1914LL*Power(r, 10LL)*Power(xj, 10LL) - 46LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  44LL*Power(xi, 20LL)*Power(xj, 6LL)*

                                  (-43375500LL - 79521750LL*r*xj - 72292500LL*Power(r, 2LL)*Power(xj, 2LL) -

            43375500LL*Power(r, 3LL)*Power(xj, 3LL) -

            19278000LL*Power(r, 4LL)*Power(xj, 4LL) -

            6747300LL*Power(r, 5LL)*Power(xj, 5LL) -

            1927800LL*Power(r, 6LL)*Power(xj, 6LL) -

            441180LL*Power(r, 7LL)*Power(xj, 7LL) - 109620LL*Power(r, 8LL)*Power(xj, 8LL) -

            14715LL*Power(r, 9LL)*Power(xj, 9LL) - 690LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  6LL*Power(xj, 26LL)*(6547290750LL + 7202019825LL*r*xj +

                                                       3790536750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       1263512250LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                       297297000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                       52026975LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                       6936930LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                       707850LL*Power(r, 7LL)*Power(xj, 7LL) + 54450LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                       3025LL*Power(r, 9LL)*Power(xj, 9LL) + 110LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                       2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  44LL*Power(xi, 6LL)*Power(xj, 20LL)*

                                  (100049928300LL - 5205782925LL*r*xj -

            25852279950LL*Power(r, 2LL)*Power(xj, 2LL) -

            8238935250LL*Power(r, 3LL)*Power(xj, 3LL) -

            784614600LL*Power(r, 4LL)*Power(xj, 4LL) +

            136745280LL*Power(r, 5LL)*Power(xj, 5LL) +

            52950240LL*Power(r, 6LL)*Power(xj, 6LL) +

            7931520LL*Power(r, 7LL)*Power(xj, 7LL) +

            685440LL*Power(r, 8LL)*Power(xj, 8LL) + 34425LL*Power(r, 9LL)*Power(xj, 9LL) +

            822LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  3LL*Power(xi, 26LL)*(935550LL + 1715175LL*r*xj +

                                                       1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                       145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                       9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                       330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                       4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  2244LL*Power(xi, 14LL)*Power(xj, 12LL)*

                                  (-15479100LL - 28676025LL*r*xj - 22821750LL*Power(r, 2LL)*Power(xj, 2LL) -

            22689450LL*Power(r, 3LL)*Power(xj, 3LL) -

            1852200LL*Power(r, 4LL)*Power(xj, 4LL) -

            2372580LL*Power(r, 5LL)*Power(xj, 5LL) -

            1252440LL*Power(r, 6LL)*Power(xj, 6LL) -

            228600LL*Power(r, 7LL)*Power(xj, 7LL) - 15000LL*Power(r, 8LL)*Power(xj, 8LL) +

            450LL*Power(r, 9LL)*Power(xj, 9LL) + 108LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  2244LL*Power(xi, 12LL)*Power(xj, 14LL)*

                                  (-27556200LL - 14104125LL*r*xj - 108438750LL*Power(r, 2LL)*Power(xj, 2LL) +

            15375150LL*Power(r, 3LL)*Power(xj, 3LL) -

            5632200LL*Power(r, 4LL)*Power(xj, 4LL) -

            8370180LL*Power(r, 5LL)*Power(xj, 5LL) -

            2119320LL*Power(r, 6LL)*Power(xj, 6LL) -

            198000LL*Power(r, 7LL)*Power(xj, 7LL) + 2400LL*Power(r, 8LL)*Power(xj, 8LL) +

            2010LL*Power(r, 9LL)*Power(xj, 9LL) + 156LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  330LL*Power(xi, 18LL)*Power(xj, 8LL)*

                                  (-20241900LL - 37110150LL*r*xj - 33736500LL*Power(r, 2LL)*Power(xj, 2LL) -

            20241900LL*Power(r, 3LL)*Power(xj, 3LL) -

            8996400LL*Power(r, 4LL)*Power(xj, 4LL) -

            3211803LL*Power(r, 5LL)*Power(xj, 5LL) -

            773514LL*Power(r, 6LL)*Power(xj, 6LL) - 263898LL*Power(r, 7LL)*Power(xj, 7LL) -

            53202LL*Power(r, 8LL)*Power(xj, 8LL) - 4393LL*Power(r, 9LL)*Power(xj, 9LL) -

            62LL*Power(r, 10LL)*Power(xj, 10LL) + 6LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  165LL*Power(xi, 8LL)*Power(xj, 18LL)*

                                  (-11754743490LL + 11330341155LL*r*xj +

            1384290810LL*Power(r, 2LL)*Power(xj, 2LL) -

            2116476810LL*Power(r, 3LL)*Power(xj, 3LL) -

            782225640LL*Power(r, 4LL)*Power(xj, 4LL) -

            97437186LL*Power(r, 5LL)*Power(xj, 5LL) +

            2679012LL*Power(r, 6LL)*Power(xj, 6LL) +

            2436804LL*Power(r, 7LL)*Power(xj, 7LL) +

            347316LL*Power(r, 8LL)*Power(xj, 8LL) + 25014LL*Power(r, 9LL)*Power(xj, 9LL) +

            916LL*Power(r, 10LL)*Power(xj, 10LL) + 12LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  4LL*Power(xi, 4LL)*Power(xj, 22LL)*

                                  (921052717200LL + 543777678675LL*r*xj +

            99905825250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10883876850LL*Power(r, 3LL)*Power(xj, 3LL) -

            9266934600LL*Power(r, 4LL)*Power(xj, 4LL) -

            2236505040LL*Power(r, 5LL)*Power(xj, 5LL) -

            316673280LL*Power(r, 6LL)*Power(xj, 6LL) -

            28779300LL*Power(r, 7LL)*Power(xj, 7LL) -

            1601820LL*Power(r, 8LL)*Power(xj, 8LL) - 40095LL*Power(r, 9LL)*Power(xj, 9LL) +

            726LL*Power(r, 10LL)*Power(xj, 10LL) + 58LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  4LL*Power(xi, 22LL)*Power(xj, 4LL)*

                                  (95426100LL + 174947850LL*r*xj + 159043500LL*Power(r, 2LL)*Power(xj, 2LL) +

            95426100LL*Power(r, 3LL)*Power(xj, 3LL) +

            42411600LL*Power(r, 4LL)*Power(xj, 4LL) +

            14844060LL*Power(r, 5LL)*Power(xj, 5LL) +

            4241160LL*Power(r, 6LL)*Power(xj, 6LL) +

            1009800LL*Power(r, 7LL)*Power(xj, 7LL) +

            201960LL*Power(r, 8LL)*Power(xj, 8LL) + 37125LL*Power(r, 9LL)*Power(xj, 9LL) +

            3102LL*Power(r, 10LL)*Power(xj, 10LL) + 58LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  66LL*Power(xi, 16LL)*Power(xj, 10LL)*

                                  (-263144700LL - 482431950LL*r*xj - 438574500LL*Power(r, 2LL)*Power(xj, 2LL) -

            259704900LL*Power(r, 3LL)*Power(xj, 3LL) -

            130712400LL*Power(r, 4LL)*Power(xj, 4LL) -

            27031095LL*Power(r, 5LL)*Power(xj, 5LL) -

            13816530LL*Power(r, 6LL)*Power(xj, 6LL) -

            4240170LL*Power(r, 7LL)*Power(xj, 7LL) -

            537330LL*Power(r, 8LL)*Power(xj, 8LL) - 20565LL*Power(r, 9LL)*Power(xj, 9LL) +

            1146LL*Power(r, 10LL)*Power(xj, 10LL) + 86LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  Power(xi, 24LL)*Power(xj, 2LL)*

                                  (47713050LL + 87473925LL*r*xj + 79521750LL*Power(r, 2LL)*Power(xj, 2LL) +

            47713050LL*Power(r, 3LL)*Power(xj, 3LL) +

            21205800LL*Power(r, 4LL)*Power(xj, 4LL) +

            7422030LL*Power(r, 5LL)*Power(xj, 5LL) +

            2120580LL*Power(r, 6LL)*Power(xj, 6LL) +

            504900LL*Power(r, 7LL)*Power(xj, 7LL) + 100980LL*Power(r, 8LL)*Power(xj, 8LL) +

            16830LL*Power(r, 9LL)*Power(xj, 9LL) + 2244LL*Power(r, 10LL)*Power(xj, 10LL) +

            92LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  33LL*Power(xi, 10LL)*Power(xj, 16LL)*

                                  (5519319750LL - 27722883825LL*r*xj +

            11646151650LL*Power(r, 2LL)*Power(xj, 2LL) +

            955234350LL*Power(r, 3LL)*Power(xj, 3LL) -

            2729953800LL*Power(r, 4LL)*Power(xj, 4LL) -

            902572650LL*Power(r, 5LL)*Power(xj, 5LL) -

            105286860LL*Power(r, 6LL)*Power(xj, 6LL) +

            622260LL*Power(r, 7LL)*Power(xj, 7LL) + 1538340LL*Power(r, 8LL)*Power(xj, 8LL) +

            178830LL*Power(r, 9LL)*Power(xj, 9LL) + 9060LL*Power(r, 10LL)*Power(xj, 10LL) +

            172LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (1.403325e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 17LL)*Power(xi + xj, 16LL))

                - (5613300LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                   Power(Power(xi, 2LL) - Power(xj, 2LL), 17LL) +

                   20790LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                   (-960LL*Power(r, 3LL)*Power(xi, 24LL) - 30LL*Power(r, 4LL)*Power(xi, 25LL) -

                    8220LL*Power(r, 3LL)*Power(xi, 22LL)*Power(xj, 2LL) -

                    160LL*Power(r, 4LL)*Power(xi, 23LL)*Power(xj, 2LL) +

                    893520LL*Power(r, 2LL)*Power(xi, 15LL)*Power(xj, 8LL) + 225LL*xi*Power(xj, 20LL) +

                    360LL*r*Power(xi, 2LL)*Power(xj, 20LL) +

                    180LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 20LL) -

                    240LL*Power(r, 2LL)*Power(xi, 23LL)*(51LL + Power(r, 2LL)*Power(xj, 2LL)) +

                    45LL*Power(xi, 3LL)*Power(xj, 18LL)*(-85LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    60LL*r*Power(xi, 22LL)*(1224LL + 137LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    3060LL*Power(xi, 15LL)*Power(xj, 6LL)*

                    (-11875LL + 146LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    2LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

                    (230520LL*r*Power(xj, 2LL) - 188LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    1020LL*r*Power(xi, 13LL)*Power(xj, 8LL)*

                    (-1482LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    6LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

                    (-510LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    30LL*Power(xi, 4LL)*Power(xj, 16LL)*

                    (-204LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    510LL*Power(xi, 6LL)*Power(xj, 14LL)*

                    (-96LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    20LL*r*Power(xi, 7LL)*Power(xj, 14LL)*

                    (2316LL*r*Power(xj, 2LL) + 16LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    510LL*Power(xi, 10LL)*Power(xj, 10LL)*

                    (-15896LL*r*Power(xj, 2LL) + 48LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    34LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

                    (6900LL*r*Power(xj, 2LL) + 64LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    90LL*Power(xi, 20LL)*(20708LL*r*Power(xj, 2LL) +

                                          116LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    1020LL*Power(xi, 12LL)*Power(xj, 8LL)*

                    (-52LL*r*Power(xj, 2LL) + 124LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    1020LL*Power(xi, 14LL)*Power(xj, 6LL)*

                    (-17192LL*r*Power(xj, 2LL) + 156LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    2LL*r*Power(xi, 21LL)*(-86250LL*r*Power(xj, 2LL) +

                                           188LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    15LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

                    (-62288LL*r*Power(xj, 2LL) + 272LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    90LL*Power(xi, 8LL)*Power(xj, 12LL)*

                    (13328LL*r*Power(xj, 2LL) + 304LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    r*Power(xi, 19LL)*Power(xj, 2LL)*

                    (-405420LL*r*Power(xj, 2LL) + 2176LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    45LL*Power(xi, 18LL)*Power(xj, 2LL)*

                    (-167352LL*r*Power(xj, 2LL) + 2720LL*Power(r, 3LL)*Power(xj, 4LL)) -

                    15LL*Power(xi, 16LL)*Power(xj, 4LL)*

                    (-83232LL*r*Power(xj, 2LL) + 3808LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    2LL*Power(xi, 9LL)*Power(xj, 12LL)*

                    (3977235LL + 115260LL*Power(r, 2LL)*Power(xj, 2LL) -

            47LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    1020LL*Power(xi, 13LL)*Power(xj, 8LL)*

                    (23775LL - 741LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    6LL*Power(xi, 5LL)*Power(xj, 16LL)*

                    (5100LL - 255LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    20LL*Power(xi, 7LL)*Power(xj, 14LL)*

                    (-1683LL + 1158LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    34LL*Power(xi, 11LL)*Power(xj, 10LL)*

                    (-1158885LL + 3450LL*Power(r, 2LL)*Power(xj, 2LL) +

            16LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    2LL*Power(xi, 21LL)*(-87210LL - 43125LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         47LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    15LL*Power(xi, 17LL)*Power(xj, 4LL)*

                    (1992273LL - 31144LL*Power(r, 2LL)*Power(xj, 2LL) +

            68LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    Power(xi, 19LL)*Power(xj, 2LL)*

                    (-5204385LL - 202710LL*Power(r, 2LL)*Power(xj, 2LL) +

            544LL*Power(r, 4LL)*Power(xj, 4LL))) +

                   41580LL*exp(2LL*r*xj)*Power(xj, 15LL)*

                   (-240LL*Power(r, 4LL)*Power(xi, 24LL) - 6LL*Power(r, 5LL)*Power(xi, 25LL) +

                    135LL*Power(xj, 20LL) + 225LL*r*xi*Power(xj, 20LL) -

                    80LL*Power(r, 3LL)*Power(xi, 23LL)*(51LL + Power(r, 2LL)*Power(xj, 2LL)) +

                    45LL*r*Power(xi, 3LL)*Power(xj, 18LL)*(-85LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    45LL*Power(xi, 2LL)*Power(xj, 18LL)*(-51LL + 4LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    30LL*Power(r, 2LL)*Power(xi, 22LL)*(1224LL + 137LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    3060LL*r*Power(xi, 15LL)*Power(xj, 6LL)*

                    (-11875LL + 146LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    2LL*r*Power(xi, 9LL)*Power(xj, 12LL)*

                    (3977235LL + 115260LL*Power(r, 2LL)*Power(xj, 2LL) -

            47LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    1020LL*r*Power(xi, 13LL)*Power(xj, 8LL)*

                    (23775LL - 741LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    6LL*r*Power(xi, 5LL)*Power(xj, 16LL)*

                    (5100LL - 255LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    30LL*Power(xi, 4LL)*Power(xj, 16LL)*

                    (612LL - 102LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) -

                    510LL*Power(xi, 6LL)*Power(xj, 14LL)*

                    (180LL - 48LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    20LL*r*Power(xi, 7LL)*Power(xj, 14LL)*

                    (-1683LL + 1158LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    510LL*Power(xi, 10LL)*Power(xj, 10LL)*

                    (-83889LL - 7948LL*Power(r, 2LL)*Power(xj, 2LL) + 12LL*Power(r, 4LL)*Power(xj, 4LL))

                    - 34LL*r*Power(xi, 11LL)*Power(xj, 10LL)*

                    (-1158885LL + 3450LL*Power(r, 2LL)*Power(xj, 2LL) +

            16LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    90LL*Power(xi, 20LL)*(3876LL + 10354LL*Power(r, 2LL)*Power(xj, 2LL) +

                                          29LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    1020LL*Power(xi, 12LL)*Power(xj, 8LL)*

                    (-172098LL - 26LL*Power(r, 2LL)*Power(xj, 2LL) + 31LL*Power(r, 4LL)*Power(xj, 4LL))

                    - 1020LL*Power(xi, 14LL)*Power(xj, 6LL)*

                    (210168LL - 8596LL*Power(r, 2LL)*Power(xj, 2LL) + 39LL*Power(r, 4LL)*Power(xj, 4LL))

                    + 2LL*r*Power(xi, 21LL)*(-87210LL - 43125LL*Power(r, 2LL)*Power(xj, 2LL) +

                                             47LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    15LL*r*Power(xi, 17LL)*Power(xj, 4LL)*

                    (1992273LL - 31144LL*Power(r, 2LL)*Power(xj, 2LL) +

            68LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    90LL*Power(xi, 8LL)*Power(xj, 12LL)*

                    (17425LL + 6664LL*Power(r, 2LL)*Power(xj, 2LL) + 76LL*Power(r, 4LL)*Power(xj, 4LL))

                    + r*Power(xi, 19LL)*Power(xj, 2LL)*(-5204385LL - 202710LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        544LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    45LL*Power(xi, 18LL)*Power(xj, 2LL)*

                    (-267615LL - 83676LL*Power(r, 2LL)*Power(xj, 2LL) +

            680LL*Power(r, 4LL)*Power(xj, 4LL)) -

                    15LL*Power(xi, 16LL)*Power(xj, 4LL)*

                    (6000651LL - 41616LL*Power(r, 2LL)*Power(xj, 2LL) +

            952LL*Power(r, 4LL)*Power(xj, 4LL))) +

                   exp(2LL*r*xi)*Power(xi, 8LL)*

                   (2LL*Power(xi, 2LL)*Power(xj, 24LL)*

                    (402658381125LL*xj + 346661815500LL*r*Power(xj, 2LL) +

            136666079550LL*Power(r, 2LL)*Power(xj, 3LL) +

            31978346400LL*Power(r, 3LL)*Power(xj, 4LL) +

            4743914175LL*Power(r, 4LL)*Power(xj, 5LL) +

            419999580LL*Power(r, 5LL)*Power(xj, 6LL) +

            12162150LL*Power(r, 6LL)*Power(xj, 7LL) -

            2035440LL*Power(r, 7LL)*Power(xj, 8LL) -

            307395LL*Power(r, 8LL)*Power(xj, 9LL) - 19140LL*Power(r, 9LL)*Power(xj, 10LL) -

            506LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    44LL*Power(xi, 20LL)*Power(xj, 6LL)*

                    (-79521750LL*xj - 144585000LL*r*Power(xj, 2LL) -

            130126500LL*Power(r, 2LL)*Power(xj, 3LL) -

            77112000LL*Power(r, 3LL)*Power(xj, 4LL) -

            33736500LL*Power(r, 4LL)*Power(xj, 5LL) -

            11566800LL*Power(r, 5LL)*Power(xj, 6LL) -

            3088260LL*Power(r, 6LL)*Power(xj, 7LL) -

            876960LL*Power(r, 7LL)*Power(xj, 8LL) - 132435LL*Power(r, 8LL)*Power(xj, 9LL) -

            6900LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    6LL*Power(xj, 26LL)*(7202019825LL*xj + 7581073500LL*r*Power(xj, 2LL) +

                                         3790536750LL*Power(r, 2LL)*Power(xj, 3LL) +

                                         1189188000LL*Power(r, 3LL)*Power(xj, 4LL) +

                                         260134875LL*Power(r, 4LL)*Power(xj, 5LL) +

                                         41621580LL*Power(r, 5LL)*Power(xj, 6LL) +

                                         4954950LL*Power(r, 6LL)*Power(xj, 7LL) +

                                         435600LL*Power(r, 7LL)*Power(xj, 8LL) + 27225LL*Power(r, 8LL)*Power(xj, 9LL) +

                                         1100LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    44LL*Power(xi, 6LL)*Power(xj, 20LL)*

                    (-5205782925LL*xj - 51704559900LL*r*Power(xj, 2LL) -

            24716805750LL*Power(r, 2LL)*Power(xj, 3LL) -

            3138458400LL*Power(r, 3LL)*Power(xj, 4LL) +

            683726400LL*Power(r, 4LL)*Power(xj, 5LL) +

            317701440LL*Power(r, 5LL)*Power(xj, 6LL) +

            55520640LL*Power(r, 6LL)*Power(xj, 7LL) +

            5483520LL*Power(r, 7LL)*Power(xj, 8LL) +

            309825LL*Power(r, 8LL)*Power(xj, 9LL) + 8220LL*Power(r, 9LL)*Power(xj, 10LL) +

            22LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    3LL*Power(xi, 26LL)*(1715175LL*xj + 3118500LL*r*Power(xj, 2LL) +

                                         2806650LL*Power(r, 2LL)*Power(xj, 3LL) +

                                         1663200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                         727650LL*Power(r, 4LL)*Power(xj, 5LL) + 249480LL*Power(r, 5LL)*Power(xj, 6LL) +

                                         69300LL*Power(r, 6LL)*Power(xj, 7LL) + 15840LL*Power(r, 7LL)*Power(xj, 8LL) +

                                         2970LL*Power(r, 8LL)*Power(xj, 9LL) + 440LL*Power(r, 9LL)*Power(xj, 10LL) +

                                         44LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    2244LL*Power(xi, 14LL)*Power(xj, 12LL)*

                    (-28676025LL*xj - 45643500LL*r*Power(xj, 2LL) -

            68068350LL*Power(r, 2LL)*Power(xj, 3LL) -

            7408800LL*Power(r, 3LL)*Power(xj, 4LL) -

            11862900LL*Power(r, 4LL)*Power(xj, 5LL) -

            7514640LL*Power(r, 5LL)*Power(xj, 6LL) -

            1600200LL*Power(r, 6LL)*Power(xj, 7LL) -

            120000LL*Power(r, 7LL)*Power(xj, 8LL) + 4050LL*Power(r, 8LL)*Power(xj, 9LL) +

            1080LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    2244LL*Power(xi, 12LL)*Power(xj, 14LL)*

                    (-14104125LL*xj - 216877500LL*r*Power(xj, 2LL) +

            46125450LL*Power(r, 2LL)*Power(xj, 3LL) -

            22528800LL*Power(r, 3LL)*Power(xj, 4LL) -

            41850900LL*Power(r, 4LL)*Power(xj, 5LL) -

            12715920LL*Power(r, 5LL)*Power(xj, 6LL) -

            1386000LL*Power(r, 6LL)*Power(xj, 7LL) + 19200LL*Power(r, 7LL)*Power(xj, 8LL) +

            18090LL*Power(r, 8LL)*Power(xj, 9LL) + 1560LL*Power(r, 9LL)*Power(xj, 10LL) +

            44LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    330LL*Power(xi, 18LL)*Power(xj, 8LL)*

                    (-37110150LL*xj - 67473000LL*r*Power(xj, 2LL) -

            60725700LL*Power(r, 2LL)*Power(xj, 3LL) -

            35985600LL*Power(r, 3LL)*Power(xj, 4LL) -

            16059015LL*Power(r, 4LL)*Power(xj, 5LL) -

            4641084LL*Power(r, 5LL)*Power(xj, 6LL) -

            1847286LL*Power(r, 6LL)*Power(xj, 7LL) -

            425616LL*Power(r, 7LL)*Power(xj, 8LL) - 39537LL*Power(r, 8LL)*Power(xj, 9LL) -

            620LL*Power(r, 9LL)*Power(xj, 10LL) + 66LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    165LL*Power(xi, 8LL)*Power(xj, 18LL)*

                    (11330341155LL*xj + 2768581620LL*r*Power(xj, 2LL) -

            6349430430LL*Power(r, 2LL)*Power(xj, 3LL) -

            3128902560LL*Power(r, 3LL)*Power(xj, 4LL) -

            487185930LL*Power(r, 4LL)*Power(xj, 5LL) +

            16074072LL*Power(r, 5LL)*Power(xj, 6LL) +

            17057628LL*Power(r, 6LL)*Power(xj, 7LL) +

            2778528LL*Power(r, 7LL)*Power(xj, 8LL) +

            225126LL*Power(r, 8LL)*Power(xj, 9LL) + 9160LL*Power(r, 9LL)*Power(xj, 10LL) +

            132LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    4LL*Power(xi, 4LL)*Power(xj, 22LL)*

                    (543777678675LL*xj + 199811650500LL*r*Power(xj, 2LL) -

            32651630550LL*Power(r, 2LL)*Power(xj, 3LL) -

            37067738400LL*Power(r, 3LL)*Power(xj, 4LL) -

            11182525200LL*Power(r, 4LL)*Power(xj, 5LL) -

            1900039680LL*Power(r, 5LL)*Power(xj, 6LL) -

            201455100LL*Power(r, 6LL)*Power(xj, 7LL) -

            12814560LL*Power(r, 7LL)*Power(xj, 8LL) -

            360855LL*Power(r, 8LL)*Power(xj, 9LL) + 7260LL*Power(r, 9LL)*Power(xj, 10LL) +

            638LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    4LL*Power(xi, 22LL)*Power(xj, 4LL)*

                    (174947850LL*xj + 318087000LL*r*Power(xj, 2LL) +

            286278300LL*Power(r, 2LL)*Power(xj, 3LL) +

            169646400LL*Power(r, 3LL)*Power(xj, 4LL) +

            74220300LL*Power(r, 4LL)*Power(xj, 5LL) +

            25446960LL*Power(r, 5LL)*Power(xj, 6LL) +

            7068600LL*Power(r, 6LL)*Power(xj, 7LL) +

            1615680LL*Power(r, 7LL)*Power(xj, 8LL) +

            334125LL*Power(r, 8LL)*Power(xj, 9LL) + 31020LL*Power(r, 9LL)*Power(xj, 10LL) +

            638LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    66LL*Power(xi, 16LL)*Power(xj, 10LL)*

                    (-482431950LL*xj - 877149000LL*r*Power(xj, 2LL) -

            779114700LL*Power(r, 2LL)*Power(xj, 3LL) -

            522849600LL*Power(r, 3LL)*Power(xj, 4LL) -

            135155475LL*Power(r, 4LL)*Power(xj, 5LL) -

            82899180LL*Power(r, 5LL)*Power(xj, 6LL) -

            29681190LL*Power(r, 6LL)*Power(xj, 7LL) -

            4298640LL*Power(r, 7LL)*Power(xj, 8LL) -

            185085LL*Power(r, 8LL)*Power(xj, 9LL) + 11460LL*Power(r, 9LL)*Power(xj, 10LL) +

            946LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    Power(xi, 24LL)*Power(xj, 2LL)*

                    (87473925LL*xj + 159043500LL*r*Power(xj, 2LL) +

            143139150LL*Power(r, 2LL)*Power(xj, 3LL) +

            84823200LL*Power(r, 3LL)*Power(xj, 4LL) +

            37110150LL*Power(r, 4LL)*Power(xj, 5LL) +

            12723480LL*Power(r, 5LL)*Power(xj, 6LL) +

            3534300LL*Power(r, 6LL)*Power(xj, 7LL) +

            807840LL*Power(r, 7LL)*Power(xj, 8LL) + 151470LL*Power(r, 8LL)*Power(xj, 9LL) +

            22440LL*Power(r, 9LL)*Power(xj, 10LL) + 1012LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    33LL*Power(xi, 10LL)*Power(xj, 16LL)*

                    (-27722883825LL*xj + 23292303300LL*r*Power(xj, 2LL) +

            2865703050LL*Power(r, 2LL)*Power(xj, 3LL) -

            10919815200LL*Power(r, 3LL)*Power(xj, 4LL) -

            4512863250LL*Power(r, 4LL)*Power(xj, 5LL) -

            631721160LL*Power(r, 5LL)*Power(xj, 6LL) +

            4355820LL*Power(r, 6LL)*Power(xj, 7LL) +

            12306720LL*Power(r, 7LL)*Power(xj, 8LL) +

            1609470LL*Power(r, 8LL)*Power(xj, 9LL) + 90600LL*Power(r, 9LL)*Power(xj, 10LL) +

            1892LL*Power(r, 10LL)*Power(xj, 11LL))) +

                   2LL*exp(2LL*r*xi)*Power(xi, 9LL)*

                   (2LL*Power(xi, 2LL)*Power(xj, 24LL)*

                    (436049563950LL + 402658381125LL*r*xj +

            173330907750LL*Power(r, 2LL)*Power(xj, 2LL) +

            45555359850LL*Power(r, 3LL)*Power(xj, 3LL) +

            7994586600LL*Power(r, 4LL)*Power(xj, 4LL) +

            948782835LL*Power(r, 5LL)*Power(xj, 5LL) +

            69999930LL*Power(r, 6LL)*Power(xj, 6LL) +

            1737450LL*Power(r, 7LL)*Power(xj, 7LL) - 254430LL*Power(r, 8LL)*Power(xj, 8LL) -

            34155LL*Power(r, 9LL)*Power(xj, 9LL) - 1914LL*Power(r, 10LL)*Power(xj, 10LL) -

            46LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    44LL*Power(xi, 20LL)*Power(xj, 6LL)*

                    (-43375500LL - 79521750LL*r*xj - 72292500LL*Power(r, 2LL)*Power(xj, 2LL) -

            43375500LL*Power(r, 3LL)*Power(xj, 3LL) -

            19278000LL*Power(r, 4LL)*Power(xj, 4LL) -

            6747300LL*Power(r, 5LL)*Power(xj, 5LL) -

            1927800LL*Power(r, 6LL)*Power(xj, 6LL) - 441180LL*Power(r, 7LL)*Power(xj, 7LL) -

            109620LL*Power(r, 8LL)*Power(xj, 8LL) - 14715LL*Power(r, 9LL)*Power(xj, 9LL) -

            690LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    6LL*Power(xj, 26LL)*(6547290750LL + 7202019825LL*r*xj +

                                         3790536750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         1263512250LL*Power(r, 3LL)*Power(xj, 3LL) +

                                         297297000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                         52026975LL*Power(r, 5LL)*Power(xj, 5LL) +

                                         6936930LL*Power(r, 6LL)*Power(xj, 6LL) + 707850LL*Power(r, 7LL)*Power(xj, 7LL) +

                                         54450LL*Power(r, 8LL)*Power(xj, 8LL) + 3025LL*Power(r, 9LL)*Power(xj, 9LL) +

                                         110LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    44LL*Power(xi, 6LL)*Power(xj, 20LL)*

                    (100049928300LL - 5205782925LL*r*xj -

            25852279950LL*Power(r, 2LL)*Power(xj, 2LL) -

            8238935250LL*Power(r, 3LL)*Power(xj, 3LL) -

            784614600LL*Power(r, 4LL)*Power(xj, 4LL) +

            136745280LL*Power(r, 5LL)*Power(xj, 5LL) +

            52950240LL*Power(r, 6LL)*Power(xj, 6LL) +

            7931520LL*Power(r, 7LL)*Power(xj, 7LL) + 685440LL*Power(r, 8LL)*Power(xj, 8LL) +

            34425LL*Power(r, 9LL)*Power(xj, 9LL) + 822LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    3LL*Power(xi, 26LL)*(935550LL + 1715175LL*r*xj +

                                         1559250LL*Power(r, 2LL)*Power(xj, 2LL) + 935550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                         415800LL*Power(r, 4LL)*Power(xj, 4LL) + 145530LL*Power(r, 5LL)*Power(xj, 5LL) +

                                         41580LL*Power(r, 6LL)*Power(xj, 6LL) + 9900LL*Power(r, 7LL)*Power(xj, 7LL) +

                                         1980LL*Power(r, 8LL)*Power(xj, 8LL) + 330LL*Power(r, 9LL)*Power(xj, 9LL) +

                                         44LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    2244LL*Power(xi, 14LL)*Power(xj, 12LL)*

                    (-15479100LL - 28676025LL*r*xj - 22821750LL*Power(r, 2LL)*Power(xj, 2LL) -

            22689450LL*Power(r, 3LL)*Power(xj, 3LL) -

            1852200LL*Power(r, 4LL)*Power(xj, 4LL) -

            2372580LL*Power(r, 5LL)*Power(xj, 5LL) -

            1252440LL*Power(r, 6LL)*Power(xj, 6LL) - 228600LL*Power(r, 7LL)*Power(xj, 7LL) -

            15000LL*Power(r, 8LL)*Power(xj, 8LL) + 450LL*Power(r, 9LL)*Power(xj, 9LL) +

            108LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    2244LL*Power(xi, 12LL)*Power(xj, 14LL)*

                    (-27556200LL - 14104125LL*r*xj - 108438750LL*Power(r, 2LL)*Power(xj, 2LL) +

            15375150LL*Power(r, 3LL)*Power(xj, 3LL) -

            5632200LL*Power(r, 4LL)*Power(xj, 4LL) -

            8370180LL*Power(r, 5LL)*Power(xj, 5LL) -

            2119320LL*Power(r, 6LL)*Power(xj, 6LL) - 198000LL*Power(r, 7LL)*Power(xj, 7LL) +

            2400LL*Power(r, 8LL)*Power(xj, 8LL) + 2010LL*Power(r, 9LL)*Power(xj, 9LL) +

            156LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    330LL*Power(xi, 18LL)*Power(xj, 8LL)*

                    (-20241900LL - 37110150LL*r*xj - 33736500LL*Power(r, 2LL)*Power(xj, 2LL) -

            20241900LL*Power(r, 3LL)*Power(xj, 3LL) -

            8996400LL*Power(r, 4LL)*Power(xj, 4LL) -

            3211803LL*Power(r, 5LL)*Power(xj, 5LL) - 773514LL*Power(r, 6LL)*Power(xj, 6LL) -

            263898LL*Power(r, 7LL)*Power(xj, 7LL) - 53202LL*Power(r, 8LL)*Power(xj, 8LL) -

            4393LL*Power(r, 9LL)*Power(xj, 9LL) - 62LL*Power(r, 10LL)*Power(xj, 10LL) +

            6LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    165LL*Power(xi, 8LL)*Power(xj, 18LL)*

                    (-11754743490LL + 11330341155LL*r*xj +

            1384290810LL*Power(r, 2LL)*Power(xj, 2LL) -

            2116476810LL*Power(r, 3LL)*Power(xj, 3LL) -

            782225640LL*Power(r, 4LL)*Power(xj, 4LL) -

            97437186LL*Power(r, 5LL)*Power(xj, 5LL) +

            2679012LL*Power(r, 6LL)*Power(xj, 6LL) +

            2436804LL*Power(r, 7LL)*Power(xj, 7LL) + 347316LL*Power(r, 8LL)*Power(xj, 8LL) +

            25014LL*Power(r, 9LL)*Power(xj, 9LL) + 916LL*Power(r, 10LL)*Power(xj, 10LL) +

            12LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    4LL*Power(xi, 4LL)*Power(xj, 22LL)*

                    (921052717200LL + 543777678675LL*r*xj +

            99905825250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10883876850LL*Power(r, 3LL)*Power(xj, 3LL) -

            9266934600LL*Power(r, 4LL)*Power(xj, 4LL) -

            2236505040LL*Power(r, 5LL)*Power(xj, 5LL) -

            316673280LL*Power(r, 6LL)*Power(xj, 6LL) -

            28779300LL*Power(r, 7LL)*Power(xj, 7LL) -

            1601820LL*Power(r, 8LL)*Power(xj, 8LL) - 40095LL*Power(r, 9LL)*Power(xj, 9LL) +

            726LL*Power(r, 10LL)*Power(xj, 10LL) + 58LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    4LL*Power(xi, 22LL)*Power(xj, 4LL)*

                    (95426100LL + 174947850LL*r*xj + 159043500LL*Power(r, 2LL)*Power(xj, 2LL) +

            95426100LL*Power(r, 3LL)*Power(xj, 3LL) +

            42411600LL*Power(r, 4LL)*Power(xj, 4LL) +

            14844060LL*Power(r, 5LL)*Power(xj, 5LL) +

            4241160LL*Power(r, 6LL)*Power(xj, 6LL) +

            1009800LL*Power(r, 7LL)*Power(xj, 7LL) + 201960LL*Power(r, 8LL)*Power(xj, 8LL) +

            37125LL*Power(r, 9LL)*Power(xj, 9LL) + 3102LL*Power(r, 10LL)*Power(xj, 10LL) +

            58LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    66LL*Power(xi, 16LL)*Power(xj, 10LL)*

                    (-263144700LL - 482431950LL*r*xj - 438574500LL*Power(r, 2LL)*Power(xj, 2LL) -

            259704900LL*Power(r, 3LL)*Power(xj, 3LL) -

            130712400LL*Power(r, 4LL)*Power(xj, 4LL) -

            27031095LL*Power(r, 5LL)*Power(xj, 5LL) -

            13816530LL*Power(r, 6LL)*Power(xj, 6LL) -

            4240170LL*Power(r, 7LL)*Power(xj, 7LL) - 537330LL*Power(r, 8LL)*Power(xj, 8LL) -

            20565LL*Power(r, 9LL)*Power(xj, 9LL) + 1146LL*Power(r, 10LL)*Power(xj, 10LL) +

            86LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    Power(xi, 24LL)*Power(xj, 2LL)*

                    (47713050LL + 87473925LL*r*xj + 79521750LL*Power(r, 2LL)*Power(xj, 2LL) +

            47713050LL*Power(r, 3LL)*Power(xj, 3LL) +

            21205800LL*Power(r, 4LL)*Power(xj, 4LL) +

            7422030LL*Power(r, 5LL)*Power(xj, 5LL) +

            2120580LL*Power(r, 6LL)*Power(xj, 6LL) + 504900LL*Power(r, 7LL)*Power(xj, 7LL) +

            100980LL*Power(r, 8LL)*Power(xj, 8LL) + 16830LL*Power(r, 9LL)*Power(xj, 9LL) +

            2244LL*Power(r, 10LL)*Power(xj, 10LL) + 92LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    33LL*Power(xi, 10LL)*Power(xj, 16LL)*

                    (5519319750LL - 27722883825LL*r*xj +

            11646151650LL*Power(r, 2LL)*Power(xj, 2LL) +

            955234350LL*Power(r, 3LL)*Power(xj, 3LL) -

            2729953800LL*Power(r, 4LL)*Power(xj, 4LL) -

            902572650LL*Power(r, 5LL)*Power(xj, 5LL) -

            105286860LL*Power(r, 6LL)*Power(xj, 6LL) +

            622260LL*Power(r, 7LL)*Power(xj, 7LL) + 1538340LL*Power(r, 8LL)*Power(xj, 8LL) +

            178830LL*Power(r, 9LL)*Power(xj, 9LL) + 9060LL*Power(r, 10LL)*Power(xj, 10LL) +

            172LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (2.80665e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 17LL)*Power(xi + xj, 17LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_6S_3S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_3S_6S(r, xj, xi);
}
