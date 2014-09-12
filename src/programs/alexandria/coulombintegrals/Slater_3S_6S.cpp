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

cl_R Slater_3S_6S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (64049LL*xi)/393216LL

            ;
        }
        else
        {
            S = (1LL/r)*((-12804747411456000LL + 12804747411456000LL*exp(2LL*rxi) -

                          23523793155237375LL*rxi - 21438091487562750LL*Power(rxi, 2LL) -

                          12909495448599750LL*Power(rxi, 3LL) - 5771367188086500LL*Power(rxi, 4LL) -

                          2040067705876200LL*Power(rxi, 5LL) - 592812380160000LL*Power(rxi, 6LL) -

                          145331660073600LL*Power(rxi, 7LL) - 30604380206400LL*Power(rxi, 8LL) -

                          5606134934400LL*Power(rxi, 9LL) - 900980720640LL*Power(rxi, 10LL) -

                          127672796160LL*Power(rxi, 11LL) - 15968010240LL*Power(rxi, 12LL) -

                          1754726400LL*Power(rxi, 13LL) - 167116800LL*Power(rxi, 14LL) -

                          13369344LL*Power(rxi, 15LL) - 835584LL*Power(rxi, 16LL) - 32768LL*Power(rxi, 17LL))/

                         (1.2804747411456e16*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 16LL) + 17LL*Power(xi, 15LL)*xj + 136LL*Power(xi, 14LL)*Power(xj, 2LL) +

                        680LL*Power(xi, 13LL)*Power(xj, 3LL) + 2380LL*Power(xi, 12LL)*Power(xj, 4LL) +

                        6188LL*Power(xi, 11LL)*Power(xj, 5LL) + 12376LL*Power(xi, 10LL)*Power(xj, 6LL) +

                        19448LL*Power(xi, 9LL)*Power(xj, 7LL) + 24310LL*Power(xi, 8LL)*Power(xj, 8LL) +

                        24310LL*Power(xi, 7LL)*Power(xj, 9LL) + 19448LL*Power(xi, 6LL)*Power(xj, 10LL) +

                        12376LL*Power(xi, 5LL)*Power(xj, 11LL) + 4760LL*Power(xi, 4LL)*Power(xj, 12LL) +

                        1360LL*Power(xi, 3LL)*Power(xj, 13LL) + 272LL*Power(xi, 2LL)*Power(xj, 14LL) +

                        34LL*xi*Power(xj, 15LL) + 2LL*Power(xj, 16LL)))/(6LL*Power(xi + xj, 17LL))

            ;
        }
        else
        {
            S = (1LL/r)*((2806650LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 17LL) +

                          20790LL*exp(2LL*rxj)*Power(rxj, 14LL)*

                          (-240LL*Power(rxi, 24LL) - 6LL*Power(rxi, 25LL) + 135LL*Power(rxj, 20LL) +

                           225LL*rxi*Power(rxj, 20LL) - 80LL*Power(rxi, 23LL)*(51LL + Power(rxj, 2LL)) +

                           45LL*Power(rxi, 3LL)*Power(rxj, 18LL)*(-85LL + 2LL*Power(rxj, 2LL)) +

                           45LL*Power(rxi, 2LL)*Power(rxj, 18LL)*(-51LL + 4LL*Power(rxj, 2LL)) -

                           30LL*Power(rxi, 22LL)*(1224LL + 137LL*Power(rxj, 2LL)) +

                           3060LL*Power(rxi, 15LL)*Power(rxj, 6LL)*(-11875LL + 146LL*Power(rxj, 2LL)) +

                           2LL*Power(rxi, 9LL)*Power(rxj, 12LL)*

                           (3977235LL + 115260LL*Power(rxj, 2LL) - 47LL*Power(rxj, 4LL)) +

                           1020LL*Power(rxi, 13LL)*Power(rxj, 8LL)*

                           (23775LL - 741LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           6LL*Power(rxi, 5LL)*Power(rxj, 16LL)*

                           (5100LL - 255LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           30LL*Power(rxi, 4LL)*Power(rxj, 16LL)*

                           (612LL - 102LL*Power(rxj, 2LL) + Power(rxj, 4LL)) -

                           510LL*Power(rxi, 6LL)*Power(rxj, 14LL)*

                           (180LL - 48LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           20LL*Power(rxi, 7LL)*Power(rxj, 14LL)*

                           (-1683LL + 1158LL*Power(rxj, 2LL) + 4LL*Power(rxj, 4LL)) +

                           510LL*Power(rxi, 10LL)*Power(rxj, 10LL)*

                           (-83889LL - 7948LL*Power(rxj, 2LL) + 12LL*Power(rxj, 4LL)) -

                           34LL*Power(rxi, 11LL)*Power(rxj, 10LL)*

                           (-1158885LL + 3450LL*Power(rxj, 2LL) + 16LL*Power(rxj, 4LL)) -

                           90LL*Power(rxi, 20LL)*(3876LL + 10354LL*Power(rxj, 2LL) + 29LL*Power(rxj, 4LL)) +

                           1020LL*Power(rxi, 12LL)*Power(rxj, 8LL)*

                           (-172098LL - 26LL*Power(rxj, 2LL) + 31LL*Power(rxj, 4LL)) -

                           1020LL*Power(rxi, 14LL)*Power(rxj, 6LL)*

                           (210168LL - 8596LL*Power(rxj, 2LL) + 39LL*Power(rxj, 4LL)) +

                           2LL*Power(rxi, 21LL)*(-87210LL - 43125LL*Power(rxj, 2LL) + 47LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 17LL)*Power(rxj, 4LL)*

                           (1992273LL - 31144LL*Power(rxj, 2LL) + 68LL*Power(rxj, 4LL)) -

                           90LL*Power(rxi, 8LL)*Power(rxj, 12LL)*

                           (17425LL + 6664LL*Power(rxj, 2LL) + 76LL*Power(rxj, 4LL)) +

                           Power(rxi, 19LL)*Power(rxj, 2LL)*

                           (-5204385LL - 202710LL*Power(rxj, 2LL) + 544LL*Power(rxj, 4LL)) +

                           45LL*Power(rxi, 18LL)*Power(rxj, 2LL)*

                           (-267615LL - 83676LL*Power(rxj, 2LL) + 680LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 16LL)*Power(rxj, 4LL)*

                           (6000651LL - 41616LL*Power(rxj, 2LL) + 952LL*Power(rxj, 4LL))) -

                          exp(2LL*rxi)*Power(rxi, 8LL)*

                          (44LL*Power(rxi, 20LL)*Power(rxj, 6LL)*

                           (-43375500LL - 79521750LL*rxj - 72292500LL*Power(rxj, 2LL) -

           43375500LL*Power(rxj, 3LL) - 19278000LL*Power(rxj, 4LL) -

           6747300LL*Power(rxj, 5LL) - 1927800LL*Power(rxj, 6LL) -

           441180LL*Power(rxj, 7LL) - 109620LL*Power(rxj, 8LL) - 14715LL*Power(rxj, 9LL) -

           690LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) -

                           6LL*Power(rxj, 26LL)*(6547290750LL + 7202019825LL*rxj +

                                                 3790536750LL*Power(rxj, 2LL) + 1263512250LL*Power(rxj, 3LL) +

                                                 297297000LL*Power(rxj, 4LL) + 52026975LL*Power(rxj, 5LL) +

                                                 6936930LL*Power(rxj, 6LL) + 707850LL*Power(rxj, 7LL) + 54450LL*Power(rxj, 8LL) +

                                                 3025LL*Power(rxj, 9LL) + 110LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) -

                           44LL*Power(rxi, 6LL)*Power(rxj, 20LL)*

                           (100049928300LL - 5205782925LL*rxj - 25852279950LL*Power(rxj, 2LL) -

           8238935250LL*Power(rxj, 3LL) - 784614600LL*Power(rxj, 4LL) +

           136745280LL*Power(rxj, 5LL) + 52950240LL*Power(rxj, 6LL) +

           7931520LL*Power(rxj, 7LL) + 685440LL*Power(rxj, 8LL) + 34425LL*Power(rxj, 9LL) +

           822LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) +

                           3LL*Power(rxi, 26LL)*(935550LL + 1715175LL*rxj + 1559250LL*Power(rxj, 2LL) +

                                                 935550LL*Power(rxj, 3LL) + 415800LL*Power(rxj, 4LL) + 145530LL*Power(rxj, 5LL) +

                                                 41580LL*Power(rxj, 6LL) + 9900LL*Power(rxj, 7LL) + 1980LL*Power(rxj, 8LL) +

                                                 330LL*Power(rxj, 9LL) + 44LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           2244LL*Power(rxi, 14LL)*Power(rxj, 12LL)*

                           (-15479100LL - 28676025LL*rxj - 22821750LL*Power(rxj, 2LL) -

           22689450LL*Power(rxj, 3LL) - 1852200LL*Power(rxj, 4LL) -

           2372580LL*Power(rxj, 5LL) - 1252440LL*Power(rxj, 6LL) -

           228600LL*Power(rxj, 7LL) - 15000LL*Power(rxj, 8LL) + 450LL*Power(rxj, 9LL) +

           108LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           2244LL*Power(rxi, 12LL)*Power(rxj, 14LL)*

                           (-27556200LL - 14104125LL*rxj - 108438750LL*Power(rxj, 2LL) +

           15375150LL*Power(rxj, 3LL) - 5632200LL*Power(rxj, 4LL) -

           8370180LL*Power(rxj, 5LL) - 2119320LL*Power(rxj, 6LL) -

           198000LL*Power(rxj, 7LL) + 2400LL*Power(rxj, 8LL) + 2010LL*Power(rxj, 9LL) +

           156LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           330LL*Power(rxi, 18LL)*Power(rxj, 8LL)*

                           (-20241900LL - 37110150LL*rxj - 33736500LL*Power(rxj, 2LL) -

           20241900LL*Power(rxj, 3LL) - 8996400LL*Power(rxj, 4LL) -

           3211803LL*Power(rxj, 5LL) - 773514LL*Power(rxj, 6LL) -

           263898LL*Power(rxj, 7LL) - 53202LL*Power(rxj, 8LL) - 4393LL*Power(rxj, 9LL) -

           62LL*Power(rxj, 10LL) + 6LL*Power(rxj, 11LL)) +

                           165LL*Power(rxi, 8LL)*Power(rxj, 18LL)*

                           (-11754743490LL + 11330341155LL*rxj + 1384290810LL*Power(rxj, 2LL) -

           2116476810LL*Power(rxj, 3LL) - 782225640LL*Power(rxj, 4LL) -

           97437186LL*Power(rxj, 5LL) + 2679012LL*Power(rxj, 6LL) +

           2436804LL*Power(rxj, 7LL) + 347316LL*Power(rxj, 8LL) + 25014LL*Power(rxj, 9LL) +

           916LL*Power(rxj, 10LL) + 12LL*Power(rxj, 11LL)) +

                           2LL*Power(rxi, 2LL)*Power(rxj, 24LL)*

                           (-436049563950LL - 402658381125LL*rxj - 173330907750LL*Power(rxj, 2LL) -

           45555359850LL*Power(rxj, 3LL) - 7994586600LL*Power(rxj, 4LL) -

           948782835LL*Power(rxj, 5LL) - 69999930LL*Power(rxj, 6LL) -

           1737450LL*Power(rxj, 7LL) + 254430LL*Power(rxj, 8LL) + 34155LL*Power(rxj, 9LL) +

           1914LL*Power(rxj, 10LL) + 46LL*Power(rxj, 11LL)) -

                           4LL*Power(rxi, 4LL)*Power(rxj, 22LL)*

                           (921052717200LL + 543777678675LL*rxj + 99905825250LL*Power(rxj, 2LL) -

           10883876850LL*Power(rxj, 3LL) - 9266934600LL*Power(rxj, 4LL) -

           2236505040LL*Power(rxj, 5LL) - 316673280LL*Power(rxj, 6LL) -

           28779300LL*Power(rxj, 7LL) - 1601820LL*Power(rxj, 8LL) -

           40095LL*Power(rxj, 9LL) + 726LL*Power(rxj, 10LL) + 58LL*Power(rxj, 11LL)) +

                           4LL*Power(rxi, 22LL)*Power(rxj, 4LL)*

                           (95426100LL + 174947850LL*rxj + 159043500LL*Power(rxj, 2LL) +

           95426100LL*Power(rxj, 3LL) + 42411600LL*Power(rxj, 4LL) +

           14844060LL*Power(rxj, 5LL) + 4241160LL*Power(rxj, 6LL) +

           1009800LL*Power(rxj, 7LL) + 201960LL*Power(rxj, 8LL) + 37125LL*Power(rxj, 9LL) +

           3102LL*Power(rxj, 10LL) + 58LL*Power(rxj, 11LL)) +

                           66LL*Power(rxi, 16LL)*Power(rxj, 10LL)*

                           (-263144700LL - 482431950LL*rxj - 438574500LL*Power(rxj, 2LL) -

           259704900LL*Power(rxj, 3LL) - 130712400LL*Power(rxj, 4LL) -

           27031095LL*Power(rxj, 5LL) - 13816530LL*Power(rxj, 6LL) -

           4240170LL*Power(rxj, 7LL) - 537330LL*Power(rxj, 8LL) - 20565LL*Power(rxj, 9LL) +

           1146LL*Power(rxj, 10LL) + 86LL*Power(rxj, 11LL)) -

                           Power(rxi, 24LL)*Power(rxj, 2LL)*

                           (47713050LL + 87473925LL*rxj + 79521750LL*Power(rxj, 2LL) +

           47713050LL*Power(rxj, 3LL) + 21205800LL*Power(rxj, 4LL) +

           7422030LL*Power(rxj, 5LL) + 2120580LL*Power(rxj, 6LL) +

           504900LL*Power(rxj, 7LL) + 100980LL*Power(rxj, 8LL) + 16830LL*Power(rxj, 9LL) +

           2244LL*Power(rxj, 10LL) + 92LL*Power(rxj, 11LL)) -

                           33LL*Power(rxi, 10LL)*Power(rxj, 16LL)*

                           (5519319750LL - 27722883825LL*rxj + 11646151650LL*Power(rxj, 2LL) +

           955234350LL*Power(rxj, 3LL) - 2729953800LL*Power(rxj, 4LL) -

           902572650LL*Power(rxj, 5LL) - 105286860LL*Power(rxj, 6LL) +

           622260LL*Power(rxj, 7LL) + 1538340LL*Power(rxj, 8LL) + 178830LL*Power(rxj, 9LL) +

           9060LL*Power(rxj, 10LL) + 172LL*Power(rxj, 11LL))))/

                         (2.80665e6*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 17LL)*Power(rxi + rxj, 17LL))

                         );
        }

    }
    return S;
}


cl_R Slater_6S_3S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_3S_6S(r, xj, xi);
}
