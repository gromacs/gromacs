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

cl_R Slater_2S_6S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (32555LL*xi)/196608LL

            ;
        }
        else
        {
            S = (1LL/r)*((-125536739328000LL + 125536739328000LL*exp(2LL*rxi) - 230286692010375LL*rxi -

                          209499905364750LL*Power(rxi, 2LL) - 125847482260500LL*Power(rxi, 3LL) -

                          56052916920000LL*Power(rxi, 4LL) - 19698207328800LL*Power(rxi, 5LL) -

                          5671583517600LL*Power(rxi, 6LL) - 1370593224000LL*Power(rxi, 7LL) -

                          282291609600LL*Power(rxi, 8LL) - 49989139200LL*Power(rxi, 9LL) -

                          7633866240LL*Power(rxi, 10LL) - 1002193920LL*Power(rxi, 11LL) -

                          111820800LL*Power(rxi, 12LL) - 10321920LL*Power(rxi, 13LL) -

                          737280LL*Power(rxi, 14LL) - 32768LL*Power(rxi, 15LL))/

                         (1.25536739328e14*exp(2LL*rxi))

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

                        6435LL*Power(xi, 7LL)*Power(xj, 7LL) + 6435LL*Power(xi, 6LL)*Power(xj, 8LL) +

                        5005LL*Power(xi, 5LL)*Power(xj, 9LL) + 3003LL*Power(xi, 4LL)*Power(xj, 10LL) +

                        1365LL*Power(xi, 3LL)*Power(xj, 11LL) + 315LL*Power(xi, 2LL)*Power(xj, 12LL) +

                        45LL*xi*Power(xj, 13LL) + 3LL*Power(xj, 14LL)))/(6LL*Power(xi + xj, 15LL))

            ;
        }
        else
        {
            S = (1LL/r)*((935550LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 15LL) +

                          51975LL*exp(2LL*rxj)*Power(rxj, 14LL)*

                          (-144LL*Power(rxi, 18LL) - 6LL*Power(rxi, 19LL) -

                           63999LL*Power(rxi, 11LL)*Power(rxj, 6LL) + 18LL*Power(rxj, 16LL) +

                           27LL*rxi*Power(rxj, 16LL) + 18LL*Power(rxi, 2LL)*Power(rxj, 14LL)*

                           (-15LL + Power(rxj, 2LL)) -

                           270LL*Power(rxi, 4LL)*Power(rxj, 12LL)*(-7LL + Power(rxj, 2LL)) +

                           3LL*Power(rxi, 3LL)*Power(rxj, 14LL)*(-135LL + 2LL*Power(rxj, 2LL)) -

                           918LL*Power(rxi, 16LL)*(4LL + 3LL*Power(rxj, 2LL)) -

                           117LL*Power(rxi, 9LL)*Power(rxj, 8LL)*(-1045LL + 4LL*Power(rxj, 2LL)) -

                           4LL*Power(rxi, 17LL)*(306LL + 23LL*Power(rxj, 2LL)) -

                           3LL*Power(rxi, 15LL)*Power(rxj, 2LL)*(9441LL + 28LL*Power(rxj, 2LL)) +

                           3LL*Power(rxi, 7LL)*Power(rxj, 10LL)*(27261LL + 28LL*Power(rxj, 2LL)) +

                           9LL*Power(rxi, 13LL)*Power(rxj, 4LL)*(-12915LL + 52LL*Power(rxj, 2LL)) +

                           234LL*Power(rxi, 10LL)*Power(rxj, 6LL)*(-4209LL + 55LL*Power(rxj, 2LL)) -

                           78LL*Power(rxi, 8LL)*Power(rxj, 8LL)*(6655LL + 69LL*Power(rxj, 2LL)) -

                           90LL*Power(rxi, 14LL)*Power(rxj, 2LL)*(1117LL + 77LL*Power(rxj, 2LL)) +

                           Power(rxi, 5LL)*Power(rxj, 12LL)*(6111LL + 92LL*Power(rxj, 2LL)) -

                           18LL*Power(rxi, 6LL)*Power(rxj, 10LL)*(3107LL + 259LL*Power(rxj, 2LL)) +

                           18LL*Power(rxi, 12LL)*Power(rxj, 4LL)*(-31885LL + 403LL*Power(rxj, 2LL))) -

                          exp(2LL*rxi)*Power(rxi, 6LL)*

                          (3465LL*Power(rxi, 12LL)*Power(rxj, 12LL)*

                           (1351350LL + 2483775LL*rxj + 2189250LL*Power(rxj, 2LL) +

           1499400LL*Power(rxj, 3LL) + 512400LL*Power(rxj, 4LL) +

           191940LL*Power(rxj, 5LL) + 73080LL*Power(rxj, 6LL) + 18200LL*Power(rxj, 7LL) +

           2680LL*Power(rxj, 8LL) + 220LL*Power(rxj, 9LL) + 8LL*Power(rxj, 10LL)) -

                           330LL*Power(rxi, 8LL)*Power(rxj, 16LL)*

                           (-2409750LL - 79762725LL*rxj - 9440550LL*Power(rxj, 2LL) -

           6036975LL*Power(rxj, 3LL) - 10098900LL*Power(rxj, 4LL) -

           4800285LL*Power(rxj, 5LL) - 1163190LL*Power(rxj, 6LL) -

           164670LL*Power(rxj, 7LL) - 13110LL*Power(rxj, 8LL) - 365LL*Power(rxj, 9LL) +

           26LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) -

                           2LL*Power(rxj, 24LL)*(1240539300LL + 1516214700LL*rxj +

                                                 891891000LL*Power(rxj, 2LL) + 334459125LL*Power(rxj, 3LL) +

                                                 89189100LL*Power(rxj, 4LL) + 17837820LL*Power(rxj, 5LL) +

                                                 2744280LL*Power(rxj, 6LL) + 326700LL*Power(rxj, 7LL) + 29700LL*Power(rxj, 8LL) +

                                                 1980LL*Power(rxj, 9LL) + 88LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) +

                           Power(rxi, 24LL)*(935550LL + 1715175LL*rxj + 1559250LL*Power(rxj, 2LL) +

                                             935550LL*Power(rxj, 3LL) + 415800LL*Power(rxj, 4LL) + 145530LL*Power(rxj, 5LL) +

                                             41580LL*Power(rxj, 6LL) + 9900LL*Power(rxj, 7LL) + 1980LL*Power(rxj, 8LL) +

                                             330LL*Power(rxj, 9LL) + 44LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           110LL*Power(rxi, 6LL)*Power(rxj, 18LL)*

                           (-313749450LL + 140006475LL*rxj + 40682250LL*Power(rxj, 2LL) -

           63603225LL*Power(rxj, 3LL) - 41107500LL*Power(rxj, 4LL) -

           11688705LL*Power(rxj, 5LL) - 1918350LL*Power(rxj, 6LL) -

           179550LL*Power(rxj, 7LL) - 5670LL*Power(rxj, 8LL) + 735LL*Power(rxj, 9LL) +

           98LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           10LL*Power(rxi, 2LL)*Power(rxj, 22LL)*

                           (-2825672850LL - 2653375725LL*rxj - 1114863750LL*Power(rxj, 2LL) -

           260134875LL*Power(rxj, 3LL) - 29729700LL*Power(rxj, 4LL) +

           1486485LL*Power(rxj, 5LL) + 1295910LL*Power(rxj, 6LL) +

           272250LL*Power(rxj, 7LL) + 34650LL*Power(rxj, 8LL) + 2915LL*Power(rxj, 9LL) +

           154LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           165LL*Power(rxi, 16LL)*Power(rxj, 8LL)*

                           (7739550LL + 14189175LL*rxj + 12899250LL*Power(rxj, 2LL) +

           7739550LL*Power(rxj, 3LL) + 3439800LL*Power(rxj, 4LL) +

           1210860LL*Power(rxj, 5LL) + 330120LL*Power(rxj, 6LL) + 86400LL*Power(rxj, 7LL) +

           18480LL*Power(rxj, 8LL) + 2460LL*Power(rxj, 9LL) + 168LL*Power(rxj, 10LL) +

           4LL*Power(rxj, 11LL)) - 5LL*Power(rxi, 22LL)*Power(rxj, 2LL)*

                           (2806650LL + 5145525LL*rxj + 4677750LL*Power(rxj, 2LL) +

           2806650LL*Power(rxj, 3LL) + 1247400LL*Power(rxj, 4LL) +

           436590LL*Power(rxj, 5LL) + 124740LL*Power(rxj, 6LL) + 29700LL*Power(rxj, 7LL) +

           5940LL*Power(rxj, 8LL) + 990LL*Power(rxj, 9LL) + 132LL*Power(rxj, 10LL) +

           8LL*Power(rxj, 11LL)) - 55LL*Power(rxi, 18LL)*Power(rxj, 6LL)*

                           (7739550LL + 14189175LL*rxj + 12899250LL*Power(rxj, 2LL) +

           7739550LL*Power(rxj, 3LL) + 3439800LL*Power(rxj, 4LL) +

           1203930LL*Power(rxj, 5LL) + 343980LL*Power(rxj, 6LL) + 80820LL*Power(rxj, 7LL) +

           17460LL*Power(rxj, 8LL) + 2790LL*Power(rxj, 9LL) + 244LL*Power(rxj, 10LL) +

           8LL*Power(rxj, 11LL)) - 22LL*Power(rxi, 4LL)*Power(rxj, 20LL)*

                           (2199137850LL + 366522975LL*rxj - 665232750LL*Power(rxj, 2LL) -

           422542575LL*Power(rxj, 3LL) - 123095700LL*Power(rxj, 4LL) -

           20724795LL*Power(rxj, 5LL) - 1838970LL*Power(rxj, 6LL) +

           12150LL*Power(rxj, 7LL) + 26910LL*Power(rxj, 8LL) + 3735LL*Power(rxj, 9LL) +

           258LL*Power(rxj, 10LL) + 8LL*Power(rxj, 11LL)) +

                           33LL*Power(rxi, 10LL)*Power(rxj, 14LL)*

                           (-188215650LL - 280764225LL*rxj - 416886750LL*Power(rxj, 2LL) -

           131922000LL*Power(rxj, 3LL) - 59043600LL*Power(rxj, 4LL) -

           34671420LL*Power(rxj, 5LL) - 11740680LL*Power(rxj, 6LL) -

           2266200LL*Power(rxj, 7LL) - 255000LL*Power(rxj, 8LL) - 15060LL*Power(rxj, 9LL) -

           216LL*Power(rxj, 10LL) + 16LL*Power(rxj, 11LL)) +

                           11LL*Power(rxi, 20LL)*Power(rxj, 4LL)*

                           (8930250LL + 16372125LL*rxj + 14883750LL*Power(rxj, 2LL) +

           8930250LL*Power(rxj, 3LL) + 3969000LL*Power(rxj, 4LL) +

           1389150LL*Power(rxj, 5LL) + 396900LL*Power(rxj, 6LL) + 94500LL*Power(rxj, 7LL) +

           18900LL*Power(rxj, 8LL) + 3290LL*Power(rxj, 9LL) + 364LL*Power(rxj, 10LL) +

           16LL*Power(rxj, 11LL)) - 33LL*Power(rxi, 14LL)*Power(rxj, 10LL)*

                           (85135050LL + 156080925LL*rxj + 141891750LL*Power(rxj, 2LL) +

           84848400LL*Power(rxj, 3LL) + 38984400LL*Power(rxj, 4LL) +

           12157740LL*Power(rxj, 5LL) + 3814440LL*Power(rxj, 6LL) +

           1072200LL*Power(rxj, 7LL) + 198120LL*Power(rxj, 8LL) + 21020LL*Power(rxj, 9LL) +

           1096LL*Power(rxj, 10LL) + 16LL*Power(rxj, 11LL))))/

                         (935550LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 15LL)*Power(rxi + rxj, 15LL))

                         );
        }

    }
    return S;
}


cl_R Slater_6S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_2S_6S(r, xj, xi);
}
