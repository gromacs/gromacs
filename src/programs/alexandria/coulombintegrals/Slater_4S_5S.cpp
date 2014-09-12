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

cl_R Slater_4S_5S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (234137LL*xi)/1.31072e6

            ;
        }
        else
        {
            S = (1LL/r)*((-14227497123840000LL + 14227497123840000LL*exp(2LL*rxi) -

                          25913502934444125LL*rxi - 23372011621208250LL*Power(rxi, 2LL) -

                          13907709869303250LL*Power(rxi, 3LL) - 6137735659555500LL*Power(rxi, 4LL) -

                          2140857388870200LL*Power(rxi, 5LL) - 614116575072000LL*Power(rxi, 6LL) -

                          148809580920000LL*Power(rxi, 7LL) - 31036639233600LL*Power(rxi, 8LL) -

                          5645342102400LL*Power(rxi, 9LL) - 903333150720LL*Power(rxi, 10LL) -

                          127744081920LL*Power(rxi, 11LL) - 15968010240LL*Power(rxi, 12LL) -

                          1754726400LL*Power(rxi, 13LL) - 167116800LL*Power(rxi, 14LL) -

                          13369344LL*Power(rxi, 15LL) - 835584LL*Power(rxi, 16LL) - 32768LL*Power(rxi, 17LL))/

                         (1.422749712384e16*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(4LL*Power(xi, 16LL) + 68LL*Power(xi, 15LL)*xj + 544LL*Power(xi, 14LL)*Power(xj, 2LL) +

                        2720LL*Power(xi, 13LL)*Power(xj, 3LL) + 9520LL*Power(xi, 12LL)*Power(xj, 4LL) +

                        24752LL*Power(xi, 11LL)*Power(xj, 5LL) + 49504LL*Power(xi, 10LL)*Power(xj, 6LL) +

                        77792LL*Power(xi, 9LL)*Power(xj, 7LL) + 97240LL*Power(xi, 8LL)*Power(xj, 8LL) +

                        97240LL*Power(xi, 7LL)*Power(xj, 9LL) + 61880LL*Power(xi, 6LL)*Power(xj, 10LL) +

                        30940LL*Power(xi, 5LL)*Power(xj, 11LL) + 11900LL*Power(xi, 4LL)*Power(xj, 12LL) +

                        3400LL*Power(xi, 3LL)*Power(xj, 13LL) + 680LL*Power(xi, 2LL)*Power(xj, 14LL) +

                        85LL*xi*Power(xj, 15LL) + 5LL*Power(xj, 16LL)))/(20LL*Power(xi + xj, 17LL))

            ;
        }
        else
        {
            S = (1LL/r)*((56700LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 17LL) +

                          9LL*exp(2LL*rxj)*Power(rxj, 12LL)*

                          (-980LL*Power(rxi, 28LL) - 20LL*Power(rxi, 29LL) + 6300LL*Power(rxj, 22LL) +

                           11025LL*rxi*Power(rxj, 22LL) - 50LL*Power(rxi, 27LL)*(441LL + 2LL*Power(rxj, 2LL)) +

                           3150LL*Power(rxi, 2LL)*Power(rxj, 20LL)*(-34LL + 3LL*Power(rxj, 2LL)) +

                           525LL*Power(rxi, 3LL)*Power(rxj, 20LL)*(-357LL + 10LL*Power(rxj, 2LL)) -

                           420LL*Power(rxi, 26LL)*(700LL + 19LL*Power(rxj, 2LL)) +

                           1050LL*Power(rxi, 4LL)*Power(rxj, 18LL)*

                           (816LL - 153LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           210LL*Power(rxi, 5LL)*Power(rxj, 18LL)*

                           (7140LL - 425LL*Power(rxj, 2LL) + 3LL*Power(rxj, 4LL)) +

                           42LL*Power(rxi, 25LL)*(-59500LL - 6035LL*Power(rxj, 2LL) + 18LL*Power(rxj, 4LL)) +

                           84LL*Power(rxi, 24LL)*(-160650LL - 52700LL*Power(rxj, 2LL) + 397LL*Power(rxj, 4LL)) -

                           28LL*Power(rxi, 12LL)*Power(rxj, 10LL)*

                           (100849950LL + 27100125LL*Power(rxj, 2LL) + 186150LL*Power(rxj, 4LL) -

           2177LL*Power(rxj, 6LL)) +

                           140LL*Power(rxi, 6LL)*Power(rxj, 16LL)*

                           (-30600LL + 9180LL*Power(rxj, 2LL) - 255LL*Power(rxj, 4LL) + Power(rxj, 6LL)) -

                           2380LL*Power(rxi, 8LL)*Power(rxj, 14LL)*

                           (-6300LL + 2700LL*Power(rxj, 2LL) - 120LL*Power(rxj, 4LL) + Power(rxj, 6LL)) +

                           10LL*Power(rxi, 7LL)*Power(rxj, 16LL)*

                           (-749700LL + 71400LL*Power(rxj, 2LL) - 1071LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL))

                           + 204LL*Power(rxi, 15LL)*Power(rxj, 8LL)*

                           (28962255LL - 1744750LL*Power(rxj, 2LL) + 9555LL*Power(rxj, 4LL) +

           6LL*Power(rxj, 6LL)) - 42LL*Power(rxi, 11LL)*Power(rxj, 12LL)*

                           (-12911925LL - 1634550LL*Power(rxj, 2LL) - 7103LL*Power(rxj, 4LL) +

           18LL*Power(rxj, 6LL)) + 2LL*Power(rxi, 9LL)*Power(rxj, 14LL)*

                           (16948575LL - 1184400LL*Power(rxj, 2LL) + 63861LL*Power(rxj, 4LL) +

           50LL*Power(rxj, 6LL)) + 28LL*Power(rxi, 22LL)*

                           (-2180250LL - 10993050LL*Power(rxj, 2LL) + 14925LL*Power(rxj, 4LL) +

           73LL*Power(rxj, 6LL)) - 952LL*Power(rxi, 14LL)*Power(rxj, 8LL)*

                           (16966215LL + 725175LL*Power(rxj, 2LL) - 36075LL*Power(rxj, 4LL) +

           79LL*Power(rxj, 6LL)) - 84LL*Power(rxi, 10LL)*Power(rxj, 12LL)*

                           (1723800LL + 279225LL*Power(rxj, 2LL) + 45600LL*Power(rxj, 4LL) +

           107LL*Power(rxj, 6LL)) - 35LL*Power(rxi, 17LL)*Power(rxj, 6LL)*

                           (132637869LL - 2205240LL*Power(rxj, 2LL) - 48348LL*Power(rxj, 4LL) +

           136LL*Power(rxj, 6LL)) - 6LL*Power(rxi, 21LL)*Power(rxj, 2LL)*

                           (192298050LL + 12644275LL*Power(rxj, 2LL) - 218029LL*Power(rxj, 4LL) +

           204LL*Power(rxj, 6LL)) + 4LL*Power(rxi, 13LL)*Power(rxj, 10LL)*

                           (1259522775LL + 15895425LL*Power(rxj, 2LL) - 493017LL*Power(rxj, 4LL) +

           263LL*Power(rxj, 6LL)) - 140LL*Power(rxi, 16LL)*Power(rxj, 6LL)*

                           (180826281LL - 15101406LL*Power(rxj, 2LL) + 160140LL*Power(rxj, 4LL) +

           442LL*Power(rxj, 6LL)) - 2LL*Power(rxi, 23LL)*

                           (21366450LL + 23526300LL*Power(rxj, 2LL) - 246729LL*Power(rxj, 4LL) +

           526LL*Power(rxj, 6LL)) + 7LL*Power(rxi, 19LL)*Power(rxj, 4LL)*

                           (-811081215LL + 39095550LL*Power(rxj, 2LL) - 515916LL*Power(rxj, 4LL) +

           680LL*Power(rxj, 6LL)) + 70LL*Power(rxi, 18LL)*Power(rxj, 4LL)*

                           (-180554454LL + 9873711LL*Power(rxj, 2LL) - 414120LL*Power(rxj, 4LL) +

           2924LL*Power(rxj, 6LL)) -

                           14LL*Power(rxi, 20LL)*Power(rxj, 2LL)*

                           (136919700LL + 71867115LL*Power(rxj, 2LL) - 2154150LL*Power(rxj, 4LL) +

           10268LL*Power(rxj, 6LL))) -

                          4LL*exp(2LL*rxi)*Power(rxi, 10LL)*

                          (-10710LL*Power(rxi, 12LL)*Power(rxj, 12LL)*

                           (-3555LL - 127008LL*rxj + 138384LL*Power(rxj, 2LL) - 74556LL*Power(rxj, 3LL) -

           22284LL*Power(rxj, 4LL) + 408LL*Power(rxj, 5LL) + 576LL*Power(rxj, 6LL) +

           60LL*Power(rxj, 7LL) + 2LL*Power(rxj, 8LL)) +

                           2LL*Power(rxi, 20LL)*Power(rxj, 4LL)*

                           (963900LL + 1735020LL*rxj + 1542240LL*Power(rxj, 2LL) +

           899640LL*Power(rxj, 3LL) + 385560LL*Power(rxj, 4LL) + 128520LL*Power(rxj, 5LL) +

           34272LL*Power(rxj, 6LL) + 9126LL*Power(rxj, 7LL) + 333LL*Power(rxj, 8LL) -

           20LL*Power(rxj, 9LL)) - 2LL*Power(rxj, 24LL)*

                           (119041650LL + 107137485LL*rxj + 45110520LL*Power(rxj, 2LL) +

           11695320LL*Power(rxj, 3LL) + 2063880LL*Power(rxj, 4LL) +

           257985LL*Power(rxj, 5LL) + 22932LL*Power(rxj, 6LL) + 1404LL*Power(rxj, 7LL) +

           54LL*Power(rxj, 8LL) + Power(rxj, 9LL)) +

                           2LL*Power(rxi, 2LL)*Power(rxj, 22LL)*

                           (-3264488325LL - 2505368880LL*rxj - 881390160LL*Power(rxj, 2LL) -

           185775660LL*Power(rxj, 3LL) - 25639740LL*Power(rxj, 4LL) -

           2361555LL*Power(rxj, 5LL) - 139356LL*Power(rxj, 6LL) - 4482LL*Power(rxj, 7LL) -

           27LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) +

                           Power(rxi, 24LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj, 2LL) +

                                             13230LL*Power(rxj, 3LL) + 5670LL*Power(rxj, 4LL) + 1890LL*Power(rxj, 5LL) +

                                             504LL*Power(rxj, 6LL) + 108LL*Power(rxj, 7LL) + 18LL*Power(rxj, 8LL) +

                                             2LL*Power(rxj, 9LL)) - 102LL*Power(rxi, 10LL)*Power(rxj, 14LL)*

                           (44986725LL - 97433280LL*rxj + 44467920LL*Power(rxj, 2LL) +

           15857100LL*Power(rxj, 3LL) - 457380LL*Power(rxj, 4LL) -

           620550LL*Power(rxj, 5LL) - 83160LL*Power(rxj, 6LL) - 4068LL*Power(rxj, 7LL) -

           6LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) +

                           102LL*Power(rxi, 14LL)*Power(rxj, 10LL)*

                           (-859950LL - 1437345LL*rxj - 2260440LL*Power(rxj, 2LL) +

           810810LL*Power(rxj, 3LL) - 1056510LL*Power(rxj, 4LL) -

           217854LL*Power(rxj, 5LL) + 6552LL*Power(rxj, 6LL) + 3852LL*Power(rxj, 7LL) +

           258LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) -

                           Power(rxi, 22LL)*Power(rxj, 2LL)*

                           (240975LL + 433755LL*rxj + 385560LL*Power(rxj, 2LL) + 224910LL*Power(rxj, 3LL) +

           96390LL*Power(rxj, 4LL) + 32130LL*Power(rxj, 5LL) + 8568LL*Power(rxj, 6LL) +

           1836LL*Power(rxj, 7LL) + 306LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) +

                           2LL*Power(rxi, 4LL)*Power(rxj, 20LL)*

                           (-18032978565LL - 9823683240LL*rxj - 2047323600LL*Power(rxj, 2LL) -

           129098340LL*Power(rxj, 3LL) + 26410860LL*Power(rxj, 4LL) +

           7094304LL*Power(rxj, 5LL) + 788256LL*Power(rxj, 6LL) + 48654LL*Power(rxj, 7LL) +

           1593LL*Power(rxj, 8LL) + 20LL*Power(rxj, 9LL)) -

                           6LL*Power(rxi, 16LL)*Power(rxj, 8LL)*

                           (-5622750LL - 10120950LL*rxj - 8996400LL*Power(rxj, 2LL) -

           5698350LL*Power(rxj, 3LL) - 897750LL*Power(rxj, 4LL) -

           1641591LL*Power(rxj, 5LL) - 211932LL*Power(rxj, 6LL) + 10224LL*Power(rxj, 7LL) +

           2364LL*Power(rxj, 8LL) + 73LL*Power(rxj, 9LL)) +

                           2LL*Power(rxi, 18LL)*Power(rxj, 6LL)*

                           (-4819500LL - 8675100LL*rxj - 7711200LL*Power(rxj, 2LL) -

           4498200LL*Power(rxj, 3LL) - 1927800LL*Power(rxj, 4LL) -

           561519LL*Power(rxj, 5LL) - 279468LL*Power(rxj, 6LL) - 20682LL*Power(rxj, 7LL) +

           1305LL*Power(rxj, 8LL) + 106LL*Power(rxj, 9LL)) +

                           3LL*Power(rxi, 8LL)*Power(rxj, 16LL)*

                           (-9364244085LL + 6940428705LL*rxj + 2117684520LL*Power(rxj, 2LL) -

           230268150LL*Power(rxj, 3LL) - 149610510LL*Power(rxj, 4LL) -

           21824334LL*Power(rxj, 5LL) - 1223208LL*Power(rxj, 6LL) +

           12708LL*Power(rxj, 7LL) + 4470LL*Power(rxj, 8LL) + 146LL*Power(rxj, 9LL)) -

                           Power(rxi, 6LL)*Power(rxj, 18LL)*

                           (57304872765LL + 7147185255LL*rxj - 5801702760LL*Power(rxj, 2LL) -

           2053388610LL*Power(rxj, 3LL) - 271655370LL*Power(rxj, 4LL) -

           10864854LL*Power(rxj, 5LL) + 1337112LL*Power(rxj, 6LL) +

           202716LL*Power(rxj, 7LL) + 10746LL*Power(rxj, 8LL) + 212LL*Power(rxj, 9LL))))/

                         (56700LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 17LL)*Power(rxi + rxj, 17LL))

                         );
        }

    }
    return S;
}


cl_R Slater_5S_4S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_4S_5S(r, xj, xi);
}
