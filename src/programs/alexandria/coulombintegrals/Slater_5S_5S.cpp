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

cl_R Slater_5S_5S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (43191LL*xi)/262144LL

            ;
        }
        else
        {
            S = (1LL/r)*(1LL + (-1LL - (481097LL*rxi)/262144LL - (218953LL*Power(rxi, 2LL))/131072LL -

                                (988003LL*Power(rxi, 3LL))/983040LL - (110459LL*Power(rxi, 4LL))/245760LL -

                                (65243LL*Power(rxi, 5LL))/409600LL - (4769LL*Power(rxi, 6LL))/102400LL -

                                (186229LL*Power(rxi, 7LL))/1.6128e7 - (713LL*Power(rxi, 8LL))/288000LL -

                                (33791LL*Power(rxi, 9LL))/7.2576e7 - (11LL*Power(rxi, 10LL))/141750LL -

                                Power(rxi, 11LL)/86625LL - (2LL*Power(rxi, 12LL))/1.299375e6 -

                                (4LL*Power(rxi, 13LL))/2.1718125e7 - Power(rxi, 14LL)/5.0675625e7 -

                                (2LL*Power(rxi, 15LL))/1.064188125e9 - Power(rxi, 16LL)/6.38512875e9 -

                                Power(rxi, 17LL)/9.0455990625e10 - Power(rxi, 18LL)/1.62820783125e12 -

                                Power(rxi, 19LL)/4.6403923190625e13)/exp(2LL*rxi)

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 18LL) + 19LL*Power(xi, 17LL)*xj + 171LL*Power(xi, 16LL)*Power(xj, 2LL) +

                        969LL*Power(xi, 15LL)*Power(xj, 3LL) + 3876LL*Power(xi, 14LL)*Power(xj, 4LL) +

                        11628LL*Power(xi, 13LL)*Power(xj, 5LL) + 27132LL*Power(xi, 12LL)*Power(xj, 6LL) +

                        50388LL*Power(xi, 11LL)*Power(xj, 7LL) + 75582LL*Power(xi, 10LL)*Power(xj, 8LL) +

                        92378LL*Power(xi, 9LL)*Power(xj, 9LL) + 75582LL*Power(xi, 8LL)*Power(xj, 10LL) +

                        50388LL*Power(xi, 7LL)*Power(xj, 11LL) + 27132LL*Power(xi, 6LL)*Power(xj, 12LL) +

                        11628LL*Power(xi, 5LL)*Power(xj, 13LL) + 3876LL*Power(xi, 4LL)*Power(xj, 14LL) +

                        969LL*Power(xi, 3LL)*Power(xj, 15LL) + 171LL*Power(xi, 2LL)*Power(xj, 16LL) +

                        19LL*xi*Power(xj, 17LL) + Power(xj, 18LL)))/(5LL*Power(xi + xj, 19LL))

            ;
        }
        else
        {
            S = (1LL/r)*((70875LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 19LL) +

                          exp(2LL*rxj)*Power(rxj, 12LL)*

                          (-630LL*Power(rxi, 34LL) - 10LL*Power(rxi, 35LL) + 70875LL*Power(rxj, 26LL) +

                           127575LL*rxi*Power(rxj, 26LL) - 30LL*Power(rxi, 33LL)*(630LL + Power(rxj, 2LL)) +

                           14175LL*Power(rxi, 2LL)*Power(rxj, 24LL)*(-95LL + 8LL*Power(rxj, 2LL)) +

                           4725LL*Power(rxi, 3LL)*Power(rxj, 24LL)*(-513LL + 14LL*Power(rxj, 2LL)) -

                           90LL*Power(rxi, 32LL)*(3920LL + 43LL*Power(rxj, 2LL)) +

                           4725LL*Power(rxi, 5LL)*Power(rxj, 22LL)*

                           (4617LL - 266LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           14175LL*Power(rxi, 4LL)*Power(rxj, 22LL)*

                           (855LL - 152LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           36LL*Power(rxi, 31LL)*(-124950LL - 4985LL*Power(rxj, 2LL) + 13LL*Power(rxj, 4LL)) +

                           36LL*Power(rxi, 30LL)*(-1124550LL - 127960LL*Power(rxj, 2LL) +

                                                  863LL*Power(rxj, 4LL)) + 135LL*Power(rxi, 7LL)*Power(rxj, 20LL)*

                           (-915705LL + 83790LL*Power(rxj, 2LL) - 1330LL*Power(rxj, 4LL) + 4LL*Power(rxj, 6LL))

                           + 315LL*Power(rxi, 6LL)*Power(rxj, 20LL)*

                           (-218025LL + 61560LL*Power(rxj, 2LL) - 1710LL*Power(rxj, 4LL) + 8LL*Power(rxj, 6LL))

                           - 36LL*Power(rxi, 29LL)*(7122150LL + 2102730LL*Power(rxj, 2LL) - 23294LL*Power(rxj, 4LL) +

                                                    37LL*Power(rxj, 6LL)) - 36LL*Power(rxi, 28LL)*

                           (30523500LL + 23401350LL*Power(rxj, 2LL) - 299250LL*Power(rxj, 4LL) +

           1297LL*Power(rxj, 6LL)) +

                           Power(rxi, 17LL)*Power(rxj, 10LL)*

                           (1073961177975LL - 21753487980LL*Power(rxj, 2LL) -

           745994340LL*Power(rxj, 4LL) + 5307156LL*Power(rxj, 6LL) - 818LL*Power(rxj, 8LL))

                           + 10LL*Power(rxi, 9LL)*Power(rxj, 18LL)*

                           (49448070LL - 6409935LL*Power(rxj, 2LL) + 161595LL*Power(rxj, 4LL) -

           1026LL*Power(rxj, 6LL) + Power(rxj, 8LL)) +

                           90LL*Power(rxi, 8LL)*Power(rxj, 18LL)*

                           (3052350LL - 1220940LL*Power(rxj, 2LL) + 53865LL*Power(rxj, 4LL) -

           532LL*Power(rxj, 6LL) + Power(rxj, 8LL)) -

                           1710LL*Power(rxi, 10LL)*Power(rxj, 16LL)*

                           (481950LL - 257040LL*Power(rxj, 2LL) + 16065LL*Power(rxj, 4LL) -

           252LL*Power(rxj, 6LL) + Power(rxj, 8LL)) +

                           6LL*Power(rxi, 11LL)*Power(rxj, 16LL)*

                           (-207559800LL + 50390550LL*Power(rxj, 2LL) - 1165815LL*Power(rxj, 4LL) +

           21396LL*Power(rxj, 6LL) + 5LL*Power(rxj, 8LL)) -

                           18LL*Power(rxi, 13LL)*Power(rxj, 14LL)*

                           (-1703720025LL - 155669850LL*Power(rxj, 2LL) - 7410270LL*Power(rxj, 4LL) -

           1532LL*Power(rxj, 6LL) + 26LL*Power(rxj, 8LL)) +

                           18LL*Power(rxi, 15LL)*Power(rxj, 12LL)*

                           (19380896325LL + 1329128850LL*Power(rxj, 2LL) - 7608930LL*Power(rxj, 4LL) -

           116238LL*Power(rxj, 6LL) + 74LL*Power(rxj, 8LL)) -

                           18LL*Power(rxi, 12LL)*Power(rxj, 14LL)*

                           (89026875LL + 179071200LL*Power(rxj, 2LL) + 1552950LL*Power(rxj, 4LL) +

           295820LL*Power(rxj, 6LL) + 146LL*Power(rxj, 8LL)) +

                           18LL*Power(rxi, 25LL)*Power(rxj, 2LL)*

                           (-5449970925LL - 1137574935LL*Power(rxj, 2LL) + 37834755LL*Power(rxj, 4LL) -

           273062LL*Power(rxj, 6LL) + 171LL*Power(rxj, 8LL)) -

                           9LL*Power(rxi, 19LL)*Power(rxj, 8LL)*

                           (-37914907275LL + 7613889570LL*Power(rxj, 2LL) - 170524620LL*Power(rxj, 4LL) +

           397936LL*Power(rxj, 6LL) + 342LL*Power(rxj, 8LL)) +

                           Power(rxi, 27LL)*(-2884470750LL - 6409935000LL*Power(rxj, 2LL) +

                                             28332990LL*Power(rxj, 4LL) + 58104LL*Power(rxj, 6LL) + 818LL*Power(rxj, 8LL)) -

                           3LL*Power(rxi, 23LL)*Power(rxj, 4LL)*

                           (219130630425LL - 11118046590LL*Power(rxj, 2LL) +

           327611970LL*Power(rxj, 4LL) - 2920908LL*Power(rxj, 6LL) + 2584LL*Power(rxj, 8LL)

                           ) + 3LL*Power(rxi, 21LL)*Power(rxj, 6LL)*

                           (-345162539925LL + 19030764690LL*Power(rxj, 2LL) -

           141976170LL*Power(rxj, 4LL) - 1441872LL*Power(rxj, 6LL) + 2584LL*Power(rxj, 8LL)

                           ) + 63LL*Power(rxi, 20LL)*Power(rxj, 6LL)*

                           (-50980542525LL + 6240202920LL*Power(rxj, 2LL) - 201314310LL*Power(rxj, 4LL) +

           956080LL*Power(rxj, 6LL) + 2584LL*Power(rxj, 8LL)) +

                           18LL*Power(rxi, 14LL)*Power(rxj, 12LL)*

                           (-7803332775LL - 2519206200LL*Power(rxj, 2LL) - 119719950LL*Power(rxj, 4LL) +

           182280LL*Power(rxj, 6LL) + 2734LL*Power(rxj, 8LL)) -

                           18LL*Power(rxi, 26LL)*(195859125LL + 1794781800LL*Power(rxj, 2LL) +

                                                  67337235LL*Power(rxj, 4LL) - 1659700LL*Power(rxj, 6LL) + 4089LL*Power(rxj, 8LL))

                           + 9LL*Power(rxi, 18LL)*Power(rxj, 8LL)*

                           (-357591274425LL + 8328390840LL*Power(rxj, 2LL) +

           912042180LL*Power(rxj, 4LL) - 12842480LL*Power(rxj, 6LL) +

           10678LL*Power(rxj, 8LL)) -

                           9LL*Power(rxi, 16LL)*Power(rxj, 10LL)*

                           (128599724925LL + 21298077360LL*Power(rxj, 2LL) -

           267928500LL*Power(rxj, 4LL) - 5458320LL*Power(rxj, 6LL) +

           14722LL*Power(rxj, 8LL)) +

                           18LL*Power(rxi, 24LL)*Power(rxj, 2LL)*

                           (-7604930025LL - 8866107180LL*Power(rxj, 2LL) + 399272265LL*Power(rxj, 4LL) -

           5925780LL*Power(rxj, 6LL) + 17651LL*Power(rxj, 8LL)) -

                           9LL*Power(rxi, 22LL)*Power(rxj, 4LL)*

                           (129194933175LL + 3909863160LL*Power(rxj, 2LL) + 91420770LL*Power(rxj, 4LL) -

           8762040LL*Power(rxj, 6LL) + 43928LL*Power(rxj, 8LL))) +

                          exp(2LL*rxi)*Power(rxi, 12LL)*

                          (Power(rxi, 8LL)*Power(rxj, 18LL)*

                           (3218321469825LL - 341234165475LL*rxj - 393132783960LL*Power(rxj, 2LL) -

           57092294070LL*Power(rxj, 3LL) + 822786930LL*Power(rxj, 4LL) +

           982835910LL*Power(rxj, 5LL) + 106664040LL*Power(rxj, 6LL) +

           4915116LL*Power(rxj, 7LL) + 73602LL*Power(rxj, 8LL) - 818LL*Power(rxj, 9LL)) +

                           10LL*Power(rxj, 26LL)*(352546425LL + 288447075LL*rxj +

                                                  109884600LL*Power(rxj, 2LL) + 25639740LL*Power(rxj, 3LL) +

                                                  4048380LL*Power(rxj, 4LL) + 449820LL*Power(rxj, 5LL) + 35280LL*Power(rxj, 6LL) +

                                                  1890LL*Power(rxj, 7LL) + 63LL*Power(rxj, 8LL) + Power(rxj, 9LL)) +

                           30LL*Power(rxi, 2LL)*Power(rxj, 24LL)*

                           (4562958015LL + 3269982555LL*rxj + 1076869080LL*Power(rxj, 2LL) +

           213664500LL*Power(rxj, 3LL) + 28081620LL*Power(rxj, 4LL) +

           2523276LL*Power(rxj, 5LL) + 153552LL*Power(rxj, 6LL) + 5982LL*Power(rxj, 7LL) +

           129LL*Power(rxj, 8LL) + Power(rxj, 9LL)) -

                           15LL*Power(rxi, 24LL)*Power(rxj, 2LL)*

                           (-89775LL - 161595LL*rxj - 143640LL*Power(rxj, 2LL) - 83790LL*Power(rxj, 3LL) -

           35910LL*Power(rxj, 4LL) - 11970LL*Power(rxj, 5LL) - 3192LL*Power(rxj, 6LL) -

           684LL*Power(rxj, 7LL) - 114LL*Power(rxj, 8LL) + 2LL*Power(rxj, 9LL)) -

                           5LL*Power(rxi, 26LL)*(14175LL + 25515LL*rxj + 22680LL*Power(rxj, 2LL) +

                                                 13230LL*Power(rxj, 3LL) + 5670LL*Power(rxj, 4LL) + 1890LL*Power(rxj, 5LL) +

                                                 504LL*Power(rxj, 6LL) + 108LL*Power(rxj, 7LL) + 18LL*Power(rxj, 8LL) +

                                                 2LL*Power(rxj, 9LL)) - 1938LL*Power(rxi, 14LL)*Power(rxj, 12LL)*

                           (-826875LL + 15824025LL*rxj - 23398200LL*Power(rxj, 2LL) +

           12344850LL*Power(rxj, 3LL) + 1244250LL*Power(rxj, 4LL) -

           384930LL*Power(rxj, 5LL) - 59640LL*Power(rxj, 6LL) - 1848LL*Power(rxj, 7LL) +

           84LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) +

                           1938LL*Power(rxi, 12LL)*Power(rxj, 14LL)*

                           (72476775LL - 180008325LL*rxj + 98907480LL*Power(rxj, 2LL) +

           11224710LL*Power(rxj, 3LL) - 4235490LL*Power(rxj, 4LL) -

           791910LL*Power(rxj, 5LL) - 31080LL*Power(rxj, 6LL) + 2232LL*Power(rxj, 7LL) +

           204LL*Power(rxj, 8LL) + 4LL*Power(rxj, 9LL)) +

                           342LL*Power(rxi, 16LL)*Power(rxj, 10LL)*

                           (2409750LL + 3641400LL*rxj + 9424800LL*Power(rxj, 2LL) -

           8193150LL*Power(rxj, 3LL) + 6301050LL*Power(rxj, 4LL) +

           400470LL*Power(rxj, 5LL) - 143640LL*Power(rxj, 6LL) - 15518LL*Power(rxj, 7LL) -

           281LL*Power(rxj, 8LL) + 9LL*Power(rxj, 9LL)) -

                           171LL*Power(rxi, 10LL)*Power(rxj, 16LL)*

                           (-6768406575LL + 6280474725LL*rxj + 438336360LL*Power(rxj, 2LL) -

           400731030LL*Power(rxj, 3LL) - 74168430LL*Power(rxj, 4LL) -

           2490810LL*Power(rxj, 5LL) + 461160LL*Power(rxj, 6LL) + 51244LL*Power(rxj, 7LL) +

           1858LL*Power(rxj, 8LL) + 18LL*Power(rxj, 9LL)) +

                           9LL*Power(rxi, 22LL)*Power(rxj, 4LL)*

                           (-1346625LL - 2423925LL*rxj - 2154600LL*Power(rxj, 2LL) -

           1256850LL*Power(rxj, 3LL) - 538650LL*Power(rxj, 4LL) -

           179550LL*Power(rxj, 5LL) - 47880LL*Power(rxj, 6LL) - 14264LL*Power(rxj, 7LL) +

           292LL*Power(rxj, 8LL) + 52LL*Power(rxj, 9LL)) -

                           9LL*Power(rxi, 4LL)*Power(rxj, 22LL)*

                           (-129194933175LL - 73043543475LL*rxj - 17732214360LL*Power(rxj, 2LL) -

           2275149870LL*Power(rxj, 3LL) - 134674470LL*Power(rxj, 4LL) +

           3148110LL*Power(rxj, 5LL) + 1197000LL*Power(rxj, 6LL) +

           93176LL*Power(rxj, 7LL) + 3452LL*Power(rxj, 8LL) + 52LL*Power(rxj, 9LL)) +

                           9LL*Power(rxi, 6LL)*Power(rxj, 20LL)*

                           (356863797675LL + 115054179975LL*rxj + 3909863160LL*Power(rxj, 2LL) -

           3706015530LL*Power(rxj, 3LL) - 798544530LL*Power(rxj, 4LL) -

           75669510LL*Power(rxj, 5LL) - 3319400LL*Power(rxj, 6LL) -

           6456LL*Power(rxj, 7LL) + 5188LL*Power(rxj, 8LL) + 148LL*Power(rxj, 9LL)) -

                           9LL*Power(rxi, 20LL)*Power(rxj, 6LL)*

                           (-7630875LL - 13735575LL*rxj - 12209400LL*Power(rxj, 2LL) -

           7122150LL*Power(rxj, 3LL) - 3052350LL*Power(rxj, 4LL) -

           777210LL*Power(rxj, 5LL) - 591640LL*Power(rxj, 6LL) + 3064LL*Power(rxj, 7LL) +

           5468LL*Power(rxj, 8LL) + 148LL*Power(rxj, 9LL)) +

                           2LL*Power(rxi, 18LL)*Power(rxj, 8LL)*

                           (-137355750LL - 247240350LL*rxj - 219769200LL*Power(rxj, 2LL) -

           151171650LL*Power(rxj, 3LL) + 13976550LL*Power(rxj, 4LL) -

           66692430LL*Power(rxj, 5LL) - 1640520LL*Power(rxj, 6LL) +

           1046142LL*Power(rxj, 7LL) + 66249LL*Power(rxj, 8LL) + 409LL*Power(rxj, 9LL))))/

                         (70875LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 19LL)*Power(rxi + rxj, 19LL))

                         );
        }

    }
    return S;
}
