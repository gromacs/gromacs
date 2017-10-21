/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifdef HAVE_LIBCLN

cl_R Slater_1S_6S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (341LL*xi)/2048LL

            ;
        }
        else
        {
            S = (1LL/r)*((-74724249600LL + 74724249600LL*exp(2LL*rxi) - 137006619750LL*rxi -

                          124564740300LL*Power(rxi, 2LL) - 74754654975LL*Power(rxi, 3LL) -

                          33239155950LL*Power(rxi, 4LL) - 11644853220LL*Power(rxi, 5LL) -

                          3334050720LL*Power(rxi, 6LL) - 797528160LL*Power(rxi, 7LL) -

                          161235360LL*Power(rxi, 8LL) - 27593280LL*Power(rxi, 9LL) - 3953664LL*Power(rxi, 10LL) -

                          459264LL*Power(rxi, 11LL) - 39936LL*Power(rxi, 12LL) - 2048LL*Power(rxi, 13LL))/

                         (7.47242496e10*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 12LL) + 13LL*Power(xi, 11LL)*xj + 78LL*Power(xi, 10LL)*Power(xj, 2LL) +

                        286LL*Power(xi, 9LL)*Power(xj, 3LL) + 715LL*Power(xi, 8LL)*Power(xj, 4LL) +

                        1287LL*Power(xi, 7LL)*Power(xj, 5LL) + 1716LL*Power(xi, 6LL)*Power(xj, 6LL) +

                        1716LL*Power(xi, 5LL)*Power(xj, 7LL) + 1287LL*Power(xi, 4LL)*Power(xj, 8LL) +

                        715LL*Power(xi, 3LL)*Power(xj, 9LL) + 286LL*Power(xi, 2LL)*Power(xj, 10LL) +

                        78LL*xi*Power(xj, 11LL) + 6LL*Power(xj, 12LL)))/(6LL*Power(xi + xj, 13LL))

            ;
        }
        else
        {
            S = (1LL/r)*((935550LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 13LL) +

                          311850LL*exp(2LL*rxj)*Power(rxj, 14LL)*

                          (-24LL*Power(rxi, 12LL) - 3LL*Power(rxi, 13LL) -

                           507LL*Power(rxi, 10LL)*Power(rxj, 2LL) - 52LL*Power(rxi, 11LL)*Power(rxj, 2LL) -

                           2145LL*Power(rxi, 8LL)*Power(rxj, 4LL) - 143LL*Power(rxi, 9LL)*Power(rxj, 4LL) -

                           2574LL*Power(rxi, 6LL)*Power(rxj, 6LL) - 858LL*Power(rxi, 4LL)*Power(rxj, 8LL) +

                           143LL*Power(rxi, 5LL)*Power(rxj, 8LL) - 39LL*Power(rxi, 2LL)*Power(rxj, 10LL) +

                           52LL*Power(rxi, 3LL)*Power(rxj, 10LL) + 3LL*Power(rxj, 12LL) + 3LL*rxi*Power(rxj, 12LL)

                          ) - exp(2LL*rxi)*Power(rxi, 4LL)*

                          (110LL*Power(rxi, 18LL)*Power(rxj, 4LL)*

                           (663390LL + 1216215LL*rxj + 1105650LL*Power(rxj, 2LL) +

                            663390LL*Power(rxj, 3LL) + 294840LL*Power(rxj, 4LL) + 103194LL*Power(rxj, 5LL) +

                            29484LL*Power(rxj, 6LL) + 7020LL*Power(rxj, 7LL) + 1404LL*Power(rxj, 8LL) +

                            237LL*Power(rxj, 9LL) + 30LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) -

                           330LL*Power(rxi, 16LL)*Power(rxj, 6LL)*

                           (810810LL + 1486485LL*rxj + 1351350LL*Power(rxj, 2LL) +

                            810810LL*Power(rxj, 3LL) + 360360LL*Power(rxj, 4LL) + 126126LL*Power(rxj, 5LL) +

                            36036LL*Power(rxj, 6LL) + 8556LL*Power(rxj, 7LL) + 1740LL*Power(rxj, 8LL) +

                            291LL*Power(rxj, 9LL) + 34LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) +

                           330LL*Power(rxi, 6LL)*Power(rxj, 16LL)*

                           (3169530LL + 7960680LL*rxj + 5798520LL*Power(rxj, 2LL) +

                            3144960LL*Power(rxj, 3LL) + 1572480LL*Power(rxj, 4LL) +

                            638001LL*Power(rxj, 5LL) + 191646LL*Power(rxj, 6LL) + 41886LL*Power(rxj, 7LL) +

                            6630LL*Power(rxj, 8LL) + 741LL*Power(rxj, 9LL) + 54LL*Power(rxj, 10LL) +

                            2LL*Power(rxj, 11LL)) - 110LL*Power(rxi, 4LL)*Power(rxj, 18LL)*

                           (12162150LL + 8108100LL*rxj + 6486480LL*Power(rxj, 2LL) +

                            5675670LL*Power(rxj, 3LL) + 3243240LL*Power(rxj, 4LL) +

                            1216215LL*Power(rxj, 5LL) + 319410LL*Power(rxj, 6LL) + 61074LL*Power(rxj, 7LL) +

                            8586LL*Power(rxj, 8LL) + 867LL*Power(rxj, 9LL) + 58LL*Power(rxj, 10LL) +

                            2LL*Power(rxj, 11LL)) + Power(rxi, 22LL)*

                           (935550LL + 1715175LL*rxj + 1559250LL*Power(rxj, 2LL) +

                            935550LL*Power(rxj, 3LL) + 415800LL*Power(rxj, 4LL) + 145530LL*Power(rxj, 5LL) +

                            41580LL*Power(rxj, 6LL) + 9900LL*Power(rxj, 7LL) + 1980LL*Power(rxj, 8LL) +

                            330LL*Power(rxj, 9LL) + 44LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           11LL*Power(rxi, 20LL)*Power(rxj, 2LL)*

                           (1105650LL + 2027025LL*rxj + 1842750LL*Power(rxj, 2LL) +

                            1105650LL*Power(rxj, 3LL) + 491400LL*Power(rxj, 4LL) +

                            171990LL*Power(rxj, 5LL) + 49140LL*Power(rxj, 6LL) + 11700LL*Power(rxj, 7LL) +

                            2340LL*Power(rxj, 8LL) + 390LL*Power(rxj, 9LL) + 52LL*Power(rxj, 10LL) +

                            4LL*Power(rxj, 11LL)) + 11LL*Power(rxi, 2LL)*Power(rxj, 20LL)*

                           (-48648600LL + 2027025LL*rxj + 44594550LL*Power(rxj, 2LL) +

                            36486450LL*Power(rxj, 3LL) + 16216200LL*Power(rxj, 4LL) +

                            4864860LL*Power(rxj, 5LL) + 1065960LL*Power(rxj, 6LL) +

                            176040LL*Power(rxj, 7LL) + 21960LL*Power(rxj, 8LL) + 2010LL*Power(rxj, 9LL) +

                            124LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           Power(rxj, 22LL)*(340540200LL + 468242775LL*rxj + 312161850LL*Power(rxj, 2LL) +

                                             133783650LL*Power(rxj, 3LL) + 41164200LL*Power(rxj, 4LL) +

                                             9604980LL*Power(rxj, 5LL) + 1746360LL*Power(rxj, 6LL) +

                                             249480LL*Power(rxj, 7LL) + 27720LL*Power(rxj, 8LL) + 2310LL*Power(rxj, 9LL) +

                                             132LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           165LL*Power(rxi, 14LL)*Power(rxj, 8LL)*

                           (4054050LL + 7432425LL*rxj + 6756750LL*Power(rxj, 2LL) +

                            4054050LL*Power(rxj, 3LL) + 1801800LL*Power(rxj, 4LL) +

                            631260LL*Power(rxj, 5LL) + 178920LL*Power(rxj, 6LL) + 43176LL*Power(rxj, 7LL) +

                            8904LL*Power(rxj, 8LL) + 1428LL*Power(rxj, 9LL) + 152LL*Power(rxj, 10LL) +

                            8LL*Power(rxj, 11LL)) - 231LL*Power(rxi, 12LL)*Power(rxj, 10LL)*

                           (5212350LL + 9555975LL*rxj + 8687250LL*Power(rxj, 2LL) +

                            5209650LL*Power(rxj, 3LL) + 2327400LL*Power(rxj, 4LL) +

                            801540LL*Power(rxj, 5LL) + 230040LL*Power(rxj, 6LL) + 57240LL*Power(rxj, 7LL) +

                            11640LL*Power(rxj, 8LL) + 1740LL*Power(rxj, 9LL) + 168LL*Power(rxj, 10LL) +

                            8LL*Power(rxj, 11LL)) + 231LL*Power(rxi, 10LL)*Power(rxj, 12LL)*

                           (6949800LL + 12746025LL*rxj + 11535750LL*Power(rxj, 2LL) +

                            7056450LL*Power(rxj, 3LL) + 3040200LL*Power(rxj, 4LL) +

                            1051920LL*Power(rxj, 5LL) + 316800LL*Power(rxj, 6LL) + 79680LL*Power(rxj, 7LL) +

                            15360LL*Power(rxj, 8LL) + 2100LL*Power(rxj, 9LL) + 184LL*Power(rxj, 10LL) +

                            8LL*Power(rxj, 11LL)) - 165LL*Power(rxi, 8LL)*Power(rxj, 14LL)*

                           (9775080LL + 17424855LL*rxj + 17019450LL*Power(rxj, 2LL) +

                            9519930LL*Power(rxj, 3LL) + 4059720LL*Power(rxj, 4LL) +

                            1519056LL*Power(rxj, 5LL) + 475776LL*Power(rxj, 6LL) + 114720LL*Power(rxj, 7LL) +

                            20256LL*Power(rxj, 8LL) + 2508LL*Power(rxj, 9LL) + 200LL*Power(rxj, 10LL) +

                            8LL*Power(rxj, 11LL))))/

                         (935550LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 13LL)*Power(rxi + rxj, 13LL))

                         );
        }

    }
    return S;
}


cl_R Slater_6S_1S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_1S_6S(r, xj, xi);
}

#else

double Slater_1S_6S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (341*xi)/2048

            ;
        }
        else
        {
            S = (1/r)*((-74724249600 + 74724249600*exp(2*rxi) - 137006619750*rxi -

                          124564740300*power(rxi, 2) - 74754654975*power(rxi, 3) -

                          33239155950*power(rxi, 4) - 11644853220*power(rxi, 5) -

                          3334050720*power(rxi, 6) - 797528160*power(rxi, 7) -

                          161235360*power(rxi, 8) - 27593280*power(rxi, 9) - 3953664*power(rxi, 10) -

                          459264*power(rxi, 11) - 39936*power(rxi, 12) - 2048*power(rxi, 13))/

                         (7.47242496e10*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(power(xi, 12) + 13*power(xi, 11)*xj + 78*power(xi, 10)*power(xj, 2) +

                        286*power(xi, 9)*power(xj, 3) + 715*power(xi, 8)*power(xj, 4) +

                        1287*power(xi, 7)*power(xj, 5) + 1716*power(xi, 6)*power(xj, 6) +

                        1716*power(xi, 5)*power(xj, 7) + 1287*power(xi, 4)*power(xj, 8) +

                        715*power(xi, 3)*power(xj, 9) + 286*power(xi, 2)*power(xj, 10) +

                        78*xi*power(xj, 11) + 6*power(xj, 12)))/(6*power(xi + xj, 13))

            ;
        }
        else
        {
            S = (1/r)*((935550*exp(2*(rxi + rxj))*power(power(rxi, 2) - power(rxj, 2), 13) +

                          311850*exp(2*rxj)*power(rxj, 14)*

                          (-24*power(rxi, 12) - 3*power(rxi, 13) -

                           507*power(rxi, 10)*power(rxj, 2) - 52*power(rxi, 11)*power(rxj, 2) -

                           2145*power(rxi, 8)*power(rxj, 4) - 143*power(rxi, 9)*power(rxj, 4) -

                           2574*power(rxi, 6)*power(rxj, 6) - 858*power(rxi, 4)*power(rxj, 8) +

                           143*power(rxi, 5)*power(rxj, 8) - 39*power(rxi, 2)*power(rxj, 10) +

                           52*power(rxi, 3)*power(rxj, 10) + 3*power(rxj, 12) + 3*rxi*power(rxj, 12)

                          ) - exp(2*rxi)*power(rxi, 4)*

                          (110*power(rxi, 18)*power(rxj, 4)*

                           (663390 + 1216215*rxj + 1105650*power(rxj, 2) +

                            663390*power(rxj, 3) + 294840*power(rxj, 4) + 103194*power(rxj, 5) +

                            29484*power(rxj, 6) + 7020*power(rxj, 7) + 1404*power(rxj, 8) +

                            237*power(rxj, 9) + 30*power(rxj, 10) + 2*power(rxj, 11)) -

                           330*power(rxi, 16)*power(rxj, 6)*

                           (810810 + 1486485*rxj + 1351350*power(rxj, 2) +

                            810810*power(rxj, 3) + 360360*power(rxj, 4) + 126126*power(rxj, 5) +

                            36036*power(rxj, 6) + 8556*power(rxj, 7) + 1740*power(rxj, 8) +

                            291*power(rxj, 9) + 34*power(rxj, 10) + 2*power(rxj, 11)) +

                           330*power(rxi, 6)*power(rxj, 16)*

                           (3169530 + 7960680*rxj + 5798520*power(rxj, 2) +

                            3144960*power(rxj, 3) + 1572480*power(rxj, 4) +

                            638001*power(rxj, 5) + 191646*power(rxj, 6) + 41886*power(rxj, 7) +

                            6630*power(rxj, 8) + 741*power(rxj, 9) + 54*power(rxj, 10) +

                            2*power(rxj, 11)) - 110*power(rxi, 4)*power(rxj, 18)*

                           (12162150 + 8108100*rxj + 6486480*power(rxj, 2) +

                            5675670*power(rxj, 3) + 3243240*power(rxj, 4) +

                            1216215*power(rxj, 5) + 319410*power(rxj, 6) + 61074*power(rxj, 7) +

                            8586*power(rxj, 8) + 867*power(rxj, 9) + 58*power(rxj, 10) +

                            2*power(rxj, 11)) + power(rxi, 22)*

                           (935550 + 1715175*rxj + 1559250*power(rxj, 2) +

                            935550*power(rxj, 3) + 415800*power(rxj, 4) + 145530*power(rxj, 5) +

                            41580*power(rxj, 6) + 9900*power(rxj, 7) + 1980*power(rxj, 8) +

                            330*power(rxj, 9) + 44*power(rxj, 10) + 4*power(rxj, 11)) -

                           11*power(rxi, 20)*power(rxj, 2)*

                           (1105650 + 2027025*rxj + 1842750*power(rxj, 2) +

                            1105650*power(rxj, 3) + 491400*power(rxj, 4) +

                            171990*power(rxj, 5) + 49140*power(rxj, 6) + 11700*power(rxj, 7) +

                            2340*power(rxj, 8) + 390*power(rxj, 9) + 52*power(rxj, 10) +

                            4*power(rxj, 11)) + 11*power(rxi, 2)*power(rxj, 20)*

                           (-48648600 + 2027025*rxj + 44594550*power(rxj, 2) +

                            36486450*power(rxj, 3) + 16216200*power(rxj, 4) +

                            4864860*power(rxj, 5) + 1065960*power(rxj, 6) +

                            176040*power(rxj, 7) + 21960*power(rxj, 8) + 2010*power(rxj, 9) +

                            124*power(rxj, 10) + 4*power(rxj, 11)) -

                           power(rxj, 22)*(340540200 + 468242775*rxj + 312161850*power(rxj, 2) +

                                             133783650*power(rxj, 3) + 41164200*power(rxj, 4) +

                                             9604980*power(rxj, 5) + 1746360*power(rxj, 6) +

                                             249480*power(rxj, 7) + 27720*power(rxj, 8) + 2310*power(rxj, 9) +

                                             132*power(rxj, 10) + 4*power(rxj, 11)) +

                           165*power(rxi, 14)*power(rxj, 8)*

                           (4054050 + 7432425*rxj + 6756750*power(rxj, 2) +

                            4054050*power(rxj, 3) + 1801800*power(rxj, 4) +

                            631260*power(rxj, 5) + 178920*power(rxj, 6) + 43176*power(rxj, 7) +

                            8904*power(rxj, 8) + 1428*power(rxj, 9) + 152*power(rxj, 10) +

                            8*power(rxj, 11)) - 231*power(rxi, 12)*power(rxj, 10)*

                           (5212350 + 9555975*rxj + 8687250*power(rxj, 2) +

                            5209650*power(rxj, 3) + 2327400*power(rxj, 4) +

                            801540*power(rxj, 5) + 230040*power(rxj, 6) + 57240*power(rxj, 7) +

                            11640*power(rxj, 8) + 1740*power(rxj, 9) + 168*power(rxj, 10) +

                            8*power(rxj, 11)) + 231*power(rxi, 10)*power(rxj, 12)*

                           (6949800 + 12746025*rxj + 11535750*power(rxj, 2) +

                            7056450*power(rxj, 3) + 3040200*power(rxj, 4) +

                            1051920*power(rxj, 5) + 316800*power(rxj, 6) + 79680*power(rxj, 7) +

                            15360*power(rxj, 8) + 2100*power(rxj, 9) + 184*power(rxj, 10) +

                            8*power(rxj, 11)) - 165*power(rxi, 8)*power(rxj, 14)*

                           (9775080 + 17424855*rxj + 17019450*power(rxj, 2) +

                            9519930*power(rxj, 3) + 4059720*power(rxj, 4) +

                            1519056*power(rxj, 5) + 475776*power(rxj, 6) + 114720*power(rxj, 7) +

                            20256*power(rxj, 8) + 2508*power(rxj, 9) + 200*power(rxj, 10) +

                            8*power(rxj, 11))))/

                         (935550*exp(2*(rxi + rxj))*power(rxi - rxj, 13)*power(rxi + rxj, 13))

                         );
        }

    }
    return S;
}


double Slater_6S_1S(double r, double xi, double xj)
{
    return Slater_1S_6S(r, xj, xi);
}

#endif
