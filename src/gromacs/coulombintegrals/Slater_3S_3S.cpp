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

#if HAVE_LIBCLN
cl_R Slater_3S_3S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (793LL*xi)/3072LL

            ;
        }
        else
        {
            S = (1LL/r)*((-1437004800LL + 1437004800LL*exp(2LL*rxi) - 2503064025LL*rxi -

                          2132118450LL*Power(rxi, 2LL) - 1180664100LL*Power(rxi, 3LL) -

                          476506800LL*Power(rxi, 4LL) - 148856400LL*Power(rxi, 5LL) -

                          37255680LL*Power(rxi, 6LL) - 7603200LL*Power(rxi, 7LL) - 1267200LL*Power(rxi, 8LL) -

                          168960LL*Power(rxi, 9LL) - 16896LL*Power(rxi, 10LL) - 1024LL*Power(rxi, 11LL))/

                         (1.4370048e9*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 10LL) + 11LL*Power(xi, 9LL)*xj + 55LL*Power(xi, 8LL)*Power(xj, 2LL) +

                        165LL*Power(xi, 7LL)*Power(xj, 3LL) + 330LL*Power(xi, 6LL)*Power(xj, 4LL) +

                        462LL*Power(xi, 5LL)*Power(xj, 5LL) + 330LL*Power(xi, 4LL)*Power(xj, 6LL) +

                        165LL*Power(xi, 3LL)*Power(xj, 7LL) + 55LL*Power(xi, 2LL)*Power(xj, 8LL) +

                        11LL*xi*Power(xj, 9LL) + Power(xj, 10LL)))/(3LL*Power(xi + xj, 11LL))

            ;
        }
        else
        {
            S = (1LL/r)*((135LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 11LL) +

                          exp(2LL*rxj)*Power(rxj, 8LL)*

                          (-150LL*Power(rxi, 18LL) - 6LL*Power(rxi, 19LL) + 135LL*Power(rxj, 14LL) +

                           225LL*rxi*Power(rxj, 14LL) + 10LL*Power(rxi, 17LL)*(-165LL + Power(rxj, 2LL)) -

                           30LL*Power(rxi, 16LL)*(330LL + Power(rxj, 2LL)) +

                           45LL*Power(rxi, 3LL)*Power(rxj, 12LL)*(-55LL + 2LL*Power(rxj, 2LL)) +

                           45LL*Power(rxi, 2LL)*Power(rxj, 12LL)*(-33LL + 4LL*Power(rxj, 2LL)) +

                           Power(rxi, 9LL)*Power(rxj, 6LL)*

                           (234135LL - 4950LL*Power(rxj, 2LL) - 34LL*Power(rxj, 4LL)) -

                           5LL*Power(rxi, 7LL)*Power(rxj, 8LL)*

                           (6237LL - 1242LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           3LL*Power(rxi, 5LL)*Power(rxj, 10LL)*

                           (4125LL - 330LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 4LL)*Power(rxj, 10LL)*

                           (495LL - 132LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) -

                           165LL*Power(rxi, 6LL)*Power(rxj, 8LL)*

                           (135LL - 60LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) -

                           5LL*Power(rxi, 13LL)*Power(rxj, 2LL)*

                           (43875LL - 3438LL*Power(rxj, 2LL) + 22LL*Power(rxj, 4LL)) +

                           5LL*Power(rxi, 11LL)*Power(rxj, 4LL)*

                           (7695LL - 2442LL*Power(rxj, 2LL) + 22LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 8LL)*Power(rxj, 6LL)*

                           (-33LL - 3564LL*Power(rxj, 2LL) + 26LL*Power(rxj, 4LL)) +

                           Power(rxi, 15LL)*(-32175LL - 3690LL*Power(rxj, 2LL) + 34LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 10LL)*Power(rxj, 4LL)*

                           (-32277LL + 1364LL*Power(rxj, 2LL) + 66LL*Power(rxj, 4LL)) +

                           15LL*Power(rxi, 14LL)*(-3003LL - 2932LL*Power(rxj, 2LL) + 94LL*Power(rxj, 4LL)) -

                           15LL*Power(rxi, 12LL)*Power(rxj, 2LL)*

                           (28119LL - 5252LL*Power(rxj, 2LL) + 154LL*Power(rxj, 4LL))) +

                          exp(2LL*rxi)*Power(rxi, 8LL)*

                          (-5LL*Power(rxi, 2LL)*Power(rxj, 12LL)*

                           (-84357LL - 43875LL*rxj - 8796LL*Power(rxj, 2LL) - 738LL*Power(rxj, 3LL) -

                            6LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) -

                           3LL*Power(rxi, 14LL)*(45LL + 75LL*rxj + 60LL*Power(rxj, 2LL) + 30LL*Power(rxj, 3LL) +

                                                 10LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) -

                           55LL*Power(rxi, 8LL)*Power(rxj, 6LL)*

                           (-405LL - 567LL*rxj - 972LL*Power(rxj, 2LL) - 90LL*Power(rxj, 3LL) +

                            18LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           55LL*Power(rxi, 6LL)*Power(rxj, 8LL)*

                           (9LL - 4257LL*rxj - 372LL*Power(rxj, 2LL) + 222LL*Power(rxj, 3LL) +

                            42LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           3LL*Power(rxj, 14LL)*(15015LL + 10725LL*rxj + 3300LL*Power(rxj, 2LL) +

                                                 550LL*Power(rxj, 3LL) + 50LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           5LL*Power(rxi, 12LL)*Power(rxj, 2LL)*

                           (297LL + 495LL*rxj + 396LL*Power(rxj, 2LL) + 198LL*Power(rxj, 3LL) +

                            66LL*Power(rxj, 4LL) + 2LL*Power(rxj, 5LL)) +

                           Power(rxi, 10LL)*Power(rxj, 4LL)*

                           (-7425LL - 12375LL*rxj - 9900LL*Power(rxj, 2LL) - 6210LL*Power(rxj, 3LL) -

                            390LL*Power(rxj, 4LL) + 34LL*Power(rxj, 5LL)) -

                           Power(rxi, 4LL)*Power(rxj, 10LL)*

                           (-484155LL + 38475LL*rxj + 78780LL*Power(rxj, 2LL) + 17190LL*Power(rxj, 3LL) +

                            1410LL*Power(rxj, 4LL) + 34LL*Power(rxj, 5LL))))/

                         (135LL*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 11LL)*Power(rxi + rxj, 11LL))

                         );
        }

    }
    return S;
}

#else

double Slater_3S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (793*xi)/3072

            ;
        }
        else
        {
            S = (1/r)*((-1437004800 + 1437004800*exp(2*rxi) - 2503064025*rxi -

                          2132118450*pow(rxi, 2) - 1180664100*pow(rxi, 3) -

                          476506800*pow(rxi, 4) - 148856400*pow(rxi, 5) -

                          37255680*pow(rxi, 6) - 7603200*pow(rxi, 7) - 1267200*pow(rxi, 8) -

                          168960*pow(rxi, 9) - 16896*pow(rxi, 10) - 1024*pow(rxi, 11))/

                         (1.4370048e9*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 10) + 11*pow(xi, 9)*xj + 55*pow(xi, 8)*pow(xj, 2) +

                        165*pow(xi, 7)*pow(xj, 3) + 330*pow(xi, 6)*pow(xj, 4) +

                        462*pow(xi, 5)*pow(xj, 5) + 330*pow(xi, 4)*pow(xj, 6) +

                        165*pow(xi, 3)*pow(xj, 7) + 55*pow(xi, 2)*pow(xj, 8) +

                        11*xi*pow(xj, 9) + pow(xj, 10)))/(3*pow(xi + xj, 11))

            ;
        }
        else
        {
            S = (1/r)*((135*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 11) +

                          exp(2*rxj)*pow(rxj, 8)*

                          (-150*pow(rxi, 18) - 6*pow(rxi, 19) + 135*pow(rxj, 14) +

                           225*rxi*pow(rxj, 14) + 10*pow(rxi, 17)*(-165 + pow(rxj, 2)) -

                           30*pow(rxi, 16)*(330 + pow(rxj, 2)) +

                           45*pow(rxi, 3)*pow(rxj, 12)*(-55 + 2*pow(rxj, 2)) +

                           45*pow(rxi, 2)*pow(rxj, 12)*(-33 + 4*pow(rxj, 2)) +

                           pow(rxi, 9)*pow(rxj, 6)*

                           (234135 - 4950*pow(rxj, 2) - 34*pow(rxj, 4)) -

                           5*pow(rxi, 7)*pow(rxj, 8)*

                           (6237 - 1242*pow(rxj, 2) + 2*pow(rxj, 4)) +

                           3*pow(rxi, 5)*pow(rxj, 10)*

                           (4125 - 330*pow(rxj, 2) + 2*pow(rxj, 4)) +

                           15*pow(rxi, 4)*pow(rxj, 10)*

                           (495 - 132*pow(rxj, 2) + 2*pow(rxj, 4)) -

                           165*pow(rxi, 6)*pow(rxj, 8)*

                           (135 - 60*pow(rxj, 2) + 2*pow(rxj, 4)) -

                           5*pow(rxi, 13)*pow(rxj, 2)*

                           (43875 - 3438*pow(rxj, 2) + 22*pow(rxj, 4)) +

                           5*pow(rxi, 11)*pow(rxj, 4)*

                           (7695 - 2442*pow(rxj, 2) + 22*pow(rxj, 4)) +

                           15*pow(rxi, 8)*pow(rxj, 6)*

                           (-33 - 3564*pow(rxj, 2) + 26*pow(rxj, 4)) +

                           pow(rxi, 15)*(-32175 - 3690*pow(rxj, 2) + 34*pow(rxj, 4)) +

                           15*pow(rxi, 10)*pow(rxj, 4)*

                           (-32277 + 1364*pow(rxj, 2) + 66*pow(rxj, 4)) +

                           15*pow(rxi, 14)*(-3003 - 2932*pow(rxj, 2) + 94*pow(rxj, 4)) -

                           15*pow(rxi, 12)*pow(rxj, 2)*

                           (28119 - 5252*pow(rxj, 2) + 154*pow(rxj, 4))) +

                          exp(2*rxi)*pow(rxi, 8)*

                          (-5*pow(rxi, 2)*pow(rxj, 12)*

                           (-84357 - 43875*rxj - 8796*pow(rxj, 2) - 738*pow(rxj, 3) -

                            6*pow(rxj, 4) + 2*pow(rxj, 5)) -

                           3*pow(rxi, 14)*(45 + 75*rxj + 60*pow(rxj, 2) + 30*pow(rxj, 3) +

                                                 10*pow(rxj, 4) + 2*pow(rxj, 5)) -

                           55*pow(rxi, 8)*pow(rxj, 6)*

                           (-405 - 567*rxj - 972*pow(rxj, 2) - 90*pow(rxj, 3) +

                            18*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           55*pow(rxi, 6)*pow(rxj, 8)*

                           (9 - 4257*rxj - 372*pow(rxj, 2) + 222*pow(rxj, 3) +

                            42*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           3*pow(rxj, 14)*(15015 + 10725*rxj + 3300*pow(rxj, 2) +

                                                 550*pow(rxj, 3) + 50*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           5*pow(rxi, 12)*pow(rxj, 2)*

                           (297 + 495*rxj + 396*pow(rxj, 2) + 198*pow(rxj, 3) +

                            66*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           pow(rxi, 10)*pow(rxj, 4)*

                           (-7425 - 12375*rxj - 9900*pow(rxj, 2) - 6210*pow(rxj, 3) -

                            390*pow(rxj, 4) + 34*pow(rxj, 5)) -

                           pow(rxi, 4)*pow(rxj, 10)*

                           (-484155 + 38475*rxj + 78780*pow(rxj, 2) + 17190*pow(rxj, 3) +

                            1410*pow(rxj, 4) + 34*pow(rxj, 5))))/

                         (135*exp(2*(rxi + rxj))*pow(rxi - rxj, 11)*pow(rxi + rxj, 11))

                         );
        }

    }
    return S;
}

#endif
