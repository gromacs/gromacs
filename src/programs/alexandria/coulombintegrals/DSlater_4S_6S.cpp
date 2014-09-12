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

cl_R DSlater_4S_6S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-(1e9*cl_float(537882537, precision)+342262612LL)*xi + (1e9*cl_float(583896481, precision)+962393600LL)*exp(2LL*r*xi)*xi -

                  (1e9*cl_float(983737185, precision)+444263250LL)*r*Power(xi, 2LL) -

                  (1e9*cl_float(892447675, precision)+429810500LL)*Power(r, 2LL)*Power(xi, 3LL) -

                  (1e9*cl_float(535089949, precision)+244650800LL)*Power(r, 3LL)*Power(xi, 4LL) -

                  (1e9*cl_float(238344161, precision)+324519250LL)*Power(r, 4LL)*Power(xi, 5LL) -

                  840485675935105200LL*Power(r, 5LL)*Power(xi, 6LL) -

                  244148817764451600LL*Power(r, 6LL)*Power(xi, 7LL) -

                  60013996461619200LL*Power(r, 7LL)*Power(xi, 8LL) -

                  12723407730633600LL*Power(r, 8LL)*Power(xi, 9LL) -

                  2358784581753600LL*Power(r, 9LL)*Power(xi, 10LL) -

                  386141399804160LL*Power(r, 10LL)*Power(xi, 11LL) -

                  56170897735680LL*Power(r, 11LL)*Power(xi, 12LL) -

                  7281412669440LL*Power(r, 12LL)*Power(xi, 13LL) -

                  840163000320LL*Power(r, 13LL)*Power(xi, 14LL) -

                  85730918400LL*Power(r, 14LL)*Power(xi, 15LL) -

                  7620526080LL*Power(r, 15LL)*Power(xi, 16LL) -

                  571539456LL*Power(r, 16LL)*Power(xi, 17LL) -

                  33619968LL*Power(r, 17LL)*Power(xi, 18LL) - 1245184LL*Power(r, 18LL)*Power(xi, 19LL))/

                (2.919482409811968e18*exp(2LL*r*xi)*r) +

                (-(1e9*cl_float(291948240, precision)+981196800LL) + (1e9*cl_float(291948240, precision)+981196800LL)*exp(2LL*r*xi) -

                 (1e9*cl_float(537882537, precision)+342262612LL)*r*xi -

                 (1e9*cl_float(491868592, precision)+722131625LL)*Power(r, 2LL)*Power(xi, 2LL) -

                 (1e9*cl_float(297482558, precision)+476603500LL)*Power(r, 3LL)*Power(xi, 3LL) -

                 (1e9*cl_float(133772487, precision)+311162700LL)*Power(r, 4LL)*Power(xi, 4LL) -

                 476688322649038500LL*Power(r, 5LL)*Power(xi, 5LL) -

                 140080945989184200LL*Power(r, 6LL)*Power(xi, 6LL) -

                 34878402537778800LL*Power(r, 7LL)*Power(xi, 7LL) -

                 7501749557702400LL*Power(r, 8LL)*Power(xi, 8LL) -

                 1413711970070400LL*Power(r, 9LL)*Power(xi, 9LL) -

                 235878458175360LL*Power(r, 10LL)*Power(xi, 10LL) -

                 35103763618560LL*Power(r, 11LL)*Power(xi, 11LL) -

                 4680908144640LL*Power(r, 12LL)*Power(xi, 12LL) -

                 560108666880LL*Power(r, 13LL)*Power(xi, 13LL) -

                 60011642880LL*Power(r, 14LL)*Power(xi, 14LL) -

                 5715394560LL*Power(r, 15LL)*Power(xi, 15LL) -

                 476282880LL*Power(r, 16LL)*Power(xi, 16LL) -

                 33619968LL*Power(r, 17LL)*Power(xi, 17LL) - 1867776LL*Power(r, 18LL)*Power(xi, 18LL) -

                 65536LL*Power(r, 19LL)*Power(xi, 19LL))/

                (2.919482409811968e18*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-(1e9*cl_float(291948240, precision)+981196800LL) + (1e9*cl_float(291948240, precision)+981196800LL)*exp(2LL*r*xi) -

                     (1e9*cl_float(537882537, precision)+342262612LL)*r*xi -

                     (1e9*cl_float(491868592, precision)+722131625LL)*Power(r, 2LL)*Power(xi, 2LL) -

                     (1e9*cl_float(297482558, precision)+476603500LL)*Power(r, 3LL)*Power(xi, 3LL) -

                     (1e9*cl_float(133772487, precision)+311162700LL)*Power(r, 4LL)*Power(xi, 4LL) -

                     476688322649038500LL*Power(r, 5LL)*Power(xi, 5LL) -

                     140080945989184200LL*Power(r, 6LL)*Power(xi, 6LL) -

                     34878402537778800LL*Power(r, 7LL)*Power(xi, 7LL) -

                     7501749557702400LL*Power(r, 8LL)*Power(xi, 8LL) -

                     1413711970070400LL*Power(r, 9LL)*Power(xi, 9LL) -

                     235878458175360LL*Power(r, 10LL)*Power(xi, 10LL) -

                     35103763618560LL*Power(r, 11LL)*Power(xi, 11LL) -

                     4680908144640LL*Power(r, 12LL)*Power(xi, 12LL) -

                     560108666880LL*Power(r, 13LL)*Power(xi, 13LL) -

                     60011642880LL*Power(r, 14LL)*Power(xi, 14LL) -

                     5715394560LL*Power(r, 15LL)*Power(xi, 15LL) -

                     476282880LL*Power(r, 16LL)*Power(xi, 16LL) -

                     33619968LL*Power(r, 17LL)*Power(xi, 17LL) -

                     1867776LL*Power(r, 18LL)*Power(xi, 18LL) - 65536LL*Power(r, 19LL)*Power(xi, 19LL)))/

                (1.459741204905984e18*exp(2LL*r*xi)*r)

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
            S = (1871100LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 19LL) +

                 495LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-672LL*Power(r, 6LL)*Power(xi, 30LL) - 12LL*Power(r, 7LL)*Power(xi, 31LL) +

                  3780LL*Power(xj, 24LL) + 6615LL*r*xi*Power(xj, 24LL) -

                  136LL*Power(r, 5LL)*Power(xi, 29LL)*(126LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  1890LL*Power(xi, 2LL)*Power(xj, 22LL)*(-38LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  315LL*r*Power(xi, 3LL)*Power(xj, 22LL)*

                  (-399LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  84LL*Power(r, 4LL)*Power(xi, 28LL)*(3060LL + 121LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  630LL*Power(xi, 4LL)*Power(xj, 20LL)*

                  (1026LL - 171LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  63LL*r*Power(xi, 5LL)*Power(xj, 20LL)*

                  (17955LL - 950LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  84LL*Power(r, 2LL)*Power(xi, 26LL)*

                  (-174420LL - 71535LL*Power(r, 2LL)*Power(xj, 2LL) +

            179LL*Power(r, 4LL)*Power(xj, 4LL)) -

                  63LL*r*Power(xi, 19LL)*Power(xj, 6LL)*

                  (468377895LL - 14898090LL*Power(r, 2LL)*Power(xj, 2LL) +

            78812LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  Power(xi, 27LL)*(-2441880LL*Power(r, 3LL) -

                                   327978LL*Power(r, 5LL)*Power(xj, 2LL) + 496LL*Power(r, 7LL)*Power(xj, 4LL)) +

                  2LL*r*Power(xi, 11LL)*Power(xj, 14LL)*

                  (613624095LL + 56366730LL*Power(r, 2LL)*Power(xj, 2LL) +

            383607LL*Power(r, 4LL)*Power(xj, 4LL) - 248LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  42LL*Power(xi, 6LL)*Power(xj, 18LL)*

                  (-87210LL + 23085LL*Power(r, 2LL)*Power(xj, 2LL) -

            570LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  798LL*Power(xi, 8LL)*Power(xj, 16LL)*

                  (-18360LL + 6885LL*Power(r, 2LL)*Power(xj, 2LL) -

            270LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  3LL*r*Power(xi, 7LL)*Power(xj, 18LL)*

                  (-2136645LL + 179550LL*Power(r, 2LL)*Power(xj, 2LL) -

            2394LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  1596LL*Power(xi, 14LL)*Power(xj, 10LL)*

                  (-34484670LL - 2408985LL*Power(r, 2LL)*Power(xj, 2LL) +

            32810LL*Power(r, 4LL)*Power(xj, 4LL) + 22LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  7980LL*Power(xi, 16LL)*Power(xj, 8LL)*

                  (15696909LL - 494343LL*Power(r, 2LL)*Power(xj, 2LL) -

            4182LL*Power(r, 4LL)*Power(xj, 4LL) + 34LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  2LL*r*Power(xi, 9LL)*Power(xj, 16LL)*

                  (19433295LL - 690795LL*Power(r, 2LL)*Power(xj, 2LL) +

            55251LL*Power(r, 4LL)*Power(xj, 4LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  6LL*r*Power(xi, 25LL)*(-8546580LL - 11329605LL*Power(r, 2LL)*Power(xj, 2LL) -

                                         24003LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  6LL*r*Power(xi, 13LL)*Power(xj, 12LL)*

                  (-2361196215LL - 54738810LL*Power(r, 2LL)*Power(xj, 2LL) +

            388626LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  38LL*r*Power(xi, 15LL)*Power(xj, 10LL)*

                  (808181955LL - 17168130LL*Power(r, 2LL)*Power(xj, 2LL) -

            32130LL*Power(r, 4LL)*Power(xj, 4LL) + 106LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  84LL*Power(xi, 10LL)*Power(xj, 14LL)*

                  (3168630LL + 683145LL*Power(r, 2LL)*Power(xj, 2LL) +

            54315LL*Power(r, 4LL)*Power(xj, 4LL) + 193LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  19LL*r*Power(xi, 17LL)*Power(xj, 8LL)*

                  (-2525985LL + 33479460LL*Power(r, 2LL)*Power(xj, 2LL) -

            406980LL*Power(r, 4LL)*Power(xj, 4LL) + 272LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  84LL*Power(xi, 12LL)*Power(xj, 12LL)*

                  (-88925130LL - 19869345LL*Power(r, 2LL)*Power(xj, 2LL) -

            235790LL*Power(r, 4LL)*Power(xj, 4LL) + 643LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  210LL*Power(xi, 18LL)*Power(xj, 6LL)*

                  (-496605582LL + 32638599LL*Power(r, 2LL)*Power(xj, 2LL) -

            564604LL*Power(r, 4LL)*Power(xj, 4LL) + 1292LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  42LL*Power(xi, 20LL)*Power(xj, 4LL)*

                  (-777723210LL - 46394505LL*Power(r, 2LL)*Power(xj, 2LL) +

            625670LL*Power(r, 4LL)*Power(xj, 4LL) + 1292LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  42LL*Power(xi, 24LL)*(-1918620LL - 11344995LL*Power(r, 2LL)*Power(xj, 2LL) -

                                        323070LL*Power(r, 4LL)*Power(xj, 4LL) + 2114LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  r*Power(xi, 23LL)*Power(xj, 2LL)*

                  (1919335635LL + 275096430LL*Power(r, 2LL)*Power(xj, 2LL) -

            3302586LL*Power(r, 4LL)*Power(xj, 4LL) + 4028LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  r*Power(xi, 21LL)*Power(xj, 4LL)*

                  (-14708379735LL + 255168270LL*Power(r, 2LL)*Power(xj, 2LL) -

            2899134LL*Power(r, 4LL)*Power(xj, 4LL) + 5168LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  42LL*Power(xi, 22LL)*Power(xj, 2LL)*

                  (81654210LL + 66273255LL*Power(r, 2LL)*Power(xj, 2LL) -

            1203870LL*Power(r, 4LL)*Power(xj, 4LL) + 5206LL*Power(r, 6LL)*Power(xj, 6LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 10LL)*

                 (21318LL*Power(xi, 14LL)*Power(xj, 14LL)*

                  (-3146850LL + 4890375LL*r*xj - 24522750LL*Power(r, 2LL)*Power(xj, 2LL) +

            12162150LL*Power(r, 3LL)*Power(xj, 3LL) -

            1549800LL*Power(r, 4LL)*Power(xj, 4LL) -

            1615950LL*Power(r, 5LL)*Power(xj, 5LL) -

            185220LL*Power(r, 6LL)*Power(xj, 6LL) + 12240LL*Power(r, 7LL)*Power(xj, 7LL) +

            3960LL*Power(r, 8LL)*Power(xj, 8LL) + 300LL*Power(r, 9LL)*Power(xj, 9LL) +

            8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  3LL*Power(xi, 24LL)*Power(xj, 4LL)*

                  (53326350LL + 97764975LL*r*xj + 88877250LL*Power(r, 2LL)*Power(xj, 2LL) +

            53326350LL*Power(r, 3LL)*Power(xj, 3LL) +

            23700600LL*Power(r, 4LL)*Power(xj, 4LL) +

            8295210LL*Power(r, 5LL)*Power(xj, 5LL) +

            2370060LL*Power(r, 6LL)*Power(xj, 6LL) +

            564300LL*Power(r, 7LL)*Power(xj, 7LL) + 112860LL*Power(r, 8LL)*Power(xj, 8LL) +

            22440LL*Power(r, 9LL)*Power(xj, 9LL) + 1056LL*Power(r, 10LL)*Power(xj, 10LL) -

            20LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  4LL*Power(xj, 28LL)*(13749310575LL + 13749310575LL*r*xj +

                                       6547290750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       1964187225LL*Power(r, 3LL)*Power(xj, 3LL) +

                                       413513100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       64324260LL*Power(r, 5LL)*Power(xj, 5LL) +

                                       7567560LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       675675LL*Power(r, 7LL)*Power(xj, 7LL) + 45045LL*Power(r, 8LL)*Power(xj, 8LL) +

                                       2145LL*Power(r, 9LL)*Power(xj, 9LL) + 66LL*Power(r, 10LL)*Power(xj, 10LL) +

                                       Power(r, 11LL)*Power(xj, 11LL)) -

                  1254LL*Power(xi, 16LL)*Power(xj, 12LL)*

                  (-20241900LL - 38315025LL*r*xj - 21687750LL*Power(r, 2LL)*Power(xj, 2LL) -

            50122800LL*Power(r, 3LL)*Power(xj, 3LL) +

            14137200LL*Power(r, 4LL)*Power(xj, 4LL) -

            5853330LL*Power(r, 5LL)*Power(xj, 5LL) -

            2687580LL*Power(r, 6LL)*Power(xj, 6LL) -

            208530LL*Power(r, 7LL)*Power(xj, 7LL) + 19530LL*Power(r, 8LL)*Power(xj, 8LL) +

            3630LL*Power(r, 9LL)*Power(xj, 9LL) + 172LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  627LL*Power(xi, 12LL)*Power(xj, 16LL)*

                  (-1240964550LL + 4740389325LL*r*xj -

            3311818650LL*Power(r, 2LL)*Power(xj, 2LL) +

            134804250LL*Power(r, 3LL)*Power(xj, 3LL) +

            407673000LL*Power(r, 4LL)*Power(xj, 4LL) +

            58641030LL*Power(r, 5LL)*Power(xj, 5LL) -

            3549420LL*Power(r, 6LL)*Power(xj, 6LL) -

            1641060LL*Power(r, 7LL)*Power(xj, 7LL) -

            167940LL*Power(r, 8LL)*Power(xj, 8LL) - 6990LL*Power(r, 9LL)*Power(xj, 9LL) -

            36LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  Power(xi, 28LL)*(935550LL + 1715175LL*r*xj +

                                   1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                   330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                   4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  2LL*Power(xi, 2LL)*Power(xj, 26LL)*

                  (-937068397650LL - 815439881025LL*r*xj -

            332904552750LL*Power(r, 2LL)*Power(xj, 2LL) -

            84006776700LL*Power(r, 3LL)*Power(xj, 3LL) -

            14504767200LL*Power(r, 4LL)*Power(xj, 4LL) -

            1786235220LL*Power(r, 5LL)*Power(xj, 5LL) -

            157754520LL*Power(r, 6LL)*Power(xj, 6LL) -

            9667350LL*Power(r, 7LL)*Power(xj, 7LL) -

            367290LL*Power(r, 8LL)*Power(xj, 8LL) - 5115LL*Power(r, 9LL)*Power(xj, 9LL) +

            198LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  6LL*Power(xi, 4LL)*Power(xj, 24LL)*

                  (-2262441500550LL - 1503711230175LL*r*xj -

            426178264050LL*Power(r, 2LL)*Power(xj, 2LL) -

            60134347350LL*Power(r, 3LL)*Power(xj, 3LL) -

            2014551000LL*Power(r, 4LL)*Power(xj, 4LL) +

            846111420LL*Power(r, 5LL)*Power(xj, 5LL) +

            184864680LL*Power(r, 6LL)*Power(xj, 6LL) +

            20183130LL*Power(r, 7LL)*Power(xj, 7LL) +

            1367190LL*Power(r, 8LL)*Power(xj, 8LL) + 57255LL*Power(r, 9LL)*Power(xj, 9LL) +

            1298LL*Power(r, 10LL)*Power(xj, 10LL) + 10LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  Power(xi, 26LL)*Power(xj, 2LL)*

                  (17775450LL + 32588325LL*r*xj + 29625750LL*Power(r, 2LL)*Power(xj, 2LL) +

            17775450LL*Power(r, 3LL)*Power(xj, 3LL) +

            7900200LL*Power(r, 4LL)*Power(xj, 4LL) +

            2765070LL*Power(r, 5LL)*Power(xj, 5LL) +

            790020LL*Power(r, 6LL)*Power(xj, 6LL) + 188100LL*Power(r, 7LL)*Power(xj, 7LL) +

            37620LL*Power(r, 8LL)*Power(xj, 8LL) + 6270LL*Power(r, 9LL)*Power(xj, 9LL) +

            836LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  9LL*Power(xi, 22LL)*Power(xj, 6LL)*

                  (-100727550LL - 184667175LL*r*xj - 167879250LL*Power(r, 2LL)*Power(xj, 2LL) -

            100727550LL*Power(r, 3LL)*Power(xj, 3LL) -

            44767800LL*Power(r, 4LL)*Power(xj, 4LL) -

            15668730LL*Power(r, 5LL)*Power(xj, 5LL) -

            4476780LL*Power(r, 6LL)*Power(xj, 6LL) -

            971520LL*Power(r, 7LL)*Power(xj, 7LL) - 307560LL*Power(r, 8LL)*Power(xj, 8LL) -

            27060LL*Power(r, 9LL)*Power(xj, 9LL) + 264LL*Power(r, 10LL)*Power(xj, 10LL) +

            64LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  9LL*Power(xi, 6LL)*Power(xj, 22LL)*

                  (3452543428950LL + 1097992509075LL*r*xj -

            101420792550LL*Power(r, 2LL)*Power(xj, 2LL) -

            110557373850LL*Power(r, 3LL)*Power(xj, 3LL) -

            24909330600LL*Power(r, 4LL)*Power(xj, 4LL) -

            2686726350LL*Power(r, 5LL)*Power(xj, 5LL) -

            93485700LL*Power(r, 6LL)*Power(xj, 6LL) +

            12941280LL*Power(r, 7LL)*Power(xj, 7LL) +

            2081640LL*Power(r, 8LL)*Power(xj, 8LL) +

            137940LL*Power(r, 9LL)*Power(xj, 9LL) + 4664LL*Power(r, 10LL)*Power(xj, 10LL) +

            64LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  22LL*Power(xi, 20LL)*Power(xj, 8LL)*

                  (-164826900LL - 302182650LL*r*xj - 274711500LL*Power(r, 2LL)*Power(xj, 2LL) -

            164826900LL*Power(r, 3LL)*Power(xj, 3LL) -

            73256400LL*Power(r, 4LL)*Power(xj, 4LL) -

            26991090LL*Power(r, 5LL)*Power(xj, 5LL) -

            4622940LL*Power(r, 6LL)*Power(xj, 6LL) -

            2941110LL*Power(r, 7LL)*Power(xj, 7LL) -

            438930LL*Power(r, 8LL)*Power(xj, 8LL) - 5505LL*Power(r, 9LL)*Power(xj, 9LL) +

            2082LL*Power(r, 10LL)*Power(xj, 10LL) + 82LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  22LL*Power(xi, 18LL)*Power(xj, 10LL)*

                  (-494480700LL - 906547950LL*r*xj - 824134500LL*Power(r, 2LL)*Power(xj, 2LL) -

            475684650LL*Power(r, 3LL)*Power(xj, 3LL) -

            294953400LL*Power(r, 4LL)*Power(xj, 4LL) +

            2663010LL*Power(r, 5LL)*Power(xj, 5LL) -

            40797540LL*Power(r, 6LL)*Power(xj, 6LL) -

            10248390LL*Power(r, 7LL)*Power(xj, 7LL) -

            434610LL*Power(r, 8LL)*Power(xj, 8LL) + 65865LL*Power(r, 9LL)*Power(xj, 9LL) +

            6366LL*Power(r, 10LL)*Power(xj, 10LL) + 136LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  11LL*Power(xi, 8LL)*Power(xj, 20LL)*

                  (-2338604626050LL + 656001834075LL*r*xj +

            504510561450LL*Power(r, 2LL)*Power(xj, 2LL) +

            51560967150LL*Power(r, 3LL)*Power(xj, 3LL) -

            15574998600LL*Power(r, 4LL)*Power(xj, 4LL) -

            5055778350LL*Power(r, 5LL)*Power(xj, 5LL) -

            626213700LL*Power(r, 6LL)*Power(xj, 6LL) -

            34768620LL*Power(r, 7LL)*Power(xj, 7LL) +

            207540LL*Power(r, 8LL)*Power(xj, 8LL) + 150240LL*Power(r, 9LL)*Power(xj, 9LL) +

            8464LL*Power(r, 10LL)*Power(xj, 10LL) + 164LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  11LL*Power(xi, 10LL)*Power(xj, 18LL)*

                  (742805182350LL - 933111659025LL*r*xj +

            57080542050LL*Power(r, 2LL)*Power(xj, 2LL) +

            129505209750LL*Power(r, 3LL)*Power(xj, 3LL) +

            19066887000LL*Power(r, 4LL)*Power(xj, 4LL) -

            1817573310LL*Power(r, 5LL)*Power(xj, 5LL) -

            810647460LL*Power(r, 6LL)*Power(xj, 6LL) -

            97669980LL*Power(r, 7LL)*Power(xj, 7LL) -

            5173020LL*Power(r, 8LL)*Power(xj, 8LL) - 37770LL*Power(r, 9LL)*Power(xj, 9LL) +

            8212LL*Power(r, 10LL)*Power(xj, 10LL) + 272LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (1.8711e6*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 19LL)*

                 Power(xi + xj, 19LL)) + (1871100LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 19LL) +

                                          495LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                                          (-672LL*Power(r, 6LL)*Power(xi, 30LL) - 12LL*Power(r, 7LL)*Power(xi, 31LL) +

                                  3780LL*Power(xj, 24LL) + 6615LL*r*xi*Power(xj, 24LL) -

                                  136LL*Power(r, 5LL)*Power(xi, 29LL)*(126LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  1890LL*Power(xi, 2LL)*Power(xj, 22LL)*(-38LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  315LL*r*Power(xi, 3LL)*Power(xj, 22LL)*

                                  (-399LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  84LL*Power(r, 4LL)*Power(xi, 28LL)*(3060LL + 121LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  630LL*Power(xi, 4LL)*Power(xj, 20LL)*

                                  (1026LL - 171LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  63LL*r*Power(xi, 5LL)*Power(xj, 20LL)*

                                  (17955LL - 950LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  84LL*Power(r, 2LL)*Power(xi, 26LL)*

                                  (-174420LL - 71535LL*Power(r, 2LL)*Power(xj, 2LL) +

            179LL*Power(r, 4LL)*Power(xj, 4LL)) -

                                  63LL*r*Power(xi, 19LL)*Power(xj, 6LL)*

                                  (468377895LL - 14898090LL*Power(r, 2LL)*Power(xj, 2LL) +

            78812LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  Power(xi, 27LL)*(-2441880LL*Power(r, 3LL) -

                                                   327978LL*Power(r, 5LL)*Power(xj, 2LL) + 496LL*Power(r, 7LL)*Power(xj, 4LL)) +

                                  2LL*r*Power(xi, 11LL)*Power(xj, 14LL)*

                                  (613624095LL + 56366730LL*Power(r, 2LL)*Power(xj, 2LL) +

            383607LL*Power(r, 4LL)*Power(xj, 4LL) - 248LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  42LL*Power(xi, 6LL)*Power(xj, 18LL)*

                                  (-87210LL + 23085LL*Power(r, 2LL)*Power(xj, 2LL) -

            570LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  798LL*Power(xi, 8LL)*Power(xj, 16LL)*

                                  (-18360LL + 6885LL*Power(r, 2LL)*Power(xj, 2LL) -

            270LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  3LL*r*Power(xi, 7LL)*Power(xj, 18LL)*

                                  (-2136645LL + 179550LL*Power(r, 2LL)*Power(xj, 2LL) -

            2394LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  1596LL*Power(xi, 14LL)*Power(xj, 10LL)*

                                  (-34484670LL - 2408985LL*Power(r, 2LL)*Power(xj, 2LL) +

            32810LL*Power(r, 4LL)*Power(xj, 4LL) + 22LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  7980LL*Power(xi, 16LL)*Power(xj, 8LL)*

                                  (15696909LL - 494343LL*Power(r, 2LL)*Power(xj, 2LL) -

            4182LL*Power(r, 4LL)*Power(xj, 4LL) + 34LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  2LL*r*Power(xi, 9LL)*Power(xj, 16LL)*

                                  (19433295LL - 690795LL*Power(r, 2LL)*Power(xj, 2LL) +

            55251LL*Power(r, 4LL)*Power(xj, 4LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  6LL*r*Power(xi, 25LL)*(-8546580LL - 11329605LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                         24003LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  6LL*r*Power(xi, 13LL)*Power(xj, 12LL)*

                                  (-2361196215LL - 54738810LL*Power(r, 2LL)*Power(xj, 2LL) +

            388626LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  38LL*r*Power(xi, 15LL)*Power(xj, 10LL)*

                                  (808181955LL - 17168130LL*Power(r, 2LL)*Power(xj, 2LL) -

            32130LL*Power(r, 4LL)*Power(xj, 4LL) + 106LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  84LL*Power(xi, 10LL)*Power(xj, 14LL)*

                                  (3168630LL + 683145LL*Power(r, 2LL)*Power(xj, 2LL) +

            54315LL*Power(r, 4LL)*Power(xj, 4LL) + 193LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  19LL*r*Power(xi, 17LL)*Power(xj, 8LL)*

                                  (-2525985LL + 33479460LL*Power(r, 2LL)*Power(xj, 2LL) -

            406980LL*Power(r, 4LL)*Power(xj, 4LL) + 272LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  84LL*Power(xi, 12LL)*Power(xj, 12LL)*

                                  (-88925130LL - 19869345LL*Power(r, 2LL)*Power(xj, 2LL) -

            235790LL*Power(r, 4LL)*Power(xj, 4LL) + 643LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  210LL*Power(xi, 18LL)*Power(xj, 6LL)*

                                  (-496605582LL + 32638599LL*Power(r, 2LL)*Power(xj, 2LL) -

            564604LL*Power(r, 4LL)*Power(xj, 4LL) + 1292LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  42LL*Power(xi, 20LL)*Power(xj, 4LL)*

                                  (-777723210LL - 46394505LL*Power(r, 2LL)*Power(xj, 2LL) +

            625670LL*Power(r, 4LL)*Power(xj, 4LL) + 1292LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  42LL*Power(xi, 24LL)*(-1918620LL - 11344995LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                        323070LL*Power(r, 4LL)*Power(xj, 4LL) + 2114LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  r*Power(xi, 23LL)*Power(xj, 2LL)*

                                  (1919335635LL + 275096430LL*Power(r, 2LL)*Power(xj, 2LL) -

            3302586LL*Power(r, 4LL)*Power(xj, 4LL) + 4028LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  r*Power(xi, 21LL)*Power(xj, 4LL)*

                                  (-14708379735LL + 255168270LL*Power(r, 2LL)*Power(xj, 2LL) -

            2899134LL*Power(r, 4LL)*Power(xj, 4LL) + 5168LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  42LL*Power(xi, 22LL)*Power(xj, 2LL)*

                                  (81654210LL + 66273255LL*Power(r, 2LL)*Power(xj, 2LL) -

            1203870LL*Power(r, 4LL)*Power(xj, 4LL) + 5206LL*Power(r, 6LL)*Power(xj, 6LL))) -

                                          2LL*exp(2LL*r*xi)*Power(xi, 10LL)*

                                          (21318LL*Power(xi, 14LL)*Power(xj, 14LL)*

                                  (-3146850LL + 4890375LL*r*xj - 24522750LL*Power(r, 2LL)*Power(xj, 2LL) +

            12162150LL*Power(r, 3LL)*Power(xj, 3LL) -

            1549800LL*Power(r, 4LL)*Power(xj, 4LL) -

            1615950LL*Power(r, 5LL)*Power(xj, 5LL) -

            185220LL*Power(r, 6LL)*Power(xj, 6LL) + 12240LL*Power(r, 7LL)*Power(xj, 7LL) +

            3960LL*Power(r, 8LL)*Power(xj, 8LL) + 300LL*Power(r, 9LL)*Power(xj, 9LL) +

            8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                  3LL*Power(xi, 24LL)*Power(xj, 4LL)*

                                  (53326350LL + 97764975LL*r*xj + 88877250LL*Power(r, 2LL)*Power(xj, 2LL) +

            53326350LL*Power(r, 3LL)*Power(xj, 3LL) +

            23700600LL*Power(r, 4LL)*Power(xj, 4LL) +

            8295210LL*Power(r, 5LL)*Power(xj, 5LL) +

            2370060LL*Power(r, 6LL)*Power(xj, 6LL) +

            564300LL*Power(r, 7LL)*Power(xj, 7LL) + 112860LL*Power(r, 8LL)*Power(xj, 8LL) +

            22440LL*Power(r, 9LL)*Power(xj, 9LL) + 1056LL*Power(r, 10LL)*Power(xj, 10LL) -

            20LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  4LL*Power(xj, 28LL)*(13749310575LL + 13749310575LL*r*xj +

                                                       6547290750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       1964187225LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                       413513100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                       64324260LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                       7567560LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                       675675LL*Power(r, 7LL)*Power(xj, 7LL) + 45045LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                       2145LL*Power(r, 9LL)*Power(xj, 9LL) + 66LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                       Power(r, 11LL)*Power(xj, 11LL)) -

                                  1254LL*Power(xi, 16LL)*Power(xj, 12LL)*

                                  (-20241900LL - 38315025LL*r*xj - 21687750LL*Power(r, 2LL)*Power(xj, 2LL) -

            50122800LL*Power(r, 3LL)*Power(xj, 3LL) +

            14137200LL*Power(r, 4LL)*Power(xj, 4LL) -

            5853330LL*Power(r, 5LL)*Power(xj, 5LL) -

            2687580LL*Power(r, 6LL)*Power(xj, 6LL) -

            208530LL*Power(r, 7LL)*Power(xj, 7LL) + 19530LL*Power(r, 8LL)*Power(xj, 8LL) +

            3630LL*Power(r, 9LL)*Power(xj, 9LL) + 172LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  627LL*Power(xi, 12LL)*Power(xj, 16LL)*

                                  (-1240964550LL + 4740389325LL*r*xj -

            3311818650LL*Power(r, 2LL)*Power(xj, 2LL) +

            134804250LL*Power(r, 3LL)*Power(xj, 3LL) +

            407673000LL*Power(r, 4LL)*Power(xj, 4LL) +

            58641030LL*Power(r, 5LL)*Power(xj, 5LL) -

            3549420LL*Power(r, 6LL)*Power(xj, 6LL) -

            1641060LL*Power(r, 7LL)*Power(xj, 7LL) -

            167940LL*Power(r, 8LL)*Power(xj, 8LL) - 6990LL*Power(r, 9LL)*Power(xj, 9LL) -

            36LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  Power(xi, 28LL)*(935550LL + 1715175LL*r*xj +

                                                   1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                   330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                   4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  2LL*Power(xi, 2LL)*Power(xj, 26LL)*

                                  (-937068397650LL - 815439881025LL*r*xj -

            332904552750LL*Power(r, 2LL)*Power(xj, 2LL) -

            84006776700LL*Power(r, 3LL)*Power(xj, 3LL) -

            14504767200LL*Power(r, 4LL)*Power(xj, 4LL) -

            1786235220LL*Power(r, 5LL)*Power(xj, 5LL) -

            157754520LL*Power(r, 6LL)*Power(xj, 6LL) -

            9667350LL*Power(r, 7LL)*Power(xj, 7LL) -

            367290LL*Power(r, 8LL)*Power(xj, 8LL) - 5115LL*Power(r, 9LL)*Power(xj, 9LL) +

            198LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  6LL*Power(xi, 4LL)*Power(xj, 24LL)*

                                  (-2262441500550LL - 1503711230175LL*r*xj -

            426178264050LL*Power(r, 2LL)*Power(xj, 2LL) -

            60134347350LL*Power(r, 3LL)*Power(xj, 3LL) -

            2014551000LL*Power(r, 4LL)*Power(xj, 4LL) +

            846111420LL*Power(r, 5LL)*Power(xj, 5LL) +

            184864680LL*Power(r, 6LL)*Power(xj, 6LL) +

            20183130LL*Power(r, 7LL)*Power(xj, 7LL) +

            1367190LL*Power(r, 8LL)*Power(xj, 8LL) + 57255LL*Power(r, 9LL)*Power(xj, 9LL) +

            1298LL*Power(r, 10LL)*Power(xj, 10LL) + 10LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  Power(xi, 26LL)*Power(xj, 2LL)*

                                  (17775450LL + 32588325LL*r*xj + 29625750LL*Power(r, 2LL)*Power(xj, 2LL) +

            17775450LL*Power(r, 3LL)*Power(xj, 3LL) +

            7900200LL*Power(r, 4LL)*Power(xj, 4LL) +

            2765070LL*Power(r, 5LL)*Power(xj, 5LL) +

            790020LL*Power(r, 6LL)*Power(xj, 6LL) + 188100LL*Power(r, 7LL)*Power(xj, 7LL) +

            37620LL*Power(r, 8LL)*Power(xj, 8LL) + 6270LL*Power(r, 9LL)*Power(xj, 9LL) +

            836LL*Power(r, 10LL)*Power(xj, 10LL) + 16LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  9LL*Power(xi, 22LL)*Power(xj, 6LL)*

                                  (-100727550LL - 184667175LL*r*xj - 167879250LL*Power(r, 2LL)*Power(xj, 2LL) -

            100727550LL*Power(r, 3LL)*Power(xj, 3LL) -

            44767800LL*Power(r, 4LL)*Power(xj, 4LL) -

            15668730LL*Power(r, 5LL)*Power(xj, 5LL) -

            4476780LL*Power(r, 6LL)*Power(xj, 6LL) -

            971520LL*Power(r, 7LL)*Power(xj, 7LL) - 307560LL*Power(r, 8LL)*Power(xj, 8LL) -

            27060LL*Power(r, 9LL)*Power(xj, 9LL) + 264LL*Power(r, 10LL)*Power(xj, 10LL) +

            64LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  9LL*Power(xi, 6LL)*Power(xj, 22LL)*

                                  (3452543428950LL + 1097992509075LL*r*xj -

            101420792550LL*Power(r, 2LL)*Power(xj, 2LL) -

            110557373850LL*Power(r, 3LL)*Power(xj, 3LL) -

            24909330600LL*Power(r, 4LL)*Power(xj, 4LL) -

            2686726350LL*Power(r, 5LL)*Power(xj, 5LL) -

            93485700LL*Power(r, 6LL)*Power(xj, 6LL) +

            12941280LL*Power(r, 7LL)*Power(xj, 7LL) +

            2081640LL*Power(r, 8LL)*Power(xj, 8LL) +

            137940LL*Power(r, 9LL)*Power(xj, 9LL) + 4664LL*Power(r, 10LL)*Power(xj, 10LL) +

            64LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  22LL*Power(xi, 20LL)*Power(xj, 8LL)*

                                  (-164826900LL - 302182650LL*r*xj - 274711500LL*Power(r, 2LL)*Power(xj, 2LL) -

            164826900LL*Power(r, 3LL)*Power(xj, 3LL) -

            73256400LL*Power(r, 4LL)*Power(xj, 4LL) -

            26991090LL*Power(r, 5LL)*Power(xj, 5LL) -

            4622940LL*Power(r, 6LL)*Power(xj, 6LL) -

            2941110LL*Power(r, 7LL)*Power(xj, 7LL) -

            438930LL*Power(r, 8LL)*Power(xj, 8LL) - 5505LL*Power(r, 9LL)*Power(xj, 9LL) +

            2082LL*Power(r, 10LL)*Power(xj, 10LL) + 82LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  22LL*Power(xi, 18LL)*Power(xj, 10LL)*

                                  (-494480700LL - 906547950LL*r*xj - 824134500LL*Power(r, 2LL)*Power(xj, 2LL) -

            475684650LL*Power(r, 3LL)*Power(xj, 3LL) -

            294953400LL*Power(r, 4LL)*Power(xj, 4LL) +

            2663010LL*Power(r, 5LL)*Power(xj, 5LL) -

            40797540LL*Power(r, 6LL)*Power(xj, 6LL) -

            10248390LL*Power(r, 7LL)*Power(xj, 7LL) -

            434610LL*Power(r, 8LL)*Power(xj, 8LL) + 65865LL*Power(r, 9LL)*Power(xj, 9LL) +

            6366LL*Power(r, 10LL)*Power(xj, 10LL) + 136LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  11LL*Power(xi, 8LL)*Power(xj, 20LL)*

                                  (-2338604626050LL + 656001834075LL*r*xj +

            504510561450LL*Power(r, 2LL)*Power(xj, 2LL) +

            51560967150LL*Power(r, 3LL)*Power(xj, 3LL) -

            15574998600LL*Power(r, 4LL)*Power(xj, 4LL) -

            5055778350LL*Power(r, 5LL)*Power(xj, 5LL) -

            626213700LL*Power(r, 6LL)*Power(xj, 6LL) -

            34768620LL*Power(r, 7LL)*Power(xj, 7LL) +

            207540LL*Power(r, 8LL)*Power(xj, 8LL) + 150240LL*Power(r, 9LL)*Power(xj, 9LL) +

            8464LL*Power(r, 10LL)*Power(xj, 10LL) + 164LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  11LL*Power(xi, 10LL)*Power(xj, 18LL)*

                                  (742805182350LL - 933111659025LL*r*xj +

            57080542050LL*Power(r, 2LL)*Power(xj, 2LL) +

            129505209750LL*Power(r, 3LL)*Power(xj, 3LL) +

            19066887000LL*Power(r, 4LL)*Power(xj, 4LL) -

            1817573310LL*Power(r, 5LL)*Power(xj, 5LL) -

            810647460LL*Power(r, 6LL)*Power(xj, 6LL) -

            97669980LL*Power(r, 7LL)*Power(xj, 7LL) -

            5173020LL*Power(r, 8LL)*Power(xj, 8LL) - 37770LL*Power(r, 9LL)*Power(xj, 9LL) +

            8212LL*Power(r, 10LL)*Power(xj, 10LL) + 272LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (935550LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 19LL)*Power(xi + xj, 18LL)) -

                (3742200LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 19LL) +

                 495LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-4032LL*Power(r, 5LL)*Power(xi, 30LL) - 84LL*Power(r, 6LL)*Power(xi, 31LL) -

        20328LL*Power(r, 5LL)*Power(xi, 28LL)*Power(xj, 2LL) -

        272LL*Power(r, 6LL)*Power(xi, 29LL)*Power(xj, 2LL) + 6615LL*xi*Power(xj, 24LL) +

        11340LL*r*Power(xi, 2LL)*Power(xj, 24LL) +

        6300LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 24LL) -

        680LL*Power(r, 4LL)*Power(xi, 29LL)*(126LL + Power(r, 2LL)*Power(xj, 2LL)) +

        315LL*Power(xi, 3LL)*Power(xj, 22LL)*(-399LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

        336LL*Power(r, 3LL)*Power(xi, 28LL)*(3060LL + 121LL*Power(r, 2LL)*Power(xj, 2LL)) +

        630LL*Power(xi, 4LL)*Power(xj, 20LL)*

        (-342LL*r*Power(xj, 2LL) + 8LL*Power(r, 3LL)*Power(xj, 4LL)) +

        63LL*r*Power(xi, 5LL)*Power(xj, 20LL)*

        (-1900LL*r*Power(xj, 2LL) + 24LL*Power(r, 3LL)*Power(xj, 4LL)) +

        84LL*Power(r, 2LL)*Power(xi, 26LL)*

        (-143070LL*r*Power(xj, 2LL) + 716LL*Power(r, 3LL)*Power(xj, 4LL)) -

        63LL*r*Power(xi, 19LL)*Power(xj, 6LL)*

        (-29796180LL*r*Power(xj, 2LL) + 315248LL*Power(r, 3LL)*Power(xj, 4LL)) +

        63LL*Power(xi, 5LL)*Power(xj, 20LL)*

        (17955LL - 950LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL)) +

        168LL*r*Power(xi, 26LL)*(-174420LL - 71535LL*Power(r, 2LL)*Power(xj, 2LL) +

                                 179LL*Power(r, 4LL)*Power(xj, 4LL)) -

        63LL*Power(xi, 19LL)*Power(xj, 6LL)*

        (468377895LL - 14898090LL*Power(r, 2LL)*Power(xj, 2LL) +

            78812LL*Power(r, 4LL)*Power(xj, 4LL)) +

        Power(xi, 27LL)*(-7325640LL*Power(r, 2LL) -

                         1639890LL*Power(r, 4LL)*Power(xj, 2LL) + 3472LL*Power(r, 6LL)*Power(xj, 4LL)) +

        2LL*r*Power(xi, 11LL)*Power(xj, 14LL)*

        (112733460LL*r*Power(xj, 2LL) + 1534428LL*Power(r, 3LL)*Power(xj, 4LL) -

            1488LL*Power(r, 5LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 6LL)*Power(xj, 18LL)*

        (46170LL*r*Power(xj, 2LL) - 2280LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) -

        798LL*Power(xi, 8LL)*Power(xj, 16LL)*

        (13770LL*r*Power(xj, 2LL) - 1080LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) +

        3LL*r*Power(xi, 7LL)*Power(xj, 18LL)*

        (359100LL*r*Power(xj, 2LL) - 9576LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) +

        1596LL*Power(xi, 14LL)*Power(xj, 10LL)*

        (-4817970LL*r*Power(xj, 2LL) + 131240LL*Power(r, 3LL)*Power(xj, 4LL) +

            132LL*Power(r, 5LL)*Power(xj, 6LL)) -

        7980LL*Power(xi, 16LL)*Power(xj, 8LL)*

        (-988686LL*r*Power(xj, 2LL) - 16728LL*Power(r, 3LL)*Power(xj, 4LL) +

            204LL*Power(r, 5LL)*Power(xj, 6LL)) +

        2LL*r*Power(xi, 9LL)*Power(xj, 16LL)*

        (-1381590LL*r*Power(xj, 2LL) + 221004LL*Power(r, 3LL)*Power(xj, 4LL) +

            408LL*Power(r, 5LL)*Power(xj, 6LL)) +

        6LL*r*Power(xi, 25LL)*(-22659210LL*r*Power(xj, 2LL) -

                               96012LL*Power(r, 3LL)*Power(xj, 4LL) + 552LL*Power(r, 5LL)*Power(xj, 6LL)) -

        6LL*r*Power(xi, 13LL)*Power(xj, 12LL)*

        (-109477620LL*r*Power(xj, 2LL) + 1554504LL*Power(r, 3LL)*Power(xj, 4LL) +

            552LL*Power(r, 5LL)*Power(xj, 6LL)) +

        38LL*r*Power(xi, 15LL)*Power(xj, 10LL)*

        (-34336260LL*r*Power(xj, 2LL) - 128520LL*Power(r, 3LL)*Power(xj, 4LL) +

            636LL*Power(r, 5LL)*Power(xj, 6LL)) -

        84LL*Power(xi, 10LL)*Power(xj, 14LL)*

        (1366290LL*r*Power(xj, 2LL) + 217260LL*Power(r, 3LL)*Power(xj, 4LL) +

            1158LL*Power(r, 5LL)*Power(xj, 6LL)) -

        19LL*r*Power(xi, 17LL)*Power(xj, 8LL)*

        (66958920LL*r*Power(xj, 2LL) - 1627920LL*Power(r, 3LL)*Power(xj, 4LL) +

            1632LL*Power(r, 5LL)*Power(xj, 6LL)) +

        84LL*Power(xi, 12LL)*Power(xj, 12LL)*

        (-39738690LL*r*Power(xj, 2LL) - 943160LL*Power(r, 3LL)*Power(xj, 4LL) +

            3858LL*Power(r, 5LL)*Power(xj, 6LL)) +

        210LL*Power(xi, 18LL)*Power(xj, 6LL)*

        (65277198LL*r*Power(xj, 2LL) - 2258416LL*Power(r, 3LL)*Power(xj, 4LL) +

            7752LL*Power(r, 5LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 20LL)*Power(xj, 4LL)*

        (-92789010LL*r*Power(xj, 2LL) + 2502680LL*Power(r, 3LL)*Power(xj, 4LL) +

            7752LL*Power(r, 5LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 24LL)*(-22689990LL*r*Power(xj, 2LL) -

                              1292280LL*Power(r, 3LL)*Power(xj, 4LL) + 12684LL*Power(r, 5LL)*Power(xj, 6LL)) -

        r*Power(xi, 23LL)*Power(xj, 2LL)*

        (550192860LL*r*Power(xj, 2LL) - 13210344LL*Power(r, 3LL)*Power(xj, 4LL) +

            24168LL*Power(r, 5LL)*Power(xj, 6LL)) +

        r*Power(xi, 21LL)*Power(xj, 4LL)*

        (510336540LL*r*Power(xj, 2LL) - 11596536LL*Power(r, 3LL)*Power(xj, 4LL) +

            31008LL*Power(r, 5LL)*Power(xj, 6LL)) -

        42LL*Power(xi, 22LL)*Power(xj, 2LL)*

        (132546510LL*r*Power(xj, 2LL) - 4815480LL*Power(r, 3LL)*Power(xj, 4LL) +

            31236LL*Power(r, 5LL)*Power(xj, 6LL)) +

        2LL*Power(xi, 11LL)*Power(xj, 14LL)*

        (613624095LL + 56366730LL*Power(r, 2LL)*Power(xj, 2LL) +

            383607LL*Power(r, 4LL)*Power(xj, 4LL) - 248LL*Power(r, 6LL)*Power(xj, 6LL)) +

        3LL*Power(xi, 7LL)*Power(xj, 18LL)*

        (-2136645LL + 179550LL*Power(r, 2LL)*Power(xj, 2LL) -

            2394LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

        2LL*Power(xi, 9LL)*Power(xj, 16LL)*

        (19433295LL - 690795LL*Power(r, 2LL)*Power(xj, 2LL) +

            55251LL*Power(r, 4LL)*Power(xj, 4LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL)) +

        6LL*Power(xi, 25LL)*(-8546580LL - 11329605LL*Power(r, 2LL)*Power(xj, 2LL) -

                             24003LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

        6LL*Power(xi, 13LL)*Power(xj, 12LL)*

        (-2361196215LL - 54738810LL*Power(r, 2LL)*Power(xj, 2LL) +

            388626LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) +

        38LL*Power(xi, 15LL)*Power(xj, 10LL)*

        (808181955LL - 17168130LL*Power(r, 2LL)*Power(xj, 2LL) -

            32130LL*Power(r, 4LL)*Power(xj, 4LL) + 106LL*Power(r, 6LL)*Power(xj, 6LL)) -

        19LL*Power(xi, 17LL)*Power(xj, 8LL)*

        (-2525985LL + 33479460LL*Power(r, 2LL)*Power(xj, 2LL) -

            406980LL*Power(r, 4LL)*Power(xj, 4LL) + 272LL*Power(r, 6LL)*Power(xj, 6LL)) -

        Power(xi, 23LL)*Power(xj, 2LL)*

        (1919335635LL + 275096430LL*Power(r, 2LL)*Power(xj, 2LL) -

            3302586LL*Power(r, 4LL)*Power(xj, 4LL) + 4028LL*Power(r, 6LL)*Power(xj, 6LL)) +

        Power(xi, 21LL)*Power(xj, 4LL)*

        (-14708379735LL + 255168270LL*Power(r, 2LL)*Power(xj, 2LL) -

            2899134LL*Power(r, 4LL)*Power(xj, 4LL) + 5168LL*Power(r, 6LL)*Power(xj, 6LL))) +

                 990LL*exp(2LL*r*xj)*Power(xj, 15LL)*

                 (-672LL*Power(r, 6LL)*Power(xi, 30LL) - 12LL*Power(r, 7LL)*Power(xi, 31LL) +

        3780LL*Power(xj, 24LL) + 6615LL*r*xi*Power(xj, 24LL) -

        136LL*Power(r, 5LL)*Power(xi, 29LL)*(126LL + Power(r, 2LL)*Power(xj, 2LL)) +

        1890LL*Power(xi, 2LL)*Power(xj, 22LL)*(-38LL + 3LL*Power(r, 2LL)*Power(xj, 2LL)) +

        315LL*r*Power(xi, 3LL)*Power(xj, 22LL)*(-399LL + 10LL*Power(r, 2LL)*Power(xj, 2LL)) -

        84LL*Power(r, 4LL)*Power(xi, 28LL)*(3060LL + 121LL*Power(r, 2LL)*Power(xj, 2LL)) +

        630LL*Power(xi, 4LL)*Power(xj, 20LL)*

        (1026LL - 171LL*Power(r, 2LL)*Power(xj, 2LL) + 2LL*Power(r, 4LL)*Power(xj, 4LL)) +

        63LL*r*Power(xi, 5LL)*Power(xj, 20LL)*

        (17955LL - 950LL*Power(r, 2LL)*Power(xj, 2LL) + 6LL*Power(r, 4LL)*Power(xj, 4LL)) +

        84LL*Power(r, 2LL)*Power(xi, 26LL)*

        (-174420LL - 71535LL*Power(r, 2LL)*Power(xj, 2LL) +

            179LL*Power(r, 4LL)*Power(xj, 4LL)) -

        63LL*r*Power(xi, 19LL)*Power(xj, 6LL)*

        (468377895LL - 14898090LL*Power(r, 2LL)*Power(xj, 2LL) +

            78812LL*Power(r, 4LL)*Power(xj, 4LL)) +

        Power(xi, 27LL)*(-2441880LL*Power(r, 3LL) -

                         327978LL*Power(r, 5LL)*Power(xj, 2LL) + 496LL*Power(r, 7LL)*Power(xj, 4LL)) +

        2LL*r*Power(xi, 11LL)*Power(xj, 14LL)*

        (613624095LL + 56366730LL*Power(r, 2LL)*Power(xj, 2LL) +

            383607LL*Power(r, 4LL)*Power(xj, 4LL) - 248LL*Power(r, 6LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 6LL)*Power(xj, 18LL)*

        (-87210LL + 23085LL*Power(r, 2LL)*Power(xj, 2LL) -

            570LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

        798LL*Power(xi, 8LL)*Power(xj, 16LL)*

        (-18360LL + 6885LL*Power(r, 2LL)*Power(xj, 2LL) -

            270LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

        3LL*r*Power(xi, 7LL)*Power(xj, 18LL)*

        (-2136645LL + 179550LL*Power(r, 2LL)*Power(xj, 2LL) -

            2394LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) +

        1596LL*Power(xi, 14LL)*Power(xj, 10LL)*

        (-34484670LL - 2408985LL*Power(r, 2LL)*Power(xj, 2LL) +

            32810LL*Power(r, 4LL)*Power(xj, 4LL) + 22LL*Power(r, 6LL)*Power(xj, 6LL)) -

        7980LL*Power(xi, 16LL)*Power(xj, 8LL)*

        (15696909LL - 494343LL*Power(r, 2LL)*Power(xj, 2LL) -

            4182LL*Power(r, 4LL)*Power(xj, 4LL) + 34LL*Power(r, 6LL)*Power(xj, 6LL)) +

        2LL*r*Power(xi, 9LL)*Power(xj, 16LL)*

        (19433295LL - 690795LL*Power(r, 2LL)*Power(xj, 2LL) +

            55251LL*Power(r, 4LL)*Power(xj, 4LL) + 68LL*Power(r, 6LL)*Power(xj, 6LL)) +

        6LL*r*Power(xi, 25LL)*(-8546580LL - 11329605LL*Power(r, 2LL)*Power(xj, 2LL) -

                               24003LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) -

        6LL*r*Power(xi, 13LL)*Power(xj, 12LL)*

        (-2361196215LL - 54738810LL*Power(r, 2LL)*Power(xj, 2LL) +

            388626LL*Power(r, 4LL)*Power(xj, 4LL) + 92LL*Power(r, 6LL)*Power(xj, 6LL)) +

        38LL*r*Power(xi, 15LL)*Power(xj, 10LL)*

        (808181955LL - 17168130LL*Power(r, 2LL)*Power(xj, 2LL) -

            32130LL*Power(r, 4LL)*Power(xj, 4LL) + 106LL*Power(r, 6LL)*Power(xj, 6LL)) -

        84LL*Power(xi, 10LL)*Power(xj, 14LL)*

        (3168630LL + 683145LL*Power(r, 2LL)*Power(xj, 2LL) +

            54315LL*Power(r, 4LL)*Power(xj, 4LL) + 193LL*Power(r, 6LL)*Power(xj, 6LL)) -

        19LL*r*Power(xi, 17LL)*Power(xj, 8LL)*

        (-2525985LL + 33479460LL*Power(r, 2LL)*Power(xj, 2LL) -

            406980LL*Power(r, 4LL)*Power(xj, 4LL) + 272LL*Power(r, 6LL)*Power(xj, 6LL)) +

        84LL*Power(xi, 12LL)*Power(xj, 12LL)*

        (-88925130LL - 19869345LL*Power(r, 2LL)*Power(xj, 2LL) -

            235790LL*Power(r, 4LL)*Power(xj, 4LL) + 643LL*Power(r, 6LL)*Power(xj, 6LL)) +

        210LL*Power(xi, 18LL)*Power(xj, 6LL)*

        (-496605582LL + 32638599LL*Power(r, 2LL)*Power(xj, 2LL) -

            564604LL*Power(r, 4LL)*Power(xj, 4LL) + 1292LL*Power(r, 6LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 20LL)*Power(xj, 4LL)*

        (-777723210LL - 46394505LL*Power(r, 2LL)*Power(xj, 2LL) +

            625670LL*Power(r, 4LL)*Power(xj, 4LL) + 1292LL*Power(r, 6LL)*Power(xj, 6LL)) +

        42LL*Power(xi, 24LL)*(-1918620LL - 11344995LL*Power(r, 2LL)*Power(xj, 2LL) -

                              323070LL*Power(r, 4LL)*Power(xj, 4LL) + 2114LL*Power(r, 6LL)*Power(xj, 6LL)) -

        r*Power(xi, 23LL)*Power(xj, 2LL)*

        (1919335635LL + 275096430LL*Power(r, 2LL)*Power(xj, 2LL) -

            3302586LL*Power(r, 4LL)*Power(xj, 4LL) + 4028LL*Power(r, 6LL)*Power(xj, 6LL)) +

        r*Power(xi, 21LL)*Power(xj, 4LL)*

        (-14708379735LL + 255168270LL*Power(r, 2LL)*Power(xj, 2LL) -

            2899134LL*Power(r, 4LL)*Power(xj, 4LL) + 5168LL*Power(r, 6LL)*Power(xj, 6LL)) -

        42LL*Power(xi, 22LL)*Power(xj, 2LL)*

        (81654210LL + 66273255LL*Power(r, 2LL)*Power(xj, 2LL) -

            1203870LL*Power(r, 4LL)*Power(xj, 4LL) + 5206LL*Power(r, 6LL)*Power(xj, 6LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 10LL)*

                 (21318LL*Power(xi, 14LL)*Power(xj, 14LL)*

        (4890375LL*xj - 49045500LL*r*Power(xj, 2LL) +

            36486450LL*Power(r, 2LL)*Power(xj, 3LL) -

            6199200LL*Power(r, 3LL)*Power(xj, 4LL) -

            8079750LL*Power(r, 4LL)*Power(xj, 5LL) -

            1111320LL*Power(r, 5LL)*Power(xj, 6LL) + 85680LL*Power(r, 6LL)*Power(xj, 7LL) +

            31680LL*Power(r, 7LL)*Power(xj, 8LL) + 2700LL*Power(r, 8LL)*Power(xj, 9LL) +

            80LL*Power(r, 9LL)*Power(xj, 10LL)) +

        3LL*Power(xi, 24LL)*Power(xj, 4LL)*

        (97764975LL*xj + 177754500LL*r*Power(xj, 2LL) +

            159979050LL*Power(r, 2LL)*Power(xj, 3LL) +

            94802400LL*Power(r, 3LL)*Power(xj, 4LL) +

            41476050LL*Power(r, 4LL)*Power(xj, 5LL) +

            14220360LL*Power(r, 5LL)*Power(xj, 6LL) +

            3950100LL*Power(r, 6LL)*Power(xj, 7LL) +

            902880LL*Power(r, 7LL)*Power(xj, 8LL) + 201960LL*Power(r, 8LL)*Power(xj, 9LL) +

            10560LL*Power(r, 9LL)*Power(xj, 10LL) - 220LL*Power(r, 10LL)*Power(xj, 11LL)) -

        4LL*Power(xj, 28LL)*(13749310575LL*xj + 13094581500LL*r*Power(xj, 2LL) +

                             5892561675LL*Power(r, 2LL)*Power(xj, 3LL) +

                             1654052400LL*Power(r, 3LL)*Power(xj, 4LL) +

                             321621300LL*Power(r, 4LL)*Power(xj, 5LL) +

                             45405360LL*Power(r, 5LL)*Power(xj, 6LL) +

                             4729725LL*Power(r, 6LL)*Power(xj, 7LL) +

                             360360LL*Power(r, 7LL)*Power(xj, 8LL) + 19305LL*Power(r, 8LL)*Power(xj, 9LL) +

                             660LL*Power(r, 9LL)*Power(xj, 10LL) + 11LL*Power(r, 10LL)*Power(xj, 11LL)) -

        1254LL*Power(xi, 16LL)*Power(xj, 12LL)*

        (-38315025LL*xj - 43375500LL*r*Power(xj, 2LL) -

            150368400LL*Power(r, 2LL)*Power(xj, 3LL) +

            56548800LL*Power(r, 3LL)*Power(xj, 4LL) -

            29266650LL*Power(r, 4LL)*Power(xj, 5LL) -

            16125480LL*Power(r, 5LL)*Power(xj, 6LL) -

            1459710LL*Power(r, 6LL)*Power(xj, 7LL) +

            156240LL*Power(r, 7LL)*Power(xj, 8LL) + 32670LL*Power(r, 8LL)*Power(xj, 9LL) +

            1720LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

        627LL*Power(xi, 12LL)*Power(xj, 16LL)*

        (4740389325LL*xj - 6623637300LL*r*Power(xj, 2LL) +

            404412750LL*Power(r, 2LL)*Power(xj, 3LL) +

            1630692000LL*Power(r, 3LL)*Power(xj, 4LL) +

            293205150LL*Power(r, 4LL)*Power(xj, 5LL) -

            21296520LL*Power(r, 5LL)*Power(xj, 6LL) -

            11487420LL*Power(r, 6LL)*Power(xj, 7LL) -

            1343520LL*Power(r, 7LL)*Power(xj, 8LL) - 62910LL*Power(r, 8LL)*Power(xj, 9LL) -

            360LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) +

        Power(xi, 28LL)*(1715175LL*xj + 3118500LL*r*Power(xj, 2LL) +

                         2806650LL*Power(r, 2LL)*Power(xj, 3LL) +

                         1663200LL*Power(r, 3LL)*Power(xj, 4LL) +

                         727650LL*Power(r, 4LL)*Power(xj, 5LL) + 249480LL*Power(r, 5LL)*Power(xj, 6LL) +

                         69300LL*Power(r, 6LL)*Power(xj, 7LL) + 15840LL*Power(r, 7LL)*Power(xj, 8LL) +

                         2970LL*Power(r, 8LL)*Power(xj, 9LL) + 440LL*Power(r, 9LL)*Power(xj, 10LL) +

                         44LL*Power(r, 10LL)*Power(xj, 11LL)) +

        2LL*Power(xi, 2LL)*Power(xj, 26LL)*

        (-815439881025LL*xj - 665809105500LL*r*Power(xj, 2LL) -

            252020330100LL*Power(r, 2LL)*Power(xj, 3LL) -

            58019068800LL*Power(r, 3LL)*Power(xj, 4LL) -

            8931176100LL*Power(r, 4LL)*Power(xj, 5LL) -

            946527120LL*Power(r, 5LL)*Power(xj, 6LL) -

            67671450LL*Power(r, 6LL)*Power(xj, 7LL) -

            2938320LL*Power(r, 7LL)*Power(xj, 8LL) - 46035LL*Power(r, 8LL)*Power(xj, 9LL) +

            1980LL*Power(r, 9LL)*Power(xj, 10LL) + 88LL*Power(r, 10LL)*Power(xj, 11LL)) +

        6LL*Power(xi, 4LL)*Power(xj, 24LL)*

        (-1503711230175LL*xj - 852356528100LL*r*Power(xj, 2LL) -

            180403042050LL*Power(r, 2LL)*Power(xj, 3LL) -

            8058204000LL*Power(r, 3LL)*Power(xj, 4LL) +

            4230557100LL*Power(r, 4LL)*Power(xj, 5LL) +

            1109188080LL*Power(r, 5LL)*Power(xj, 6LL) +

            141281910LL*Power(r, 6LL)*Power(xj, 7LL) +

            10937520LL*Power(r, 7LL)*Power(xj, 8LL) +

            515295LL*Power(r, 8LL)*Power(xj, 9LL) + 12980LL*Power(r, 9LL)*Power(xj, 10LL) +

            110LL*Power(r, 10LL)*Power(xj, 11LL)) -

        Power(xi, 26LL)*Power(xj, 2LL)*

        (32588325LL*xj + 59251500LL*r*Power(xj, 2LL) +

            53326350LL*Power(r, 2LL)*Power(xj, 3LL) +

            31600800LL*Power(r, 3LL)*Power(xj, 4LL) +

            13825350LL*Power(r, 4LL)*Power(xj, 5LL) +

            4740120LL*Power(r, 5LL)*Power(xj, 6LL) +

            1316700LL*Power(r, 6LL)*Power(xj, 7LL) +

            300960LL*Power(r, 7LL)*Power(xj, 8LL) + 56430LL*Power(r, 8LL)*Power(xj, 9LL) +

            8360LL*Power(r, 9LL)*Power(xj, 10LL) + 176LL*Power(r, 10LL)*Power(xj, 11LL)) +

        9LL*Power(xi, 22LL)*Power(xj, 6LL)*

        (-184667175LL*xj - 335758500LL*r*Power(xj, 2LL) -

            302182650LL*Power(r, 2LL)*Power(xj, 3LL) -

            179071200LL*Power(r, 3LL)*Power(xj, 4LL) -

            78343650LL*Power(r, 4LL)*Power(xj, 5LL) -

            26860680LL*Power(r, 5LL)*Power(xj, 6LL) -

            6800640LL*Power(r, 6LL)*Power(xj, 7LL) -

            2460480LL*Power(r, 7LL)*Power(xj, 8LL) -

            243540LL*Power(r, 8LL)*Power(xj, 9LL) + 2640LL*Power(r, 9LL)*Power(xj, 10LL) +

            704LL*Power(r, 10LL)*Power(xj, 11LL)) -

        9LL*Power(xi, 6LL)*Power(xj, 22LL)*

        (1097992509075LL*xj - 202841585100LL*r*Power(xj, 2LL) -

            331672121550LL*Power(r, 2LL)*Power(xj, 3LL) -

            99637322400LL*Power(r, 3LL)*Power(xj, 4LL) -

            13433631750LL*Power(r, 4LL)*Power(xj, 5LL) -

            560914200LL*Power(r, 5LL)*Power(xj, 6LL) +

            90588960LL*Power(r, 6LL)*Power(xj, 7LL) +

            16653120LL*Power(r, 7LL)*Power(xj, 8LL) +

            1241460LL*Power(r, 8LL)*Power(xj, 9LL) +

            46640LL*Power(r, 9LL)*Power(xj, 10LL) + 704LL*Power(r, 10LL)*Power(xj, 11LL)) -

        22LL*Power(xi, 20LL)*Power(xj, 8LL)*

        (-302182650LL*xj - 549423000LL*r*Power(xj, 2LL) -

            494480700LL*Power(r, 2LL)*Power(xj, 3LL) -

            293025600LL*Power(r, 3LL)*Power(xj, 4LL) -

            134955450LL*Power(r, 4LL)*Power(xj, 5LL) -

            27737640LL*Power(r, 5LL)*Power(xj, 6LL) -

            20587770LL*Power(r, 6LL)*Power(xj, 7LL) -

            3511440LL*Power(r, 7LL)*Power(xj, 8LL) - 49545LL*Power(r, 8LL)*Power(xj, 9LL) +

            20820LL*Power(r, 9LL)*Power(xj, 10LL) + 902LL*Power(r, 10LL)*Power(xj, 11LL)) +

        22LL*Power(xi, 18LL)*Power(xj, 10LL)*

        (-906547950LL*xj - 1648269000LL*r*Power(xj, 2LL) -

            1427053950LL*Power(r, 2LL)*Power(xj, 3LL) -

            1179813600LL*Power(r, 3LL)*Power(xj, 4LL) +

            13315050LL*Power(r, 4LL)*Power(xj, 5LL) -

            244785240LL*Power(r, 5LL)*Power(xj, 6LL) -

            71738730LL*Power(r, 6LL)*Power(xj, 7LL) -

            3476880LL*Power(r, 7LL)*Power(xj, 8LL) +

            592785LL*Power(r, 8LL)*Power(xj, 9LL) + 63660LL*Power(r, 9LL)*Power(xj, 10LL) +

            1496LL*Power(r, 10LL)*Power(xj, 11LL)) +

        11LL*Power(xi, 8LL)*Power(xj, 20LL)*

        (656001834075LL*xj + 1009021122900LL*r*Power(xj, 2LL) +

            154682901450LL*Power(r, 2LL)*Power(xj, 3LL) -

            62299994400LL*Power(r, 3LL)*Power(xj, 4LL) -

            25278891750LL*Power(r, 4LL)*Power(xj, 5LL) -

            3757282200LL*Power(r, 5LL)*Power(xj, 6LL) -

            243380340LL*Power(r, 6LL)*Power(xj, 7LL) +

            1660320LL*Power(r, 7LL)*Power(xj, 8LL) +

            1352160LL*Power(r, 8LL)*Power(xj, 9LL) +

            84640LL*Power(r, 9LL)*Power(xj, 10LL) + 1804LL*Power(r, 10LL)*Power(xj, 11LL)) -

        11LL*Power(xi, 10LL)*Power(xj, 18LL)*

        (-933111659025LL*xj + 114161084100LL*r*Power(xj, 2LL) +

            388515629250LL*Power(r, 2LL)*Power(xj, 3LL) +

            76267548000LL*Power(r, 3LL)*Power(xj, 4LL) -

            9087866550LL*Power(r, 4LL)*Power(xj, 5LL) -

            4863884760LL*Power(r, 5LL)*Power(xj, 6LL) -

            683689860LL*Power(r, 6LL)*Power(xj, 7LL) -

            41384160LL*Power(r, 7LL)*Power(xj, 8LL) -

            339930LL*Power(r, 8LL)*Power(xj, 9LL) + 82120LL*Power(r, 9LL)*Power(xj, 10LL) +

            2992LL*Power(r, 10LL)*Power(xj, 11LL))) -

                 4LL*exp(2LL*r*xi)*Power(xi, 11LL)*

                 (21318LL*Power(xi, 14LL)*Power(xj, 14LL)*

        (-3146850LL + 4890375LL*r*xj - 24522750LL*Power(r, 2LL)*Power(xj, 2LL) +

            12162150LL*Power(r, 3LL)*Power(xj, 3LL) -

            1549800LL*Power(r, 4LL)*Power(xj, 4LL) -

            1615950LL*Power(r, 5LL)*Power(xj, 5LL) - 185220LL*Power(r, 6LL)*Power(xj, 6LL) +

            12240LL*Power(r, 7LL)*Power(xj, 7LL) + 3960LL*Power(r, 8LL)*Power(xj, 8LL) +

            300LL*Power(r, 9LL)*Power(xj, 9LL) + 8LL*Power(r, 10LL)*Power(xj, 10LL)) +

        3LL*Power(xi, 24LL)*Power(xj, 4LL)*

        (53326350LL + 97764975LL*r*xj + 88877250LL*Power(r, 2LL)*Power(xj, 2LL) +

            53326350LL*Power(r, 3LL)*Power(xj, 3LL) +

            23700600LL*Power(r, 4LL)*Power(xj, 4LL) +

            8295210LL*Power(r, 5LL)*Power(xj, 5LL) +

            2370060LL*Power(r, 6LL)*Power(xj, 6LL) + 564300LL*Power(r, 7LL)*Power(xj, 7LL) +

            112860LL*Power(r, 8LL)*Power(xj, 8LL) + 22440LL*Power(r, 9LL)*Power(xj, 9LL) +

            1056LL*Power(r, 10LL)*Power(xj, 10LL) - 20LL*Power(r, 11LL)*Power(xj, 11LL)) -

        4LL*Power(xj, 28LL)*(13749310575LL + 13749310575LL*r*xj +

                             6547290750LL*Power(r, 2LL)*Power(xj, 2LL) +

                             1964187225LL*Power(r, 3LL)*Power(xj, 3LL) +

                             413513100LL*Power(r, 4LL)*Power(xj, 4LL) +

                             64324260LL*Power(r, 5LL)*Power(xj, 5LL) +

                             7567560LL*Power(r, 6LL)*Power(xj, 6LL) + 675675LL*Power(r, 7LL)*Power(xj, 7LL) +

                             45045LL*Power(r, 8LL)*Power(xj, 8LL) + 2145LL*Power(r, 9LL)*Power(xj, 9LL) +

                             66LL*Power(r, 10LL)*Power(xj, 10LL) + Power(r, 11LL)*Power(xj, 11LL)) -

        1254LL*Power(xi, 16LL)*Power(xj, 12LL)*

        (-20241900LL - 38315025LL*r*xj - 21687750LL*Power(r, 2LL)*Power(xj, 2LL) -

            50122800LL*Power(r, 3LL)*Power(xj, 3LL) +

            14137200LL*Power(r, 4LL)*Power(xj, 4LL) -

            5853330LL*Power(r, 5LL)*Power(xj, 5LL) -

            2687580LL*Power(r, 6LL)*Power(xj, 6LL) - 208530LL*Power(r, 7LL)*Power(xj, 7LL) +

            19530LL*Power(r, 8LL)*Power(xj, 8LL) + 3630LL*Power(r, 9LL)*Power(xj, 9LL) +

            172LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

        627LL*Power(xi, 12LL)*Power(xj, 16LL)*

        (-1240964550LL + 4740389325LL*r*xj -

            3311818650LL*Power(r, 2LL)*Power(xj, 2LL) +

            134804250LL*Power(r, 3LL)*Power(xj, 3LL) +

            407673000LL*Power(r, 4LL)*Power(xj, 4LL) +

            58641030LL*Power(r, 5LL)*Power(xj, 5LL) -

            3549420LL*Power(r, 6LL)*Power(xj, 6LL) -

            1641060LL*Power(r, 7LL)*Power(xj, 7LL) - 167940LL*Power(r, 8LL)*Power(xj, 8LL) -

            6990LL*Power(r, 9LL)*Power(xj, 9LL) - 36LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

        Power(xi, 28LL)*(935550LL + 1715175LL*r*xj +

                         1559250LL*Power(r, 2LL)*Power(xj, 2LL) + 935550LL*Power(r, 3LL)*Power(xj, 3LL) +

                         415800LL*Power(r, 4LL)*Power(xj, 4LL) + 145530LL*Power(r, 5LL)*Power(xj, 5LL) +

                         41580LL*Power(r, 6LL)*Power(xj, 6LL) + 9900LL*Power(r, 7LL)*Power(xj, 7LL) +

                         1980LL*Power(r, 8LL)*Power(xj, 8LL) + 330LL*Power(r, 9LL)*Power(xj, 9LL) +

                         44LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

        2LL*Power(xi, 2LL)*Power(xj, 26LL)*

        (-937068397650LL - 815439881025LL*r*xj -

            332904552750LL*Power(r, 2LL)*Power(xj, 2LL) -

            84006776700LL*Power(r, 3LL)*Power(xj, 3LL) -

            14504767200LL*Power(r, 4LL)*Power(xj, 4LL) -

            1786235220LL*Power(r, 5LL)*Power(xj, 5LL) -

            157754520LL*Power(r, 6LL)*Power(xj, 6LL) -

            9667350LL*Power(r, 7LL)*Power(xj, 7LL) - 367290LL*Power(r, 8LL)*Power(xj, 8LL) -

            5115LL*Power(r, 9LL)*Power(xj, 9LL) + 198LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

        6LL*Power(xi, 4LL)*Power(xj, 24LL)*

        (-2262441500550LL - 1503711230175LL*r*xj -

            426178264050LL*Power(r, 2LL)*Power(xj, 2LL) -

            60134347350LL*Power(r, 3LL)*Power(xj, 3LL) -

            2014551000LL*Power(r, 4LL)*Power(xj, 4LL) +

            846111420LL*Power(r, 5LL)*Power(xj, 5LL) +

            184864680LL*Power(r, 6LL)*Power(xj, 6LL) +

            20183130LL*Power(r, 7LL)*Power(xj, 7LL) +

            1367190LL*Power(r, 8LL)*Power(xj, 8LL) + 57255LL*Power(r, 9LL)*Power(xj, 9LL) +

            1298LL*Power(r, 10LL)*Power(xj, 10LL) + 10LL*Power(r, 11LL)*Power(xj, 11LL)) -

        Power(xi, 26LL)*Power(xj, 2LL)*

        (17775450LL + 32588325LL*r*xj + 29625750LL*Power(r, 2LL)*Power(xj, 2LL) +

            17775450LL*Power(r, 3LL)*Power(xj, 3LL) +

            7900200LL*Power(r, 4LL)*Power(xj, 4LL) +

            2765070LL*Power(r, 5LL)*Power(xj, 5LL) + 790020LL*Power(r, 6LL)*Power(xj, 6LL) +

            188100LL*Power(r, 7LL)*Power(xj, 7LL) + 37620LL*Power(r, 8LL)*Power(xj, 8LL) +

            6270LL*Power(r, 9LL)*Power(xj, 9LL) + 836LL*Power(r, 10LL)*Power(xj, 10LL) +

            16LL*Power(r, 11LL)*Power(xj, 11LL)) +

        9LL*Power(xi, 22LL)*Power(xj, 6LL)*

        (-100727550LL - 184667175LL*r*xj - 167879250LL*Power(r, 2LL)*Power(xj, 2LL) -

            100727550LL*Power(r, 3LL)*Power(xj, 3LL) -

            44767800LL*Power(r, 4LL)*Power(xj, 4LL) -

            15668730LL*Power(r, 5LL)*Power(xj, 5LL) -

            4476780LL*Power(r, 6LL)*Power(xj, 6LL) - 971520LL*Power(r, 7LL)*Power(xj, 7LL) -

            307560LL*Power(r, 8LL)*Power(xj, 8LL) - 27060LL*Power(r, 9LL)*Power(xj, 9LL) +

            264LL*Power(r, 10LL)*Power(xj, 10LL) + 64LL*Power(r, 11LL)*Power(xj, 11LL)) -

        9LL*Power(xi, 6LL)*Power(xj, 22LL)*

        (3452543428950LL + 1097992509075LL*r*xj -

            101420792550LL*Power(r, 2LL)*Power(xj, 2LL) -

            110557373850LL*Power(r, 3LL)*Power(xj, 3LL) -

            24909330600LL*Power(r, 4LL)*Power(xj, 4LL) -

            2686726350LL*Power(r, 5LL)*Power(xj, 5LL) -

            93485700LL*Power(r, 6LL)*Power(xj, 6LL) +

            12941280LL*Power(r, 7LL)*Power(xj, 7LL) +

            2081640LL*Power(r, 8LL)*Power(xj, 8LL) + 137940LL*Power(r, 9LL)*Power(xj, 9LL) +

            4664LL*Power(r, 10LL)*Power(xj, 10LL) + 64LL*Power(r, 11LL)*Power(xj, 11LL)) -

        22LL*Power(xi, 20LL)*Power(xj, 8LL)*

        (-164826900LL - 302182650LL*r*xj - 274711500LL*Power(r, 2LL)*Power(xj, 2LL) -

            164826900LL*Power(r, 3LL)*Power(xj, 3LL) -

            73256400LL*Power(r, 4LL)*Power(xj, 4LL) -

            26991090LL*Power(r, 5LL)*Power(xj, 5LL) -

            4622940LL*Power(r, 6LL)*Power(xj, 6LL) -

            2941110LL*Power(r, 7LL)*Power(xj, 7LL) - 438930LL*Power(r, 8LL)*Power(xj, 8LL) -

            5505LL*Power(r, 9LL)*Power(xj, 9LL) + 2082LL*Power(r, 10LL)*Power(xj, 10LL) +

            82LL*Power(r, 11LL)*Power(xj, 11LL)) +

        22LL*Power(xi, 18LL)*Power(xj, 10LL)*

        (-494480700LL - 906547950LL*r*xj - 824134500LL*Power(r, 2LL)*Power(xj, 2LL) -

            475684650LL*Power(r, 3LL)*Power(xj, 3LL) -

            294953400LL*Power(r, 4LL)*Power(xj, 4LL) +

            2663010LL*Power(r, 5LL)*Power(xj, 5LL) -

            40797540LL*Power(r, 6LL)*Power(xj, 6LL) -

            10248390LL*Power(r, 7LL)*Power(xj, 7LL) -

            434610LL*Power(r, 8LL)*Power(xj, 8LL) + 65865LL*Power(r, 9LL)*Power(xj, 9LL) +

            6366LL*Power(r, 10LL)*Power(xj, 10LL) + 136LL*Power(r, 11LL)*Power(xj, 11LL)) +

        11LL*Power(xi, 8LL)*Power(xj, 20LL)*

        (-2338604626050LL + 656001834075LL*r*xj +

            504510561450LL*Power(r, 2LL)*Power(xj, 2LL) +

            51560967150LL*Power(r, 3LL)*Power(xj, 3LL) -

            15574998600LL*Power(r, 4LL)*Power(xj, 4LL) -

            5055778350LL*Power(r, 5LL)*Power(xj, 5LL) -

            626213700LL*Power(r, 6LL)*Power(xj, 6LL) -

            34768620LL*Power(r, 7LL)*Power(xj, 7LL) +

            207540LL*Power(r, 8LL)*Power(xj, 8LL) + 150240LL*Power(r, 9LL)*Power(xj, 9LL) +

            8464LL*Power(r, 10LL)*Power(xj, 10LL) + 164LL*Power(r, 11LL)*Power(xj, 11LL)) -

        11LL*Power(xi, 10LL)*Power(xj, 18LL)*

        (742805182350LL - 933111659025LL*r*xj +

            57080542050LL*Power(r, 2LL)*Power(xj, 2LL) +

            129505209750LL*Power(r, 3LL)*Power(xj, 3LL) +

            19066887000LL*Power(r, 4LL)*Power(xj, 4LL) -

            1817573310LL*Power(r, 5LL)*Power(xj, 5LL) -

            810647460LL*Power(r, 6LL)*Power(xj, 6LL) -

            97669980LL*Power(r, 7LL)*Power(xj, 7LL) -

            5173020LL*Power(r, 8LL)*Power(xj, 8LL) - 37770LL*Power(r, 9LL)*Power(xj, 9LL) +

            8212LL*Power(r, 10LL)*Power(xj, 10LL) + 272LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (1.8711e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 19LL)*Power(xi + xj, 19LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_6S_4S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_4S_6S(r, xj, xi);
}
