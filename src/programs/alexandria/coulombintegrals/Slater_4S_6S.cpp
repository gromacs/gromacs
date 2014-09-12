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

cl_R Slater_4S_6S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (82633LL*xi)/524288LL

            ;
        }
        else
        {
            S = (1LL/r)*((-(1e9*cl_float(291948240, precision)+981196800LL) + (1e9*cl_float(291948240, precision)+981196800LL)*exp(2LL*rxi) -

                          (1e9*cl_float(537882537, precision)+342262612LL)*rxi - (1e9*cl_float(491868592, precision)+722131625LL)*Power(rxi, 2LL) -

                          (1e9*cl_float(297482558, precision)+476603500LL)*Power(rxi, 3LL) - (1e9*cl_float(133772487, precision)+311162700LL)*Power(rxi, 4LL) -

                          476688322649038500LL*Power(rxi, 5LL) - 140080945989184200LL*Power(rxi, 6LL) -

                          34878402537778800LL*Power(rxi, 7LL) - 7501749557702400LL*Power(rxi, 8LL) -

                          1413711970070400LL*Power(rxi, 9LL) - 235878458175360LL*Power(rxi, 10LL) -

                          35103763618560LL*Power(rxi, 11LL) - 4680908144640LL*Power(rxi, 12LL) -

                          560108666880LL*Power(rxi, 13LL) - 60011642880LL*Power(rxi, 14LL) -

                          5715394560LL*Power(rxi, 15LL) - 476282880LL*Power(rxi, 16LL) -

                          33619968LL*Power(rxi, 17LL) - 1867776LL*Power(rxi, 18LL) - 65536LL*Power(rxi, 19LL))/

                         (2.919482409811968e18*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(2LL*Power(xi, 18LL) + 38LL*Power(xi, 17LL)*xj + 342LL*Power(xi, 16LL)*Power(xj, 2LL) +

                        1938LL*Power(xi, 15LL)*Power(xj, 3LL) + 7752LL*Power(xi, 14LL)*Power(xj, 4LL) +

                        23256LL*Power(xi, 13LL)*Power(xj, 5LL) + 54264LL*Power(xi, 12LL)*Power(xj, 6LL) +

                        100776LL*Power(xi, 11LL)*Power(xj, 7LL) + 151164LL*Power(xi, 10LL)*Power(xj, 8LL) +

                        184756LL*Power(xi, 9LL)*Power(xj, 9LL) + 184756LL*Power(xi, 8LL)*Power(xj, 10LL) +

                        151164LL*Power(xi, 7LL)*Power(xj, 11LL) + 81396LL*Power(xi, 6LL)*Power(xj, 12LL) +

                        34884LL*Power(xi, 5LL)*Power(xj, 13LL) + 11628LL*Power(xi, 4LL)*Power(xj, 14LL) +

                        2907LL*Power(xi, 3LL)*Power(xj, 15LL) + 513LL*Power(xi, 2LL)*Power(xj, 16LL) +

                        57LL*xi*Power(xj, 17LL) + 3LL*Power(xj, 18LL)))/(12LL*Power(xi + xj, 19LL))

            ;
        }
        else
        {
            S = (1LL/r)*((1871100LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 19LL) +

                          495LL*exp(2LL*rxj)*Power(rxj, 14LL)*

                          (-672LL*Power(rxi, 30LL) - 12LL*Power(rxi, 31LL) + 3780LL*Power(rxj, 24LL) +

                           6615LL*rxi*Power(rxj, 24LL) - 136LL*Power(rxi, 29LL)*(126LL + Power(rxj, 2LL)) +

                           1890LL*Power(rxi, 2LL)*Power(rxj, 22LL)*(-38LL + 3LL*Power(rxj, 2LL)) +

                           315LL*Power(rxi, 3LL)*Power(rxj, 22LL)*(-399LL + 10LL*Power(rxj, 2LL)) -

                           84LL*Power(rxi, 28LL)*(3060LL + 121LL*Power(rxj, 2LL)) +

                           630LL*Power(rxi, 4LL)*Power(rxj, 20LL)*

                           (1026LL - 171LL*Power(rxj, 2LL) + 2LL*Power(rxj, 4LL)) +

                           63LL*Power(rxi, 5LL)*Power(rxj, 20LL)*

                           (17955LL - 950LL*Power(rxj, 2LL) + 6LL*Power(rxj, 4LL)) +

                           84LL*Power(rxi, 26LL)*(-174420LL - 71535LL*Power(rxj, 2LL) + 179LL*Power(rxj, 4LL)) +

                           Power(rxi, 27LL)*(-2441880LL - 327978LL*Power(rxj, 2LL) + 496LL*Power(rxj, 4LL)) -

                           63LL*Power(rxi, 19LL)*Power(rxj, 6LL)*

                           (468377895LL - 14898090LL*Power(rxj, 2LL) + 78812LL*Power(rxj, 4LL)) +

                           2LL*Power(rxi, 11LL)*Power(rxj, 14LL)*

                           (613624095LL + 56366730LL*Power(rxj, 2LL) + 383607LL*Power(rxj, 4LL) -

           248LL*Power(rxj, 6LL)) + 42LL*Power(rxi, 6LL)*Power(rxj, 18LL)*

                           (-87210LL + 23085LL*Power(rxj, 2LL) - 570LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL)) -

                           798LL*Power(rxi, 8LL)*Power(rxj, 16LL)*

                           (-18360LL + 6885LL*Power(rxj, 2LL) - 270LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL)) +

                           3LL*Power(rxi, 7LL)*Power(rxj, 18LL)*

                           (-2136645LL + 179550LL*Power(rxj, 2LL) - 2394LL*Power(rxj, 4LL) +

           4LL*Power(rxj, 6LL)) + 1596LL*Power(rxi, 14LL)*Power(rxj, 10LL)*

                           (-34484670LL - 2408985LL*Power(rxj, 2LL) + 32810LL*Power(rxj, 4LL) +

           22LL*Power(rxj, 6LL)) - 7980LL*Power(rxi, 16LL)*Power(rxj, 8LL)*

                           (15696909LL - 494343LL*Power(rxj, 2LL) - 4182LL*Power(rxj, 4LL) +

           34LL*Power(rxj, 6LL)) + 2LL*Power(rxi, 9LL)*Power(rxj, 16LL)*

                           (19433295LL - 690795LL*Power(rxj, 2LL) + 55251LL*Power(rxj, 4LL) +

           68LL*Power(rxj, 6LL)) + 6LL*Power(rxi, 25LL)*

                           (-8546580LL - 11329605LL*Power(rxj, 2LL) - 24003LL*Power(rxj, 4LL) +

           92LL*Power(rxj, 6LL)) - 6LL*Power(rxi, 13LL)*Power(rxj, 12LL)*

                           (-2361196215LL - 54738810LL*Power(rxj, 2LL) + 388626LL*Power(rxj, 4LL) +

           92LL*Power(rxj, 6LL)) + 38LL*Power(rxi, 15LL)*Power(rxj, 10LL)*

                           (808181955LL - 17168130LL*Power(rxj, 2LL) - 32130LL*Power(rxj, 4LL) +

           106LL*Power(rxj, 6LL)) - 84LL*Power(rxi, 10LL)*Power(rxj, 14LL)*

                           (3168630LL + 683145LL*Power(rxj, 2LL) + 54315LL*Power(rxj, 4LL) +

           193LL*Power(rxj, 6LL)) - 19LL*Power(rxi, 17LL)*Power(rxj, 8LL)*

                           (-2525985LL + 33479460LL*Power(rxj, 2LL) - 406980LL*Power(rxj, 4LL) +

           272LL*Power(rxj, 6LL)) + 84LL*Power(rxi, 12LL)*Power(rxj, 12LL)*

                           (-88925130LL - 19869345LL*Power(rxj, 2LL) - 235790LL*Power(rxj, 4LL) +

           643LL*Power(rxj, 6LL)) + 210LL*Power(rxi, 18LL)*Power(rxj, 6LL)*

                           (-496605582LL + 32638599LL*Power(rxj, 2LL) - 564604LL*Power(rxj, 4LL) +

           1292LL*Power(rxj, 6LL)) +

                           42LL*Power(rxi, 20LL)*Power(rxj, 4LL)*

                           (-777723210LL - 46394505LL*Power(rxj, 2LL) + 625670LL*Power(rxj, 4LL) +

           1292LL*Power(rxj, 6LL)) +

                           42LL*Power(rxi, 24LL)*(-1918620LL - 11344995LL*Power(rxj, 2LL) -

                                                  323070LL*Power(rxj, 4LL) + 2114LL*Power(rxj, 6LL)) -

                           Power(rxi, 23LL)*Power(rxj, 2LL)*

                           (1919335635LL + 275096430LL*Power(rxj, 2LL) - 3302586LL*Power(rxj, 4LL) +

           4028LL*Power(rxj, 6LL)) +

                           Power(rxi, 21LL)*Power(rxj, 4LL)*

                           (-14708379735LL + 255168270LL*Power(rxj, 2LL) - 2899134LL*Power(rxj, 4LL) +

           5168LL*Power(rxj, 6LL)) -

                           42LL*Power(rxi, 22LL)*Power(rxj, 2LL)*

                           (81654210LL + 66273255LL*Power(rxj, 2LL) - 1203870LL*Power(rxj, 4LL) +

           5206LL*Power(rxj, 6LL))) -

                          2LL*exp(2LL*rxi)*Power(rxi, 10LL)*

                          (21318LL*Power(rxi, 14LL)*Power(rxj, 14LL)*

                           (-3146850LL + 4890375LL*rxj - 24522750LL*Power(rxj, 2LL) +

           12162150LL*Power(rxj, 3LL) - 1549800LL*Power(rxj, 4LL) -

           1615950LL*Power(rxj, 5LL) - 185220LL*Power(rxj, 6LL) + 12240LL*Power(rxj, 7LL) +

           3960LL*Power(rxj, 8LL) + 300LL*Power(rxj, 9LL) + 8LL*Power(rxj, 10LL)) +

                           3LL*Power(rxi, 24LL)*Power(rxj, 4LL)*

                           (53326350LL + 97764975LL*rxj + 88877250LL*Power(rxj, 2LL) +

           53326350LL*Power(rxj, 3LL) + 23700600LL*Power(rxj, 4LL) +

           8295210LL*Power(rxj, 5LL) + 2370060LL*Power(rxj, 6LL) +

           564300LL*Power(rxj, 7LL) + 112860LL*Power(rxj, 8LL) + 22440LL*Power(rxj, 9LL) +

           1056LL*Power(rxj, 10LL) - 20LL*Power(rxj, 11LL)) -

                           4LL*Power(rxj, 28LL)*(13749310575LL + 13749310575LL*rxj +

                                                 6547290750LL*Power(rxj, 2LL) + 1964187225LL*Power(rxj, 3LL) +

                                                 413513100LL*Power(rxj, 4LL) + 64324260LL*Power(rxj, 5LL) +

                                                 7567560LL*Power(rxj, 6LL) + 675675LL*Power(rxj, 7LL) + 45045LL*Power(rxj, 8LL) +

                                                 2145LL*Power(rxj, 9LL) + 66LL*Power(rxj, 10LL) + Power(rxj, 11LL)) -

                           1254LL*Power(rxi, 16LL)*Power(rxj, 12LL)*

                           (-20241900LL - 38315025LL*rxj - 21687750LL*Power(rxj, 2LL) -

           50122800LL*Power(rxj, 3LL) + 14137200LL*Power(rxj, 4LL) -

           5853330LL*Power(rxj, 5LL) - 2687580LL*Power(rxj, 6LL) -

           208530LL*Power(rxj, 7LL) + 19530LL*Power(rxj, 8LL) + 3630LL*Power(rxj, 9LL) +

           172LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) +

                           627LL*Power(rxi, 12LL)*Power(rxj, 16LL)*

                           (-1240964550LL + 4740389325LL*rxj - 3311818650LL*Power(rxj, 2LL) +

           134804250LL*Power(rxj, 3LL) + 407673000LL*Power(rxj, 4LL) +

           58641030LL*Power(rxj, 5LL) - 3549420LL*Power(rxj, 6LL) -

           1641060LL*Power(rxj, 7LL) - 167940LL*Power(rxj, 8LL) - 6990LL*Power(rxj, 9LL) -

           36LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           Power(rxi, 28LL)*(935550LL + 1715175LL*rxj + 1559250LL*Power(rxj, 2LL) +

                                             935550LL*Power(rxj, 3LL) + 415800LL*Power(rxj, 4LL) + 145530LL*Power(rxj, 5LL) +

                                             41580LL*Power(rxj, 6LL) + 9900LL*Power(rxj, 7LL) + 1980LL*Power(rxj, 8LL) +

                                             330LL*Power(rxj, 9LL) + 44LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) +

                           2LL*Power(rxi, 2LL)*Power(rxj, 26LL)*

                           (-937068397650LL - 815439881025LL*rxj - 332904552750LL*Power(rxj, 2LL) -

           84006776700LL*Power(rxj, 3LL) - 14504767200LL*Power(rxj, 4LL) -

           1786235220LL*Power(rxj, 5LL) - 157754520LL*Power(rxj, 6LL) -

           9667350LL*Power(rxj, 7LL) - 367290LL*Power(rxj, 8LL) - 5115LL*Power(rxj, 9LL) +

           198LL*Power(rxj, 10LL) + 8LL*Power(rxj, 11LL)) +

                           6LL*Power(rxi, 4LL)*Power(rxj, 24LL)*

                           (-2262441500550LL - 1503711230175LL*rxj - 426178264050LL*Power(rxj, 2LL) -

           60134347350LL*Power(rxj, 3LL) - 2014551000LL*Power(rxj, 4LL) +

           846111420LL*Power(rxj, 5LL) + 184864680LL*Power(rxj, 6LL) +

           20183130LL*Power(rxj, 7LL) + 1367190LL*Power(rxj, 8LL) +

           57255LL*Power(rxj, 9LL) + 1298LL*Power(rxj, 10LL) + 10LL*Power(rxj, 11LL)) -

                           Power(rxi, 26LL)*Power(rxj, 2LL)*

                           (17775450LL + 32588325LL*rxj + 29625750LL*Power(rxj, 2LL) +

           17775450LL*Power(rxj, 3LL) + 7900200LL*Power(rxj, 4LL) +

           2765070LL*Power(rxj, 5LL) + 790020LL*Power(rxj, 6LL) +

           188100LL*Power(rxj, 7LL) + 37620LL*Power(rxj, 8LL) + 6270LL*Power(rxj, 9LL) +

           836LL*Power(rxj, 10LL) + 16LL*Power(rxj, 11LL)) +

                           9LL*Power(rxi, 22LL)*Power(rxj, 6LL)*

                           (-100727550LL - 184667175LL*rxj - 167879250LL*Power(rxj, 2LL) -

           100727550LL*Power(rxj, 3LL) - 44767800LL*Power(rxj, 4LL) -

           15668730LL*Power(rxj, 5LL) - 4476780LL*Power(rxj, 6LL) -

           971520LL*Power(rxj, 7LL) - 307560LL*Power(rxj, 8LL) - 27060LL*Power(rxj, 9LL) +

           264LL*Power(rxj, 10LL) + 64LL*Power(rxj, 11LL)) -

                           9LL*Power(rxi, 6LL)*Power(rxj, 22LL)*

                           (3452543428950LL + 1097992509075LL*rxj - 101420792550LL*Power(rxj, 2LL) -

           110557373850LL*Power(rxj, 3LL) - 24909330600LL*Power(rxj, 4LL) -

           2686726350LL*Power(rxj, 5LL) - 93485700LL*Power(rxj, 6LL) +

           12941280LL*Power(rxj, 7LL) + 2081640LL*Power(rxj, 8LL) +

           137940LL*Power(rxj, 9LL) + 4664LL*Power(rxj, 10LL) + 64LL*Power(rxj, 11LL)) -

                           22LL*Power(rxi, 20LL)*Power(rxj, 8LL)*

                           (-164826900LL - 302182650LL*rxj - 274711500LL*Power(rxj, 2LL) -

           164826900LL*Power(rxj, 3LL) - 73256400LL*Power(rxj, 4LL) -

           26991090LL*Power(rxj, 5LL) - 4622940LL*Power(rxj, 6LL) -

           2941110LL*Power(rxj, 7LL) - 438930LL*Power(rxj, 8LL) - 5505LL*Power(rxj, 9LL) +

           2082LL*Power(rxj, 10LL) + 82LL*Power(rxj, 11LL)) +

                           22LL*Power(rxi, 18LL)*Power(rxj, 10LL)*

                           (-494480700LL - 906547950LL*rxj - 824134500LL*Power(rxj, 2LL) -

           475684650LL*Power(rxj, 3LL) - 294953400LL*Power(rxj, 4LL) +

           2663010LL*Power(rxj, 5LL) - 40797540LL*Power(rxj, 6LL) -

           10248390LL*Power(rxj, 7LL) - 434610LL*Power(rxj, 8LL) +

           65865LL*Power(rxj, 9LL) + 6366LL*Power(rxj, 10LL) + 136LL*Power(rxj, 11LL)) +

                           11LL*Power(rxi, 8LL)*Power(rxj, 20LL)*

                           (-2338604626050LL + 656001834075LL*rxj + 504510561450LL*Power(rxj, 2LL) +

           51560967150LL*Power(rxj, 3LL) - 15574998600LL*Power(rxj, 4LL) -

           5055778350LL*Power(rxj, 5LL) - 626213700LL*Power(rxj, 6LL) -

           34768620LL*Power(rxj, 7LL) + 207540LL*Power(rxj, 8LL) +

           150240LL*Power(rxj, 9LL) + 8464LL*Power(rxj, 10LL) + 164LL*Power(rxj, 11LL)) -

                           11LL*Power(rxi, 10LL)*Power(rxj, 18LL)*

                           (742805182350LL - 933111659025LL*rxj + 57080542050LL*Power(rxj, 2LL) +

           129505209750LL*Power(rxj, 3LL) + 19066887000LL*Power(rxj, 4LL) -

           1817573310LL*Power(rxj, 5LL) - 810647460LL*Power(rxj, 6LL) -

           97669980LL*Power(rxj, 7LL) - 5173020LL*Power(rxj, 8LL) -

           37770LL*Power(rxj, 9LL) + 8212LL*Power(rxj, 10LL) + 272LL*Power(rxj, 11LL))))/

                         (1.8711e6*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 19LL)*Power(rxi + rxj, 19LL))

                         );
        }

    }
    return S;
}


cl_R Slater_6S_4S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_4S_6S(r, xj, xi);
}
