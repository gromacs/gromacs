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

cl_R Slater_5S_6S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (235451LL*xi)/1.572864e6

            ;
        }
        else
        {
            S = (1LL/r)*((-(1e9*cl_float(153272826, precision)+515128320000LL) + (1e9*cl_float(153272826, precision)+515128320000LL)*exp(2LL*rxi) -

                          (1e9*cl_float(283601367, precision)+780029361562LL)*rxi - (1e9*cl_float(260657082, precision)+529802083125LL)*Power(rxi, 2LL) -

                          (1e9*cl_float(158611593, precision)+237817407187LL)*Power(rxi, 3LL) -

                          (1e9*cl_float(718622941, precision)+12650916875LL)*Power(rxi, 4LL) -

                          (1e9*cl_float(258482050, precision)+83510960150LL)*Power(rxi, 5LL) - (1e9*cl_float(768533806, precision)+7827219800LL)*Power(rxi, 6LL) -

                          (1e9*cl_float(194179852, precision)+3340075400LL)*Power(rxi, 7LL) - (1e9*cl_float(425318313, precision)+470450400LL)*Power(rxi, 8LL) -

                          819670099680432000LL*Power(rxi, 9LL) - 140553592289510400LL*Power(rxi, 10LL) -

                          21625475644281600LL*Power(rxi, 11LL) - 3003582726144000LL*Power(rxi, 12LL) -

                          378073350144000LL*Power(rxi, 13LL) - 43208382873600LL*Power(rxi, 14LL) -

                          4480869335040LL*Power(rxi, 15LL) - 420081500160LL*Power(rxi, 16LL) -

                          35300966400LL*Power(rxi, 17LL) - 2614886400LL*Power(rxi, 18LL) -

                          165150720LL*Power(rxi, 19LL) - 8257536LL*Power(rxi, 20LL) - 262144LL*Power(rxi, 21LL))/

                         (1.5327282651512832e21*exp(2LL*rxi))

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(5LL*Power(xi, 20LL) + 105LL*Power(xi, 19LL)*xj +

                        1050LL*Power(xi, 18LL)*Power(xj, 2LL) + 6650LL*Power(xi, 17LL)*Power(xj, 3LL) +

                        29925LL*Power(xi, 16LL)*Power(xj, 4LL) + 101745LL*Power(xi, 15LL)*Power(xj, 5LL) +

                        271320LL*Power(xi, 14LL)*Power(xj, 6LL) + 581400LL*Power(xi, 13LL)*Power(xj, 7LL) +

                        1017450LL*Power(xi, 12LL)*Power(xj, 8LL) + 1469650LL*Power(xi, 11LL)*Power(xj, 9LL) +

                        1763580LL*Power(xi, 10LL)*Power(xj, 10LL) + 1763580LL*Power(xi, 9LL)*Power(xj, 11LL) +

                        1220940LL*Power(xi, 8LL)*Power(xj, 12LL) + 697680LL*Power(xi, 7LL)*Power(xj, 13LL) +

                        325584LL*Power(xi, 6LL)*Power(xj, 14LL) + 122094LL*Power(xi, 5LL)*Power(xj, 15LL) +

                        35910LL*Power(xi, 4LL)*Power(xj, 16LL) + 7980LL*Power(xi, 3LL)*Power(xj, 17LL) +

                        1260LL*Power(xi, 2LL)*Power(xj, 18LL) + 126LL*xi*Power(xj, 19LL) + 6LL*Power(xj, 20LL)))/

                (30LL*Power(xi + xj, 21LL))

            ;
        }
        else
        {
            S = (1LL/r)*((4677750LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 21LL) +

                          110LL*exp(2LL*rxj)*Power(rxj, 14LL)*

                          (-432LL*Power(rxi, 36LL) - 6LL*Power(rxi, 37LL) + 42525LL*Power(rxj, 28LL) +

                           76545LL*rxi*Power(rxj, 28LL) +

                           19845LL*Power(rxi, 3LL)*Power(rxj, 26LL)*(-81LL + 2LL*Power(rxj, 2LL)) -

                           1134LL*Power(rxi, 34LL)*(272LL + 5LL*Power(rxj, 2LL)) -

                           8LL*Power(rxi, 35LL)*(1836LL + 7LL*Power(rxj, 2LL)) +

                           8505LL*Power(rxi, 2LL)*Power(rxj, 26LL)*(-105LL + 8LL*Power(rxj, 2LL)) +

                           378LL*Power(rxi, 33LL)*(-11628LL - 666LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           5670LL*Power(rxi, 5LL)*Power(rxj, 24LL)*

                           (2835LL - 147LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           17010LL*Power(rxi, 4LL)*Power(rxj, 24LL)*

                           (525LL - 84LL*Power(rxj, 2LL) + Power(rxj, 4LL)) +

                           378LL*Power(rxi, 32LL)*(-116280LL - 17444LL*Power(rxj, 2LL) + 59LL*Power(rxj, 4LL)) +

                           162LL*Power(rxi, 7LL)*Power(rxj, 22LL)*

                           (-628425LL + 51450LL*Power(rxj, 2LL) - 735LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL))

                           + 378LL*Power(rxi, 6LL)*Power(rxj, 22LL)*

                           (-149625LL + 37800LL*Power(rxj, 2LL) - 945LL*Power(rxj, 4LL) + 4LL*Power(rxj, 6LL))

                           - 18LL*Power(rxi, 31LL)*(17093160LL + 6309387LL*Power(rxj, 2LL) - 23562LL*Power(rxj, 4LL) +

                                                    16LL*Power(rxj, 6LL)) + 54LL*Power(rxi, 30LL)*

                           (-26860680LL - 24843735LL*Power(rxj, 2LL) - 40180LL*Power(rxj, 4LL) +

           578LL*Power(rxj, 6LL)) + 378LL*Power(rxi, 23LL)*Power(rxj, 6LL)*

                           (-14625683325LL + 704051250LL*Power(rxj, 2LL) - 10752861LL*Power(rxj, 4LL) +

           33478LL*Power(rxj, 6LL)) +

                           3LL*Power(rxi, 9LL)*Power(rxj, 20LL)*

                           (152707275LL - 17595900LL*Power(rxj, 2LL) + 396900LL*Power(rxj, 4LL) -

           2268LL*Power(rxj, 6LL) + 2LL*Power(rxj, 8LL)) +

                           27LL*Power(rxi, 8LL)*Power(rxj, 20LL)*

                           (9426375LL - 3351600LL*Power(rxj, 2LL) + 132300LL*Power(rxj, 4LL) -

           1176LL*Power(rxj, 6LL) + 2LL*Power(rxj, 8LL)) -

                           567LL*Power(rxi, 10LL)*Power(rxj, 18LL)*

                           (1526175LL - 718200LL*Power(rxj, 2LL) + 39900LL*Power(rxj, 4LL) -

           560LL*Power(rxj, 6LL) + 2LL*Power(rxj, 8LL)) -

                           54LL*Power(rxi, 13LL)*Power(rxj, 16LL)*

                           (-1356769575LL - 127011675LL*Power(rxj, 2LL) - 3867843LL*Power(rxj, 4LL) -

           8556LL*Power(rxj, 6LL) + 7LL*Power(rxj, 8LL)) +

                           7LL*Power(rxi, 11LL)*Power(rxj, 18LL)*

                           (-151091325LL + 45272250LL*Power(rxj, 2LL) - 647676LL*Power(rxj, 4LL) +

           15336LL*Power(rxj, 6LL) + 8LL*Power(rxj, 8LL)) +

                           18LL*Power(rxi, 15LL)*Power(rxj, 14LL)*

                           (63046289250LL + 3917182500LL*Power(rxj, 2LL) + 10158435LL*Power(rxj, 4LL) -

           178842LL*Power(rxj, 6LL) + 16LL*Power(rxj, 8LL)) +

                           378LL*Power(rxi, 21LL)*Power(rxj, 8LL)*

                           (-8559820125LL + 17573325LL*Power(rxj, 2LL) + 7421001LL*Power(rxj, 4LL) -

           49096LL*Power(rxj, 6LL) + 19LL*Power(rxj, 8LL)) -

                           378LL*Power(rxi, 12LL)*Power(rxj, 16LL)*

                           (17296650LL + 14244300LL*Power(rxj, 2LL) + 360525LL*Power(rxj, 4LL) +

           15928LL*Power(rxj, 6LL) + 22LL*Power(rxj, 8LL)) -

                           189LL*Power(rxi, 25LL)*Power(rxj, 4LL)*

                           (9994948425LL + 63821700LL*Power(rxj, 2LL) - 1458540LL*Power(rxj, 4LL) -

           18756LL*Power(rxj, 6LL) + 38LL*Power(rxj, 8LL)) -

                           189LL*Power(rxi, 24LL)*Power(rxj, 4LL)*

                           (17962854525LL + 4036942800LL*Power(rxj, 2LL) - 126472500LL*Power(rxj, 4LL) +

           765464LL*Power(rxj, 6LL) + 190LL*Power(rxj, 8LL)) -

                           21LL*Power(rxi, 19LL)*Power(rxj, 10LL)*

                           (-228066210225LL + 13487616450LL*Power(rxj, 2LL) -

           85465800LL*Power(rxj, 4LL) - 320112LL*Power(rxj, 6LL) + 328LL*Power(rxj, 8LL)) -

                           189LL*Power(rxi, 18LL)*Power(rxj, 10LL)*

                           (86069971575LL + 2157712200LL*Power(rxj, 2LL) - 158179560LL*Power(rxj, 4LL) +

           578816LL*Power(rxj, 6LL) + 978LL*Power(rxj, 8LL)) -

                           2LL*Power(rxi, 29LL)*(2085060285LL + 5450330025LL*Power(rxj, 2LL) +

                                                 127424745LL*Power(rxj, 4LL) - 1398276LL*Power(rxj, 6LL) + 1159LL*Power(rxj, 8LL)

                                                 ) - 378LL*Power(rxi, 22LL)*Power(rxj, 6LL)*

                           (37244490525LL - 2411839800LL*Power(rxj, 2LL) + 92951775LL*Power(rxj, 4LL) -

           942172LL*Power(rxj, 6LL) + 1292LL*Power(rxj, 8LL)) -

                           27LL*Power(rxi, 16LL)*Power(rxj, 12LL)*

                           (164245367475LL + 26909517600LL*Power(rxj, 2LL) + 62674920LL*Power(rxj, 4LL) -

           3885112LL*Power(rxj, 6LL) + 2122LL*Power(rxj, 8LL)) +

                           3LL*Power(rxi, 27LL)*Power(rxj, 2LL)*

                           (-63819198135LL - 21841975890LL*Power(rxj, 2LL) +

           442430100LL*Power(rxj, 4LL) - 2756664LL*Power(rxj, 6LL) + 2296LL*Power(rxj, 8LL)

                           ) + Power(rxi, 17LL)*Power(rxj, 12LL)*

                           (4851990871875LL + 21622847400LL*Power(rxj, 2LL) -

           2153738160LL*Power(rxj, 4LL) + 3608388LL*Power(rxj, 6LL) +

           2318LL*Power(rxj, 8LL)) +

                           18LL*Power(rxi, 14LL)*Power(rxj, 14LL)*

                           (-23418646650LL - 6922729800LL*Power(rxj, 2LL) - 259958475LL*Power(rxj, 4LL) -

           697732LL*Power(rxj, 6LL) + 3030LL*Power(rxj, 8LL)) +

                           126LL*Power(rxi, 20LL)*Power(rxj, 8LL)*

                           (-186637212225LL + 13028280300LL*Power(rxj, 2LL) -

           116198775LL*Power(rxj, 4LL) - 1266160LL*Power(rxj, 6LL) + 4332LL*Power(rxj, 8LL)

                           ) - 54LL*Power(rxi, 28LL)*(102965940LL + 1089437580LL*Power(rxj, 2LL) +

                                                      102508245LL*Power(rxj, 4LL) - 1593144LL*Power(rxj, 6LL) + 4538LL*Power(rxj, 8LL)

                                                      ) + 63LL*Power(rxi, 26LL)*Power(rxj, 2LL)*

                           (-4544129205LL - 7396000920LL*Power(rxj, 2LL) + 149614020LL*Power(rxj, 4LL) -

           1684112LL*Power(rxj, 6LL) + 5922LL*Power(rxj, 8LL))) +

                          exp(2LL*rxi)*Power(rxi, 12LL)*

                          (6LL*Power(rxi, 24LL)*Power(rxj, 6LL)*

                           (1036901250LL + 1900985625LL*rxj + 1728168750LL*Power(rxj, 2LL) +

           1036901250LL*Power(rxj, 3LL) + 460845000LL*Power(rxj, 4LL) +

           161295750LL*Power(rxj, 5LL) + 46084500LL*Power(rxj, 6LL) +

           9084900LL*Power(rxj, 7LL) + 4082100LL*Power(rxj, 8LL) +

           121935LL*Power(rxj, 9LL) - 21494LL*Power(rxj, 10LL) - 766LL*Power(rxj, 11LL)) +

                           5LL*Power(rxi, 28LL)*Power(rxj, 2LL)*

                           (19646550LL + 36018675LL*rxj + 32744250LL*Power(rxj, 2LL) +

           19646550LL*Power(rxj, 3LL) + 8731800LL*Power(rxj, 4LL) +

           3056130LL*Power(rxj, 5LL) + 873180LL*Power(rxj, 6LL) +

           207900LL*Power(rxj, 7LL) + 41580LL*Power(rxj, 8LL) + 6930LL*Power(rxj, 9LL) +

           924LL*Power(rxj, 10LL) - 4LL*Power(rxj, 11LL)) +

                           26334LL*Power(rxi, 16LL)*Power(rxj, 14LL)*

                           (43880400LL - 186686775LL*rxj + 576771750LL*Power(rxj, 2LL) -

           398603250LL*Power(rxj, 3LL) + 72552600LL*Power(rxj, 4LL) +

           27903120LL*Power(rxj, 5LL) - 342720LL*Power(rxj, 6LL) -

           574800LL*Power(rxj, 7LL) - 50800LL*Power(rxj, 8LL) - 945LL*Power(rxj, 9LL) +

           58LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) +

                           10LL*Power(rxj, 30LL)*(97302813300LL + 89194245525LL*rxj +

                                                  38780106750LL*Power(rxj, 2LL) + 10576392750LL*Power(rxj, 3LL) +

                                                  2014551000LL*Power(rxj, 4LL) + 282037140LL*Power(rxj, 5LL) +

                                                  29688120LL*Power(rxj, 6LL) + 2356200LL*Power(rxj, 7LL) +

                                                  138600LL*Power(rxj, 8LL) + 5775LL*Power(rxj, 9LL) + 154LL*Power(rxj, 10LL) +

                                                  2LL*Power(rxj, 11LL)) + 10LL*Power(rxi, 2LL)*Power(rxj, 28LL)*

                           (4582499159700LL + 3733416276975LL*rxj + 1428215931450LL*Power(rxj, 2LL) +

           338545295550LL*Power(rxj, 3LL) + 55198697400LL*Power(rxj, 4LL) +

           6486854220LL*Power(rxj, 5LL) + 558419400LL*Power(rxj, 6LL) +

           34939080LL*Power(rxj, 7LL) + 1532520LL*Power(rxj, 8LL) +

           43285LL*Power(rxj, 9LL) + 638LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) -

                           110LL*Power(rxi, 10LL)*Power(rxj, 20LL)*

                           (-14063418170550LL + 6795156458475LL*rxj + 2067471236250LL*Power(rxj, 2LL) -

           214664924250LL*Power(rxj, 3LL) - 124416469800LL*Power(rxj, 4LL) -

           14935545450LL*Power(rxj, 5LL) - 256688460LL*Power(rxj, 6LL) +

           105750900LL*Power(rxj, 7LL) + 11502180LL*Power(rxj, 8LL) +

           518085LL*Power(rxj, 9LL) + 9294LL*Power(rxj, 10LL) + 2LL*Power(rxj, 11LL)) +

                           55LL*Power(rxi, 20LL)*Power(rxj, 10LL)*

                           (1730682450LL + 3172917825LL*rxj + 2884470750LL*Power(rxj, 2LL) +

           1571960250LL*Power(rxj, 3LL) + 1404081000LL*Power(rxj, 4LL) -

           426654270LL*Power(rxj, 5LL) + 283536540LL*Power(rxj, 6LL) +

           39116700LL*Power(rxj, 7LL) - 2659860LL*Power(rxj, 8LL) -

           528850LL*Power(rxj, 9LL) - 18236LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           5LL*Power(rxi, 30LL)*(935550LL + 1715175LL*rxj + 1559250LL*Power(rxj, 2LL) +

                                                 935550LL*Power(rxj, 3LL) + 415800LL*Power(rxj, 4LL) + 145530LL*Power(rxj, 5LL) +

                                                 41580LL*Power(rxj, 6LL) + 9900LL*Power(rxj, 7LL) + 1980LL*Power(rxj, 8LL) +

                                                 330LL*Power(rxj, 9LL) + 44LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           13167LL*Power(rxi, 14LL)*Power(rxj, 16LL)*

                           (-2319354450LL + 8540029575LL*rxj - 7335672750LL*Power(rxj, 2LL) +

           1133154750LL*Power(rxj, 3LL) + 575014200LL*Power(rxj, 4LL) -

           913710LL*Power(rxj, 5LL) - 14863940LL*Power(rxj, 6LL) -

           1687300LL*Power(rxj, 7LL) - 46900LL*Power(rxj, 8LL) + 3210LL*Power(rxj, 9LL) +

           236LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           770LL*Power(rxi, 18LL)*Power(rxj, 12LL)*

                           (329653800LL + 654729075LL*rxj + 45785250LL*Power(rxj, 2LL) +

           1602483750LL*Power(rxj, 3LL) - 915705000LL*Power(rxj, 4LL) +

           266036400LL*Power(rxj, 5LL) + 63745920LL*Power(rxj, 6LL) -

           2304000LL*Power(rxj, 7LL) - 1074240LL*Power(rxj, 8LL) -

           64635LL*Power(rxj, 9LL) - 514LL*Power(rxj, 10LL) + 34LL*Power(rxj, 11LL)) +

                           385LL*Power(rxi, 12LL)*Power(rxj, 18LL)*

                           (973565393850LL - 1429122323475LL*rxj + 298281831750LL*Power(rxj, 2LL) +

           138841148250LL*Power(rxj, 3LL) - 2454240600LL*Power(rxj, 4LL) -

           4925394810LL*Power(rxj, 5LL) - 623832300LL*Power(rxj, 6LL) -

           19098540LL*Power(rxj, 7LL) + 2083140LL*Power(rxj, 8LL) +

           212430LL*Power(rxj, 9LL) + 7012LL*Power(rxj, 10LL) + 68LL*Power(rxj, 11LL)) +

                           14LL*Power(rxi, 26LL)*Power(rxj, 4LL)*

                           (-70166250LL - 128638125LL*rxj - 116943750LL*Power(rxj, 2LL) -

           70166250LL*Power(rxj, 3LL) - 31185000LL*Power(rxj, 4LL) -

           10914750LL*Power(rxj, 5LL) - 3118500LL*Power(rxj, 6LL) -

           742500LL*Power(rxj, 7LL) - 148500LL*Power(rxj, 8LL) - 32615LL*Power(rxj, 9LL) -

           154LL*Power(rxj, 10LL) + 74LL*Power(rxj, 11LL)) -

                           7LL*Power(rxi, 4LL)*Power(rxj, 26LL)*

                           (-69822945249750LL - 46669577290875LL*rxj -

           14025037430250LL*Power(rxj, 2LL) - 2430881664750LL*Power(rxj, 3LL) -

           251629270200LL*Power(rxj, 4LL) - 12434519790LL*Power(rxj, 5LL) +

           452930940LL*Power(rxj, 6LL) + 131125500LL*Power(rxj, 7LL) +

           11018700LL*Power(rxj, 8LL) + 514470LL*Power(rxj, 9LL) +

           13332LL*Power(rxj, 10LL) + 148LL*Power(rxj, 11LL)) -

                           50LL*Power(rxi, 8LL)*Power(rxj, 22LL)*

                           (-51768833574150LL - 5003280391725LL*rxj + 4493439477450LL*Power(rxj, 2LL) +

           1286866176750LL*Power(rxj, 3LL) + 111437476920LL*Power(rxj, 4LL) -

           6620313546LL*Power(rxj, 5LL) - 2406603276LL*Power(rxj, 6LL) -

           242686620LL*Power(rxj, 7LL) - 12228876LL*Power(rxj, 8LL) -

           256223LL*Power(rxj, 9LL) + 2486LL*Power(rxj, 10LL) + 158LL*Power(rxj, 11LL)) +

                           25LL*Power(rxi, 22LL)*Power(rxj, 8LL)*

                           (-1119853350LL - 2053064475LL*rxj - 1866422250LL*Power(rxj, 2LL) -

           1119853350LL*Power(rxj, 3LL) - 497712600LL*Power(rxj, 4LL) -

           194415606LL*Power(rxj, 5LL) - 9338868LL*Power(rxj, 6LL) -

           31217076LL*Power(rxj, 7LL) - 2256804LL*Power(rxj, 8LL) +

           246774LL*Power(rxj, 9LL) + 22836LL*Power(rxj, 10LL) + 316LL*Power(rxj, 11LL)) +

                           3LL*Power(rxi, 6LL)*Power(rxj, 24LL)*

                           (596006592662250LL + 266778699697125LL*rxj +

           37515651153750LL*Power(rxj, 2LL) - 2214626163750LL*Power(rxj, 3LL) -

           1538075107800LL*Power(rxj, 4LL) - 248955308910LL*Power(rxj, 5LL) -

           21434337540LL*Power(rxj, 6LL) - 957980100LL*Power(rxj, 7LL) -

           4874100LL*Power(rxj, 8LL) + 1831830LL*Power(rxj, 9LL) +

           91828LL*Power(rxj, 10LL) + 1532LL*Power(rxj, 11LL))))/

                         (4.67775e6*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 21LL)*Power(rxi + rxj, 21LL))

                         );
        }

    }
    return S;
}


cl_R Slater_6S_5S(cl_R r, cl_R xi, cl_R xj)
{
    return Slater_5S_6S(r, xj, xi);
}
