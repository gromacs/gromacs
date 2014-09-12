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

cl_R Slater_6S_6S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = (1172755LL*xi)/8.388608e6

            ;
        }
        else
        {
            S = (1LL/r)*(1LL + (-1LL - (15604461LL*rxi)/8.388608e6 - (7215853LL*Power(rxi, 2LL))/4.194304e6 -

                                (19903123LL*Power(rxi, 3LL))/1.8874368e7 -

                                (1136755LL*Power(rxi, 4LL))/2.359296e6 - (309683LL*Power(rxi, 5LL))/1.769472e6 -

                                (13797LL*Power(rxi, 6LL))/262144LL - (238423LL*Power(rxi, 7LL))/1.769472e7 -

                                (139189LL*Power(rxi, 8LL))/4.644864e7 - (163811LL*Power(rxi, 9LL))/2.7869184e8 -

                                (71677LL*Power(rxi, 10LL))/6.967296e8 -

                                (186367LL*Power(rxi, 11LL))/1.14960384e10 - (13LL*Power(rxi, 12LL))/5.6133e6 -

                                Power(rxi, 13LL)/3.31695e6 - Power(rxi, 14LL)/2.786238e7 -

                                Power(rxi, 15LL)/2.5540515e8 - Power(rxi, 16LL)/2.5540515e9 -

                                Power(rxi, 17LL)/2.791213425e10 - Power(rxi, 18LL)/3.34945611e11 -

                                Power(rxi, 19LL)/4.4547766263e12 - Power(rxi, 20LL)/6.68216493945e13 -

                                Power(rxi, 21LL)/1.16937886440375e15 - Power(rxi, 22LL)/2.57263350168825e16 -

                                Power(rxi, 23LL)/8.875585580824462e17)/exp(2LL*rxi)

                         );
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = (xi*xj*(Power(xi, 22LL) + 23LL*Power(xi, 21LL)*xj + 253LL*Power(xi, 20LL)*Power(xj, 2LL) +

                        1771LL*Power(xi, 19LL)*Power(xj, 3LL) + 8855LL*Power(xi, 18LL)*Power(xj, 4LL) +

                        33649LL*Power(xi, 17LL)*Power(xj, 5LL) + 100947LL*Power(xi, 16LL)*Power(xj, 6LL) +

                        245157LL*Power(xi, 15LL)*Power(xj, 7LL) + 490314LL*Power(xi, 14LL)*Power(xj, 8LL) +

                        817190LL*Power(xi, 13LL)*Power(xj, 9LL) + 1144066LL*Power(xi, 12LL)*Power(xj, 10LL) +

                        1352078LL*Power(xi, 11LL)*Power(xj, 11LL) +

                        1144066LL*Power(xi, 10LL)*Power(xj, 12LL) + 817190LL*Power(xi, 9LL)*Power(xj, 13LL) +

                        490314LL*Power(xi, 8LL)*Power(xj, 14LL) + 245157LL*Power(xi, 7LL)*Power(xj, 15LL) +

                        100947LL*Power(xi, 6LL)*Power(xj, 16LL) + 33649LL*Power(xi, 5LL)*Power(xj, 17LL) +

                        8855LL*Power(xi, 4LL)*Power(xj, 18LL) + 1771LL*Power(xi, 3LL)*Power(xj, 19LL) +

                        253LL*Power(xi, 2LL)*Power(xj, 20LL) + 23LL*xi*Power(xj, 21LL) + Power(xj, 22LL)))/

                (6LL*Power(xi + xj, 23LL))

            ;
        }
        else
        {
            S = (1LL/r)*((2806650LL*exp(2LL*(rxi + rxj))*Power(Power(rxi, 2LL) - Power(rxj, 2LL), 23LL) +

                          exp(2LL*rxj)*Power(rxj, 14LL)*

                          (-1056LL*Power(rxi, 42LL) - 12LL*Power(rxi, 43LL) + 2806650LL*Power(rxj, 32LL) +

                           5145525LL*rxi*Power(rxj, 32LL) - 88LL*Power(rxi, 41LL)*(510LL + Power(rxj, 2LL)) +

                           935550LL*Power(rxi, 2LL)*Power(rxj, 30LL)*(-69LL + 5LL*Power(rxj, 2LL)) +

                           467775LL*Power(rxi, 3LL)*Power(rxj, 30LL)*(-253LL + 6LL*Power(rxj, 2LL)) -

                           132LL*Power(rxi, 40LL)*(9180LL + 89LL*Power(rxj, 2LL)) +

                           311850LL*Power(rxi, 4LL)*Power(rxj, 28LL)*

                           (2277LL - 345LL*Power(rxj, 2LL) + 4LL*Power(rxj, 4LL)) +

                           31185LL*Power(rxi, 5LL)*Power(rxj, 28LL)*

                           (41745LL - 2070LL*Power(rxj, 2LL) + 14LL*Power(rxj, 4LL)) +

                           1980LL*Power(rxi, 38LL)*(-162792LL - 11859LL*Power(rxj, 2LL) +

                                                    41LL*Power(rxj, 4LL)) + 22LL*Power(rxi, 39LL)*

                           (-1046520LL - 30885LL*Power(rxj, 2LL) + 44LL*Power(rxj, 4LL)) +

                           62370LL*Power(rxi, 6LL)*Power(rxj, 26LL)*

                           (-79695LL + 18975LL*Power(rxj, 2LL) - 460LL*Power(rxj, 4LL) + 2LL*Power(rxj, 6LL)) -

                           110LL*Power(rxi, 37LL)*(30767688LL + 4989438LL*Power(rxj, 2LL) -

                                                   25359LL*Power(rxj, 4LL) + 20LL*Power(rxj, 6LL)) +

                           1485LL*Power(rxi, 7LL)*Power(rxj, 26LL)*

                           (-6136515LL + 478170LL*Power(rxj, 2LL) - 6762LL*Power(rxj, 4LL) +

           20LL*Power(rxj, 6LL)) - 132LL*Power(rxi, 36LL)*

                           (201455100LL + 69647445LL*Power(rxj, 2LL) - 318735LL*Power(rxj, 4LL) +

           353LL*Power(rxj, 6LL)) + 495LL*Power(rxi, 9LL)*Power(rxj, 24LL)*

                           (92047725LL - 10041570LL*Power(rxj, 2LL) + 223146LL*Power(rxj, 4LL) -

           1380LL*Power(rxj, 6LL) + 2LL*Power(rxj, 8LL)) +

                           2970LL*Power(rxi, 8LL)*Power(rxj, 24LL)*

                           (8367975LL - 2789325LL*Power(rxj, 2LL) + 106260LL*Power(rxj, 4LL) -

           966LL*Power(rxj, 6LL) + 2LL*Power(rxj, 8LL)) -

                           22LL*Power(rxi, 35LL)*(6950200950LL + 5142653145LL*Power(rxj, 2LL) +

                                                  7644510LL*Power(rxj, 4LL) - 235635LL*Power(rxj, 6LL) + 124LL*Power(rxj, 8LL)) -

                           132LL*Power(rxi, 34LL)*(4633467300LL + 7767871650LL*Power(rxj, 2LL) +

                                                   160904205LL*Power(rxj, 4LL) - 2493315LL*Power(rxj, 6LL) + 5281LL*Power(rxj, 8LL)

                                                   ) - 495LL*Power(rxi, 27LL)*Power(rxj, 6LL)*

                           (8395934795325LL - 439434024750LL*Power(rxj, 2LL) +

           11948496210LL*Power(rxj, 4LL) - 118623972LL*Power(rxj, 6LL) +

           248906LL*Power(rxj, 8LL)) +

                           11LL*Power(rxi, 15LL)*Power(rxj, 18LL)*

                           (1488922594425LL + 252796524750LL*Power(rxj, 2LL) +

           6172031250LL*Power(rxj, 4LL) + 104343660LL*Power(rxj, 6LL) +

           66810LL*Power(rxj, 8LL) - 88LL*Power(rxj, 10LL)) +

                           66LL*Power(rxi, 10LL)*Power(rxj, 22LL)*

                           (-1430923725LL + 627598125LL*Power(rxj, 2LL) - 33471900LL*Power(rxj, 4LL) +

           478170LL*Power(rxj, 6LL) - 2070LL*Power(rxj, 8LL) + 2LL*Power(rxj, 10LL)) -

                           1518LL*Power(rxi, 12LL)*Power(rxj, 20LL)*

                           (-186642225LL + 103690125LL*Power(rxj, 2LL) - 7276500LL*Power(rxj, 4LL) +

           145530LL*Power(rxj, 6LL) - 990LL*Power(rxj, 8LL) + 2LL*Power(rxj, 10LL)) +

                           3LL*Power(rxi, 11LL)*Power(rxj, 22LL)*

                           (-57713923575LL + 8284295250LL*Power(rxj, 2LL) - 257733630LL*Power(rxj, 4LL) +

           2504700LL*Power(rxj, 6LL) - 7590LL*Power(rxj, 8LL) + 4LL*Power(rxj, 10LL)) +

                           11LL*Power(rxi, 13LL)*Power(rxj, 20LL)*

                           (56066193225LL - 6918959250LL*Power(rxj, 2LL) + 430816050LL*Power(rxj, 4LL) -

           3349620LL*Power(rxj, 6LL) + 33690LL*Power(rxj, 8LL) + 8LL*Power(rxj, 10LL)) +

                           55LL*Power(rxi, 17LL)*Power(rxj, 16LL)*

                           (7416068831325LL + 658162968750LL*Power(rxj, 2LL) +

           11421785970LL*Power(rxj, 4LL) - 22800852LL*Power(rxj, 6LL) -

           224214LL*Power(rxj, 8LL) + 40LL*Power(rxj, 10LL)) -

                           198LL*Power(rxi, 14LL)*Power(rxj, 18LL)*

                           (12601626975LL + 2529410625LL*Power(rxj, 2LL) + 582340500LL*Power(rxj, 4LL) +

           3239250LL*Power(rxj, 6LL) + 132690LL*Power(rxj, 8LL) + 74LL*Power(rxj, 10LL)) -

                           231LL*Power(rxi, 25LL)*Power(rxj, 8LL)*

                           (21444497452125LL - 909858116250LL*Power(rxj, 2LL) +

           1447333650LL*Power(rxj, 4LL) + 178686540LL*Power(rxj, 6LL) -

           747270LL*Power(rxj, 8LL) + 184LL*Power(rxj, 10LL)) -

                           198LL*Power(rxi, 20LL)*Power(rxj, 12LL)*

                           (42449899182075LL + 4344172457625LL*Power(rxj, 2LL) -

           85249741500LL*Power(rxj, 4LL) - 1059301110LL*Power(rxj, 6LL) +

           6582370LL*Power(rxj, 8LL) + 194LL*Power(rxj, 10LL)) +

                           11LL*Power(rxi, 19LL)*Power(rxj, 14LL)*

                           (239338679943825LL + 8851966719750LL*Power(rxj, 2LL) -

           112537092150LL*Power(rxj, 4LL) - 1100275380LL*Power(rxj, 6LL) +

           2919090LL*Power(rxj, 8LL) + 248LL*Power(rxj, 10LL)) -

                           330LL*Power(rxi, 28LL)*Power(rxj, 4LL)*

                           (4860066085875LL + 2524912849305LL*Power(rxj, 2LL) -

           109538431380LL*Power(rxj, 4LL) + 1633704282LL*Power(rxj, 6LL) -

           6421278LL*Power(rxj, 8LL) + 322LL*Power(rxj, 10LL)) +

                           33LL*Power(rxi, 29LL)*Power(rxj, 4LL)*

                           (-31641507079875LL - 2157639318450LL*Power(rxj, 2LL) +

           74910015810LL*Power(rxj, 4LL) - 522003060LL*Power(rxj, 6LL) -

           250470LL*Power(rxj, 8LL) + 1288LL*Power(rxj, 10LL)) -

                           330LL*Power(rxi, 18LL)*Power(rxj, 14LL)*

                           (4867016286825LL + 1199363925375LL*Power(rxj, 2LL) +

           26817947100LL*Power(rxj, 4LL) - 167333418LL*Power(rxj, 6LL) -

           1476138LL*Power(rxj, 8LL) + 1294LL*Power(rxj, 10LL)) +

                           66LL*Power(rxi, 16LL)*Power(rxj, 16LL)*

                           (-1657759205025LL - 682207855875LL*Power(rxj, 2LL) -

           31509229500LL*Power(rxj, 4LL) - 492146550LL*Power(rxj, 6LL) -

           11910LL*Power(rxj, 8LL) + 2594LL*Power(rxj, 10LL)) +

                           1386LL*Power(rxi, 26LL)*Power(rxj, 6LL)*

                           (-6066588045375LL + 98854491375LL*Power(rxj, 2LL) -

           12496954500LL*Power(rxj, 4LL) + 420813750LL*Power(rxj, 6LL) -

           2881210LL*Power(rxj, 8LL) + 2622LL*Power(rxj, 10LL)) +

                           11LL*Power(rxi, 23LL)*Power(rxj, 10LL)*

                           (149900659402725LL - 26541339882750LL*Power(rxj, 2LL) +

           594745455150LL*Power(rxj, 4LL) - 1399125420LL*Power(rxj, 6LL) -

           7887390LL*Power(rxj, 8LL) + 4232LL*Power(rxj, 10LL)) -

                           11LL*Power(rxi, 31LL)*Power(rxj, 2LL)*

                           (7685082491625LL + 5034333946950LL*Power(rxj, 2LL) -

           108088893990LL*Power(rxj, 4LL) + 1254174300LL*Power(rxj, 6LL) -

           6355950LL*Power(rxj, 8LL) + 4232LL*Power(rxj, 10LL)) -

                           462LL*Power(rxi, 24LL)*Power(rxj, 8LL)*

                           (40495013164125LL - 3973079865375LL*Power(rxj, 2LL) +

           110288047500LL*Power(rxj, 4LL) - 381623130LL*Power(rxj, 6LL) -

           4811370LL*Power(rxj, 8LL) + 9338LL*Power(rxj, 10LL)) +

                           198LL*Power(rxi, 32LL)*(-9126526500LL - 152866565775LL*Power(rxj, 2LL) -

                                                   32383266300LL*Power(rxj, 4LL) + 709444890LL*Power(rxj, 6LL) -

                                                   5562070LL*Power(rxj, 8LL) + 11042LL*Power(rxj, 10LL)) +

                           2LL*Power(rxi, 33LL)*(-764522104500LL - 3357151476525LL*Power(rxj, 2LL) -

                                                 242177564475LL*Power(rxj, 4LL) + 4513719870LL*Power(rxj, 6LL) -

                                                 20531775LL*Power(rxj, 8LL) + 11236LL*Power(rxj, 10LL)) -

                           Power(rxi, 21LL)*Power(rxj, 12LL)*

                           (-5533525427435775LL + 138591131159250LL*Power(rxj, 2LL) +

           2815739907750LL*Power(rxj, 4LL) - 32922004500LL*Power(rxj, 6LL) +

           11347050LL*Power(rxj, 8LL) + 22472LL*Power(rxj, 10LL)) +

                           66LL*Power(rxi, 22LL)*Power(rxj, 10LL)*

                           (-283522589265825LL + 7639225988625LL*Power(rxj, 2LL) +

           480728209500LL*Power(rxj, 4LL) - 8458349130LL*Power(rxj, 6LL) +

           9771090LL*Power(rxj, 8LL) + 31786LL*Power(rxj, 10LL)) -

                           66LL*Power(rxi, 30LL)*Power(rxj, 2LL)*

                           (1678609807875LL + 4713298976925LL*Power(rxj, 2LL) -

           30578971500LL*Power(rxj, 4LL) + 53723250LL*Power(rxj, 6LL) -

           9140190LL*Power(rxj, 8LL) + 38042LL*Power(rxj, 10LL))) +

                          exp(2LL*rxi)*Power(rxi, 14LL)*

                          (-302841LL*Power(rxi, 16LL)*Power(rxj, 16LL)*

                           (-361285650LL + 1346857875LL*rxj - 1306923750LL*Power(rxj, 2LL) +

           321527250LL*Power(rxj, 3LL) + 55737000LL*Power(rxj, 4LL) -

           9297750LL*Power(rxj, 5LL) - 1843380LL*Power(rxj, 6LL) -

           50820LL*Power(rxj, 7LL) + 7340LL*Power(rxj, 8LL) + 570LL*Power(rxj, 9LL) +

           12LL*Power(rxj, 10LL)) + 12LL*Power(rxj, 32LL)*

                           (150587687250LL + 127420350750LL*rxj + 50968140300LL*Power(rxj, 2LL) +

           12742035075LL*Power(rxj, 3LL) + 2216006100LL*Power(rxj, 4LL) +

           282037140LL*Power(rxj, 5LL) + 26860680LL*Power(rxj, 6LL) +

           1918620LL*Power(rxj, 7LL) + 100980LL*Power(rxj, 8LL) + 3740LL*Power(rxj, 9LL) +

           88LL*Power(rxj, 10LL) + Power(rxj, 11LL)) -

                           3LL*Power(rxi, 32LL)*(935550LL + 1715175LL*rxj + 1559250LL*Power(rxj, 2LL) +

                                                 935550LL*Power(rxj, 3LL) + 415800LL*Power(rxj, 4LL) + 145530LL*Power(rxj, 5LL) +

                                                 41580LL*Power(rxj, 6LL) + 9900LL*Power(rxj, 7LL) + 1980LL*Power(rxj, 8LL) +

                                                 330LL*Power(rxj, 9LL) + 44LL*Power(rxj, 10LL) + 4LL*Power(rxj, 11LL)) -

                           11LL*Power(rxi, 30LL)*Power(rxj, 2LL)*

                           (-5868450LL - 10758825LL*rxj - 9780750LL*Power(rxj, 2LL) -

           5868450LL*Power(rxj, 3LL) - 2608200LL*Power(rxj, 4LL) -

           912870LL*Power(rxj, 5LL) - 260820LL*Power(rxj, 6LL) - 62100LL*Power(rxj, 7LL) -

           12420LL*Power(rxj, 8LL) - 2070LL*Power(rxj, 9LL) - 276LL*Power(rxj, 10LL) +

           8LL*Power(rxj, 11LL)) - 5313LL*Power(rxi, 14LL)*Power(rxj, 18LL)*

                           (-302299148250LL + 495525217275LL*rxj - 161894625750LL*Power(rxj, 2LL) -

           26085287250LL*Power(rxj, 3LL) + 5971779000LL*Power(rxj, 4LL) +

           1231357050LL*Power(rxj, 5LL) + 33184620LL*Power(rxj, 6LL) -

           7768980LL*Power(rxj, 7LL) - 751620LL*Power(rxj, 8LL) - 23190LL*Power(rxj, 9LL) -

           20LL*Power(rxj, 10LL) + 8LL*Power(rxj, 11LL)) +

                           5313LL*Power(rxi, 18LL)*Power(rxj, 14LL)*

                           (469625850LL - 3082655475LL*rxj + 8474631750LL*Power(rxj, 2LL) -

           6813281250LL*Power(rxj, 3LL) + 1665711000LL*Power(rxj, 4LL) +

           232996050LL*Power(rxj, 5LL) - 39477060LL*Power(rxj, 6LL) -

           6196500LL*Power(rxj, 7LL) - 121380LL*Power(rxj, 8LL) + 16330LL*Power(rxj, 9LL) +

           812LL*Power(rxj, 10LL) + 8LL*Power(rxj, 11LL)) +

                           11LL*Power(rxi, 2LL)*Power(rxj, 30LL)*

                           (10071658847250LL + 7685082491625LL*rxj + 2751598183950LL*Power(rxj, 2LL) +

           610391177550LL*Power(rxj, 3LL) + 93214459800LL*Power(rxj, 4LL) +

           10285306290LL*Power(rxj, 5LL) + 835769340LL*Power(rxj, 6LL) +

           49894380LL*Power(rxj, 7LL) + 2134620LL*Power(rxj, 8LL) +

           61770LL*Power(rxj, 9LL) + 1068LL*Power(rxj, 10LL) + 8LL*Power(rxj, 11LL)) +

                           11LL*Power(rxi, 28LL)*Power(rxj, 4LL)*

                           (-64552950LL - 118347075LL*rxj - 107588250LL*Power(rxj, 2LL) -

           64552950LL*Power(rxj, 3LL) - 28690200LL*Power(rxj, 4LL) -

           10041570LL*Power(rxj, 5LL) - 2869020LL*Power(rxj, 6LL) -

           683100LL*Power(rxj, 7LL) - 136620LL*Power(rxj, 8LL) - 33690LL*Power(rxj, 9LL) +

           1332LL*Power(rxj, 10LL) + 88LL*Power(rxj, 11LL)) -

                           11LL*Power(rxi, 4LL)*Power(rxj, 28LL)*

                           (-145801982576250LL - 94924521239625LL*rxj -

           28279793861550LL*Power(rxj, 2LL) - 5034333946950LL*Power(rxj, 3LL) -

           582898793400LL*Power(rxj, 4LL) - 44032284450LL*Power(rxj, 5LL) -

           1930850460LL*Power(rxj, 6LL) - 15289020LL*Power(rxj, 7LL) +

           3824820LL*Power(rxj, 8LL) + 253590LL*Power(rxj, 9LL) + 7380LL*Power(rxj, 10LL) +

           88LL*Power(rxj, 11LL)) - 253LL*Power(rxi, 20LL)*Power(rxj, 12LL)*

                           (1119853350LL + 2437660575LL*rxj - 1979538750LL*Power(rxj, 2LL) +

           10991153250LL*Power(rxj, 3LL) - 8219799000LL*Power(rxj, 4LL) +

           2482996950LL*Power(rxj, 5LL) + 218260980LL*Power(rxj, 6LL) -

           47838060LL*Power(rxj, 7LL) - 5151420LL*Power(rxj, 8LL) -

           44850LL*Power(rxj, 9LL) + 8292LL*Power(rxj, 10LL) + 184LL*Power(rxj, 11LL)) +

                           253LL*Power(rxi, 12LL)*Power(rxj, 20LL)*

                           (33221660229450LL - 21871642005675LL*rxj - 1992841562250LL*Power(rxj, 2LL) +

           1153971299250LL*Power(rxj, 3LL) + 201395565000LL*Power(rxj, 4LL) +

           1321478550LL*Power(rxj, 5LL) - 2305327500LL*Power(rxj, 6LL) -

           232090380LL*Power(rxj, 7LL) - 8375580LL*Power(rxj, 8LL) +

           32670LL*Power(rxj, 9LL) + 9924LL*Power(rxj, 10LL) + 184LL*Power(rxj, 11LL)) +

                           11LL*Power(rxi, 6LL)*Power(rxj, 26LL)*

                           (764390093717250LL + 377817065789625LL*rxj +

           75747385479150LL*Power(rxj, 2LL) + 6472917955350LL*Power(rxj, 3LL) -

           183473829000LL*Power(rxj, 4LL) - 108088893990LL*Power(rxj, 5LL) -

           12770008020LL*Power(rxj, 6LL) - 820676340LL*Power(rxj, 7LL) -

           29919780LL*Power(rxj, 8LL) - 471270LL*Power(rxj, 9LL) +

           4236LL*Power(rxj, 10LL) + 200LL*Power(rxj, 11LL)) -

                           11LL*Power(rxi, 26LL)*Power(rxj, 6LL)*

                           (-451870650LL - 828429525LL*rxj - 753117750LL*Power(rxj, 2LL) -

           451870650LL*Power(rxj, 3LL) - 200831400LL*Power(rxj, 4LL) -

           70290990LL*Power(rxj, 5LL) - 20083140LL*Power(rxj, 6LL) -

           3349620LL*Power(rxj, 7LL) - 2388420LL*Power(rxj, 8LL) +

           66810LL*Power(rxj, 9LL) + 15564LL*Power(rxj, 10LL) + 200LL*Power(rxj, 11LL)) -

                           11LL*Power(rxi, 24LL)*Power(rxj, 8LL)*

                           (2259353250LL + 4142147625LL*rxj + 3765588750LL*Power(rxj, 2LL) +

           2259353250LL*Power(rxj, 3LL) + 1004157000LL*Power(rxj, 4LL) +

           430816050LL*Power(rxj, 5LL) - 58306500LL*Power(rxj, 6LL) +

           104343660LL*Power(rxj, 7LL) - 71460LL*Power(rxj, 8LL) -

           1121070LL*Power(rxj, 9LL) - 38820LL*Power(rxj, 10LL) + 248LL*Power(rxj, 11LL)) +

                           11LL*Power(rxi, 8LL)*Power(rxj, 24LL)*

                           (1700790552893250LL + 450334446494625LL*rxj -

           12455665913250LL*Power(rxj, 2LL) - 19774531113750LL*Power(rxj, 3LL) -

           3286152941400LL*Power(rxj, 4LL) - 224730047430LL*Power(rxj, 5LL) +

           322339500LL*Power(rxj, 6LL) + 1254174300LL*Power(rxj, 7LL) +

           100117260LL*Power(rxj, 8LL) + 3733050LL*Power(rxj, 9LL) +

           63372LL*Power(rxj, 10LL) + 248LL*Power(rxj, 11LL)) +

                           Power(rxi, 22LL)*Power(rxj, 10LL)*

                           (94440965850LL + 173141770725LL*rxj + 157401609750LL*Power(rxj, 2LL) +

           76108551750LL*Power(rxj, 3LL) + 115303419000LL*Power(rxj, 4LL) -

           67892343750LL*Power(rxj, 5LL) + 32481672300LL*Power(rxj, 6LL) +

           1254046860LL*Power(rxj, 7LL) - 487125540LL*Power(rxj, 8LL) -

           32109990LL*Power(rxj, 9LL) + 38412LL*Power(rxj, 10LL) + 22472LL*Power(rxj, 11LL))

                           - Power(rxi, 10LL)*Power(rxj, 22LL)*(-18712490891544450LL + 1648907253429975LL*rxj +

                                                                1835562897803250LL*Power(rxj, 2LL) + 210177224853750LL*Power(rxj, 3LL) -

                                                                17320778937000LL*Power(rxj, 4LL) - 5914505623950LL*Power(rxj, 5LL) -

                                                                539122413060LL*Power(rxj, 6LL) - 17226100980LL*Power(rxj, 7LL) +

                                                                603252540LL*Power(rxj, 8LL) + 69915450LL*Power(rxj, 9LL) +

                                                                2186316LL*Power(rxj, 10LL) + 22472LL*Power(rxj, 11LL))))/

                         (2.80665e6*exp(2LL*(rxi + rxj))*Power(rxi - rxj, 23LL)*Power(rxi + rxj, 23LL))

                         );
        }

    }
    return S;
}
