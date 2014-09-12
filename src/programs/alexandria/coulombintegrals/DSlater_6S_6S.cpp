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

cl_R DSlater_6S_6S(cl_R r, cl_R xi, cl_R xj)
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
            S = (1LL + (-1LL - (15604461LL*r*xi)/8.388608e6 -

                        (7215853LL*Power(r, 2LL)*Power(xi, 2LL))/4.194304e6 -

                        (19903123LL*Power(r, 3LL)*Power(xi, 3LL))/1.8874368e7 -

                        (1136755LL*Power(r, 4LL)*Power(xi, 4LL))/2.359296e6 -

                        (309683LL*Power(r, 5LL)*Power(xi, 5LL))/1.769472e6 -

                        (13797LL*Power(r, 6LL)*Power(xi, 6LL))/262144LL -

                        (238423LL*Power(r, 7LL)*Power(xi, 7LL))/1.769472e7 -

                        (139189LL*Power(r, 8LL)*Power(xi, 8LL))/4.644864e7 -

                        (163811LL*Power(r, 9LL)*Power(xi, 9LL))/2.7869184e8 -

                        (71677LL*Power(r, 10LL)*Power(xi, 10LL))/6.967296e8 -

                        (186367LL*Power(r, 11LL)*Power(xi, 11LL))/1.14960384e10 -

                        (13LL*Power(r, 12LL)*Power(xi, 12LL))/5.6133e6 -

                        (Power(r, 13LL)*Power(xi, 13LL))/3.31695e6 -

                        (Power(r, 14LL)*Power(xi, 14LL))/2.786238e7 -

                        (Power(r, 15LL)*Power(xi, 15LL))/2.5540515e8 -

                        (Power(r, 16LL)*Power(xi, 16LL))/2.5540515e9 -

                        (Power(r, 17LL)*Power(xi, 17LL))/2.791213425e10 -

                        (Power(r, 18LL)*Power(xi, 18LL))/3.34945611e11 -

                        (Power(r, 19LL)*Power(xi, 19LL))/4.4547766263e12 -

                        (Power(r, 20LL)*Power(xi, 20LL))/6.68216493945e13 -

                        (Power(r, 21LL)*Power(xi, 21LL))/1.16937886440375e15 -

                        (Power(r, 22LL)*Power(xi, 22LL))/2.57263350168825e16 -

                        (Power(r, 23LL)*Power(xi, 23LL))/8.875585580824462e17)/exp(2LL*r*xi))/

                Power(r, 2LL) - (((-15604461LL*xi)/8.388608e6 -

                                  (7215853LL*r*Power(xi, 2LL))/2.097152e6 -

                                  (19903123LL*Power(r, 2LL)*Power(xi, 3LL))/6.291456e6 -

                                  (1136755LL*Power(r, 3LL)*Power(xi, 4LL))/589824LL -

                                  (1548415LL*Power(r, 4LL)*Power(xi, 5LL))/1.769472e6 -

                                  (41391LL*Power(r, 5LL)*Power(xi, 6LL))/131072LL -

                                  (1668961LL*Power(r, 6LL)*Power(xi, 7LL))/1.769472e7 -

                                  (139189LL*Power(r, 7LL)*Power(xi, 8LL))/5.80608e6 -

                                  (163811LL*Power(r, 8LL)*Power(xi, 9LL))/3.096576e7 -

                                  (71677LL*Power(r, 9LL)*Power(xi, 10LL))/6.967296e7 -

                                  (186367LL*Power(r, 10LL)*Power(xi, 11LL))/1.0450944e9 -

                                  (13LL*Power(r, 11LL)*Power(xi, 12LL))/467775LL -

                                  (Power(r, 12LL)*Power(xi, 13LL))/255150LL -

                                  (Power(r, 13LL)*Power(xi, 14LL))/1.99017e6 -

                                  (Power(r, 14LL)*Power(xi, 15LL))/1.702701e7 -

                                  (4LL*Power(r, 15LL)*Power(xi, 16LL))/6.38512875e8 -

                                  (Power(r, 16LL)*Power(xi, 17LL))/1.64189025e9 -

                                  (Power(r, 17LL)*Power(xi, 18LL))/1.86080895e10 -

                                  (Power(r, 18LL)*Power(xi, 19LL))/2.344619277e11 -

                                  (Power(r, 19LL)*Power(xi, 20LL))/3.341082469725e12 -

                                  (Power(r, 20LL)*Power(xi, 21LL))/5.568470782875e13 -

                                  (Power(r, 21LL)*Power(xi, 22LL))/1.16937886440375e15 -

                                  (Power(r, 22LL)*Power(xi, 23LL))/3.858950252532375e16)/exp(2LL*r*xi) -

                                 (2LL*xi*(-1LL - (15604461LL*r*xi)/8.388608e6 -

                                          (7215853LL*Power(r, 2LL)*Power(xi, 2LL))/4.194304e6 -

                                          (19903123LL*Power(r, 3LL)*Power(xi, 3LL))/1.8874368e7 -

                                          (1136755LL*Power(r, 4LL)*Power(xi, 4LL))/2.359296e6 -

                                          (309683LL*Power(r, 5LL)*Power(xi, 5LL))/1.769472e6 -

                                          (13797LL*Power(r, 6LL)*Power(xi, 6LL))/262144LL -

                                          (238423LL*Power(r, 7LL)*Power(xi, 7LL))/1.769472e7 -

                                          (139189LL*Power(r, 8LL)*Power(xi, 8LL))/4.644864e7 -

                                          (163811LL*Power(r, 9LL)*Power(xi, 9LL))/2.7869184e8 -

                                          (71677LL*Power(r, 10LL)*Power(xi, 10LL))/6.967296e8 -

                                          (186367LL*Power(r, 11LL)*Power(xi, 11LL))/1.14960384e10 -

                                          (13LL*Power(r, 12LL)*Power(xi, 12LL))/5.6133e6 -

                                          (Power(r, 13LL)*Power(xi, 13LL))/3.31695e6 -

                                          (Power(r, 14LL)*Power(xi, 14LL))/2.786238e7 -

                                          (Power(r, 15LL)*Power(xi, 15LL))/2.5540515e8 -

                                          (Power(r, 16LL)*Power(xi, 16LL))/2.5540515e9 -

                                          (Power(r, 17LL)*Power(xi, 17LL))/2.791213425e10 -

                                          (Power(r, 18LL)*Power(xi, 18LL))/3.34945611e11 -

                                          (Power(r, 19LL)*Power(xi, 19LL))/4.4547766263e12 -

                                          (Power(r, 20LL)*Power(xi, 20LL))/6.68216493945e13 -

                                          (Power(r, 21LL)*Power(xi, 21LL))/1.16937886440375e15 -

                                          (Power(r, 22LL)*Power(xi, 22LL))/2.57263350168825e16 -

                                          (Power(r, 23LL)*Power(xi, 23LL))/8.875585580824462e17))/exp(2LL*r*xi))/r

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
            S = (2806650LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 23LL) +

                 exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-1056LL*Power(r, 10LL)*Power(xi, 42LL) - 12LL*Power(r, 11LL)*Power(xi, 43LL) +

                  2806650LL*Power(xj, 32LL) + 5145525LL*r*xi*Power(xj, 32LL) -

                  88LL*Power(r, 9LL)*Power(xi, 41LL)*(510LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  935550LL*Power(xi, 2LL)*Power(xj, 30LL)*(-69LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  467775LL*r*Power(xi, 3LL)*Power(xj, 30LL)*

                  (-253LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  132LL*Power(r, 8LL)*Power(xi, 40LL)*(9180LL + 89LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  311850LL*Power(xi, 4LL)*Power(xj, 28LL)*

                  (2277LL - 345LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  31185LL*r*Power(xi, 5LL)*Power(xj, 28LL)*

                  (41745LL - 2070LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

                  + 1980LL*Power(r, 6LL)*Power(xi, 38LL)*

                  (-162792LL - 11859LL*Power(r, 2LL)*Power(xj, 2LL) +

            41LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  22LL*Power(r, 7LL)*Power(xi, 39LL)*

                  (-1046520LL - 30885LL*Power(r, 2LL)*Power(xj, 2LL) +

            44LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  62370LL*Power(xi, 6LL)*Power(xj, 26LL)*

                  (-79695LL + 18975LL*Power(r, 2LL)*Power(xj, 2LL) -

            460LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  110LL*Power(r, 5LL)*Power(xi, 37LL)*

                  (30767688LL + 4989438LL*Power(r, 2LL)*Power(xj, 2LL) -

            25359LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  1485LL*r*Power(xi, 7LL)*Power(xj, 26LL)*

                  (-6136515LL + 478170LL*Power(r, 2LL)*Power(xj, 2LL) -

            6762LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  132LL*Power(r, 4LL)*Power(xi, 36LL)*

                  (201455100LL + 69647445LL*Power(r, 2LL)*Power(xj, 2LL) -

            318735LL*Power(r, 4LL)*Power(xj, 4LL) + 353LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  495LL*r*Power(xi, 9LL)*Power(xj, 24LL)*

                  (92047725LL - 10041570LL*Power(r, 2LL)*Power(xj, 2LL) +

            223146LL*Power(r, 4LL)*Power(xj, 4LL) - 1380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  2970LL*Power(xi, 8LL)*Power(xj, 24LL)*

                  (8367975LL - 2789325LL*Power(r, 2LL)*Power(xj, 2LL) +

            106260LL*Power(r, 4LL)*Power(xj, 4LL) - 966LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  22LL*Power(r, 3LL)*Power(xi, 35LL)*

                  (6950200950LL + 5142653145LL*Power(r, 2LL)*Power(xj, 2LL) +

            7644510LL*Power(r, 4LL)*Power(xj, 4LL) -

            235635LL*Power(r, 6LL)*Power(xj, 6LL) + 124LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  132LL*Power(r, 2LL)*Power(xi, 34LL)*

                  (4633467300LL + 7767871650LL*Power(r, 2LL)*Power(xj, 2LL) +

            160904205LL*Power(r, 4LL)*Power(xj, 4LL) -

            2493315LL*Power(r, 6LL)*Power(xj, 6LL) + 5281LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  495LL*r*Power(xi, 27LL)*Power(xj, 6LL)*

                  (8395934795325LL - 439434024750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11948496210LL*Power(r, 4LL)*Power(xj, 4LL) -

            118623972LL*Power(r, 6LL)*Power(xj, 6LL) +

            248906LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  11LL*r*Power(xi, 15LL)*Power(xj, 18LL)*

                  (1488922594425LL + 252796524750LL*Power(r, 2LL)*Power(xj, 2LL) +

            6172031250LL*Power(r, 4LL)*Power(xj, 4LL) +

            104343660LL*Power(r, 6LL)*Power(xj, 6LL) +

            66810LL*Power(r, 8LL)*Power(xj, 8LL) - 88LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  66LL*Power(xi, 10LL)*Power(xj, 22LL)*

                  (-1430923725LL + 627598125LL*Power(r, 2LL)*Power(xj, 2LL) -

            33471900LL*Power(r, 4LL)*Power(xj, 4LL) +

            478170LL*Power(r, 6LL)*Power(xj, 6LL) - 2070LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 10LL)*Power(xj, 10LL)) -

                  1518LL*Power(xi, 12LL)*Power(xj, 20LL)*

                  (-186642225LL + 103690125LL*Power(r, 2LL)*Power(xj, 2LL) -

            7276500LL*Power(r, 4LL)*Power(xj, 4LL) +

            145530LL*Power(r, 6LL)*Power(xj, 6LL) - 990LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  3LL*r*Power(xi, 11LL)*Power(xj, 22LL)*

                  (-57713923575LL + 8284295250LL*Power(r, 2LL)*Power(xj, 2LL) -

            257733630LL*Power(r, 4LL)*Power(xj, 4LL) +

            2504700LL*Power(r, 6LL)*Power(xj, 6LL) - 7590LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  11LL*r*Power(xi, 13LL)*Power(xj, 20LL)*

                  (56066193225LL - 6918959250LL*Power(r, 2LL)*Power(xj, 2LL) +

            430816050LL*Power(r, 4LL)*Power(xj, 4LL) -

            3349620LL*Power(r, 6LL)*Power(xj, 6LL) +

            33690LL*Power(r, 8LL)*Power(xj, 8LL) + 8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  55LL*r*Power(xi, 17LL)*Power(xj, 16LL)*

                  (7416068831325LL + 658162968750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11421785970LL*Power(r, 4LL)*Power(xj, 4LL) -

            22800852LL*Power(r, 6LL)*Power(xj, 6LL) -

            224214LL*Power(r, 8LL)*Power(xj, 8LL) + 40LL*Power(r, 10LL)*Power(xj, 10LL)) -

                  198LL*Power(xi, 14LL)*Power(xj, 18LL)*

                  (12601626975LL + 2529410625LL*Power(r, 2LL)*Power(xj, 2LL) +

            582340500LL*Power(r, 4LL)*Power(xj, 4LL) +

            3239250LL*Power(r, 6LL)*Power(xj, 6LL) +

            132690LL*Power(r, 8LL)*Power(xj, 8LL) + 74LL*Power(r, 10LL)*Power(xj, 10LL)) -

                  231LL*r*Power(xi, 25LL)*Power(xj, 8LL)*

                  (21444497452125LL - 909858116250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1447333650LL*Power(r, 4LL)*Power(xj, 4LL) +

            178686540LL*Power(r, 6LL)*Power(xj, 6LL) -

            747270LL*Power(r, 8LL)*Power(xj, 8LL) + 184LL*Power(r, 10LL)*Power(xj, 10LL)) -

                  198LL*Power(xi, 20LL)*Power(xj, 12LL)*

                  (42449899182075LL + 4344172457625LL*Power(r, 2LL)*Power(xj, 2LL) -

            85249741500LL*Power(r, 4LL)*Power(xj, 4LL) -

            1059301110LL*Power(r, 6LL)*Power(xj, 6LL) +

            6582370LL*Power(r, 8LL)*Power(xj, 8LL) + 194LL*Power(r, 10LL)*Power(xj, 10LL))

                  + 11LL*r*Power(xi, 19LL)*Power(xj, 14LL)*

                  (239338679943825LL + 8851966719750LL*Power(r, 2LL)*Power(xj, 2LL) -

            112537092150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1100275380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2919090LL*Power(r, 8LL)*Power(xj, 8LL) + 248LL*Power(r, 10LL)*Power(xj, 10LL))

                  - 330LL*Power(xi, 28LL)*Power(xj, 4LL)*

                  (4860066085875LL + 2524912849305LL*Power(r, 2LL)*Power(xj, 2LL) -

            109538431380LL*Power(r, 4LL)*Power(xj, 4LL) +

            1633704282LL*Power(r, 6LL)*Power(xj, 6LL) -

            6421278LL*Power(r, 8LL)*Power(xj, 8LL) + 322LL*Power(r, 10LL)*Power(xj, 10LL))

                  + 33LL*r*Power(xi, 29LL)*Power(xj, 4LL)*

                  (-31641507079875LL - 2157639318450LL*Power(r, 2LL)*Power(xj, 2LL) +

            74910015810LL*Power(r, 4LL)*Power(xj, 4LL) -

            522003060LL*Power(r, 6LL)*Power(xj, 6LL) -

            250470LL*Power(r, 8LL)*Power(xj, 8LL) + 1288LL*Power(r, 10LL)*Power(xj, 10LL))

                  - 330LL*Power(xi, 18LL)*Power(xj, 14LL)*

                  (4867016286825LL + 1199363925375LL*Power(r, 2LL)*Power(xj, 2LL) +

            26817947100LL*Power(r, 4LL)*Power(xj, 4LL) -

            167333418LL*Power(r, 6LL)*Power(xj, 6LL) -

            1476138LL*Power(r, 8LL)*Power(xj, 8LL) + 1294LL*Power(r, 10LL)*Power(xj, 10LL))

                  + 66LL*Power(xi, 16LL)*Power(xj, 16LL)*

                  (-1657759205025LL - 682207855875LL*Power(r, 2LL)*Power(xj, 2LL) -

            31509229500LL*Power(r, 4LL)*Power(xj, 4LL) -

            492146550LL*Power(r, 6LL)*Power(xj, 6LL) -

            11910LL*Power(r, 8LL)*Power(xj, 8LL) + 2594LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  1386LL*Power(xi, 26LL)*Power(xj, 6LL)*

                  (-6066588045375LL + 98854491375LL*Power(r, 2LL)*Power(xj, 2LL) -

            12496954500LL*Power(r, 4LL)*Power(xj, 4LL) +

            420813750LL*Power(r, 6LL)*Power(xj, 6LL) -

            2881210LL*Power(r, 8LL)*Power(xj, 8LL) + 2622LL*Power(r, 10LL)*Power(xj, 10LL))

                  + 11LL*r*Power(xi, 23LL)*Power(xj, 10LL)*

                  (149900659402725LL - 26541339882750LL*Power(r, 2LL)*Power(xj, 2LL) +

            594745455150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1399125420LL*Power(r, 6LL)*Power(xj, 6LL) -

            7887390LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                  - 11LL*r*Power(xi, 31LL)*Power(xj, 2LL)*

                  (7685082491625LL + 5034333946950LL*Power(r, 2LL)*Power(xj, 2LL) -

            108088893990LL*Power(r, 4LL)*Power(xj, 4LL) +

            1254174300LL*Power(r, 6LL)*Power(xj, 6LL) -

            6355950LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                  - 462LL*Power(xi, 24LL)*Power(xj, 8LL)*

                  (40495013164125LL - 3973079865375LL*Power(r, 2LL)*Power(xj, 2LL) +

            110288047500LL*Power(r, 4LL)*Power(xj, 4LL) -

            381623130LL*Power(r, 6LL)*Power(xj, 6LL) -

            4811370LL*Power(r, 8LL)*Power(xj, 8LL) + 9338LL*Power(r, 10LL)*Power(xj, 10LL))

                  + 198LL*Power(xi, 32LL)*(-9126526500LL - 152866565775LL*Power(r, 2LL)*Power(xj, 2LL) -

                                           32383266300LL*Power(r, 4LL)*Power(xj, 4LL) +

                                           709444890LL*Power(r, 6LL)*Power(xj, 6LL) -

                                           5562070LL*Power(r, 8LL)*Power(xj, 8LL) + 11042LL*Power(r, 10LL)*Power(xj, 10LL)

                                           ) + 2LL*r*Power(xi, 33LL)*(-764522104500LL -

                                                                      3357151476525LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                      242177564475LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                      4513719870LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                      20531775LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                      11236LL*Power(r, 10LL)*Power(xj, 10LL)) -

                  r*Power(xi, 21LL)*Power(xj, 12LL)*

                  (-5533525427435775LL + 138591131159250LL*Power(r, 2LL)*Power(xj, 2LL) +

            2815739907750LL*Power(r, 4LL)*Power(xj, 4LL) -

            32922004500LL*Power(r, 6LL)*Power(xj, 6LL) +

            11347050LL*Power(r, 8LL)*Power(xj, 8LL) +

            22472LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  66LL*Power(xi, 22LL)*Power(xj, 10LL)*

                  (-283522589265825LL + 7639225988625LL*Power(r, 2LL)*Power(xj, 2LL) +

            480728209500LL*Power(r, 4LL)*Power(xj, 4LL) -

            8458349130LL*Power(r, 6LL)*Power(xj, 6LL) +

            9771090LL*Power(r, 8LL)*Power(xj, 8LL) + 31786LL*Power(r, 10LL)*Power(xj, 10LL)

                  ) - 66LL*Power(xi, 30LL)*Power(xj, 2LL)*

                  (1678609807875LL + 4713298976925LL*Power(r, 2LL)*Power(xj, 2LL) -

            30578971500LL*Power(r, 4LL)*Power(xj, 4LL) +

            53723250LL*Power(r, 6LL)*Power(xj, 6LL) -

            9140190LL*Power(r, 8LL)*Power(xj, 8LL) + 38042LL*Power(r, 10LL)*Power(xj, 10LL))

                 ) + exp(2LL*r*xi)*Power(xi, 14LL)*

                 (-302841LL*Power(xi, 16LL)*Power(xj, 16LL)*

                  (-361285650LL + 1346857875LL*r*xj -

            1306923750LL*Power(r, 2LL)*Power(xj, 2LL) +

            321527250LL*Power(r, 3LL)*Power(xj, 3LL) +

            55737000LL*Power(r, 4LL)*Power(xj, 4LL) -

            9297750LL*Power(r, 5LL)*Power(xj, 5LL) -

            1843380LL*Power(r, 6LL)*Power(xj, 6LL) - 50820LL*Power(r, 7LL)*Power(xj, 7LL) +

            7340LL*Power(r, 8LL)*Power(xj, 8LL) + 570LL*Power(r, 9LL)*Power(xj, 9LL) +

            12LL*Power(r, 10LL)*Power(xj, 10LL)) +

                  12LL*Power(xj, 32LL)*(150587687250LL + 127420350750LL*r*xj +

                                        50968140300LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        12742035075LL*Power(r, 3LL)*Power(xj, 3LL) +

                                        2216006100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                        282037140LL*Power(r, 5LL)*Power(xj, 5LL) +

                                        26860680LL*Power(r, 6LL)*Power(xj, 6LL) +

                                        1918620LL*Power(r, 7LL)*Power(xj, 7LL) +

                                        100980LL*Power(r, 8LL)*Power(xj, 8LL) + 3740LL*Power(r, 9LL)*Power(xj, 9LL) +

                                        88LL*Power(r, 10LL)*Power(xj, 10LL) + Power(r, 11LL)*Power(xj, 11LL)) -

                  3LL*Power(xi, 32LL)*(935550LL + 1715175LL*r*xj +

                                       1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                       330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                       4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  11LL*Power(xi, 30LL)*Power(xj, 2LL)*

                  (-5868450LL - 10758825LL*r*xj - 9780750LL*Power(r, 2LL)*Power(xj, 2LL) -

            5868450LL*Power(r, 3LL)*Power(xj, 3LL) -

            2608200LL*Power(r, 4LL)*Power(xj, 4LL) -

            912870LL*Power(r, 5LL)*Power(xj, 5LL) - 260820LL*Power(r, 6LL)*Power(xj, 6LL) -

            62100LL*Power(r, 7LL)*Power(xj, 7LL) - 12420LL*Power(r, 8LL)*Power(xj, 8LL) -

            2070LL*Power(r, 9LL)*Power(xj, 9LL) - 276LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  5313LL*Power(xi, 14LL)*Power(xj, 18LL)*

                  (-302299148250LL + 495525217275LL*r*xj -

            161894625750LL*Power(r, 2LL)*Power(xj, 2LL) -

            26085287250LL*Power(r, 3LL)*Power(xj, 3LL) +

            5971779000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1231357050LL*Power(r, 5LL)*Power(xj, 5LL) +

            33184620LL*Power(r, 6LL)*Power(xj, 6LL) -

            7768980LL*Power(r, 7LL)*Power(xj, 7LL) -

            751620LL*Power(r, 8LL)*Power(xj, 8LL) - 23190LL*Power(r, 9LL)*Power(xj, 9LL) -

            20LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  5313LL*Power(xi, 18LL)*Power(xj, 14LL)*

                  (469625850LL - 3082655475LL*r*xj +

            8474631750LL*Power(r, 2LL)*Power(xj, 2LL) -

            6813281250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1665711000LL*Power(r, 4LL)*Power(xj, 4LL) +

            232996050LL*Power(r, 5LL)*Power(xj, 5LL) -

            39477060LL*Power(r, 6LL)*Power(xj, 6LL) -

            6196500LL*Power(r, 7LL)*Power(xj, 7LL) -

            121380LL*Power(r, 8LL)*Power(xj, 8LL) + 16330LL*Power(r, 9LL)*Power(xj, 9LL) +

            812LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  11LL*Power(xi, 2LL)*Power(xj, 30LL)*

                  (10071658847250LL + 7685082491625LL*r*xj +

            2751598183950LL*Power(r, 2LL)*Power(xj, 2LL) +

            610391177550LL*Power(r, 3LL)*Power(xj, 3LL) +

            93214459800LL*Power(r, 4LL)*Power(xj, 4LL) +

            10285306290LL*Power(r, 5LL)*Power(xj, 5LL) +

            835769340LL*Power(r, 6LL)*Power(xj, 6LL) +

            49894380LL*Power(r, 7LL)*Power(xj, 7LL) +

            2134620LL*Power(r, 8LL)*Power(xj, 8LL) + 61770LL*Power(r, 9LL)*Power(xj, 9LL) +

            1068LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  11LL*Power(xi, 28LL)*Power(xj, 4LL)*

                  (-64552950LL - 118347075LL*r*xj - 107588250LL*Power(r, 2LL)*Power(xj, 2LL) -

            64552950LL*Power(r, 3LL)*Power(xj, 3LL) -

            28690200LL*Power(r, 4LL)*Power(xj, 4LL) -

            10041570LL*Power(r, 5LL)*Power(xj, 5LL) -

            2869020LL*Power(r, 6LL)*Power(xj, 6LL) -

            683100LL*Power(r, 7LL)*Power(xj, 7LL) - 136620LL*Power(r, 8LL)*Power(xj, 8LL) -

            33690LL*Power(r, 9LL)*Power(xj, 9LL) + 1332LL*Power(r, 10LL)*Power(xj, 10LL) +

            88LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  11LL*Power(xi, 4LL)*Power(xj, 28LL)*

                  (-145801982576250LL - 94924521239625LL*r*xj -

            28279793861550LL*Power(r, 2LL)*Power(xj, 2LL) -

            5034333946950LL*Power(r, 3LL)*Power(xj, 3LL) -

            582898793400LL*Power(r, 4LL)*Power(xj, 4LL) -

            44032284450LL*Power(r, 5LL)*Power(xj, 5LL) -

            1930850460LL*Power(r, 6LL)*Power(xj, 6LL) -

            15289020LL*Power(r, 7LL)*Power(xj, 7LL) +

            3824820LL*Power(r, 8LL)*Power(xj, 8LL) +

            253590LL*Power(r, 9LL)*Power(xj, 9LL) + 7380LL*Power(r, 10LL)*Power(xj, 10LL) +

            88LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  253LL*Power(xi, 20LL)*Power(xj, 12LL)*

                  (1119853350LL + 2437660575LL*r*xj -

            1979538750LL*Power(r, 2LL)*Power(xj, 2LL) +

            10991153250LL*Power(r, 3LL)*Power(xj, 3LL) -

            8219799000LL*Power(r, 4LL)*Power(xj, 4LL) +

            2482996950LL*Power(r, 5LL)*Power(xj, 5LL) +

            218260980LL*Power(r, 6LL)*Power(xj, 6LL) -

            47838060LL*Power(r, 7LL)*Power(xj, 7LL) -

            5151420LL*Power(r, 8LL)*Power(xj, 8LL) - 44850LL*Power(r, 9LL)*Power(xj, 9LL) +

            8292LL*Power(r, 10LL)*Power(xj, 10LL) + 184LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  253LL*Power(xi, 12LL)*Power(xj, 20LL)*

                  (33221660229450LL - 21871642005675LL*r*xj -

            1992841562250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1153971299250LL*Power(r, 3LL)*Power(xj, 3LL) +

            201395565000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1321478550LL*Power(r, 5LL)*Power(xj, 5LL) -

            2305327500LL*Power(r, 6LL)*Power(xj, 6LL) -

            232090380LL*Power(r, 7LL)*Power(xj, 7LL) -

            8375580LL*Power(r, 8LL)*Power(xj, 8LL) + 32670LL*Power(r, 9LL)*Power(xj, 9LL) +

            9924LL*Power(r, 10LL)*Power(xj, 10LL) + 184LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  11LL*Power(xi, 6LL)*Power(xj, 26LL)*

                  (764390093717250LL + 377817065789625LL*r*xj +

            75747385479150LL*Power(r, 2LL)*Power(xj, 2LL) +

            6472917955350LL*Power(r, 3LL)*Power(xj, 3LL) -

            183473829000LL*Power(r, 4LL)*Power(xj, 4LL) -

            108088893990LL*Power(r, 5LL)*Power(xj, 5LL) -

            12770008020LL*Power(r, 6LL)*Power(xj, 6LL) -

            820676340LL*Power(r, 7LL)*Power(xj, 7LL) -

            29919780LL*Power(r, 8LL)*Power(xj, 8LL) -

            471270LL*Power(r, 9LL)*Power(xj, 9LL) + 4236LL*Power(r, 10LL)*Power(xj, 10LL) +

            200LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  11LL*Power(xi, 26LL)*Power(xj, 6LL)*

                  (-451870650LL - 828429525LL*r*xj - 753117750LL*Power(r, 2LL)*Power(xj, 2LL) -

            451870650LL*Power(r, 3LL)*Power(xj, 3LL) -

            200831400LL*Power(r, 4LL)*Power(xj, 4LL) -

            70290990LL*Power(r, 5LL)*Power(xj, 5LL) -

            20083140LL*Power(r, 6LL)*Power(xj, 6LL) -

            3349620LL*Power(r, 7LL)*Power(xj, 7LL) -

            2388420LL*Power(r, 8LL)*Power(xj, 8LL) + 66810LL*Power(r, 9LL)*Power(xj, 9LL) +

            15564LL*Power(r, 10LL)*Power(xj, 10LL) + 200LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  11LL*Power(xi, 24LL)*Power(xj, 8LL)*

                  (2259353250LL + 4142147625LL*r*xj +

            3765588750LL*Power(r, 2LL)*Power(xj, 2LL) +

            2259353250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1004157000LL*Power(r, 4LL)*Power(xj, 4LL) +

            430816050LL*Power(r, 5LL)*Power(xj, 5LL) -

            58306500LL*Power(r, 6LL)*Power(xj, 6LL) +

            104343660LL*Power(r, 7LL)*Power(xj, 7LL) -

            71460LL*Power(r, 8LL)*Power(xj, 8LL) - 1121070LL*Power(r, 9LL)*Power(xj, 9LL) -

            38820LL*Power(r, 10LL)*Power(xj, 10LL) + 248LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  11LL*Power(xi, 8LL)*Power(xj, 24LL)*

                  (1700790552893250LL + 450334446494625LL*r*xj -

            12455665913250LL*Power(r, 2LL)*Power(xj, 2LL) -

            19774531113750LL*Power(r, 3LL)*Power(xj, 3LL) -

            3286152941400LL*Power(r, 4LL)*Power(xj, 4LL) -

            224730047430LL*Power(r, 5LL)*Power(xj, 5LL) +

            322339500LL*Power(r, 6LL)*Power(xj, 6LL) +

            1254174300LL*Power(r, 7LL)*Power(xj, 7LL) +

            100117260LL*Power(r, 8LL)*Power(xj, 8LL) +

            3733050LL*Power(r, 9LL)*Power(xj, 9LL) +

            63372LL*Power(r, 10LL)*Power(xj, 10LL) + 248LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  Power(xi, 22LL)*Power(xj, 10LL)*

                  (94440965850LL + 173141770725LL*r*xj +

            157401609750LL*Power(r, 2LL)*Power(xj, 2LL) +

            76108551750LL*Power(r, 3LL)*Power(xj, 3LL) +

            115303419000LL*Power(r, 4LL)*Power(xj, 4LL) -

            67892343750LL*Power(r, 5LL)*Power(xj, 5LL) +

            32481672300LL*Power(r, 6LL)*Power(xj, 6LL) +

            1254046860LL*Power(r, 7LL)*Power(xj, 7LL) -

            487125540LL*Power(r, 8LL)*Power(xj, 8LL) -

            32109990LL*Power(r, 9LL)*Power(xj, 9LL) +

            38412LL*Power(r, 10LL)*Power(xj, 10LL) + 22472LL*Power(r, 11LL)*Power(xj, 11LL))

                  - Power(xi, 10LL)*Power(xj, 22LL)*(-18712490891544450LL + 1648907253429975LL*r*xj +

                                                     1835562897803250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                     210177224853750LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                     17320778937000LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                     5914505623950LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                     539122413060LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                     17226100980LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                     603252540LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                     69915450LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                     2186316LL*Power(r, 10LL)*Power(xj, 10LL) + 22472LL*Power(r, 11LL)*Power(xj, 11LL)

                                                     )))/(2.80665e6*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 23LL)*

                                                          Power(xi + xj, 23LL)) + (2806650LL*exp(2LL*r*(xi + xj))*

                                                                                   Power(Power(xi, 2LL) - Power(xj, 2LL), 23LL) +

                                                                                   exp(2LL*r*xj)*Power(xj, 14LL)*

                                                                                   (-1056LL*Power(r, 10LL)*Power(xi, 42LL) - 12LL*Power(r, 11LL)*Power(xi, 43LL) +

                                                                             2806650LL*Power(xj, 32LL) + 5145525LL*r*xi*Power(xj, 32LL) -

                                                                             88LL*Power(r, 9LL)*Power(xi, 41LL)*(510LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                                                             935550LL*Power(xi, 2LL)*Power(xj, 30LL)*(-69LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                                                             467775LL*r*Power(xi, 3LL)*Power(xj, 30LL)*

                                                                             (-253LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                                                             132LL*Power(r, 8LL)*Power(xi, 40LL)*(9180LL + 89LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                                                             311850LL*Power(xi, 4LL)*Power(xj, 28LL)*

                                                                             (2277LL - 345LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                                             31185LL*r*Power(xi, 5LL)*Power(xj, 28LL)*

                                                                             (41745LL - 2070LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

                                                                             + 1980LL*Power(r, 6LL)*Power(xi, 38LL)*

                                                                             (-162792LL - 11859LL*Power(r, 2LL)*Power(xj, 2LL) +

            41LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                                             22LL*Power(r, 7LL)*Power(xi, 39LL)*

                                                                             (-1046520LL - 30885LL*Power(r, 2LL)*Power(xj, 2LL) +

            44LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                                                             62370LL*Power(xi, 6LL)*Power(xj, 26LL)*

                                                                             (-79695LL + 18975LL*Power(r, 2LL)*Power(xj, 2LL) -

            460LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                                                             110LL*Power(r, 5LL)*Power(xi, 37LL)*

                                                                             (30767688LL + 4989438LL*Power(r, 2LL)*Power(xj, 2LL) -

            25359LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                                                             1485LL*r*Power(xi, 7LL)*Power(xj, 26LL)*

                                                                             (-6136515LL + 478170LL*Power(r, 2LL)*Power(xj, 2LL) -

            6762LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                                                             132LL*Power(r, 4LL)*Power(xi, 36LL)*

                                                                             (201455100LL + 69647445LL*Power(r, 2LL)*Power(xj, 2LL) -

            318735LL*Power(r, 4LL)*Power(xj, 4LL) + 353LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                                                             495LL*r*Power(xi, 9LL)*Power(xj, 24LL)*

                                                                             (92047725LL - 10041570LL*Power(r, 2LL)*Power(xj, 2LL) +

            223146LL*Power(r, 4LL)*Power(xj, 4LL) - 1380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                                                             2970LL*Power(xi, 8LL)*Power(xj, 24LL)*

                                                                             (8367975LL - 2789325LL*Power(r, 2LL)*Power(xj, 2LL) +

            106260LL*Power(r, 4LL)*Power(xj, 4LL) - 966LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                                                             22LL*Power(r, 3LL)*Power(xi, 35LL)*

                                                                             (6950200950LL + 5142653145LL*Power(r, 2LL)*Power(xj, 2LL) +

            7644510LL*Power(r, 4LL)*Power(xj, 4LL) -

            235635LL*Power(r, 6LL)*Power(xj, 6LL) + 124LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                                                             132LL*Power(r, 2LL)*Power(xi, 34LL)*

                                                                             (4633467300LL + 7767871650LL*Power(r, 2LL)*Power(xj, 2LL) +

            160904205LL*Power(r, 4LL)*Power(xj, 4LL) -

            2493315LL*Power(r, 6LL)*Power(xj, 6LL) + 5281LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                                                             495LL*r*Power(xi, 27LL)*Power(xj, 6LL)*

                                                                             (8395934795325LL - 439434024750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11948496210LL*Power(r, 4LL)*Power(xj, 4LL) -

            118623972LL*Power(r, 6LL)*Power(xj, 6LL) +

            248906LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                                                             11LL*r*Power(xi, 15LL)*Power(xj, 18LL)*

                                                                             (1488922594425LL + 252796524750LL*Power(r, 2LL)*Power(xj, 2LL) +

            6172031250LL*Power(r, 4LL)*Power(xj, 4LL) +

            104343660LL*Power(r, 6LL)*Power(xj, 6LL) +

            66810LL*Power(r, 8LL)*Power(xj, 8LL) - 88LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             66LL*Power(xi, 10LL)*Power(xj, 22LL)*

                                                                             (-1430923725LL + 627598125LL*Power(r, 2LL)*Power(xj, 2LL) -

            33471900LL*Power(r, 4LL)*Power(xj, 4LL) +

            478170LL*Power(r, 6LL)*Power(xj, 6LL) - 2070LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 10LL)*Power(xj, 10LL)) -

                                                                             1518LL*Power(xi, 12LL)*Power(xj, 20LL)*

                                                                             (-186642225LL + 103690125LL*Power(r, 2LL)*Power(xj, 2LL) -

            7276500LL*Power(r, 4LL)*Power(xj, 4LL) +

            145530LL*Power(r, 6LL)*Power(xj, 6LL) - 990LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             3LL*r*Power(xi, 11LL)*Power(xj, 22LL)*

                                                                             (-57713923575LL + 8284295250LL*Power(r, 2LL)*Power(xj, 2LL) -

            257733630LL*Power(r, 4LL)*Power(xj, 4LL) +

            2504700LL*Power(r, 6LL)*Power(xj, 6LL) - 7590LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             11LL*r*Power(xi, 13LL)*Power(xj, 20LL)*

                                                                             (56066193225LL - 6918959250LL*Power(r, 2LL)*Power(xj, 2LL) +

            430816050LL*Power(r, 4LL)*Power(xj, 4LL) -

            3349620LL*Power(r, 6LL)*Power(xj, 6LL) +

            33690LL*Power(r, 8LL)*Power(xj, 8LL) + 8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             55LL*r*Power(xi, 17LL)*Power(xj, 16LL)*

                                                                             (7416068831325LL + 658162968750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11421785970LL*Power(r, 4LL)*Power(xj, 4LL) -

            22800852LL*Power(r, 6LL)*Power(xj, 6LL) -

            224214LL*Power(r, 8LL)*Power(xj, 8LL) + 40LL*Power(r, 10LL)*Power(xj, 10LL)) -

                                                                             198LL*Power(xi, 14LL)*Power(xj, 18LL)*

                                                                             (12601626975LL + 2529410625LL*Power(r, 2LL)*Power(xj, 2LL) +

            582340500LL*Power(r, 4LL)*Power(xj, 4LL) +

            3239250LL*Power(r, 6LL)*Power(xj, 6LL) +

            132690LL*Power(r, 8LL)*Power(xj, 8LL) + 74LL*Power(r, 10LL)*Power(xj, 10LL)) -

                                                                             231LL*r*Power(xi, 25LL)*Power(xj, 8LL)*

                                                                             (21444497452125LL - 909858116250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1447333650LL*Power(r, 4LL)*Power(xj, 4LL) +

            178686540LL*Power(r, 6LL)*Power(xj, 6LL) -

            747270LL*Power(r, 8LL)*Power(xj, 8LL) + 184LL*Power(r, 10LL)*Power(xj, 10LL)) -

                                                                             198LL*Power(xi, 20LL)*Power(xj, 12LL)*

                                                                             (42449899182075LL + 4344172457625LL*Power(r, 2LL)*Power(xj, 2LL) -

            85249741500LL*Power(r, 4LL)*Power(xj, 4LL) -

            1059301110LL*Power(r, 6LL)*Power(xj, 6LL) +

            6582370LL*Power(r, 8LL)*Power(xj, 8LL) + 194LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             + 11LL*r*Power(xi, 19LL)*Power(xj, 14LL)*

                                                                             (239338679943825LL + 8851966719750LL*Power(r, 2LL)*Power(xj, 2LL) -

            112537092150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1100275380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2919090LL*Power(r, 8LL)*Power(xj, 8LL) + 248LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             - 330LL*Power(xi, 28LL)*Power(xj, 4LL)*

                                                                             (4860066085875LL + 2524912849305LL*Power(r, 2LL)*Power(xj, 2LL) -

            109538431380LL*Power(r, 4LL)*Power(xj, 4LL) +

            1633704282LL*Power(r, 6LL)*Power(xj, 6LL) -

            6421278LL*Power(r, 8LL)*Power(xj, 8LL) + 322LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             + 33LL*r*Power(xi, 29LL)*Power(xj, 4LL)*

                                                                             (-31641507079875LL - 2157639318450LL*Power(r, 2LL)*Power(xj, 2LL) +

            74910015810LL*Power(r, 4LL)*Power(xj, 4LL) -

            522003060LL*Power(r, 6LL)*Power(xj, 6LL) -

            250470LL*Power(r, 8LL)*Power(xj, 8LL) + 1288LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             - 330LL*Power(xi, 18LL)*Power(xj, 14LL)*

                                                                             (4867016286825LL + 1199363925375LL*Power(r, 2LL)*Power(xj, 2LL) +

            26817947100LL*Power(r, 4LL)*Power(xj, 4LL) -

            167333418LL*Power(r, 6LL)*Power(xj, 6LL) -

            1476138LL*Power(r, 8LL)*Power(xj, 8LL) + 1294LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             + 66LL*Power(xi, 16LL)*Power(xj, 16LL)*

                                                                             (-1657759205025LL - 682207855875LL*Power(r, 2LL)*Power(xj, 2LL) -

            31509229500LL*Power(r, 4LL)*Power(xj, 4LL) -

            492146550LL*Power(r, 6LL)*Power(xj, 6LL) -

            11910LL*Power(r, 8LL)*Power(xj, 8LL) + 2594LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             1386LL*Power(xi, 26LL)*Power(xj, 6LL)*

                                                                             (-6066588045375LL + 98854491375LL*Power(r, 2LL)*Power(xj, 2LL) -

            12496954500LL*Power(r, 4LL)*Power(xj, 4LL) +

            420813750LL*Power(r, 6LL)*Power(xj, 6LL) -

            2881210LL*Power(r, 8LL)*Power(xj, 8LL) + 2622LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             + 11LL*r*Power(xi, 23LL)*Power(xj, 10LL)*

                                                                             (149900659402725LL - 26541339882750LL*Power(r, 2LL)*Power(xj, 2LL) +

            594745455150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1399125420LL*Power(r, 6LL)*Power(xj, 6LL) -

            7887390LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             - 11LL*r*Power(xi, 31LL)*Power(xj, 2LL)*

                                                                             (7685082491625LL + 5034333946950LL*Power(r, 2LL)*Power(xj, 2LL) -

            108088893990LL*Power(r, 4LL)*Power(xj, 4LL) +

            1254174300LL*Power(r, 6LL)*Power(xj, 6LL) -

            6355950LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             - 462LL*Power(xi, 24LL)*Power(xj, 8LL)*

                                                                             (40495013164125LL - 3973079865375LL*Power(r, 2LL)*Power(xj, 2LL) +

            110288047500LL*Power(r, 4LL)*Power(xj, 4LL) -

            381623130LL*Power(r, 6LL)*Power(xj, 6LL) -

            4811370LL*Power(r, 8LL)*Power(xj, 8LL) + 9338LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                             + 198LL*Power(xi, 32LL)*(-9126526500LL - 152866565775LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                                                      32383266300LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                                                      709444890LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                                                      5562070LL*Power(r, 8LL)*Power(xj, 8LL) + 11042LL*Power(r, 10LL)*Power(xj, 10LL)

                                                                                                      ) + 2LL*r*Power(xi, 33LL)*(-764522104500LL -

                                                                                                                                 3357151476525LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                                                                                                 242177564475LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                                                                                 4513719870LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                                                                                 20531775LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                                                                                 11236LL*Power(r, 10LL)*Power(xj, 10LL)) -

                                                                             r*Power(xi, 21LL)*Power(xj, 12LL)*

                                                                             (-5533525427435775LL + 138591131159250LL*Power(r, 2LL)*Power(xj, 2LL) +

            2815739907750LL*Power(r, 4LL)*Power(xj, 4LL) -

            32922004500LL*Power(r, 6LL)*Power(xj, 6LL) +

            11347050LL*Power(r, 8LL)*Power(xj, 8LL) +

            22472LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             66LL*Power(xi, 22LL)*Power(xj, 10LL)*

                                                                             (-283522589265825LL + 7639225988625LL*Power(r, 2LL)*Power(xj, 2LL) +

            480728209500LL*Power(r, 4LL)*Power(xj, 4LL) -

            8458349130LL*Power(r, 6LL)*Power(xj, 6LL) +

            9771090LL*Power(r, 8LL)*Power(xj, 8LL) + 31786LL*Power(r, 10LL)*Power(xj, 10LL)

                                                                             ) - 66LL*Power(xi, 30LL)*Power(xj, 2LL)*

                                                                             (1678609807875LL + 4713298976925LL*Power(r, 2LL)*Power(xj, 2LL) -

            30578971500LL*Power(r, 4LL)*Power(xj, 4LL) +

            53723250LL*Power(r, 6LL)*Power(xj, 6LL) -

            9140190LL*Power(r, 8LL)*Power(xj, 8LL) + 38042LL*Power(r, 10LL)*Power(xj, 10LL))

                                                                                   ) + exp(2LL*r*xi)*Power(xi, 14LL)*

                                                                                   (-302841LL*Power(xi, 16LL)*Power(xj, 16LL)*

                                                                             (-361285650LL + 1346857875LL*r*xj -

            1306923750LL*Power(r, 2LL)*Power(xj, 2LL) +

            321527250LL*Power(r, 3LL)*Power(xj, 3LL) +

            55737000LL*Power(r, 4LL)*Power(xj, 4LL) -

            9297750LL*Power(r, 5LL)*Power(xj, 5LL) -

            1843380LL*Power(r, 6LL)*Power(xj, 6LL) - 50820LL*Power(r, 7LL)*Power(xj, 7LL) +

            7340LL*Power(r, 8LL)*Power(xj, 8LL) + 570LL*Power(r, 9LL)*Power(xj, 9LL) +

            12LL*Power(r, 10LL)*Power(xj, 10LL)) +

                                                                             12LL*Power(xj, 32LL)*(150587687250LL + 127420350750LL*r*xj +

                                                                                                   50968140300LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                                                   12742035075LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                                                                   2216006100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                                                   282037140LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                                                                   26860680LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                                                   1918620LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                                                   100980LL*Power(r, 8LL)*Power(xj, 8LL) + 3740LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                                                   88LL*Power(r, 10LL)*Power(xj, 10LL) + Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             3LL*Power(xi, 32LL)*(935550LL + 1715175LL*r*xj +

                                                                                                  1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                                                  935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                                                                  145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                                                  9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                                                  330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                                                                  4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             11LL*Power(xi, 30LL)*Power(xj, 2LL)*

                                                                             (-5868450LL - 10758825LL*r*xj - 9780750LL*Power(r, 2LL)*Power(xj, 2LL) -

            5868450LL*Power(r, 3LL)*Power(xj, 3LL) -

            2608200LL*Power(r, 4LL)*Power(xj, 4LL) -

            912870LL*Power(r, 5LL)*Power(xj, 5LL) - 260820LL*Power(r, 6LL)*Power(xj, 6LL) -

            62100LL*Power(r, 7LL)*Power(xj, 7LL) - 12420LL*Power(r, 8LL)*Power(xj, 8LL) -

            2070LL*Power(r, 9LL)*Power(xj, 9LL) - 276LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             5313LL*Power(xi, 14LL)*Power(xj, 18LL)*

                                                                             (-302299148250LL + 495525217275LL*r*xj -

            161894625750LL*Power(r, 2LL)*Power(xj, 2LL) -

            26085287250LL*Power(r, 3LL)*Power(xj, 3LL) +

            5971779000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1231357050LL*Power(r, 5LL)*Power(xj, 5LL) +

            33184620LL*Power(r, 6LL)*Power(xj, 6LL) -

            7768980LL*Power(r, 7LL)*Power(xj, 7LL) -

            751620LL*Power(r, 8LL)*Power(xj, 8LL) - 23190LL*Power(r, 9LL)*Power(xj, 9LL) -

            20LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             5313LL*Power(xi, 18LL)*Power(xj, 14LL)*

                                                                             (469625850LL - 3082655475LL*r*xj +

            8474631750LL*Power(r, 2LL)*Power(xj, 2LL) -

            6813281250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1665711000LL*Power(r, 4LL)*Power(xj, 4LL) +

            232996050LL*Power(r, 5LL)*Power(xj, 5LL) -

            39477060LL*Power(r, 6LL)*Power(xj, 6LL) -

            6196500LL*Power(r, 7LL)*Power(xj, 7LL) -

            121380LL*Power(r, 8LL)*Power(xj, 8LL) + 16330LL*Power(r, 9LL)*Power(xj, 9LL) +

            812LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             11LL*Power(xi, 2LL)*Power(xj, 30LL)*

                                                                             (10071658847250LL + 7685082491625LL*r*xj +

            2751598183950LL*Power(r, 2LL)*Power(xj, 2LL) +

            610391177550LL*Power(r, 3LL)*Power(xj, 3LL) +

            93214459800LL*Power(r, 4LL)*Power(xj, 4LL) +

            10285306290LL*Power(r, 5LL)*Power(xj, 5LL) +

            835769340LL*Power(r, 6LL)*Power(xj, 6LL) +

            49894380LL*Power(r, 7LL)*Power(xj, 7LL) +

            2134620LL*Power(r, 8LL)*Power(xj, 8LL) + 61770LL*Power(r, 9LL)*Power(xj, 9LL) +

            1068LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             11LL*Power(xi, 28LL)*Power(xj, 4LL)*

                                                                             (-64552950LL - 118347075LL*r*xj - 107588250LL*Power(r, 2LL)*Power(xj, 2LL) -

            64552950LL*Power(r, 3LL)*Power(xj, 3LL) -

            28690200LL*Power(r, 4LL)*Power(xj, 4LL) -

            10041570LL*Power(r, 5LL)*Power(xj, 5LL) -

            2869020LL*Power(r, 6LL)*Power(xj, 6LL) -

            683100LL*Power(r, 7LL)*Power(xj, 7LL) - 136620LL*Power(r, 8LL)*Power(xj, 8LL) -

            33690LL*Power(r, 9LL)*Power(xj, 9LL) + 1332LL*Power(r, 10LL)*Power(xj, 10LL) +

            88LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             11LL*Power(xi, 4LL)*Power(xj, 28LL)*

                                                                             (-145801982576250LL - 94924521239625LL*r*xj -

            28279793861550LL*Power(r, 2LL)*Power(xj, 2LL) -

            5034333946950LL*Power(r, 3LL)*Power(xj, 3LL) -

            582898793400LL*Power(r, 4LL)*Power(xj, 4LL) -

            44032284450LL*Power(r, 5LL)*Power(xj, 5LL) -

            1930850460LL*Power(r, 6LL)*Power(xj, 6LL) -

            15289020LL*Power(r, 7LL)*Power(xj, 7LL) +

            3824820LL*Power(r, 8LL)*Power(xj, 8LL) +

            253590LL*Power(r, 9LL)*Power(xj, 9LL) + 7380LL*Power(r, 10LL)*Power(xj, 10LL) +

            88LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             253LL*Power(xi, 20LL)*Power(xj, 12LL)*

                                                                             (1119853350LL + 2437660575LL*r*xj -

            1979538750LL*Power(r, 2LL)*Power(xj, 2LL) +

            10991153250LL*Power(r, 3LL)*Power(xj, 3LL) -

            8219799000LL*Power(r, 4LL)*Power(xj, 4LL) +

            2482996950LL*Power(r, 5LL)*Power(xj, 5LL) +

            218260980LL*Power(r, 6LL)*Power(xj, 6LL) -

            47838060LL*Power(r, 7LL)*Power(xj, 7LL) -

            5151420LL*Power(r, 8LL)*Power(xj, 8LL) - 44850LL*Power(r, 9LL)*Power(xj, 9LL) +

            8292LL*Power(r, 10LL)*Power(xj, 10LL) + 184LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             253LL*Power(xi, 12LL)*Power(xj, 20LL)*

                                                                             (33221660229450LL - 21871642005675LL*r*xj -

            1992841562250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1153971299250LL*Power(r, 3LL)*Power(xj, 3LL) +

            201395565000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1321478550LL*Power(r, 5LL)*Power(xj, 5LL) -

            2305327500LL*Power(r, 6LL)*Power(xj, 6LL) -

            232090380LL*Power(r, 7LL)*Power(xj, 7LL) -

            8375580LL*Power(r, 8LL)*Power(xj, 8LL) + 32670LL*Power(r, 9LL)*Power(xj, 9LL) +

            9924LL*Power(r, 10LL)*Power(xj, 10LL) + 184LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             11LL*Power(xi, 6LL)*Power(xj, 26LL)*

                                                                             (764390093717250LL + 377817065789625LL*r*xj +

            75747385479150LL*Power(r, 2LL)*Power(xj, 2LL) +

            6472917955350LL*Power(r, 3LL)*Power(xj, 3LL) -

            183473829000LL*Power(r, 4LL)*Power(xj, 4LL) -

            108088893990LL*Power(r, 5LL)*Power(xj, 5LL) -

            12770008020LL*Power(r, 6LL)*Power(xj, 6LL) -

            820676340LL*Power(r, 7LL)*Power(xj, 7LL) -

            29919780LL*Power(r, 8LL)*Power(xj, 8LL) -

            471270LL*Power(r, 9LL)*Power(xj, 9LL) + 4236LL*Power(r, 10LL)*Power(xj, 10LL) +

            200LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             11LL*Power(xi, 26LL)*Power(xj, 6LL)*

                                                                             (-451870650LL - 828429525LL*r*xj - 753117750LL*Power(r, 2LL)*Power(xj, 2LL) -

            451870650LL*Power(r, 3LL)*Power(xj, 3LL) -

            200831400LL*Power(r, 4LL)*Power(xj, 4LL) -

            70290990LL*Power(r, 5LL)*Power(xj, 5LL) -

            20083140LL*Power(r, 6LL)*Power(xj, 6LL) -

            3349620LL*Power(r, 7LL)*Power(xj, 7LL) -

            2388420LL*Power(r, 8LL)*Power(xj, 8LL) + 66810LL*Power(r, 9LL)*Power(xj, 9LL) +

            15564LL*Power(r, 10LL)*Power(xj, 10LL) + 200LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                                                             11LL*Power(xi, 24LL)*Power(xj, 8LL)*

                                                                             (2259353250LL + 4142147625LL*r*xj +

            3765588750LL*Power(r, 2LL)*Power(xj, 2LL) +

            2259353250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1004157000LL*Power(r, 4LL)*Power(xj, 4LL) +

            430816050LL*Power(r, 5LL)*Power(xj, 5LL) -

            58306500LL*Power(r, 6LL)*Power(xj, 6LL) +

            104343660LL*Power(r, 7LL)*Power(xj, 7LL) -

            71460LL*Power(r, 8LL)*Power(xj, 8LL) - 1121070LL*Power(r, 9LL)*Power(xj, 9LL) -

            38820LL*Power(r, 10LL)*Power(xj, 10LL) + 248LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             11LL*Power(xi, 8LL)*Power(xj, 24LL)*

                                                                             (1700790552893250LL + 450334446494625LL*r*xj -

            12455665913250LL*Power(r, 2LL)*Power(xj, 2LL) -

            19774531113750LL*Power(r, 3LL)*Power(xj, 3LL) -

            3286152941400LL*Power(r, 4LL)*Power(xj, 4LL) -

            224730047430LL*Power(r, 5LL)*Power(xj, 5LL) +

            322339500LL*Power(r, 6LL)*Power(xj, 6LL) +

            1254174300LL*Power(r, 7LL)*Power(xj, 7LL) +

            100117260LL*Power(r, 8LL)*Power(xj, 8LL) +

            3733050LL*Power(r, 9LL)*Power(xj, 9LL) +

            63372LL*Power(r, 10LL)*Power(xj, 10LL) + 248LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                                                             Power(xi, 22LL)*Power(xj, 10LL)*

                                                                             (94440965850LL + 173141770725LL*r*xj +

            157401609750LL*Power(r, 2LL)*Power(xj, 2LL) +

            76108551750LL*Power(r, 3LL)*Power(xj, 3LL) +

            115303419000LL*Power(r, 4LL)*Power(xj, 4LL) -

            67892343750LL*Power(r, 5LL)*Power(xj, 5LL) +

            32481672300LL*Power(r, 6LL)*Power(xj, 6LL) +

            1254046860LL*Power(r, 7LL)*Power(xj, 7LL) -

            487125540LL*Power(r, 8LL)*Power(xj, 8LL) -

            32109990LL*Power(r, 9LL)*Power(xj, 9LL) +

            38412LL*Power(r, 10LL)*Power(xj, 10LL) + 22472LL*Power(r, 11LL)*Power(xj, 11LL))

                                                                             - Power(xi, 10LL)*Power(xj, 22LL)*(-18712490891544450LL + 1648907253429975LL*r*xj +

                                                                                                                1835562897803250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                                                                210177224853750LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                                                                                17320778937000LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                                                                5914505623950LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                                                                                539122413060LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                                                                                17226100980LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                                                                                603252540LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                                                                                69915450LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                                                                                2186316LL*Power(r, 10LL)*Power(xj, 10LL) + 22472LL*Power(r, 11LL)*Power(xj, 11LL)

                                                                                                                )))/(1.403325e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 23LL)*Power(xi + xj, 22LL))

                - (5613300LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                   Power(Power(xi, 2LL) - Power(xj, 2LL), 23LL) +

                   exp(2LL*r*xj)*Power(xj, 14LL)*

                   (-10560LL*Power(r, 9LL)*Power(xi, 42LL) - 132LL*Power(r, 10LL)*Power(xi, 43LL) -

                    23496LL*Power(r, 9LL)*Power(xi, 40LL)*Power(xj, 2LL) -

                    176LL*Power(r, 10LL)*Power(xi, 41LL)*Power(xj, 2LL) +

                    5145525LL*xi*Power(xj, 32LL) + 9355500LL*r*Power(xi, 2LL)*Power(xj, 32LL) +

                    5613300LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 32LL) -

                    792LL*Power(r, 8LL)*Power(xi, 41LL)*(510LL + Power(r, 2LL)*Power(xj, 2LL)) +

                    467775LL*Power(xi, 3LL)*Power(xj, 30LL)*(-253LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    1056LL*Power(r, 7LL)*Power(xi, 40LL)*(9180LL + 89LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    311850LL*Power(xi, 4LL)*Power(xj, 28LL)*

                    (-690LL*r*Power(xj, 2LL) + 16LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    31185LL*r*Power(xi, 5LL)*Power(xj, 28LL)*

                    (-4140LL*r*Power(xj, 2LL) + 56LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    1980LL*Power(r, 6LL)*Power(xi, 38LL)*

                    (-23718LL*r*Power(xj, 2LL) + 164LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    22LL*Power(r, 7LL)*Power(xi, 39LL)*

                    (-61770LL*r*Power(xj, 2LL) + 176LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    31185LL*Power(xi, 5LL)*Power(xj, 28LL)*

                    (41745LL - 2070LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

                    + 11880LL*Power(r, 5LL)*Power(xi, 38LL)*

                    (-162792LL - 11859LL*Power(r, 2LL)*Power(xj, 2LL) +

            41LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    154LL*Power(r, 6LL)*Power(xi, 39LL)*

                    (-1046520LL - 30885LL*Power(r, 2LL)*Power(xj, 2LL) +

            44LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    62370LL*Power(xi, 6LL)*Power(xj, 26LL)*

                    (37950LL*r*Power(xj, 2LL) - 1840LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) -

                    110LL*Power(r, 5LL)*Power(xi, 37LL)*

                    (9978876LL*r*Power(xj, 2LL) - 101436LL*Power(r, 3LL)*Power(xj, 4LL) +

            120LL*Power(r, 5LL)*Power(xj, 6LL)) +

                    1485LL*r*Power(xi, 7LL)*Power(xj, 26LL)*

                    (956340LL*r*Power(xj, 2LL) - 27048LL*Power(r, 3LL)*Power(xj, 4LL) +

            120LL*Power(r, 5LL)*Power(xj, 6LL)) -

                    132LL*Power(r, 4LL)*Power(xi, 36LL)*

                    (139294890LL*r*Power(xj, 2LL) - 1274940LL*Power(r, 3LL)*Power(xj, 4LL) +

            2118LL*Power(r, 5LL)*Power(xj, 6LL)) -

                    550LL*Power(r, 4LL)*Power(xi, 37LL)*

                    (30767688LL + 4989438LL*Power(r, 2LL)*Power(xj, 2LL) -

            25359LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    1485LL*Power(xi, 7LL)*Power(xj, 26LL)*

                    (-6136515LL + 478170LL*Power(r, 2LL)*Power(xj, 2LL) -

            6762LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

                    528LL*Power(r, 3LL)*Power(xi, 36LL)*

                    (201455100LL + 69647445LL*Power(r, 2LL)*Power(xj, 2LL) -

            318735LL*Power(r, 4LL)*Power(xj, 4LL) + 353LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    495LL*r*Power(xi, 9LL)*Power(xj, 24LL)*

                    (-20083140LL*r*Power(xj, 2LL) + 892584LL*Power(r, 3LL)*Power(xj, 4LL) -

            8280LL*Power(r, 5LL)*Power(xj, 6LL) + 16LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    2970LL*Power(xi, 8LL)*Power(xj, 24LL)*

                    (-5578650LL*r*Power(xj, 2LL) + 425040LL*Power(r, 3LL)*Power(xj, 4LL) -

            5796LL*Power(r, 5LL)*Power(xj, 6LL) + 16LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    22LL*Power(r, 3LL)*Power(xi, 35LL)*

                    (10285306290LL*r*Power(xj, 2LL) + 30578040LL*Power(r, 3LL)*Power(xj, 4LL) -

            1413810LL*Power(r, 5LL)*Power(xj, 6LL) + 992LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    132LL*Power(r, 2LL)*Power(xi, 34LL)*

                    (15535743300LL*r*Power(xj, 2LL) + 643616820LL*Power(r, 3LL)*Power(xj, 4LL) -

            14959890LL*Power(r, 5LL)*Power(xj, 6LL) + 42248LL*Power(r, 7LL)*Power(xj, 8LL))

                    - 495LL*r*Power(xi, 27LL)*Power(xj, 6LL)*

                    (-878868049500LL*r*Power(xj, 2LL) +

            47793984840LL*Power(r, 3LL)*Power(xj, 4LL) -

            711743832LL*Power(r, 5LL)*Power(xj, 6LL) +

            1991248LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    495LL*Power(xi, 9LL)*Power(xj, 24LL)*

                    (92047725LL - 10041570LL*Power(r, 2LL)*Power(xj, 2LL) +

            223146LL*Power(r, 4LL)*Power(xj, 4LL) - 1380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    66LL*Power(r, 2LL)*Power(xi, 35LL)*

                    (6950200950LL + 5142653145LL*Power(r, 2LL)*Power(xj, 2LL) +

            7644510LL*Power(r, 4LL)*Power(xj, 4LL) -

            235635LL*Power(r, 6LL)*Power(xj, 6LL) + 124LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    264LL*r*Power(xi, 34LL)*(4633467300LL +

                                             7767871650LL*Power(r, 2LL)*Power(xj, 2LL) +

                                             160904205LL*Power(r, 4LL)*Power(xj, 4LL) -

                                             2493315LL*Power(r, 6LL)*Power(xj, 6LL) + 5281LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    495LL*Power(xi, 27LL)*Power(xj, 6LL)*

                    (8395934795325LL - 439434024750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11948496210LL*Power(r, 4LL)*Power(xj, 4LL) -

            118623972LL*Power(r, 6LL)*Power(xj, 6LL) + 248906LL*Power(r, 8LL)*Power(xj, 8LL)

                    ) + 11LL*r*Power(xi, 15LL)*Power(xj, 18LL)*

                    (505593049500LL*r*Power(xj, 2LL) +

            24688125000LL*Power(r, 3LL)*Power(xj, 4LL) +

            626061960LL*Power(r, 5LL)*Power(xj, 6LL) +

            534480LL*Power(r, 7LL)*Power(xj, 8LL) - 880LL*Power(r, 9LL)*Power(xj, 10LL)) +

                    66LL*Power(xi, 10LL)*Power(xj, 22LL)*

                    (1255196250LL*r*Power(xj, 2LL) - 133887600LL*Power(r, 3LL)*Power(xj, 4LL) +

            2869020LL*Power(r, 5LL)*Power(xj, 6LL) - 16560LL*Power(r, 7LL)*Power(xj, 8LL) +

            20LL*Power(r, 9LL)*Power(xj, 10LL)) -

                    1518LL*Power(xi, 12LL)*Power(xj, 20LL)*

                    (207380250LL*r*Power(xj, 2LL) - 29106000LL*Power(r, 3LL)*Power(xj, 4LL) +

            873180LL*Power(r, 5LL)*Power(xj, 6LL) - 7920LL*Power(r, 7LL)*Power(xj, 8LL) +

            20LL*Power(r, 9LL)*Power(xj, 10LL)) +

                    3LL*r*Power(xi, 11LL)*Power(xj, 22LL)*

                    (16568590500LL*r*Power(xj, 2LL) - 1030934520LL*Power(r, 3LL)*Power(xj, 4LL) +

            15028200LL*Power(r, 5LL)*Power(xj, 6LL) -

            60720LL*Power(r, 7LL)*Power(xj, 8LL) + 40LL*Power(r, 9LL)*Power(xj, 10LL)) +

                    11LL*r*Power(xi, 13LL)*Power(xj, 20LL)*

                    (-13837918500LL*r*Power(xj, 2LL) + 1723264200LL*Power(r, 3LL)*Power(xj, 4LL) -

            20097720LL*Power(r, 5LL)*Power(xj, 6LL) +

            269520LL*Power(r, 7LL)*Power(xj, 8LL) + 80LL*Power(r, 9LL)*Power(xj, 10LL)) +

                    55LL*r*Power(xi, 17LL)*Power(xj, 16LL)*

                    (1316325937500LL*r*Power(xj, 2LL) +

            45687143880LL*Power(r, 3LL)*Power(xj, 4LL) -

            136805112LL*Power(r, 5LL)*Power(xj, 6LL) -

            1793712LL*Power(r, 7LL)*Power(xj, 8LL) + 400LL*Power(r, 9LL)*Power(xj, 10LL)) -

                    198LL*Power(xi, 14LL)*Power(xj, 18LL)*

                    (5058821250LL*r*Power(xj, 2LL) + 2329362000LL*Power(r, 3LL)*Power(xj, 4LL) +

            19435500LL*Power(r, 5LL)*Power(xj, 6LL) +

            1061520LL*Power(r, 7LL)*Power(xj, 8LL) + 740LL*Power(r, 9LL)*Power(xj, 10LL)) -

                    231LL*r*Power(xi, 25LL)*Power(xj, 8LL)*

                    (-1819716232500LL*r*Power(xj, 2LL) +

            5789334600LL*Power(r, 3LL)*Power(xj, 4LL) +

            1072119240LL*Power(r, 5LL)*Power(xj, 6LL) -

            5978160LL*Power(r, 7LL)*Power(xj, 8LL) + 1840LL*Power(r, 9LL)*Power(xj, 10LL)) -

                    198LL*Power(xi, 20LL)*Power(xj, 12LL)*

                    (8688344915250LL*r*Power(xj, 2LL) -

            340998966000LL*Power(r, 3LL)*Power(xj, 4LL) -

            6355806660LL*Power(r, 5LL)*Power(xj, 6LL) +

            52658960LL*Power(r, 7LL)*Power(xj, 8LL) + 1940LL*Power(r, 9LL)*Power(xj, 10LL))

                    + 11LL*r*Power(xi, 19LL)*Power(xj, 14LL)*

                    (17703933439500LL*r*Power(xj, 2LL) -

            450148368600LL*Power(r, 3LL)*Power(xj, 4LL) -

            6601652280LL*Power(r, 5LL)*Power(xj, 6LL) +

            23352720LL*Power(r, 7LL)*Power(xj, 8LL) + 2480LL*Power(r, 9LL)*Power(xj, 10LL))

                    - 330LL*Power(xi, 28LL)*Power(xj, 4LL)*

                    (5049825698610LL*r*Power(xj, 2LL) -

            438153725520LL*Power(r, 3LL)*Power(xj, 4LL) +

            9802225692LL*Power(r, 5LL)*Power(xj, 6LL) -

            51370224LL*Power(r, 7LL)*Power(xj, 8LL) + 3220LL*Power(r, 9LL)*Power(xj, 10LL))

                    + 33LL*r*Power(xi, 29LL)*Power(xj, 4LL)*

                    (-4315278636900LL*r*Power(xj, 2LL) +

            299640063240LL*Power(r, 3LL)*Power(xj, 4LL) -

            3132018360LL*Power(r, 5LL)*Power(xj, 6LL) -

            2003760LL*Power(r, 7LL)*Power(xj, 8LL) + 12880LL*Power(r, 9LL)*Power(xj, 10LL))

                    - 330LL*Power(xi, 18LL)*Power(xj, 14LL)*

                    (2398727850750LL*r*Power(xj, 2LL) +

            107271788400LL*Power(r, 3LL)*Power(xj, 4LL) -

            1004000508LL*Power(r, 5LL)*Power(xj, 6LL) -

            11809104LL*Power(r, 7LL)*Power(xj, 8LL) + 12940LL*Power(r, 9LL)*Power(xj, 10LL))

                    + 66LL*Power(xi, 16LL)*Power(xj, 16LL)*

                    (-1364415711750LL*r*Power(xj, 2LL) -

            126036918000LL*Power(r, 3LL)*Power(xj, 4LL) -

            2952879300LL*Power(r, 5LL)*Power(xj, 6LL) -

            95280LL*Power(r, 7LL)*Power(xj, 8LL) + 25940LL*Power(r, 9LL)*Power(xj, 10LL)) +

                    1386LL*Power(xi, 26LL)*Power(xj, 6LL)*

                    (197708982750LL*r*Power(xj, 2LL) -

            49987818000LL*Power(r, 3LL)*Power(xj, 4LL) +

            2524882500LL*Power(r, 5LL)*Power(xj, 6LL) -

            23049680LL*Power(r, 7LL)*Power(xj, 8LL) + 26220LL*Power(r, 9LL)*Power(xj, 10LL))

                    + 11LL*r*Power(xi, 23LL)*Power(xj, 10LL)*

                    (-53082679765500LL*r*Power(xj, 2LL) +

            2378981820600LL*Power(r, 3LL)*Power(xj, 4LL) -

            8394752520LL*Power(r, 5LL)*Power(xj, 6LL) -

            63099120LL*Power(r, 7LL)*Power(xj, 8LL) + 42320LL*Power(r, 9LL)*Power(xj, 10LL))

                    - 11LL*r*Power(xi, 31LL)*Power(xj, 2LL)*

                    (10068667893900LL*r*Power(xj, 2LL) -

            432355575960LL*Power(r, 3LL)*Power(xj, 4LL) +

            7525045800LL*Power(r, 5LL)*Power(xj, 6LL) -

            50847600LL*Power(r, 7LL)*Power(xj, 8LL) + 42320LL*Power(r, 9LL)*Power(xj, 10LL))

                    - 462LL*Power(xi, 24LL)*Power(xj, 8LL)*

                    (-7946159730750LL*r*Power(xj, 2LL) +

            441152190000LL*Power(r, 3LL)*Power(xj, 4LL) -

            2289738780LL*Power(r, 5LL)*Power(xj, 6LL) -

            38490960LL*Power(r, 7LL)*Power(xj, 8LL) + 93380LL*Power(r, 9LL)*Power(xj, 10LL))

                    + 198LL*Power(xi, 32LL)*(-305733131550LL*r*Power(xj, 2LL) -

                                             129533065200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                             4256669340LL*Power(r, 5LL)*Power(xj, 6LL) -

                                             44496560LL*Power(r, 7LL)*Power(xj, 8LL) + 110420LL*Power(r, 9LL)*Power(xj, 10LL)

                                             ) + 2LL*r*Power(xi, 33LL)*(-6714302953050LL*r*Power(xj, 2LL) -

                                                                        968710257900LL*Power(r, 3LL)*Power(xj, 4LL) +

                                                                        27082319220LL*Power(r, 5LL)*Power(xj, 6LL) -

                                                                        164254200LL*Power(r, 7LL)*Power(xj, 8LL) +

                                                                        112360LL*Power(r, 9LL)*Power(xj, 10LL)) -

                    r*Power(xi, 21LL)*Power(xj, 12LL)*

                    (277182262318500LL*r*Power(xj, 2LL) +

            11262959631000LL*Power(r, 3LL)*Power(xj, 4LL) -

            197532027000LL*Power(r, 5LL)*Power(xj, 6LL) +

            90776400LL*Power(r, 7LL)*Power(xj, 8LL) + 224720LL*Power(r, 9LL)*Power(xj, 10LL)

                    ) + 66LL*Power(xi, 22LL)*Power(xj, 10LL)*

                    (15278451977250LL*r*Power(xj, 2LL) +

            1922912838000LL*Power(r, 3LL)*Power(xj, 4LL) -

            50750094780LL*Power(r, 5LL)*Power(xj, 6LL) +

            78168720LL*Power(r, 7LL)*Power(xj, 8LL) + 317860LL*Power(r, 9LL)*Power(xj, 10LL)

                    ) - 66LL*Power(xi, 30LL)*Power(xj, 2LL)*

                    (9426597953850LL*r*Power(xj, 2LL) -

            122315886000LL*Power(r, 3LL)*Power(xj, 4LL) +

            322339500LL*Power(r, 5LL)*Power(xj, 6LL) -

            73121520LL*Power(r, 7LL)*Power(xj, 8LL) + 380420LL*Power(r, 9LL)*Power(xj, 10LL)

                    ) + 11LL*Power(xi, 15LL)*Power(xj, 18LL)*

                    (1488922594425LL + 252796524750LL*Power(r, 2LL)*Power(xj, 2LL) +

            6172031250LL*Power(r, 4LL)*Power(xj, 4LL) +

            104343660LL*Power(r, 6LL)*Power(xj, 6LL) +

            66810LL*Power(r, 8LL)*Power(xj, 8LL) - 88LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    3LL*Power(xi, 11LL)*Power(xj, 22LL)*

                    (-57713923575LL + 8284295250LL*Power(r, 2LL)*Power(xj, 2LL) -

            257733630LL*Power(r, 4LL)*Power(xj, 4LL) +

            2504700LL*Power(r, 6LL)*Power(xj, 6LL) - 7590LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    11LL*Power(xi, 13LL)*Power(xj, 20LL)*

                    (56066193225LL - 6918959250LL*Power(r, 2LL)*Power(xj, 2LL) +

            430816050LL*Power(r, 4LL)*Power(xj, 4LL) -

            3349620LL*Power(r, 6LL)*Power(xj, 6LL) + 33690LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    55LL*Power(xi, 17LL)*Power(xj, 16LL)*

                    (7416068831325LL + 658162968750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11421785970LL*Power(r, 4LL)*Power(xj, 4LL) -

            22800852LL*Power(r, 6LL)*Power(xj, 6LL) -

            224214LL*Power(r, 8LL)*Power(xj, 8LL) + 40LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    231LL*Power(xi, 25LL)*Power(xj, 8LL)*

                    (21444497452125LL - 909858116250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1447333650LL*Power(r, 4LL)*Power(xj, 4LL) +

            178686540LL*Power(r, 6LL)*Power(xj, 6LL) -

            747270LL*Power(r, 8LL)*Power(xj, 8LL) + 184LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    11LL*Power(xi, 19LL)*Power(xj, 14LL)*

                    (239338679943825LL + 8851966719750LL*Power(r, 2LL)*Power(xj, 2LL) -

            112537092150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1100275380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2919090LL*Power(r, 8LL)*Power(xj, 8LL) + 248LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    33LL*Power(xi, 29LL)*Power(xj, 4LL)*

                    (-31641507079875LL - 2157639318450LL*Power(r, 2LL)*Power(xj, 2LL) +

            74910015810LL*Power(r, 4LL)*Power(xj, 4LL) -

            522003060LL*Power(r, 6LL)*Power(xj, 6LL) -

            250470LL*Power(r, 8LL)*Power(xj, 8LL) + 1288LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    11LL*Power(xi, 23LL)*Power(xj, 10LL)*

                    (149900659402725LL - 26541339882750LL*Power(r, 2LL)*Power(xj, 2LL) +

            594745455150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1399125420LL*Power(r, 6LL)*Power(xj, 6LL) -

            7887390LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                    - 11LL*Power(xi, 31LL)*Power(xj, 2LL)*(7685082491625LL +

                                                           5034333946950LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                           108088893990LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                           1254174300LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                           6355950LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                    + 2LL*Power(xi, 33LL)*(-764522104500LL - 3357151476525LL*Power(r, 2LL)*Power(xj, 2LL) -

                                           242177564475LL*Power(r, 4LL)*Power(xj, 4LL) +

                                           4513719870LL*Power(r, 6LL)*Power(xj, 6LL) -

                                           20531775LL*Power(r, 8LL)*Power(xj, 8LL) + 11236LL*Power(r, 10LL)*Power(xj, 10LL)

                                           ) - Power(xi, 21LL)*Power(xj, 12LL)*(-5533525427435775LL +

                                                                                138591131159250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                                                2815739907750LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                                                32922004500LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                                                11347050LL*Power(r, 8LL)*Power(xj, 8LL) + 22472LL*Power(r, 10LL)*Power(xj, 10LL))

                   ) + 2LL*exp(2LL*r*xj)*Power(xj, 15LL)*

                   (-1056LL*Power(r, 10LL)*Power(xi, 42LL) - 12LL*Power(r, 11LL)*Power(xi, 43LL) +

                    2806650LL*Power(xj, 32LL) + 5145525LL*r*xi*Power(xj, 32LL) -

                    88LL*Power(r, 9LL)*Power(xi, 41LL)*(510LL + Power(r, 2LL)*Power(xj, 2LL)) +

                    935550LL*Power(xi, 2LL)*Power(xj, 30LL)*(-69LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    467775LL*r*Power(xi, 3LL)*Power(xj, 30LL)*

                    (-253LL + 6LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    132LL*Power(r, 8LL)*Power(xi, 40LL)*(9180LL + 89LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    311850LL*Power(xi, 4LL)*Power(xj, 28LL)*

                    (2277LL - 345LL*Power(r, 2LL)*Power(xj, 2LL) + 4LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    31185LL*r*Power(xi, 5LL)*Power(xj, 28LL)*

                    (41745LL - 2070LL*Power(r, 2LL)*Power(xj, 2LL) + 14LL*Power(r, 4LL)*Power(xj, 4LL))

                    + 1980LL*Power(r, 6LL)*Power(xi, 38LL)*

                    (-162792LL - 11859LL*Power(r, 2LL)*Power(xj, 2LL) +

            41LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    22LL*Power(r, 7LL)*Power(xi, 39LL)*

                    (-1046520LL - 30885LL*Power(r, 2LL)*Power(xj, 2LL) +

            44LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    62370LL*Power(xi, 6LL)*Power(xj, 26LL)*

                    (-79695LL + 18975LL*Power(r, 2LL)*Power(xj, 2LL) -

            460LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                    110LL*Power(r, 5LL)*Power(xi, 37LL)*

                    (30767688LL + 4989438LL*Power(r, 2LL)*Power(xj, 2LL) -

            25359LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    1485LL*r*Power(xi, 7LL)*Power(xj, 26LL)*

                    (-6136515LL + 478170LL*Power(r, 2LL)*Power(xj, 2LL) -

            6762LL*Power(r, 4LL)*Power(xj, 4LL) + 20LL*Power(r, 6LL)*Power(xj, 6LL)) -

                    132LL*Power(r, 4LL)*Power(xi, 36LL)*

                    (201455100LL + 69647445LL*Power(r, 2LL)*Power(xj, 2LL) -

            318735LL*Power(r, 4LL)*Power(xj, 4LL) + 353LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    495LL*r*Power(xi, 9LL)*Power(xj, 24LL)*

                    (92047725LL - 10041570LL*Power(r, 2LL)*Power(xj, 2LL) +

            223146LL*Power(r, 4LL)*Power(xj, 4LL) - 1380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    2970LL*Power(xi, 8LL)*Power(xj, 24LL)*

                    (8367975LL - 2789325LL*Power(r, 2LL)*Power(xj, 2LL) +

            106260LL*Power(r, 4LL)*Power(xj, 4LL) - 966LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    22LL*Power(r, 3LL)*Power(xi, 35LL)*

                    (6950200950LL + 5142653145LL*Power(r, 2LL)*Power(xj, 2LL) +

            7644510LL*Power(r, 4LL)*Power(xj, 4LL) -

            235635LL*Power(r, 6LL)*Power(xj, 6LL) + 124LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    132LL*Power(r, 2LL)*Power(xi, 34LL)*

                    (4633467300LL + 7767871650LL*Power(r, 2LL)*Power(xj, 2LL) +

            160904205LL*Power(r, 4LL)*Power(xj, 4LL) -

            2493315LL*Power(r, 6LL)*Power(xj, 6LL) + 5281LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    495LL*r*Power(xi, 27LL)*Power(xj, 6LL)*

                    (8395934795325LL - 439434024750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11948496210LL*Power(r, 4LL)*Power(xj, 4LL) -

            118623972LL*Power(r, 6LL)*Power(xj, 6LL) + 248906LL*Power(r, 8LL)*Power(xj, 8LL)

                    ) + 11LL*r*Power(xi, 15LL)*Power(xj, 18LL)*

                    (1488922594425LL + 252796524750LL*Power(r, 2LL)*Power(xj, 2LL) +

            6172031250LL*Power(r, 4LL)*Power(xj, 4LL) +

            104343660LL*Power(r, 6LL)*Power(xj, 6LL) +

            66810LL*Power(r, 8LL)*Power(xj, 8LL) - 88LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    66LL*Power(xi, 10LL)*Power(xj, 22LL)*

                    (-1430923725LL + 627598125LL*Power(r, 2LL)*Power(xj, 2LL) -

            33471900LL*Power(r, 4LL)*Power(xj, 4LL) +

            478170LL*Power(r, 6LL)*Power(xj, 6LL) - 2070LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    1518LL*Power(xi, 12LL)*Power(xj, 20LL)*

                    (-186642225LL + 103690125LL*Power(r, 2LL)*Power(xj, 2LL) -

            7276500LL*Power(r, 4LL)*Power(xj, 4LL) +

            145530LL*Power(r, 6LL)*Power(xj, 6LL) - 990LL*Power(r, 8LL)*Power(xj, 8LL) +

            2LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    3LL*r*Power(xi, 11LL)*Power(xj, 22LL)*

                    (-57713923575LL + 8284295250LL*Power(r, 2LL)*Power(xj, 2LL) -

            257733630LL*Power(r, 4LL)*Power(xj, 4LL) +

            2504700LL*Power(r, 6LL)*Power(xj, 6LL) - 7590LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    11LL*r*Power(xi, 13LL)*Power(xj, 20LL)*

                    (56066193225LL - 6918959250LL*Power(r, 2LL)*Power(xj, 2LL) +

            430816050LL*Power(r, 4LL)*Power(xj, 4LL) -

            3349620LL*Power(r, 6LL)*Power(xj, 6LL) + 33690LL*Power(r, 8LL)*Power(xj, 8LL) +

            8LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    55LL*r*Power(xi, 17LL)*Power(xj, 16LL)*

                    (7416068831325LL + 658162968750LL*Power(r, 2LL)*Power(xj, 2LL) +

            11421785970LL*Power(r, 4LL)*Power(xj, 4LL) -

            22800852LL*Power(r, 6LL)*Power(xj, 6LL) -

            224214LL*Power(r, 8LL)*Power(xj, 8LL) + 40LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    198LL*Power(xi, 14LL)*Power(xj, 18LL)*

                    (12601626975LL + 2529410625LL*Power(r, 2LL)*Power(xj, 2LL) +

            582340500LL*Power(r, 4LL)*Power(xj, 4LL) +

            3239250LL*Power(r, 6LL)*Power(xj, 6LL) +

            132690LL*Power(r, 8LL)*Power(xj, 8LL) + 74LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    231LL*r*Power(xi, 25LL)*Power(xj, 8LL)*

                    (21444497452125LL - 909858116250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1447333650LL*Power(r, 4LL)*Power(xj, 4LL) +

            178686540LL*Power(r, 6LL)*Power(xj, 6LL) -

            747270LL*Power(r, 8LL)*Power(xj, 8LL) + 184LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    198LL*Power(xi, 20LL)*Power(xj, 12LL)*

                    (42449899182075LL + 4344172457625LL*Power(r, 2LL)*Power(xj, 2LL) -

            85249741500LL*Power(r, 4LL)*Power(xj, 4LL) -

            1059301110LL*Power(r, 6LL)*Power(xj, 6LL) +

            6582370LL*Power(r, 8LL)*Power(xj, 8LL) + 194LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    11LL*r*Power(xi, 19LL)*Power(xj, 14LL)*

                    (239338679943825LL + 8851966719750LL*Power(r, 2LL)*Power(xj, 2LL) -

            112537092150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1100275380LL*Power(r, 6LL)*Power(xj, 6LL) +

            2919090LL*Power(r, 8LL)*Power(xj, 8LL) + 248LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    330LL*Power(xi, 28LL)*Power(xj, 4LL)*

                    (4860066085875LL + 2524912849305LL*Power(r, 2LL)*Power(xj, 2LL) -

            109538431380LL*Power(r, 4LL)*Power(xj, 4LL) +

            1633704282LL*Power(r, 6LL)*Power(xj, 6LL) -

            6421278LL*Power(r, 8LL)*Power(xj, 8LL) + 322LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    33LL*r*Power(xi, 29LL)*Power(xj, 4LL)*

                    (-31641507079875LL - 2157639318450LL*Power(r, 2LL)*Power(xj, 2LL) +

            74910015810LL*Power(r, 4LL)*Power(xj, 4LL) -

            522003060LL*Power(r, 6LL)*Power(xj, 6LL) -

            250470LL*Power(r, 8LL)*Power(xj, 8LL) + 1288LL*Power(r, 10LL)*Power(xj, 10LL)) -

                    330LL*Power(xi, 18LL)*Power(xj, 14LL)*

                    (4867016286825LL + 1199363925375LL*Power(r, 2LL)*Power(xj, 2LL) +

            26817947100LL*Power(r, 4LL)*Power(xj, 4LL) -

            167333418LL*Power(r, 6LL)*Power(xj, 6LL) -

            1476138LL*Power(r, 8LL)*Power(xj, 8LL) + 1294LL*Power(r, 10LL)*Power(xj, 10LL))

                    + 66LL*Power(xi, 16LL)*Power(xj, 16LL)*

                    (-1657759205025LL - 682207855875LL*Power(r, 2LL)*Power(xj, 2LL) -

            31509229500LL*Power(r, 4LL)*Power(xj, 4LL) -

            492146550LL*Power(r, 6LL)*Power(xj, 6LL) -

            11910LL*Power(r, 8LL)*Power(xj, 8LL) + 2594LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    1386LL*Power(xi, 26LL)*Power(xj, 6LL)*

                    (-6066588045375LL + 98854491375LL*Power(r, 2LL)*Power(xj, 2LL) -

            12496954500LL*Power(r, 4LL)*Power(xj, 4LL) +

            420813750LL*Power(r, 6LL)*Power(xj, 6LL) -

            2881210LL*Power(r, 8LL)*Power(xj, 8LL) + 2622LL*Power(r, 10LL)*Power(xj, 10LL))

                    + 11LL*r*Power(xi, 23LL)*Power(xj, 10LL)*

                    (149900659402725LL - 26541339882750LL*Power(r, 2LL)*Power(xj, 2LL) +

            594745455150LL*Power(r, 4LL)*Power(xj, 4LL) -

            1399125420LL*Power(r, 6LL)*Power(xj, 6LL) -

            7887390LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                    - 11LL*r*Power(xi, 31LL)*Power(xj, 2LL)*

                    (7685082491625LL + 5034333946950LL*Power(r, 2LL)*Power(xj, 2LL) -

            108088893990LL*Power(r, 4LL)*Power(xj, 4LL) +

            1254174300LL*Power(r, 6LL)*Power(xj, 6LL) -

            6355950LL*Power(r, 8LL)*Power(xj, 8LL) + 4232LL*Power(r, 10LL)*Power(xj, 10LL))

                    - 462LL*Power(xi, 24LL)*Power(xj, 8LL)*

                    (40495013164125LL - 3973079865375LL*Power(r, 2LL)*Power(xj, 2LL) +

            110288047500LL*Power(r, 4LL)*Power(xj, 4LL) -

            381623130LL*Power(r, 6LL)*Power(xj, 6LL) -

            4811370LL*Power(r, 8LL)*Power(xj, 8LL) + 9338LL*Power(r, 10LL)*Power(xj, 10LL))

                    + 198LL*Power(xi, 32LL)*(-9126526500LL - 152866565775LL*Power(r, 2LL)*Power(xj, 2LL) -

                                             32383266300LL*Power(r, 4LL)*Power(xj, 4LL) +

                                             709444890LL*Power(r, 6LL)*Power(xj, 6LL) -

                                             5562070LL*Power(r, 8LL)*Power(xj, 8LL) + 11042LL*Power(r, 10LL)*Power(xj, 10LL))

                    + 2LL*r*Power(xi, 33LL)*(-764522104500LL - 3357151476525LL*Power(r, 2LL)*Power(xj, 2LL) -

                                             242177564475LL*Power(r, 4LL)*Power(xj, 4LL) +

                                             4513719870LL*Power(r, 6LL)*Power(xj, 6LL) -

                                             20531775LL*Power(r, 8LL)*Power(xj, 8LL) + 11236LL*Power(r, 10LL)*Power(xj, 10LL)

                                             ) - r*Power(xi, 21LL)*Power(xj, 12LL)*

                    (-5533525427435775LL + 138591131159250LL*Power(r, 2LL)*Power(xj, 2LL) +

            2815739907750LL*Power(r, 4LL)*Power(xj, 4LL) -

            32922004500LL*Power(r, 6LL)*Power(xj, 6LL) +

            11347050LL*Power(r, 8LL)*Power(xj, 8LL) + 22472LL*Power(r, 10LL)*Power(xj, 10LL)

                    ) + 66LL*Power(xi, 22LL)*Power(xj, 10LL)*

                    (-283522589265825LL + 7639225988625LL*Power(r, 2LL)*Power(xj, 2LL) +

            480728209500LL*Power(r, 4LL)*Power(xj, 4LL) -

            8458349130LL*Power(r, 6LL)*Power(xj, 6LL) +

            9771090LL*Power(r, 8LL)*Power(xj, 8LL) + 31786LL*Power(r, 10LL)*Power(xj, 10LL))

                    - 66LL*Power(xi, 30LL)*Power(xj, 2LL)*(1678609807875LL +

                                                           4713298976925LL*Power(r, 2LL)*Power(xj, 2LL) -

                                                           30578971500LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                           53723250LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                           9140190LL*Power(r, 8LL)*Power(xj, 8LL) + 38042LL*Power(r, 10LL)*Power(xj, 10LL)))

                   + exp(2LL*r*xi)*Power(xi, 14LL)*

                   (-302841LL*Power(xi, 16LL)*Power(xj, 16LL)*

                    (1346857875LL*xj - 2613847500LL*r*Power(xj, 2LL) +

            964581750LL*Power(r, 2LL)*Power(xj, 3LL) +

            222948000LL*Power(r, 3LL)*Power(xj, 4LL) -

            46488750LL*Power(r, 4LL)*Power(xj, 5LL) -

            11060280LL*Power(r, 5LL)*Power(xj, 6LL) -

            355740LL*Power(r, 6LL)*Power(xj, 7LL) + 58720LL*Power(r, 7LL)*Power(xj, 8LL) +

            5130LL*Power(r, 8LL)*Power(xj, 9LL) + 120LL*Power(r, 9LL)*Power(xj, 10LL)) +

                    12LL*Power(xj, 32LL)*(127420350750LL*xj + 101936280600LL*r*Power(xj, 2LL) +

                                          38226105225LL*Power(r, 2LL)*Power(xj, 3LL) +

                                          8864024400LL*Power(r, 3LL)*Power(xj, 4LL) +

                                          1410185700LL*Power(r, 4LL)*Power(xj, 5LL) +

                                          161164080LL*Power(r, 5LL)*Power(xj, 6LL) +

                                          13430340LL*Power(r, 6LL)*Power(xj, 7LL) +

                                          807840LL*Power(r, 7LL)*Power(xj, 8LL) + 33660LL*Power(r, 8LL)*Power(xj, 9LL) +

                                          880LL*Power(r, 9LL)*Power(xj, 10LL) + 11LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    3LL*Power(xi, 32LL)*(1715175LL*xj + 3118500LL*r*Power(xj, 2LL) +

                                         2806650LL*Power(r, 2LL)*Power(xj, 3LL) +

                                         1663200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                         727650LL*Power(r, 4LL)*Power(xj, 5LL) + 249480LL*Power(r, 5LL)*Power(xj, 6LL) +

                                         69300LL*Power(r, 6LL)*Power(xj, 7LL) + 15840LL*Power(r, 7LL)*Power(xj, 8LL) +

                                         2970LL*Power(r, 8LL)*Power(xj, 9LL) + 440LL*Power(r, 9LL)*Power(xj, 10LL) +

                                         44LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 30LL)*Power(xj, 2LL)*

                    (-10758825LL*xj - 19561500LL*r*Power(xj, 2LL) -

            17605350LL*Power(r, 2LL)*Power(xj, 3LL) -

            10432800LL*Power(r, 3LL)*Power(xj, 4LL) -

            4564350LL*Power(r, 4LL)*Power(xj, 5LL) -

            1564920LL*Power(r, 5LL)*Power(xj, 6LL) -

            434700LL*Power(r, 6LL)*Power(xj, 7LL) - 99360LL*Power(r, 7LL)*Power(xj, 8LL) -

            18630LL*Power(r, 8LL)*Power(xj, 9LL) - 2760LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    5313LL*Power(xi, 14LL)*Power(xj, 18LL)*

                    (495525217275LL*xj - 323789251500LL*r*Power(xj, 2LL) -

            78255861750LL*Power(r, 2LL)*Power(xj, 3LL) +

            23887116000LL*Power(r, 3LL)*Power(xj, 4LL) +

            6156785250LL*Power(r, 4LL)*Power(xj, 5LL) +

            199107720LL*Power(r, 5LL)*Power(xj, 6LL) -

            54382860LL*Power(r, 6LL)*Power(xj, 7LL) -

            6012960LL*Power(r, 7LL)*Power(xj, 8LL) -

            208710LL*Power(r, 8LL)*Power(xj, 9LL) - 200LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    5313LL*Power(xi, 18LL)*Power(xj, 14LL)*

                    (-3082655475LL*xj + 16949263500LL*r*Power(xj, 2LL) -

            20439843750LL*Power(r, 2LL)*Power(xj, 3LL) +

            6662844000LL*Power(r, 3LL)*Power(xj, 4LL) +

            1164980250LL*Power(r, 4LL)*Power(xj, 5LL) -

            236862360LL*Power(r, 5LL)*Power(xj, 6LL) -

            43375500LL*Power(r, 6LL)*Power(xj, 7LL) -

            971040LL*Power(r, 7LL)*Power(xj, 8LL) + 146970LL*Power(r, 8LL)*Power(xj, 9LL) +

            8120LL*Power(r, 9LL)*Power(xj, 10LL) + 88LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 2LL)*Power(xj, 30LL)*

                    (7685082491625LL*xj + 5503196367900LL*r*Power(xj, 2LL) +

            1831173532650LL*Power(r, 2LL)*Power(xj, 3LL) +

            372857839200LL*Power(r, 3LL)*Power(xj, 4LL) +

            51426531450LL*Power(r, 4LL)*Power(xj, 5LL) +

            5014616040LL*Power(r, 5LL)*Power(xj, 6LL) +

            349260660LL*Power(r, 6LL)*Power(xj, 7LL) +

            17076960LL*Power(r, 7LL)*Power(xj, 8LL) +

            555930LL*Power(r, 8LL)*Power(xj, 9LL) + 10680LL*Power(r, 9LL)*Power(xj, 10LL) +

            88LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 28LL)*Power(xj, 4LL)*

                    (-118347075LL*xj - 215176500LL*r*Power(xj, 2LL) -

            193658850LL*Power(r, 2LL)*Power(xj, 3LL) -

            114760800LL*Power(r, 3LL)*Power(xj, 4LL) -

            50207850LL*Power(r, 4LL)*Power(xj, 5LL) -

            17214120LL*Power(r, 5LL)*Power(xj, 6LL) -

            4781700LL*Power(r, 6LL)*Power(xj, 7LL) -

            1092960LL*Power(r, 7LL)*Power(xj, 8LL) -

            303210LL*Power(r, 8LL)*Power(xj, 9LL) + 13320LL*Power(r, 9LL)*Power(xj, 10LL) +

            968LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 4LL)*Power(xj, 28LL)*

                    (-94924521239625LL*xj - 56559587723100LL*r*Power(xj, 2LL) -

            15103001840850LL*Power(r, 2LL)*Power(xj, 3LL) -

            2331595173600LL*Power(r, 3LL)*Power(xj, 4LL) -

            220161422250LL*Power(r, 4LL)*Power(xj, 5LL) -

            11585102760LL*Power(r, 5LL)*Power(xj, 6LL) -

            107023140LL*Power(r, 6LL)*Power(xj, 7LL) +

            30598560LL*Power(r, 7LL)*Power(xj, 8LL) +

            2282310LL*Power(r, 8LL)*Power(xj, 9LL) +

            73800LL*Power(r, 9LL)*Power(xj, 10LL) + 968LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    253LL*Power(xi, 20LL)*Power(xj, 12LL)*

                    (2437660575LL*xj - 3959077500LL*r*Power(xj, 2LL) +

            32973459750LL*Power(r, 2LL)*Power(xj, 3LL) -

            32879196000LL*Power(r, 3LL)*Power(xj, 4LL) +

            12414984750LL*Power(r, 4LL)*Power(xj, 5LL) +

            1309565880LL*Power(r, 5LL)*Power(xj, 6LL) -

            334866420LL*Power(r, 6LL)*Power(xj, 7LL) -

            41211360LL*Power(r, 7LL)*Power(xj, 8LL) -

            403650LL*Power(r, 8LL)*Power(xj, 9LL) + 82920LL*Power(r, 9LL)*Power(xj, 10LL) +

            2024LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    253LL*Power(xi, 12LL)*Power(xj, 20LL)*

                    (-21871642005675LL*xj - 3985683124500LL*r*Power(xj, 2LL) +

            3461913897750LL*Power(r, 2LL)*Power(xj, 3LL) +

            805582260000LL*Power(r, 3LL)*Power(xj, 4LL) +

            6607392750LL*Power(r, 4LL)*Power(xj, 5LL) -

            13831965000LL*Power(r, 5LL)*Power(xj, 6LL) -

            1624632660LL*Power(r, 6LL)*Power(xj, 7LL) -

            67004640LL*Power(r, 7LL)*Power(xj, 8LL) +

            294030LL*Power(r, 8LL)*Power(xj, 9LL) + 99240LL*Power(r, 9LL)*Power(xj, 10LL) +

            2024LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 6LL)*Power(xj, 26LL)*

                    (377817065789625LL*xj + 151494770958300LL*r*Power(xj, 2LL) +

            19418753866050LL*Power(r, 2LL)*Power(xj, 3LL) -

            733895316000LL*Power(r, 3LL)*Power(xj, 4LL) -

            540444469950LL*Power(r, 4LL)*Power(xj, 5LL) -

            76620048120LL*Power(r, 5LL)*Power(xj, 6LL) -

            5744734380LL*Power(r, 6LL)*Power(xj, 7LL) -

            239358240LL*Power(r, 7LL)*Power(xj, 8LL) -

            4241430LL*Power(r, 8LL)*Power(xj, 9LL) +

            42360LL*Power(r, 9LL)*Power(xj, 10LL) + 2200LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 26LL)*Power(xj, 6LL)*

                    (-828429525LL*xj - 1506235500LL*r*Power(xj, 2LL) -

            1355611950LL*Power(r, 2LL)*Power(xj, 3LL) -

            803325600LL*Power(r, 3LL)*Power(xj, 4LL) -

            351454950LL*Power(r, 4LL)*Power(xj, 5LL) -

            120498840LL*Power(r, 5LL)*Power(xj, 6LL) -

            23447340LL*Power(r, 6LL)*Power(xj, 7LL) -

            19107360LL*Power(r, 7LL)*Power(xj, 8LL) +

            601290LL*Power(r, 8LL)*Power(xj, 9LL) +

            155640LL*Power(r, 9LL)*Power(xj, 10LL) + 2200LL*Power(r, 10LL)*Power(xj, 11LL))

                    - 11LL*Power(xi, 24LL)*Power(xj, 8LL)*(4142147625LL*xj + 7531177500LL*r*Power(xj, 2LL) +

                                                           6778059750LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                           4016628000LL*Power(r, 3LL)*Power(xj, 4LL) +

                                                           2154080250LL*Power(r, 4LL)*Power(xj, 5LL) -

                                                           349839000LL*Power(r, 5LL)*Power(xj, 6LL) +

                                                           730405620LL*Power(r, 6LL)*Power(xj, 7LL) -

                                                           571680LL*Power(r, 7LL)*Power(xj, 8LL) -

                                                           10089630LL*Power(r, 8LL)*Power(xj, 9LL) -

                                                           388200LL*Power(r, 9LL)*Power(xj, 10LL) + 2728LL*Power(r, 10LL)*Power(xj, 11LL))

                    + 11LL*Power(xi, 8LL)*Power(xj, 24LL)*(450334446494625LL*xj -

                                                           24911331826500LL*r*Power(xj, 2LL) -

                                                           59323593341250LL*Power(r, 2LL)*Power(xj, 3LL) -

                                                           13144611765600LL*Power(r, 3LL)*Power(xj, 4LL) -

                                                           1123650237150LL*Power(r, 4LL)*Power(xj, 5LL) +

                                                           1934037000LL*Power(r, 5LL)*Power(xj, 6LL) +

                                                           8779220100LL*Power(r, 6LL)*Power(xj, 7LL) +

                                                           800938080LL*Power(r, 7LL)*Power(xj, 8LL) +

                                                           33597450LL*Power(r, 8LL)*Power(xj, 9LL) +

                                                           633720LL*Power(r, 9LL)*Power(xj, 10LL) + 2728LL*Power(r, 10LL)*Power(xj, 11LL))

                    + Power(xi, 22LL)*Power(xj, 10LL)*(173141770725LL*xj +

                                                       314803219500LL*r*Power(xj, 2LL) +

                                                       228325655250LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                       461213676000LL*Power(r, 3LL)*Power(xj, 4LL) -

                                                       339461718750LL*Power(r, 4LL)*Power(xj, 5LL) +

                                                       194890033800LL*Power(r, 5LL)*Power(xj, 6LL) +

                                                       8778328020LL*Power(r, 6LL)*Power(xj, 7LL) -

                                                       3897004320LL*Power(r, 7LL)*Power(xj, 8LL) -

                                                       288989910LL*Power(r, 8LL)*Power(xj, 9LL) +

                                                       384120LL*Power(r, 9LL)*Power(xj, 10LL) + 247192LL*Power(r, 10LL)*Power(xj, 11LL)

                                                       ) - Power(xi, 10LL)*Power(xj, 22LL)*(1648907253429975LL*xj +

                                                                                            3671125795606500LL*r*Power(xj, 2LL) +

                                                                                            630531674561250LL*Power(r, 2LL)*Power(xj, 3LL) -

                                                                                            69283115748000LL*Power(r, 3LL)*Power(xj, 4LL) -

                                                                                            29572528119750LL*Power(r, 4LL)*Power(xj, 5LL) -

                                                                                            3234734478360LL*Power(r, 5LL)*Power(xj, 6LL) -

                                                                                            120582706860LL*Power(r, 6LL)*Power(xj, 7LL) +

                                                                                            4826020320LL*Power(r, 7LL)*Power(xj, 8LL) +

                                                                                            629239050LL*Power(r, 8LL)*Power(xj, 9LL) +

                                                                                            21863160LL*Power(r, 9LL)*Power(xj, 10LL) +

                                                                                            247192LL*Power(r, 10LL)*Power(xj, 11LL))) +

                   2LL*exp(2LL*r*xi)*Power(xi, 15LL)*

                   (-302841LL*Power(xi, 16LL)*Power(xj, 16LL)*

                    (-361285650LL + 1346857875LL*r*xj -

            1306923750LL*Power(r, 2LL)*Power(xj, 2LL) +

            321527250LL*Power(r, 3LL)*Power(xj, 3LL) +

            55737000LL*Power(r, 4LL)*Power(xj, 4LL) -

            9297750LL*Power(r, 5LL)*Power(xj, 5LL) -

            1843380LL*Power(r, 6LL)*Power(xj, 6LL) - 50820LL*Power(r, 7LL)*Power(xj, 7LL) +

            7340LL*Power(r, 8LL)*Power(xj, 8LL) + 570LL*Power(r, 9LL)*Power(xj, 9LL) +

            12LL*Power(r, 10LL)*Power(xj, 10LL)) +

                    12LL*Power(xj, 32LL)*(150587687250LL + 127420350750LL*r*xj +

                                          50968140300LL*Power(r, 2LL)*Power(xj, 2LL) +

                                          12742035075LL*Power(r, 3LL)*Power(xj, 3LL) +

                                          2216006100LL*Power(r, 4LL)*Power(xj, 4LL) +

                                          282037140LL*Power(r, 5LL)*Power(xj, 5LL) +

                                          26860680LL*Power(r, 6LL)*Power(xj, 6LL) +

                                          1918620LL*Power(r, 7LL)*Power(xj, 7LL) + 100980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                          3740LL*Power(r, 9LL)*Power(xj, 9LL) + 88LL*Power(r, 10LL)*Power(xj, 10LL) +

                                          Power(r, 11LL)*Power(xj, 11LL)) -

                    3LL*Power(xi, 32LL)*(935550LL + 1715175LL*r*xj +

                                         1559250LL*Power(r, 2LL)*Power(xj, 2LL) + 935550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                         415800LL*Power(r, 4LL)*Power(xj, 4LL) + 145530LL*Power(r, 5LL)*Power(xj, 5LL) +

                                         41580LL*Power(r, 6LL)*Power(xj, 6LL) + 9900LL*Power(r, 7LL)*Power(xj, 7LL) +

                                         1980LL*Power(r, 8LL)*Power(xj, 8LL) + 330LL*Power(r, 9LL)*Power(xj, 9LL) +

                                         44LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 30LL)*Power(xj, 2LL)*

                    (-5868450LL - 10758825LL*r*xj - 9780750LL*Power(r, 2LL)*Power(xj, 2LL) -

            5868450LL*Power(r, 3LL)*Power(xj, 3LL) -

            2608200LL*Power(r, 4LL)*Power(xj, 4LL) - 912870LL*Power(r, 5LL)*Power(xj, 5LL) -

            260820LL*Power(r, 6LL)*Power(xj, 6LL) - 62100LL*Power(r, 7LL)*Power(xj, 7LL) -

            12420LL*Power(r, 8LL)*Power(xj, 8LL) - 2070LL*Power(r, 9LL)*Power(xj, 9LL) -

            276LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    5313LL*Power(xi, 14LL)*Power(xj, 18LL)*

                    (-302299148250LL + 495525217275LL*r*xj -

            161894625750LL*Power(r, 2LL)*Power(xj, 2LL) -

            26085287250LL*Power(r, 3LL)*Power(xj, 3LL) +

            5971779000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1231357050LL*Power(r, 5LL)*Power(xj, 5LL) +

            33184620LL*Power(r, 6LL)*Power(xj, 6LL) -

            7768980LL*Power(r, 7LL)*Power(xj, 7LL) - 751620LL*Power(r, 8LL)*Power(xj, 8LL) -

            23190LL*Power(r, 9LL)*Power(xj, 9LL) - 20LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    5313LL*Power(xi, 18LL)*Power(xj, 14LL)*

                    (469625850LL - 3082655475LL*r*xj + 8474631750LL*Power(r, 2LL)*Power(xj, 2LL) -

            6813281250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1665711000LL*Power(r, 4LL)*Power(xj, 4LL) +

            232996050LL*Power(r, 5LL)*Power(xj, 5LL) -

            39477060LL*Power(r, 6LL)*Power(xj, 6LL) -

            6196500LL*Power(r, 7LL)*Power(xj, 7LL) - 121380LL*Power(r, 8LL)*Power(xj, 8LL) +

            16330LL*Power(r, 9LL)*Power(xj, 9LL) + 812LL*Power(r, 10LL)*Power(xj, 10LL) +

            8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 2LL)*Power(xj, 30LL)*

                    (10071658847250LL + 7685082491625LL*r*xj +

            2751598183950LL*Power(r, 2LL)*Power(xj, 2LL) +

            610391177550LL*Power(r, 3LL)*Power(xj, 3LL) +

            93214459800LL*Power(r, 4LL)*Power(xj, 4LL) +

            10285306290LL*Power(r, 5LL)*Power(xj, 5LL) +

            835769340LL*Power(r, 6LL)*Power(xj, 6LL) +

            49894380LL*Power(r, 7LL)*Power(xj, 7LL) +

            2134620LL*Power(r, 8LL)*Power(xj, 8LL) + 61770LL*Power(r, 9LL)*Power(xj, 9LL) +

            1068LL*Power(r, 10LL)*Power(xj, 10LL) + 8LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 28LL)*Power(xj, 4LL)*

                    (-64552950LL - 118347075LL*r*xj - 107588250LL*Power(r, 2LL)*Power(xj, 2LL) -

            64552950LL*Power(r, 3LL)*Power(xj, 3LL) -

            28690200LL*Power(r, 4LL)*Power(xj, 4LL) -

            10041570LL*Power(r, 5LL)*Power(xj, 5LL) -

            2869020LL*Power(r, 6LL)*Power(xj, 6LL) - 683100LL*Power(r, 7LL)*Power(xj, 7LL) -

            136620LL*Power(r, 8LL)*Power(xj, 8LL) - 33690LL*Power(r, 9LL)*Power(xj, 9LL) +

            1332LL*Power(r, 10LL)*Power(xj, 10LL) + 88LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 4LL)*Power(xj, 28LL)*

                    (-145801982576250LL - 94924521239625LL*r*xj -

            28279793861550LL*Power(r, 2LL)*Power(xj, 2LL) -

            5034333946950LL*Power(r, 3LL)*Power(xj, 3LL) -

            582898793400LL*Power(r, 4LL)*Power(xj, 4LL) -

            44032284450LL*Power(r, 5LL)*Power(xj, 5LL) -

            1930850460LL*Power(r, 6LL)*Power(xj, 6LL) -

            15289020LL*Power(r, 7LL)*Power(xj, 7LL) +

            3824820LL*Power(r, 8LL)*Power(xj, 8LL) + 253590LL*Power(r, 9LL)*Power(xj, 9LL) +

            7380LL*Power(r, 10LL)*Power(xj, 10LL) + 88LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    253LL*Power(xi, 20LL)*Power(xj, 12LL)*

                    (1119853350LL + 2437660575LL*r*xj -

            1979538750LL*Power(r, 2LL)*Power(xj, 2LL) +

            10991153250LL*Power(r, 3LL)*Power(xj, 3LL) -

            8219799000LL*Power(r, 4LL)*Power(xj, 4LL) +

            2482996950LL*Power(r, 5LL)*Power(xj, 5LL) +

            218260980LL*Power(r, 6LL)*Power(xj, 6LL) -

            47838060LL*Power(r, 7LL)*Power(xj, 7LL) -

            5151420LL*Power(r, 8LL)*Power(xj, 8LL) - 44850LL*Power(r, 9LL)*Power(xj, 9LL) +

            8292LL*Power(r, 10LL)*Power(xj, 10LL) + 184LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    253LL*Power(xi, 12LL)*Power(xj, 20LL)*

                    (33221660229450LL - 21871642005675LL*r*xj -

            1992841562250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1153971299250LL*Power(r, 3LL)*Power(xj, 3LL) +

            201395565000LL*Power(r, 4LL)*Power(xj, 4LL) +

            1321478550LL*Power(r, 5LL)*Power(xj, 5LL) -

            2305327500LL*Power(r, 6LL)*Power(xj, 6LL) -

            232090380LL*Power(r, 7LL)*Power(xj, 7LL) -

            8375580LL*Power(r, 8LL)*Power(xj, 8LL) + 32670LL*Power(r, 9LL)*Power(xj, 9LL) +

            9924LL*Power(r, 10LL)*Power(xj, 10LL) + 184LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 6LL)*Power(xj, 26LL)*

                    (764390093717250LL + 377817065789625LL*r*xj +

            75747385479150LL*Power(r, 2LL)*Power(xj, 2LL) +

            6472917955350LL*Power(r, 3LL)*Power(xj, 3LL) -

            183473829000LL*Power(r, 4LL)*Power(xj, 4LL) -

            108088893990LL*Power(r, 5LL)*Power(xj, 5LL) -

            12770008020LL*Power(r, 6LL)*Power(xj, 6LL) -

            820676340LL*Power(r, 7LL)*Power(xj, 7LL) -

            29919780LL*Power(r, 8LL)*Power(xj, 8LL) -

            471270LL*Power(r, 9LL)*Power(xj, 9LL) + 4236LL*Power(r, 10LL)*Power(xj, 10LL) +

            200LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 26LL)*Power(xj, 6LL)*

                    (-451870650LL - 828429525LL*r*xj - 753117750LL*Power(r, 2LL)*Power(xj, 2LL) -

            451870650LL*Power(r, 3LL)*Power(xj, 3LL) -

            200831400LL*Power(r, 4LL)*Power(xj, 4LL) -

            70290990LL*Power(r, 5LL)*Power(xj, 5LL) -

            20083140LL*Power(r, 6LL)*Power(xj, 6LL) -

            3349620LL*Power(r, 7LL)*Power(xj, 7LL) -

            2388420LL*Power(r, 8LL)*Power(xj, 8LL) + 66810LL*Power(r, 9LL)*Power(xj, 9LL) +

            15564LL*Power(r, 10LL)*Power(xj, 10LL) + 200LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    11LL*Power(xi, 24LL)*Power(xj, 8LL)*

                    (2259353250LL + 4142147625LL*r*xj +

            3765588750LL*Power(r, 2LL)*Power(xj, 2LL) +

            2259353250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1004157000LL*Power(r, 4LL)*Power(xj, 4LL) +

            430816050LL*Power(r, 5LL)*Power(xj, 5LL) -

            58306500LL*Power(r, 6LL)*Power(xj, 6LL) +

            104343660LL*Power(r, 7LL)*Power(xj, 7LL) -

            71460LL*Power(r, 8LL)*Power(xj, 8LL) - 1121070LL*Power(r, 9LL)*Power(xj, 9LL) -

            38820LL*Power(r, 10LL)*Power(xj, 10LL) + 248LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    11LL*Power(xi, 8LL)*Power(xj, 24LL)*

                    (1700790552893250LL + 450334446494625LL*r*xj -

            12455665913250LL*Power(r, 2LL)*Power(xj, 2LL) -

            19774531113750LL*Power(r, 3LL)*Power(xj, 3LL) -

            3286152941400LL*Power(r, 4LL)*Power(xj, 4LL) -

            224730047430LL*Power(r, 5LL)*Power(xj, 5LL) +

            322339500LL*Power(r, 6LL)*Power(xj, 6LL) +

            1254174300LL*Power(r, 7LL)*Power(xj, 7LL) +

            100117260LL*Power(r, 8LL)*Power(xj, 8LL) +

            3733050LL*Power(r, 9LL)*Power(xj, 9LL) +

            63372LL*Power(r, 10LL)*Power(xj, 10LL) + 248LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    Power(xi, 22LL)*Power(xj, 10LL)*

                    (94440965850LL + 173141770725LL*r*xj +

            157401609750LL*Power(r, 2LL)*Power(xj, 2LL) +

            76108551750LL*Power(r, 3LL)*Power(xj, 3LL) +

            115303419000LL*Power(r, 4LL)*Power(xj, 4LL) -

            67892343750LL*Power(r, 5LL)*Power(xj, 5LL) +

            32481672300LL*Power(r, 6LL)*Power(xj, 6LL) +

            1254046860LL*Power(r, 7LL)*Power(xj, 7LL) -

            487125540LL*Power(r, 8LL)*Power(xj, 8LL) -

            32109990LL*Power(r, 9LL)*Power(xj, 9LL) +

            38412LL*Power(r, 10LL)*Power(xj, 10LL) + 22472LL*Power(r, 11LL)*Power(xj, 11LL))

                    - Power(xi, 10LL)*Power(xj, 22LL)*(-18712490891544450LL + 1648907253429975LL*r*xj +

                                                       1835562897803250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       210177224853750LL*Power(r, 3LL)*Power(xj, 3LL) -

                                                       17320778937000LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                       5914505623950LL*Power(r, 5LL)*Power(xj, 5LL) -

                                                       539122413060LL*Power(r, 6LL)*Power(xj, 6LL) -

                                                       17226100980LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                       603252540LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                       69915450LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                       2186316LL*Power(r, 10LL)*Power(xj, 10LL) + 22472LL*Power(r, 11LL)*Power(xj, 11LL))

                   ))/(2.80665e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 23LL)*Power(xi + xj, 23LL))

            ;
        }

    }
    return S;
}
