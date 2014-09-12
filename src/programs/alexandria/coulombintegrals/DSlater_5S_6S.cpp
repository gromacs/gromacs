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

cl_R DSlater_5S_6S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-(1e9*cl_float(283601367, precision)+780029361562LL)*xi + (1e9*cl_float(306545653, precision)+30256640000LL)*exp(2LL*r*xi)*xi -

                  (1e9*cl_float(521314165, precision)+59604166250LL)*r*Power(xi, 2LL) -

                  (1e9*cl_float(475834779, precision)+713452221562LL)*Power(r, 2LL)*Power(xi, 3LL) -

                  (1e9*cl_float(287449176, precision)+450603667500LL)*Power(r, 3LL)*Power(xi, 4LL) -

                  (1e9*cl_float(129241025, precision)+417554800750LL)*Power(r, 4LL)*Power(xi, 5LL) -

                  (1e9*cl_float(461120284, precision)+6963318800LL)*Power(r, 5LL)*Power(xi, 6LL) -

                  (1e9*cl_float(135925896, precision)+63380527800LL)*Power(r, 6LL)*Power(xi, 7LL) -

                  (1e9*cl_float(340254650, precision)+7763603200LL)*Power(r, 7LL)*Power(xi, 8LL) -

                  (1e9*cl_float(737703089, precision)+712388800LL)*Power(r, 8LL)*Power(xi, 9LL) -

                  (1e9*cl_float(140553592, precision)+289510400LL)*Power(r, 9LL)*Power(xi, 10LL) -

                  237880232087097600LL*Power(r, 10LL)*Power(xi, 11LL) -

                  36042992713728000LL*Power(r, 11LL)*Power(xi, 12LL) -

                  4914953551872000LL*Power(r, 12LL)*Power(xi, 13LL) -

                  604917360230400LL*Power(r, 13LL)*Power(xi, 14LL) -

                  67213040025600LL*Power(r, 14LL)*Power(xi, 15LL) -

                  6721304002560LL*Power(r, 15LL)*Power(xi, 16LL) -

                  600116428800LL*Power(r, 16LL)*Power(xi, 17LL) -

                  47067955200LL*Power(r, 17LL)*Power(xi, 18LL) -

                  3137863680LL*Power(r, 18LL)*Power(xi, 19LL) -

                  165150720LL*Power(r, 19LL)*Power(xi, 20LL) - 5505024LL*Power(r, 20LL)*Power(xi, 21LL))/

                (1.5327282651512832e21*exp(2LL*r*xi)*r) +

                (-(1e9*cl_float(153272826, precision)+515128320000LL) + (1e9*cl_float(153272826, precision)+515128320000LL)*exp(2LL*r*xi) -

                 (1e9*cl_float(283601367, precision)+780029361562LL)*r*xi -

                 (1e9*cl_float(260657082, precision)+529802083125LL)*Power(r, 2LL)*Power(xi, 2LL) -

                 (1e9*cl_float(158611593, precision)+237817407187LL)*Power(r, 3LL)*Power(xi, 3LL) -

                 (1e9*cl_float(718622941, precision)+12650916875LL)*Power(r, 4LL)*Power(xi, 4LL) -

                 (1e9*cl_float(258482050, precision)+83510960150LL)*Power(r, 5LL)*Power(xi, 5LL) -

                 (1e9*cl_float(768533806, precision)+7827219800LL)*Power(r, 6LL)*Power(xi, 6LL) -

                 (1e9*cl_float(194179852, precision)+3340075400LL)*Power(r, 7LL)*Power(xi, 7LL) -

                 (1e9*cl_float(425318313, precision)+470450400LL)*Power(r, 8LL)*Power(xi, 8LL) -

                 819670099680432000LL*Power(r, 9LL)*Power(xi, 9LL) -

                 140553592289510400LL*Power(r, 10LL)*Power(xi, 10LL) -

                 21625475644281600LL*Power(r, 11LL)*Power(xi, 11LL) -

                 3003582726144000LL*Power(r, 12LL)*Power(xi, 12LL) -

                 378073350144000LL*Power(r, 13LL)*Power(xi, 13LL) -

                 43208382873600LL*Power(r, 14LL)*Power(xi, 14LL) -

                 4480869335040LL*Power(r, 15LL)*Power(xi, 15LL) -

                 420081500160LL*Power(r, 16LL)*Power(xi, 16LL) -

                 35300966400LL*Power(r, 17LL)*Power(xi, 17LL) -

                 2614886400LL*Power(r, 18LL)*Power(xi, 18LL) -

                 165150720LL*Power(r, 19LL)*Power(xi, 19LL) -

                 8257536LL*Power(r, 20LL)*Power(xi, 20LL) - 262144LL*Power(r, 21LL)*Power(xi, 21LL))/

                (1.5327282651512832e21*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-(1e9*cl_float(153272826, precision)+515128320000LL) + (1e9*cl_float(153272826, precision)+515128320000LL)*exp(2LL*r*xi) -

                     (1e9*cl_float(283601367, precision)+780029361562LL)*r*xi -

                     (1e9*cl_float(260657082, precision)+529802083125LL)*Power(r, 2LL)*Power(xi, 2LL) -

                     (1e9*cl_float(158611593, precision)+237817407187LL)*Power(r, 3LL)*Power(xi, 3LL) -

                     (1e9*cl_float(718622941, precision)+12650916875LL)*Power(r, 4LL)*Power(xi, 4LL) -

                     (1e9*cl_float(258482050, precision)+83510960150LL)*Power(r, 5LL)*Power(xi, 5LL) -

                     (1e9*cl_float(768533806, precision)+7827219800LL)*Power(r, 6LL)*Power(xi, 6LL) -

                     (1e9*cl_float(194179852, precision)+3340075400LL)*Power(r, 7LL)*Power(xi, 7LL) -

                     (1e9*cl_float(425318313, precision)+470450400LL)*Power(r, 8LL)*Power(xi, 8LL) -

                     819670099680432000LL*Power(r, 9LL)*Power(xi, 9LL) -

                     140553592289510400LL*Power(r, 10LL)*Power(xi, 10LL) -

                     21625475644281600LL*Power(r, 11LL)*Power(xi, 11LL) -

                     3003582726144000LL*Power(r, 12LL)*Power(xi, 12LL) -

                     378073350144000LL*Power(r, 13LL)*Power(xi, 13LL) -

                     43208382873600LL*Power(r, 14LL)*Power(xi, 14LL) -

                     4480869335040LL*Power(r, 15LL)*Power(xi, 15LL) -

                     420081500160LL*Power(r, 16LL)*Power(xi, 16LL) -

                     35300966400LL*Power(r, 17LL)*Power(xi, 17LL) -

                     2614886400LL*Power(r, 18LL)*Power(xi, 18LL) -

                     165150720LL*Power(r, 19LL)*Power(xi, 19LL) -

                     8257536LL*Power(r, 20LL)*Power(xi, 20LL) - 262144LL*Power(r, 21LL)*Power(xi, 21LL)))/

                (7.663641325756416e20*exp(2LL*r*xi)*r)

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
            S = (4677750LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 21LL) +

                 110LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                 (-432LL*Power(r, 8LL)*Power(xi, 36LL) - 6LL*Power(r, 9LL)*Power(xi, 37LL) +

                  42525LL*Power(xj, 28LL) + 76545LL*r*xi*Power(xj, 28LL) +

                  19845LL*r*Power(xi, 3LL)*Power(xj, 26LL)*

                  (-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  1134LL*Power(r, 6LL)*Power(xi, 34LL)*(272LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  8LL*Power(r, 7LL)*Power(xi, 35LL)*(1836LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  8505LL*Power(xi, 2LL)*Power(xj, 26LL)*(-105LL + 8LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  378LL*Power(r, 5LL)*Power(xi, 33LL)*

                  (-11628LL - 666LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  5670LL*r*Power(xi, 5LL)*Power(xj, 24LL)*

                  (2835LL - 147LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  17010LL*Power(xi, 4LL)*Power(xj, 24LL)*

                  (525LL - 84LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                  378LL*Power(r, 4LL)*Power(xi, 32LL)*

                  (-116280LL - 17444LL*Power(r, 2LL)*Power(xj, 2LL) +

            59LL*Power(r, 4LL)*Power(xj, 4LL)) +

                  162LL*r*Power(xi, 7LL)*Power(xj, 22LL)*

                  (-628425LL + 51450LL*Power(r, 2LL)*Power(xj, 2LL) -

            735LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  378LL*Power(xi, 6LL)*Power(xj, 22LL)*

                  (-149625LL + 37800LL*Power(r, 2LL)*Power(xj, 2LL) -

            945LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                  18LL*Power(r, 3LL)*Power(xi, 31LL)*

                  (17093160LL + 6309387LL*Power(r, 2LL)*Power(xj, 2LL) -

            23562LL*Power(r, 4LL)*Power(xj, 4LL) + 16LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  54LL*Power(r, 2LL)*Power(xi, 30LL)*

                  (-26860680LL - 24843735LL*Power(r, 2LL)*Power(xj, 2LL) -

            40180LL*Power(r, 4LL)*Power(xj, 4LL) + 578LL*Power(r, 6LL)*Power(xj, 6LL)) +

                  378LL*r*Power(xi, 23LL)*Power(xj, 6LL)*

                  (-14625683325LL + 704051250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10752861LL*Power(r, 4LL)*Power(xj, 4LL) + 33478LL*Power(r, 6LL)*Power(xj, 6LL))

                  + 3LL*r*Power(xi, 9LL)*Power(xj, 20LL)*

                  (152707275LL - 17595900LL*Power(r, 2LL)*Power(xj, 2LL) +

            396900LL*Power(r, 4LL)*Power(xj, 4LL) - 2268LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  27LL*Power(xi, 8LL)*Power(xj, 20LL)*

                  (9426375LL - 3351600LL*Power(r, 2LL)*Power(xj, 2LL) +

            132300LL*Power(r, 4LL)*Power(xj, 4LL) - 1176LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  567LL*Power(xi, 10LL)*Power(xj, 18LL)*

                  (1526175LL - 718200LL*Power(r, 2LL)*Power(xj, 2LL) +

            39900LL*Power(r, 4LL)*Power(xj, 4LL) - 560LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  54LL*r*Power(xi, 13LL)*Power(xj, 16LL)*

                  (-1356769575LL - 127011675LL*Power(r, 2LL)*Power(xj, 2LL) -

            3867843LL*Power(r, 4LL)*Power(xj, 4LL) - 8556LL*Power(r, 6LL)*Power(xj, 6LL) +

            7LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  7LL*r*Power(xi, 11LL)*Power(xj, 18LL)*

                  (-151091325LL + 45272250LL*Power(r, 2LL)*Power(xj, 2LL) -

            647676LL*Power(r, 4LL)*Power(xj, 4LL) + 15336LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  18LL*r*Power(xi, 15LL)*Power(xj, 14LL)*

                  (63046289250LL + 3917182500LL*Power(r, 2LL)*Power(xj, 2LL) +

            10158435LL*Power(r, 4LL)*Power(xj, 4LL) -

            178842LL*Power(r, 6LL)*Power(xj, 6LL) + 16LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  378LL*r*Power(xi, 21LL)*Power(xj, 8LL)*

                  (-8559820125LL + 17573325LL*Power(r, 2LL)*Power(xj, 2LL) +

            7421001LL*Power(r, 4LL)*Power(xj, 4LL) -

            49096LL*Power(r, 6LL)*Power(xj, 6LL) + 19LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  378LL*Power(xi, 12LL)*Power(xj, 16LL)*

                  (17296650LL + 14244300LL*Power(r, 2LL)*Power(xj, 2LL) +

            360525LL*Power(r, 4LL)*Power(xj, 4LL) + 15928LL*Power(r, 6LL)*Power(xj, 6LL) +

            22LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  189LL*r*Power(xi, 25LL)*Power(xj, 4LL)*

                  (9994948425LL + 63821700LL*Power(r, 2LL)*Power(xj, 2LL) -

            1458540LL*Power(r, 4LL)*Power(xj, 4LL) -

            18756LL*Power(r, 6LL)*Power(xj, 6LL) + 38LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  189LL*Power(xi, 24LL)*Power(xj, 4LL)*

                  (17962854525LL + 4036942800LL*Power(r, 2LL)*Power(xj, 2LL) -

            126472500LL*Power(r, 4LL)*Power(xj, 4LL) +

            765464LL*Power(r, 6LL)*Power(xj, 6LL) + 190LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  21LL*r*Power(xi, 19LL)*Power(xj, 10LL)*

                  (-228066210225LL + 13487616450LL*Power(r, 2LL)*Power(xj, 2LL) -

            85465800LL*Power(r, 4LL)*Power(xj, 4LL) -

            320112LL*Power(r, 6LL)*Power(xj, 6LL) + 328LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  189LL*Power(xi, 18LL)*Power(xj, 10LL)*

                  (86069971575LL + 2157712200LL*Power(r, 2LL)*Power(xj, 2LL) -

            158179560LL*Power(r, 4LL)*Power(xj, 4LL) +

            578816LL*Power(r, 6LL)*Power(xj, 6LL) + 978LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  2LL*r*Power(xi, 29LL)*(2085060285LL +

                                         5450330025LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         127424745LL*Power(r, 4LL)*Power(xj, 4LL) -

                                         1398276LL*Power(r, 6LL)*Power(xj, 6LL) + 1159LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  378LL*Power(xi, 22LL)*Power(xj, 6LL)*

                  (37244490525LL - 2411839800LL*Power(r, 2LL)*Power(xj, 2LL) +

            92951775LL*Power(r, 4LL)*Power(xj, 4LL) -

            942172LL*Power(r, 6LL)*Power(xj, 6LL) + 1292LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  27LL*Power(xi, 16LL)*Power(xj, 12LL)*

                  (164245367475LL + 26909517600LL*Power(r, 2LL)*Power(xj, 2LL) +

            62674920LL*Power(r, 4LL)*Power(xj, 4LL) -

            3885112LL*Power(r, 6LL)*Power(xj, 6LL) + 2122LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  3LL*r*Power(xi, 27LL)*Power(xj, 2LL)*

                  (-63819198135LL - 21841975890LL*Power(r, 2LL)*Power(xj, 2LL) +

            442430100LL*Power(r, 4LL)*Power(xj, 4LL) -

            2756664LL*Power(r, 6LL)*Power(xj, 6LL) + 2296LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  r*Power(xi, 17LL)*Power(xj, 12LL)*

                  (4851990871875LL + 21622847400LL*Power(r, 2LL)*Power(xj, 2LL) -

            2153738160LL*Power(r, 4LL)*Power(xj, 4LL) +

            3608388LL*Power(r, 6LL)*Power(xj, 6LL) + 2318LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  18LL*Power(xi, 14LL)*Power(xj, 14LL)*

                  (-23418646650LL - 6922729800LL*Power(r, 2LL)*Power(xj, 2LL) -

            259958475LL*Power(r, 4LL)*Power(xj, 4LL) -

            697732LL*Power(r, 6LL)*Power(xj, 6LL) + 3030LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  126LL*Power(xi, 20LL)*Power(xj, 8LL)*

                  (-186637212225LL + 13028280300LL*Power(r, 2LL)*Power(xj, 2LL) -

            116198775LL*Power(r, 4LL)*Power(xj, 4LL) -

            1266160LL*Power(r, 6LL)*Power(xj, 6LL) + 4332LL*Power(r, 8LL)*Power(xj, 8LL)) -

                  54LL*Power(xi, 28LL)*(102965940LL + 1089437580LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        102508245LL*Power(r, 4LL)*Power(xj, 4LL) -

                                        1593144LL*Power(r, 6LL)*Power(xj, 6LL) + 4538LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  63LL*Power(xi, 26LL)*Power(xj, 2LL)*

                  (-4544129205LL - 7396000920LL*Power(r, 2LL)*Power(xj, 2LL) +

            149614020LL*Power(r, 4LL)*Power(xj, 4LL) -

            1684112LL*Power(r, 6LL)*Power(xj, 6LL) + 5922LL*Power(r, 8LL)*Power(xj, 8LL))) +

                 exp(2LL*r*xi)*Power(xi, 12LL)*

                 (6LL*Power(xi, 24LL)*Power(xj, 6LL)*

                  (1036901250LL + 1900985625LL*r*xj +

            1728168750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1036901250LL*Power(r, 3LL)*Power(xj, 3LL) +

            460845000LL*Power(r, 4LL)*Power(xj, 4LL) +

            161295750LL*Power(r, 5LL)*Power(xj, 5LL) +

            46084500LL*Power(r, 6LL)*Power(xj, 6LL) +

            9084900LL*Power(r, 7LL)*Power(xj, 7LL) +

            4082100LL*Power(r, 8LL)*Power(xj, 8LL) +

            121935LL*Power(r, 9LL)*Power(xj, 9LL) -

            21494LL*Power(r, 10LL)*Power(xj, 10LL) - 766LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  5LL*Power(xi, 28LL)*Power(xj, 2LL)*

                  (19646550LL + 36018675LL*r*xj + 32744250LL*Power(r, 2LL)*Power(xj, 2LL) +

            19646550LL*Power(r, 3LL)*Power(xj, 3LL) +

            8731800LL*Power(r, 4LL)*Power(xj, 4LL) +

            3056130LL*Power(r, 5LL)*Power(xj, 5LL) +

            873180LL*Power(r, 6LL)*Power(xj, 6LL) + 207900LL*Power(r, 7LL)*Power(xj, 7LL) +

            41580LL*Power(r, 8LL)*Power(xj, 8LL) + 6930LL*Power(r, 9LL)*Power(xj, 9LL) +

            924LL*Power(r, 10LL)*Power(xj, 10LL) - 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  26334LL*Power(xi, 16LL)*Power(xj, 14LL)*

                  (43880400LL - 186686775LL*r*xj + 576771750LL*Power(r, 2LL)*Power(xj, 2LL) -

            398603250LL*Power(r, 3LL)*Power(xj, 3LL) +

            72552600LL*Power(r, 4LL)*Power(xj, 4LL) +

            27903120LL*Power(r, 5LL)*Power(xj, 5LL) -

            342720LL*Power(r, 6LL)*Power(xj, 6LL) - 574800LL*Power(r, 7LL)*Power(xj, 7LL) -

            50800LL*Power(r, 8LL)*Power(xj, 8LL) - 945LL*Power(r, 9LL)*Power(xj, 9LL) +

            58LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  10LL*Power(xj, 30LL)*(97302813300LL + 89194245525LL*r*xj +

                                        38780106750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                        10576392750LL*Power(r, 3LL)*Power(xj, 3LL) +

                                        2014551000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                        282037140LL*Power(r, 5LL)*Power(xj, 5LL) +

                                        29688120LL*Power(r, 6LL)*Power(xj, 6LL) +

                                        2356200LL*Power(r, 7LL)*Power(xj, 7LL) +

                                        138600LL*Power(r, 8LL)*Power(xj, 8LL) + 5775LL*Power(r, 9LL)*Power(xj, 9LL) +

                                        154LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  10LL*Power(xi, 2LL)*Power(xj, 28LL)*

                  (4582499159700LL + 3733416276975LL*r*xj +

            1428215931450LL*Power(r, 2LL)*Power(xj, 2LL) +

            338545295550LL*Power(r, 3LL)*Power(xj, 3LL) +

            55198697400LL*Power(r, 4LL)*Power(xj, 4LL) +

            6486854220LL*Power(r, 5LL)*Power(xj, 5LL) +

            558419400LL*Power(r, 6LL)*Power(xj, 6LL) +

            34939080LL*Power(r, 7LL)*Power(xj, 7LL) +

            1532520LL*Power(r, 8LL)*Power(xj, 8LL) + 43285LL*Power(r, 9LL)*Power(xj, 9LL) +

            638LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  110LL*Power(xi, 10LL)*Power(xj, 20LL)*

                  (-14063418170550LL + 6795156458475LL*r*xj +

            2067471236250LL*Power(r, 2LL)*Power(xj, 2LL) -

            214664924250LL*Power(r, 3LL)*Power(xj, 3LL) -

            124416469800LL*Power(r, 4LL)*Power(xj, 4LL) -

            14935545450LL*Power(r, 5LL)*Power(xj, 5LL) -

            256688460LL*Power(r, 6LL)*Power(xj, 6LL) +

            105750900LL*Power(r, 7LL)*Power(xj, 7LL) +

            11502180LL*Power(r, 8LL)*Power(xj, 8LL) +

            518085LL*Power(r, 9LL)*Power(xj, 9LL) + 9294LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  55LL*Power(xi, 20LL)*Power(xj, 10LL)*

                  (1730682450LL + 3172917825LL*r*xj +

            2884470750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1571960250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1404081000LL*Power(r, 4LL)*Power(xj, 4LL) -

            426654270LL*Power(r, 5LL)*Power(xj, 5LL) +

            283536540LL*Power(r, 6LL)*Power(xj, 6LL) +

            39116700LL*Power(r, 7LL)*Power(xj, 7LL) -

            2659860LL*Power(r, 8LL)*Power(xj, 8LL) -

            528850LL*Power(r, 9LL)*Power(xj, 9LL) -

            18236LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  5LL*Power(xi, 30LL)*(935550LL + 1715175LL*r*xj +

                                       1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                       935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                       145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                       9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                       330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                       4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  13167LL*Power(xi, 14LL)*Power(xj, 16LL)*

                  (-2319354450LL + 8540029575LL*r*xj -

            7335672750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1133154750LL*Power(r, 3LL)*Power(xj, 3LL) +

            575014200LL*Power(r, 4LL)*Power(xj, 4LL) -

            913710LL*Power(r, 5LL)*Power(xj, 5LL) -

            14863940LL*Power(r, 6LL)*Power(xj, 6LL) -

            1687300LL*Power(r, 7LL)*Power(xj, 7LL) - 46900LL*Power(r, 8LL)*Power(xj, 8LL) +

            3210LL*Power(r, 9LL)*Power(xj, 9LL) + 236LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  770LL*Power(xi, 18LL)*Power(xj, 12LL)*

                  (329653800LL + 654729075LL*r*xj + 45785250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1602483750LL*Power(r, 3LL)*Power(xj, 3LL) -

            915705000LL*Power(r, 4LL)*Power(xj, 4LL) +

            266036400LL*Power(r, 5LL)*Power(xj, 5LL) +

            63745920LL*Power(r, 6LL)*Power(xj, 6LL) -

            2304000LL*Power(r, 7LL)*Power(xj, 7LL) -

            1074240LL*Power(r, 8LL)*Power(xj, 8LL) - 64635LL*Power(r, 9LL)*Power(xj, 9LL) -

            514LL*Power(r, 10LL)*Power(xj, 10LL) + 34LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  385LL*Power(xi, 12LL)*Power(xj, 18LL)*

                  (973565393850LL - 1429122323475LL*r*xj +

            298281831750LL*Power(r, 2LL)*Power(xj, 2LL) +

            138841148250LL*Power(r, 3LL)*Power(xj, 3LL) -

            2454240600LL*Power(r, 4LL)*Power(xj, 4LL) -

            4925394810LL*Power(r, 5LL)*Power(xj, 5LL) -

            623832300LL*Power(r, 6LL)*Power(xj, 6LL) -

            19098540LL*Power(r, 7LL)*Power(xj, 7LL) +

            2083140LL*Power(r, 8LL)*Power(xj, 8LL) +

            212430LL*Power(r, 9LL)*Power(xj, 9LL) + 7012LL*Power(r, 10LL)*Power(xj, 10LL) +

            68LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  14LL*Power(xi, 26LL)*Power(xj, 4LL)*

                  (-70166250LL - 128638125LL*r*xj - 116943750LL*Power(r, 2LL)*Power(xj, 2LL) -

            70166250LL*Power(r, 3LL)*Power(xj, 3LL) -

            31185000LL*Power(r, 4LL)*Power(xj, 4LL) -

            10914750LL*Power(r, 5LL)*Power(xj, 5LL) -

            3118500LL*Power(r, 6LL)*Power(xj, 6LL) -

            742500LL*Power(r, 7LL)*Power(xj, 7LL) - 148500LL*Power(r, 8LL)*Power(xj, 8LL) -

            32615LL*Power(r, 9LL)*Power(xj, 9LL) - 154LL*Power(r, 10LL)*Power(xj, 10LL) +

            74LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  7LL*Power(xi, 4LL)*Power(xj, 26LL)*

                  (-69822945249750LL - 46669577290875LL*r*xj -

            14025037430250LL*Power(r, 2LL)*Power(xj, 2LL) -

            2430881664750LL*Power(r, 3LL)*Power(xj, 3LL) -

            251629270200LL*Power(r, 4LL)*Power(xj, 4LL) -

            12434519790LL*Power(r, 5LL)*Power(xj, 5LL) +

            452930940LL*Power(r, 6LL)*Power(xj, 6LL) +

            131125500LL*Power(r, 7LL)*Power(xj, 7LL) +

            11018700LL*Power(r, 8LL)*Power(xj, 8LL) +

            514470LL*Power(r, 9LL)*Power(xj, 9LL) +

            13332LL*Power(r, 10LL)*Power(xj, 10LL) + 148LL*Power(r, 11LL)*Power(xj, 11LL)) -

                  50LL*Power(xi, 8LL)*Power(xj, 22LL)*

                  (-51768833574150LL - 5003280391725LL*r*xj +

            4493439477450LL*Power(r, 2LL)*Power(xj, 2LL) +

            1286866176750LL*Power(r, 3LL)*Power(xj, 3LL) +

            111437476920LL*Power(r, 4LL)*Power(xj, 4LL) -

            6620313546LL*Power(r, 5LL)*Power(xj, 5LL) -

            2406603276LL*Power(r, 6LL)*Power(xj, 6LL) -

            242686620LL*Power(r, 7LL)*Power(xj, 7LL) -

            12228876LL*Power(r, 8LL)*Power(xj, 8LL) -

            256223LL*Power(r, 9LL)*Power(xj, 9LL) + 2486LL*Power(r, 10LL)*Power(xj, 10LL) +

            158LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  25LL*Power(xi, 22LL)*Power(xj, 8LL)*

                  (-1119853350LL - 2053064475LL*r*xj -

            1866422250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1119853350LL*Power(r, 3LL)*Power(xj, 3LL) -

            497712600LL*Power(r, 4LL)*Power(xj, 4LL) -

            194415606LL*Power(r, 5LL)*Power(xj, 5LL) -

            9338868LL*Power(r, 6LL)*Power(xj, 6LL) -

            31217076LL*Power(r, 7LL)*Power(xj, 7LL) -

            2256804LL*Power(r, 8LL)*Power(xj, 8LL) +

            246774LL*Power(r, 9LL)*Power(xj, 9LL) +

            22836LL*Power(r, 10LL)*Power(xj, 10LL) + 316LL*Power(r, 11LL)*Power(xj, 11LL)) +

                  3LL*Power(xi, 6LL)*Power(xj, 24LL)*

                  (596006592662250LL + 266778699697125LL*r*xj +

            37515651153750LL*Power(r, 2LL)*Power(xj, 2LL) -

            2214626163750LL*Power(r, 3LL)*Power(xj, 3LL) -

            1538075107800LL*Power(r, 4LL)*Power(xj, 4LL) -

            248955308910LL*Power(r, 5LL)*Power(xj, 5LL) -

            21434337540LL*Power(r, 6LL)*Power(xj, 6LL) -

            957980100LL*Power(r, 7LL)*Power(xj, 7LL) -

            4874100LL*Power(r, 8LL)*Power(xj, 8LL) +

            1831830LL*Power(r, 9LL)*Power(xj, 9LL) +

            91828LL*Power(r, 10LL)*Power(xj, 10LL) + 1532LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (4.67775e6*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 21LL)*

                 Power(xi + xj, 21LL)) + (4677750LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 21LL) +

                                          110LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                                          (-432LL*Power(r, 8LL)*Power(xi, 36LL) - 6LL*Power(r, 9LL)*Power(xi, 37LL) +

                                  42525LL*Power(xj, 28LL) + 76545LL*r*xi*Power(xj, 28LL) +

                                  19845LL*r*Power(xi, 3LL)*Power(xj, 26LL)*

                                  (-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  1134LL*Power(r, 6LL)*Power(xi, 34LL)*(272LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  8LL*Power(r, 7LL)*Power(xi, 35LL)*(1836LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  8505LL*Power(xi, 2LL)*Power(xj, 26LL)*(-105LL + 8LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  378LL*Power(r, 5LL)*Power(xi, 33LL)*

                                  (-11628LL - 666LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  5670LL*r*Power(xi, 5LL)*Power(xj, 24LL)*

                                  (2835LL - 147LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  17010LL*Power(xi, 4LL)*Power(xj, 24LL)*

                                  (525LL - 84LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                                  378LL*Power(r, 4LL)*Power(xi, 32LL)*

                                  (-116280LL - 17444LL*Power(r, 2LL)*Power(xj, 2LL) +

            59LL*Power(r, 4LL)*Power(xj, 4LL)) +

                                  162LL*r*Power(xi, 7LL)*Power(xj, 22LL)*

                                  (-628425LL + 51450LL*Power(r, 2LL)*Power(xj, 2LL) -

            735LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  378LL*Power(xi, 6LL)*Power(xj, 22LL)*

                                  (-149625LL + 37800LL*Power(r, 2LL)*Power(xj, 2LL) -

            945LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                                  18LL*Power(r, 3LL)*Power(xi, 31LL)*

                                  (17093160LL + 6309387LL*Power(r, 2LL)*Power(xj, 2LL) -

            23562LL*Power(r, 4LL)*Power(xj, 4LL) + 16LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  54LL*Power(r, 2LL)*Power(xi, 30LL)*

                                  (-26860680LL - 24843735LL*Power(r, 2LL)*Power(xj, 2LL) -

            40180LL*Power(r, 4LL)*Power(xj, 4LL) + 578LL*Power(r, 6LL)*Power(xj, 6LL)) +

                                  378LL*r*Power(xi, 23LL)*Power(xj, 6LL)*

                                  (-14625683325LL + 704051250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10752861LL*Power(r, 4LL)*Power(xj, 4LL) + 33478LL*Power(r, 6LL)*Power(xj, 6LL))

                                  + 3LL*r*Power(xi, 9LL)*Power(xj, 20LL)*

                                  (152707275LL - 17595900LL*Power(r, 2LL)*Power(xj, 2LL) +

            396900LL*Power(r, 4LL)*Power(xj, 4LL) - 2268LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  27LL*Power(xi, 8LL)*Power(xj, 20LL)*

                                  (9426375LL - 3351600LL*Power(r, 2LL)*Power(xj, 2LL) +

            132300LL*Power(r, 4LL)*Power(xj, 4LL) - 1176LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  567LL*Power(xi, 10LL)*Power(xj, 18LL)*

                                  (1526175LL - 718200LL*Power(r, 2LL)*Power(xj, 2LL) +

            39900LL*Power(r, 4LL)*Power(xj, 4LL) - 560LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  54LL*r*Power(xi, 13LL)*Power(xj, 16LL)*

                                  (-1356769575LL - 127011675LL*Power(r, 2LL)*Power(xj, 2LL) -

            3867843LL*Power(r, 4LL)*Power(xj, 4LL) - 8556LL*Power(r, 6LL)*Power(xj, 6LL) +

            7LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  7LL*r*Power(xi, 11LL)*Power(xj, 18LL)*

                                  (-151091325LL + 45272250LL*Power(r, 2LL)*Power(xj, 2LL) -

            647676LL*Power(r, 4LL)*Power(xj, 4LL) + 15336LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  18LL*r*Power(xi, 15LL)*Power(xj, 14LL)*

                                  (63046289250LL + 3917182500LL*Power(r, 2LL)*Power(xj, 2LL) +

            10158435LL*Power(r, 4LL)*Power(xj, 4LL) -

            178842LL*Power(r, 6LL)*Power(xj, 6LL) + 16LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  378LL*r*Power(xi, 21LL)*Power(xj, 8LL)*

                                  (-8559820125LL + 17573325LL*Power(r, 2LL)*Power(xj, 2LL) +

            7421001LL*Power(r, 4LL)*Power(xj, 4LL) -

            49096LL*Power(r, 6LL)*Power(xj, 6LL) + 19LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  378LL*Power(xi, 12LL)*Power(xj, 16LL)*

                                  (17296650LL + 14244300LL*Power(r, 2LL)*Power(xj, 2LL) +

            360525LL*Power(r, 4LL)*Power(xj, 4LL) + 15928LL*Power(r, 6LL)*Power(xj, 6LL) +

            22LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  189LL*r*Power(xi, 25LL)*Power(xj, 4LL)*

                                  (9994948425LL + 63821700LL*Power(r, 2LL)*Power(xj, 2LL) -

            1458540LL*Power(r, 4LL)*Power(xj, 4LL) -

            18756LL*Power(r, 6LL)*Power(xj, 6LL) + 38LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  189LL*Power(xi, 24LL)*Power(xj, 4LL)*

                                  (17962854525LL + 4036942800LL*Power(r, 2LL)*Power(xj, 2LL) -

            126472500LL*Power(r, 4LL)*Power(xj, 4LL) +

            765464LL*Power(r, 6LL)*Power(xj, 6LL) + 190LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  21LL*r*Power(xi, 19LL)*Power(xj, 10LL)*

                                  (-228066210225LL + 13487616450LL*Power(r, 2LL)*Power(xj, 2LL) -

            85465800LL*Power(r, 4LL)*Power(xj, 4LL) -

            320112LL*Power(r, 6LL)*Power(xj, 6LL) + 328LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  189LL*Power(xi, 18LL)*Power(xj, 10LL)*

                                  (86069971575LL + 2157712200LL*Power(r, 2LL)*Power(xj, 2LL) -

            158179560LL*Power(r, 4LL)*Power(xj, 4LL) +

            578816LL*Power(r, 6LL)*Power(xj, 6LL) + 978LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  2LL*r*Power(xi, 29LL)*(2085060285LL +

                                                         5450330025LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                         127424745LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                         1398276LL*Power(r, 6LL)*Power(xj, 6LL) + 1159LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  378LL*Power(xi, 22LL)*Power(xj, 6LL)*

                                  (37244490525LL - 2411839800LL*Power(r, 2LL)*Power(xj, 2LL) +

            92951775LL*Power(r, 4LL)*Power(xj, 4LL) -

            942172LL*Power(r, 6LL)*Power(xj, 6LL) + 1292LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  27LL*Power(xi, 16LL)*Power(xj, 12LL)*

                                  (164245367475LL + 26909517600LL*Power(r, 2LL)*Power(xj, 2LL) +

            62674920LL*Power(r, 4LL)*Power(xj, 4LL) -

            3885112LL*Power(r, 6LL)*Power(xj, 6LL) + 2122LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  3LL*r*Power(xi, 27LL)*Power(xj, 2LL)*

                                  (-63819198135LL - 21841975890LL*Power(r, 2LL)*Power(xj, 2LL) +

            442430100LL*Power(r, 4LL)*Power(xj, 4LL) -

            2756664LL*Power(r, 6LL)*Power(xj, 6LL) + 2296LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  r*Power(xi, 17LL)*Power(xj, 12LL)*

                                  (4851990871875LL + 21622847400LL*Power(r, 2LL)*Power(xj, 2LL) -

            2153738160LL*Power(r, 4LL)*Power(xj, 4LL) +

            3608388LL*Power(r, 6LL)*Power(xj, 6LL) + 2318LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  18LL*Power(xi, 14LL)*Power(xj, 14LL)*

                                  (-23418646650LL - 6922729800LL*Power(r, 2LL)*Power(xj, 2LL) -

            259958475LL*Power(r, 4LL)*Power(xj, 4LL) -

            697732LL*Power(r, 6LL)*Power(xj, 6LL) + 3030LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  126LL*Power(xi, 20LL)*Power(xj, 8LL)*

                                  (-186637212225LL + 13028280300LL*Power(r, 2LL)*Power(xj, 2LL) -

            116198775LL*Power(r, 4LL)*Power(xj, 4LL) -

            1266160LL*Power(r, 6LL)*Power(xj, 6LL) + 4332LL*Power(r, 8LL)*Power(xj, 8LL)) -

                                  54LL*Power(xi, 28LL)*(102965940LL + 1089437580LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        102508245LL*Power(r, 4LL)*Power(xj, 4LL) -

                                                        1593144LL*Power(r, 6LL)*Power(xj, 6LL) + 4538LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  63LL*Power(xi, 26LL)*Power(xj, 2LL)*

                                  (-4544129205LL - 7396000920LL*Power(r, 2LL)*Power(xj, 2LL) +

            149614020LL*Power(r, 4LL)*Power(xj, 4LL) -

            1684112LL*Power(r, 6LL)*Power(xj, 6LL) + 5922LL*Power(r, 8LL)*Power(xj, 8LL))) +

                                          exp(2LL*r*xi)*Power(xi, 12LL)*

                                          (6LL*Power(xi, 24LL)*Power(xj, 6LL)*

                                  (1036901250LL + 1900985625LL*r*xj +

            1728168750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1036901250LL*Power(r, 3LL)*Power(xj, 3LL) +

            460845000LL*Power(r, 4LL)*Power(xj, 4LL) +

            161295750LL*Power(r, 5LL)*Power(xj, 5LL) +

            46084500LL*Power(r, 6LL)*Power(xj, 6LL) +

            9084900LL*Power(r, 7LL)*Power(xj, 7LL) +

            4082100LL*Power(r, 8LL)*Power(xj, 8LL) +

            121935LL*Power(r, 9LL)*Power(xj, 9LL) -

            21494LL*Power(r, 10LL)*Power(xj, 10LL) - 766LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  5LL*Power(xi, 28LL)*Power(xj, 2LL)*

                                  (19646550LL + 36018675LL*r*xj + 32744250LL*Power(r, 2LL)*Power(xj, 2LL) +

            19646550LL*Power(r, 3LL)*Power(xj, 3LL) +

            8731800LL*Power(r, 4LL)*Power(xj, 4LL) +

            3056130LL*Power(r, 5LL)*Power(xj, 5LL) +

            873180LL*Power(r, 6LL)*Power(xj, 6LL) + 207900LL*Power(r, 7LL)*Power(xj, 7LL) +

            41580LL*Power(r, 8LL)*Power(xj, 8LL) + 6930LL*Power(r, 9LL)*Power(xj, 9LL) +

            924LL*Power(r, 10LL)*Power(xj, 10LL) - 4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  26334LL*Power(xi, 16LL)*Power(xj, 14LL)*

                                  (43880400LL - 186686775LL*r*xj + 576771750LL*Power(r, 2LL)*Power(xj, 2LL) -

            398603250LL*Power(r, 3LL)*Power(xj, 3LL) +

            72552600LL*Power(r, 4LL)*Power(xj, 4LL) +

            27903120LL*Power(r, 5LL)*Power(xj, 5LL) -

            342720LL*Power(r, 6LL)*Power(xj, 6LL) - 574800LL*Power(r, 7LL)*Power(xj, 7LL) -

            50800LL*Power(r, 8LL)*Power(xj, 8LL) - 945LL*Power(r, 9LL)*Power(xj, 9LL) +

            58LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  10LL*Power(xj, 30LL)*(97302813300LL + 89194245525LL*r*xj +

                                                        38780106750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                        10576392750LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                        2014551000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                        282037140LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                        29688120LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                        2356200LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                        138600LL*Power(r, 8LL)*Power(xj, 8LL) + 5775LL*Power(r, 9LL)*Power(xj, 9LL) +

                                                        154LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  10LL*Power(xi, 2LL)*Power(xj, 28LL)*

                                  (4582499159700LL + 3733416276975LL*r*xj +

            1428215931450LL*Power(r, 2LL)*Power(xj, 2LL) +

            338545295550LL*Power(r, 3LL)*Power(xj, 3LL) +

            55198697400LL*Power(r, 4LL)*Power(xj, 4LL) +

            6486854220LL*Power(r, 5LL)*Power(xj, 5LL) +

            558419400LL*Power(r, 6LL)*Power(xj, 6LL) +

            34939080LL*Power(r, 7LL)*Power(xj, 7LL) +

            1532520LL*Power(r, 8LL)*Power(xj, 8LL) + 43285LL*Power(r, 9LL)*Power(xj, 9LL) +

            638LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  110LL*Power(xi, 10LL)*Power(xj, 20LL)*

                                  (-14063418170550LL + 6795156458475LL*r*xj +

            2067471236250LL*Power(r, 2LL)*Power(xj, 2LL) -

            214664924250LL*Power(r, 3LL)*Power(xj, 3LL) -

            124416469800LL*Power(r, 4LL)*Power(xj, 4LL) -

            14935545450LL*Power(r, 5LL)*Power(xj, 5LL) -

            256688460LL*Power(r, 6LL)*Power(xj, 6LL) +

            105750900LL*Power(r, 7LL)*Power(xj, 7LL) +

            11502180LL*Power(r, 8LL)*Power(xj, 8LL) +

            518085LL*Power(r, 9LL)*Power(xj, 9LL) + 9294LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  55LL*Power(xi, 20LL)*Power(xj, 10LL)*

                                  (1730682450LL + 3172917825LL*r*xj +

            2884470750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1571960250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1404081000LL*Power(r, 4LL)*Power(xj, 4LL) -

            426654270LL*Power(r, 5LL)*Power(xj, 5LL) +

            283536540LL*Power(r, 6LL)*Power(xj, 6LL) +

            39116700LL*Power(r, 7LL)*Power(xj, 7LL) -

            2659860LL*Power(r, 8LL)*Power(xj, 8LL) -

            528850LL*Power(r, 9LL)*Power(xj, 9LL) -

            18236LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  5LL*Power(xi, 30LL)*(935550LL + 1715175LL*r*xj +

                                                       1559250LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                       935550LL*Power(r, 3LL)*Power(xj, 3LL) + 415800LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                       145530LL*Power(r, 5LL)*Power(xj, 5LL) + 41580LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                       9900LL*Power(r, 7LL)*Power(xj, 7LL) + 1980LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                       330LL*Power(r, 9LL)*Power(xj, 9LL) + 44LL*Power(r, 10LL)*Power(xj, 10LL) +

                                                       4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  13167LL*Power(xi, 14LL)*Power(xj, 16LL)*

                                  (-2319354450LL + 8540029575LL*r*xj -

            7335672750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1133154750LL*Power(r, 3LL)*Power(xj, 3LL) +

            575014200LL*Power(r, 4LL)*Power(xj, 4LL) -

            913710LL*Power(r, 5LL)*Power(xj, 5LL) -

            14863940LL*Power(r, 6LL)*Power(xj, 6LL) -

            1687300LL*Power(r, 7LL)*Power(xj, 7LL) - 46900LL*Power(r, 8LL)*Power(xj, 8LL) +

            3210LL*Power(r, 9LL)*Power(xj, 9LL) + 236LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  770LL*Power(xi, 18LL)*Power(xj, 12LL)*

                                  (329653800LL + 654729075LL*r*xj + 45785250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1602483750LL*Power(r, 3LL)*Power(xj, 3LL) -

            915705000LL*Power(r, 4LL)*Power(xj, 4LL) +

            266036400LL*Power(r, 5LL)*Power(xj, 5LL) +

            63745920LL*Power(r, 6LL)*Power(xj, 6LL) -

            2304000LL*Power(r, 7LL)*Power(xj, 7LL) -

            1074240LL*Power(r, 8LL)*Power(xj, 8LL) - 64635LL*Power(r, 9LL)*Power(xj, 9LL) -

            514LL*Power(r, 10LL)*Power(xj, 10LL) + 34LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  385LL*Power(xi, 12LL)*Power(xj, 18LL)*

                                  (973565393850LL - 1429122323475LL*r*xj +

            298281831750LL*Power(r, 2LL)*Power(xj, 2LL) +

            138841148250LL*Power(r, 3LL)*Power(xj, 3LL) -

            2454240600LL*Power(r, 4LL)*Power(xj, 4LL) -

            4925394810LL*Power(r, 5LL)*Power(xj, 5LL) -

            623832300LL*Power(r, 6LL)*Power(xj, 6LL) -

            19098540LL*Power(r, 7LL)*Power(xj, 7LL) +

            2083140LL*Power(r, 8LL)*Power(xj, 8LL) +

            212430LL*Power(r, 9LL)*Power(xj, 9LL) + 7012LL*Power(r, 10LL)*Power(xj, 10LL) +

            68LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  14LL*Power(xi, 26LL)*Power(xj, 4LL)*

                                  (-70166250LL - 128638125LL*r*xj - 116943750LL*Power(r, 2LL)*Power(xj, 2LL) -

            70166250LL*Power(r, 3LL)*Power(xj, 3LL) -

            31185000LL*Power(r, 4LL)*Power(xj, 4LL) -

            10914750LL*Power(r, 5LL)*Power(xj, 5LL) -

            3118500LL*Power(r, 6LL)*Power(xj, 6LL) -

            742500LL*Power(r, 7LL)*Power(xj, 7LL) - 148500LL*Power(r, 8LL)*Power(xj, 8LL) -

            32615LL*Power(r, 9LL)*Power(xj, 9LL) - 154LL*Power(r, 10LL)*Power(xj, 10LL) +

            74LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  7LL*Power(xi, 4LL)*Power(xj, 26LL)*

                                  (-69822945249750LL - 46669577290875LL*r*xj -

            14025037430250LL*Power(r, 2LL)*Power(xj, 2LL) -

            2430881664750LL*Power(r, 3LL)*Power(xj, 3LL) -

            251629270200LL*Power(r, 4LL)*Power(xj, 4LL) -

            12434519790LL*Power(r, 5LL)*Power(xj, 5LL) +

            452930940LL*Power(r, 6LL)*Power(xj, 6LL) +

            131125500LL*Power(r, 7LL)*Power(xj, 7LL) +

            11018700LL*Power(r, 8LL)*Power(xj, 8LL) +

            514470LL*Power(r, 9LL)*Power(xj, 9LL) +

            13332LL*Power(r, 10LL)*Power(xj, 10LL) + 148LL*Power(r, 11LL)*Power(xj, 11LL)) -

                                  50LL*Power(xi, 8LL)*Power(xj, 22LL)*

                                  (-51768833574150LL - 5003280391725LL*r*xj +

            4493439477450LL*Power(r, 2LL)*Power(xj, 2LL) +

            1286866176750LL*Power(r, 3LL)*Power(xj, 3LL) +

            111437476920LL*Power(r, 4LL)*Power(xj, 4LL) -

            6620313546LL*Power(r, 5LL)*Power(xj, 5LL) -

            2406603276LL*Power(r, 6LL)*Power(xj, 6LL) -

            242686620LL*Power(r, 7LL)*Power(xj, 7LL) -

            12228876LL*Power(r, 8LL)*Power(xj, 8LL) -

            256223LL*Power(r, 9LL)*Power(xj, 9LL) + 2486LL*Power(r, 10LL)*Power(xj, 10LL) +

            158LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  25LL*Power(xi, 22LL)*Power(xj, 8LL)*

                                  (-1119853350LL - 2053064475LL*r*xj -

            1866422250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1119853350LL*Power(r, 3LL)*Power(xj, 3LL) -

            497712600LL*Power(r, 4LL)*Power(xj, 4LL) -

            194415606LL*Power(r, 5LL)*Power(xj, 5LL) -

            9338868LL*Power(r, 6LL)*Power(xj, 6LL) -

            31217076LL*Power(r, 7LL)*Power(xj, 7LL) -

            2256804LL*Power(r, 8LL)*Power(xj, 8LL) +

            246774LL*Power(r, 9LL)*Power(xj, 9LL) +

            22836LL*Power(r, 10LL)*Power(xj, 10LL) + 316LL*Power(r, 11LL)*Power(xj, 11LL)) +

                                  3LL*Power(xi, 6LL)*Power(xj, 24LL)*

                                  (596006592662250LL + 266778699697125LL*r*xj +

            37515651153750LL*Power(r, 2LL)*Power(xj, 2LL) -

            2214626163750LL*Power(r, 3LL)*Power(xj, 3LL) -

            1538075107800LL*Power(r, 4LL)*Power(xj, 4LL) -

            248955308910LL*Power(r, 5LL)*Power(xj, 5LL) -

            21434337540LL*Power(r, 6LL)*Power(xj, 6LL) -

            957980100LL*Power(r, 7LL)*Power(xj, 7LL) -

            4874100LL*Power(r, 8LL)*Power(xj, 8LL) +

            1831830LL*Power(r, 9LL)*Power(xj, 9LL) +

            91828LL*Power(r, 10LL)*Power(xj, 10LL) + 1532LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (2.338875e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 21LL)*Power(xi + xj, 20LL))

                - (9355500LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                   Power(Power(xi, 2LL) - Power(xj, 2LL), 21LL) +

                   110LL*exp(2LL*r*xj)*Power(xj, 14LL)*

                   (-3456LL*Power(r, 7LL)*Power(xi, 36LL) - 54LL*Power(r, 8LL)*Power(xi, 37LL) -

                    11340LL*Power(r, 7LL)*Power(xi, 34LL)*Power(xj, 2LL) -

                    112LL*Power(r, 8LL)*Power(xi, 35LL)*Power(xj, 2LL) + 76545LL*xi*Power(xj, 28LL) +

                    136080LL*r*Power(xi, 2LL)*Power(xj, 28LL) +

                    79380LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 28LL) +

                    19845LL*Power(xi, 3LL)*Power(xj, 26LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    6804LL*Power(r, 5LL)*Power(xi, 34LL)*(272LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    56LL*Power(r, 6LL)*Power(xi, 35LL)*(1836LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    378LL*Power(r, 5LL)*Power(xi, 33LL)*

                    (-1332LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    5670LL*r*Power(xi, 5LL)*Power(xj, 24LL)*

                    (-294LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    17010LL*Power(xi, 4LL)*Power(xj, 24LL)*

                    (-168LL*r*Power(xj, 2LL) + 4LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    378LL*Power(r, 4LL)*Power(xi, 32LL)*

                    (-34888LL*r*Power(xj, 2LL) + 236LL*Power(r, 3LL)*Power(xj, 4LL)) +

                    1890LL*Power(r, 4LL)*Power(xi, 33LL)*

                    (-11628LL - 666LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    5670LL*Power(xi, 5LL)*Power(xj, 24LL)*

                    (2835LL - 147LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    1512LL*Power(r, 3LL)*Power(xi, 32LL)*

                    (-116280LL - 17444LL*Power(r, 2LL)*Power(xj, 2LL) +

            59LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    162LL*r*Power(xi, 7LL)*Power(xj, 22LL)*

                    (102900LL*r*Power(xj, 2LL) - 2940LL*Power(r, 3LL)*Power(xj, 4LL) +

            12LL*Power(r, 5LL)*Power(xj, 6LL)) +

                    378LL*Power(xi, 6LL)*Power(xj, 22LL)*

                    (75600LL*r*Power(xj, 2LL) - 3780LL*Power(r, 3LL)*Power(xj, 4LL) +

            24LL*Power(r, 5LL)*Power(xj, 6LL)) -

                    18LL*Power(r, 3LL)*Power(xi, 31LL)*

                    (12618774LL*r*Power(xj, 2LL) - 94248LL*Power(r, 3LL)*Power(xj, 4LL) +

            96LL*Power(r, 5LL)*Power(xj, 6LL)) +

                    54LL*Power(r, 2LL)*Power(xi, 30LL)*

                    (-49687470LL*r*Power(xj, 2LL) - 160720LL*Power(r, 3LL)*Power(xj, 4LL) +

            3468LL*Power(r, 5LL)*Power(xj, 6LL)) +

                    378LL*r*Power(xi, 23LL)*Power(xj, 6LL)*

                    (1408102500LL*r*Power(xj, 2LL) - 43011444LL*Power(r, 3LL)*Power(xj, 4LL) +

            200868LL*Power(r, 5LL)*Power(xj, 6LL)) +

                    162LL*Power(xi, 7LL)*Power(xj, 22LL)*

                    (-628425LL + 51450LL*Power(r, 2LL)*Power(xj, 2LL) -

            735LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) -

                    54LL*Power(r, 2LL)*Power(xi, 31LL)*

                    (17093160LL + 6309387LL*Power(r, 2LL)*Power(xj, 2LL) -

            23562LL*Power(r, 4LL)*Power(xj, 4LL) + 16LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    108LL*r*Power(xi, 30LL)*(-26860680LL - 24843735LL*Power(r, 2LL)*Power(xj, 2LL) -

                                             40180LL*Power(r, 4LL)*Power(xj, 4LL) + 578LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    378LL*Power(xi, 23LL)*Power(xj, 6LL)*

                    (-14625683325LL + 704051250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10752861LL*Power(r, 4LL)*Power(xj, 4LL) + 33478LL*Power(r, 6LL)*Power(xj, 6LL))

                    + 3LL*r*Power(xi, 9LL)*Power(xj, 20LL)*

                    (-35191800LL*r*Power(xj, 2LL) + 1587600LL*Power(r, 3LL)*Power(xj, 4LL) -

            13608LL*Power(r, 5LL)*Power(xj, 6LL) + 16LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    27LL*Power(xi, 8LL)*Power(xj, 20LL)*

                    (-6703200LL*r*Power(xj, 2LL) + 529200LL*Power(r, 3LL)*Power(xj, 4LL) -

            7056LL*Power(r, 5LL)*Power(xj, 6LL) + 16LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    567LL*Power(xi, 10LL)*Power(xj, 18LL)*

                    (-1436400LL*r*Power(xj, 2LL) + 159600LL*Power(r, 3LL)*Power(xj, 4LL) -

            3360LL*Power(r, 5LL)*Power(xj, 6LL) + 16LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    54LL*r*Power(xi, 13LL)*Power(xj, 16LL)*

                    (-254023350LL*r*Power(xj, 2LL) - 15471372LL*Power(r, 3LL)*Power(xj, 4LL) -

            51336LL*Power(r, 5LL)*Power(xj, 6LL) + 56LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    7LL*r*Power(xi, 11LL)*Power(xj, 18LL)*

                    (90544500LL*r*Power(xj, 2LL) - 2590704LL*Power(r, 3LL)*Power(xj, 4LL) +

            92016LL*Power(r, 5LL)*Power(xj, 6LL) + 64LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    18LL*r*Power(xi, 15LL)*Power(xj, 14LL)*

                    (7834365000LL*r*Power(xj, 2LL) + 40633740LL*Power(r, 3LL)*Power(xj, 4LL) -

            1073052LL*Power(r, 5LL)*Power(xj, 6LL) + 128LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    378LL*r*Power(xi, 21LL)*Power(xj, 8LL)*

                    (35146650LL*r*Power(xj, 2LL) + 29684004LL*Power(r, 3LL)*Power(xj, 4LL) -

            294576LL*Power(r, 5LL)*Power(xj, 6LL) + 152LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    378LL*Power(xi, 12LL)*Power(xj, 16LL)*

                    (28488600LL*r*Power(xj, 2LL) + 1442100LL*Power(r, 3LL)*Power(xj, 4LL) +

            95568LL*Power(r, 5LL)*Power(xj, 6LL) + 176LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    189LL*r*Power(xi, 25LL)*Power(xj, 4LL)*

                    (127643400LL*r*Power(xj, 2LL) - 5834160LL*Power(r, 3LL)*Power(xj, 4LL) -

            112536LL*Power(r, 5LL)*Power(xj, 6LL) + 304LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    189LL*Power(xi, 24LL)*Power(xj, 4LL)*

                    (8073885600LL*r*Power(xj, 2LL) - 505890000LL*Power(r, 3LL)*Power(xj, 4LL) +

            4592784LL*Power(r, 5LL)*Power(xj, 6LL) + 1520LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    21LL*r*Power(xi, 19LL)*Power(xj, 10LL)*

                    (26975232900LL*r*Power(xj, 2LL) - 341863200LL*Power(r, 3LL)*Power(xj, 4LL) -

            1920672LL*Power(r, 5LL)*Power(xj, 6LL) + 2624LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    189LL*Power(xi, 18LL)*Power(xj, 10LL)*

                    (4315424400LL*r*Power(xj, 2LL) - 632718240LL*Power(r, 3LL)*Power(xj, 4LL) +

            3472896LL*Power(r, 5LL)*Power(xj, 6LL) + 7824LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    2LL*r*Power(xi, 29LL)*(10900660050LL*r*Power(xj, 2LL) +

                                           509698980LL*Power(r, 3LL)*Power(xj, 4LL) -

                                           8389656LL*Power(r, 5LL)*Power(xj, 6LL) + 9272LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    378LL*Power(xi, 22LL)*Power(xj, 6LL)*

                    (-4823679600LL*r*Power(xj, 2LL) + 371807100LL*Power(r, 3LL)*Power(xj, 4LL) -

            5653032LL*Power(r, 5LL)*Power(xj, 6LL) + 10336LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    27LL*Power(xi, 16LL)*Power(xj, 12LL)*

                    (53819035200LL*r*Power(xj, 2LL) + 250699680LL*Power(r, 3LL)*Power(xj, 4LL) -

            23310672LL*Power(r, 5LL)*Power(xj, 6LL) + 16976LL*Power(r, 7LL)*Power(xj, 8LL))

                    + 3LL*r*Power(xi, 27LL)*Power(xj, 2LL)*

                    (-43683951780LL*r*Power(xj, 2LL) + 1769720400LL*Power(r, 3LL)*Power(xj, 4LL) -

            16539984LL*Power(r, 5LL)*Power(xj, 6LL) + 18368LL*Power(r, 7LL)*Power(xj, 8LL))

                    + r*Power(xi, 17LL)*Power(xj, 12LL)*(43245694800LL*r*Power(xj, 2LL) -

                                                         8614952640LL*Power(r, 3LL)*Power(xj, 4LL) +

                                                         21650328LL*Power(r, 5LL)*Power(xj, 6LL) + 18544LL*Power(r, 7LL)*Power(xj, 8LL))

                    + 18LL*Power(xi, 14LL)*Power(xj, 14LL)*

                    (-13845459600LL*r*Power(xj, 2LL) - 1039833900LL*Power(r, 3LL)*Power(xj, 4LL) -

            4186392LL*Power(r, 5LL)*Power(xj, 6LL) + 24240LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    126LL*Power(xi, 20LL)*Power(xj, 8LL)*

                    (26056560600LL*r*Power(xj, 2LL) - 464795100LL*Power(r, 3LL)*Power(xj, 4LL) -

            7596960LL*Power(r, 5LL)*Power(xj, 6LL) + 34656LL*Power(r, 7LL)*Power(xj, 8LL)) -

                    54LL*Power(xi, 28LL)*(2178875160LL*r*Power(xj, 2LL) +

                                          410032980LL*Power(r, 3LL)*Power(xj, 4LL) -

                                          9558864LL*Power(r, 5LL)*Power(xj, 6LL) + 36304LL*Power(r, 7LL)*Power(xj, 8LL)) +

                    63LL*Power(xi, 26LL)*Power(xj, 2LL)*

                    (-14792001840LL*r*Power(xj, 2LL) + 598456080LL*Power(r, 3LL)*Power(xj, 4LL) -

            10104672LL*Power(r, 5LL)*Power(xj, 6LL) + 47376LL*Power(r, 7LL)*Power(xj, 8LL))

                    + 3LL*Power(xi, 9LL)*Power(xj, 20LL)*(152707275LL -

                                                          17595900LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                          396900LL*Power(r, 4LL)*Power(xj, 4LL) - 2268LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                          2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    54LL*Power(xi, 13LL)*Power(xj, 16LL)*

                    (-1356769575LL - 127011675LL*Power(r, 2LL)*Power(xj, 2LL) -

            3867843LL*Power(r, 4LL)*Power(xj, 4LL) - 8556LL*Power(r, 6LL)*Power(xj, 6LL) +

            7LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    7LL*Power(xi, 11LL)*Power(xj, 18LL)*

                    (-151091325LL + 45272250LL*Power(r, 2LL)*Power(xj, 2LL) -

            647676LL*Power(r, 4LL)*Power(xj, 4LL) + 15336LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    18LL*Power(xi, 15LL)*Power(xj, 14LL)*

                    (63046289250LL + 3917182500LL*Power(r, 2LL)*Power(xj, 2LL) +

            10158435LL*Power(r, 4LL)*Power(xj, 4LL) -

            178842LL*Power(r, 6LL)*Power(xj, 6LL) + 16LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    378LL*Power(xi, 21LL)*Power(xj, 8LL)*

                    (-8559820125LL + 17573325LL*Power(r, 2LL)*Power(xj, 2LL) +

            7421001LL*Power(r, 4LL)*Power(xj, 4LL) - 49096LL*Power(r, 6LL)*Power(xj, 6LL) +

            19LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    189LL*Power(xi, 25LL)*Power(xj, 4LL)*

                    (9994948425LL + 63821700LL*Power(r, 2LL)*Power(xj, 2LL) -

            1458540LL*Power(r, 4LL)*Power(xj, 4LL) - 18756LL*Power(r, 6LL)*Power(xj, 6LL) +

            38LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    21LL*Power(xi, 19LL)*Power(xj, 10LL)*

                    (-228066210225LL + 13487616450LL*Power(r, 2LL)*Power(xj, 2LL) -

            85465800LL*Power(r, 4LL)*Power(xj, 4LL) -

            320112LL*Power(r, 6LL)*Power(xj, 6LL) + 328LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    2LL*Power(xi, 29LL)*(2085060285LL + 5450330025LL*Power(r, 2LL)*Power(xj, 2LL) +

                                         127424745LL*Power(r, 4LL)*Power(xj, 4LL) -

                                         1398276LL*Power(r, 6LL)*Power(xj, 6LL) + 1159LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    3LL*Power(xi, 27LL)*Power(xj, 2LL)*

                    (-63819198135LL - 21841975890LL*Power(r, 2LL)*Power(xj, 2LL) +

            442430100LL*Power(r, 4LL)*Power(xj, 4LL) -

            2756664LL*Power(r, 6LL)*Power(xj, 6LL) + 2296LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    Power(xi, 17LL)*Power(xj, 12LL)*

                    (4851990871875LL + 21622847400LL*Power(r, 2LL)*Power(xj, 2LL) -

            2153738160LL*Power(r, 4LL)*Power(xj, 4LL) +

            3608388LL*Power(r, 6LL)*Power(xj, 6LL) + 2318LL*Power(r, 8LL)*Power(xj, 8LL))) +

                   220LL*exp(2LL*r*xj)*Power(xj, 15LL)*

                   (-432LL*Power(r, 8LL)*Power(xi, 36LL) - 6LL*Power(r, 9LL)*Power(xi, 37LL) +

                    42525LL*Power(xj, 28LL) + 76545LL*r*xi*Power(xj, 28LL) +

                    19845LL*r*Power(xi, 3LL)*Power(xj, 26LL)*(-81LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    1134LL*Power(r, 6LL)*Power(xi, 34LL)*(272LL + 5LL*Power(r, 2LL)*Power(xj, 2LL)) -

                    8LL*Power(r, 7LL)*Power(xi, 35LL)*(1836LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    8505LL*Power(xi, 2LL)*Power(xj, 26LL)*(-105LL + 8LL*Power(r, 2LL)*Power(xj, 2LL)) +

                    378LL*Power(r, 5LL)*Power(xi, 33LL)*

                    (-11628LL - 666LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    5670LL*r*Power(xi, 5LL)*Power(xj, 24LL)*

                    (2835LL - 147LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    17010LL*Power(xi, 4LL)*Power(xj, 24LL)*

                    (525LL - 84LL*Power(r, 2LL)*Power(xj, 2LL) + Power(r, 4LL)*Power(xj, 4LL)) +

                    378LL*Power(r, 4LL)*Power(xi, 32LL)*

                    (-116280LL - 17444LL*Power(r, 2LL)*Power(xj, 2LL) +

            59LL*Power(r, 4LL)*Power(xj, 4LL)) +

                    162LL*r*Power(xi, 7LL)*Power(xj, 22LL)*

                    (-628425LL + 51450LL*Power(r, 2LL)*Power(xj, 2LL) -

            735LL*Power(r, 4LL)*Power(xj, 4LL) + 2LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    378LL*Power(xi, 6LL)*Power(xj, 22LL)*

                    (-149625LL + 37800LL*Power(r, 2LL)*Power(xj, 2LL) -

            945LL*Power(r, 4LL)*Power(xj, 4LL) + 4LL*Power(r, 6LL)*Power(xj, 6LL)) -

                    18LL*Power(r, 3LL)*Power(xi, 31LL)*

                    (17093160LL + 6309387LL*Power(r, 2LL)*Power(xj, 2LL) -

            23562LL*Power(r, 4LL)*Power(xj, 4LL) + 16LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    54LL*Power(r, 2LL)*Power(xi, 30LL)*

                    (-26860680LL - 24843735LL*Power(r, 2LL)*Power(xj, 2LL) -

            40180LL*Power(r, 4LL)*Power(xj, 4LL) + 578LL*Power(r, 6LL)*Power(xj, 6LL)) +

                    378LL*r*Power(xi, 23LL)*Power(xj, 6LL)*

                    (-14625683325LL + 704051250LL*Power(r, 2LL)*Power(xj, 2LL) -

            10752861LL*Power(r, 4LL)*Power(xj, 4LL) + 33478LL*Power(r, 6LL)*Power(xj, 6LL))

                    + 3LL*r*Power(xi, 9LL)*Power(xj, 20LL)*

                    (152707275LL - 17595900LL*Power(r, 2LL)*Power(xj, 2LL) +

            396900LL*Power(r, 4LL)*Power(xj, 4LL) - 2268LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    27LL*Power(xi, 8LL)*Power(xj, 20LL)*

                    (9426375LL - 3351600LL*Power(r, 2LL)*Power(xj, 2LL) +

            132300LL*Power(r, 4LL)*Power(xj, 4LL) - 1176LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    567LL*Power(xi, 10LL)*Power(xj, 18LL)*

                    (1526175LL - 718200LL*Power(r, 2LL)*Power(xj, 2LL) +

            39900LL*Power(r, 4LL)*Power(xj, 4LL) - 560LL*Power(r, 6LL)*Power(xj, 6LL) +

            2LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    54LL*r*Power(xi, 13LL)*Power(xj, 16LL)*

                    (-1356769575LL - 127011675LL*Power(r, 2LL)*Power(xj, 2LL) -

            3867843LL*Power(r, 4LL)*Power(xj, 4LL) - 8556LL*Power(r, 6LL)*Power(xj, 6LL) +

            7LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    7LL*r*Power(xi, 11LL)*Power(xj, 18LL)*

                    (-151091325LL + 45272250LL*Power(r, 2LL)*Power(xj, 2LL) -

            647676LL*Power(r, 4LL)*Power(xj, 4LL) + 15336LL*Power(r, 6LL)*Power(xj, 6LL) +

            8LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    18LL*r*Power(xi, 15LL)*Power(xj, 14LL)*

                    (63046289250LL + 3917182500LL*Power(r, 2LL)*Power(xj, 2LL) +

            10158435LL*Power(r, 4LL)*Power(xj, 4LL) -

            178842LL*Power(r, 6LL)*Power(xj, 6LL) + 16LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    378LL*r*Power(xi, 21LL)*Power(xj, 8LL)*

                    (-8559820125LL + 17573325LL*Power(r, 2LL)*Power(xj, 2LL) +

            7421001LL*Power(r, 4LL)*Power(xj, 4LL) - 49096LL*Power(r, 6LL)*Power(xj, 6LL) +

            19LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    378LL*Power(xi, 12LL)*Power(xj, 16LL)*

                    (17296650LL + 14244300LL*Power(r, 2LL)*Power(xj, 2LL) +

            360525LL*Power(r, 4LL)*Power(xj, 4LL) + 15928LL*Power(r, 6LL)*Power(xj, 6LL) +

            22LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    189LL*r*Power(xi, 25LL)*Power(xj, 4LL)*

                    (9994948425LL + 63821700LL*Power(r, 2LL)*Power(xj, 2LL) -

            1458540LL*Power(r, 4LL)*Power(xj, 4LL) - 18756LL*Power(r, 6LL)*Power(xj, 6LL) +

            38LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    189LL*Power(xi, 24LL)*Power(xj, 4LL)*

                    (17962854525LL + 4036942800LL*Power(r, 2LL)*Power(xj, 2LL) -

            126472500LL*Power(r, 4LL)*Power(xj, 4LL) +

            765464LL*Power(r, 6LL)*Power(xj, 6LL) + 190LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    21LL*r*Power(xi, 19LL)*Power(xj, 10LL)*

                    (-228066210225LL + 13487616450LL*Power(r, 2LL)*Power(xj, 2LL) -

            85465800LL*Power(r, 4LL)*Power(xj, 4LL) -

            320112LL*Power(r, 6LL)*Power(xj, 6LL) + 328LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    189LL*Power(xi, 18LL)*Power(xj, 10LL)*

                    (86069971575LL + 2157712200LL*Power(r, 2LL)*Power(xj, 2LL) -

            158179560LL*Power(r, 4LL)*Power(xj, 4LL) +

            578816LL*Power(r, 6LL)*Power(xj, 6LL) + 978LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    2LL*r*Power(xi, 29LL)*(2085060285LL + 5450330025LL*Power(r, 2LL)*Power(xj, 2LL) +

                                           127424745LL*Power(r, 4LL)*Power(xj, 4LL) -

                                           1398276LL*Power(r, 6LL)*Power(xj, 6LL) + 1159LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    378LL*Power(xi, 22LL)*Power(xj, 6LL)*

                    (37244490525LL - 2411839800LL*Power(r, 2LL)*Power(xj, 2LL) +

            92951775LL*Power(r, 4LL)*Power(xj, 4LL) -

            942172LL*Power(r, 6LL)*Power(xj, 6LL) + 1292LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    27LL*Power(xi, 16LL)*Power(xj, 12LL)*

                    (164245367475LL + 26909517600LL*Power(r, 2LL)*Power(xj, 2LL) +

            62674920LL*Power(r, 4LL)*Power(xj, 4LL) -

            3885112LL*Power(r, 6LL)*Power(xj, 6LL) + 2122LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    3LL*r*Power(xi, 27LL)*Power(xj, 2LL)*

                    (-63819198135LL - 21841975890LL*Power(r, 2LL)*Power(xj, 2LL) +

            442430100LL*Power(r, 4LL)*Power(xj, 4LL) -

            2756664LL*Power(r, 6LL)*Power(xj, 6LL) + 2296LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    r*Power(xi, 17LL)*Power(xj, 12LL)*

                    (4851990871875LL + 21622847400LL*Power(r, 2LL)*Power(xj, 2LL) -

            2153738160LL*Power(r, 4LL)*Power(xj, 4LL) +

            3608388LL*Power(r, 6LL)*Power(xj, 6LL) + 2318LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    18LL*Power(xi, 14LL)*Power(xj, 14LL)*

                    (-23418646650LL - 6922729800LL*Power(r, 2LL)*Power(xj, 2LL) -

            259958475LL*Power(r, 4LL)*Power(xj, 4LL) -

            697732LL*Power(r, 6LL)*Power(xj, 6LL) + 3030LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    126LL*Power(xi, 20LL)*Power(xj, 8LL)*

                    (-186637212225LL + 13028280300LL*Power(r, 2LL)*Power(xj, 2LL) -

            116198775LL*Power(r, 4LL)*Power(xj, 4LL) -

            1266160LL*Power(r, 6LL)*Power(xj, 6LL) + 4332LL*Power(r, 8LL)*Power(xj, 8LL)) -

                    54LL*Power(xi, 28LL)*(102965940LL + 1089437580LL*Power(r, 2LL)*Power(xj, 2LL) +

                                          102508245LL*Power(r, 4LL)*Power(xj, 4LL) -

                                          1593144LL*Power(r, 6LL)*Power(xj, 6LL) + 4538LL*Power(r, 8LL)*Power(xj, 8LL)) +

                    63LL*Power(xi, 26LL)*Power(xj, 2LL)*

                    (-4544129205LL - 7396000920LL*Power(r, 2LL)*Power(xj, 2LL) +

            149614020LL*Power(r, 4LL)*Power(xj, 4LL) -

            1684112LL*Power(r, 6LL)*Power(xj, 6LL) + 5922LL*Power(r, 8LL)*Power(xj, 8LL))) +

                   exp(2LL*r*xi)*Power(xi, 12LL)*

                   (6LL*Power(xi, 24LL)*Power(xj, 6LL)*

                    (1900985625LL*xj + 3456337500LL*r*Power(xj, 2LL) +

            3110703750LL*Power(r, 2LL)*Power(xj, 3LL) +

            1843380000LL*Power(r, 3LL)*Power(xj, 4LL) +

            806478750LL*Power(r, 4LL)*Power(xj, 5LL) +

            276507000LL*Power(r, 5LL)*Power(xj, 6LL) +

            63594300LL*Power(r, 6LL)*Power(xj, 7LL) +

            32656800LL*Power(r, 7LL)*Power(xj, 8LL) +

            1097415LL*Power(r, 8LL)*Power(xj, 9LL) -

            214940LL*Power(r, 9LL)*Power(xj, 10LL) - 8426LL*Power(r, 10LL)*Power(xj, 11LL))

                    + 5LL*Power(xi, 28LL)*Power(xj, 2LL)*(36018675LL*xj + 65488500LL*r*Power(xj, 2LL) +

                                                          58939650LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                          34927200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                                          15280650LL*Power(r, 4LL)*Power(xj, 5LL) +

                                                          5239080LL*Power(r, 5LL)*Power(xj, 6LL) +

                                                          1455300LL*Power(r, 6LL)*Power(xj, 7LL) +

                                                          332640LL*Power(r, 7LL)*Power(xj, 8LL) + 62370LL*Power(r, 8LL)*Power(xj, 9LL) +

                                                          9240LL*Power(r, 9LL)*Power(xj, 10LL) - 44LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    26334LL*Power(xi, 16LL)*Power(xj, 14LL)*

                    (-186686775LL*xj + 1153543500LL*r*Power(xj, 2LL) -

            1195809750LL*Power(r, 2LL)*Power(xj, 3LL) +

            290210400LL*Power(r, 3LL)*Power(xj, 4LL) +

            139515600LL*Power(r, 4LL)*Power(xj, 5LL) -

            2056320LL*Power(r, 5LL)*Power(xj, 6LL) -

            4023600LL*Power(r, 6LL)*Power(xj, 7LL) -

            406400LL*Power(r, 7LL)*Power(xj, 8LL) - 8505LL*Power(r, 8LL)*Power(xj, 9LL) +

            580LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    10LL*Power(xj, 30LL)*(89194245525LL*xj + 77560213500LL*r*Power(xj, 2LL) +

                                          31729178250LL*Power(r, 2LL)*Power(xj, 3LL) +

                                          8058204000LL*Power(r, 3LL)*Power(xj, 4LL) +

                                          1410185700LL*Power(r, 4LL)*Power(xj, 5LL) +

                                          178128720LL*Power(r, 5LL)*Power(xj, 6LL) +

                                          16493400LL*Power(r, 6LL)*Power(xj, 7LL) +

                                          1108800LL*Power(r, 7LL)*Power(xj, 8LL) + 51975LL*Power(r, 8LL)*Power(xj, 9LL) +

                                          1540LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    10LL*Power(xi, 2LL)*Power(xj, 28LL)*

                    (3733416276975LL*xj + 2856431862900LL*r*Power(xj, 2LL) +

            1015635886650LL*Power(r, 2LL)*Power(xj, 3LL) +

            220794789600LL*Power(r, 3LL)*Power(xj, 4LL) +

            32434271100LL*Power(r, 4LL)*Power(xj, 5LL) +

            3350516400LL*Power(r, 5LL)*Power(xj, 6LL) +

            244573560LL*Power(r, 6LL)*Power(xj, 7LL) +

            12260160LL*Power(r, 7LL)*Power(xj, 8LL) +

            389565LL*Power(r, 8LL)*Power(xj, 9LL) + 6380LL*Power(r, 9LL)*Power(xj, 10LL) +

            22LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    110LL*Power(xi, 10LL)*Power(xj, 20LL)*

                    (6795156458475LL*xj + 4134942472500LL*r*Power(xj, 2LL) -

            643994772750LL*Power(r, 2LL)*Power(xj, 3LL) -

            497665879200LL*Power(r, 3LL)*Power(xj, 4LL) -

            74677727250LL*Power(r, 4LL)*Power(xj, 5LL) -

            1540130760LL*Power(r, 5LL)*Power(xj, 6LL) +

            740256300LL*Power(r, 6LL)*Power(xj, 7LL) +

            92017440LL*Power(r, 7LL)*Power(xj, 8LL) +

            4662765LL*Power(r, 8LL)*Power(xj, 9LL) +

            92940LL*Power(r, 9LL)*Power(xj, 10LL) + 22LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    55LL*Power(xi, 20LL)*Power(xj, 10LL)*

                    (3172917825LL*xj + 5768941500LL*r*Power(xj, 2LL) +

            4715880750LL*Power(r, 2LL)*Power(xj, 3LL) +

            5616324000LL*Power(r, 3LL)*Power(xj, 4LL) -

            2133271350LL*Power(r, 4LL)*Power(xj, 5LL) +

            1701219240LL*Power(r, 5LL)*Power(xj, 6LL) +

            273816900LL*Power(r, 6LL)*Power(xj, 7LL) -

            21278880LL*Power(r, 7LL)*Power(xj, 8LL) -

            4759650LL*Power(r, 8LL)*Power(xj, 9LL) -

            182360LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    5LL*Power(xi, 30LL)*(1715175LL*xj + 3118500LL*r*Power(xj, 2LL) +

                                         2806650LL*Power(r, 2LL)*Power(xj, 3LL) +

                                         1663200LL*Power(r, 3LL)*Power(xj, 4LL) +

                                         727650LL*Power(r, 4LL)*Power(xj, 5LL) + 249480LL*Power(r, 5LL)*Power(xj, 6LL) +

                                         69300LL*Power(r, 6LL)*Power(xj, 7LL) + 15840LL*Power(r, 7LL)*Power(xj, 8LL) +

                                         2970LL*Power(r, 8LL)*Power(xj, 9LL) + 440LL*Power(r, 9LL)*Power(xj, 10LL) +

                                         44LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    13167LL*Power(xi, 14LL)*Power(xj, 16LL)*

                    (8540029575LL*xj - 14671345500LL*r*Power(xj, 2LL) +

            3399464250LL*Power(r, 2LL)*Power(xj, 3LL) +

            2300056800LL*Power(r, 3LL)*Power(xj, 4LL) -

            4568550LL*Power(r, 4LL)*Power(xj, 5LL) -

            89183640LL*Power(r, 5LL)*Power(xj, 6LL) -

            11811100LL*Power(r, 6LL)*Power(xj, 7LL) -

            375200LL*Power(r, 7LL)*Power(xj, 8LL) + 28890LL*Power(r, 8LL)*Power(xj, 9LL) +

            2360LL*Power(r, 9LL)*Power(xj, 10LL) + 44LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    770LL*Power(xi, 18LL)*Power(xj, 12LL)*

                    (654729075LL*xj + 91570500LL*r*Power(xj, 2LL) +

            4807451250LL*Power(r, 2LL)*Power(xj, 3LL) -

            3662820000LL*Power(r, 3LL)*Power(xj, 4LL) +

            1330182000LL*Power(r, 4LL)*Power(xj, 5LL) +

            382475520LL*Power(r, 5LL)*Power(xj, 6LL) -

            16128000LL*Power(r, 6LL)*Power(xj, 7LL) -

            8593920LL*Power(r, 7LL)*Power(xj, 8LL) -

            581715LL*Power(r, 8LL)*Power(xj, 9LL) - 5140LL*Power(r, 9LL)*Power(xj, 10LL) +

            374LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    385LL*Power(xi, 12LL)*Power(xj, 18LL)*

                    (-1429122323475LL*xj + 596563663500LL*r*Power(xj, 2LL) +

            416523444750LL*Power(r, 2LL)*Power(xj, 3LL) -

            9816962400LL*Power(r, 3LL)*Power(xj, 4LL) -

            24626974050LL*Power(r, 4LL)*Power(xj, 5LL) -

            3742993800LL*Power(r, 5LL)*Power(xj, 6LL) -

            133689780LL*Power(r, 6LL)*Power(xj, 7LL) +

            16665120LL*Power(r, 7LL)*Power(xj, 8LL) +

            1911870LL*Power(r, 8LL)*Power(xj, 9LL) +

            70120LL*Power(r, 9LL)*Power(xj, 10LL) + 748LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    14LL*Power(xi, 26LL)*Power(xj, 4LL)*

                    (-128638125LL*xj - 233887500LL*r*Power(xj, 2LL) -

            210498750LL*Power(r, 2LL)*Power(xj, 3LL) -

            124740000LL*Power(r, 3LL)*Power(xj, 4LL) -

            54573750LL*Power(r, 4LL)*Power(xj, 5LL) -

            18711000LL*Power(r, 5LL)*Power(xj, 6LL) -

            5197500LL*Power(r, 6LL)*Power(xj, 7LL) -

            1188000LL*Power(r, 7LL)*Power(xj, 8LL) -

            293535LL*Power(r, 8LL)*Power(xj, 9LL) - 1540LL*Power(r, 9LL)*Power(xj, 10LL) +

            814LL*Power(r, 10LL)*Power(xj, 11LL)) -

                    7LL*Power(xi, 4LL)*Power(xj, 26LL)*

                    (-46669577290875LL*xj - 28050074860500LL*r*Power(xj, 2LL) -

            7292644994250LL*Power(r, 2LL)*Power(xj, 3LL) -

            1006517080800LL*Power(r, 3LL)*Power(xj, 4LL) -

            62172598950LL*Power(r, 4LL)*Power(xj, 5LL) +

            2717585640LL*Power(r, 5LL)*Power(xj, 6LL) +

            917878500LL*Power(r, 6LL)*Power(xj, 7LL) +

            88149600LL*Power(r, 7LL)*Power(xj, 8LL) +

            4630230LL*Power(r, 8LL)*Power(xj, 9LL) +

            133320LL*Power(r, 9LL)*Power(xj, 10LL) + 1628LL*Power(r, 10LL)*Power(xj, 11LL))

                    - 50LL*Power(xi, 8LL)*Power(xj, 22LL)*(-5003280391725LL*xj +

                                                           8986878954900LL*r*Power(xj, 2LL) +

                                                           3860598530250LL*Power(r, 2LL)*Power(xj, 3LL) +

                                                           445749907680LL*Power(r, 3LL)*Power(xj, 4LL) -

                                                           33101567730LL*Power(r, 4LL)*Power(xj, 5LL) -

                                                           14439619656LL*Power(r, 5LL)*Power(xj, 6LL) -

                                                           1698806340LL*Power(r, 6LL)*Power(xj, 7LL) -

                                                           97831008LL*Power(r, 7LL)*Power(xj, 8LL) -

                                                           2306007LL*Power(r, 8LL)*Power(xj, 9LL) +

                                                           24860LL*Power(r, 9LL)*Power(xj, 10LL) + 1738LL*Power(r, 10LL)*Power(xj, 11LL)) +

                    25LL*Power(xi, 22LL)*Power(xj, 8LL)*

                    (-2053064475LL*xj - 3732844500LL*r*Power(xj, 2LL) -

            3359560050LL*Power(r, 2LL)*Power(xj, 3LL) -

            1990850400LL*Power(r, 3LL)*Power(xj, 4LL) -

            972078030LL*Power(r, 4LL)*Power(xj, 5LL) -

            56033208LL*Power(r, 5LL)*Power(xj, 6LL) -

            218519532LL*Power(r, 6LL)*Power(xj, 7LL) -

            18054432LL*Power(r, 7LL)*Power(xj, 8LL) +

            2220966LL*Power(r, 8LL)*Power(xj, 9LL) +

            228360LL*Power(r, 9LL)*Power(xj, 10LL) + 3476LL*Power(r, 10LL)*Power(xj, 11LL))

                    + 3LL*Power(xi, 6LL)*Power(xj, 24LL)*(266778699697125LL*xj +

                                                          75031302307500LL*r*Power(xj, 2LL) -

                                                          6643878491250LL*Power(r, 2LL)*Power(xj, 3LL) -

                                                          6152300431200LL*Power(r, 3LL)*Power(xj, 4LL) -

                                                          1244776544550LL*Power(r, 4LL)*Power(xj, 5LL) -

                                                          128606025240LL*Power(r, 5LL)*Power(xj, 6LL) -

                                                          6705860700LL*Power(r, 6LL)*Power(xj, 7LL) -

                                                          38992800LL*Power(r, 7LL)*Power(xj, 8LL) +

                                                          16486470LL*Power(r, 8LL)*Power(xj, 9LL) +

                                                          918280LL*Power(r, 9LL)*Power(xj, 10LL) + 16852LL*Power(r, 10LL)*Power(xj, 11LL)))

                   + 2LL*exp(2LL*r*xi)*Power(xi, 13LL)*

                   (6LL*Power(xi, 24LL)*Power(xj, 6LL)*

                    (1036901250LL + 1900985625LL*r*xj +

            1728168750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1036901250LL*Power(r, 3LL)*Power(xj, 3LL) +

            460845000LL*Power(r, 4LL)*Power(xj, 4LL) +

            161295750LL*Power(r, 5LL)*Power(xj, 5LL) +

            46084500LL*Power(r, 6LL)*Power(xj, 6LL) +

            9084900LL*Power(r, 7LL)*Power(xj, 7LL) +

            4082100LL*Power(r, 8LL)*Power(xj, 8LL) + 121935LL*Power(r, 9LL)*Power(xj, 9LL) -

            21494LL*Power(r, 10LL)*Power(xj, 10LL) - 766LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    5LL*Power(xi, 28LL)*Power(xj, 2LL)*

                    (19646550LL + 36018675LL*r*xj + 32744250LL*Power(r, 2LL)*Power(xj, 2LL) +

            19646550LL*Power(r, 3LL)*Power(xj, 3LL) +

            8731800LL*Power(r, 4LL)*Power(xj, 4LL) +

            3056130LL*Power(r, 5LL)*Power(xj, 5LL) + 873180LL*Power(r, 6LL)*Power(xj, 6LL) +

            207900LL*Power(r, 7LL)*Power(xj, 7LL) + 41580LL*Power(r, 8LL)*Power(xj, 8LL) +

            6930LL*Power(r, 9LL)*Power(xj, 9LL) + 924LL*Power(r, 10LL)*Power(xj, 10LL) -

            4LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    26334LL*Power(xi, 16LL)*Power(xj, 14LL)*

                    (43880400LL - 186686775LL*r*xj + 576771750LL*Power(r, 2LL)*Power(xj, 2LL) -

            398603250LL*Power(r, 3LL)*Power(xj, 3LL) +

            72552600LL*Power(r, 4LL)*Power(xj, 4LL) +

            27903120LL*Power(r, 5LL)*Power(xj, 5LL) -

            342720LL*Power(r, 6LL)*Power(xj, 6LL) - 574800LL*Power(r, 7LL)*Power(xj, 7LL) -

            50800LL*Power(r, 8LL)*Power(xj, 8LL) - 945LL*Power(r, 9LL)*Power(xj, 9LL) +

            58LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    10LL*Power(xj, 30LL)*(97302813300LL + 89194245525LL*r*xj +

                                          38780106750LL*Power(r, 2LL)*Power(xj, 2LL) +

                                          10576392750LL*Power(r, 3LL)*Power(xj, 3LL) +

                                          2014551000LL*Power(r, 4LL)*Power(xj, 4LL) +

                                          282037140LL*Power(r, 5LL)*Power(xj, 5LL) +

                                          29688120LL*Power(r, 6LL)*Power(xj, 6LL) +

                                          2356200LL*Power(r, 7LL)*Power(xj, 7LL) + 138600LL*Power(r, 8LL)*Power(xj, 8LL) +

                                          5775LL*Power(r, 9LL)*Power(xj, 9LL) + 154LL*Power(r, 10LL)*Power(xj, 10LL) +

                                          2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    10LL*Power(xi, 2LL)*Power(xj, 28LL)*

                    (4582499159700LL + 3733416276975LL*r*xj +

            1428215931450LL*Power(r, 2LL)*Power(xj, 2LL) +

            338545295550LL*Power(r, 3LL)*Power(xj, 3LL) +

            55198697400LL*Power(r, 4LL)*Power(xj, 4LL) +

            6486854220LL*Power(r, 5LL)*Power(xj, 5LL) +

            558419400LL*Power(r, 6LL)*Power(xj, 6LL) +

            34939080LL*Power(r, 7LL)*Power(xj, 7LL) +

            1532520LL*Power(r, 8LL)*Power(xj, 8LL) + 43285LL*Power(r, 9LL)*Power(xj, 9LL) +

            638LL*Power(r, 10LL)*Power(xj, 10LL) + 2LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    110LL*Power(xi, 10LL)*Power(xj, 20LL)*

                    (-14063418170550LL + 6795156458475LL*r*xj +

            2067471236250LL*Power(r, 2LL)*Power(xj, 2LL) -

            214664924250LL*Power(r, 3LL)*Power(xj, 3LL) -

            124416469800LL*Power(r, 4LL)*Power(xj, 4LL) -

            14935545450LL*Power(r, 5LL)*Power(xj, 5LL) -

            256688460LL*Power(r, 6LL)*Power(xj, 6LL) +

            105750900LL*Power(r, 7LL)*Power(xj, 7LL) +

            11502180LL*Power(r, 8LL)*Power(xj, 8LL) +

            518085LL*Power(r, 9LL)*Power(xj, 9LL) + 9294LL*Power(r, 10LL)*Power(xj, 10LL) +

            2LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    55LL*Power(xi, 20LL)*Power(xj, 10LL)*

                    (1730682450LL + 3172917825LL*r*xj +

            2884470750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1571960250LL*Power(r, 3LL)*Power(xj, 3LL) +

            1404081000LL*Power(r, 4LL)*Power(xj, 4LL) -

            426654270LL*Power(r, 5LL)*Power(xj, 5LL) +

            283536540LL*Power(r, 6LL)*Power(xj, 6LL) +

            39116700LL*Power(r, 7LL)*Power(xj, 7LL) -

            2659860LL*Power(r, 8LL)*Power(xj, 8LL) - 528850LL*Power(r, 9LL)*Power(xj, 9LL) -

            18236LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    5LL*Power(xi, 30LL)*(935550LL + 1715175LL*r*xj +

                                         1559250LL*Power(r, 2LL)*Power(xj, 2LL) + 935550LL*Power(r, 3LL)*Power(xj, 3LL) +

                                         415800LL*Power(r, 4LL)*Power(xj, 4LL) + 145530LL*Power(r, 5LL)*Power(xj, 5LL) +

                                         41580LL*Power(r, 6LL)*Power(xj, 6LL) + 9900LL*Power(r, 7LL)*Power(xj, 7LL) +

                                         1980LL*Power(r, 8LL)*Power(xj, 8LL) + 330LL*Power(r, 9LL)*Power(xj, 9LL) +

                                         44LL*Power(r, 10LL)*Power(xj, 10LL) + 4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    13167LL*Power(xi, 14LL)*Power(xj, 16LL)*

                    (-2319354450LL + 8540029575LL*r*xj -

            7335672750LL*Power(r, 2LL)*Power(xj, 2LL) +

            1133154750LL*Power(r, 3LL)*Power(xj, 3LL) +

            575014200LL*Power(r, 4LL)*Power(xj, 4LL) -

            913710LL*Power(r, 5LL)*Power(xj, 5LL) -

            14863940LL*Power(r, 6LL)*Power(xj, 6LL) -

            1687300LL*Power(r, 7LL)*Power(xj, 7LL) - 46900LL*Power(r, 8LL)*Power(xj, 8LL) +

            3210LL*Power(r, 9LL)*Power(xj, 9LL) + 236LL*Power(r, 10LL)*Power(xj, 10LL) +

            4LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    770LL*Power(xi, 18LL)*Power(xj, 12LL)*

                    (329653800LL + 654729075LL*r*xj + 45785250LL*Power(r, 2LL)*Power(xj, 2LL) +

            1602483750LL*Power(r, 3LL)*Power(xj, 3LL) -

            915705000LL*Power(r, 4LL)*Power(xj, 4LL) +

            266036400LL*Power(r, 5LL)*Power(xj, 5LL) +

            63745920LL*Power(r, 6LL)*Power(xj, 6LL) -

            2304000LL*Power(r, 7LL)*Power(xj, 7LL) -

            1074240LL*Power(r, 8LL)*Power(xj, 8LL) - 64635LL*Power(r, 9LL)*Power(xj, 9LL) -

            514LL*Power(r, 10LL)*Power(xj, 10LL) + 34LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    385LL*Power(xi, 12LL)*Power(xj, 18LL)*

                    (973565393850LL - 1429122323475LL*r*xj +

            298281831750LL*Power(r, 2LL)*Power(xj, 2LL) +

            138841148250LL*Power(r, 3LL)*Power(xj, 3LL) -

            2454240600LL*Power(r, 4LL)*Power(xj, 4LL) -

            4925394810LL*Power(r, 5LL)*Power(xj, 5LL) -

            623832300LL*Power(r, 6LL)*Power(xj, 6LL) -

            19098540LL*Power(r, 7LL)*Power(xj, 7LL) +

            2083140LL*Power(r, 8LL)*Power(xj, 8LL) + 212430LL*Power(r, 9LL)*Power(xj, 9LL) +

            7012LL*Power(r, 10LL)*Power(xj, 10LL) + 68LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    14LL*Power(xi, 26LL)*Power(xj, 4LL)*

                    (-70166250LL - 128638125LL*r*xj - 116943750LL*Power(r, 2LL)*Power(xj, 2LL) -

            70166250LL*Power(r, 3LL)*Power(xj, 3LL) -

            31185000LL*Power(r, 4LL)*Power(xj, 4LL) -

            10914750LL*Power(r, 5LL)*Power(xj, 5LL) -

            3118500LL*Power(r, 6LL)*Power(xj, 6LL) - 742500LL*Power(r, 7LL)*Power(xj, 7LL) -

            148500LL*Power(r, 8LL)*Power(xj, 8LL) - 32615LL*Power(r, 9LL)*Power(xj, 9LL) -

            154LL*Power(r, 10LL)*Power(xj, 10LL) + 74LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    7LL*Power(xi, 4LL)*Power(xj, 26LL)*

                    (-69822945249750LL - 46669577290875LL*r*xj -

            14025037430250LL*Power(r, 2LL)*Power(xj, 2LL) -

            2430881664750LL*Power(r, 3LL)*Power(xj, 3LL) -

            251629270200LL*Power(r, 4LL)*Power(xj, 4LL) -

            12434519790LL*Power(r, 5LL)*Power(xj, 5LL) +

            452930940LL*Power(r, 6LL)*Power(xj, 6LL) +

            131125500LL*Power(r, 7LL)*Power(xj, 7LL) +

            11018700LL*Power(r, 8LL)*Power(xj, 8LL) +

            514470LL*Power(r, 9LL)*Power(xj, 9LL) + 13332LL*Power(r, 10LL)*Power(xj, 10LL) +

            148LL*Power(r, 11LL)*Power(xj, 11LL)) -

                    50LL*Power(xi, 8LL)*Power(xj, 22LL)*

                    (-51768833574150LL - 5003280391725LL*r*xj +

            4493439477450LL*Power(r, 2LL)*Power(xj, 2LL) +

            1286866176750LL*Power(r, 3LL)*Power(xj, 3LL) +

            111437476920LL*Power(r, 4LL)*Power(xj, 4LL) -

            6620313546LL*Power(r, 5LL)*Power(xj, 5LL) -

            2406603276LL*Power(r, 6LL)*Power(xj, 6LL) -

            242686620LL*Power(r, 7LL)*Power(xj, 7LL) -

            12228876LL*Power(r, 8LL)*Power(xj, 8LL) -

            256223LL*Power(r, 9LL)*Power(xj, 9LL) + 2486LL*Power(r, 10LL)*Power(xj, 10LL) +

            158LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    25LL*Power(xi, 22LL)*Power(xj, 8LL)*

                    (-1119853350LL - 2053064475LL*r*xj -

            1866422250LL*Power(r, 2LL)*Power(xj, 2LL) -

            1119853350LL*Power(r, 3LL)*Power(xj, 3LL) -

            497712600LL*Power(r, 4LL)*Power(xj, 4LL) -

            194415606LL*Power(r, 5LL)*Power(xj, 5LL) -

            9338868LL*Power(r, 6LL)*Power(xj, 6LL) -

            31217076LL*Power(r, 7LL)*Power(xj, 7LL) -

            2256804LL*Power(r, 8LL)*Power(xj, 8LL) + 246774LL*Power(r, 9LL)*Power(xj, 9LL) +

            22836LL*Power(r, 10LL)*Power(xj, 10LL) + 316LL*Power(r, 11LL)*Power(xj, 11LL)) +

                    3LL*Power(xi, 6LL)*Power(xj, 24LL)*

                    (596006592662250LL + 266778699697125LL*r*xj +

            37515651153750LL*Power(r, 2LL)*Power(xj, 2LL) -

            2214626163750LL*Power(r, 3LL)*Power(xj, 3LL) -

            1538075107800LL*Power(r, 4LL)*Power(xj, 4LL) -

            248955308910LL*Power(r, 5LL)*Power(xj, 5LL) -

            21434337540LL*Power(r, 6LL)*Power(xj, 6LL) -

            957980100LL*Power(r, 7LL)*Power(xj, 7LL) -

            4874100LL*Power(r, 8LL)*Power(xj, 8LL) + 1831830LL*Power(r, 9LL)*Power(xj, 9LL) +

            91828LL*Power(r, 10LL)*Power(xj, 10LL) + 1532LL*Power(r, 11LL)*Power(xj, 11LL))))/

                (4.67775e6*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 21LL)*Power(xi + xj, 21LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_6S_5S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_5S_6S(r, xj, xi);
}
