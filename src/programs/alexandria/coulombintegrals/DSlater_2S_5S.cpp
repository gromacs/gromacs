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

cl_R DSlater_2S_5S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-224622748350LL*xi + 249080832000LL*exp(2LL*r*xi)*xi -

                  400329329400LL*r*Power(xi, 2LL) - 351747621225LL*Power(r, 2LL)*Power(xi, 3LL) -

                  202556554200LL*Power(r, 3LL)*Power(xi, 4LL) -

                  85662076500LL*Power(r, 4LL)*Power(xi, 5LL) -

                  28229160960LL*Power(r, 5LL)*Power(xi, 6LL) -

                  7498370880LL*Power(r, 6LL)*Power(xi, 7LL) -

                  1635828480LL*Power(r, 7LL)*Power(xi, 8LL) -

                  295289280LL*Power(r, 8LL)*Power(xi, 9LL) - 43929600LL*Power(r, 9LL)*Power(xi, 10LL) -

                  5271552LL*Power(r, 10LL)*Power(xi, 11LL) - 479232LL*Power(r, 11LL)*Power(xi, 12LL) -

                  26624LL*Power(r, 12LL)*Power(xi, 13LL))/(1.24540416e11*exp(2LL*r*xi)*r) +

                (-124540416000LL + 124540416000LL*exp(2LL*r*xi) - 224622748350LL*r*xi -

                 200164664700LL*Power(r, 2LL)*Power(xi, 2LL) -

                 117249207075LL*Power(r, 3LL)*Power(xi, 3LL) -

                 50639138550LL*Power(r, 4LL)*Power(xi, 4LL) -

                 17132415300LL*Power(r, 5LL)*Power(xi, 5LL) -

                 4704860160LL*Power(r, 6LL)*Power(xi, 6LL) -

                 1071195840LL*Power(r, 7LL)*Power(xi, 7LL) - 204478560LL*Power(r, 8LL)*Power(xi, 8LL) -

                 32809920LL*Power(r, 9LL)*Power(xi, 9LL) - 4392960LL*Power(r, 10LL)*Power(xi, 10LL) -

                 479232LL*Power(r, 11LL)*Power(xi, 11LL) - 39936LL*Power(r, 12LL)*Power(xi, 12LL) -

                 2048LL*Power(r, 13LL)*Power(xi, 13LL))/

                (1.24540416e11*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-124540416000LL + 124540416000LL*exp(2LL*r*xi) - 224622748350LL*r*xi -

                     200164664700LL*Power(r, 2LL)*Power(xi, 2LL) -

                     117249207075LL*Power(r, 3LL)*Power(xi, 3LL) -

                     50639138550LL*Power(r, 4LL)*Power(xi, 4LL) -

                     17132415300LL*Power(r, 5LL)*Power(xi, 5LL) -

                     4704860160LL*Power(r, 6LL)*Power(xi, 6LL) -

                     1071195840LL*Power(r, 7LL)*Power(xi, 7LL) -

                     204478560LL*Power(r, 8LL)*Power(xi, 8LL) - 32809920LL*Power(r, 9LL)*Power(xi, 9LL) -

                     4392960LL*Power(r, 10LL)*Power(xi, 10LL) - 479232LL*Power(r, 11LL)*Power(xi, 11LL) -

                     39936LL*Power(r, 12LL)*Power(xi, 12LL) - 2048LL*Power(r, 13LL)*Power(xi, 13LL)))/

                (6.2270208e10*exp(2LL*r*xi)*r)

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
            S = (28350LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                 945LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-210LL*Power(r, 2LL)*Power(xi, 16LL) - 10LL*Power(r, 3LL)*Power(xi, 17LL) +

                  30LL*Power(xj, 14LL) + 45LL*r*xi*Power(xj, 14LL) +

                  39LL*r*Power(xi, 7LL)*Power(xj, 8LL)*(1309LL - 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  858LL*Power(xi, 8LL)*Power(xj, 6LL)*(-305LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  30LL*Power(xi, 2LL)*Power(xj, 12LL)*(-13LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  390LL*Power(xi, 4LL)*Power(xj, 10LL)*(-6LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  143LL*r*Power(xi, 9LL)*Power(xj, 6LL)*(-153LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  5LL*r*Power(xi, 3LL)*Power(xj, 12LL)*(-117LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  45LL*r*Power(xi, 15LL)*(35LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  138LL*Power(xi, 12LL)*Power(xj, 2LL)*(580LL + 13LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  150LL*Power(xi, 14LL)*(28LL + 17LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  13LL*r*Power(xi, 11LL)*Power(xj, 4LL)*

                  (-4071LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  3LL*r*Power(xi, 13LL)*Power(xj, 2LL)*(-8135LL + 26LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*(2171LL + 30LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  234LL*Power(xi, 10LL)*Power(xj, 4LL)*(-1235LL + 33LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  78LL*Power(xi, 6LL)*Power(xj, 8LL)*(550LL + 47LL*Power(r, 2LL)*Power(xj, 2LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-819LL*Power(xi, 10LL)*Power(xj, 10LL)*

                  (22275LL + 39780LL*r*xj + 38160LL*Power(r, 2LL)*Power(xj, 2LL) +

            16560LL*Power(r, 3LL)*Power(xj, 3LL) + 9840LL*Power(r, 4LL)*Power(xj, 4LL) +

            3900LL*Power(r, 5LL)*Power(xj, 5LL) + 816LL*Power(r, 6LL)*Power(xj, 6LL) +

            88LL*Power(r, 7LL)*Power(xj, 7LL) + 4LL*Power(r, 8LL)*Power(xj, 8LL)) +

                  Power(xi, 20LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                                   1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                                   108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                                   2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  Power(xj, 20LL)*(16216200LL + 18243225LL*r*xj +

                                   9729720LL*Power(r, 2LL)*Power(xj, 2LL) +

                                   3243240LL*Power(r, 3LL)*Power(xj, 3LL) +

                                   748440LL*Power(r, 4LL)*Power(xj, 4LL) + 124740LL*Power(r, 5LL)*Power(xj, 5LL) +

                                   15120LL*Power(r, 6LL)*Power(xj, 6LL) + 1296LL*Power(r, 7LL)*Power(xj, 7LL) +

                                   72LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  18LL*Power(xi, 16LL)*Power(xj, 4LL)*

                  (61425LL + 110565LL*r*xj + 98280LL*Power(r, 2LL)*Power(xj, 2LL) +

            57330LL*Power(r, 3LL)*Power(xj, 3LL) + 24570LL*Power(r, 4LL)*Power(xj, 4LL) +

            8190LL*Power(r, 5LL)*Power(xj, 5LL) + 2184LL*Power(r, 6LL)*Power(xj, 6LL) +

            496LL*Power(r, 7LL)*Power(xj, 7LL) + 64LL*Power(r, 8LL)*Power(xj, 8LL) +

            3LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  18LL*Power(xi, 4LL)*Power(xj, 16LL)*

                  (6572475LL - 3161340LL*r*xj - 4782960LL*Power(r, 2LL)*Power(xj, 2LL) -

            1912365LL*Power(r, 3LL)*Power(xj, 3LL) -

            378105LL*Power(r, 4LL)*Power(xj, 4LL) - 34125LL*Power(r, 5LL)*Power(xj, 5LL) +

            1092LL*Power(r, 6LL)*Power(xj, 6LL) + 650LL*Power(r, 7LL)*Power(xj, 7LL) +

            71LL*Power(r, 8LL)*Power(xj, 8LL) + 3LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  21LL*Power(xi, 8LL)*Power(xj, 12LL)*

                  (-1063800LL - 2775735LL*r*xj - 862920LL*Power(r, 2LL)*Power(xj, 2LL) -

            1132020LL*Power(r, 3LL)*Power(xj, 3LL) -

            698580LL*Power(r, 4LL)*Power(xj, 4LL) - 196920LL*Power(r, 5LL)*Power(xj, 5LL) -

            28992LL*Power(r, 6LL)*Power(xj, 6LL) - 2064LL*Power(r, 7LL)*Power(xj, 7LL) -

            24LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  21LL*Power(xi, 12LL)*Power(xj, 8LL)*

                  (482625LL + 868725LL*r*xj + 772200LL*Power(r, 2LL)*Power(xj, 2LL) +

            455400LL*Power(r, 3LL)*Power(xj, 3LL) + 178200LL*Power(r, 4LL)*Power(xj, 4LL) +

            72180LL*Power(r, 5LL)*Power(xj, 5LL) + 19920LL*Power(r, 6LL)*Power(xj, 6LL) +

            2952LL*Power(r, 7LL)*Power(xj, 7LL) + 204LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  6LL*Power(xi, 6LL)*Power(xj, 14LL)*

                  (-10357200LL + 5071815LL*r*xj - 6463800LL*Power(r, 2LL)*Power(xj, 2LL) -

            7151130LL*Power(r, 3LL)*Power(xj, 3LL) -

            2572290LL*Power(r, 4LL)*Power(xj, 4LL) -

            468720LL*Power(r, 5LL)*Power(xj, 5LL) - 42672LL*Power(r, 6LL)*Power(xj, 6LL) -

            648LL*Power(r, 7LL)*Power(xj, 7LL) + 228LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  Power(xi, 18LL)*Power(xj, 2LL)*

                  (184275LL + 331695LL*r*xj + 294840LL*Power(r, 2LL)*Power(xj, 2LL) +

            171990LL*Power(r, 3LL)*Power(xj, 3LL) + 73710LL*Power(r, 4LL)*Power(xj, 4LL) +

            24570LL*Power(r, 5LL)*Power(xj, 5LL) + 6552LL*Power(r, 6LL)*Power(xj, 6LL) +

            1404LL*Power(r, 7LL)*Power(xj, 7LL) + 234LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL)) +

                  Power(xi, 2LL)*Power(xj, 18LL)*

                  (-133783650LL - 107432325LL*r*xj - 35675640LL*Power(r, 2LL)*Power(xj, 2LL) -

            5135130LL*Power(r, 3LL)*Power(xj, 3LL) +

            270270LL*Power(r, 4LL)*Power(xj, 4LL) + 270270LL*Power(r, 5LL)*Power(xj, 5LL) +

            57960LL*Power(r, 6LL)*Power(xj, 6LL) + 6948LL*Power(r, 7LL)*Power(xj, 7LL) +

            486LL*Power(r, 8LL)*Power(xj, 8LL) + 16LL*Power(r, 9LL)*Power(xj, 9LL)) -

                  6LL*Power(xi, 14LL)*Power(xj, 6LL)*

                  (675675LL + 1216215LL*r*xj + 1081080LL*Power(r, 2LL)*Power(xj, 2LL) +

            630630LL*Power(r, 3LL)*Power(xj, 3LL) + 270270LL*Power(r, 4LL)*Power(xj, 4LL) +

            88200LL*Power(r, 5LL)*Power(xj, 5LL) + 26544LL*Power(r, 6LL)*Power(xj, 6LL) +

            5160LL*Power(r, 7LL)*Power(xj, 7LL) + 492LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (28350LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 13LL)*

                 Power(xi + xj, 13LL)) + (28350LL*exp(2LL*r*(xi + xj))*

                                          Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                                          945LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                                          (-210LL*Power(r, 2LL)*Power(xi, 16LL) - 10LL*Power(r, 3LL)*Power(xi, 17LL) +

                                  30LL*Power(xj, 14LL) + 45LL*r*xi*Power(xj, 14LL) +

                                  39LL*r*Power(xi, 7LL)*Power(xj, 8LL)*(1309LL - 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  858LL*Power(xi, 8LL)*Power(xj, 6LL)*(-305LL + Power(r, 2LL)*Power(xj, 2LL)) +

                                  30LL*Power(xi, 2LL)*Power(xj, 12LL)*(-13LL + Power(r, 2LL)*Power(xj, 2LL)) -

                                  390LL*Power(xi, 4LL)*Power(xj, 10LL)*(-6LL + Power(r, 2LL)*Power(xj, 2LL)) -

                                  143LL*r*Power(xi, 9LL)*Power(xj, 6LL)*(-153LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  5LL*r*Power(xi, 3LL)*Power(xj, 12LL)*(-117LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  45LL*r*Power(xi, 15LL)*(35LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  138LL*Power(xi, 12LL)*Power(xj, 2LL)*(580LL + 13LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  150LL*Power(xi, 14LL)*(28LL + 17LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  13LL*r*Power(xi, 11LL)*Power(xj, 4LL)*

                                  (-4071LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  3LL*r*Power(xi, 13LL)*Power(xj, 2LL)*(-8135LL + 26LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*(2171LL + 30LL*Power(r, 2LL)*Power(xj, 2LL)) +

                                  234LL*Power(xi, 10LL)*Power(xj, 4LL)*(-1235LL + 33LL*Power(r, 2LL)*Power(xj, 2LL)) -

                                  78LL*Power(xi, 6LL)*Power(xj, 8LL)*(550LL + 47LL*Power(r, 2LL)*Power(xj, 2LL))) -

                                          2LL*exp(2LL*r*xi)*Power(xi, 6LL)*

                                          (-819LL*Power(xi, 10LL)*Power(xj, 10LL)*

                                  (22275LL + 39780LL*r*xj + 38160LL*Power(r, 2LL)*Power(xj, 2LL) +

            16560LL*Power(r, 3LL)*Power(xj, 3LL) + 9840LL*Power(r, 4LL)*Power(xj, 4LL) +

            3900LL*Power(r, 5LL)*Power(xj, 5LL) + 816LL*Power(r, 6LL)*Power(xj, 6LL) +

            88LL*Power(r, 7LL)*Power(xj, 7LL) + 4LL*Power(r, 8LL)*Power(xj, 8LL)) +

                                  Power(xi, 20LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                                                   1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                                                   108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                                                   2LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  Power(xj, 20LL)*(16216200LL + 18243225LL*r*xj +

                                                   9729720LL*Power(r, 2LL)*Power(xj, 2LL) +

                                                   3243240LL*Power(r, 3LL)*Power(xj, 3LL) +

                                                   748440LL*Power(r, 4LL)*Power(xj, 4LL) + 124740LL*Power(r, 5LL)*Power(xj, 5LL) +

                                                   15120LL*Power(r, 6LL)*Power(xj, 6LL) + 1296LL*Power(r, 7LL)*Power(xj, 7LL) +

                                                   72LL*Power(r, 8LL)*Power(xj, 8LL) + 2LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  18LL*Power(xi, 16LL)*Power(xj, 4LL)*

                                  (61425LL + 110565LL*r*xj + 98280LL*Power(r, 2LL)*Power(xj, 2LL) +

            57330LL*Power(r, 3LL)*Power(xj, 3LL) + 24570LL*Power(r, 4LL)*Power(xj, 4LL) +

            8190LL*Power(r, 5LL)*Power(xj, 5LL) + 2184LL*Power(r, 6LL)*Power(xj, 6LL) +

            496LL*Power(r, 7LL)*Power(xj, 7LL) + 64LL*Power(r, 8LL)*Power(xj, 8LL) +

            3LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  18LL*Power(xi, 4LL)*Power(xj, 16LL)*

                                  (6572475LL - 3161340LL*r*xj - 4782960LL*Power(r, 2LL)*Power(xj, 2LL) -

            1912365LL*Power(r, 3LL)*Power(xj, 3LL) -

            378105LL*Power(r, 4LL)*Power(xj, 4LL) - 34125LL*Power(r, 5LL)*Power(xj, 5LL) +

            1092LL*Power(r, 6LL)*Power(xj, 6LL) + 650LL*Power(r, 7LL)*Power(xj, 7LL) +

            71LL*Power(r, 8LL)*Power(xj, 8LL) + 3LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  21LL*Power(xi, 8LL)*Power(xj, 12LL)*

                                  (-1063800LL - 2775735LL*r*xj - 862920LL*Power(r, 2LL)*Power(xj, 2LL) -

            1132020LL*Power(r, 3LL)*Power(xj, 3LL) -

            698580LL*Power(r, 4LL)*Power(xj, 4LL) - 196920LL*Power(r, 5LL)*Power(xj, 5LL) -

            28992LL*Power(r, 6LL)*Power(xj, 6LL) - 2064LL*Power(r, 7LL)*Power(xj, 7LL) -

            24LL*Power(r, 8LL)*Power(xj, 8LL) + 4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  21LL*Power(xi, 12LL)*Power(xj, 8LL)*

                                  (482625LL + 868725LL*r*xj + 772200LL*Power(r, 2LL)*Power(xj, 2LL) +

            455400LL*Power(r, 3LL)*Power(xj, 3LL) + 178200LL*Power(r, 4LL)*Power(xj, 4LL) +

            72180LL*Power(r, 5LL)*Power(xj, 5LL) + 19920LL*Power(r, 6LL)*Power(xj, 6LL) +

            2952LL*Power(r, 7LL)*Power(xj, 7LL) + 204LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  6LL*Power(xi, 6LL)*Power(xj, 14LL)*

                                  (-10357200LL + 5071815LL*r*xj - 6463800LL*Power(r, 2LL)*Power(xj, 2LL) -

            7151130LL*Power(r, 3LL)*Power(xj, 3LL) -

            2572290LL*Power(r, 4LL)*Power(xj, 4LL) -

            468720LL*Power(r, 5LL)*Power(xj, 5LL) - 42672LL*Power(r, 6LL)*Power(xj, 6LL) -

            648LL*Power(r, 7LL)*Power(xj, 7LL) + 228LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  Power(xi, 18LL)*Power(xj, 2LL)*

                                  (184275LL + 331695LL*r*xj + 294840LL*Power(r, 2LL)*Power(xj, 2LL) +

            171990LL*Power(r, 3LL)*Power(xj, 3LL) + 73710LL*Power(r, 4LL)*Power(xj, 4LL) +

            24570LL*Power(r, 5LL)*Power(xj, 5LL) + 6552LL*Power(r, 6LL)*Power(xj, 6LL) +

            1404LL*Power(r, 7LL)*Power(xj, 7LL) + 234LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL)) +

                                  Power(xi, 2LL)*Power(xj, 18LL)*

                                  (-133783650LL - 107432325LL*r*xj - 35675640LL*Power(r, 2LL)*Power(xj, 2LL) -

            5135130LL*Power(r, 3LL)*Power(xj, 3LL) +

            270270LL*Power(r, 4LL)*Power(xj, 4LL) + 270270LL*Power(r, 5LL)*Power(xj, 5LL) +

            57960LL*Power(r, 6LL)*Power(xj, 6LL) + 6948LL*Power(r, 7LL)*Power(xj, 7LL) +

            486LL*Power(r, 8LL)*Power(xj, 8LL) + 16LL*Power(r, 9LL)*Power(xj, 9LL)) -

                                  6LL*Power(xi, 14LL)*Power(xj, 6LL)*

                                  (675675LL + 1216215LL*r*xj + 1081080LL*Power(r, 2LL)*Power(xj, 2LL) +

            630630LL*Power(r, 3LL)*Power(xj, 3LL) + 270270LL*Power(r, 4LL)*Power(xj, 4LL) +

            88200LL*Power(r, 5LL)*Power(xj, 5LL) + 26544LL*Power(r, 6LL)*Power(xj, 6LL) +

            5160LL*Power(r, 7LL)*Power(xj, 7LL) + 492LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (14175LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 13LL)*Power(xi + xj, 12LL)) -

                (56700LL*exp(2LL*r*(xi + xj))*(xi + xj)*

                 Power(Power(xi, 2LL) - Power(xj, 2LL), 13LL) +

                 945LL*exp(2LL*r*xj)*Power(xj, 12LL)*

                 (-420LL*r*Power(xi, 16LL) - 30LL*Power(r, 2LL)*Power(xi, 17LL) -

        5100LL*r*Power(xi, 14LL)*Power(xj, 2LL) -

        180LL*Power(r, 2LL)*Power(xi, 15LL)*Power(xj, 2LL) -

        3588LL*r*Power(xi, 12LL)*Power(xj, 4LL) +

        156LL*Power(r, 2LL)*Power(xi, 13LL)*Power(xj, 4LL) +

        15444LL*r*Power(xi, 10LL)*Power(xj, 6LL) +

        572LL*Power(r, 2LL)*Power(xi, 11LL)*Power(xj, 6LL) +

        1716LL*r*Power(xi, 8LL)*Power(xj, 8LL) -

        572LL*Power(r, 2LL)*Power(xi, 9LL)*Power(xj, 8LL) -

        7332LL*r*Power(xi, 6LL)*Power(xj, 10LL) -

        156LL*Power(r, 2LL)*Power(xi, 7LL)*Power(xj, 10LL) -

        780LL*r*Power(xi, 4LL)*Power(xj, 12LL) +

        180LL*Power(r, 2LL)*Power(xi, 5LL)*Power(xj, 12LL) + 45LL*xi*Power(xj, 14LL) +

        60LL*r*Power(xi, 2LL)*Power(xj, 14LL) +

        20LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 14LL) +

        39LL*Power(xi, 7LL)*Power(xj, 8LL)*(1309LL - 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        143LL*Power(xi, 9LL)*Power(xj, 6LL)*(-153LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        5LL*Power(xi, 3LL)*Power(xj, 12LL)*(-117LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        45LL*Power(xi, 15LL)*(35LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        13LL*Power(xi, 11LL)*Power(xj, 4LL)*(-4071LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*Power(xi, 13LL)*Power(xj, 2LL)*(-8135LL + 26LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*Power(xi, 5LL)*Power(xj, 10LL)*(2171LL + 30LL*Power(r, 2LL)*Power(xj, 2LL))) +

                 1890LL*exp(2LL*r*xj)*Power(xj, 13LL)*

                 (-210LL*Power(r, 2LL)*Power(xi, 16LL) - 10LL*Power(r, 3LL)*Power(xi, 17LL) +

        30LL*Power(xj, 14LL) + 45LL*r*xi*Power(xj, 14LL) +

        39LL*r*Power(xi, 7LL)*Power(xj, 8LL)*(1309LL - 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        858LL*Power(xi, 8LL)*Power(xj, 6LL)*(-305LL + Power(r, 2LL)*Power(xj, 2LL)) +

        30LL*Power(xi, 2LL)*Power(xj, 12LL)*(-13LL + Power(r, 2LL)*Power(xj, 2LL)) -

        390LL*Power(xi, 4LL)*Power(xj, 10LL)*(-6LL + Power(r, 2LL)*Power(xj, 2LL)) -

        143LL*r*Power(xi, 9LL)*Power(xj, 6LL)*(-153LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        5LL*r*Power(xi, 3LL)*Power(xj, 12LL)*(-117LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        45LL*r*Power(xi, 15LL)*(35LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        138LL*Power(xi, 12LL)*Power(xj, 2LL)*(580LL + 13LL*Power(r, 2LL)*Power(xj, 2LL)) -

        150LL*Power(xi, 14LL)*(28LL + 17LL*Power(r, 2LL)*Power(xj, 2LL)) +

        13LL*r*Power(xi, 11LL)*Power(xj, 4LL)*(-4071LL + 22LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*r*Power(xi, 13LL)*Power(xj, 2LL)*(-8135LL + 26LL*Power(r, 2LL)*Power(xj, 2LL)) +

        3LL*r*Power(xi, 5LL)*Power(xj, 10LL)*(2171LL + 30LL*Power(r, 2LL)*Power(xj, 2LL)) +

        234LL*Power(xi, 10LL)*Power(xj, 4LL)*(-1235LL + 33LL*Power(r, 2LL)*Power(xj, 2LL)) -

        78LL*Power(xi, 6LL)*Power(xj, 8LL)*(550LL + 47LL*Power(r, 2LL)*Power(xj, 2LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 6LL)*

                 (-819LL*Power(xi, 10LL)*Power(xj, 10LL)*

        (39780LL*xj + 76320LL*r*Power(xj, 2LL) + 49680LL*Power(r, 2LL)*Power(xj, 3LL) +

            39360LL*Power(r, 3LL)*Power(xj, 4LL) + 19500LL*Power(r, 4LL)*Power(xj, 5LL) +

            4896LL*Power(r, 5LL)*Power(xj, 6LL) + 616LL*Power(r, 6LL)*Power(xj, 7LL) +

            32LL*Power(r, 7LL)*Power(xj, 8LL)) +

        Power(xi, 20LL)*(25515LL*xj + 45360LL*r*Power(xj, 2LL) +

                         39690LL*Power(r, 2LL)*Power(xj, 3LL) + 22680LL*Power(r, 3LL)*Power(xj, 4LL) +

                         9450LL*Power(r, 4LL)*Power(xj, 5LL) + 3024LL*Power(r, 5LL)*Power(xj, 6LL) +

                         756LL*Power(r, 6LL)*Power(xj, 7LL) + 144LL*Power(r, 7LL)*Power(xj, 8LL) +

                         18LL*Power(r, 8LL)*Power(xj, 9LL)) -

        Power(xj, 20LL)*(18243225LL*xj + 19459440LL*r*Power(xj, 2LL) +

                         9729720LL*Power(r, 2LL)*Power(xj, 3LL) +

                         2993760LL*Power(r, 3LL)*Power(xj, 4LL) +

                         623700LL*Power(r, 4LL)*Power(xj, 5LL) + 90720LL*Power(r, 5LL)*Power(xj, 6LL) +

                         9072LL*Power(r, 6LL)*Power(xj, 7LL) + 576LL*Power(r, 7LL)*Power(xj, 8LL) +

                         18LL*Power(r, 8LL)*Power(xj, 9LL)) +

        18LL*Power(xi, 16LL)*Power(xj, 4LL)*

        (110565LL*xj + 196560LL*r*Power(xj, 2LL) +

            171990LL*Power(r, 2LL)*Power(xj, 3LL) + 98280LL*Power(r, 3LL)*Power(xj, 4LL) +

            40950LL*Power(r, 4LL)*Power(xj, 5LL) + 13104LL*Power(r, 5LL)*Power(xj, 6LL) +

            3472LL*Power(r, 6LL)*Power(xj, 7LL) + 512LL*Power(r, 7LL)*Power(xj, 8LL) +

            27LL*Power(r, 8LL)*Power(xj, 9LL)) -

        18LL*Power(xi, 4LL)*Power(xj, 16LL)*

        (-3161340LL*xj - 9565920LL*r*Power(xj, 2LL) -

            5737095LL*Power(r, 2LL)*Power(xj, 3LL) -

            1512420LL*Power(r, 3LL)*Power(xj, 4LL) -

            170625LL*Power(r, 4LL)*Power(xj, 5LL) + 6552LL*Power(r, 5LL)*Power(xj, 6LL) +

            4550LL*Power(r, 6LL)*Power(xj, 7LL) + 568LL*Power(r, 7LL)*Power(xj, 8LL) +

            27LL*Power(r, 8LL)*Power(xj, 9LL)) -

        21LL*Power(xi, 8LL)*Power(xj, 12LL)*

        (-2775735LL*xj - 1725840LL*r*Power(xj, 2LL) -

            3396060LL*Power(r, 2LL)*Power(xj, 3LL) -

            2794320LL*Power(r, 3LL)*Power(xj, 4LL) -

            984600LL*Power(r, 4LL)*Power(xj, 5LL) - 173952LL*Power(r, 5LL)*Power(xj, 6LL) -

            14448LL*Power(r, 6LL)*Power(xj, 7LL) - 192LL*Power(r, 7LL)*Power(xj, 8LL) +

            36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        21LL*Power(xi, 12LL)*Power(xj, 8LL)*

        (868725LL*xj + 1544400LL*r*Power(xj, 2LL) +

            1366200LL*Power(r, 2LL)*Power(xj, 3LL) +

            712800LL*Power(r, 3LL)*Power(xj, 4LL) + 360900LL*Power(r, 4LL)*Power(xj, 5LL) +

            119520LL*Power(r, 5LL)*Power(xj, 6LL) + 20664LL*Power(r, 6LL)*Power(xj, 7LL) +

            1632LL*Power(r, 7LL)*Power(xj, 8LL) + 36LL*Power(r, 8LL)*Power(xj, 9LL)) +

        6LL*Power(xi, 6LL)*Power(xj, 14LL)*

        (5071815LL*xj - 12927600LL*r*Power(xj, 2LL) -

            21453390LL*Power(r, 2LL)*Power(xj, 3LL) -

            10289160LL*Power(r, 3LL)*Power(xj, 4LL) -

            2343600LL*Power(r, 4LL)*Power(xj, 5LL) -

            256032LL*Power(r, 5LL)*Power(xj, 6LL) - 4536LL*Power(r, 6LL)*Power(xj, 7LL) +

            1824LL*Power(r, 7LL)*Power(xj, 8LL) + 144LL*Power(r, 8LL)*Power(xj, 9LL)) -

        Power(xi, 18LL)*Power(xj, 2LL)*

        (331695LL*xj + 589680LL*r*Power(xj, 2LL) +

            515970LL*Power(r, 2LL)*Power(xj, 3LL) + 294840LL*Power(r, 3LL)*Power(xj, 4LL) +

            122850LL*Power(r, 4LL)*Power(xj, 5LL) + 39312LL*Power(r, 5LL)*Power(xj, 6LL) +

            9828LL*Power(r, 6LL)*Power(xj, 7LL) + 1872LL*Power(r, 7LL)*Power(xj, 8LL) +

            144LL*Power(r, 8LL)*Power(xj, 9LL)) +

        Power(xi, 2LL)*Power(xj, 18LL)*

        (-107432325LL*xj - 71351280LL*r*Power(xj, 2LL) -

            15405390LL*Power(r, 2LL)*Power(xj, 3LL) +

            1081080LL*Power(r, 3LL)*Power(xj, 4LL) +

            1351350LL*Power(r, 4LL)*Power(xj, 5LL) +

            347760LL*Power(r, 5LL)*Power(xj, 6LL) + 48636LL*Power(r, 6LL)*Power(xj, 7LL) +

            3888LL*Power(r, 7LL)*Power(xj, 8LL) + 144LL*Power(r, 8LL)*Power(xj, 9LL)) -

        6LL*Power(xi, 14LL)*Power(xj, 6LL)*

        (1216215LL*xj + 2162160LL*r*Power(xj, 2LL) +

            1891890LL*Power(r, 2LL)*Power(xj, 3LL) +

            1081080LL*Power(r, 3LL)*Power(xj, 4LL) + 441000LL*Power(r, 4LL)*Power(xj, 5LL) +

            159264LL*Power(r, 5LL)*Power(xj, 6LL) + 36120LL*Power(r, 6LL)*Power(xj, 7LL) +

            3936LL*Power(r, 7LL)*Power(xj, 8LL) + 144LL*Power(r, 8LL)*Power(xj, 9LL))) -

                 4LL*exp(2LL*r*xi)*Power(xi, 7LL)*

                 (-819LL*Power(xi, 10LL)*Power(xj, 10LL)*

        (22275LL + 39780LL*r*xj + 38160LL*Power(r, 2LL)*Power(xj, 2LL) +

            16560LL*Power(r, 3LL)*Power(xj, 3LL) + 9840LL*Power(r, 4LL)*Power(xj, 4LL) +

            3900LL*Power(r, 5LL)*Power(xj, 5LL) + 816LL*Power(r, 6LL)*Power(xj, 6LL) +

            88LL*Power(r, 7LL)*Power(xj, 7LL) + 4LL*Power(r, 8LL)*Power(xj, 8LL)) +

        Power(xi, 20LL)*(14175LL + 25515LL*r*xj + 22680LL*Power(r, 2LL)*Power(xj, 2LL) +

                         13230LL*Power(r, 3LL)*Power(xj, 3LL) + 5670LL*Power(r, 4LL)*Power(xj, 4LL) +

                         1890LL*Power(r, 5LL)*Power(xj, 5LL) + 504LL*Power(r, 6LL)*Power(xj, 6LL) +

                         108LL*Power(r, 7LL)*Power(xj, 7LL) + 18LL*Power(r, 8LL)*Power(xj, 8LL) +

                         2LL*Power(r, 9LL)*Power(xj, 9LL)) -

        Power(xj, 20LL)*(16216200LL + 18243225LL*r*xj +

                         9729720LL*Power(r, 2LL)*Power(xj, 2LL) +

                         3243240LL*Power(r, 3LL)*Power(xj, 3LL) + 748440LL*Power(r, 4LL)*Power(xj, 4LL) +

                         124740LL*Power(r, 5LL)*Power(xj, 5LL) + 15120LL*Power(r, 6LL)*Power(xj, 6LL) +

                         1296LL*Power(r, 7LL)*Power(xj, 7LL) + 72LL*Power(r, 8LL)*Power(xj, 8LL) +

                         2LL*Power(r, 9LL)*Power(xj, 9LL)) +

        18LL*Power(xi, 16LL)*Power(xj, 4LL)*

        (61425LL + 110565LL*r*xj + 98280LL*Power(r, 2LL)*Power(xj, 2LL) +

            57330LL*Power(r, 3LL)*Power(xj, 3LL) + 24570LL*Power(r, 4LL)*Power(xj, 4LL) +

            8190LL*Power(r, 5LL)*Power(xj, 5LL) + 2184LL*Power(r, 6LL)*Power(xj, 6LL) +

            496LL*Power(r, 7LL)*Power(xj, 7LL) + 64LL*Power(r, 8LL)*Power(xj, 8LL) +

            3LL*Power(r, 9LL)*Power(xj, 9LL)) -

        18LL*Power(xi, 4LL)*Power(xj, 16LL)*

        (6572475LL - 3161340LL*r*xj - 4782960LL*Power(r, 2LL)*Power(xj, 2LL) -

            1912365LL*Power(r, 3LL)*Power(xj, 3LL) - 378105LL*Power(r, 4LL)*Power(xj, 4LL) -

            34125LL*Power(r, 5LL)*Power(xj, 5LL) + 1092LL*Power(r, 6LL)*Power(xj, 6LL) +

            650LL*Power(r, 7LL)*Power(xj, 7LL) + 71LL*Power(r, 8LL)*Power(xj, 8LL) +

            3LL*Power(r, 9LL)*Power(xj, 9LL)) -

        21LL*Power(xi, 8LL)*Power(xj, 12LL)*

        (-1063800LL - 2775735LL*r*xj - 862920LL*Power(r, 2LL)*Power(xj, 2LL) -

            1132020LL*Power(r, 3LL)*Power(xj, 3LL) - 698580LL*Power(r, 4LL)*Power(xj, 4LL) -

            196920LL*Power(r, 5LL)*Power(xj, 5LL) - 28992LL*Power(r, 6LL)*Power(xj, 6LL) -

            2064LL*Power(r, 7LL)*Power(xj, 7LL) - 24LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        21LL*Power(xi, 12LL)*Power(xj, 8LL)*

        (482625LL + 868725LL*r*xj + 772200LL*Power(r, 2LL)*Power(xj, 2LL) +

            455400LL*Power(r, 3LL)*Power(xj, 3LL) + 178200LL*Power(r, 4LL)*Power(xj, 4LL) +

            72180LL*Power(r, 5LL)*Power(xj, 5LL) + 19920LL*Power(r, 6LL)*Power(xj, 6LL) +

            2952LL*Power(r, 7LL)*Power(xj, 7LL) + 204LL*Power(r, 8LL)*Power(xj, 8LL) +

            4LL*Power(r, 9LL)*Power(xj, 9LL)) +

        6LL*Power(xi, 6LL)*Power(xj, 14LL)*

        (-10357200LL + 5071815LL*r*xj - 6463800LL*Power(r, 2LL)*Power(xj, 2LL) -

            7151130LL*Power(r, 3LL)*Power(xj, 3LL) -

            2572290LL*Power(r, 4LL)*Power(xj, 4LL) - 468720LL*Power(r, 5LL)*Power(xj, 5LL) -

            42672LL*Power(r, 6LL)*Power(xj, 6LL) - 648LL*Power(r, 7LL)*Power(xj, 7LL) +

            228LL*Power(r, 8LL)*Power(xj, 8LL) + 16LL*Power(r, 9LL)*Power(xj, 9LL)) -

        Power(xi, 18LL)*Power(xj, 2LL)*

        (184275LL + 331695LL*r*xj + 294840LL*Power(r, 2LL)*Power(xj, 2LL) +

            171990LL*Power(r, 3LL)*Power(xj, 3LL) + 73710LL*Power(r, 4LL)*Power(xj, 4LL) +

            24570LL*Power(r, 5LL)*Power(xj, 5LL) + 6552LL*Power(r, 6LL)*Power(xj, 6LL) +

            1404LL*Power(r, 7LL)*Power(xj, 7LL) + 234LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL)) +

        Power(xi, 2LL)*Power(xj, 18LL)*

        (-133783650LL - 107432325LL*r*xj - 35675640LL*Power(r, 2LL)*Power(xj, 2LL) -

            5135130LL*Power(r, 3LL)*Power(xj, 3LL) + 270270LL*Power(r, 4LL)*Power(xj, 4LL) +

            270270LL*Power(r, 5LL)*Power(xj, 5LL) + 57960LL*Power(r, 6LL)*Power(xj, 6LL) +

            6948LL*Power(r, 7LL)*Power(xj, 7LL) + 486LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL)) -

        6LL*Power(xi, 14LL)*Power(xj, 6LL)*

        (675675LL + 1216215LL*r*xj + 1081080LL*Power(r, 2LL)*Power(xj, 2LL) +

            630630LL*Power(r, 3LL)*Power(xj, 3LL) + 270270LL*Power(r, 4LL)*Power(xj, 4LL) +

            88200LL*Power(r, 5LL)*Power(xj, 5LL) + 26544LL*Power(r, 6LL)*Power(xj, 6LL) +

            5160LL*Power(r, 7LL)*Power(xj, 7LL) + 492LL*Power(r, 8LL)*Power(xj, 8LL) +

            16LL*Power(r, 9LL)*Power(xj, 9LL))))/

                (28350LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 13LL)*Power(xi + xj, 13LL))

            ;
        }

    }
    return S;
}


cl_R DSlater_5S_2S(cl_R r, cl_R xi, cl_R xj)
{
    return DSlater_2S_5S(r, xj, xi);
}
