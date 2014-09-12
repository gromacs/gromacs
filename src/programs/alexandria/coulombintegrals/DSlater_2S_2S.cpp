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

cl_R DSlater_2S_2S(cl_R r, cl_R xi, cl_R xj)
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
            S = -(-131985LL*xi + 161280LL*exp(2LL*r*xi)*xi - 205380LL*r*Power(xi, 2LL) -

                  149940LL*Power(r, 2LL)*Power(xi, 3LL) - 67200LL*Power(r, 3LL)*Power(xi, 4LL) -

                  20160LL*Power(r, 4LL)*Power(xi, 5LL) - 4032LL*Power(r, 5LL)*Power(xi, 6LL) -

                  448LL*Power(r, 6LL)*Power(xi, 7LL))/(80640LL*exp(2LL*r*xi)*r) +

                (-80640LL + 80640LL*exp(2LL*r*xi) - 131985LL*r*xi -

                 102690LL*Power(r, 2LL)*Power(xi, 2LL) - 49980LL*Power(r, 3LL)*Power(xi, 3LL) -

                 16800LL*Power(r, 4LL)*Power(xi, 4LL) - 4032LL*Power(r, 5LL)*Power(xi, 5LL) -

                 672LL*Power(r, 6LL)*Power(xi, 6LL) - 64LL*Power(r, 7LL)*Power(xi, 7LL))/

                (80640LL*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-80640LL + 80640LL*exp(2LL*r*xi) - 131985LL*r*xi -

                     102690LL*Power(r, 2LL)*Power(xi, 2LL) - 49980LL*Power(r, 3LL)*Power(xi, 3LL) -

                     16800LL*Power(r, 4LL)*Power(xi, 4LL) - 4032LL*Power(r, 5LL)*Power(xi, 5LL) -

                     672LL*Power(r, 6LL)*Power(xi, 6LL) - 64LL*Power(r, 7LL)*Power(xi, 7LL)))/

                (40320LL*exp(2LL*r*xi)*r)

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
            S = (6LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 7LL) -

                 exp(2LL*r*xi)*Power(xi, 6LL)*

                 (21LL*Power(xi, 4LL)*Power(xj, 4LL)*

                  (6LL + 11LL*r*xj + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

                  2LL*Power(xj, 8LL)*(90LL + 54LL*r*xj + 12LL*Power(r, 2LL)*Power(xj, 2LL) +

                                      Power(r, 3LL)*Power(xj, 3LL)) +

                  Power(xi, 8LL)*(6LL + 9LL*r*xj + 6LL*Power(r, 2LL)*Power(xj, 2LL) +

                                  2LL*Power(r, 3LL)*Power(xj, 3LL)) +

                  Power(xi, 2LL)*Power(xj, 6LL)*

                  (-390LL - 69LL*r*xj + 18LL*Power(r, 2LL)*Power(xj, 2LL) +

            4LL*Power(r, 3LL)*Power(xj, 3LL)) -

                  Power(xi, 6LL)*Power(xj, 2LL)*

                  (42LL + 63LL*r*xj + 42LL*Power(r, 2LL)*Power(xj, 2LL) +

            4LL*Power(r, 3LL)*Power(xj, 3LL))) +

                 exp(2LL*r*xj)*Power(xj, 6LL)*

                 (-24LL*Power(r, 2LL)*Power(xi, 10LL) - 2LL*Power(r, 3LL)*Power(xi, 11LL) -

                  69LL*r*Power(xi, 7LL)*Power(xj, 2LL) + 6LL*Power(xj, 8LL) + 9LL*r*xi*Power(xj, 8LL) +

                  4LL*r*Power(xi, 9LL)*(-27LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  18LL*Power(xi, 8LL)*(-10LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  6LL*Power(xi, 2LL)*Power(xj, 6LL)*(-7LL + Power(r, 2LL)*Power(xj, 2LL)) -

                  42LL*Power(xi, 4LL)*Power(xj, 4LL)*(-3LL + Power(r, 2LL)*Power(xj, 2LL)) +

                  r*Power(xi, 3LL)*Power(xj, 6LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  6LL*Power(xi, 6LL)*Power(xj, 2LL)*(-65LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

                  Power(xi, 5LL)*(231LL*r*Power(xj, 4LL) - 4LL*Power(r, 3LL)*Power(xj, 6LL))))/

                (6LL*exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 7LL)*Power(xi + xj, 7LL)) +

                (6LL*exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 7LL) -

                 exp(2LL*r*xi)*Power(xi, 6LL)*

                 (21LL*Power(xi, 4LL)*Power(xj, 4LL)*

        (6LL + 11LL*r*xj + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        2LL*Power(xj, 8LL)*(90LL + 54LL*r*xj + 12LL*Power(r, 2LL)*Power(xj, 2LL) +

                            Power(r, 3LL)*Power(xj, 3LL)) +

        Power(xi, 8LL)*(6LL + 9LL*r*xj + 6LL*Power(r, 2LL)*Power(xj, 2LL) +

                        2LL*Power(r, 3LL)*Power(xj, 3LL)) +

        Power(xi, 2LL)*Power(xj, 6LL)*

        (-390LL - 69LL*r*xj + 18LL*Power(r, 2LL)*Power(xj, 2LL) +

            4LL*Power(r, 3LL)*Power(xj, 3LL)) -

        Power(xi, 6LL)*Power(xj, 2LL)*

        (42LL + 63LL*r*xj + 42LL*Power(r, 2LL)*Power(xj, 2LL) +

            4LL*Power(r, 3LL)*Power(xj, 3LL))) +

                 exp(2LL*r*xj)*Power(xj, 6LL)*

                 (-24LL*Power(r, 2LL)*Power(xi, 10LL) - 2LL*Power(r, 3LL)*Power(xi, 11LL) -

        69LL*r*Power(xi, 7LL)*Power(xj, 2LL) + 6LL*Power(xj, 8LL) + 9LL*r*xi*Power(xj, 8LL) +

        4LL*r*Power(xi, 9LL)*(-27LL + Power(r, 2LL)*Power(xj, 2LL)) +

        18LL*Power(xi, 8LL)*(-10LL + Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 2LL)*Power(xj, 6LL)*(-7LL + Power(r, 2LL)*Power(xj, 2LL)) -

        42LL*Power(xi, 4LL)*Power(xj, 4LL)*(-3LL + Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 3LL)*Power(xj, 6LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 6LL)*Power(xj, 2LL)*(-65LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 5LL)*(231LL*r*Power(xj, 4LL) - 4LL*Power(r, 3LL)*Power(xj, 6LL))))/

                (3LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 7LL)*Power(xi + xj, 6LL)) -

                (12LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi, 2LL) - Power(xj, 2LL), 7LL) -

                 exp(2LL*r*xi)*Power(xi, 6LL)*

                 (21LL*Power(xi, 4LL)*Power(xj, 4LL)*(11LL*xj + 4LL*r*Power(xj, 2LL)) -

        2LL*Power(xj, 8LL)*(54LL*xj + 24LL*r*Power(xj, 2LL) +

                            3LL*Power(r, 2LL)*Power(xj, 3LL)) +

        Power(xi, 8LL)*(9LL*xj + 12LL*r*Power(xj, 2LL) + 6LL*Power(r, 2LL)*Power(xj, 3LL)) +

        Power(xi, 2LL)*Power(xj, 6LL)*

        (-69LL*xj + 36LL*r*Power(xj, 2LL) + 12LL*Power(r, 2LL)*Power(xj, 3LL)) -

        Power(xi, 6LL)*Power(xj, 2LL)*

        (63LL*xj + 84LL*r*Power(xj, 2LL) + 12LL*Power(r, 2LL)*Power(xj, 3LL))) -

                 2LL*exp(2LL*r*xi)*Power(xi, 7LL)*

                 (21LL*Power(xi, 4LL)*Power(xj, 4LL)*

        (6LL + 11LL*r*xj + 2LL*Power(r, 2LL)*Power(xj, 2LL)) -

        2LL*Power(xj, 8LL)*(90LL + 54LL*r*xj + 12LL*Power(r, 2LL)*Power(xj, 2LL) +

                            Power(r, 3LL)*Power(xj, 3LL)) +

        Power(xi, 8LL)*(6LL + 9LL*r*xj + 6LL*Power(r, 2LL)*Power(xj, 2LL) +

                        2LL*Power(r, 3LL)*Power(xj, 3LL)) +

        Power(xi, 2LL)*Power(xj, 6LL)*

        (-390LL - 69LL*r*xj + 18LL*Power(r, 2LL)*Power(xj, 2LL) +

            4LL*Power(r, 3LL)*Power(xj, 3LL)) -

        Power(xi, 6LL)*Power(xj, 2LL)*

        (42LL + 63LL*r*xj + 42LL*Power(r, 2LL)*Power(xj, 2LL) +

            4LL*Power(r, 3LL)*Power(xj, 3LL))) +

                 exp(2LL*r*xj)*Power(xj, 6LL)*

                 (-48LL*r*Power(xi, 10LL) - 6LL*Power(r, 2LL)*Power(xi, 11LL) -

        69LL*Power(xi, 7LL)*Power(xj, 2LL) + 36LL*r*Power(xi, 8LL)*Power(xj, 2LL) +

        8LL*Power(r, 2LL)*Power(xi, 9LL)*Power(xj, 2LL) +

        84LL*r*Power(xi, 6LL)*Power(xj, 4LL) - 84LL*r*Power(xi, 4LL)*Power(xj, 6LL) +

        9LL*xi*Power(xj, 8LL) + 12LL*r*Power(xi, 2LL)*Power(xj, 8LL) +

        4LL*Power(r, 2LL)*Power(xi, 3LL)*Power(xj, 8LL) +

        4LL*Power(xi, 9LL)*(-27LL + Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 3LL)*Power(xj, 6LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 5LL)*(231LL*Power(xj, 4LL) - 12LL*Power(r, 2LL)*Power(xj, 6LL))) +

                 2LL*exp(2LL*r*xj)*Power(xj, 7LL)*

                 (-24LL*Power(r, 2LL)*Power(xi, 10LL) - 2LL*Power(r, 3LL)*Power(xi, 11LL) -

        69LL*r*Power(xi, 7LL)*Power(xj, 2LL) + 6LL*Power(xj, 8LL) + 9LL*r*xi*Power(xj, 8LL) +

        4LL*r*Power(xi, 9LL)*(-27LL + Power(r, 2LL)*Power(xj, 2LL)) +

        18LL*Power(xi, 8LL)*(-10LL + Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 2LL)*Power(xj, 6LL)*(-7LL + Power(r, 2LL)*Power(xj, 2LL)) -

        42LL*Power(xi, 4LL)*Power(xj, 4LL)*(-3LL + Power(r, 2LL)*Power(xj, 2LL)) +

        r*Power(xi, 3LL)*Power(xj, 6LL)*(-63LL + 2LL*Power(r, 2LL)*Power(xj, 2LL)) +

        6LL*Power(xi, 6LL)*Power(xj, 2LL)*(-65LL + 7LL*Power(r, 2LL)*Power(xj, 2LL)) +

        Power(xi, 5LL)*(231LL*r*Power(xj, 4LL) - 4LL*Power(r, 3LL)*Power(xj, 6LL))))/

                (6LL*exp(2LL*r*(xi + xj))*r*Power(xi - xj, 7LL)*Power(xi + xj, 7LL))

            ;
        }

    }
    return S;
}
