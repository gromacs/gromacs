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
            S = -(-131985LL*xi + 161280LL*exp(2LL*r*xi)*xi - 205380LL*r*Pow(xi, 2LL) -

                  149940LL*Pow(r, 2LL)*Pow(xi, 3LL) - 67200LL*Pow(r, 3LL)*Pow(xi, 4LL) -

                  20160LL*Pow(r, 4LL)*Pow(xi, 5LL) - 4032LL*Pow(r, 5LL)*Pow(xi, 6LL) -

                  448LL*Pow(r, 6LL)*Pow(xi, 7LL))/(80640LL*exp(2LL*r*xi)*r) +

                (-80640LL + 80640LL*exp(2LL*r*xi) - 131985LL*r*xi -

                 102690LL*Pow(r, 2LL)*Pow(xi, 2LL) - 49980LL*Pow(r, 3LL)*Pow(xi, 3LL) -

                 16800LL*Pow(r, 4LL)*Pow(xi, 4LL) - 4032LL*Pow(r, 5LL)*Pow(xi, 5LL) -

                 672LL*Pow(r, 6LL)*Pow(xi, 6LL) - 64LL*Pow(r, 7LL)*Pow(xi, 7LL))/

                (80640LL*exp(2LL*r*xi)*Pow(r, 2LL)) +

                (xi*(-80640LL + 80640LL*exp(2LL*r*xi) - 131985LL*r*xi -

                     102690LL*Pow(r, 2LL)*Pow(xi, 2LL) - 49980LL*Pow(r, 3LL)*Pow(xi, 3LL) -

                     16800LL*Pow(r, 4LL)*Pow(xi, 4LL) - 4032LL*Pow(r, 5LL)*Pow(xi, 5LL) -

                     672LL*Pow(r, 6LL)*Pow(xi, 6LL) - 64LL*Pow(r, 7LL)*Pow(xi, 7LL)))/

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
            S = (6LL*exp(2LL*r*(xi + xj))*Pow(Pow(xi, 2LL) - Pow(xj, 2LL), 7LL) -

                 exp(2LL*r*xi)*Pow(xi, 6LL)*

                 (21LL*Pow(xi, 4LL)*Pow(xj, 4LL)*

                  (6LL + 11LL*r*xj + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) -

                  2LL*Pow(xj, 8LL)*(90LL + 54LL*r*xj + 12LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                      Pow(r, 3LL)*Pow(xj, 3LL)) +

                  Pow(xi, 8LL)*(6LL + 9LL*r*xj + 6LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                  2LL*Pow(r, 3LL)*Pow(xj, 3LL)) +

                  Pow(xi, 2LL)*Pow(xj, 6LL)*

                  (-390LL - 69LL*r*xj + 18LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   4LL*Pow(r, 3LL)*Pow(xj, 3LL)) -

                  Pow(xi, 6LL)*Pow(xj, 2LL)*

                  (42LL + 63LL*r*xj + 42LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   4LL*Pow(r, 3LL)*Pow(xj, 3LL))) +

                 exp(2LL*r*xj)*Pow(xj, 6LL)*

                 (-24LL*Pow(r, 2LL)*Pow(xi, 10LL) - 2LL*Pow(r, 3LL)*Pow(xi, 11LL) -

                  69LL*r*Pow(xi, 7LL)*Pow(xj, 2LL) + 6LL*Pow(xj, 8LL) + 9LL*r*xi*Pow(xj, 8LL) +

                  4LL*r*Pow(xi, 9LL)*(-27LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  18LL*Pow(xi, 8LL)*(-10LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  6LL*Pow(xi, 2LL)*Pow(xj, 6LL)*(-7LL + Pow(r, 2LL)*Pow(xj, 2LL)) -

                  42LL*Pow(xi, 4LL)*Pow(xj, 4LL)*(-3LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  r*Pow(xi, 3LL)*Pow(xj, 6LL)*(-63LL + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  6LL*Pow(xi, 6LL)*Pow(xj, 2LL)*(-65LL + 7LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  Pow(xi, 5LL)*(231LL*r*Pow(xj, 4LL) - 4LL*Pow(r, 3LL)*Pow(xj, 6LL))))/

                (6LL*exp(2LL*r*(xi + xj))*Pow(r, 2LL)*Pow(xi - xj, 7LL)*Pow(xi + xj, 7LL)) +

                (6LL*exp(2LL*r*(xi + xj))*Pow(Pow(xi, 2LL) - Pow(xj, 2LL), 7LL) -

                 exp(2LL*r*xi)*Pow(xi, 6LL)*

                 (21LL*Pow(xi, 4LL)*Pow(xj, 4LL)*

                  (6LL + 11LL*r*xj + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) -

                  2LL*Pow(xj, 8LL)*(90LL + 54LL*r*xj + 12LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                      Pow(r, 3LL)*Pow(xj, 3LL)) +

                  Pow(xi, 8LL)*(6LL + 9LL*r*xj + 6LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                  2LL*Pow(r, 3LL)*Pow(xj, 3LL)) +

                  Pow(xi, 2LL)*Pow(xj, 6LL)*

                  (-390LL - 69LL*r*xj + 18LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   4LL*Pow(r, 3LL)*Pow(xj, 3LL)) -

                  Pow(xi, 6LL)*Pow(xj, 2LL)*

                  (42LL + 63LL*r*xj + 42LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   4LL*Pow(r, 3LL)*Pow(xj, 3LL))) +

                 exp(2LL*r*xj)*Pow(xj, 6LL)*

                 (-24LL*Pow(r, 2LL)*Pow(xi, 10LL) - 2LL*Pow(r, 3LL)*Pow(xi, 11LL) -

                  69LL*r*Pow(xi, 7LL)*Pow(xj, 2LL) + 6LL*Pow(xj, 8LL) + 9LL*r*xi*Pow(xj, 8LL) +

                  4LL*r*Pow(xi, 9LL)*(-27LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  18LL*Pow(xi, 8LL)*(-10LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  6LL*Pow(xi, 2LL)*Pow(xj, 6LL)*(-7LL + Pow(r, 2LL)*Pow(xj, 2LL)) -

                  42LL*Pow(xi, 4LL)*Pow(xj, 4LL)*(-3LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  r*Pow(xi, 3LL)*Pow(xj, 6LL)*(-63LL + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  6LL*Pow(xi, 6LL)*Pow(xj, 2LL)*(-65LL + 7LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  Pow(xi, 5LL)*(231LL*r*Pow(xj, 4LL) - 4LL*Pow(r, 3LL)*Pow(xj, 6LL))))/

                (3LL*exp(2LL*r*(xi + xj))*r*Pow(xi - xj, 7LL)*Pow(xi + xj, 6LL)) -

                (12LL*exp(2LL*r*(xi + xj))*(xi + xj)*Pow(Pow(xi, 2LL) - Pow(xj, 2LL), 7LL) -

                 exp(2LL*r*xi)*Pow(xi, 6LL)*

                 (21LL*Pow(xi, 4LL)*Pow(xj, 4LL)*(11LL*xj + 4LL*r*Pow(xj, 2LL)) -

                  2LL*Pow(xj, 8LL)*(54LL*xj + 24LL*r*Pow(xj, 2LL) +

                                      3LL*Pow(r, 2LL)*Pow(xj, 3LL)) +

                  Pow(xi, 8LL)*(9LL*xj + 12LL*r*Pow(xj, 2LL) + 6LL*Pow(r, 2LL)*Pow(xj, 3LL)) +

                  Pow(xi, 2LL)*Pow(xj, 6LL)*

                  (-69LL*xj + 36LL*r*Pow(xj, 2LL) + 12LL*Pow(r, 2LL)*Pow(xj, 3LL)) -

                  Pow(xi, 6LL)*Pow(xj, 2LL)*

                  (63LL*xj + 84LL*r*Pow(xj, 2LL) + 12LL*Pow(r, 2LL)*Pow(xj, 3LL))) -

                 2LL*exp(2LL*r*xi)*Pow(xi, 7LL)*

                 (21LL*Pow(xi, 4LL)*Pow(xj, 4LL)*

                  (6LL + 11LL*r*xj + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) -

                  2LL*Pow(xj, 8LL)*(90LL + 54LL*r*xj + 12LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                      Pow(r, 3LL)*Pow(xj, 3LL)) +

                  Pow(xi, 8LL)*(6LL + 9LL*r*xj + 6LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                                  2LL*Pow(r, 3LL)*Pow(xj, 3LL)) +

                  Pow(xi, 2LL)*Pow(xj, 6LL)*

                  (-390LL - 69LL*r*xj + 18LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   4LL*Pow(r, 3LL)*Pow(xj, 3LL)) -

                  Pow(xi, 6LL)*Pow(xj, 2LL)*

                  (42LL + 63LL*r*xj + 42LL*Pow(r, 2LL)*Pow(xj, 2LL) +

                   4LL*Pow(r, 3LL)*Pow(xj, 3LL))) +

                 exp(2LL*r*xj)*Pow(xj, 6LL)*

                 (-48LL*r*Pow(xi, 10LL) - 6LL*Pow(r, 2LL)*Pow(xi, 11LL) -

                  69LL*Pow(xi, 7LL)*Pow(xj, 2LL) + 36LL*r*Pow(xi, 8LL)*Pow(xj, 2LL) +

                  8LL*Pow(r, 2LL)*Pow(xi, 9LL)*Pow(xj, 2LL) +

                  84LL*r*Pow(xi, 6LL)*Pow(xj, 4LL) - 84LL*r*Pow(xi, 4LL)*Pow(xj, 6LL) +

                  9LL*xi*Pow(xj, 8LL) + 12LL*r*Pow(xi, 2LL)*Pow(xj, 8LL) +

                  4LL*Pow(r, 2LL)*Pow(xi, 3LL)*Pow(xj, 8LL) +

                  4LL*Pow(xi, 9LL)*(-27LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  Pow(xi, 3LL)*Pow(xj, 6LL)*(-63LL + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  Pow(xi, 5LL)*(231LL*Pow(xj, 4LL) - 12LL*Pow(r, 2LL)*Pow(xj, 6LL))) +

                 2LL*exp(2LL*r*xj)*Pow(xj, 7LL)*

                 (-24LL*Pow(r, 2LL)*Pow(xi, 10LL) - 2LL*Pow(r, 3LL)*Pow(xi, 11LL) -

                  69LL*r*Pow(xi, 7LL)*Pow(xj, 2LL) + 6LL*Pow(xj, 8LL) + 9LL*r*xi*Pow(xj, 8LL) +

                  4LL*r*Pow(xi, 9LL)*(-27LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  18LL*Pow(xi, 8LL)*(-10LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  6LL*Pow(xi, 2LL)*Pow(xj, 6LL)*(-7LL + Pow(r, 2LL)*Pow(xj, 2LL)) -

                  42LL*Pow(xi, 4LL)*Pow(xj, 4LL)*(-3LL + Pow(r, 2LL)*Pow(xj, 2LL)) +

                  r*Pow(xi, 3LL)*Pow(xj, 6LL)*(-63LL + 2LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  6LL*Pow(xi, 6LL)*Pow(xj, 2LL)*(-65LL + 7LL*Pow(r, 2LL)*Pow(xj, 2LL)) +

                  Pow(xi, 5LL)*(231LL*r*Pow(xj, 4LL) - 4LL*Pow(r, 3LL)*Pow(xj, 6LL))))/

                (6LL*exp(2LL*r*(xi + xj))*r*Pow(xi - xj, 7LL)*Pow(xi + xj, 7LL))

            ;
        }

    }
    return S;
}

#else

double DSlater_2S_2S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = -(-131985*xi + 161280*exp(2*r*xi)*xi - 205380*r*pow(xi, 2) -

                  149940*pow(r, 2)*pow(xi, 3) - 67200*pow(r, 3)*pow(xi, 4) -

                  20160*pow(r, 4)*pow(xi, 5) - 4032*pow(r, 5)*pow(xi, 6) -

                  448*pow(r, 6)*pow(xi, 7))/(80640*exp(2*r*xi)*r) +

                (-80640 + 80640*exp(2*r*xi) - 131985*r*xi -

                 102690*pow(r, 2)*pow(xi, 2) - 49980*pow(r, 3)*pow(xi, 3) -

                 16800*pow(r, 4)*pow(xi, 4) - 4032*pow(r, 5)*pow(xi, 5) -

                 672*pow(r, 6)*pow(xi, 6) - 64*pow(r, 7)*pow(xi, 7))/

                (80640*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-80640 + 80640*exp(2*r*xi) - 131985*r*xi -

                     102690*pow(r, 2)*pow(xi, 2) - 49980*pow(r, 3)*pow(xi, 3) -

                     16800*pow(r, 4)*pow(xi, 4) - 4032*pow(r, 5)*pow(xi, 5) -

                     672*pow(r, 6)*pow(xi, 6) - 64*pow(r, 7)*pow(xi, 7)))/

                (40320*exp(2*r*xi)*r)

            ;
        }

    }
    else
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = (6*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) -

                 exp(2*r*xi)*pow(xi, 6)*

                 (21*pow(xi, 4)*pow(xj, 4)*

                  (6 + 11*r*xj + 2*pow(r, 2)*pow(xj, 2)) -

                  2*pow(xj, 8)*(90 + 54*r*xj + 12*pow(r, 2)*pow(xj, 2) +

                                      pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 8)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-390 - 69*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (42 + 63*r*xj + 42*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3))) +

                 exp(2*r*xj)*pow(xj, 6)*

                 (-24*pow(r, 2)*pow(xi, 10) - 2*pow(r, 3)*pow(xi, 11) -

                  69*r*pow(xi, 7)*pow(xj, 2) + 6*pow(xj, 8) + 9*r*xi*pow(xj, 8) +

                  4*r*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  18*pow(xi, 8)*(-10 + pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 2)*pow(xj, 6)*(-7 + pow(r, 2)*pow(xj, 2)) -

                  42*pow(xi, 4)*pow(xj, 4)*(-3 + pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 6)*pow(xj, 2)*(-65 + 7*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*r*pow(xj, 4) - 4*pow(r, 3)*pow(xj, 6))))/

                (6*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 7)*pow(xi + xj, 7)) +

                (6*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) -

                 exp(2*r*xi)*pow(xi, 6)*

                 (21*pow(xi, 4)*pow(xj, 4)*

                  (6 + 11*r*xj + 2*pow(r, 2)*pow(xj, 2)) -

                  2*pow(xj, 8)*(90 + 54*r*xj + 12*pow(r, 2)*pow(xj, 2) +

                                      pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 8)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-390 - 69*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (42 + 63*r*xj + 42*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3))) +

                 exp(2*r*xj)*pow(xj, 6)*

                 (-24*pow(r, 2)*pow(xi, 10) - 2*pow(r, 3)*pow(xi, 11) -

                  69*r*pow(xi, 7)*pow(xj, 2) + 6*pow(xj, 8) + 9*r*xi*pow(xj, 8) +

                  4*r*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  18*pow(xi, 8)*(-10 + pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 2)*pow(xj, 6)*(-7 + pow(r, 2)*pow(xj, 2)) -

                  42*pow(xi, 4)*pow(xj, 4)*(-3 + pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 6)*pow(xj, 2)*(-65 + 7*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*r*pow(xj, 4) - 4*pow(r, 3)*pow(xj, 6))))/

                (3*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 6)) -

                (12*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 7) -

                 exp(2*r*xi)*pow(xi, 6)*

                 (21*pow(xi, 4)*pow(xj, 4)*(11*xj + 4*r*pow(xj, 2)) -

                  2*pow(xj, 8)*(54*xj + 24*r*pow(xj, 2) +

                                      3*pow(r, 2)*pow(xj, 3)) +

                  pow(xi, 8)*(9*xj + 12*r*pow(xj, 2) + 6*pow(r, 2)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-69*xj + 36*r*pow(xj, 2) + 12*pow(r, 2)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (63*xj + 84*r*pow(xj, 2) + 12*pow(r, 2)*pow(xj, 3))) -

                 2*exp(2*r*xi)*pow(xi, 7)*

                 (21*pow(xi, 4)*pow(xj, 4)*

                  (6 + 11*r*xj + 2*pow(r, 2)*pow(xj, 2)) -

                  2*pow(xj, 8)*(90 + 54*r*xj + 12*pow(r, 2)*pow(xj, 2) +

                                      pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 8)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-390 - 69*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (42 + 63*r*xj + 42*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3))) +

                 exp(2*r*xj)*pow(xj, 6)*

                 (-48*r*pow(xi, 10) - 6*pow(r, 2)*pow(xi, 11) -

                  69*pow(xi, 7)*pow(xj, 2) + 36*r*pow(xi, 8)*pow(xj, 2) +

                  8*pow(r, 2)*pow(xi, 9)*pow(xj, 2) +

                  84*r*pow(xi, 6)*pow(xj, 4) - 84*r*pow(xi, 4)*pow(xj, 6) +

                  9*xi*pow(xj, 8) + 12*r*pow(xi, 2)*pow(xj, 8) +

                  4*pow(r, 2)*pow(xi, 3)*pow(xj, 8) +

                  4*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*pow(xj, 4) - 12*pow(r, 2)*pow(xj, 6))) +

                 2*exp(2*r*xj)*pow(xj, 7)*

                 (-24*pow(r, 2)*pow(xi, 10) - 2*pow(r, 3)*pow(xi, 11) -

                  69*r*pow(xi, 7)*pow(xj, 2) + 6*pow(xj, 8) + 9*r*xi*pow(xj, 8) +

                  4*r*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  18*pow(xi, 8)*(-10 + pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 2)*pow(xj, 6)*(-7 + pow(r, 2)*pow(xj, 2)) -

                  42*pow(xi, 4)*pow(xj, 4)*(-3 + pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 6)*pow(xj, 2)*(-65 + 7*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*r*pow(xj, 4) - 4*pow(r, 3)*pow(xj, 6))))/

                (6*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 7))

            ;
        }

    }
    return S;
}

#endif
