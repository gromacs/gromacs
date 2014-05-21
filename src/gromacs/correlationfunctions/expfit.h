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
#ifndef GMX_EXPFIT_H
#define GMX_EXPFIT_H

#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Enum to select fitting functions
 */
enum {
    effnNONE, effnEXP1, effnEXP2, effnEXP3,   effnVAC,
    effnEXP5, effnEXP7, effnEXP9, effnERF, effnERREST, effnNR
};

extern const int   nfp_ffn[effnNR];

extern const char *s_ffn[effnNR+2];

extern const char *longs_ffn[effnNR];


/*! \brief
 * Returns  corresponding to the selected enum option in sffn
 * \param[in] sffn Two dimensional string array coming from parse_common_args
 * \return the ffn enum
 */
int sffn2effn(const char **sffn);

/*! \brief
 * Fitting to exponential?
 *
 ***********************************************/
/*void do_expfit(int ndata, real c1[], real dt,
               real begintimefit, real endtimefit);
 */
/*! \brief
 * This procedure fits y=exp(a+bx) for n (x,y) pairs to determine a and b.
 * The standard deviations of a and b, sa and sb, are also calculated.
 *
 * Routine from Computers in physics, 7(3) (1993), p. 280-285.
 * \param[in] n Number of data points in the x, y, dy arrays
 * \param[in] x The x-axis
 * \param[in] y The y-axis
 * \param[in] Dy the uncertainties in the y values.
 * \param[out] a The coefficient a in the eqn y = exp(a+bx)
 * \param[out] sa The uncertainty in a
 * \param[out] b The coefficient b in the eqn y = exp(a+bx)
 * \param[out] sb The uncertainty in b
 */
void expfit(int n, real x[], real y[], real Dy[],
            real *a, real *sa,
            real *b, real *sb);


#ifdef __cplusplus
}
#endif

#endif
