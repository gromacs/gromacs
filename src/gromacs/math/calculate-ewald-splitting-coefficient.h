/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \file
 *
 * \brief Declares functions for computing Ewald splitting coefficients
 *
 * These belong in the maths module because they do simple maths and
 * are used many parts of Gromacs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 */
#ifndef GMX_MATH_CALCULATE_EWALD_SPLITTING_COEFFICIENT_H
#define GMX_MATH_CALCULATE_EWALD_SPLITTING_COEFFICIENT_H

#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Computes the Ewald splitting coefficient for Coulomb
 *
 * Returns a value of beta that satisfies rtol > erfc(beta * rc)
 * (and is very close to equality). That value is used the same way in
 * all Coulomb-based Ewald methods.
 *
 * \param[in] rc    Cutoff radius
 * \param[in] rtol  Required maximum value of the short-ranged
 *                  potential at the cutoff (ie. ewald-rtol)
 * \return          The value of the splitting coefficient that
 *                  produces the required dtol at rc.
 */
real
calc_ewaldcoeff_q(real rc, real rtol);

/*! \brief Computes the Ewald splitting coefficient for LJ
 *
 * Returns a value of beta that satisfies dtol > erfc(beta * rc) * (1
 * + beta^2 * rc^2 + 0.5 * beta^4 * rc^4) (and is very close to
 * equality), which is used in LJ-PME.
 *
 * \param[in] rc    Cutoff radius
 * \param[in] rtol  Required maximum value of the short-ranged
 *                  potential at the cutoff (ie. ewald-rtol-lj)
 * \return          The value of the splitting coefficient that
 *                  produces the required dtol at rc.
 */
real
calc_ewaldcoeff_lj(real rc, real rtol);

#ifdef __cplusplus
}
#endif

#endif
