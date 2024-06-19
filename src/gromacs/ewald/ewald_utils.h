/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 *
 * \brief Declares utility functions related to Ewald.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_UTILS_H
#define GMX_EWALD_UTILS_H

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

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
real calc_ewaldcoeff_q(real rc, real rtol);

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
real calc_ewaldcoeff_lj(real rc, real rtol);


/*! \libinternal \brief Class to handle box scaling for Ewald and PME.
 *
 * At construction contents of inputrec determine whether scaling is necessary
 * as well as the scaling factor used. Later, the scaleBox method can be used
 * to apply the appropriate scaling (if needed) for Ewald-based methods.
 *
 */
class EwaldBoxZScaler
{

private:
    bool scaleWithWalls_; /**< True if the simulation uses two walls and the box needs to be scaled in PME */
    real scalingFactor_; /**< Box The scaling factor PME uses with walls */

public:
    EwaldBoxZScaler() = delete;

    /*! \brief Constructor that takes the input record to initialize Ewald box scaling appropriately. */
    EwaldBoxZScaler(bool havePbcXY2Walls, real wallEwaldZfac)
    {
        if (havePbcXY2Walls)
        {
            scaleWithWalls_ = true;
            scalingFactor_  = wallEwaldZfac;
        }
        else
        {
            scaleWithWalls_ = false;
            scalingFactor_  = 1;
        }
    }

    /*! \brief Copy and scale the box for PME.
     *
     * When PME is used with 2D periodicity and two walls, the
     * copy of the \p box passed is scaled with the Z scaling factor.
     *
     * \param[in] box        The current box matrix
     * \param[out] scaledBox Scaled copy of the box matrix.
     */
    void scaleBox(const matrix box, matrix scaledBox) const
    {
        GMX_ASSERT(box, "invalid source box pointer");
        GMX_ASSERT(scaledBox, "invalid target box pointer");

        copy_mat(box, scaledBox);
        if (scaleWithWalls_)
        {
            svmul(scalingFactor_, scaledBox[ZZ], scaledBox[ZZ]);
        }
    }
};

#endif
