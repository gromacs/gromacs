/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2017, by the GROMACS development team, led by
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

#include <assert.h>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
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
        real scalingFactor_;  /**< Box The scaling factor PME uses with walls */

    public:
        EwaldBoxZScaler() = delete;

        /*! \brief Constructor that takes the input record to initialize Ewald box scaling appropriately. */
        EwaldBoxZScaler(const t_inputrec &ir)
        {
            if (inputrecPbcXY2Walls(&ir))
            {
                scaleWithWalls_ = true;
                scalingFactor_  = ir.wall_ewald_zfac;
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
        void scaleBox(const matrix box,
                      matrix       scaledBox)
        {
            assert(box);
            assert(scaledBox);

            copy_mat(box, scaledBox);
            if (scaleWithWalls_)
            {
                svmul(scalingFactor_, scaledBox[ZZ], scaledBox[ZZ]);
            }
        }
};

#endif
