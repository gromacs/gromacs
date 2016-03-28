/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief
 * This file contains the definition of a container for force and virial
 * output.
 *
 * Currently the only container defined here is one used in algorithms
 * that provide their own virial tensor contribution.
 * We can consider adding another containter for forces and shift forces.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_FORCEOUTPUT_H
#define GMX_MDTYPES_FORCEOUTPUT_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal \brief Container for force and virial for algorithms that provide their own virial tensor contribution
 *
 * \note The \p force_ data member is a reference to an external force buffer.
 */
class ForceWithVirial
{
    public:
        /*! \brief Constructor
         *
         * \param[in] force          A force buffer that will be used for storing forces
         * \param[in] computeVirial  True when algorithms are required to provide their virial contribution (for the current force evaluation)
         */
        ForceWithVirial(ArrayRef<RVec> force, bool computeVirial) :
            force_(force),
            computeVirial_(computeVirial)
        {
            for (int dim1 = 0; dim1 < DIM; dim1++)
            {
                for (int dim2 = 0; dim2 < DIM; dim2++)
                {
                    virial_[dim1][dim2] = 0;
                }
            }
        }

        /*! \brief Adds a virial contribution
         *
         * \note Can be called with \p computeVirial=false.
         * \note It is recommended to accumulate the virial contributions
         *       of a module internally before calling this method, as that
         *       will reduce rounding errors.
         *
         * \param[in] virial  The virial contribution to add
         */
        void addVirialContribution(const matrix virial)
        {
            if (computeVirial_)
            {
                for (int dim1 = 0; dim1 < DIM; dim1++)
                {
                    for (int dim2 = 0; dim2 < DIM; dim2++)
                    {
                        virial_[dim1][dim2] += virial[dim1][dim2];
                    }
                }
            }
        }

        /*! \brief Adds a virial diagonal contribution
         *
         * \note Can be called with \p computeVirial=false.
         * \note It is recommended to accumulate the virial contributions
         *       of a module internally before calling this method, as that
         *       will reduce rounding errors.
         *
         * \param[in] virial  The virial contribution to add
         */
        void addVirialContribution(const RVec virial)
        {
            if (computeVirial_)
            {
                for (int dim = 0; dim < DIM; dim++)
                {
                    virial_[dim][dim] += virial[dim];
                }
            }
        }

        /*! \brief Returns the accumulated virial contributions
         */
        const matrix &getVirial() const
        {
            return virial_;
        }

        const ArrayRef<RVec> force_;         //!< Force accumulation buffer reference
        const bool           computeVirial_; //!< True when algorithms are required to provide their virial contribution (for the current force evaluation)
    private:
        matrix               virial_;        //!< Virial accumulation buffer
};

}

#endif
