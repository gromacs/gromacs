/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal \brief Container for force and virial for algorithms that compute shift forces for virial calculation
 *
 * This force output class should be used when computing forces whos virial contribution
 * is computed using the "single sum virial" algorithm (see the reference manual for
 * details). To handle the virial contributions of forces working between periodic
 * images correctly, so-called "shift forces" need to be accumulated for the different
 * periodic images.
 */
class ForceWithShiftForces
{
public:
    /*! \brief Constructor
     *
     * \param[in] force          A force buffer that will be used for storing forces
     * \param[in] computeVirial  True when algorithms are required to provide their virial contribution (for the current force evaluation)
     * \param[in] shiftForces    A shift forces buffer of size SHIFTS, only used with \p computeVirial = true
     */
    ForceWithShiftForces(const gmx::ArrayRefWithPadding<gmx::RVec>& force,
                         const bool                                 computeVirial,
                         const gmx::ArrayRef<gmx::RVec>&            shiftForces) :
        force_(force),
        computeVirial_(computeVirial),
        shiftForces_(computeVirial ? shiftForces : gmx::ArrayRef<gmx::RVec>())
    {
    }

    //! Returns an arrayref to the force buffer without padding
    gmx::ArrayRef<gmx::RVec> force() { return force_.unpaddedArrayRef(); }

    //! Returns a const arrayref to the force buffer without padding
    gmx::ArrayRef<const gmx::RVec> force() const { return force_.unpaddedConstArrayRef(); }

    //! Returns whether the virial needs to be computed
    bool computeVirial() const { return computeVirial_; }

    //! Returns the shift forces buffer
    gmx::ArrayRef<gmx::RVec> shiftForces() { return shiftForces_; }

    //! Returns a const shift forces buffer
    gmx::ArrayRef<const gmx::RVec> shiftForces() const { return shiftForces_; }

private:
    //! The force buffer
    gmx::ArrayRefWithPadding<gmx::RVec> force_;
    //! True when virial computation is requested
    bool computeVirial_;
    //! A buffer for storing the shift forces, size SHIFTS
    gmx::ArrayRef<gmx::RVec> shiftForces_;
};

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
    ForceWithVirial(const ArrayRef<RVec>& force, const bool computeVirial) :
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
    const matrix& getVirial() const { return virial_; }

    const ArrayRef<RVec> force_;         //!< Force accumulation buffer reference
    const bool           computeVirial_; //!< True when algorithms are required to provide their virial contribution (for the current force evaluation)
private:
    matrix virial_; //!< Virial accumulation buffer
};

/*! \libinternal \brief Force and virial output buffers for use in force computation
 */
class ForceOutputs
{
public:
    //! Constructor
    ForceOutputs(const ForceWithShiftForces& forceWithShiftForces, const ForceWithVirial& forceWithVirial) :
        forceWithShiftForces_(forceWithShiftForces),
        forceWithVirial_(forceWithVirial)
    {
    }

    //! Returns a reference to the force with shift forces object
    ForceWithShiftForces& forceWithShiftForces() { return forceWithShiftForces_; }

    //! Returns a reference to the force with virial object
    ForceWithVirial& forceWithVirial() { return forceWithVirial_; }

private:
    //! Force output buffer used by legacy modules (without SIMD padding)
    ForceWithShiftForces forceWithShiftForces_;
    //! Force with direct virial contribution (if there are any; without SIMD padding)
    ForceWithVirial forceWithVirial_;
};

} // namespace gmx

#endif
