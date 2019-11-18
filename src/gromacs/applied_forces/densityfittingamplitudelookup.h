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
/*! \internal \file
 * \brief
 * Declares amplitude lookup for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_DENSITYFITTINGAMPLITUDELOOKUP_H
#define GMX_APPLIED_FORCES_DENSITYFITTINGAMPLITUDELOOKUP_H

#include <map>
#include <memory>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

struct t_mdatoms;

namespace gmx
{

/*! \brief
 * The methods that determine how amplitudes are spread on a grid in density guided simulations.
 */
enum class DensityFittingAmplitudeMethod
{
    Unity,  //!< same spread amplitude, unity, for all atoms
    Mass,   //!< atom mass is the spread amplitude
    Charge, //!< partial charge determines the spread amplitude
    Count
};

//! The names of the methods to determine the amplitude of the atoms to be spread on a grid
const EnumerationArray<DensityFittingAmplitudeMethod, const char* const> c_densityFittingAmplitudeMethodNames = {
    { "unity", "mass", "charge" }
};

class DensityFittingAmplitudeLookupImpl;

/*! \internal \brief Class that translates atom properties into amplitudes.
 *
 */
class DensityFittingAmplitudeLookup
{
public:
    //! Construct force provider for density fitting from its parameters
    explicit DensityFittingAmplitudeLookup(const DensityFittingAmplitudeMethod& method);
    ~DensityFittingAmplitudeLookup();
    //! Copy constructor
    DensityFittingAmplitudeLookup(const DensityFittingAmplitudeLookup& other);
    //! Copy assignment
    DensityFittingAmplitudeLookup& operator=(const DensityFittingAmplitudeLookup& other);
    //! Move constructor
    DensityFittingAmplitudeLookup(DensityFittingAmplitudeLookup&& other) noexcept;
    //! Move assignment
    DensityFittingAmplitudeLookup& operator=(DensityFittingAmplitudeLookup&& other) noexcept;
    /*! \brief Return the amplitudes for spreading atoms of a given local index.
     * \param[in] atoms the atom information
     * \param[in] localIndex the local atom indices
     * \returns amplitudes
     */
    const std::vector<real>& operator()(const t_mdatoms& atoms, ArrayRef<const int> localIndex);

private:
    std::unique_ptr<DensityFittingAmplitudeLookupImpl> impl_;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_DENSITYFITTINGAMPLITUDELOOKUP_H
