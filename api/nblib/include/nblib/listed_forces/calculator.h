/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \inpublicapi \file
 * \brief
 * Implements a force calculator based on GROMACS data structures.
 *
 * Intended for internal use inside the ForceCalculator.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_LISTEDFORCES_CALCULATOR_H
#define NBLIB_LISTEDFORCES_CALCULATOR_H

#include <cstddef>

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

#include "nblib/basicdefinitions.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/vector.h"

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx

namespace nblib
{
class Box;
class PbcHolder;
template<class T>
class ForceBufferProxy;

/*! \internal \brief Object to calculate forces and energies of listed forces
 *
 */
class ListedForceCalculator
{
public:
    using EnergyType = std::array<real, std::tuple_size<ListedInteractionData>::value>;

    ListedForceCalculator(const ListedInteractionData& interactions,
                          size_t                       bufferSize,
                          int                          numThreads,
                          const Box&                   box);

    /*! \brief Dispatch the listed force kernels and reduce the forces
     *
     * This function adds the computed listed forces to all values in the passed in forces buffer,
     * so it can be regarded as an output only param. In case this is being used in a simulation
     * that uses the same force buffer for both non-bonded and listed forces, this call should be
     * made only after the compute() call from the non-bonded ForceCalculator
     *
     * This function also stores the forces and energies from listed interactions in the internal
     * buffer of the ListedForceCalculator object
     *
     * \param[in]  coordinates     input coordinates for the force calculation
     * \param[inout] forces        output for adding the forces
     * \param[inout] shiftForces   output for adding shift forces
     * \param[out] energies        output for potential energies
     * \param[in]  usePbc          whether or not to consider periodic boundary conditions
     */
    void compute(gmx::ArrayRef<const Vec3> coordinates,
                 gmx::ArrayRef<Vec3>       forces,
                 gmx::ArrayRef<Vec3>       shiftForces,
                 gmx::ArrayRef<real>       energies,
                 bool                      usePbc = false);

    //! \brief Alternative overload without shift forces
    void compute(gmx::ArrayRef<const Vec3> coordinates,
                 gmx::ArrayRef<Vec3>       forces,
                 gmx::ArrayRef<real>       energies,
                 bool                      usePbc = false);

    //! \brief default, but moved to separate compilation unit
    ~ListedForceCalculator();

private:
    int numThreads;

    //! holds the array of energies computed
    EnergyType energyBuffer_;

    //! holds the listed interactions split into groups for multithreading
    std::vector<ListedInteractionData> threadedInteractions_;

    //! reduction force buffers
    std::vector<ForceBufferProxy<Vec3>> threadedForceBuffers_;

    //! reduction shift force buffers
    std::vector<std::vector<Vec3>> threadedShiftForceBuffers_;

    //! PBC objects
    std::unique_ptr<PbcHolder> pbcHolder_;

    //! compute listed forces and energies, overwrites the internal buffers
    template<class ShiftForce>
    void computeForcesAndEnergies(gmx::ArrayRef<const Vec3>                  x,
                                  gmx::ArrayRef<Vec3>                        forces,
                                  [[maybe_unused]] gmx::ArrayRef<ShiftForce> shiftForces,
                                  bool                                       usePbc = false);
};

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CALCULATOR_H
