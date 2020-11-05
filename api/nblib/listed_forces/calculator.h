/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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

#include <memory>
#include <unordered_map>

#include "nblib/listed_forces/definitions.h"

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
class ForceBuffer;

/*! \internal \brief object to calculate listed forces
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
     * \param[in] coordinates to be used for the force calculation
     * \param[out] forces buffer to store the output forces
     */
    void compute(gmx::ArrayRef<const Vec3> coordinates, gmx::ArrayRef<Vec3> forces, bool usePbc = false);

    //! Alternative overload with the energies in an output buffer
    void compute(gmx::ArrayRef<const Vec3> coordinates,
                 gmx::ArrayRef<Vec3>       forces,
                 EnergyType&               energies,
                 bool                      usePbc = false);

    /*! \brief We need to declare the destructor here to move the (still default) implementation
     *  to the .cpp file. Omitting this declaration would mean an inline destructor
     *  which can't compile because the unique_ptr dtor needs ~ForceBuffer, which is not available
     * here because it's incomplete.
     */
    ~ListedForceCalculator();

private:
    int numThreads;

    //! the main buffer to hold the final listed forces
    std::vector<gmx::RVec> masterForceBuffer_;

    //! holds the array of energies computed
    EnergyType energyBuffer_;

    //! holds the listed interactions split into groups for multithreading
    std::vector<ListedInteractionData> threadedInteractions_;

    //! reduction force buffers
    std::vector<std::unique_ptr<ForceBuffer<gmx::RVec>>> threadedForceBuffers_;

    //! PBC objects
    std::unique_ptr<PbcHolder> pbcHolder_;

    //! compute listed forces and energies, overwrites the internal buffers
    void computeForcesAndEnergies(gmx::ArrayRef<const Vec3> x, bool usePbc = false);
};

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CALCULATOR_H
