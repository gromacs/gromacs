/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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

/*! \internal \file
 *
 * \brief
 * Declares functions to check bias sharing properties.
 *
 * This actual sharing of biases is currently implemeted in BiasState.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_BIASSHARING_H
#define GMX_AWH_BIASSHARING_H

#include <cstddef>

#include <memory>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxmpi.h"

struct gmx_multisim_t;
struct t_commrec;

namespace gmx
{

template<typename>
class ArrayRef;
class AwhParams;

class BiasSharing
{
public:
    /*! \brief Constructor
     *
     * \param[in] awhParams              Parameters for all biases in this simulation
     * \param[in] commRecord             Intra-simulation communication record
     * \param[in] simulationMainComm     MPI communicator for all main ranks of all simulations that share this bias
     */
    BiasSharing(const AwhParams& awhParams, const t_commrec& commRecord, MPI_Comm simulationMainComm);

    ~BiasSharing();

    //! Returns the number of simulations sharing bias \p biasIndex
    int numSharingSimulations(int biasIndex) const { return numSharingSimulations_[biasIndex]; }

    //! Returns the index of our simulation in the simulations sharing bias \p biasIndex
    int sharingSimulationIndex(int biasIndex) const { return sharingSimulationIndices_[biasIndex]; }

    //! Sums data of type int over the main ranks of all simulations sharing bias \p biasIndex
    void sumOverSharingMainRanks(ArrayRef<int> data, int biasIndex) const;

    //! Sums data of type long over the main ranks of all simulations sharing bias \p biasIndex
    void sumOverSharingMainRanks(ArrayRef<long> data, int biasIndex) const;

    //! Sums data of type double over the main ranks of all simulations sharing bias \p biasIndex
    void sumOverSharingMainRanks(ArrayRef<double> data, int biasIndex) const;

    /*! \brief Sums data of type int over all simulations sharing bias \p biasIndex
     *
     * The summing is performed over the main ranks of the simulations sharing bias \p biasIndex
     * and the result is broadcasted to the other ranks within each simulation.
     */
    void sumOverSharingSimulations(ArrayRef<int> data, int biasIndex) const;

    /*! \brief Sums data of type long over all simulations sharing bias \p biasIndex
     *
     * The summing is performed over the main ranks of the simulations sharing bias \p biasIndex
     * and the result is broadcasted to the other ranks within each simulation.
     */
    void sumOverSharingSimulations(ArrayRef<double> data, int biasIndex) const;

private:
    //! The number of simulations sharing for each bias
    std::vector<int> numSharingSimulations_;
    //! The index of our simulations in the simulations for each bias
    std::vector<int> sharingSimulationIndices_;
    //! Reference to the intra-simulation communication record
    const t_commrec& commRecord_;

    //! Communicator between main ranks sharing a bias, for each bias
    std::vector<MPI_Comm> multiSimCommPerBias_;
    //! List of MPI communicators created by this object so we can destroy them on distruction
    std::vector<MPI_Comm> createdCommList_;

    GMX_DISALLOW_COPY_AND_ASSIGN(BiasSharing);
};

/*! \brief Returns if any bias is sharing within a simulation.
 *
 * \param[in] awhParams  The AWH parameters.
 */
bool haveBiasSharingWithinSimulation(const AwhParams& awhParams);

/*! \brief Checks whether biases are compatible for sharing between simulations, throws when not.
 *
 * Should be called simultaneously on the main rank of every simulation.
 * Note that this only checks for technical compatibility. It is up to
 * the user to check that the sharing physically makes sense.
 * Throws an exception when shared biases are not compatible.
 *
 * \param[in] awhParams     The AWH parameters.
 * \param[in] pointSize     Vector of grid-point sizes for each bias.
 * \param[in] biasSharing   Object for communication for sharing bias data over simulations.
 */
void biasesAreCompatibleForSharingBetweenSimulations(const AwhParams&       awhParams,
                                                     ArrayRef<const size_t> pointSize,
                                                     const BiasSharing&     biasSharing);

} // namespace gmx

#endif /* GMX_AWH_BIASSHARING_H */
