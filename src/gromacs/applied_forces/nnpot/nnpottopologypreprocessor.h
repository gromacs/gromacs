/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * NNPotTopologyPrepocessor class responsible for
 * all modificatios of the topology during input pre-processing.
 * Inherits from QMMMTopologyPreprocessor.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GROMACS_APPLIED_FORCES_NNPOT_NNPOTTOPOLGYPREPROCESSOR_H
#define GROMACS_APPLIED_FORCES_NNPOT_NNPOTTOPOLGYPREPROCESSOR_H

#include "gromacs/applied_forces/qmmm/qmmmtopologypreprocessor.h"
#include "gromacs/utility/arrayref.h"

struct gmx_mtop_t;

namespace gmx
{

/*! \brief
 * Class responsible for all modifications of the topology during input pre-processing.
 *
 * Inherits from QMMMTopologyPreprocessor and performs the following modifications
 * for a scheme analoguous to mechanical embedding in QMMM:
 * 1) Split molecules containing NNP input atoms from other molecules in blocks
 * 2) Exclude non-bonded interactions between NNP atoms.
 * 3) Build vector with atomic numbers of all atoms.
 * 4) Make F_CONNBOND between atoms within NNP region.
 * 5) Remove angles and settles containing 2 or more NNP atoms.
 * 6) Remove dihedrals containing 3 or more NNP atoms.
 */
class NNPotTopologyPreprocessor : public QMMMTopologyPreprocessor
{
public:
    /*! \brief Constructor for NNPotTopologyPreprocessor from its parameters
     *
     * \param[in] inputIndices Array with global indices of NN input atoms
     */
    NNPotTopologyPreprocessor(ArrayRef<const Index> inputIndices);

    /*! \brief Process mtop topology and builds topInfo_
     * containing information about topology modifications.
     *
     * \param[in,out] mtop Topology that needs to be modified
     */
    void preprocess(gmx_mtop_t* mtop);
};

} // namespace gmx

#endif
