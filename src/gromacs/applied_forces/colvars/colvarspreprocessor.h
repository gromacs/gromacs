/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Declares the Colvars GROMACS proxy class during pre-processing.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_COLVARSPREPROCESSOR_H
#define GMX_APPLIED_FORCES_COLVARSPREPROCESSOR_H


#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/real.h"

#include "colvarproxygromacs.h"

enum class PbcType : int;


namespace gmx
{
class MDLogger;

/*! \internal \brief
 * Class that read a colvars configuration file during pre-processing and
 * retrieve the colvars atoms coordinates to be stored in tpr KVT.
 */
class ColvarsPreProcessor : public ColvarProxyGromacs
{
public:
    /*! \brief Construct ColvarsPreProcessor from its parameters
     *

     * \param[in] colvarsConfigString Content of the colvars input file.
     * \param[in] atoms Atoms topology
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] logger GROMACS logger instance
     * \param[in] ensembleTemperature the constant ensemble temperature
     * \param[in] seed the colvars seed for random number generator
     * \param[in] box Matrix with full box of the system
     * \param[in] x Coordinates of each atom in the system
     */
    ColvarsPreProcessor(const std::string&   colvarsConfigString,
                        t_atoms              atoms,
                        PbcType              pbcType,
                        const MDLogger*      logger,
                        real                 ensembleTemperature,
                        int                  seed,
                        const matrix         box,
                        ArrayRef<const RVec> x);


    //! Return a vector of the colvars atoms coordinates
    std::vector<RVec> getColvarsCoords();

    //! Save all input files of colvars (outside the config file) in the tpr file through the key-value-tree
    bool inputStreamsToKVT(KeyValueTreeObjectBuilder treeBuilder, const std::string& tag);

private:
    //! Atoms coordinates of the whole system
    ArrayRef<const RVec> x_;
};


} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARSPREPROCESSOR_H
