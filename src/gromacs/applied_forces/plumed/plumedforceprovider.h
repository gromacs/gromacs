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
 * Declares Plumed force provider class
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_PLUMEDFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_PLUMEDFORCEPROVIDER_H

#include <memory>

#include "gromacs/mdtypes/iforceprovider.h"

namespace PLMD
{
class Plumed;
}
namespace gmx
{
struct PlumedOptions;
/*! \internal \brief
 * Implements IForceProvider for PLUMED.
 */
class PlumedForceProvider final : public IForceProvider
{
public:
    /*! \brief Initialize the PLUMED interface with the given options
     *
     * \param options PLUMED options
     */
    PlumedForceProvider(const PlumedOptions& options);
    ~PlumedForceProvider();
    /*! \brief Tells PLUMED to output the checkpoint data
     *
     * If the PLUMED API version is not greater than 3 it will do nothing.
     */
    void writeCheckpointData();
    /*! \brief Calculate the forces with PLUMED
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

private:
    std::unique_ptr<PLMD::Plumed> plumed_;
    int                           plumedAPIversion_;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_PLUMEDFORCEPROVIDER_H
