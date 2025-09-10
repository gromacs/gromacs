/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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

/*! \libinternal \file
 * \brief Declares the builder class for constructing an FmmForceProvider.
 *
 * \author Muhammad Umair Sadiq <mumairsadiq1@gmail.com>
 */

#ifndef GMX_FMM_FORCEPROVIDER_BUILDER_H
#define GMX_FMM_FORCEPROVIDER_BUILDER_H

#include "fmmforceprovider.h"

namespace gmx
{

struct MDModulesNotifiers;
struct IFmmOptions;

/*! \brief Helper class to configure and construct FMM Force Provider.
 *
 * Initializes FMM inputs (e.g., topology, PBC type, and logger)
 * from simulation setup and constructs a fully configured FmmForceProvider.
 */
class FmmForceProviderBuilder
{
public:
    FmmForceProviderBuilder& subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers);
    FmmForceProviderBuilder& setFmmOptions(const IFmmOptions* fmmOptions);

    std::unique_ptr<FmmForceProvider> build();

private:
    const IFmmOptions* fmmOptions_ = nullptr;
    const gmx_mtop_t*  mtop_       = nullptr;
    const PbcType*     pbcType_    = nullptr;
    const MDLogger*    logger_     = nullptr;
};

} // namespace gmx

#endif // GMX_FMM_FORCEPROVIDER_BUILDER_H
