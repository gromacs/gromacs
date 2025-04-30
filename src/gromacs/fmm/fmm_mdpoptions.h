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
#ifndef GMX_FMM_MDP_OPTIONS_H
#define GMX_FMM_MDP_OPTIONS_H

#include <vector>

#include "gromacs/mdtypes/imdpoptionprovider.h"

#include "fmmoptions.h"

namespace gmx
{

/**
 * @brief MDP option provider that manages all FMM backends.
 *
 * Unlike IFmmOptions, which encapsulates backend-specific logic,
 * this class centralizes MDP declaration, transformation, and output,
 * and delegates to the currently active backend (e.g., ExaFMM, FMSolvr).
 * Designed to support multiple FMM implementations in a unified way.
 */
class FmmMdpOptions : public IMdpOptionProvider
{
public:
    /**
     * @brief Constructs the FMM MDP options provider.
     */
    FmmMdpOptions();

    /**
     * @brief Declares FMM options in the MDP system.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    /**
     * @brief Registers MDP transformation rules for FMM options.
     */
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /**
     * @brief Adds FMM options to the MDP output.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /**
     * @brief Returns the currently selected FMM backend based on MDP options.
     * Returns ActiveFmmBackend::Inactive if no backend is enabled.
     */
    ActiveFmmBackend activeFmmBackend() const;

    /**
     * @brief Returns the ExaFMM options if ExaFMM is the active backend.
     *
     * @return Reference to the ExaFmmOptions instance.
     * @throws gmx::InternalError if the active backend is not ExaFMM.
     */
    const ExaFmmOptions& exaFmmOptions() const;

    /**
     * @brief Returns the FMSolvr options if FMSolvr is the active backend.
     *
     * @return Reference to the FMSolvrOptions instance.
     * @throws gmx::InternalError if the active backend is not FMSolvr.
     */
    const FMSolvrOptions& fmSolvrOptions() const;

private:
    ExaFmmOptions  exaFmmOptions_;  ///< Options specific to the ExaFMM backend
    FMSolvrOptions fmSolvrOptions_; ///< Options specific to the FMSolvr backend
    ActiveFmmBackend activeFmmBackend_ = ActiveFmmBackend::Inactive; ///< Currently selected FMM backend
};

} // namespace gmx

#endif // GMX_FMM_MDP_OPTIONS_H
