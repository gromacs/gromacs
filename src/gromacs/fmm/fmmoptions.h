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
#ifndef GMX_FMM_OPTIONS_H
#define GMX_FMM_OPTIONS_H

#include <string>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

namespace gmx
{
class OptionSectionHandle;
class KeyValueTreeObjectBuilder;

/**
 * @brief Enumeration of possible direct interaction providers.
 */
enum class FmmDirectProvider
{
    Gromacs, //!< Use GROMACS direct interactions
    Fmm,     //!< Use direct interactions from the FMM backend
    Count
};

/**
 * @brief Indicates which FMM backend is active based on MDP configuration.
 *
 * Only one FMM backend should be active at a time. This enum reflects the
 * active backend as determined by user-specified MDP options.
 */
enum class ActiveFmmBackend
{
    Inactive, //!< No FMM backend is active.
    ExaFmm,   //!< ExaFMM backend is active.
    FMSolvr,  //!< FMSolvr backend is active.
    Count
};

//! String names corresponding to ActiveFmmBackend enum values.
static const EnumerationArray<ActiveFmmBackend, const char*> c_activeFmmBackendNames = {
    { "inactive", "exafmm", "fmsolvr" }
};

//! MDP option name to enable one of the FMM backends (e.g., ExaFMM).
const std::string c_fmmActiveOptionName = "backend";

//! MDP option name to set the multipole expansion order for ExaFMM
const std::string c_fmmExaFmmOrderOptionName = "exafmm-order";

//! MDP option name to set the multipole expansion order for FMSolvr
const std::string c_fmmFMSolvrOrderOptionName = "fmsolvr-order";

/**
 * @brief Interface for FMM option sets used in MDP handling.
 *
 * This interface ensures that all FMM backend option types (e.g., ExaFMM, FMSolvr)
 * provide consistent implementations of the three MDP-related methods:
 * - initMdpOptionsFmm(): declares MDP options
 * - initMdpTransformFmm(): sets up option transformations
 * - buildMdpOutputFmm(): outputs values to the MDP file
 *
 * These methods follow the same structure as IMdpOptionProvider, but IFmmOptions does not
 * inherit from it because:
 * - Each FMM backend (e.g., ExaFMM, FMSolvr) is part of a larger option group,
 *   not a standalone MDP option provider.
 * - Only FmmMdpOptions is registered as the provider and delegates to the backends.
 *
 * This keeps the logic for each backend modular while centralizing integration.
 *
 */
struct IFmmOptions
{
    virtual void initMdpOptionsFmm(OptionSectionHandle& section)             = 0;
    virtual void initMdpTransformFmm(IKeyValueTreeTransformRules* rules)     = 0;
    virtual void buildMdpOutputFmm(KeyValueTreeObjectBuilder* builder) const = 0;

    virtual ~IFmmOptions() = default;
    GMX_DEFAULT_CONSTRUCTORS(IFmmOptions);
};

struct ExaFmmOptions : IFmmOptions
{
    int               order          = 6;
    int               directRange    = 2;
    FmmDirectProvider directProvider = FmmDirectProvider::Gromacs;

    void initMdpOptionsFmm(OptionSectionHandle& section) override;
    void initMdpTransformFmm(IKeyValueTreeTransformRules* rules) override;
    void buildMdpOutputFmm(KeyValueTreeObjectBuilder* builder) const override;
};

struct FMSolvrOptions : IFmmOptions
{
    int               order              = 8;
    int               directRange        = 1;
    FmmDirectProvider directProvider     = FmmDirectProvider::Fmm;
    bool              dipoleCompensation = true;
    int               treeDepth          = 3;
    bool              sparse             = false;

    void initMdpOptionsFmm(OptionSectionHandle& section) override;
    void initMdpTransformFmm(IKeyValueTreeTransformRules* rules) override;
    void buildMdpOutputFmm(KeyValueTreeObjectBuilder* builder) const override;
};

} // namespace gmx

#endif // GMX_FMM_OPTIONS_H
