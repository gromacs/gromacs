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
#include "gmxpre.h"

#include "fmm_mdpoptions.h"

#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

#include "fmm_mdmodule.h"

namespace gmx
{

namespace
{

//! Helper function to make a std::string containing the module name
std::string moduleName()
{
    return std::string(FmmModuleInfo::sc_name);
}

} // namespace

FmmMdpOptions::FmmMdpOptions()
{
    exaFmmOptions_  = ExaFmmOptions{};
    fmSolvrOptions_ = FMSolvrOptions{};
}

void FmmMdpOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    OptionSectionHandle fmmSectionMdp = options->addSection(OptionSection(moduleName().c_str()));
    fmmSectionMdp.addOption(EnumOption<ActiveFmmBackend>(c_fmmActiveOptionName.c_str())
                                    .enumValue(c_activeFmmBackendNames)
                                    .store(&activeFmmBackend_));
    exaFmmOptions_.initMdpOptionsFmm(fmmSectionMdp);
    fmSolvrOptions_.initMdpOptionsFmm(fmmSectionMdp);
}

void FmmMdpOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, FmmModuleInfo::sc_name, c_fmmActiveOptionName);
    exaFmmOptions_.initMdpTransformFmm(rules);
    fmSolvrOptions_.initMdpTransformFmm(rules);
}

void FmmMdpOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, FmmModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(builder, FmmModuleInfo::sc_name, "module", "; Fast Multipole Method");
    addMdpOutputValue<std::string>(builder,
                                   FmmModuleInfo::sc_name,
                                   c_fmmActiveOptionName,
                                   c_activeFmmBackendNames[activeFmmBackend_]);

    if (activeFmmBackend_ == ActiveFmmBackend::FMSolvr)
    {
        fmSolvrOptions_.buildMdpOutputFmm(builder);
    }
    else if (activeFmmBackend_ == ActiveFmmBackend::ExaFmm)
    {
        exaFmmOptions_.buildMdpOutputFmm(builder);
    }
}

ActiveFmmBackend FmmMdpOptions::activeFmmBackend() const
{
    return activeFmmBackend_;
}

const ExaFmmOptions& FmmMdpOptions::exaFmmOptions() const
{
    if (activeFmmBackend() != ActiveFmmBackend::ExaFmm)
    {
        GMX_THROW(gmx::InternalError("ExaFmmOptions requested, but active backend is not ExaFmm."));
    }
    return exaFmmOptions_;
}

const FMSolvrOptions& FmmMdpOptions::fmSolvrOptions() const
{
    if (activeFmmBackend() != ActiveFmmBackend::FMSolvr)
    {
        GMX_THROW(
                gmx::InternalError("FMSolvrOptions requested, but active backend is not FMSolvr."));
    }
    return fmSolvrOptions_;
}

} // namespace gmx
