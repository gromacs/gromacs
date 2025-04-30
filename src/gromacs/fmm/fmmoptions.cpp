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
#include "fmmoptions.h"

#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"

#include "fmm_mdmodule.h"

namespace gmx
{

//! String names corresponding to FmmDirectProvider enum values.
static const EnumerationArray<FmmDirectProvider, const char*> c_fmmDirectProviderNames = {
    { "GROMACS", "FMM" }
};

// ExaFMM MDP option names

//! MDP option name to configure the direct interaction range for ExaFMM (1 or 2)
const std::string c_fmmExaFmmDirectRangeOptionName = "exafmm-direct-range";
//! MDP option name to select the direct interaction provider for ExaFMM (GROMACS or FMM)
const std::string c_fmmExaFmmDirectProviderOptionName = "exafmm-direct-provider";
//! MDP option name to set the maximum number of particles per cell in the ExaFMM tree
const std::string c_fmmExaFmmMaxParticlesPerCellOptionName = "exafmm-max-particles-per-cell";


// FMSolvr MDP option names

//! MDP option name to configure the direct interaction range for FMSolvr (1 or 2)
const std::string c_fmmFMSolvrDirectRangeOptionName = "fmsolvr-direct-range";
//! MDP option name to select the direct interaction provider for FMSolvr (GROMACS or FMM).
const std::string c_fmmFMSolvrDirectProviderOptionName = "fmsolvr-direct-provider";
//! MDP option name to enable dipole compensation for improved accuracy in FMSolvr
const std::string c_fmmFMSolvrDipoleCompensationOptionName = "fmsolvr-dipole-compensation";
//! MDP option name to set the tree depth for FMSolvr (controls spatial subdivision granularity)
const std::string c_fmmFMSolvrTreeDepthOptionName = "fmsolvr-tree-depth";
//! MDP option name to enable sparse tree representation in FMSolvr for memory optimization
const std::string c_fmmFMSolvrSparseOptionName = "fmsolvr-sparse";


void ExaFmmOptions::initMdpOptionsFmm(OptionSectionHandle& section)
{
    section.addOption(IntegerOption(c_fmmExaFmmOrderOptionName.c_str()).store(&order));
    section.addOption(IntegerOption(c_fmmExaFmmDirectRangeOptionName.c_str()).store(&directRange));
    section.addOption(EnumOption<FmmDirectProvider>(c_fmmExaFmmDirectProviderOptionName.c_str())
                              .enumValue(c_fmmDirectProviderNames)
                              .store(&directProvider));
}

void ExaFmmOptions::initMdpTransformFmm(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<int>(
            rules, &fromStdString<int>, FmmModuleInfo::sc_name, c_fmmExaFmmOrderOptionName);
    addMdpTransformFromString<int>(
            rules, &fromStdString<int>, FmmModuleInfo::sc_name, c_fmmExaFmmDirectRangeOptionName);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, FmmModuleInfo::sc_name, c_fmmExaFmmDirectProviderOptionName);
    addMdpTransformFromString<int>(
            rules, &fromStdString<int>, FmmModuleInfo::sc_name, c_fmmExaFmmMaxParticlesPerCellOptionName);
}

void ExaFmmOptions::buildMdpOutputFmm(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputValue<int>(builder, FmmModuleInfo::sc_name, c_fmmExaFmmOrderOptionName, order);
    addMdpOutputValue<int>(builder, FmmModuleInfo::sc_name, c_fmmExaFmmDirectRangeOptionName, directRange);
    addMdpOutputValue<std::string>(builder,
                                   FmmModuleInfo::sc_name,
                                   c_fmmExaFmmDirectProviderOptionName,
                                   c_fmmDirectProviderNames[directProvider]);
}

void FMSolvrOptions::initMdpOptionsFmm(OptionSectionHandle& section)
{
    section.addOption(IntegerOption(c_fmmFMSolvrOrderOptionName.c_str()).store(&order));
    section.addOption(IntegerOption(c_fmmFMSolvrDirectRangeOptionName.c_str()).store(&directRange));
    section.addOption(EnumOption<FmmDirectProvider>(c_fmmFMSolvrDirectProviderOptionName.c_str())
                              .enumValue(c_fmmDirectProviderNames)
                              .store(&directProvider));
    section.addOption(
            BooleanOption(c_fmmFMSolvrDipoleCompensationOptionName.c_str()).store(&dipoleCompensation));
    section.addOption(IntegerOption(c_fmmFMSolvrTreeDepthOptionName.c_str()).store(&treeDepth));
    section.addOption(BooleanOption(c_fmmFMSolvrSparseOptionName.c_str()).store(&sparse));
}

void FMSolvrOptions::initMdpTransformFmm(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<int>(
            rules, &fromStdString<int>, FmmModuleInfo::sc_name, c_fmmFMSolvrOrderOptionName);
    addMdpTransformFromString<int>(
            rules, &fromStdString<int>, FmmModuleInfo::sc_name, c_fmmFMSolvrDirectRangeOptionName);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, FmmModuleInfo::sc_name, c_fmmFMSolvrDirectProviderOptionName);
    addMdpTransformFromString<bool>(
            rules, &fromStdString<bool>, FmmModuleInfo::sc_name, c_fmmFMSolvrDipoleCompensationOptionName);
    addMdpTransformFromString<int>(
            rules, &fromStdString<int>, FmmModuleInfo::sc_name, c_fmmFMSolvrTreeDepthOptionName);
    addMdpTransformFromString<bool>(
            rules, &fromStdString<bool>, FmmModuleInfo::sc_name, c_fmmFMSolvrSparseOptionName);
}

void FMSolvrOptions::buildMdpOutputFmm(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputValue<int>(builder, FmmModuleInfo::sc_name, c_fmmFMSolvrOrderOptionName, order);
    addMdpOutputValue<int>(builder, FmmModuleInfo::sc_name, c_fmmFMSolvrDirectRangeOptionName, directRange);
    addMdpOutputValue<std::string>(builder,
                                   FmmModuleInfo::sc_name,
                                   c_fmmFMSolvrDirectProviderOptionName,
                                   c_fmmDirectProviderNames[directProvider]);
    addMdpOutputValue<bool>(
            builder, FmmModuleInfo::sc_name, c_fmmFMSolvrDipoleCompensationOptionName, dipoleCompensation);
    addMdpOutputValue<int>(builder, FmmModuleInfo::sc_name, c_fmmFMSolvrTreeDepthOptionName, treeDepth);
    addMdpOutputValue<bool>(builder, FmmModuleInfo::sc_name, c_fmmFMSolvrSparseOptionName, sparse);
}

} // namespace gmx
