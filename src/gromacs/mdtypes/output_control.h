/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief
 * Declares module for output control parameters
 *
 * It refers to an OutputControl object that is embedded in t_inputrec
 * and populated by the OutputControlModule during MDP parsing. Access
 * it through inputrec:
 *
 * Usage example for MDP parsing:
 * \code
 * t_inputrec ir;
 * gmx::MDModules mdModules;
 * gmx::KeyValueTreeTransformer transform;
 * mdModules.initMdpTransform(transform.rules());
 * // ... apply transform to populate ir.params ...
 * mdModules.adjustInputrecBasedOnModules(&ir);  // Configures modules and normalizes params
 * mdModules.assignOptionsToModules(*ir.params, &errorHandler, &ir);  // Assigns values
 *
 * // Access OutputControl from inputrec
 * const gmx::OutputControl& outputControl = ir.outputControl;
 * \endcode
 *
 * Usage example for reading from TPR:
 * \code
 * t_inputrec ir;
 * gmx::MDModules mdModules;
 * // ... read ir from TPR file ...
 * mdModules.adjustInputrecBasedOnModules(&ir);  // Optional: normalize params
 * mdModules.assignOptionsToModules(*ir.params, nullptr, &ir);  // Configures and assigns
 * \endcode
 *
 * \ingroup module_mdtypes
 * \inlibraryapi
 */
#ifndef GMX_MDTYPES_OUTPUTCONTROL_H
#define GMX_MDTYPES_OUTPUTCONTROL_H

#include <memory>
#include <string_view>

//! Forward declaration for preprocessing-only string storage (local to grompp)
struct gmx_inputrec_strings;

namespace gmx
{

struct OutputControl;
class IMDModule;
class KeyValueTreeObject;

/*! \internal
    \brief Information about the output-control module.
 *
 * Provides name and method to create an output-control module.
 */
struct OutputControlModuleInfo
{
    /*! \brief
     * Creates a module for output control parameters.
     *
     * The returned class manages output frequency parameters that control
     * how often coordinates, velocities, forces, energies, and compressed
     * trajectories are written during simulation.
     */
    static std::unique_ptr<IMDModule> create();
    //! The name of the module
    static constexpr std::string_view sc_name = "output-control";
};

/*! \brief Set the target OutputControl for an OutputControlModule
 *
 * Configures the OutputControlModule to write to an external OutputControl
 * (typically embedded in t_inputrec). This must be called before MDP options
 * are parsed.
 *
 * \param[in] module The module (must be an OutputControlModule)
 * \param[in] outputControl Pointer to the OutputControl to write to
 * \throws InternalError if module is not an OutputControlModule
 */
void setOutputControlTarget(IMDModule* module, OutputControl* outputControl);

/*! \brief Set preprocessing-only string storage for an OutputControlModule
 *
 * Configures where the OutputControlModule stores group specifications
 * (compressed-x-grps, energygrps) that are only used during grompp preprocessing.
 * These strings are NOT part of OutputControl and are NOT serialized to TPR.
 *
 * This must be called before assignOptionsToModules() when parsing MDP files in grompp.
 * It is not needed when reading TPR files in mdrun/tools.
 *
 * \param[in] module              The module (must be an OutputControlModule)
 * \param[in] preprocessingStrings Pointer to storage for preprocessing-only strings
 * \throws InternalError if module is not an OutputControlModule
 */
void setOutputControlPreprocessingStrings(IMDModule* module, gmx_inputrec_strings* preprocessingStrings);

/*! \brief Enable MDP output writing for OutputControl module (for testing only)
 *
 * By default, OutputControlModule::buildMdpOutput() does nothing to avoid duplication
 * with the inp vector mechanism. This function enables output writing so tests can
 * verify the module produces correct MDP output.
 *
 * \param[in] module  The OutputControl module (must be an OutputControlModule)
 * \param[in] enable  Whether to enable MDP output writing
 * \throws InternalError if module is not an OutputControlModule
 */
void setOutputControlWriteToMdpOutput(IMDModule* module, bool enable);

} // namespace gmx

#endif
