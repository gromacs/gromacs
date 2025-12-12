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
 *
 * \brief Declares functions that write a message about the validation state of features
 *
 * \ingroup module_mdrunutility
 * \inlibraryapi
 */
#ifndef GMX_MDRUNUTILITY_PRINT_VALIDATION_H
#define GMX_MDRUNUTILITY_PRINT_VALIDATION_H

struct DeviceInformation;
enum class IntegrationAlgorithm;

namespace gmx
{

class SimulationWorkload;
class MDLogger;

/*! \brief Write a messages for lists of features with status "experimental" and "validation pending" to the log
 *
 * \param[in] mdlog                            The logger
 * \param[in] simulationWorkload               Describes the simulation workload
 * \param[in] haveFillerParticlesInLocalState  Whether direct-halo communication is in use
 * \param[in] useH5mdOutputFile                Whether a trajectory output file uses H5md format.
 * \param[in] useModularSimulator              Whether the modular simulator is in use
 * \param[in] integrationAlgorithm             Whether integration algorithm is in use
 * \param[in] deviceInfo                       Handle to device information object, can be nullptr
 */
void logValidationMessages(const MDLogger&           mdlog,
                           const SimulationWorkload& simulationWorkload,
                           bool                      haveFillerParticlesInLocalState,
                           bool                      useH5mdOutputFile,
                           bool                      useModularSimulator,
                           IntegrationAlgorithm      integrationAlgorithm,
                           const DeviceInformation*  deviceInfo);

} // namespace gmx

#endif
