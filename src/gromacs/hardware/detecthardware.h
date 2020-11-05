/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_HARDWARE_DETECTHARDWARE_H
#define GMX_HARDWARE_DETECTHARDWARE_H

#include <memory>

struct gmx_hw_info_t;

namespace gmx
{
class HardwareTopology;
class MDLogger;
class PhysicalNodeCommunicator;

/*! \brief Run detection and make correct and consistent
 * hardware information available on all ranks.
 *
 * May do communication on MPI_COMM_WORLD when compiled with real MPI.
 *
 * This routine is designed to be called once on each process.  In a
 * thread-MPI configuration, it may only be called before the threads
 * are spawned. With real MPI, communication is needed to coordinate
 * the results. In all cases, any thread within a process may use the
 * returned handle.
 *
 * \todo Replace the use of MPI_COMM_WORLD e.g. by using a libraryCommWorld
 * argument. See https://gitlab.com/gromacs/gromacs/-/issues/3650
 */
std::unique_ptr<gmx_hw_info_t> gmx_detect_hardware(const PhysicalNodeCommunicator& physicalNodeComm);

/*! \brief Sanity check hardware topology and print some notes to log
 *
 *  \param mdlog            Logger.
 *  \param hardwareTopology Reference to hardwareTopology object.
 */
void hardwareTopologyDoubleCheckDetection(const gmx::MDLogger&         mdlog,
                                          const gmx::HardwareTopology& hardwareTopology);

/*! \brief Issue warnings to mdlog that were decided during detection
 *
 * \param[in] mdlog                Logger
 * \param[in] hardwareInformation  The hardwareInformation */
void logHardwareDetectionWarnings(const gmx::MDLogger& mdlog, const gmx_hw_info_t& hardwareInformation);

} // namespace gmx

#endif
