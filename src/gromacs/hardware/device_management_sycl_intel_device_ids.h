/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 *  \brief Lookup Intel hardware version from PCI Express ID.
 *
 *  Extracted into a separate file because it contains a huge data table.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \internal
 * \ingroup module_hardware
 */
#ifndef GMX_HARDWARE_DEVICE_MANAGEMENT_SYCL_INTEL_DEVICE_IDS_H
#define GMX_HARDWARE_DEVICE_MANAGEMENT_SYCL_INTEL_DEVICE_IDS_H

#include <optional>
#include <tuple>

/*! \brief Look up Intel hardware version from device's PCI Express ID.
 *
 * The returned values correspond to the ones \c ocloc uses.
 *
 * \param[in] pciExpressID Device ID reported in the device name.
 * \returns major.minor.revision if device is found in the database, \c std::nullopt otherwise.
 */
std::optional<std::tuple<int, int, int>> getIntelHardwareVersionFromPciExpressID(unsigned int pciExpressID);

#endif // GMX_HARDWARE_DEVICE_MANAGEMENT_SYCL_INTEL_DEVICE_IDS_H
