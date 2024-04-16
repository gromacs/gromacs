/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief This file declares a manager of GPU context and streams needed for
 * running workloads on GPUs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_gpu_utils
 */
#ifndef GMX_GPU_UTILS_GPUSTREAMMANAGER_H
#define GMX_GPU_UTILS_GPUSTREAMMANAGER_H

#include <memory>

class DeviceContext;
struct DeviceInformation;
class DeviceStream;

namespace gmx
{

class SimulationWorkload;

/*! \brief Class enum to describe the different logical streams used
 * for GPU work.
 *
 * Whether the actual streams differ is an implementation detail of
 * the manager class.
 */
enum class DeviceStreamType : int
{
    //! Stream primarily for short-ranged local nonbonded work.
    NonBondedLocal,
    //! Stream primarily for short-ranged nonlocal nonbonded work.
    NonBondedNonLocal,
    //! Stream primarily for PME work.
    Pme,
    //! Stream primarily for data exchange between PME and PP ranks.
    PmePpTransfer,
    //! Stream primarily for update and constraints.
    UpdateAndConstraints,
    //! Conventional termination of the enumeration.
    Count
};

/*! \libinternal
 * \brief Device stream and context manager.
 *
 * Given information on which device to use and what to compute with
 * it, objects of this class manage the lifetime of the GPU context
 * and streams.
 *
 * If supported by the GPU API, the available runtime and the
 * indicated device, some streams will be configured at high
 * priority. Otherwise, all streams will share the default priority
 * appropriate to the situation.
 */
class DeviceStreamManager
{
public:
    /*! \brief Constructor.
     *
     * \throws InternalError  If any of the required resources could not be initialized.
     */
    DeviceStreamManager(const DeviceInformation& deviceInfo, SimulationWorkload simulationWork, bool useTiming);
    ~DeviceStreamManager();

    /*! \brief Get the device information object of the associated device.
     *
     * \returns reference to device info.
     */
    const DeviceInformation& deviceInfo() const;

    /*! \brief Returns a handle to the device context. */
    const DeviceContext& context() const;

    /*! \brief Returns a handle to the requested GPU stream.
     *
     * \param[in] streamToGet Which stream to get.
     */
    const DeviceStream& stream(DeviceStreamType streamToGet) const;

    //! \brief Returns a handle to the GPU stream to compute bonded forces in.
    const DeviceStream& bondedStream() const;

    /*! \brief Return whether the requested GPU stream is valid for use.
     *
     * \param[in] streamToCheck Which stream to check.
     *
     * \returns Whether the stream was initialized.
     */
    bool streamIsValid(DeviceStreamType streamToCheck) const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
