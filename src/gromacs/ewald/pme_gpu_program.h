/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

/*! \libinternal \file
 * \brief
 * Declares PmeGpuProgram, which wrap arounds PmeGpuProgramImpl
 * to store permanent PME GPU context-derived data,
 * such as (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 * \inlibraryapi
 */

#ifndef GMX_EWALD_PME_PME_GPU_PROGRAM_H
#define GMX_EWALD_PME_PME_GPU_PROGRAM_H

#include <memory>

struct PmeGpuProgramImpl;
struct gmx_device_info_t;

class PmeGpuProgram
{
public:
    explicit PmeGpuProgram(const gmx_device_info_t* deviceInfo);
    ~PmeGpuProgram();

    // TODO: design getters for information inside, if needed for PME, and make this private?
    std::unique_ptr<PmeGpuProgramImpl> impl_;
};

/*! \brief This is an owning handle for the compiled PME GPU kernels.
 */
using PmeGpuProgramStorage = std::unique_ptr<PmeGpuProgram>;

/*! \brief This is a handle for passing references to PME GPU program data.
 * TODO: it should be a const reference, but for that the PmeGpu types need to be C++
 */
using PmeGpuProgramHandle = const PmeGpuProgram*;

/*! \brief
 * Factory function used to build persistent PME GPU program for the device at once.
 */
PmeGpuProgramStorage buildPmeGpuProgram(const gmx_device_info_t* /*deviceInfo*/);

#endif
