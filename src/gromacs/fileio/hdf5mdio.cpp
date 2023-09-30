/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
/* This file was inspired by ch5md by Pierre de Buyl (BSD license). */

#include "gmxpre.h"

#include "hdf5mdio.h"

#include "config.h"

#include <string>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

// FIXME: TEMPORARY FOR EASIER EDITIING:
#define GMX_USE_HDF5 1

#if GMX_USE_HDF5
#    include <h5xx/h5xx.hpp>
#endif

static unsigned convertModeStringToInt(std::string modeString)
{
#if GMX_USE_HDF5
    const gmx_unused int maxStringLength = 3;
    GMX_ASSERT(modeString.length() <= maxStringLength, "The mode string is too long");
    std::string modeStringLower = gmx::toLowerCase(modeString);
    unsigned mode = h5xx::file::in; // Reading is always enabled.

    std::size_t found = modeStringLower.find("w");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::out;
    }
    found = modeStringLower.find("t");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::trunc;
    }
    found = modeStringLower.find("e");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::excl;
    }

    return mode;
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::GmxHdf5MdIo()
{
#ifdef GMX_USE_HDF5
    file_ = new h5xx::file();
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::GmxHdf5MdIo(const std::string &fileName, const std::string &modeString)
{
#ifdef GMX_USE_HDF5
    openFile(fileName, modeString);
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::~GmxHdf5MdIo()
{
#ifdef GMX_USE_HDF5

#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::openFile(const std::string &fileName, const std::string &modeString)
{
#ifdef GMX_USE_HDF5
    unsigned mode = convertModeStringToInt(modeString);
    file_->open(fileName, mode);
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::closeFile()
{
#ifdef GMX_USE_HDF5
    file_->flush();
    file_->close();
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}
