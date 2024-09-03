/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief I/o interface to H5MD HDF5 files.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#include "gmxpre.h"

#include "h5md.h"

#include "config.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#if GMX_USE_HDF5
#    include <hdf5.h>
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")
#else
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")
#endif // GMX_USE_HDF5

namespace gmx
{
H5md::H5md(const std::filesystem::path& fileName, const H5mdFileMode mode)
{
#if GMX_USE_HDF5
    /* Disable automatic HDF5 error output, e.g. when items are not found. Explicit H5EPrint2() will
     * still print error messages. */
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    switch (mode)
    {
        case H5mdFileMode::Write:
            file_ = H5Fcreate(
                    fileName.string().c_str(), H5F_ACC_TRUNC, H5Pcreate(H5P_FILE_CREATE), H5P_DEFAULT);
            break;
        case H5mdFileMode::Read:
            file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
        default: throw NotImplementedError("Appending to H5MD is not implemented yet.");
    }
    if (file_ == H5I_INVALID_HID)
    {
#    ifdef NDEBUG
        if (debug)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
        }
#    else
        H5Eprint2(H5E_DEFAULT, nullptr);
#    endif
        throw FileIOError("Cannot open H5MD file.");
    }
    filemode_ = mode;

#else
    GMX_UNUSED_VALUE(fileName);
    GMX_UNUSED_VALUE(mode);
    throw FileIOError("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

H5md::~H5md()
{
#if GMX_USE_HDF5
    if (file_ != H5I_INVALID_HID)
    {
        H5Fclose(file_);
    }

    /* Do not throw, if GMX_USE_HDF5 is false, in the destructor. */

#endif
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
