/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
/*! \internal \file
 *  \brief Define infrastructure for OpenCL JIT compilation for Gromacs
 *
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "ocl_caching.h"

#include <assert.h>

#include <cctype>
#include <cstdio>

#include <algorithm>
#include <array>
#include <iterator>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{
namespace ocl
{

std::string makeBinaryCacheFilename(const std::string& kernelFilename, cl_device_id deviceId)
{
    // Note that the OpenCL API is defined in terms of bytes, and we
    // assume that sizeof(char) is one byte.
    std::array<char, 1024> deviceName;
    size_t                 deviceNameLength;
    cl_int                 cl_error = clGetDeviceInfo(deviceId, CL_DEVICE_NAME, deviceName.size(),
                                      deviceName.data(), &deviceNameLength);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError(formatString("Could not get OpenCL device name, error was %s",
                                             ocl_get_error_string(cl_error).c_str())));
    }

    std::string cacheFilename = "OCL-cache";
    /* remove the kernel source suffix */
    cacheFilename += "_" + stripSuffixIfPresent(kernelFilename, ".cl") + "_";
    /* We want a cache filename that's somewhat human readable, and
       describes the device because it's based on the vendor's
       information, but also always works as a filename. So we remove
       characters that are commonly illegal in filenames (dot, slash),
       or sometimes inconvenient (whitespace), or perhaps problematic
       (symbols), by permitting only alphanumeric characters from the
       current locale. We assume these work well enough in a
       filename. */
    std::copy_if(deviceName.begin(), deviceName.begin() + deviceNameLength,
                 std::back_inserter(cacheFilename), isalnum);
    cacheFilename += ".bin";

    return cacheFilename;
}

cl_program makeProgramFromCache(const std::string& filename, cl_context context, cl_device_id deviceId)
{
    // TODO all this file reading stuff should become gmx::BinaryReader
    const auto f = create_unique_with_deleter(fopen(filename.c_str(), "rb"), fclose);
    if (!f)
    {
        GMX_THROW(FileIOError("Failed to open binary cache file " + filename));
    }

    // TODO more stdio error handling
    fseek(f.get(), 0, SEEK_END);
    unsigned char*             binary;
    unique_cptr<unsigned char> binaryGuard;
    size_t                     fileSize = ftell(f.get());
    snew(binary, fileSize);
    binaryGuard.reset(binary);
    fseek(f.get(), 0, SEEK_SET);
    size_t readCount = fread(binary, 1, fileSize, f.get());

    if (readCount != fileSize)
    {
        GMX_THROW(FileIOError("Failed to read binary cache file " + filename));
    }

    /* TODO If/when caching is re-enabled, compare current build
     * options and code against the build options and the code
     * corresponding to the cache. If any change is detected then the
     * cache cannot be used.
     *
     * Also caching functionality will need full re-testing. */

    /* Create program from pre-built binary */
    cl_int     cl_error;
    cl_program program = clCreateProgramWithBinary(context, 1, &deviceId, &fileSize,
                                                   const_cast<const unsigned char**>(&binary),
                                                   nullptr, &cl_error);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not create OpenCL program, error was "
                                + ocl_get_error_string(cl_error)));
    }

    return program;
}

void writeBinaryToCache(cl_program program, const std::string& filename)
{
    size_t fileSize;
    cl_int cl_error =
            clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(fileSize), &fileSize, nullptr);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not get OpenCL program binary size, error was "
                                + ocl_get_error_string(cl_error)));
    }

    // TODO all this file writing stuff should become gmx::BinaryWriter
    unsigned char* binary;
    snew(binary, fileSize);
    const unique_cptr<unsigned char> binaryGuard(binary);

    cl_error = clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(binary), &binary, nullptr);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not get OpenCL program binary, error was "
                                + ocl_get_error_string(cl_error)));
    }

    const auto f = create_unique_with_deleter(fopen(filename.c_str(), "wb"), fclose);
    if (!f)
    {
        GMX_THROW(FileIOError("Failed to open binary cache file " + filename));
    }

    fwrite(binary, 1, fileSize, f.get());
}

} // namespace ocl
} // namespace gmx
