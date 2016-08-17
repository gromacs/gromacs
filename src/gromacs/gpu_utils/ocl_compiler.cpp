/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

#include "ocl_compiler.h"

#include "config.h"

#include <cstdio>

#include <string>
#include <vector>

#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "ocl_caching.h"

namespace gmx
{
namespace ocl
{

/*! \brief True if OpenCL binary caching is enabled.
 *
 *  Currently caching is disabled by default unless the env var override
 *  is used until we resolve concurrency issues. */
static bool useBuildCache = getenv("GMX_OCL_GENCACHE"); // (NULL == getenv("GMX_OCL_NOGENCACHE"));

/*! \brief Handles writing the OpenCL JIT compilation log to \c fplog.
 *
 * If \c fplog is non-null and either the GMX_OCL_DUMP_LOG environment
 * variable is set or the compilation failed, then the OpenCL
 * compilation log is written.
 *
 * \param fplog               Open file pointer to log file
 * \param program             OpenCL program that was compiled
 * \param deviceId            Id of the device for which compilation took place
 * \param kernelFilename      File name containing the kernel
 * \param preprocessorOptions String containing the preprocessor command-line options used for the build
 * \param buildFailed         Whether the OpenCL build succeeded
 *
 * \throws std::bad_alloc if out of memory */
static void
writeOclBuildLog(FILE              *fplog,
                 cl_program         program,
                 cl_device_id       deviceId,
                 const std::string &kernelFilename,
                 const std::string &preprocessorOptions,
                 bool               buildFailed)
{
    bool writeOutput = ((fplog != nullptr) &&
                        (buildFailed || (getenv("GMX_OCL_DUMP_LOG") != nullptr)));

    if (!writeOutput)
    {
        return;
    }

    // Get build log string size
    size_t buildLogSize;
    cl_int cl_error = clGetProgramBuildInfo(program,
                                            deviceId,
                                            CL_PROGRAM_BUILD_LOG,
                                            0,
                                            NULL,
                                            &buildLogSize);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not get OpenCL program build log size, error was " + ocl_get_error_string(cl_error)));
    }

    char             *buildLog = nullptr;
    scoped_cptr<char> buildLogGuard;
    if (buildLogSize != 0)
    {
        /* Allocate memory to fit the build log,
           it can be very large in case of errors */
        snew(buildLog, buildLogSize);
        buildLogGuard.reset(buildLog);

        /* Get the actual compilation log */
        cl_error = clGetProgramBuildInfo(program,
                                         deviceId,
                                         CL_PROGRAM_BUILD_LOG,
                                         buildLogSize,
                                         buildLog,
                                         NULL);
        if (cl_error != CL_SUCCESS)
        {
            GMX_THROW(InternalError("Could not get OpenCL program build log, error was " + ocl_get_error_string(cl_error)));
        }
    }

    std::string message;
    if (buildFailed)
    {
        message += "Compilation of source file " + kernelFilename + " failed!\n";
    }
    else
    {
        message += "Compilation of source file " + kernelFilename + " was successful!\n";
    }
    message += "-- Used build options: " + preprocessorOptions + "\n";
    message += "--------------LOG START---------------\n";
    message += buildLog;
    message += "---------------LOG END----------------\n";;

    fputs(message.c_str(), fplog);
}

/*! \brief Construct compiler options string
 *
 * \param deviceVendorId  Device vendor id. Used to
 *          automatically enable some vendor-specific options
 * \return The string with the compiler options
 */
static std::string
selectCompilerOptions(ocl_vendor_id_t deviceVendorId)
{
    std::string compilerOptions;

    if (getenv("GMX_OCL_NOOPT") )
    {
        compilerOptions += " -cl-opt-disable";
    }

    if (getenv("GMX_OCL_FASTMATH") )
    {
        compilerOptions += " -cl-fast-relaxed-math";
    }

    if ((deviceVendorId == OCL_VENDOR_NVIDIA) && getenv("GMX_OCL_VERBOSE"))
    {
        compilerOptions += " -cl-nv-verbose";
    }

    if ((deviceVendorId == OCL_VENDOR_AMD) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* To dump OpenCL build intermediate files, caching must be off */
        if (!useBuildCache)
        {
            compilerOptions += " -save-temps";
        }
    }

    if ( ( deviceVendorId == OCL_VENDOR_AMD ) && getenv("GMX_OCL_DEBUG"))
    {
        compilerOptions += " -g";
    }

    return compilerOptions;
}

/*! \brief Get the path to the main folder storing OpenCL kernels.
 *
 * By default, this function constructs the full path to the OpenCL from
 * the known location of the binary that is running, so that we handle
 * both in-source and installed builds. The user can override this
 * behavior by defining GMX_OCL_FILE_PATH environment variable.
 *
 * \return OS-normalized path string to the main folder storing OpenCL kernels
 *
 * \throws std::bad_alloc    if out of memory.
 *         FileIOError  if GMX_OCL_FILE_PATH does not specify a readable path
 */
static std::string
getKernelRootPath()
{
    std::string kernelRootPath;
    /* Use GMX_OCL_FILE_PATH if the user has defined it */
    const char *gmxOclFilePath = getenv("GMX_OCL_FILE_PATH");

    if (gmxOclFilePath == nullptr)
    {
        /* Normal way of getting ocl_root_dir. First get the right
           root path from the path to the binary that is running. */
        InstallationPrefixInfo      info           = getProgramContext().installationPrefix();
        std::string                 dataPathSuffix = (info.bSourceLayout ?
                                                      "src/gromacs/mdlib/nbnxn_ocl" :
                                                      OCL_INSTALL_DIR);
        kernelRootPath = Path::join(info.path, dataPathSuffix);
    }
    else
    {
        if (!Directory::exists(gmxOclFilePath))
        {
            GMX_THROW(FileIOError(formatString("GMX_OCL_FILE_PATH must point to the directory where OpenCL"
                                               "kernels are found, but '%s' does not exist", gmxOclFilePath)));
        }
        kernelRootPath = gmxOclFilePath;
    }

    // Make sure we return an OS-correct path format
    return Path::normalize(kernelRootPath);
}

/*!  \brief Get the warp size reported by device
 *
 *  This is platform implementation dependant and seems to only work on the Nvidia and AMD platforms!
 *  Nvidia reports 32, AMD for GPU 64. Ignore the rest
 *
 *  \param  context   Current OpenCL context
 *  \param  deviceId OpenCL device with the context
 *  \return cl_int value of the warp size
 *
 * \throws InternalError if an OpenCL error was encountered
 */
static size_t
getWarpSize(cl_context context, cl_device_id deviceId)
{
    cl_int      cl_error;
    const char *warpSizeKernel = "__kernel void test(__global int* test){test[get_local_id(0)] = 0;}";
    cl_program  program        = clCreateProgramWithSource(context, 1, (const char**)&warpSizeKernel, NULL, &cl_error);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not create OpenCL program to determine warp size, error was " + ocl_get_error_string(cl_error)));
    }

    cl_error = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not build OpenCL program to determine warp size, error was " + ocl_get_error_string(cl_error)));
    }

    cl_kernel kernel = clCreateKernel(program, "test", &cl_error);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not create OpenCL kernel to determine warp size, error was " + ocl_get_error_string(cl_error)));
    }

    size_t warpSize = 0;
    cl_error = clGetKernelWorkGroupInfo(kernel, deviceId, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                        sizeof(warpSize), &warpSize, NULL);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not measure OpenCL warp size, error was " + ocl_get_error_string(cl_error)));
    }
    if (warpSize == 0)
    {
        GMX_THROW(InternalError(formatString("Did not measure a valid OpenCL warp size")));
    }

    cl_error = clReleaseKernel(kernel);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not release OpenCL warp-size kernel, error was " + ocl_get_error_string(cl_error)));
    }
    cl_error = clReleaseProgram(program);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not release OpenCL warp-size program, error was " + ocl_get_error_string(cl_error)));
    }

    return warpSize;
}

/*! \brief Select a compilation-line define for a vendor-specific kernel choice from vendor id
 *
 * \param[in] vendorId Vendor id enumerator
 *
 * \return The appropriate compilation-line define
 */
static const char *
makeVendorFlavorChoice(ocl_vendor_id_t vendorId)
{
    const char *choice;
    switch (vendorId)
    {
        case OCL_VENDOR_AMD:
            choice = "-D_AMD_SOURCE_";
            break;
        case OCL_VENDOR_NVIDIA:
            choice = "-D_NVIDIA_SOURCE_";
            break;
        default:
            choice = "-D_WARPLESS_SOURCE_";
            break;
    }
    return choice;
}

/*! \brief Create include paths for kernel sources.
 *
 * All OpenCL kernel files are expected to be stored in one single folder.
 *
 * \throws std::bad_alloc  if out of memory.
 */
static std::string makeKernelIncludePathOption(const std::string &unescapedKernelRootPath)
{
    std::string includePathOption;

    /* Apple does not seem to accept the quoted include paths other
     * OpenCL implementations are happy with. Since the standard still says
     * it should be quoted, we handle Apple as a special case.
     */
#ifdef __APPLE__
    includePathOption += "-I";

    // Prepend all the spaces with a backslash
    for (std::string::size_type i = 0; i < unescapedKernelRootPath.length(); i++)
    {
        if (unescapedKernelRootPath[i] == ' ')
        {
            includePathOption.push_back('\\');
        }
        includePathOption.push_back(unescapedKernelRootPath[i]);
    }
#else
    includePathOption += "-I\"" + unescapedKernelRootPath + "\"";
#endif

    return includePathOption;
}

/*! \brief Builds a string with build options for the OpenCL kernels
 *
 * \throws std::bad_alloc  if out of memory. */
std::string
makePreprocessorOptions(const std::string   &kernelRootPath,
                        size_t               warpSize,
                        ocl_vendor_id_t      deviceVendorId,
                        const std::string   &extraDefines)
{
    std::string preprocessorOptions;

    /* Compose the complete build options */
    preprocessorOptions  = formatString("-DWARP_SIZE_TEST=%d", static_cast<int>(warpSize));
    preprocessorOptions += ' ';
    preprocessorOptions += makeVendorFlavorChoice(deviceVendorId);
    preprocessorOptions += ' ';
    preprocessorOptions += extraDefines;
    preprocessorOptions += ' ';
    preprocessorOptions += selectCompilerOptions(deviceVendorId);
    preprocessorOptions += ' ';
    preprocessorOptions += makeKernelIncludePathOption(kernelRootPath);

    return preprocessorOptions;
}

cl_program
compileProgram(FILE              *fplog,
               const std::string &kernelBaseFilename,
               const std::string &extraDefines,
               cl_context         context,
               cl_device_id       deviceId,
               ocl_vendor_id_t    deviceVendorId)
{
    cl_int      cl_error;
    std::string kernelRootPath = getKernelRootPath();

    GMX_RELEASE_ASSERT(fplog != nullptr, "Need a valid log file for building OpenCL programs");

    /* Load OpenCL source files */
    std::string kernelFilename = Path::join(kernelRootPath,
                                            kernelBaseFilename);

    /* Make the build options */
    std::string preprocessorOptions = makePreprocessorOptions(kernelRootPath,
                                                              getWarpSize(context, deviceId),
                                                              deviceVendorId,
                                                              extraDefines);

    bool        buildCacheWasRead = false;

    std::string cacheFilename;
    if (useBuildCache)
    {
        cacheFilename = makeBinaryCacheFilename(kernelBaseFilename, deviceId);
    }

    /* Create OpenCL program */
    cl_program program = nullptr;
    if (useBuildCache)
    {
        if (File::exists(cacheFilename, File::returnFalseOnError))
        {
            /* Check if there's a valid cache available */
            try
            {
                program           = makeProgramFromCache(cacheFilename, context, deviceId);
                buildCacheWasRead = true;
            }
            catch (FileIOError &e)
            {
                // Failing to read from the cache is not a critical error
                formatExceptionMessageToFile(fplog, e);
            }
        }
        else
        {
            fprintf(fplog, "No OpenCL binary cache file was present, so will compile kernels normally.\n");
        }
    }
    if (program == nullptr)
    {
        // Compile OpenCL program from source
        std::string kernelSource = TextReader::readFileToString(kernelFilename);
        if (kernelSource.empty())
        {
            GMX_THROW(FileIOError("Error loading OpenCL code " + kernelFilename));
        }
        const char *kernelSourcePtr  = kernelSource.c_str();
        size_t      kernelSourceSize = kernelSource.size();
        /* Create program from source code */
        program = clCreateProgramWithSource(context,
                                            1,
                                            &kernelSourcePtr,
                                            &kernelSourceSize,
                                            &cl_error);
        if (cl_error != CL_SUCCESS)
        {
            GMX_THROW(InternalError("Could not create OpenCL program, error was " + ocl_get_error_string(cl_error)));
        }
    }

    /* Build the OpenCL program, keeping the status to potentially
       write to the simulation log file. */
    cl_int buildStatus = clBuildProgram(program, 0, NULL, preprocessorOptions.c_str(), NULL, NULL);

    /* Write log first, and then throw exception that the user know what is
       the issue even if the build fails. */
    writeOclBuildLog(fplog,
                     program,
                     deviceId,
                     kernelFilename,
                     preprocessorOptions,
                     buildStatus != CL_SUCCESS);

    if (buildStatus != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not build OpenCL program, error was " + ocl_get_error_string(buildStatus)));
    }

    if (useBuildCache)
    {
        if (!buildCacheWasRead)
        {
            /* If OpenCL caching is ON, but the current cache is not
               valid => update it */
            try
            {
                writeBinaryToCache(program, cacheFilename);
            }
            catch (GromacsException &e)
            {
                // Failing to write the cache is not a critical error
                formatExceptionMessageToFile(fplog, e);
            }
        }
    }
    if ((OCL_VENDOR_NVIDIA == deviceVendorId) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* If dumping intermediate files has been requested and this is an NVIDIA card
           => write PTX to file */
        char buffer[STRLEN];

        cl_error = clGetDeviceInfo(deviceId, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
        if (cl_error != CL_SUCCESS)
        {
            GMX_THROW(InternalError("Could not get OpenCL device info, error was " + ocl_get_error_string(cl_error)));
        }
        std::string ptxFilename = buffer;
        ptxFilename += ".ptx";

        try
        {
            writeBinaryToCache(program, ptxFilename);
        }
        catch (GromacsException &e)
        {
            // Failing to write the cache is not a critical error
            formatExceptionMessageToFile(fplog, e);
        }
    }

    return program;
}

} // namespace
} // namespace
