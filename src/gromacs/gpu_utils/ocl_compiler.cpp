/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/unique_cptr.h"

#include "ocl_caching.h"

namespace gmx
{
namespace ocl
{

/*! \brief True if OpenCL binary caching is enabled.
 *
 *  Currently caching is disabled by default unless the env var override
 *  is used until we resolve concurrency issues. */
static bool useBuildCache = getenv("GMX_OCL_GENCACHE") != nullptr;

/*! \brief Handles writing the OpenCL JIT compilation log to \c fplog.
 *
 * If \c fplog is non-null and either the \c GMX_OCL_DUMP_LOG environment
 * variable is set or the compilation failed, then the OpenCL
 * compilation log is written.
 *
 * \param fplog               Open file pointer to log file
 * \param program             OpenCL program that was compiled
 * \param deviceId            Id of the device for which compilation took place
 * \param kernelFilename      File name containing the kernel
 * \param preprocessorOptions String containing the preprocessor command-line options used for the
 *                            build
 * \param buildFailed         Whether the OpenCL build succeeded
 *
 * \throws std::bad_alloc if out of memory */
static void writeOclBuildLog(FILE*              fplog,
                             cl_program         program,
                             cl_device_id       deviceId,
                             const std::string& kernelFilename,
                             const std::string& preprocessorOptions,
                             bool               buildFailed)
{
    bool writeOutput = ((fplog != nullptr) && (buildFailed || (getenv("GMX_OCL_DUMP_LOG") != nullptr)));

    if (!writeOutput)
    {
        return;
    }

    // Get build log string size
    size_t buildLogSize;
    cl_int cl_error =
            clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, 0, nullptr, &buildLogSize);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not get OpenCL program build log size, error was "
                                + ocl_get_error_string(cl_error)));
    }

    char*             buildLog = nullptr;
    unique_cptr<char> buildLogGuard;
    if (buildLogSize != 0)
    {
        /* Allocate memory to fit the build log,
           it can be very large in case of errors */
        snew(buildLog, buildLogSize);
        buildLogGuard.reset(buildLog);

        /* Get the actual compilation log */
        cl_error = clGetProgramBuildInfo(
                program, deviceId, CL_PROGRAM_BUILD_LOG, buildLogSize, buildLog, nullptr);
        if (cl_error != CL_SUCCESS)
        {
            GMX_THROW(InternalError("Could not get OpenCL program build log, error was "
                                    + ocl_get_error_string(cl_error)));
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
    message += "---------------LOG END----------------\n";
    ;

    fputs(message.c_str(), fplog);
}

/*! \brief Construct compiler options string
 *
 * \param deviceVendor  Device vendor. Used to automatically enable some
 *                      vendor-specific options.
 * \return The string with the compiler options
 */
static std::string selectCompilerOptions(DeviceVendor deviceVendor)
{
    std::string compilerOptions;

    if (getenv("GMX_OCL_NOOPT"))
    {
        compilerOptions += " -cl-opt-disable";
    }

    /* Fastmath improves performance on all supported arch,
     * but is tends to cause problems on Intel (Issue #3898) */
    if ((deviceVendor != DeviceVendor::Intel) && (getenv("GMX_OCL_DISABLE_FASTMATH") == nullptr))
    {
        compilerOptions += " -cl-fast-relaxed-math";

        // Hint to the compiler that it can flush denorms to zero.
        // In CUDA this is triggered by the -use_fast_math flag, equivalent with
        // -cl-fast-relaxed-math, hence the inclusion on the conditional block.
        compilerOptions += " -cl-denorms-are-zero";
    }

    if ((deviceVendor == DeviceVendor::Nvidia) && getenv("GMX_OCL_VERBOSE"))
    {
        compilerOptions += " -cl-nv-verbose";
    }

    if ((deviceVendor == DeviceVendor::Amd) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* To dump OpenCL build intermediate files, caching must be off */
        if (!useBuildCache)
        {
            compilerOptions += " -save-temps";
        }
    }

    if (getenv("GMX_OCL_DEBUG"))
    {
        compilerOptions += " -g";
    }

    if (deviceVendor == DeviceVendor::Amd && getenv("GMX_OCL_FORCE_AMD_WAVEFRONT64"))
    {
        compilerOptions += " -Wf,-mwavefrontsize64";
    }

    return compilerOptions;
}

/*! \brief Get the path to the folder storing an OpenCL source file.
 *
 * By default, this function constructs the full path to the OpenCL from
 * the known location of the binary that is running, so that we handle
 * both in-source and installed builds. The user can override this
 * behavior by defining GMX_OCL_FILE_PATH environment variable.
 *
 * \param[in] sourceRelativePath    Relative path to the kernel or other file in the source tree,
 *                                  from src, e.g. "gromacs/mdlib/nbnxn_ocl" for NB kernels.
 * \return OS-normalized path string to the folder storing OpenCL source file
 *
 * \throws std::bad_alloc    if out of memory.
 *         FileIOError  if GMX_OCL_FILE_PATH does not specify a readable path
 */
static std::filesystem::path getSourceRootPath(const std::string& sourceRelativePath)
{
    std::filesystem::path sourceRootPath;
    /* Use GMX_OCL_FILE_PATH if the user has defined it */
    const char* gmxOclFilePath = getenv("GMX_OCL_FILE_PATH");

    if (gmxOclFilePath == nullptr)
    {
        /* Normal way of getting ocl_root_dir. First get the right
           root path from the path to the binary that is running. */
        InstallationPrefixInfo info = getProgramContext().installationPrefix();
        std::string dataPathSuffix  = (info.sourceLayoutTreeLike_ ? "src" : GMX_INSTALL_OCLDIR);
        sourceRootPath =
                std::filesystem::path(info.path_).append(dataPathSuffix).append(sourceRelativePath);
    }
    else
    {
        if (!std::filesystem::is_directory(gmxOclFilePath))
        {
            GMX_THROW(FileIOError(
                    formatString("GMX_OCL_FILE_PATH must point to the directory where OpenCL"
                                 "kernels are found, but '%s' does not exist",
                                 gmxOclFilePath)));
        }
        sourceRootPath = std::filesystem::path(gmxOclFilePath).append(sourceRelativePath);
    }

    // Make sure we return an OS-correct path format
    return sourceRootPath.make_preferred();
}

size_t getKernelWarpSize(cl_kernel kernel, cl_device_id deviceId)
{
    size_t warpSize = 0;
    cl_int cl_error = clGetKernelWorkGroupInfo(
            kernel, deviceId, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(warpSize), &warpSize, nullptr);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not query OpenCL preferred workgroup size, error was "
                                + ocl_get_error_string(cl_error)));
    }
    if (warpSize == 0)
    {
        GMX_THROW(InternalError(formatString("Invalid OpenCL warp size encountered")));
    }
    return warpSize;
}

size_t getDeviceWarpSize(cl_context context, cl_device_id deviceId)
{
    cl_int      cl_error;
    const char* warpSizeKernel =
            "__kernel void test(__global int* test){test[get_local_id(0)] = 0;}";
    cl_program program = clCreateProgramWithSource(context, 1, &warpSizeKernel, nullptr, &cl_error);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not create OpenCL program to determine warp size, error was "
                                + ocl_get_error_string(cl_error)));
    }

    cl_error = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not build OpenCL program to determine warp size, error was "
                                + ocl_get_error_string(cl_error)));
    }

    cl_kernel kernel = clCreateKernel(program, "test", &cl_error);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not create OpenCL kernel to determine warp size, error was "
                                + ocl_get_error_string(cl_error)));
    }

    size_t warpSize = getKernelWarpSize(kernel, deviceId);

    cl_error = clReleaseKernel(kernel);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not release OpenCL warp-size kernel, error was "
                                + ocl_get_error_string(cl_error)));
    }
    cl_error = clReleaseProgram(program);
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not release OpenCL warp-size program, error was "
                                + ocl_get_error_string(cl_error)));
    }

    return warpSize;
}

/*! \brief Select a compilation-line define for a vendor-specific kernel choice from vendor id
 *
 * \param[in] deviceVendor Vendor id enumerator
 *
 * \return The appropriate compilation-line define
 */
static std::string makeVendorFlavorChoice(DeviceVendor deviceVendor)
{
    switch (deviceVendor)
    {
        case DeviceVendor::Amd: return "-D_AMD_SOURCE_";
        case DeviceVendor::Nvidia: return "-D_NVIDIA_SOURCE_";
        case DeviceVendor::Intel: return "-D_INTEL_SOURCE_";
        case DeviceVendor::Apple: return "-D_APPLE_SOURCE_";
        default: return "";
    }
}

/*! \brief Create include paths for kernel sources.
 *
 * All OpenCL kernel files are expected to be stored in one single folder.
 *
 * \throws std::bad_alloc  if out of memory.
 */
static std::string makeKernelIncludePathOption(const std::string& unescapedKernelRootPath)
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

/*! \brief Replace duplicated spaces with a single one in string
 *
 * Only the first character will be kept for multiple adjacent characters that
 * are both identical and where the first one returns true for isspace().
 *
 * \param str String that will be modified.
 */
static void removeExtraSpaces(std::string* str)
{
    GMX_RELEASE_ASSERT(str != nullptr, "A pointer to an actual string must be provided");
    std::string::iterator newEnd = std::unique(
            str->begin(), str->end(), [=](char a, char b) { return isspace(a) != 0 && (a == b); });
    str->erase(newEnd, str->end());
}

/*! \brief Builds a string with build options for the OpenCL kernels
 *
 * \throws std::bad_alloc  if out of memory. */
static std::string makePreprocessorOptions(const std::filesystem::path& kernelRootPath,
                                           const std::filesystem::path& includeRootPath,
                                           size_t                       warpSize,
                                           DeviceVendor                 deviceVendor,
                                           const std::string&           extraDefines)
{
    std::string preprocessorOptions;

    /* Compose the complete build options */
    preprocessorOptions = formatString("-DWARP_SIZE_TEST=%d", static_cast<int>(warpSize));
    preprocessorOptions += ' ';
    preprocessorOptions += makeVendorFlavorChoice(deviceVendor);
    preprocessorOptions += ' ';
    preprocessorOptions += extraDefines;
    preprocessorOptions += ' ';
    preprocessorOptions += selectCompilerOptions(deviceVendor);
    preprocessorOptions += ' ';
    preprocessorOptions += makeKernelIncludePathOption(kernelRootPath.generic_string());
    preprocessorOptions += ' ';
    preprocessorOptions += makeKernelIncludePathOption(includeRootPath.generic_string());

    // Mac OS (and maybe some other implementations) does not accept double spaces in options
    removeExtraSpaces(&preprocessorOptions);

    return preprocessorOptions;
}

cl_program compileProgram(FILE*              fplog,
                          const std::string& kernelRelativePath,
                          const std::string& kernelBaseFilename,
                          const std::string& extraDefines,
                          cl_context         context,
                          cl_device_id       deviceId,
                          DeviceVendor       deviceVendor)
{
    cl_int cl_error;
    // Let the kernel find include files from its module.
    auto kernelRootPath = getSourceRootPath(kernelRelativePath);
    // Let the kernel find include files from other modules.
    auto rootPath = getSourceRootPath("");

    GMX_RELEASE_ASSERT(fplog != nullptr, "Need a valid log file for building OpenCL programs");

    /* Load OpenCL source files */
    auto kernelFilename = std::filesystem::path(kernelRootPath).append(kernelBaseFilename);

    /* Make the build options */
    std::string preprocessorOptions = makePreprocessorOptions(
            kernelRootPath, rootPath, getDeviceWarpSize(context, deviceId), deviceVendor, extraDefines);

    bool buildCacheWasRead = false;

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
            catch (FileIOError& e)
            {
                // Failing to read from the cache is not a critical error
                formatExceptionMessageToFile(fplog, e);
            }
            fprintf(fplog,
                    "OpenCL binary cache file %s is present, will load kernels.\n",
                    cacheFilename.c_str());
        }
        else
        {
            fprintf(fplog,
                    "No OpenCL binary cache file was present for %s, so will compile kernels "
                    "normally.\n",
                    kernelBaseFilename.c_str());
        }
    }
    if (program == nullptr)
    {
        // Compile OpenCL program from source
        std::string kernelSource = TextReader::readFileToString(kernelFilename.string());
        if (kernelSource.empty())
        {
            GMX_THROW(FileIOError(gmx::formatString("Error loading OpenCL code %s",
                                                    kernelFilename.string().c_str())));
        }
        const char* kernelSourcePtr  = kernelSource.c_str();
        size_t      kernelSourceSize = kernelSource.size();
        /* Create program from source code */
        program = clCreateProgramWithSource(context, 1, &kernelSourcePtr, &kernelSourceSize, &cl_error);
        if (cl_error != CL_SUCCESS)
        {
            GMX_THROW(InternalError("Could not create OpenCL program, error was "
                                    + ocl_get_error_string(cl_error)));
        }
    }

    /* Build the OpenCL program, keeping the status to potentially
       write to the simulation log file. */
    cl_int buildStatus =
            clBuildProgram(program, 0, nullptr, preprocessorOptions.c_str(), nullptr, nullptr);

    /* Write log first, and then throw exception that the user know what is
       the issue even if the build fails. */
    writeOclBuildLog(
            fplog, program, deviceId, kernelFilename.string(), preprocessorOptions, buildStatus != CL_SUCCESS);

    if (buildStatus != CL_SUCCESS)
    {
        GMX_THROW(InternalError("Could not build OpenCL program, error was "
                                + ocl_get_error_string(buildStatus)));
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
            catch (GromacsException& e)
            {
                // Failing to write the cache is not a critical error
                formatExceptionMessageToFile(fplog, e);
            }
        }
    }
    if ((deviceVendor == DeviceVendor::Nvidia) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* If dumping intermediate files has been requested and this is an NVIDIA card
           => write PTX to file */
        char buffer[STRLEN];

        cl_error = clGetDeviceInfo(deviceId, CL_DEVICE_NAME, sizeof(buffer), buffer, nullptr);
        if (cl_error != CL_SUCCESS)
        {
            GMX_THROW(InternalError("Could not get OpenCL device info, error was "
                                    + ocl_get_error_string(cl_error)));
        }
        std::string ptxFilename = buffer;
        ptxFilename += ".ptx";

        try
        {
            writeBinaryToCache(program, ptxFilename);
        }
        catch (GromacsException& e)
        {
            // Failing to write the cache is not a critical error
            formatExceptionMessageToFile(fplog, e);
        }
    }

    return program;
}

} // namespace ocl
} // namespace gmx
