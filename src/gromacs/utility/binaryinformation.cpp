/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * \brief Implements functionality for printing information about the
 * currently running binary
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "binaryinformation.h"

#include "config.h"

#if GMX_FFT_FFTW3
// Needed for construction of the FFT library description string
#include <fftw3.h>
#endif

#ifdef HAVE_LIBMKL
#include <mkl.h>
#endif

#if HAVE_EXTRAE
#include <extrae_user_events.h>
#endif

#if GMX_HWLOC
#include <hwloc.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

/* This file is completely threadsafe - keep it that way! */

#include "buildinfo.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "cuda_version_information.h"

namespace
{

using gmx::formatString;

//! \cond Doxygen does not need to care about most of this stuff, and the macro usage is painful to document

int centeringOffset(int width, int length)
{
    return std::max(width - length, 0) / 2;
}

std::string formatCentered(int width, const char *text)
{
    const int offset = centeringOffset(width, std::strlen(text));
    return formatString("%*s%s", offset, "", text);
}

void printCopyright(gmx::TextWriter *writer)
{
    static const char * const Contributors[] = {
        "Emile Apol",
        "Rossen Apostolov",
        "Herman J.C. Berendsen",
        "Par Bjelkmar",
        "Aldert van Buuren",
        "Rudi van Drunen",
        "Anton Feenstra",
        "Gerrit Groenhof",
        "Christoph Junghans",
        "Anca Hamuraru",
        "Vincent Hindriksen",
        "Dimitrios Karkoulis",
        "Peter Kasson",
        "Jiri Kraus",
        "Carsten Kutzner",
        "Per Larsson",
        "Justin A. Lemkul",
        "Viveca Lindahl",
        "Magnus Lundborg",
        "Pieter Meulenhoff",
        "Erik Marklund",
        "Teemu Murtola",
        "Szilard Pall",
        "Sander Pronk",
        "Roland Schulz",
        "Alexey Shvetsov",
        "Michael Shirts",
        "Alfons Sijbers",
        "Peter Tieleman",
        "Teemu Virolainen",
        "Christian Wennberg",
        "Maarten Wolf"
    };
    static const char * const CopyrightText[] = {
        "Copyright (c) 1991-2000, University of Groningen, The Netherlands.",
        "Copyright (c) 2001-2017, The GROMACS development team at",
        "Uppsala University, Stockholm University and",
        "the Royal Institute of Technology, Sweden.",
        "check out http://www.gromacs.org for more information."
    };
    static const char * const LicenseText[] = {
        "GROMACS is free software; you can redistribute it and/or modify it",
        "under the terms of the GNU Lesser General Public License",
        "as published by the Free Software Foundation; either version 2.1",
        "of the License, or (at your option) any later version."
    };

#define NCONTRIBUTORS (int)asize(Contributors)
#define NCR (int)asize(CopyrightText)

// FAH has an exception permission from LGPL to allow digital signatures in Gromacs.
#ifdef GMX_FAHCORE
#define NLICENSE 0
#else
#define NLICENSE (int)asize(LicenseText)
#endif

    // TODO a centering behaviour of TextWriter could be useful here
    writer->writeLine(formatCentered(78, "GROMACS is written by:"));
    for (int i = 0; i < NCONTRIBUTORS; )
    {
        for (int j = 0; j < 4 && i < NCONTRIBUTORS; ++j, ++i)
        {
            const int width = 18;
            char      buf[30];
            const int offset = centeringOffset(width, strlen(Contributors[i]));
            GMX_RELEASE_ASSERT(strlen(Contributors[i]) + offset < asize(buf),
                               "Formatting buffer is not long enough");
            std::fill(buf, buf+width, ' ');
            std::strcpy(buf+offset, Contributors[i]);
            writer->writeString(formatString(" %-*s", width, buf));
        }
        writer->ensureLineBreak();
    }
    writer->writeLine(formatCentered(78, "and the project leaders:"));
    writer->writeLine(formatCentered(78, "Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel"));
    writer->ensureEmptyLine();
    for (int i = 0; i < NCR; ++i)
    {
        writer->writeLine(CopyrightText[i]);
    }
    writer->ensureEmptyLine();
    for (int i = 0; i < NLICENSE; ++i)
    {
        writer->writeLine(LicenseText[i]);
    }
}

// Construct a string that describes the library that provides FFT support to this build
const char *getFftDescriptionString()
{
// Define the FFT description string
#if GMX_FFT_FFTW3
#  if GMX_NATIVE_WINDOWS
    // Don't buy trouble
    return "fftw3";
#  else
    // Use the version string provided by libfftw3
#    if GMX_DOUBLE
    return fftw_version;
#    else
    return fftwf_version;
#    endif
#  endif
#endif
#if GMX_FFT_MKL
    return "Intel MKL";
#endif
#if GMX_FFT_FFTPACK
    return "fftpack (built-in)";
#endif
};

void gmx_print_version_info(gmx::TextWriter *writer)
{
    writer->writeLine(formatString("GROMACS version:    %s", gmx_version()));
    const char *const git_hash = gmx_version_git_full_hash();
    if (git_hash[0] != '\0')
    {
        writer->writeLine(formatString("GIT SHA1 hash:      %s", git_hash));
    }
    const char *const base_hash = gmx_version_git_central_base_hash();
    if (base_hash[0] != '\0')
    {
        writer->writeLine(formatString("Branched from:      %s", base_hash));
    }

#if GMX_DOUBLE
    writer->writeLine("Precision:          double");
#else
    writer->writeLine("Precision:          single");
#endif
    writer->writeLine(formatString("Memory model:       %u bit", (unsigned)(8*sizeof(void *))));

#if GMX_THREAD_MPI
    writer->writeLine("MPI library:        thread_mpi");
#elif GMX_MPI
    writer->writeLine("MPI library:        MPI");
#else
    writer->writeLine("MPI library:        none");
#endif
#if GMX_OPENMP
    writer->writeLine(formatString("OpenMP support:     enabled (GMX_OPENMP_MAX_THREADS = %d)", GMX_OPENMP_MAX_THREADS));
#else
    writer->writeLine("OpenMP support:     disabled");
#endif
    writer->writeLine(formatString("GPU support:        %s", getGpuImplementationString()));
    writer->writeLine(formatString("SIMD instructions:  %s", GMX_SIMD_STRING));
    writer->writeLine(formatString("FFT library:        %s", getFftDescriptionString()));
    writer->writeLine(formatString("RDTSCP usage:       %s", HAVE_RDTSCP ? "enabled" : "disabled"));
#ifdef GMX_USE_TNG
    writer->writeLine("TNG support:        enabled");
#else
    writer->writeLine("TNG support:        disabled");
#endif
#if GMX_HWLOC
    writer->writeLine(formatString("Hwloc support:      hwloc-%d.%d.%d",
                                   HWLOC_API_VERSION>>16,
                                   (HWLOC_API_VERSION>>8) & 0xFF,
                                   HWLOC_API_VERSION & 0xFF));
#else
    writer->writeLine("Hwloc support:      disabled");
#endif
#if HAVE_EXTRAE
    unsigned major, minor, revision;
    Extrae_get_version(&major, &minor, &revision);
    writer->writeLine(formatString("Tracing support:    enabled. Using Extrae-%d.%d.%d", major, minor, revision));
#else
    writer->writeLine("Tracing support:    disabled");
#endif


    writer->writeLine(formatString("Built on:           %s", BUILD_TIME));
    writer->writeLine(formatString("Built by:           %s", BUILD_USER));
    writer->writeLine(formatString("Build OS/arch:      %s", BUILD_HOST));
    writer->writeLine(formatString("Build CPU vendor:   %s", BUILD_CPU_VENDOR));
    writer->writeLine(formatString("Build CPU brand:    %s", BUILD_CPU_BRAND));
    writer->writeLine(formatString("Build CPU family:   %d   Model: %d   Stepping: %d",
                                   BUILD_CPU_FAMILY, BUILD_CPU_MODEL, BUILD_CPU_STEPPING));
    /* TODO: The below strings can be quite long, so it would be nice to wrap
     * them. Can wait for later, as the master branch has ready code to do all
     * that. */
    writer->writeLine(formatString("Build CPU features: %s", BUILD_CPU_FEATURES));
    writer->writeLine(formatString("C compiler:         %s", BUILD_C_COMPILER));
    writer->writeLine(formatString("C compiler flags:   %s", BUILD_CFLAGS));
    writer->writeLine(formatString("C++ compiler:       %s", BUILD_CXX_COMPILER));
    writer->writeLine(formatString("C++ compiler flags: %s", BUILD_CXXFLAGS));
#ifdef HAVE_LIBMKL
    /* MKL might be used for LAPACK/BLAS even if FFTs use FFTW, so keep it separate */
    writer->writeLine(formatString("Linked with Intel MKL version %d.%d.%d.",
                                   __INTEL_MKL__, __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__));
#endif
#if GMX_GPU == GMX_GPU_OPENCL
    writer->writeLine(formatString("OpenCL include dir: %s", OPENCL_INCLUDE_DIR));
    writer->writeLine(formatString("OpenCL library:     %s", OPENCL_LIBRARY));
    writer->writeLine(formatString("OpenCL version:     %s", OPENCL_VERSION_STRING));
#endif
#if GMX_GPU == GMX_GPU_CUDA
    writer->writeLine(formatString("CUDA compiler:      %s\n", CUDA_COMPILER_INFO));
    writer->writeLine(formatString("CUDA compiler flags:%s\n", CUDA_COMPILER_FLAGS));
    auto driverVersion = gmx::getCudaDriverVersion();
    writer->writeLine(formatString("CUDA driver:        %d.%d\n", driverVersion.first, driverVersion.second));
    auto runtimeVersion = gmx::getCudaRuntimeVersion();
    writer->writeLine(formatString("CUDA runtime:       %d.%d\n", runtimeVersion.first, runtimeVersion.second));
#endif
}

//! \endcond

} // namespace

namespace gmx
{

BinaryInformationSettings::BinaryInformationSettings()
    : bExtendedInfo_(false), bCopyright_(false),
      bGeneratedByHeader_(false), prefix_(""), suffix_("")
{
}

void printBinaryInformation(FILE                  *fp,
                            const IProgramContext &programContext)
{
    TextWriter writer(fp);
    printBinaryInformation(&writer, programContext, BinaryInformationSettings());
}

void printBinaryInformation(FILE                            *fp,
                            const IProgramContext           &programContext,
                            const BinaryInformationSettings &settings)
{
    try
    {
        TextWriter writer(fp);
        printBinaryInformation(&writer, programContext, settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void printBinaryInformation(TextWriter                      *writer,
                            const IProgramContext           &programContext,
                            const BinaryInformationSettings &settings)
{
    // TODO Perhaps the writer could be configured with the prefix and
    // suffix strings from the settings?
    const char *prefix          = settings.prefix_;
    const char *suffix          = settings.suffix_;
    const char *precisionString = "";
#if GMX_DOUBLE
    precisionString = " (double precision)";
#endif
    const char *const name = programContext.displayName();
    if (settings.bGeneratedByHeader_)
    {
        writer->writeLine(formatString("%sCreated by:%s", prefix, suffix));
    }
    // TODO: It would be nice to know here whether we are really running a
    // Gromacs binary or some other binary that is calling Gromacs; we
    // could then print "%s is part of GROMACS" or some alternative text.
    std::string title
        = formatString(":-) GROMACS - %s, %s%s (-:", name, gmx_version(), precisionString);
    const int   indent
        = centeringOffset(78 - std::strlen(prefix) - std::strlen(suffix), title.length()) + 1;
    writer->writeLine(formatString("%s%*c%s%s", prefix, indent, ' ', title.c_str(), suffix));
    writer->writeLine(formatString("%s%s", prefix, suffix));
    if (settings.bCopyright_)
    {
        GMX_RELEASE_ASSERT(prefix[0] == '\0' && suffix[0] == '\0',
                           "Prefix/suffix not supported with copyright");
        printCopyright(writer);
        writer->ensureEmptyLine();
        // This line is printed again after the copyright notice to make it
        // appear together with all the other information, so that it is not
        // necessary to read stuff above the copyright notice.
        // The line above the copyright notice puts the copyright notice is
        // context, though.
        writer->writeLine(formatString("%sGROMACS:      %s, version %s%s%s", prefix, name,
                                       gmx_version(), precisionString, suffix));
    }
    const char *const binaryPath = programContext.fullBinaryPath();
    if (!gmx::isNullOrEmpty(binaryPath))
    {
        writer->writeLine(formatString("%sExecutable:   %s%s", prefix, binaryPath, suffix));
    }
    const gmx::InstallationPrefixInfo installPrefix = programContext.installationPrefix();
    if (!gmx::isNullOrEmpty(installPrefix.path))
    {
        writer->writeLine(formatString("%sData prefix:  %s%s%s", prefix, installPrefix.path,
                                       installPrefix.bSourceLayout ? " (source tree)" : "", suffix));
    }
    const std::string workingDir = Path::getWorkingDirectory();
    if (!workingDir.empty())
    {
        writer->writeLine(formatString("%sWorking dir:  %s%s", prefix, workingDir.c_str(), suffix));
    }
    const char *const commandLine = programContext.commandLine();
    if (!gmx::isNullOrEmpty(commandLine))
    {
        writer->writeLine(formatString("%sCommand line:%s\n%s  %s%s",
                                       prefix, suffix, prefix, commandLine, suffix));
    }
    if (settings.bExtendedInfo_)
    {
        GMX_RELEASE_ASSERT(prefix[0] == '\0' && suffix[0] == '\0',
                           "Prefix/suffix not supported with extended info");
        writer->ensureEmptyLine();
        gmx_print_version_info(writer);
    }
}

} // namespace gmx
