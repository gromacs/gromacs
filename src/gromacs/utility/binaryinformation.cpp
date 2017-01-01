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

std::string printCentered(int width, const char *text)
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
        "Copyright (c) 2001-2015, The GROMACS development team at",
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

    writer->writeString(printCentered(78, "GROMACS is written by:"));
    writer->writeString(formatString("\n"));
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
        writer->writeString(formatString("\n"));
    }
    writer->writeString(printCentered(78, "and the project leaders:"));
    writer->writeString(formatString("\n"));
    writer->writeString(printCentered(78, "Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel"));
    writer->writeString(formatString("\n\n"));
    for (int i = 0; i < NCR; ++i)
    {
        writer->writeString(formatString("%s\n", CopyrightText[i]));
    }
    writer->writeString(formatString("\n"));
    for (int i = 0; i < NLICENSE; ++i)
    {
        writer->writeString(formatString("%s\n", LicenseText[i]));
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
    writer->writeString(formatString("GROMACS version:    %s\n", gmx_version()));
    const char *const git_hash = gmx_version_git_full_hash();
    if (git_hash[0] != '\0')
    {
        writer->writeString(formatString("GIT SHA1 hash:      %s\n", git_hash));
    }
    const char *const base_hash = gmx_version_git_central_base_hash();
    if (base_hash[0] != '\0')
    {
        writer->writeString(formatString("Branched from:      %s\n", base_hash));
    }

#if GMX_DOUBLE
    writer->writeString(formatString("Precision:          double\n"));
#else
    writer->writeString(formatString("Precision:          single\n"));
#endif
    writer->writeString(formatString("Memory model:       %u bit\n", (unsigned)(8*sizeof(void *))));

#if GMX_THREAD_MPI
    writer->writeString(formatString("MPI library:        thread_mpi\n"));
#elif GMX_MPI
    writer->writeString(formatString("MPI library:        MPI\n"));
#else
    writer->writeString(formatString("MPI library:        none\n"));
#endif
#if GMX_OPENMP
    writer->writeString(formatString("OpenMP support:     enabled (GMX_OPENMP_MAX_THREADS = %d)\n", GMX_OPENMP_MAX_THREADS));
#else
    writer->writeString(formatString("OpenMP support:     disabled\n"));
#endif
    writer->writeString(formatString("GPU support:        %s\n", getGpuImplementationString()));
    writer->writeString(formatString("SIMD instructions:  %s\n", GMX_SIMD_STRING));
    writer->writeString(formatString("FFT library:        %s\n", getFftDescriptionString()));
#ifdef HAVE_RDTSCP
    writer->writeString(formatString("RDTSCP usage:       enabled\n"));
#else
    writer->writeString(formatString("RDTSCP usage:       disabled\n"));
#endif
#ifdef GMX_USE_TNG
    writer->writeString(formatString("TNG support:        enabled\n"));
#else
    writer->writeString(formatString("TNG support:        disabled\n"));
#endif
#if GMX_HWLOC
    writer->writeString(formatString("Hwloc support:      hwloc-%d.%d.%d\n",
                                     HWLOC_API_VERSION>>16,
                                     (HWLOC_API_VERSION>>8) & 0xFF,
                                     HWLOC_API_VERSION & 0xFF));
#else
    writer->writeString(formatString("Hwloc support:      disabled\n"));
#endif
#if HAVE_EXTRAE
    unsigned major, minor, revision;
    Extrae_get_version(&major, &minor, &revision);
    writer->writeString(formatString("Tracing support:    enabled. Using Extrae-%d.%d.%d\n", major, minor, revision));
#else
    writer->writeString(formatString("Tracing support:    disabled\n"));
#endif


    writer->writeString(formatString("Built on:           %s\n", BUILD_TIME));
    writer->writeString(formatString("Built by:           %s\n", BUILD_USER));
    writer->writeString(formatString("Build OS/arch:      %s\n", BUILD_HOST));
    writer->writeString(formatString("Build CPU vendor:   %s\n", BUILD_CPU_VENDOR));
    writer->writeString(formatString("Build CPU brand:    %s\n", BUILD_CPU_BRAND));
    writer->writeString(formatString("Build CPU family:   %d   Model: %d   Stepping: %d\n",
                                     BUILD_CPU_FAMILY, BUILD_CPU_MODEL, BUILD_CPU_STEPPING));
    /* TODO: The below strings can be quite long, so it would be nice to wrap
     * them. Can wait for later, as the master branch has ready code to do all
     * that. */
    writer->writeString(formatString("Build CPU features: %s\n", BUILD_CPU_FEATURES));
    writer->writeString(formatString("C compiler:         %s\n", BUILD_C_COMPILER));
    writer->writeString(formatString("C compiler flags:   %s\n", BUILD_CFLAGS));
    writer->writeString(formatString("C++ compiler:       %s\n", BUILD_CXX_COMPILER));
    writer->writeString(formatString("C++ compiler flags: %s\n", BUILD_CXXFLAGS));
#ifdef HAVE_LIBMKL
    /* MKL might be used for LAPACK/BLAS even if FFTs use FFTW, so keep it separate */
    writer->writeString(formatString("Linked with Intel MKL version %d.%d.%d.\n",
                                     __INTEL_MKL__, __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__));
#endif
#if GMX_GPU == GMX_GPU_OPENCL
    writer->writeString(formatString("OpenCL include dir: %s\n", OPENCL_INCLUDE_DIR));
    writer->writeString(formatString("OpenCL library:     %s\n", OPENCL_LIBRARY));
    writer->writeString(formatString("OpenCL version:     %s\n", OPENCL_VERSION_STRING));
#endif
#if GMX_GPU == GMX_GPU_CUDA
    writer->writeString(gmx_print_version_info_cuda_gpu());
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
    const char *prefix          = settings.prefix_;
    const char *suffix          = settings.suffix_;
    const char *precisionString = "";
#if GMX_DOUBLE
    precisionString = " (double precision)";
#endif
    const char *const name = programContext.displayName();
    if (settings.bGeneratedByHeader_)
    {
        writer->writeString(formatString("%sCreated by:%s\n", prefix, suffix));
    }
    // TODO: It would be nice to know here whether we are really running a
    // Gromacs binary or some other binary that is calling Gromacs; we
    // could then print "%s is part of GROMACS" or some alternative text.
    std::string title
        = formatString(":-) GROMACS - %s, %s%s (-:", name, gmx_version(), precisionString);
    const int   indent
        = centeringOffset(78 - std::strlen(prefix) - std::strlen(suffix), title.length()) + 1;
    writer->writeString(formatString("%s%*c%s%s\n", prefix, indent, ' ', title.c_str(), suffix));
    writer->writeString(formatString("%s%s\n", prefix, suffix));
    if (settings.bCopyright_)
    {
        GMX_RELEASE_ASSERT(prefix[0] == '\0' && suffix[0] == '\0',
                           "Prefix/suffix not supported with copyright");
        printCopyright(writer);
        writer->writeString(formatString("\n"));
        // This line is printed again after the copyright notice to make it
        // appear together with all the other information, so that it is not
        // necessary to read stuff above the copyright notice.
        // The line above the copyright notice puts the copyright notice is
        // context, though.
        writer->writeString(formatString("%sGROMACS:      %s, version %s%s%s\n", prefix, name,
                                         gmx_version(), precisionString, suffix));
    }
    const char *const binaryPath = programContext.fullBinaryPath();
    if (!gmx::isNullOrEmpty(binaryPath))
    {
        writer->writeString(formatString("%sExecutable:   %s%s\n", prefix, binaryPath, suffix));
    }
    const gmx::InstallationPrefixInfo installPrefix = programContext.installationPrefix();
    if (!gmx::isNullOrEmpty(installPrefix.path))
    {
        writer->writeString(formatString("%sData prefix:  %s%s%s\n", prefix, installPrefix.path,
                                         installPrefix.bSourceLayout ? " (source tree)" : "", suffix));
    }
    const std::string workingDir = Path::getWorkingDirectory();
    if (!workingDir.empty())
    {
        writer->writeString(formatString("%sWorking dir:  %s%s\n", prefix, workingDir.c_str(), suffix));
    }
    const char *const commandLine = programContext.commandLine();
    if (!gmx::isNullOrEmpty(commandLine))
    {
        writer->writeString(formatString("%sCommand line:%s\n%s  %s%s\n",
                                         prefix, suffix, prefix, commandLine, suffix));
    }
    if (settings.bExtendedInfo_)
    {
        GMX_RELEASE_ASSERT(prefix[0] == '\0' && suffix[0] == '\0',
                           "Prefix/suffix not supported with extended info");
        writer->writeString(formatString("\n"));
        gmx_print_version_info(writer);
    }
}

} // namespace gmx
