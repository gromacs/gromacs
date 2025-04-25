/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * \brief Implements functionality for printing information about the
 * currently running binary
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/mdrun/binary_information.h"

#include "config.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

/* This file is completely threadsafe - keep it that way! */

#include "gromacs/applied_forces/colvars/colvarsMDModule.h"
#include "gromacs/applied_forces/nnpot/binary_information.h"
#include "gromacs/applied_forces/plumed/plumedMDModule.h"
#include "gromacs/applied_forces/qmmm/qmmmforceprovider.h"
#include "gromacs/fft/binary_information.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/gpu_utils/binary_information.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/simd_support.h"
#include "gromacs/linearalgebra/binary_information.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"
#include "gromacs/utility/textwriter.h"

#include "buildinfo.h"
#include "contributors.h"

namespace
{

using gmx::formatString;

int centeringOffset(int width, int length)
{
    return std::max(width - length, 0) / 2;
}

std::string formatCentered(int width, const char* text)
{
    const int offset = centeringOffset(width, std::strlen(text));
    return formatString("%*s%s", offset, "", text);
}

void writeVectorAsColumns(gmx::TextWriter*                writer,
                          const std::string&              header,
                          const std::vector<std::string>& v,
                          std::size_t                     outputWidth = 80)
{
    if (v.empty())
    {
        return;
    }

    writer->writeLine(formatCentered(outputWidth, header.c_str()));

    const std::size_t maxWidth =
            std::accumulate(v.begin(),
                            v.end(),
                            std::size_t{ 0 },
                            [](const auto a, const auto& s) { return std::max(a, s.length()); });

    const int columns     = outputWidth / (maxWidth + 1);
    const int columnWidth = outputWidth / columns;

    for (std::size_t i = 0; i < v.size(); i++)
    {
        const std::size_t padLeft = (columnWidth - v[i].length()) / 2;
        const std::string paddedString{ std::string(padLeft, ' ') + v[i] };

        // Using maxWidth+1 when calculating #columns above means we will always have spaces
        writer->writeString(formatString("%-*s", columnWidth, paddedString.c_str()));
        if ((i + 1) % columns == 0)
        {
            writer->ensureLineBreak();
        }
    }
    writer->ensureEmptyLine();
}

void writeVectorAsSingleLine(gmx::TextWriter*                writer,
                             const std::string&              header,
                             const std::vector<std::string>& v,
                             std::size_t                     outputWidth = 80)
{
    if (v.empty())
    {
        return;
    }

    writer->writeLine(formatCentered(outputWidth, header.c_str()));

    std::string s;
    for (std::size_t i = 0; i < v.size(); i++)
    {
        s += v[i];
        if (i < v.size() - 2)
        {
            s += ", ";
        }
        else if (i == v.size() - 2)
        {
            s += (v.size() == 2) ? " and " : ", and "; // Oxford comma for >2 names...
        }
    }
    writer->writeLine(formatCentered(outputWidth, s.c_str()));
    writer->ensureEmptyLine();
}

void printCopyright(gmx::TextWriter* writer)
{
    writer->writeLine(gmx::copyrightText);
    if (!GMX_FAHCORE)
    {
        // Folding At Home has different licence to allow digital
        // signatures in GROMACS, so does not need to show the normal
        // license statement.
        writer->writeLine("GROMACS is free software; you can redistribute it and/or modify it");
        writer->writeLine("under the terms of the GNU Lesser General Public License");
        writer->writeLine("as published by the Free Software Foundation; either version 2.1");
        writer->writeLine("of the License, or (at your option) any later version.");
    }
    writer->ensureEmptyLine();

    writeVectorAsColumns(writer, "Current GROMACS contributors:", gmx::currentContributors);
    writeVectorAsColumns(writer, "Previous GROMACS contributors:", gmx::previousContributors);
    writeVectorAsSingleLine(
            writer, "Coordinated by the GROMACS project leaders:", gmx::currentProjectLeaders);
}

} // namespace

namespace gmx
{

namespace
{

//! Helper function for writing binary information descriptions
void writeLine(TextWriter* writer, const std::string& label, const std::string& description)
{
    static const int sc_maxWidthOfLabelColumn = 20;
    if (label.length() > sc_maxWidthOfLabelColumn)
    {
        GMX_THROW(APIError(formatString(
                "Label '%s' has too many characters (max %d)", label.c_str(), sc_maxWidthOfLabelColumn)));
    }
    const std::string labelPlusColon = label + ':';
    // Produce lines like "SomeLabel:     Description"
    // so that the description column is left-aligned
    writer->writeLineFormatted(
            "%-*s%s", sc_maxWidthOfLabelColumn + 1, labelPlusColon.c_str(), description.c_str());
}

//! Helper function for writing optional binary information from a map of labelled descriptions
void writeOptionalLineFromMap(TextWriter*                                         writer,
                              const std::string&                                  label,
                              const std::unordered_map<std::string, std::string>& descriptionMap)
{
    if (const auto descriptionIt = descriptionMap.find(label); descriptionIt != descriptionMap.end())
    {
        writeLine(writer, label, descriptionIt->second);
    }
}

void writeExtendedInfo(gmx::TextWriter* writer)
{
    {
        std::unordered_map<std::string, std::string> descriptions = versionDescriptions();
        writeLine(writer, "GROMACS version", descriptions.at("GROMACS version"));
        writeOptionalLineFromMap(writer, "GIT SHA1 hash", descriptions);
        writeOptionalLineFromMap(writer, "Branched from", descriptions);
    }
    writeLine(writer, "Precision", GMX_DOUBLE ? "double" : "mixed");
    writeLine(writer,
              "Memory model",
              gmx::formatString("%u bit", static_cast<unsigned>(CHAR_BIT * sizeof(void*))));
    {
        std::unordered_map<std::string, std::string> descriptions = mpiDescriptions();
        writeLine(writer, "MPI library", descriptions.at("MPI library"));
        writeLine(writer, "MPI version", descriptions.at("MPI version"));
    }
    writeLine(writer, "OpenMP support", openmpDescription());
    std::unordered_map<std::string, std::string> gpuDescriptionStrings = gpuDescriptions();
    writeLine(writer, "GPU support", gpuDescriptionStrings.at("GPU support"));
    if (const std::optional<std::string> description = nbnxmGpuClusteringDescription())
    {
        writeLine(writer, "NBNxM GPU setup", description.value());
    }
    writeLine(writer, "SIMD instructions", simdDescription());
    writeLine(writer, "CPU FFT library", cpuFftDescription());
    writeLine(writer, "GPU FFT library", gpuFftDescription());
    writeLine(writer, "Multi-GPU FFT", multiGpuFftDescription());
    if (const std::optional<std::string> description = rdtscpDescription())
    {
        writeLine(writer, "RDTSCP", description.value());
    }
    writeLine(writer, "TNG support", tngDescription());
    writeLine(writer, "Hwloc support", hwlocDescription());
    writeLine(writer, "Tracing support", "disabled");

    // MDModules
    writeLine(writer, "Colvars support", colvarsDescription());
    writeLine(writer, "CP2K support", qmmmDescription());
    writeLine(writer, "Torch support", torchDescription());
    writeLine(writer, "Plumed support", plumedDescription());

    if (const std::optional<std::string> description = instrumentationApiDescription())
    {
        writeLine(writer, "Instrumentation API", description.value());
    }

    /* TODO: The below strings can be quite long, so it would be nice to wrap
     * them. Can wait for later, as the main branch has ready code to do all
     * that. */
    writeLine(writer, "C compiler", BUILD_C_COMPILER);
    writeLine(writer, "C compiler flags", std::string(BUILD_CFLAGS) + " " + CMAKE_BUILD_CONFIGURATION_C_FLAGS);
    writeLine(writer, "C++ compiler", BUILD_CXX_COMPILER);
    writeLine(writer, "C++ compiler flags", std::string(BUILD_CXXFLAGS) + " " + CMAKE_BUILD_CONFIGURATION_CXX_FLAGS);

    writeLine(writer, "BLAS library", blasDescription());
    writeLine(writer, "LAPACK library", lapackDescription());
    // Each GPU SDK populates only fields relevant to it
    writeOptionalLineFromMap(writer, "OpenCL include dir", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "OpenCL library", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "OpenCL version", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "CUDA compiler", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "CUDA compiler flags", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "CUDA driver", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "CUDA runtime", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "SYCL version", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "SYCL compiler", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "SYCL compiler flags", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "SYCL linker flags", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "SYCL GPU flags", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "SYCL targets", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "HIP compiler", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "HIP compiler flags", gpuDescriptionStrings);
    writeOptionalLineFromMap(writer, "HIP driver/runtime", gpuDescriptionStrings);
}

} // namespace

BinaryInformationSettings::BinaryInformationSettings() :
    bExtendedInfo_(false),
    bCopyright_(false),
    bProcessId_(false),
    bGeneratedByHeader_(false),
    prefix_(""),
    suffix_("")
{
}

void printBinaryInformation(FILE*                            fp,
                            const IProgramContext&           programContext,
                            const BinaryInformationSettings& settings)
{
    try
    {
        TextWriter writer(fp);
        printBinaryInformation(&writer, programContext, settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}

void printBinaryInformation(TextWriter*                      writer,
                            const IProgramContext&           programContext,
                            const BinaryInformationSettings& settings)
{
    // TODO Perhaps the writer could be configured with the prefix and
    // suffix strings from the settings?
    const char* prefix          = settings.prefix_;
    const char* suffix          = settings.suffix_;
    const char* precisionString = "";
#if GMX_DOUBLE
    precisionString = " (double precision)";
#endif
    const char* const name = programContext.displayName();
    if (settings.bGeneratedByHeader_)
    {
        writer->writeLine(formatString("%sCreated by:%s", prefix, suffix));
    }
    // TODO: It would be nice to know here whether we are really running a
    // Gromacs binary or some other binary that is calling Gromacs; we
    // could then print "%s is part of GROMACS" or some alternative text.
    std::string title = formatString(":-) GROMACS - %s, %s%s (-:", name, gmx_version(), precisionString);
    const int indent =
            centeringOffset(78 - std::strlen(prefix) - std::strlen(suffix), title.length()) + 1;
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
        writer->writeLine(formatString(
                "%sGROMACS:      %s, version %s%s%s", prefix, name, gmx_version(), precisionString, suffix));
    }
    const auto& binaryPath = programContext.fullBinaryPath();
    if (!binaryPath.empty())
    {
        writer->writeLine(formatString("%sExecutable:   %s%s", prefix, binaryPath.string().c_str(), suffix));
    }
    const gmx::InstallationPrefixInfo installPrefix = programContext.installationPrefix();
    if (!installPrefix.path_.empty())
    {
        writer->writeLine(formatString("%sData prefix:  %s%s%s",
                                       prefix,
                                       installPrefix.path_.string().c_str(),
                                       installPrefix.sourceLayoutTreeLike_ ? " (source tree)" : "",
                                       suffix));
    }
    const auto workingDir = std::filesystem::current_path();
    if (!workingDir.empty())
    {
        writer->writeLine(formatString("%sWorking dir:  %s%s", prefix, workingDir.string().c_str(), suffix));
    }
    if (settings.bProcessId_)
    {
        writer->writeLine(formatString("%sProcess ID:   %d%s", prefix, gmx_getpid(), suffix));
    }
    const char* const commandLine = programContext.commandLine();
    if (!gmx::isNullOrEmpty(commandLine))
    {
        writer->writeLine(formatString(
                "%sCommand line:%s\n%s  %s%s", prefix, suffix, prefix, commandLine, suffix));
    }
    if (settings.bExtendedInfo_)
    {
        GMX_RELEASE_ASSERT(prefix[0] == '\0' && suffix[0] == '\0',
                           "Prefix/suffix not supported with extended info");
        writer->ensureEmptyLine();
        writeExtendedInfo(writer);
    }
}

} // namespace gmx
