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
#include "gmxpre.h"

#include "warninp.h"

#include <cstdio>
#include <cstring>

#include <string>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

static const char* warningTypeString(const WarningType type)
{
    static constexpr gmx::EnumerationArray<WarningType, const char*> warningTypeName = { "NOTE",
                                                                                         "WARNING",
                                                                                         "ERROR" };
    return warningTypeName[type];
}

void WarningHandler::setFileAndLineNumber(const std::filesystem::path& fileName, int lineNumber)
{
    if (!fileName.empty())
    {
        fileName_ = fileName;
    }
    lineNumber_ = lineNumber;
}

std::filesystem::path WarningHandler::getFileName() const
{
    return fileName_;
}

void WarningHandler::addLowLevel(std::string_view message, const WarningType type)
{
    static constexpr int indent = 2;
    if (message.empty())
    {
        message = "Empty error message.";
    }

    gmx::TextLineWrapperSettings lineWrapperSettings{};
    lineWrapperSettings.setIndent(indent);
    lineWrapperSettings.setLineLength(77 - indent);
    lineWrapperSettings.setFirstLineIndent(2);
    gmx::TextLineWrapper wrapper(lineWrapperSettings);
    const std::string    wrapped = wrapper.wrapToString(std::string(message));
    if (!fileName_.empty())
    {
        if (lineNumber_ != -1)
        {
            fprintf(stderr,
                    "\n%s %d [file %s, line %d]:\n%s\n\n",
                    warningTypeString(type),
                    numberOfEntries_[type],
                    fileName_.string().c_str(),
                    lineNumber_,
                    wrapped.c_str());
        }
        else
        {
            fprintf(stderr,
                    "\n%s %d [file %s]:\n%s\n\n",
                    warningTypeString(type),
                    numberOfEntries_[type],
                    fileName_.string().c_str(),
                    wrapped.c_str());
        }
    }
    else
    {
        fprintf(stderr, "\n%s %d:\n%s\n\n", warningTypeString(type), numberOfEntries_[type], wrapped.c_str());
    }
}

void WarningHandler::addWarning(std::string_view message)
{
    if (allowWarnings_)
    {
        numberOfEntries_[WarningType::Warning]++;
        addLowLevel(message, WarningType::Warning);
    }
    else
    {
        addError(message);
    }
}

void WarningHandler::addError(std::string_view message)
{
    numberOfEntries_[WarningType::Error]++;
    addLowLevel(message, WarningType::Error);
}

void WarningHandler::addNote(std::string_view message)
{
    numberOfEntries_[WarningType::Note]++;
    addLowLevel(message, WarningType::Note);
}

static void printWarningCount(const WarningType type, int n)
{
    if (n > 0)
    {
        fprintf(stderr,
                "\nThere %s %d %s%s\n",
                (n == 1) ? "was" : "were",
                n,
                warningTypeString(type),
                (n == 1) ? "" : "s");
    }
}

// Note it is the caller's responsibility to ensure that exiting is correct behaviour
[[noreturn]] static void check_warning_error_impl(const WarningHandler&        wi,
                                                  int                          f_errno,
                                                  const std::filesystem::path& file,
                                                  int                          line)
{
    printWarningCount(WarningType::Note, wi.noteCount());
    printWarningCount(WarningType::Warning, wi.warningCount());

    gmx_fatal(f_errno,
              file.c_str(),
              line,
              "There %s %d error%s in input file(s)",
              (wi.errorCount() == 1) ? "was" : "were",
              wi.errorCount(),
              (wi.errorCount() == 1) ? "" : "s");
}

void check_warning_error(const WarningHandler& wi, int f_errno, const std::filesystem::path& file, int line)
{
    if (wi.errorCount() > 0)
    {
        check_warning_error_impl(wi, f_errno, file, line);
    }
}

void warning_error_and_exit(WarningHandler* wi, const char* s, int f_errno, const std::filesystem::path& file, int line)
{
    wi->addError(s);
    check_warning_error_impl(*wi, f_errno, file, line);
}

void warning_error_and_exit(WarningHandler*              wi,
                            const std::string&           s,
                            int                          f_errno,
                            const std::filesystem::path& file,
                            int                          line)
{
    warning_error_and_exit(wi, s.c_str(), f_errno, file, line);
}

bool warning_errors_exist(const WarningHandler& wi)
{
    return (wi.errorCount() > 0);
}

void done_warning(const WarningHandler& wi, int f_errno, const std::filesystem::path& file, int line)
{
    // If we've had an error, then this will report the number of
    // notes and warnings, and then exit.
    check_warning_error(wi, f_errno, file, line);

    // Otherwise, we report the number of notes and warnings.
    printWarningCount(WarningType::Note, wi.noteCount());
    printWarningCount(WarningType::Warning, wi.warningCount());

    if (wi.warningCount() > wi.maxWarningCount())
    {
        gmx_fatal(f_errno,
                  file.c_str(),
                  line,
                  "Too many warnings (%d).\n"
                  "If you are sure all warnings are harmless, use the -maxwarn option.",
                  wi.warningCount());
    }
}

void too_few_function(WarningHandler* wi, const std::filesystem::path& fn, int line)
{
    wi->addWarning(gmx::formatString(
            "Too few parameters on line (source file %s, line %d)", fn.string().c_str(), line));
}

void incorrect_n_param_function(WarningHandler* wi, const std::filesystem::path& fn, int line)
{
    wi->addWarning(gmx::formatString(
            "Incorrect number of parameters on line (source file %s, line %d)", fn.string().c_str(), line));
}
