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
#ifndef GMX_FILEIO_WARNINP_H
#define GMX_FILEIO_WARNINP_H

#include <filesystem>
#include <string>
#include <string_view>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"

//! Different type of warnings
enum class WarningType : int
{
    //! Useful note to the user.
    Note,
    //! Warning about potential issue.
    Warning,
    //! Error that can't be recovered.
    Error,
    Count
};

/*! \brief
 * General warning handling object.
 *
 * If \p allowWarnings is false, all warnings (calls to warning()) will be
 * transformed into errors, calls to warning_note still produce notes.
 * \p maxNumberWarnings determines the maximum number of warnings that are allowed
 * for proceeding. When this number is exceeded check_warning_error
 * and done_warning will generate a fatal error.
 * allowWarnings should only be true in programs that have
 * a -maxwarn command line option.
 */
class WarningHandler
{
public:
    WarningHandler(bool allowWarnings, int maxNumberWarnings) :
        allowWarnings_(allowWarnings), maxNumberWarnings_(maxNumberWarnings)
    {
        if (maxNumberWarnings_ < 0)
        {
            GMX_THROW(gmx::InconsistentInputError(
                    "Max number of warnings need to be a positive integer"));
        }
    }
    /*!\brief Set \p fileName and \p lineNumber for the warning
     *
     * Note that \p fileName can be nullptr, leading to no file
     * information being printed.
     */
    void setFileAndLineNumber(const std::filesystem::path& fileName, int lineNumber);
    //! Get filename for the warning */
    std::filesystem::path getFileName() const;

    int errorCount() const { return numberOfEntries_[WarningType::Error]; }

    int warningCount() const { return numberOfEntries_[WarningType::Warning]; }

    int noteCount() const { return numberOfEntries_[WarningType::Note]; }

    int maxWarningCount() const { return maxNumberWarnings_; }

    /*! \brief
     * Issue a warning, with the string \p message.
     *
     * If \p message == nullptr, then warn_buf
     * will be printed instead. The file and line set by setLine
     * are printed, and the number of warnings is incremented.
     * A fatal error will be generated after processing the input
     * when number of warnings is larger than maxNumberWarning passed during construction.
     * So addWarning should only be called for issues that should be resolved,
     * otherwise addNote should be called.
     */
    void addWarning(std::string_view message);
    /*! \brief
     *  Issue a note, with the string \p message.
     *
     *   If \p message == nullptr, then warn_buf
     * will be printed instead. The file and line set by setLine
     * are printed, number of notes is incremented.
     * This is for issues which could be a problem for some systems,
     * but 100% ok for other systems.
     */
    void addNote(std::string_view message);
    /*! \brief
     *  Issue an error, with the string \p message.
     *
     *   If \p message == nullptr, then warn_buf
     * will be printed instead. The file and line set by setLine
     * are printed, number of errors is incremented.
     */
    void addError(std::string_view message);

private:
    //! Internal method for adding specific warning \p type
    void addLowLevel(std::string_view message, WarningType type);
    //! Whether warnings are distinct from errors or not.
    bool allowWarnings_;
    //! Storage for number of different entries.
    gmx::EnumerationArray<WarningType, int> numberOfEntries_ = { 0, 0, 0 };
    //! How many warnings can be tolerated.
    int maxNumberWarnings_;
    //! Which line number corresponds to the message.
    int lineNumber_ = -1;
    //! Which file the message is coming from.
    std::filesystem::path fileName_ = "unknown";
};

/*! \brief Issue an error with warning_error() and prevent further
 * processing by calling check_warning_error().
 *
 * This is intended for use where there is no way to produce a data
 * structure that would prevent execution from segfaulting. */
[[noreturn]] void warning_error_and_exit(WarningHandler*              wi,
                                         const char*                  s,
                                         int                          f_errno,
                                         const std::filesystem::path& file,
                                         int                          line);
//! \copydoc warning_error_and_exit(WarningHandler*,const char *,int,const char *,int)
[[noreturn]] void warning_error_and_exit(WarningHandler*              wi,
                                         const std::string&           s,
                                         int                          f_errno,
                                         const std::filesystem::path& file,
                                         int                          line);

//! Return whether any error-level warnings were issued to \p wi. */
bool warning_errors_exist(const WarningHandler& wi);

void check_warning_error(const WarningHandler& wi, int f_errno, const std::filesystem::path& file, int line);
/* When warning_error has been called at least once gmx_fatal is called,
 * otherwise does nothing.
 */

void done_warning(const WarningHandler& wi, int f_errno, const std::filesystem::path& file, int line);
/* Should be called when finished processing the input file.
 * Prints the number of notes and warnings
 * and generates a fatal error when errors were found or too many
 * warnings were generatesd.
 * Frees the data structure pointed to by wi.
 */


void too_few_function(WarningHandler* wi, const std::filesystem::path& fn, int line);
#define too_few(wi) too_few_function(wi, __FILE__, __LINE__)
/* Issue a warning stating 'Too few parameters' */

void incorrect_n_param_function(WarningHandler* wi, const std::filesystem::path& fn, int line);
#define incorrect_n_param(wi) incorrect_n_param_function(wi, __FILE__, __LINE__)
/* Issue a warning stating 'Incorrect number of parameters' */

#endif
