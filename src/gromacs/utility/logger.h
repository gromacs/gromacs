/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares functionality for logging.
 *
 * See \ref page_logging for an overview of the functionality.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_LOGGER_H
#define GMX_UTILITY_LOGGER_H

#include <string>

#include "gromacs/utility/stringutil.h"

namespace gmx
{

struct LogEntry
{
    LogEntry() : asParagraph(false) {}

    std::string text;
    bool        asParagraph;
};

/*! \libinternal \brief
 * Target where log output can be written.
 *
 * \ingroup module_utility
 */
class ILogTarget
{
public:
    virtual ~ILogTarget();

    //! Writes a log entry to this target.
    virtual void writeEntry(const LogEntry& entry) = 0;
};

/*! \libinternal \brief
 * Helper class for creating log entries with ::GMX_LOG.
 *
 * \ingroup module_utility
 */
class LogEntryWriter
{
public:
    //! Appends given text as a line in the log entry.
    LogEntryWriter& appendText(const char* text)
    {
        entry_.text.append(text);
        return *this;
    }
    //! Appends given text as a line in the log entry.
    LogEntryWriter& appendText(const std::string& text)
    {
        entry_.text.append(text);
        return *this;
    }
    //! Appends given text as a line in the log entry, with printf-style formatting.
    LogEntryWriter& appendTextFormatted(gmx_fmtstr const char* fmt, ...) gmx_format(printf, 2, 3);
    //! Writes the log entry with empty lines before and after.
    LogEntryWriter& asParagraph()
    {
        entry_.asParagraph = true;
        return *this;
    }

private:
    LogEntry entry_;

    friend class LogWriteHelper;
};

/*! \internal \brief
 * Helper class for implementing ::GMX_LOG.
 *
 * \ingroup module_utility
 */
class LogWriteHelper
{
public:
    //! Initializes a helper for writing to the given target.
    explicit LogWriteHelper(ILogTarget* target) : target_(target) {}

    // Should be explicit, once that works in CUDA.
    /*! \brief
     * Returns whether anything needs to be written.
     *
     * Note that the return value is unintuitively `false` when the target
     * is active, to allow implementing ::GMX_LOG like it is now.
     */
    operator bool() const { return target_ == nullptr; }

    /*! \brief
     * Writes the entry from the given writer to the log target.
     *
     * This is implemented as an assignment operator to get proper
     * precedence for operations for the ::GMX_LOG macro; this is a common
     * technique for implementing macros that allow streming information to
     * them (see, e.g., Google Test).
     */
    LogWriteHelper& operator=(const LogEntryWriter& entryWriter)
    {
        target_->writeEntry(entryWriter.entry_);
        return *this;
    }

private:
    ILogTarget* target_;
};

/*! \libinternal \brief
 * Represents a single logging level.
 *
 * Typically this type is not used directly, but instances in MDLogger are
 * simply accessed through ::GMX_LOG in code that writes to the log.
 *
 * \ingroup module_utility
 */
class LogLevelHelper
{
public:
    //! Initializes a helper for writing to the given target.
    explicit LogLevelHelper(ILogTarget* target) : target_(target) {}

    // Both of the below should be explicit, once that works in CUDA.
    //! Returns whether the output for this log level goes anywhere.
    operator bool() const { return target_ != nullptr; }

    //! Creates a helper for ::GMX_LOG.
    operator LogWriteHelper() const { return LogWriteHelper(target_); }

private:
    ILogTarget* target_;
};

/*! \libinternal \brief
 * Declares a logging interface.
 *
 * Typically, this object is not created directly, but instead through
 * LoggerBuilder.
 *
 * For now, this is named MDLogger, since it is used only there, and it is not
 * clear whether the logging levels can be the same throughout the code.  It
 * should be relatively straightforward to split this into multiple classes
 * with different supported logging levels without changing calling code, or to
 * rename it to Logger if we do not need any specialization.
 *
 * \ingroup module_utility
 */
class MDLogger
{
public:
    //! Supported logging levels.
    enum class LogLevel
    {
        Error,
        Warning,
        Info,
        Debug,
        VerboseDebug,
        Count
    };
    //! Number of logging levels.
    static const int LogLevelCount = static_cast<int>(LogLevel::Count);

    MDLogger();
    //! Creates a logger with the given targets.
    explicit MDLogger(ILogTarget* targets[LogLevelCount]);

    //! For writing at LogLevel::Warning level.
    LogLevelHelper warning;
    //! For writing at LogLevel::Error level.
    LogLevelHelper error;
    //! For writing at LogLevel::Debug level.
    LogLevelHelper debug;
    //! For writing at LogLevel::VerboseDebug level.
    LogLevelHelper verboseDebug;
    //! For writing at LogLevel::Info level.
    LogLevelHelper info;
};

/*! \brief
 * Helper to log information using gmx::MDLogger.
 *
 * \param  logger  LogLevelHelper instance to use for logging.
 *
 * Used as
 * \code
   GMX_LOG(logger.warning).appendText(...);
   \endcode
 * and ensures that the code to format the output is only executed when the
 * output goes somewhere.
 *
 * See LogEntryWriter for functions that can be used with the macro (such as
 * the appendText() in the example).
 *
 * \ingroup module_utility
 */
#define GMX_LOG(logger)                                                  \
    if (::gmx::LogWriteHelper helper = ::gmx::LogWriteHelper(logger)) {} \
    else                                                                 \
        helper = ::gmx::LogEntryWriter()

} // namespace gmx

#endif
