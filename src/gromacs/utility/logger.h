/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares functionality for logging.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_LOGGER_H
#define GMX_UTILITY_LOGGER_H

#include <string>

namespace gmx
{

struct LogEntry
{
    LogEntry() : asParagraph(false) {}

    std::string text;
    bool        asParagraph;
};

class ILogTarget
{
    public:
        virtual ~ILogTarget();

        //! Writes a log entry to this target.
        virtual void writeEntry(const LogEntry &entry) = 0;
};

class LogEntryWriter
{
    public:
        LogEntryWriter &appendText(const char *text)
        {
            entry_.text.append(text);
            return *this;
        }
        LogEntryWriter &appendTextFormatted(const char *fmt, ...);
        LogEntryWriter &asParagraph()
        {
            entry_.asParagraph = true;
            return *this;
        }

    private:
        LogEntry    entry_;

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
        explicit LogWriteHelper(ILogTarget *target) : target_(target) {}

        /*! \brief
         * Returhs whether anything needs to be written.
         *
         * Note that the return value is unintuitively `false` when the target
         * is active, to allow implementing ::GMX_LOG like it is now.
         */
        explicit operator bool() const { return target_ == NULL; }

        /*! \brief
         * Writes the entry from the given writer to the log target.
         *
         * This is implemented as an assignment operator to get proper
         * precedence for operations for the ::GMX_LOG macro; this is a common
         * technique for implementing macros that allow streming information to
         * them (see, e.g., Google Test).
         */
        void operator=(const LogEntryWriter &entryWriter)
        {
            target_->writeEntry(entryWriter.entry_);
        }

    private:
        ILogTarget *target_;
};

/*! \libinternal \brief
 * Declares a logging interface.
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
        enum class LogLevel
        {
            Warning,
            Info
        };
        //! Number of logging levels.
        static const int LogLevelCount = static_cast<int>(LogLevel::Info) + 1;

        MDLogger();
        //! Creates a logger with the given targets.
        explicit MDLogger(ILogTarget *targets[LogLevelCount]);

        //! Returns a writer for writing an entry at LogLevel::Warning level.
        LogWriteHelper warning() const { return getWriteHelper(LogLevel::Warning); }
        //! Returns a writer for writing an entry at LogLevel::Info level.
        LogWriteHelper info() const { return getWriteHelper(LogLevel::Info); }

    private:
        LogWriteHelper getWriteHelper(LogLevel level) const
        {
            return LogWriteHelper(dest_[static_cast<int>(level)]);
        }

        ILogTarget *dest_[LogLevelCount];
};

/*! \brief
 * Helper to log information using gmx::MDLogger.
 *
 * \param  logger  LogWriteHelper instance to use for logging.
 *
 * Used as
 * \code
   GMX_LOG(logger.warning()).appendText(...);
   \endcode
 * and ensures that the code to format the output is only executed when the
 * output goes somewhere.
 *
 * See LogEntryWriter for functions that can be used with the macro (such as
 * the appendText() in the example).
 *
 * \ingroup module_utility
 */
#define GMX_LOG(logger) \
    if (::gmx::LogWriteHelper helper = (logger)) { } else \
        helper = ::gmx::LogEntryWriter()

} // namespace gmx

#endif
