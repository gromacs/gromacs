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

class LogWriteHelper
{
    public:
        explicit LogWriteHelper(ILogTarget *target) : target_(target) {}

        explicit operator bool() const { return target_ == NULL; }

        void operator=(const LogEntryWriter &entryWriter)
        {
            target_->writeEntry(entryWriter.entry_);
        }

    private:
        ILogTarget *target_;
};

class MDLogger
{
    public:
        enum class LogLevel
        {
            Warning,
            Info
        };
        static const int LogLevelCount = static_cast<int>(LogLevel::Info) + 1;

        MDLogger();
        explicit MDLogger(ILogTarget *targets[LogLevelCount]);

        LogWriteHelper warning() const { return getWriteHelper(LogLevel::Warning); }
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
 * \ingroup module_utility
 */
#define GMX_LOG(logger) \
    if (::gmx::LogWriteHelper helper = (logger)) { } else \
        helper = ::gmx::LogEntryWriter()

} // namespace gmx

#endif
