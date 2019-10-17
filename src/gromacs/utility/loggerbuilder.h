/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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
 * Declares functionality for initializing logging.
 *
 * See \ref page_logging for an overview of the functionality.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_LOGGERBUILDER_H
#define GMX_UTILITY_LOGGERBUILDER_H

#include <memory>
#include <string>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

class TextOutputStream;

class LoggerFormatterBuilder;
class LoggerOwner;

/*! \libinternal \brief
 * Initializes loggers.
 *
 * This class provides methods for specifying logging targets for a logger and
 * building the logger after all targets have been specified.  Having this
 * separate from the logger allows using different internal data structures
 * during initialization and operation, and simplifies the responsibilities of
 * the involved classes.
 *
 * \ingroup module_utility
 */
class LoggerBuilder
{
public:
    LoggerBuilder();
    ~LoggerBuilder();

    /*! \brief
     * Adds a stream to which log output is written.
     *
     * All output at level \p level or above it is written to \p stream.
     * The caller is responsible of closing and freeing \p stream once the
     * logger is discarded.
     */
    void addTargetStream(MDLogger::LogLevel level, TextOutputStream* stream);
    /*! \brief
     * Adds a file to which log output is written.
     *
     * All output at level \p level or above it is written to \p fp.
     * The caller is responsible of closing \p fp once the logger is
     * discarded.
     */
    void addTargetFile(MDLogger::LogLevel level, FILE* fp);

    /*! \brief
     * Builds the logger with the targets set for this builder.
     *
     * After this function has been called, the builder can (and should) be
     * discarded.
     */
    LoggerOwner build();

private:
    class Impl;

    PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief
 * Manages memory for a logger built with LoggerBuilder.
 *
 * This class is responsible of managing all memory allocated by LoggerBuilder
 * that is needed for operation of the actual logger.  Also the actual logger
 * instance is owned by this class.  This allows keeing the actual logger simple
 * and streamlined.
 *
 * This class supports move construction and assignment, which allows
 * initializing it on the stack and assigning a new instance if the targets
 * need to be changed.
 *
 * \ingroup module_utility
 */
class LoggerOwner
{
public:
    //! Move-constructs the owner.
    LoggerOwner(LoggerOwner&& other) noexcept;
    ~LoggerOwner();

    //! Move-assings the owner.
    LoggerOwner& operator=(LoggerOwner&& other) noexcept;

    //! Returns the logger for writing the logs.
    const MDLogger& logger() const { return *logger_; }

private:
    class Impl;

    LoggerOwner(std::unique_ptr<Impl> impl);

    PrivateImplPointer<Impl> impl_;
    const MDLogger*          logger_;

    friend class LoggerBuilder;
};

} // namespace gmx

#endif
