/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Provides helper class for constructing error messages for file I/O
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_FILEIO_ERRORMESSAGEMAKER_H
#define GMX_FILEIO_ERRORMESSAGEMAKER_H

#include <string>

#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \brief Helps to construct an `InvalidInputException` with an error
 * message, and adds custom context information useful for
 * trouble-shooting issues with reading line-based input from files.
 *
 * Example of normal usage:
 *
 * \code
   std::string       line;
   ErrorMessageMaker errorMessage(filename);
   while(getNewLine(&line))
   {
       errorMessage.setCurrentLine(&line)
       if(isCommentLine(line))
       {
           continue;
       }
       if(doesntParse(line))
       {
           GMX_THROW(errorMessage.make("Reading failed because <reasons>"));
       }
       doWork(line);
   }
   \endcode
 *
 * This will produce an error message that reads
 * \verbatim
   Reading failed because <reasons>
    on line <number> which read
    '<the actual line text>'
   \endverbatim
 */
class ErrorMessageMaker
{
    public:
        //! Constructor
        explicit ErrorMessageMaker();

        /*! \brief As each line is read and parsing begins, call this
         * function to prepare for an error message that might be
         * given later.
         *
         * The caller is responsible for ensuring that the lifetime of
         * \c line is longer than the lifetime of this object, but
         * this is normally fine. */
        void setCurrentLine(const std::string *line);
        /*! \brief Call to construct an exception with a contextualized error message
         *
         * Normally, \c setCurrentLine should already have been called
         * for each line read, so that the line-index counter is
         * accurate. */
        InvalidInputError make(const std::string &message) const;
    private:
        int                lineIndex_;
        const std::string *currentLine_;
};

} // namespace

#endif
