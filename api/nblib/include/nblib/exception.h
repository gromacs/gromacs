/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \inpublicapi \file
 * \brief
 * Implements nblib exception class
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_EXCEPTION_H
#define NBLIB_EXCEPTION_H

#include <exception>
#include <string>

namespace nblib
{

/*! \brief Base nblib exception class
 *
 * All nblib exceptions derive from this class and simply forward their message. This allows
 * exceptions to be handled uniformly across different exception types.
 */
class NbLibException : public std::exception
{
public:
    [[maybe_unused]] explicit NbLibException(const std::string& message) :
        message_("NbLib Exception: " + message)
    {
    }

    //! Overrides the what() in std::exception
    [[nodiscard]] const char* what() const noexcept override { return message_.c_str(); }

    //! Convenience call in case a string is wanted instead of a const char*
    [[nodiscard]] const std::string& reason() const& { return message_; }

private:
    std::string message_;
};

/*! \brief The exception type for user input errors
 *
 * The message should give users some hint as to how to remedy the error.
 */
class InputException final : public NbLibException
{
public:
    using NbLibException::NbLibException;
};

} // namespace nblib
#endif // NBLIB_EXCEPTION_H
