/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Declares gmx::test::CommandLine.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_CMDLINETEST_H
#define GMX_TESTUTILS_CMDLINETEST_H

#include <string>

#include "gromacs/utility/common.h"

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Helper class for tests that check command-line handling.
 *
 * This class helps in writing tests for command-line handling.
 * The create() method takes an array of const char pointers, specifying the
 * command-line arguments, each as one array element.
 * The argc() and argv() methods can then be used to obtain the argc and argv
 * (non-const char pointers) arrays for passing into methods that expect these.
 *
 * Note that although the interface allows passing the argc and argv pointers
 * to methods that modify them (typically as \p f(&argc(), argv())), currently
 * the CommandLine object is not in a consistent state internally if the
 * parameters are actually modified.
 * Currently, the C++ methods with this signature do not modify their
 * parameters, so this is not yet a problem.
 *
 * All constructors and methods that modify this class may throw an
 * std::bad_alloc.  Const methods and accessors do not throw.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class CommandLine
{
    public:
        /*! \brief
         * Initializes a command-line object from a const C array.
         *
         * \param[in] cmdline  Array of command-line arguments.
         * \tparam    count    Deduced number of elements in \p cmdline.
         *
         * \p cmdline should include the binary name as the first element if
         * that is desired in the output.
         *
         * This is not a constructor, because template constructors are not
         * possible with a private implementation class.
         */
        template <size_t count> static
        CommandLine create(const char *const (&cmdline)[count])
        {
            return CommandLine(cmdline, count);
        }

        //! Initializes an empty command-line object.
        CommandLine();
        /*! \brief
         * Initializes a command-line object.
         *
         * \param[in] cmdline  Array of command-line arguments.
         * \param[in] count    Number of elements in \p cmdline.
         *
         * \p cmdline should include the binary name as the first element if
         * that is desired in the output.
         */
        CommandLine(const char *const cmdline[], size_t count);
        //! Creates a deep copy of a command-line object.
        CommandLine(const CommandLine &other);
        ~CommandLine();

        /*! \brief
         * Append an argument to the command line.
         *
         * \param[in] arg  Argument to append.
         *
         * Strong exception safety.
         */
        void append(const char *arg);
        //! Convenience overload taking a std::string.
        void append(const std::string &arg) { append(arg.c_str()); }

        //! Returns argc for passing into C-style command-line handling.
        int &argc();
        //! Returns argv for passing into C-style command-line handling.
        char **argv();
        //! Returns argc for passing into C-style command-line handling.
        int argc() const;
        //! Returns argv for passing into C-style command-line handling.
        const char *const *argv() const;
        //! Returns a single argument.
        const char *arg(int i) const;

        //! Returns the command line formatted as a single string.
        std::string toString() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
