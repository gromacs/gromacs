/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Declares gmx::test::CommandLine.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_CMDLINETEST_H
#define GMX_TESTUTILS_CMDLINETEST_H

#include <string>

// arrayref.h is not strictly necessary for this header, but nearly all
// callers will need it to use the constructor that takes ConstArrayRef.
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/common.h"

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Helper class for tests that check command-line handling.
 *
 * This class helps in writing tests for command-line handling.
 * The constructor method takes an array of const char pointers, specifying the
 * command-line arguments, each as one array element.  It is also possible to
 * construct the command line by adding individual arguments with append() and
 * addOption().
 * The argc() and argv() methods can then be used to obtain `argc` and `argv`
 * (non-const char pointers) arrays for passing into methods that expect these.
 *
 * Note that although the interface allows passing the argc and argv pointers
 * to methods that modify them (typically as \p f(&argc(), argv())), currently
 * the CommandLine object is not in a consistent state internally if the
 * parameters are actually modified.  Reading the command line is possible
 * afterwards, but modification is not.
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
        //! Initializes an empty command-line object.
        CommandLine();
        /*! \brief
         * Initializes a command-line object from an array.
         *
         * \param[in] cmdline  Array of command-line arguments.
         *
         * \p cmdline should include the binary name as the first element if
         * that is desired in the output.
         *
         * This constructor is not explicit to make it possible to create a
         * CommandLine object directly from a C array.
         */
        CommandLine(const ConstArrayRef<const char *> &cmdline);
        //! Creates a deep copy of a command-line object.
        CommandLine(const CommandLine &other);
        ~CommandLine();

        /*! \brief
         * Initializes a command-line object in-place from an array.
         *
         * \param[in] cmdline  Array of command-line arguments.
         *
         * \p cmdline should include the binary name as the first element if
         * that is desired in the output.
         *
         * This function does the same as the constructor that takes a
         * ConstArrayRef.  Any earlier contents of the object are discarded.
         *
         * Strong exception safety.
         */
        void initFromArray(const ConstArrayRef<const char *> &cmdline);

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
        /*! \brief
         * Add an option-value pair to the command line.
         *
         * \param[in] name   Name of the option to append, which
         *                   should start with "-".
         * \param[in] value  Value of the argument to append.
         */
        void addOption(const char *name, const char *value);
        //! Convenience overload taking a std::string.
        void addOption(const char *name, const std::string &value);
        //! Overload taking an int.
        void addOption(const char *name, int value);
        //! Overload taking a double.
        void addOption(const char *name, double value);
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
