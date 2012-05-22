/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \libinternal \file
 * \brief
 * Declares gmx::test::CommandLine.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_CMDLINETEST_H
#define GMX_TESTUTILS_CMDLINETEST_H

#include <cstddef>

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * Helper class for tests that check command-line handling.
 *
 * This class helps in writing tests for command-line handling.
 * The constructor takes an array of const char pointers, specifying the
 * command-line arguments, each as one array element.
 * The argc() and argv() methods can then be used to obtain the argc and argv
 * (non-const char pointers) arrays for passing into methods that expect these.
 *
 * Note that although the interface allows passing the argc and argv pointers
 * to methods that modify them (typically as \p f(&argc(), argv())), currently
 * memory leaks occur if the parameters are actually modified.
 * Currently, the C++ methods with this signature do not modify their
 * parameters, so this is not yet a problem.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class CommandLine
{
    public:
        /*! \brief
         * Initializes a command-line object.
         *
         * \param[in] cmdline  Array of command-line arguments.
         * \tparam    count Deduced number of elements in \p cmdline.
         * \throws    std::bad_alloc if out of memory.
         *
         * \p cmdline should include the binary name as the first element if
         * that is desired in the output.
         */
        template <size_t count>
        explicit CommandLine(const char *const (&cmdline)[count])
        {
            initCommandLine(cmdline, count);
        }
        ~CommandLine();

        //! Returns argc for passing into C-style command-line handling.
        int &argc() { return argc_; }
        //! Returns argv for passing into C-style command-line handling.
        char **argv() { return argv_; }
        //! Returns argc for passing into C-style command-line handling.
        int argc() const { return argc_; }
        //! Returns argv for passing into C-style command-line handling.
        const char *const *argv() const { return argv_; }

    private:
        //! Internal helper method used to implement the constructor.
        void initCommandLine(const char *const cmdline[], size_t count);

        int                     argc_;
        char                  **argv_;
};

} // namespace test
} // namespace gmx

#endif
