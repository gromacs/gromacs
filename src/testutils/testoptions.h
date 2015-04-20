/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Functions for accessing test command-line options.
 *
 * Functions in this header allow accessing command-line options passed to the
 * test executable from tests.  This can be used to, e.g., enable additional
 * output for debugging purposes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTOPTIONS_H
#define GMX_TESTUTILS_TESTOPTIONS_H

namespace gmx
{

class Options;

namespace test
{

/*! \libinternal \brief
 * Provides additional options for the test executable.
 *
 * Typically not used directly in test code, but through the
 * #GMX_TEST_OPTIONS macro.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class TestOptionsProvider
{
    public:
        /*! \brief
         * Initializes the options from this provider.
         *
         * \param   options  The options need to be added here.
         */
        virtual void initOptions(Options *options) = 0;

    protected:
        virtual ~TestOptionsProvider() {}
};

/*! \libinternal \brief
 * Registers a test option provider with the test framework.
 *
 * \param[in] name     Name of the options provider (for ordering).
 * \param[in] provider The provider to register.
 * \throws  std::bad_alloc     if out of memory.
 * \throws  tMPI::system_error on mutex failures.
 *
 * Typically not used directly in test code, but through the
 * #GMX_TEST_OPTIONS macro.
 *
 * This gets called from constructors for global variables, so ideally
 * it would not throw to avoid unhandled exceptions.  But since this
 * is only test code, it is not worth the effort to try to remove those
 * rare exceptions (mutex failures and out-of-memory from STL).
 *
 * \ingroup module_testutils
 */
void registerTestOptions(const char *name, TestOptionsProvider *provider);
/*! \libinternal \brief
 * Initializes the options from all registered test providers.
 *
 * \param   options  The options are added here.
 *
 * This is called automatically by initTestUtils().
 *
 * \ingroup module_testutils
 */
void initTestOptions(Options *options);

// Uncrustify screws up the indentation for the example otherwise.
/* *INDENT-OFF* */
/*! \libinternal \brief
 * Macro to add additional command-line options for the test binary.
 *
 * \param  name    Unique name for the set of options.
 * \param  options Placeholder name for an gmx::Options object for adding options.
 *
 * The macro should be followed by a block that adds the desired command-line
 * options to `options` using gmx::Options::addOption().  \ref module_options
 * provides an overview of the options machinery.
 *
 * `name` must be unique within the executable to which the options are added.
 * If the macro is within an unnamed namespace, then it is sufficient that it
 * is unique within the file.
 *
 * Typical usage:
 * \code
   #include "gromacs/options/basicoptions.h"
   #include "gromacs/options/options.h"

   #include "testutils/testoptions.h"

   namespace gmx
   {
   namespace
   {

   bool g_optionValue = false;

   //! \cond
   GMX_TEST_OPTIONS(MyTestOptions, options)
   {
       options->addOption(BooleanOption("flag").store(&g_optionValue)
                              .description("My description"));
   }
   //! \endcond

   } // namespace
   } // namespace gmx
   \endcode
 *
 * \c \\cond and \c \\endcond statements are necessary around the macro to avoid
 * Doxygen warnings.
 *
 * One macro invocation per an added option, with more of the implementation
 * details hidden inside the macro, could be nicer.  But that requires more
 * elaborate macro machinery, so it is probably not worth the effort and
 * complexity.
 *
 * \ingroup module_testutils
 * \hideinitializer
 */
/* *INDENT-ON* */
#define GMX_TEST_OPTIONS(name, options) \
    class name : public ::gmx::test::TestOptionsProvider \
    { \
        public: \
            name() \
            { \
                ::gmx::test::registerTestOptions(#name, this); \
            } \
            virtual void initOptions(::gmx::Options *options); \
    }; \
    \
    static name s_ ## name ## Instance; \
    \
    void name::initOptions(::gmx::Options *options)

} // namespace test
} // namespace gmx

#endif
