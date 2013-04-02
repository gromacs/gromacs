/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares pimpl classes from integrationtests.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_integration_tests
 */
#include "integrationtests.h"

#include <boost/shared_ptr.hpp>
#include <map>

namespace gmx
{
namespace test
{
 
/********************************************************************
 * IntegrationTestFixture::Impl
 */

class IntegrationTestFixture::Impl
{
public:
    explicit Impl(IntegrationTestFixture *parent);

    /*! \brief Return the object that manaages the call to a GROMACS program.
     *
     * Implements lazy initialization of such objects.
     *
     * \param[in] programName  Name of the GROMACS program
     */
    boost::shared_ptr<ProgramCaller> getProgramCaller(std::string const& programName);
    /*! \brief Factory method to build an object that can manage
     *  the options for and call to a GROMACS program via its
     *  cmain function.
     *
     * This function flavour does not define any options for any
     * tools.
     *
     * \param[in] programName  Name of the GROMACS program
     * \param[in] cmain        The cmain function of that GROMACS program
     */
    void createProgramCaller(std::string const& programName, ProgramCaller::cmain_func_ptr cmain);
    /*! /brief Map from a string to a smart pointer to Options */
    typedef std::map<std::string, boost::shared_ptr<ProgramCaller> > ProgramCallerMap;
    /*! /brief Map containing the options from which to build
     * arguments to pass to some GROMACS program. The key is the
     * program name, and the value is an object that manages the
     * call to that program.
     *
     * There's no copy constructor for Options, so we have to
     * build a (smart) pointer to one to use it in a container.
     */
    ProgramCallerMap programCallers;
    /*! \brief Defines an option that will be used with the named program.
     *
     * The OptionType must implement the AbstractOption
     * interface. Intended to be used for pre-configuring
     * "standard options" for commonly-used programs.
     */
    template <class OptionType> void
    defineProgramOption(std::string const& programName, const char* optionName, char const* description);
    template <class OptionType> void
    setProgramOption(std::string const& programName, const char *optionName, std::string const& value);
    void setProgramOption(std::string const& programName, const char *optionName, std::string const& value);

    IntegrationTestFixture& parent_;
};

} // namespace test
} // namespace gmx
