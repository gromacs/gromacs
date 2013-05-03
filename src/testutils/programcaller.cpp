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
 * Implements classes in programcaller.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gromacs/options/optionsassigner.h"
#include "programcaller.h"
#include "argsbuilder.h"

namespace gmx
{
namespace test
{

ProgramCaller::ProgramCaller(std::string const &programName)
    : programName_(programName),
      options_(programName.c_str(), "options to use when constructing a call to a GROMACS tool")
{
}

gmx::OptionInfo*
ProgramCaller::addOption(const AbstractOption &settings)
{
    return options_.addOption(settings);
}

bool
ProgramCaller::isDefined(const char *name) const
{
    return options_.isDefined(name);
}

void
ProgramCaller::setOption(const char *optionName, std::string const& value)
{
    gmx::OptionsAssigner assigner(&options_);

    assigner.start();
    assigner.startOption(optionName);
    assigner.appendValue(value);
    assigner.finishOption();
    assigner.finish();
}

gmx::Options&
ProgramCaller::getOptions()
{
    return options_;
}

std::string const &
ProgramCaller::getProgramName()
{
    return programName_;
}

} // namespace test
} // namespace gmx
