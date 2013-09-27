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
#include "programcaller.h"
#include "gromacs/utility/uniqueptr.h"
#include "gromacs/options/optionsvisitor.h"

namespace gmx
{
namespace test
{

/*! \brief Helper class to visit the options that were added to a
    ProgramCaller and build from those options a command line that
    will run the program with those options and their values. */
class ProgramCallerOptionsVisitor : public gmx::OptionsVisitor
{
    public:
        //! Constructor
        ProgramCallerOptionsVisitor(CommandLine *commandLine)
            : commandLine_(commandLine)
        {
        }
        //! Visit any (sub)sections of the gmx::Options
        void visitSubSection(const Options &section)
        {
            gmx::OptionsIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }
        //! Visit an option and pass it to commandLine_
        void visitOption(const OptionInfo &option)
        {
            commandLine_->append("-" + option.name());
            commandLine_->append(option.formatValue(0));
        }
        //! Object to fill during the visiting process.
        CommandLine *commandLine_;
};

ProgramCaller::ProgramCaller(const char *programName,
                             CMainFunction cmain)
    : options_(programName, "options to use when constructing a call to a GROMACS tool"),
      commandLine_(new CommandLine),
      bDoneOptionsVisiting_(false),
      cmainToRun_(cmain)
{
    commandLine_->append(programName);
}

void
ProgramCaller::addOption(const AbstractOption &settings)
{
    GMX_RELEASE_ASSERT(!bDoneOptionsVisiting_, "cannot add options after visiting");
    options_.addOption(settings);
}

int
ProgramCaller::run()
{
    if (!bDoneOptionsVisiting_)
    {
        ProgramCallerOptionsVisitor theVisitor(commandLine_.get());
        theVisitor.visitSubSection(options_);
        /* _options could be thrown away at this point, since it has
           served its purpose. */
        bDoneOptionsVisiting_ = true;
    }
    return cmainToRun_(commandLine_->argc(), commandLine_->argv());
}


} // namespace test
} // namespace gmx
