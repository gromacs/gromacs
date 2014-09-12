/*
 * alexandria - driver program for tools that are part of the
 * Alexandria force field.
 * Output files are compatible with GROMACS 5.0 or newer.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements the alexandria wrapper binary.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/selection/selhelp.h"
#include "gromacs/utility/exceptions.h"

#include "alex_modules.h"

int
main(int argc, char *argv[])
{
    gmx::CommandLineProgramContext &context = gmx::initForCommandLine(&argc, &argv);
    try
    {
        gmx::CommandLineModuleManager manager("alexandria", &context);
        registerAlexandriaModules(&manager);
        manager.addHelpTopic(gmx::createSelectionHelpTopic());
        manager.setQuiet(true);
        printf("\n                   Welcome to Alexandria\n\n");
        printf("Copyright (c) 2014, David van der Spoel and Paul J. van Maaren\n");
        printf("See http://folding.bmc.uu.se/ for details.\n\n");
        printf("Alexandria is free software under the Gnu Public License v 2.\n");
        printf("Read more at http://www.gnu.org/licenses/gpl-2.0.html\n\n");
        int rc = manager.run(argc, argv);
        gmx::finalizeForCommandLine();
        return rc;
    }
    catch (const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return gmx::processExceptionAtExit(ex);
    }
}
