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
