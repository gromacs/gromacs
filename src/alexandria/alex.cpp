/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements the alex wrapper binary.
 * Copied straight from the gmx wrapper binary
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/trajectoryanalysis/modules.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/init.h"
#include "gromacs/onlinehelp/helptopicinterface.h"

#include "alex_modules.h"

namespace gmx {

class AlexHelpTopicInterface : HelpTopicInterface
{
    public:
        ~AlexHelpTopicInterface() {}

        /*! \brief
         * Returns the name of the topic.
         *
         * This should be a single lowercase word, used to identify the topic.
         * It is not used for the root of the help topic tree.
         */
        const char *name() const { return "alexandria"; }
        /*! \brief
         * Returns a title for the topic.
         *
         * May return NULL, in which case the topic is omitted from normal
         * subtopic lists and no title is printed by the methods provided in
         * helptopic.h.
         */
        const char *title() const { return "alexandria"; }

        //! Returns whether the topic has any subtopics.
        bool hasSubTopics() const { return false; }
        /*! \brief
         * Finds a subtopic by name.
         *
         * \param[in] name  Name of subtopic to find.
         * \returns   Pointer to the found subtopic, or NULL if matching topic
         *      is not found.
         */
        const HelpTopicInterface *findSubTopic(const char *name) { return NULL; }

        /*! \brief
         * Prints the help text for this topic.
         *
         * \param[in] context  Context object for writing the help.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        void writeHelp(const HelpWriterContext &context);
};

void AlexHelpTopicInterface::writeHelp(const HelpWriterContext &context)
{
    printf("BOE\n");
}

};

int
main(int argc, char *argv[])
{
    gmx::ProgramInfo &info = gmx::init("gmx", &argc, &argv);
    try
    {
        gmx::CommandLineModuleManager manager(&info);
        registerAlexandriaModules(&manager);
        
        manager.addHelpTopic(gmx::SelectionCollection::createDefaultHelpTopic());
        int rc = manager.run(argc, argv);
        gmx::finalize();
        return rc;
    }
    catch (const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return gmx::processExceptionAtExit(ex);
    }
}
