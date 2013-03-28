/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
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
 * Implements registerTrajectoryAnalysisModules().
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gromacs/trajectoryanalysis/modules.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"

#include "modules/angle.h"
#include "modules/distance.h"
#include "modules/freevolume.h"
#include "modules/select.h"

namespace gmx
{

namespace
{

/*! \internal \brief
 * Base class for trajectory analysis module command-line wrappers.
 *
 * This class exists to avoid multiple identical copies of the \p run() method
 * to be generated for the TrajAnalysisCmdLineWrapper template classes.
 *
 * \ingroup module_trajectoryanalysis
 */
class AbstractTrajAnalysisCmdLineWrapper : public CommandLineModuleInterface
{
    public:
        virtual const char *name() const             = 0;
        virtual const char *shortDescription() const = 0;

        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const HelpWriterContext &context) const;

    protected:
        //! Creates the analysis module for this wrapper.
        virtual TrajectoryAnalysisModulePointer createModule() const = 0;
};

int AbstractTrajAnalysisCmdLineWrapper::run(int argc, char *argv[])
{
    TrajectoryAnalysisModulePointer     module(createModule());
    TrajectoryAnalysisCommandLineRunner runner(module.get());
    runner.setPrintCopyright(false);
    return runner.run(argc, argv);
}

void AbstractTrajAnalysisCmdLineWrapper::writeHelp(const HelpWriterContext &context) const
{
    TrajectoryAnalysisModulePointer     module(createModule());
    TrajectoryAnalysisCommandLineRunner runner(module.get());
    runner.writeHelp(context);
}

/*! \internal \brief
 * Template for command-line wrapper of a trajectory analysis module.
 *
 * \tparam Module  Trajectory analysis module to wrap.
 *
 * \p Module should be default-constructible, derive from
 * TrajectoryAnalysisModule, and have a static public members
 * \c "const char name[]" and \c "const char shortDescription[]".
 *
 * \ingroup module_trajectoryanalysis
 */
template <class Module>
class TrajAnalysisCmdLineWrapper : public AbstractTrajAnalysisCmdLineWrapper
{
    public:
        virtual const char *name() const
        {
            return Module::name;
        }
        virtual const char *shortDescription() const
        {
            return Module::shortDescription;
        }
        virtual TrajectoryAnalysisModulePointer createModule() const
        {
            return TrajectoryAnalysisModulePointer(new Module);
        }
};

}   // namespace

void registerTrajectoryAnalysisModules(CommandLineModuleManager *manager)
{
    using namespace gmx::analysismodules;
    manager->registerModule<TrajAnalysisCmdLineWrapper<Angle> >();
    manager->registerModule<TrajAnalysisCmdLineWrapper<Distance> >();
    manager->registerModule<TrajAnalysisCmdLineWrapper<FreeVolume> >();
    manager->registerModule<TrajAnalysisCmdLineWrapper<Select> >();
}

} // namespace gmx
