/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/// \cond internal
/*! \internal \file
 * \brief Declares Python runner
 *
 * \ingroup module_python
 */
#ifndef PYGMX_RUNNER_H
#define PYGMX_RUNNER_H


#include "gmxpre.h"

#include <memory>

#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/runner.h"

namespace gmx
{

namespace pyapi
{
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

class PyOptions;

/*! \brief Wraps Trajectory Analyis Runner for Python interface.
 *
 * Exported to Python as gmx.core.TafRunner
 * \internal \ingroup module_python
 */
class PyRunner
{
    public:
        /// Empty constructor not yet used.
        PyRunner() = delete;

        /*! \brief Construct runner with a single bound module.
         *
         * \param module existing module object to register with a new runner.
         */
        PyRunner(shared_ptr<gmx::TrajectoryAnalysisModule> module);

        virtual ~PyRunner();

        /*! \brief Process options.
         *
         * Allows options to be registered by the runner and bound modules.
         * \param options existing object to be provided by calling code.
         */
        void initialize(PyOptions &options);

        /*! \brief Advance the current frame one step.
         *
         * Returns when data dependencies on the next trajectory frame have been
         * satisfied. Updates internal state of runner to begin handling the next
         * frame of input, if any.
         * \return true while there are remaining frames to be handled, otherwise false.
         */
        bool next();

    private:
        /// has a common runner for most behavior
        gmx::trajectoryanalysis::Runner runner_;

        /// binds to one analysis module
        shared_ptr<gmx::TrajectoryAnalysisModule> module_;
};

};     //end namespace pyapi
};     // end namespace gmx
#endif // header guard
/// \endcond
