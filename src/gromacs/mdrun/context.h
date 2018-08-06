/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
//
// Created by Eric Irrgang on 5/16/18.
//

#ifndef GROMACS_CONTEXT_H
#define GROMACS_CONTEXT_H

// Ugh... avoiding this header dependency is as ugly as not...
#include "gromacs/mdlib/simulationsignal.h"

namespace gmx
{

class Mdrunner;

namespace md
{

/*!
 * \brief Encapsulate some runtime context for sharing in the mdrun call stack.
 *
 * In the future, this functionality can be moved to an updated ProgramContext and
 * the Context should only provide high-level or external information directly. Its
 * primary purpose will be to register and hold factory function pointers with which
 * callers can get handles to the resources they need.
 *
 * Since those modules and resources don't exist yet, we're providing a few directly.
 */
class Context
{
    public:
        /*!
         * \brief Construct with the runner's one resource: a pointer to the owning runner.
         *
         * The Context should be owned by a runner and the Context lifetime should be entirely
         * within the Runner's life.
         *
         * \param runner non-owning pointer to the runner that owns the Context object.
         */
        explicit Context(const Mdrunner &runner);

        /*!
         * \brief Get a reference to the current array of signal flags.
         *
         * There is no guarantee that the flags have been initialized yet.
         *
         * \return pointer to signals array.
         */
        SimulationSignals * simulationSignals() const;

    private:
        const gmx::Mdrunner* runner_ {nullptr};
};

}      // end namespace md
}      // end namespace gmx

#endif //GROMACS_CONTEXT_H
