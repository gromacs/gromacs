/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
/*! \internal
 * \brief Declares the simulator interface for mdrun
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_LEGACYSIMULATOR_H
#define GMX_MDRUN_LEGACYSIMULATOR_H

#include <cstdio>

#include <memory>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "isimulator.h"

namespace gmx
{

//! Function type for simulator code.
using SimulatorFunctionType = void();

/*! \internal
 * \brief Struct to handle setting up and running the different simulation types.
 *
 * This struct is a mere aggregate of parameters to pass to run a
 * simulation, so that future changes to names and types of them consume
 * less time when refactoring other code.
 *
 * Having multiple simulation types as member functions isn't a good
 * design, and we definitely only intend one to be called, but the
 * goal is to make it easy to change the names and types of members
 * without having to make identical changes in several places in the
 * code. Once many of them have become modules, we should change this
 * approach.
 */
class LegacySimulator : public ISimulator
{
private:
    //! Implements the normal MD simulations.
    SimulatorFunctionType do_md;
    //! Implements the rerun functionality.
    SimulatorFunctionType do_rerun;
    //! Implements steepest descent EM.
    SimulatorFunctionType do_steep;
    //! Implements conjugate gradient energy minimization
    SimulatorFunctionType do_cg;
    //! Implements onjugate gradient energy minimization using the L-BFGS algorithm
    SimulatorFunctionType do_lbfgs;
    //! Implements normal mode analysis
    SimulatorFunctionType do_nm;
    //! Implements test particle insertion
    SimulatorFunctionType do_tpi;
    //! Implements MiMiC QM/MM workflow
    SimulatorFunctionType do_mimic;
    // Use the constructor of the base class
    using ISimulator::ISimulator;

public:
    // Only builder can construct
    friend class SimulatorBuilder;

    /*! \brief Function to run the correct SimulatorFunctionType,
     * based on the .mdp integrator field. */
    void run() override;
};

} // namespace gmx

#endif // GMX_MDRUN_LEGACYSIMULATOR_H
