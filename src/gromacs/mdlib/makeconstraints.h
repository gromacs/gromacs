/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief Declares and implements factory function for Constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_MAKECONSTRAINTS_H
#define GMX_MDLIB_MAKECONSTRAINTS_H

#include <memory>

#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"

struct gmx_mtop_t;

namespace gmx
{

/*! \internal
 * \brief Support type to help implement makeConstraints().
 *
 * This member type of Constraints also inherits from it, so that it
 * can access the private constructor of Constraints to support the
 * implementation of the factory function. This approach avoids having
 * to declare makeConstraints() as a template friend function. */
struct Constraints::CreationHelper : public Constraints
{
public:
    /*! \brief Constructor that can call the private constructor
     * of Constraints.
     *
     * The parameter pack insulates this helper type from changes
     * to the arguments to the constructor. */
    template<typename... Args>
    CreationHelper(Args&&... args) : Constraints(std::forward<Args>(args)...)
    {
    }
};

/*! \brief Factory function for Constraints.
 *
 * We only want an object to manage computing constraints when the
 * simulation requires one. Checking for whether the object was made
 * adds overhead to simulations that use constraints, while avoiding
 * overhead on those that do not, so is a design trade-off we might
 * reconsider some time.
 *
 * Using a private constructor and a factory function ensures that we
 * can only make a Constraints object when the prerequisites are
 * satisfied, ie. that something needs them and if necessary has
 * already been initialized.
 *
 * Using the parameter pack insulates the factory function from
 * changes to the type signature of the constructor that don't
 * affect the logic here. */
template<typename... Args>
std::unique_ptr<Constraints> makeConstraints(const gmx_mtop_t& mtop,
                                             const t_inputrec& ir,
                                             pull_t*           pull_work,
                                             bool              havePullConstraintsWork,
                                             bool              doEssentialDynamics,
                                             Args&&... args)
{
    int numConstraints = (gmx_mtop_ftype_count(mtop, F_CONSTR) + gmx_mtop_ftype_count(mtop, F_CONSTRNC));
    int numSettles = gmx_mtop_ftype_count(mtop, F_SETTLE);
    GMX_RELEASE_ASSERT(!ir.bPull || pull_work != nullptr,
                       "When COM pulling is active, it must be initialized before constraints are "
                       "initialized");
    bool doPullingWithConstraints = ir.bPull && havePullConstraintsWork;
    if (numConstraints + numSettles == 0 && !doPullingWithConstraints && !doEssentialDynamics)
    {
        // No work, so don't make a Constraints object.
        return nullptr;
    }
    return std::make_unique<Constraints::CreationHelper>(
            mtop, ir, pull_work, std::forward<Args>(args)..., numConstraints, numSettles);
}

} // namespace gmx

#endif
