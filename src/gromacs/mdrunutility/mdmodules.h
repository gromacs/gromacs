/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::MDModules.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_MDMODULES_H
#define GMX_MDRUNUTILITY_MDMODULES_H

#include "gromacs/utility/classhelpers.h"

struct t_inputrec;

namespace gmx
{

/*! \libinternal \brief
 * Factory for t_inputrec.
 *
 * This class acts as a central place for constructing t_inputrec (and possibly
 * other mdrun structures in the future) and wiring up dependencies between
 * modules that are referenced from these structures.
 *
 * The general idea is that each module takes care of its own data rather than
 * mdrun having to know about all the details of each type of force calculation.
 * Initially this is applied for simple things like electric field calculations
 * but later more complex forces will be supported too.
 *
 * In the simplest form, this module delivers the inputrec. The inputrec holds the
 * information on modules inside (efield for now only). The efield (being an
 * IInputRecExtension) has a method to initialize a forcerec which in turn has an
 * IForceProvider that actually computes forces and energies.
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
class MDModules
{
    public:
        MDModules();
        ~MDModules();

        /*! \brief
         * Returns an initialized t_inputrec structure.
         *
         * The inputrec structure is owned by MDModules and will be destroyed
         * with it.
         */
        t_inputrec *inputrec();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
