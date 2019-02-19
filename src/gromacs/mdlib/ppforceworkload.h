/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Declares force calculation workload manager.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_PPFORCEWORKLOAD_H
#define GMX_MDLIB_PPFORCEWORKLOAD_H

namespace gmx
{

/*! \libinternal
 * \brief Manage what force calculation work is required each step.
 *
 * An object of this type is updated every neighbour search stage to
 * reflect what work is required during normal MD steps, e.g. whether
 * there are bonded interactions in this PP task.
 *
 * This will remove the desire for inline getters from modules that
 * describe whether they have work to do, because that can be set up
 * once per simulation or neighborlist lifetime and not changed
 * thereafter.
 *
 * \todo Add more responsibilities, including whether GPUs are in use,
 * whether there is PME work, whether DD is active, whether NB
 * local/nonlocal regions have work, whether forces/virial/energy are
 * required. */
class PpForceWorkload
{
    public:
        //! Whether this MD step has bonded work to run on a GPU.
        bool haveGpuBondedWork = false;
        //! Whether this MD step has bonded work to run on he CPU.
        bool haveCpuBondedWork = false;
        //! Whether this MD step has listed forces bonded work to run on he CPU.
        bool haveRestraintsWork = false;
        //! Whether this MD step has listed forces work to run on he CPU.
        //  Note: currently this is haveCpuBondedWork | haveRestraintsWork
        bool haveCpuListedForceWork = false;
        //! Whether this MD step has special forces on the CPU.
        bool haveSpecialForces = false;
};

} // namespace gmx

#endif
