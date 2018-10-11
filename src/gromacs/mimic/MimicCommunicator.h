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

#ifndef GMX_MIMIC_COMMUNICATOR_H
#define GMX_MIMIC_COMMUNICATOR_H

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"

namespace gmx
{
/**
 * \inlibraryapi
 * \internal \brief
 * Class-wrapper around MiMiC communication library
 * It uses GROMACS' unit conversion to switch from GROMACS' units to a.u.
 *
 * \author Viacheslav Bolnykh <v.bolnykh@hpc-leap.eu>
 * \ingroup module_mimic
 */
class MimicCommunicator
{

    public:
        /*! \brief
         * Initializes the communicator
         */
        void init();

        /*! \brief
         * Sends the data needed for MiMiC initialization
         *
         * That includes number of atoms, element numbers, charges, masses,
         * maximal order of multipoles (0 for point-charge forcefields),
         * number of molecules, number of atoms per each molecule,
         * bond constraints data
         *
         * @param mtop global topology data
         * @param coords coordinates of all atoms
         */
        void sendInitData(gmx_mtop_t                *mtop,
                          HostVector<gmx::RVec>      coords);

        /*! \brief
         * Gets the number of MD steps to perform from MiMiC
         *
         * @return nsteps the number of MD steps to perform
         */
        int64_t getStepNumber();

        /*! \brief
         * Receive and array of updated atomic coordinates from MiMiC
         *
         * @param x array of coordinates to fill
         * @param natoms number of atoms in the system
         */
        void getCoords(HostVector<RVec> *x, int natoms);

        /*! \brief
         * Send the potential energy value to MiMiC
         *
         * @param energy energy value to send
         */
        void sendEnergies(real energy);

        /*! \brief
         * Send classical forces acting on all atoms in the system
         * to MiMiC.
         *
         * @param forces array of forces to send
         * @param natoms number of atoms in the system
         */
        void sendForces(ArrayRef<gmx::RVec> forces, int natoms);

        /*! \brief
         * Finish communications and disconnect from the server
         */
        void finalize();

};

}      // namespace gmx

#endif //GMX_MIMIC_COMMUNICATOR_H
