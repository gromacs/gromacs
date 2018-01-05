/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#ifndef GROMACS_MIMIC_COMMUNICATOR_H
#define GROMACS_MIMIC_COMMUNICATOR_H

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
 * \ingroup module_mdlib
 */
class MimicCommunicator
{

    private:
        /**
         * Conversion factor of nm to Bohr
         */
        static constexpr double bohr = 0.0529177210859;

        /**
         * Conversion factor of kJ to Hartree
         */
        static constexpr double hartree = 0.00038;

    public:
        /**
         * Initializes the communicator
         */
        void init();

        /**
         * \brief
         * Sends the data needed for MiMiC initialization
         * @param mtop global topology data
         * @param coords coordinates of all atoms
         */
        void sendInitData(gmx_mtop_t                *mtop,
                          gmx::HostVector<gmx::RVec> coords,
                          gmx::Constraints          *constr);

        /**
         * \brief
         * Gets the number of MD steps
         * @param nsteps the number of MD steps to perform
         */
        void getStepNumber(gmx_int64_t &nsteps);

        /**
         * \brief
         * Fills the array with updated atomic coordinates
         * @param x array of coordinates to fill
         * @param natoms number of atoms in the system
         */
        void getCoords(rvec *x, int natoms);

        /**
         * \brief
         * Send the energy value to MiMiC
         * @param energy energy value to send
         */
        void sendEnergies(real energy);

        /**
         * \brief
         * Send computed forces to MiMiC
         * @param forces array of forces to send
         * @param natoms number of atoms in the system
         */
        void sendForces(gmx::ArrayRef<gmx::RVec> forces, int natoms);

        /**
         * \brief
         * Finish communications and disconnect from the server
         */
        void finalize();

};

}

#endif //GROMACS_MIMIC_COMMUNICATOR_H
