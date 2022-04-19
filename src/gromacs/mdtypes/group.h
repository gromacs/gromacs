/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#ifndef GMX_MDTYPES_GROUP_H
#define GMX_MDTYPES_GROUP_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct t_grp_tcstat
{
    //! Temperature at half step
    real Th = 0;
    //! Temperature at full step
    real T = 0;
    //! Kinetic energy at half step
    tensor ekinh = { { 0 } };
    //! Kinetic energy at old half step
    tensor ekinh_old = { { 0 } };
    //! Kinetic energy at full step
    tensor ekinf = { { 0 } };
    //! Berendsen coupling lambda
    real lambda = 0;
    //! Scaling factor for NHC- full step
    double ekinscalef_nhc = 0;
    //! Scaling factor for NHC- half step
    double ekinscaleh_nhc = 0;
    //! Scaling factor for NHC- velocity
    double vscale_nhc = 0;
};

struct t_cos_acc
{
    //! The acceleration for the cosine profile
    real cos_accel = 0;
    //! The cos momenta of home particles
    real mvcos = 0;
    //! The velocity of the cosine profile
    real vcos = 0;
};

class gmx_ekindata_t
{
public:
    gmx_ekindata_t(int numTempCoupleGroups, real cos_accel, int numThreads);
    //! The number of T-coupling groups
    int ngtc;
    //! T-coupling data
    std::vector<t_grp_tcstat> tcstat;
    //! Allocated locations for *_work members
    tensor** ekin_work_alloc = nullptr;
    //! Work arrays for tcstat per thread
    tensor** ekin_work = nullptr;
    //! Work location for dekindl per thread
    real** dekindl_work = nullptr;
    //! overall kinetic energy
    tensor ekin = { { 0 } };
    //! overall 1/2 step kinetic energy
    tensor ekinh = { { 0 } };
    //! dEkin/dlambda at half step
    real dekindl = 0;
    //! dEkin/dlambda at old half step
    real dekindl_old = 0;
    //! Cosine acceleration data
    t_cos_acc cosacc;

    ~gmx_ekindata_t();

private:
    //! For size of ekin_work
    int nthreads_ = 0;
};

#define GID(igid, jgid, gnr) \
    (((igid) < (jgid)) ? ((igid) * (gnr) + (jgid)) : ((jgid) * (gnr) + (igid)))

#endif
