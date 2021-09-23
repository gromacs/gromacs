/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_VCM_H
#define GMX_MDLIB_VCM_H

#include <stdio.h>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/real.h"

struct SimulationGroups;
struct t_inputrec;
struct t_mdatoms;

namespace gmx
{
template<typename T>
class ArrayRef;
}

struct t_vcm_thread
{
    //! Linear momentum
    rvec p = { 0 };
    //! Center of mass
    rvec x = { 0 };
    //! Angular momentum
    rvec j = { 0 };
    //! Moment of inertia
    tensor i = { { 0 } };
    //! Mass
    real mass = 0;
};

struct t_vcm
{
    //! Number of groups
    int nr = 0;
    //! Size of group arrays
    int size = 0;
    //! Stride for thread data
    int stride = 0;
    //! One of the enums above
    ComRemovalAlgorithm mode = ComRemovalAlgorithm::Linear;
    //! The number of dimensions for corr.
    int ndim = 0;
    //! The time step for COMM removal
    real timeStep = 0;
    //! Number of degrees of freedom
    std::vector<real> group_ndf;
    //! Mass per group
    std::vector<real> group_mass;
    //! Linear momentum per group
    std::vector<gmx::RVec> group_p;
    //! Linear velocity per group
    std::vector<gmx::RVec> group_v;
    //! Center of mass per group
    std::vector<gmx::RVec> group_x;
    //! Angular momentum per group
    std::vector<gmx::RVec> group_j;
    //! Angular velocity (omega)
    std::vector<gmx::RVec> group_w;
    //! Moment of inertia per group
    tensor* group_i = nullptr;
    //! These two are copies to pointers in
    std::vector<char*> group_name;
    //! Tells whether dimensions are frozen per freeze group
    ivec* nFreeze = nullptr;
    //! Temporary data per thread and group
    std::vector<t_vcm_thread> thread_vcm;

    //! Tell whether the integrator conserves momentum
    bool integratorConservesMomentum = false;

    t_vcm(const SimulationGroups& groups, const t_inputrec& ir);
    ~t_vcm();
};

/* print COM removal info to log */
void reportComRemovalInfo(FILE* fp, const t_vcm& vcm);


/* Do a per group center of mass things */
void calc_vcm_grp(const t_mdatoms&               md,
                  gmx::ArrayRef<const gmx::RVec> x,
                  gmx::ArrayRef<const gmx::RVec> v,
                  t_vcm*                         vcm);

/* Set the COM velocity to zero and potentially correct the COM position.
 *
 * Processes the kinetic energy reduced over MPI before removing COM motion.
 * With mode linear, nullptr can be passed for x.
 * With acceleration correction nullptr should be passed for x at initialization
 * and a pointer to the coordinates at normal MD steps.
 * When fplog != nullptr, a warning is printed to fplog with high COM velocity.
 */
void process_and_stopcm_grp(FILE*                    fplog,
                            t_vcm*                   vcm,
                            const t_mdatoms&         mdatoms,
                            gmx::ArrayRef<gmx::RVec> x,
                            gmx::ArrayRef<gmx::RVec> v);

#endif
