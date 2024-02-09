/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
#ifndef GMX_EWALD_PME_SOLVE_H
#define GMX_EWALD_PME_SOLVE_H

#include <memory>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_pme_t;
struct PmeAndFftGrids;
struct pme_solve_work_t;
struct PmeOutput;

namespace gmx
{
template<typename>
class ArrayRef;
}

//! Class for solving PME for Coulomb and LJ
class PmeSolve
{
public:
    /*! Constructor
     *
     * \param numThreads  The number of OpenMP threads used during solve
     * \param nkx         The number of PME grid points along dimension X
     */
    PmeSolve(int numThreads, int nkx);

    ~PmeSolve();

    /*! \brief Solves PME for Coulomb
     *
     * \returns the number of grid elements solved
     */
    int solveCoulombYZX(const gmx_pme_t& pme, t_complex* grid, real vol, bool computeEnergyAndVirial, int thread);

    /*! \brief Solves PME for LJ
     *
     * \returns the number of grid elements solved
     */
    int solveLJYZX(const gmx_pme_t&              pme,
                   gmx::ArrayRef<PmeAndFftGrids> grids,
                   bool                          useLBCombinationRule,
                   real                          vol,
                   bool                          computeEnergyAndVirial,
                   int                           thread);

    //! Get Coulomb energy and virial
    void getCoulombEnergyAndVirial(PmeOutput* output) const;

    //! Get LJ energy and virial
    void getLJEnergyAndVirial(PmeOutput* output) const;

private:
    //! Returns the number of threads used for solve
    int numThreads() const { return gmx::ssize(workData_); }

    //! Returns the work data for thread \p thread
    pme_solve_work_t& workData(int thread) { return *workData_[thread]; }

    //! Returns the work data for thread \p thread
    const pme_solve_work_t& workData(int thread) const { return *workData_[thread]; }

    //! Work data for the threads, stored with unique_ptr for thread-local memory
    std::vector<std::unique_ptr<pme_solve_work_t>> workData_;
};

#endif
