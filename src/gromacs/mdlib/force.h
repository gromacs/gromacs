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
#ifndef GMX_MDLIB_FORCE_H
#define GMX_MDLIB_FORCE_H

#include <cstdint>
#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class DDBalanceRegionHandler;
struct gmx_edsam;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct SimulationGroups;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_wallcycle;
struct gmx_pme_t;
class history_t;
class InteractionDefinitions;
struct pull_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_lambda;
struct t_mdatoms;
struct t_nrnb;
struct gmx_ewald_tab_t;
class CpuPpLongRangeNonbondeds;

namespace gmx
{
template<typename>
class ArrayRefWithPadding;
class Awh;
class ForceBuffersView;
class ForceWithVirial;
class ImdSession;
struct MDModulesNotifiers;
class MdrunScheduleWorkload;
class MDLogger;
class StepWorkload;
class VirtualSitesHandler;
} // namespace gmx

struct ewald_corr_thread_t
{
    real                                                            Vcorr_q;
    real                                                            Vcorr_lj;
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real> dvdl;
    tensor                                                          vir_q;
    tensor                                                          vir_lj;
};

namespace gmx
{

/* Perform the force and, if requested, energy computation
 *
 * Without multiple time stepping the force is returned in force->force().
 *
 * With multiple time stepping the behavior depends on the integration step.
 * At fast steps (step % mtsFactor != 0), the fast force is returned in
 * force->force(). The force->forceMtsCombined() buffer is unused.
 * At slow steps, the normal force is returned in force->force(),
 * unless the \p runScheduleWork.stepWork.useOnlyMtsCombinedForceBuffer==true.
 * A MTS-combined force, F_fast + mtsFactor*F_slow, is always returned in
 * force->forceMtsCombined(). This forceMts can be used directly in a standard
 * leap-frog integrator to do multiple time stepping.
 */
void do_force(FILE*                         log,
              const t_commrec*              cr,
              const gmx_multisim_t*         ms,
              const t_inputrec&             inputrec,
              const MDModulesNotifiers&     mdModulesNotifiers,
              Awh*                          awh,
              gmx_enfrot*                   enforcedRotation,
              ImdSession*                   imdSession,
              pull_t*                       pull_work,
              int64_t                       step,
              t_nrnb*                       nrnb,
              gmx_wallcycle*                wcycle,
              const gmx_localtop_t*         top,
              const matrix                  box,
              ArrayRefWithPadding<RVec>     coordinates,
              ArrayRef<RVec>                velocities,
              const history_t*              hist,
              ForceBuffersView*             force,
              tensor                        vir_force,
              const t_mdatoms*              mdatoms,
              gmx_enerdata_t*               enerd,
              ArrayRef<const real>          lambda,
              t_forcerec*                   fr,
              const MdrunScheduleWorkload&  runScheduleWork,
              VirtualSitesHandler*          vsite,
              rvec                          mu_tot,
              double                        t,
              gmx_edsam*                    ed,
              CpuPpLongRangeNonbondeds*     longRangeNonbondeds,
              const DDBalanceRegionHandler& ddBalanceRegionHandler);

} // namespace gmx

/* Communicate coordinates (if parallel).
 * Do neighbor searching (if necessary).
 * Calculate forces.
 * Communicate forces (if parallel).
 * Spread forces for vsites (if present).
 *
 * f is always required.
 */
class CpuPpLongRangeNonbondeds
{
public:
    /* \brief Constructor
     *
     * Should be called after init_forcerec if params come from a populated forcerec
     */
    CpuPpLongRangeNonbondeds(int                         numberOfTestPaticles,
                             real                        ewaldCoeffQ,
                             real                        epsilonR,
                             gmx::ArrayRef<const double> chargeC6Sum,
                             CoulombInteractionType      eeltype,
                             VanDerWaalsType             vdwtype,
                             const t_inputrec&           inputrec,
                             t_nrnb*                     nrnb,
                             gmx_wallcycle*              wcycle,
                             FILE*                       fplog);

    ~CpuPpLongRangeNonbondeds();

    void updateAfterPartition(const t_mdatoms& md);

    /* Calculate CPU Ewald or PME-mesh forces when done on this rank and Ewald corrections, when used
     *
     * Note that Ewald dipole and net charge corrections are always computed here, independently
     * of whether the PME-mesh contribution is computed on a separate PME rank or on a GPU.
     */
    void calculate(gmx_pme_t*                     pmedata,
                   const t_commrec*               commrec,
                   gmx::ArrayRef<const gmx::RVec> coordinates,
                   gmx::ForceWithVirial*          forceWithVirial,
                   gmx_enerdata_t*                enerd,
                   const matrix                   box,
                   gmx::ArrayRef<const real>      lambda,
                   gmx::ArrayRef<const gmx::RVec> mu_tot,
                   const gmx::StepWorkload&       stepWork,
                   const DDBalanceRegionHandler&  ddBalanceRegionHandler);

private:
    //! Number of particles for test particle insertion
    int numTpiAtoms_;
    //! Ewald charge coefficient
    real ewaldCoeffQ_;
    //! Dielectric constant
    real epsilonR_;
    //! [0]: sum of charges; [1]: sum of C6's
    gmx::ArrayRef<const double> chargeC6Sum_;
    //! Cut-off treatment for Coulomb
    CoulombInteractionType coulombInteractionType_;
    //! Van der Waals interaction treatment
    VanDerWaalsType vanDerWaalsType_;
    //! Ewald geometry
    EwaldGeometry ewaldGeometry_;
    //! Epsilon for PME dipole correction
    real epsilonSurface_;
    //! Whether a long range correction is used
    bool haveEwaldSurfaceTerm_;
    //! Scaling factor for the box for Ewald
    real wallEwaldZfac_;
    //! Whether the simulation is 2D periodic with two walls
    bool havePbcXY2Walls_;
    //! Free energy perturbation type
    FreeEnergyPerturbationType freeEnergyPerturbationType_;
    //! Number of atoms on this node
    int homenr_;
    //! Whether there are perturbed interactions
    bool havePerturbed_;
    //! State A charge
    gmx::ArrayRef<const real> chargeA_;
    //! State B charge
    gmx::ArrayRef<const real> chargeB_;
    //! State A LJ c6
    gmx::ArrayRef<const real> sqrt_c6A_;
    //! State B LJ c6
    gmx::ArrayRef<const real> sqrt_c6B_;
    //! State A LJ sigma
    gmx::ArrayRef<const real> sigmaA_;
    //! State B LJ sigma
    gmx::ArrayRef<const real> sigmaB_;
    //! Ewald correction thread local virial and energy data
    std::vector<ewald_corr_thread_t> outputPerThread_;
    //! Ewald table
    std::unique_ptr<gmx_ewald_tab_t> ewaldTable_;
    //! Non bonded kernel flop counters
    t_nrnb* nrnb_;
    //! Wall cycle counters
    gmx_wallcycle* wcycle_;
};

#endif
