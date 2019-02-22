/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal \defgroup module_forceschedules Schedules used for calculating forces.
 * \ingroup group_mdrun
 *
 * \brief Implements families of force-calculation schedules that will
 * cater to special cases of simulation inputs, hardware contexts, or
 * degrees of parallelism. The IForceSchedule seam that it contains
 * is likely to be a key seam required for the NB-LIB project to use.
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 */

/*! \libinternal \file
 *
 * \brief Describes the Schedule abstraction for Force Calculations
 * in each MD step. Basically allows one to hand-craft a sequence of
 * operations to optimally overlap computation, communication and
 * reduction tasks.
 *
 * Concrete implementations of the Schedule interface can be specific
 * to the algorithm, hardware architecture, quantities of interest, etc.
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 *
 * \ingroup module_forceschedules
 * \inlibraryapi
 */

#ifndef GMX_MDTYPES_IFORCESCHEDULE_H
#define GMX_MDTYPES_IFORCESCHEDULE_H

#include <cstdio>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"

class DDBalanceRegionHandler;
class t_state;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_graph;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
struct gmx_enerdata_t;
struct gmx_edsam;
struct gmx_enfrot;
struct gmx_groups_t;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_wallcycle;
struct gmx_vsite_t;

namespace gmx
{
class Awh;
class PpForceWorkload;

/*! \libinternal \brief Struct with the values that can change during the simulation
 *
 * \todo Some users of IForceSchedule will not have the notion of step or time,
 * so we should refactor the force calculation modules so that modules that
 * need this handle their own details. If the use is just for logging, then
 * we should add that to the logger rather than propagate it through the
 * interface of the force calculation calls. */
struct MDLoopSharedPrimitives
{
    //! MD step number
    int64_t *step_;
    //! MD time
    double  *t_;

    //! Initializing with references to vars
    MDLoopSharedPrimitives(int64_t *step, double *t) :
        step_(step), t_(t)
    {}
};

/*! \libinternal \brief
 * Interface for programming static force calculation schedules
 *
 * Allows one to hand-craft sequences of operations in specific cases
 * of force calculations. The objective is to build a collection of
 * such schedules that optimally overlap computation, communication
 * (to accelerators or message passing in the network) and reduction steps.
 *
 * These schedules are envisioned to be tailor made for hardware. network/NUMA
 * awareness, quantities of interest, physics at play and so on. Different
 * concrete implementations would be utilized depending on the workloads of a
 * given simulation/rank.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IForceSchedule
{

    protected:
        //! Handles logging.
        FILE                           *log_;
        //! Handles communication.
        t_commrec                      *cr_;
        //! Coordinates multi-simulations.
        const gmx_multisim_t           *ms_;
        //! Contains user input mdp options.
        t_inputrec                     *inputrec_;

        //! Manages flop accounting.
        t_nrnb                         *nrnb_;
        //! Manages wall cycle accounting.
        gmx_wallcycle                  *wcycle_;
        /*! \brief Contains topology info for a domain over its
         * lifetime. */
        gmx_localtop_t                 *top_;

        //! Stores data and names of Groups
        gmx_groups_t                   *groups_;
        //! Contains the microstate of the system.
        t_state                        *state_;

        //! PaddedVector of  RVecs storing the forces.
        PaddedVector<gmx::RVec>        *force_;
        //! Stores Force Virials
        tensor                         *vir_force_;
        //! Atom parameters for this domain.
        t_mdatoms                      *mdatoms_;
        //! Contains energy data of various interactions
        gmx_enerdata_t                 *enerd_;
        //! Helper struct for force calculations.
        t_fcdata                       *fcd_;

        /*! \brief Describes the connectivity between the internally
         * chemically bonded parts. */
        t_graph                        *graph_;
        //! Parameters for force calculations.
        t_forcerec                     *fr_;
        //! Handles virtual sites.
        gmx_vsite_t                    *vsite_;
        //! System dipole moment
        rvec                           *mu_tot_;

        //! Essential dynamics
        gmx_edsam                      *ed_;
        //! Variables changing over the MD loop
        MDLoopSharedPrimitives         *sp_;
        //! Accelerated Weight Histogram
        Awh                            *awh_;

        //! Handles enforced rotation.
        gmx_enfrot                     *enforcedRotation_;
        //! Aggregates force-calculation work flags
        PpForceWorkload                *ppForceWorkload_;
        //! Handles DD balancing.
        const DDBalanceRegionHandler   *ddBalanceRegionHandler_;

    public:
        //! Initializes a schedule with the contextual arguments
        IForceSchedule(FILE* fplog, t_commrec* cr, const gmx_multisim_t* ms, t_inputrec* inputrec,
                       t_nrnb* nrnb, gmx_wallcycle* wcycle, gmx_localtop_t* top, gmx_groups_t* groups,
                       t_state* state, PaddedVector<gmx::RVec>* force, tensor* vir_force, t_mdatoms* mdatoms,
                       gmx_enerdata_t* enerd, t_fcdata* fcd, t_graph* graph, t_forcerec* fr,
                       gmx_vsite_t* vsite, rvec* mu_tot, gmx_edsam* ed, MDLoopSharedPrimitives* sp,
                       Awh* awh, gmx_enfrot *enforcedRotation, PpForceWorkload *ppForceWorkload,
                       const DDBalanceRegionHandler *ddBalanceRegionHandler) :
            log_(fplog),
            cr_(cr),
            ms_(ms),
            inputrec_(inputrec),
            nrnb_(nrnb),
            wcycle_(wcycle),
            top_(top),
            groups_(groups),
            state_(state),
            force_(force),
            vir_force_(vir_force),
            mdatoms_(mdatoms),
            enerd_(enerd),
            fcd_(fcd),
            graph_(graph),
            fr_(fr),
            vsite_(vsite),
            mu_tot_(mu_tot),
            ed_(ed),
            sp_(sp),
            awh_(awh),
            enforcedRotation_(enforcedRotation),
            ppForceWorkload_(ppForceWorkload),
            ddBalanceRegionHandler_(ddBalanceRegionHandler)
        {}

        virtual ~IForceSchedule() = default;

        //! Computes forces and other relevant quantities for a single MD step
        virtual void computeForces(int flags) = 0;

        friend class ScheduleBuilder;

};

}      // namespace gmx

#endif //GMX_FORCESCHEDULES_IFORCESCHEDULE_H
