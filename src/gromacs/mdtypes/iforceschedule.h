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
/*! \internal \file
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
 * \inlibraryapi
 * \ingroup module_mdlib, gmxlib
 */

#ifndef GMX_MDTYPES_IFORCESCHEDULE_H
#define GMX_MDTYPES_IFORCESCHEDULE_H

#include "gromacs/math/vectypes.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"
#include "commrec.h"
#include "forcerec.h"
#include "md_enums.h"

enum class DdOpenBalanceRegionBeforeForceComputation;
enum class DdCloseBalanceRegionAfterForceComputation;

class t_state;
struct t_forcerec;
struct gmx_groups_t;
struct gmx_localtop_t;
struct t_mdatoms;
struct t_graph;
struct gmx_enerdata_t;
struct gmx_edsam;
struct gmx_multisim_t;
struct t_commrec;
struct gmx_wallcycle;
struct t_fcdata;
struct gmx_vsite_t;
struct t_nrnb;

namespace gmx
{
class Awh;
class PpForceWorkload;

//! Struct with the values that can change during the simulation

struct MDLoopSharedPrimitives
{
    int64_t                                   *step_;
    double                                    *t_; //time

    DdOpenBalanceRegionBeforeForceComputation *ddOpenBalanceRegion_;
    DdCloseBalanceRegionAfterForceComputation *ddCloseBalanceRegion_;

    //! Initializing with references to vars

    MDLoopSharedPrimitives(int64_t *step, double *t,
                           DdOpenBalanceRegionBeforeForceComputation *ddOpenBalanceRegion,
                           DdCloseBalanceRegionAfterForceComputation *ddCloseBalanceRegion) :

        step_(step), t_(t), ddOpenBalanceRegion_(ddOpenBalanceRegion),
        ddCloseBalanceRegion_(ddCloseBalanceRegion)
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
        // Members
        FILE                           *log_;                       //! Handles logging.
        t_commrec                      *cr_;                        //! Handles communication.
        const gmx_multisim_t           *ms_;                        //! Coordinates multi-simulations.
        t_inputrec                     *inputrec_;                  //! Contains user input mdp options.

        t_nrnb                         *nrnb_;                      //! Manages flop accounting.
        gmx_wallcycle                  *wcycle_;                    //! Manages wall cycle accounting.
        gmx_localtop_t                 *top_;                       //! Contains topology info for a
                                                                    //! domain over its lifetime.
        gmx_groups_t                   *groups_;                    //! Stores data and names of Groups
        t_state                        *state_;                     //! Contains the microstate of the system.

        PaddedVector<gmx::RVec>        *force_;                     //! PaddedVector of  RVecs storing the forces.
        tensor                         *vir_force_;                 //! Stores Force Virials
        t_mdatoms                      *mdatoms_;                   //! Atom parameters for this domain.
        gmx_enerdata_t                 *enerd_;                     //! Contains energy data of various interactions
        t_fcdata                       *fcd_;                       //! Helper struct for force calculations.

        t_graph                        *graph_;                     //! Describes the connectivity between the
                                                                    //! internally chemically bonded parts.
        t_forcerec                     *fr_;                        //! Parameters for force calculations.
        gmx_vsite_t                    *vsite_;                     //! Handles virtual sites.
        rvec                           *mu_tot_;                    //! System dipole moment

        gmx_edsam                      *ed_;                        //! Essential dynamics
        MDLoopSharedPrimitives         *sp_;                        //! Variables changing over the MD loop
        Awh                            *awh_;                       //! Accelerated Weight Histogram

        gmx_enfrot                     *enforcedRotation_;          //! Handles enforced rotation.
        PpForceWorkload                *ppForceWorkload_;           //! Aggregates force-calculation work flags

    protected:

        explicit IForceSchedule() :
            log_(nullptr), cr_(nullptr), ms_(nullptr), inputrec_(nullptr), nrnb_(nullptr),
            wcycle_(nullptr), top_(nullptr), groups_(nullptr), state_(nullptr), force_(nullptr),
            vir_force_(nullptr), mdatoms_(nullptr), enerd_(nullptr), fcd_(nullptr), graph_(nullptr),
            fr_(nullptr), vsite_(nullptr), mu_tot_(nullptr), ed_(nullptr), sp_(nullptr),
            awh_(nullptr), enforcedRotation_(nullptr), ppForceWorkload_(nullptr)
        {}

    public:
        //! Initializes a schedule with the contextual arguments
        explicit IForceSchedule(FILE* fplog, t_commrec* cr, const gmx_multisim_t* ms, t_inputrec* inputrec,
                                t_nrnb* nrnb, gmx_wallcycle* wcycle, gmx_localtop_t* top, gmx_groups_t* groups,
                                t_state* state, PaddedVector<gmx::RVec>* force, tensor* vir_force, t_mdatoms* mdatoms,
                                gmx_enerdata_t* enerd, t_fcdata* fcd, t_graph* graph, t_forcerec* fr,
                                gmx_vsite_t* vsite, rvec* mu_tot, gmx_edsam* ed, MDLoopSharedPrimitives* sp,
                                Awh* awh, gmx_enfrot *enforcedRotation, PpForceWorkload *ppForceWorkload)
        {
            log_              = fplog;
            cr_               = cr;
            ms_               = ms;
            inputrec_         = inputrec;
            nrnb_             = nrnb;
            wcycle_           = wcycle;
            top_              = top;
            groups_           = groups;
            state_            = state;
            force_            = force;
            vir_force_        = vir_force;
            mdatoms_          = mdatoms;
            enerd_            = enerd;
            fcd_              = fcd;
            graph_            = graph;
            fr_               = fr;
            vsite_            = vsite;
            mu_tot_           = mu_tot;
            ed_               = ed;
            sp_               = sp;
            awh_              = awh;
            enforcedRotation_ = enforcedRotation;
            ppForceWorkload_  = ppForceWorkload;
        }

        virtual ~IForceSchedule() = default;

        //! Computes forces and other relevant quantities for a single MD step
        virtual void computeForces(int flags) = 0;

        friend class ScheduleBuilder;

};

/*! \brief Concrete Implementations of Force Schedules
 *
 * For use in functions that can change the schedule type
 * during the simulation
 */

enum ScheduleType
{
    Verlet, Group, None
};

}      // namespace gmx

#endif //GMX_MDTYPES_IFORCESCHEDULE_H
