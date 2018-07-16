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
#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/inputrec.h"
#include "md_enums.h"
#include <memory>
#include <gromacs/mdlib/nb_verlet.h>
#include "gromacs/utility/fatalerror.h"
#include "commrec.h"
#include "forcerec.h"

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

struct MDLoopSharedPrimitives {
    int64_t                                   *step_;
    double                                    *t_; //time

    DdOpenBalanceRegionBeforeForceComputation *ddOpenBalanceRegion_;
    DdCloseBalanceRegionAfterForceComputation *ddCloseBalanceRegion_;

    MDLoopSharedPrimitives(int64_t *step, double *t,
                           DdOpenBalanceRegionBeforeForceComputation *ddOpenBalanceRegion,
                           DdCloseBalanceRegionAfterForceComputation *ddCloseBalanceRegion) :

        step_(step), t_(t), ddOpenBalanceRegion_(ddOpenBalanceRegion),
        ddCloseBalanceRegion_(ddCloseBalanceRegion)
    {
        // Initialize the struct with the values that can change during the simulation
    }
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
        FILE                           *log_;
        t_commrec                      *cr_;
        const gmx_multisim_t           *ms_;
        t_inputrec                     *inputrec_;

        t_nrnb                         *nrnb_;
        gmx_wallcycle                  *wcycle_;
        gmx_localtop_t                 *top_;
        gmx_groups_t                   *groups_;
        t_state                        *state_;

        PaddedArrayRef<RVec>            force_;
        tensor                         *vir_force_;
        t_mdatoms                      *mdatoms_;
        gmx_enerdata_t                 *enerd_;
        t_fcdata                       *fcd_;

        t_graph                        *graph_;
        t_forcerec                     *fr_;
        gmx_vsite_t                    *vsite_;
        rvec                           *mu_tot_;

        gmx_edsam                      *ed_;
        MDLoopSharedPrimitives         *sp_;
        Awh                            *awh_;

        gmx_enfrot                     *enforcedRotation_;

    protected:

        IForceSchedule();

        virtual ~IForceSchedule() = 0;

        void init(FILE* fplog, t_commrec* cr, const gmx_multisim_t* ms, t_inputrec* inputrec,
                  t_nrnb* nrnb, gmx_wallcycle* wcycle, gmx_localtop_t* top, gmx_groups_t* groups,
                  t_state* state, PaddedArrayRef<RVec> force, tensor* vir_force, t_mdatoms* mdatoms,
                  gmx_enerdata_t* enerd, t_fcdata* fcd, t_graph* graph, t_forcerec* fr,
                  gmx_vsite_t* vsite, rvec* mu_tot, gmx_edsam* ed, MDLoopSharedPrimitives* sp,
                  Awh* awh, gmx_enfrot *enforcedRotation);

    public:

        virtual void computeStep(int flags) = 0;

        friend class ScheduleBuilder;

};

class DefaultVerletSchedule : public IForceSchedule
{

    public:

        void computeStep(int flags) override;

};

class DefaultGroupSchedule : public IForceSchedule
{

    public:

        void computeStep(int flags) override;

};

class SRCPUSingleNodeSchedule : public IForceSchedule
{

    public:

        void computeStep(int flags) override;
};

class DummySchedule : public IForceSchedule
{

    public:

        void computeStep(int flags) override;

};

/*! \brief Concrete Implementations of Force Schedules
 *
 * For use in functions that can change the schedule type
 * during the simulation
 */

enum ScheduleType
{
    Dummy, Verlet, Group, SRCPU_verlet
};

/*
 * TODO: Implement a schedule selector function. Also find a way to wrap the
 *       execution context and problem details into nice, readable objects
 */

/*! \libinternal \brief
 * Class to encapsulate the initialization of concrete schedules
 *
 * The class has the responsibility of initializing schedules either
 * by using a selector function to process the execution context or
 * manually specifying an identifier for the concrete schedules.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

class ScheduleBuilder
{

    private:
        std::shared_ptr<IForceSchedule> root_;
        ScheduleType                    type_;

        void setSchedule()
        {
            if (type_ == Verlet)
            {
                root_ = std::make_shared<DefaultVerletSchedule>();
            }
            else if (type_ == Group)
            {
                root_ = std::make_shared<DefaultGroupSchedule>();
            }
            else if (type_ == Dummy)
            {
                root_ = std::make_shared<DummySchedule>();
            }
            else if (type_ == SRCPU_verlet)
            {
                root_ = std::make_shared<SRCPUSingleNodeSchedule>();
            }
        }

    public:

        ScheduleBuilder() = default;

        void selectSchedule(t_inputrec *inputrec, t_commrec *cr, t_forcerec *fr
                            /* execution context and problem definition */)
        {
            if (inputrec->coulombtype == eelCUT && !EVDW_PME(inputrec->vdwtype) && !DOMAINDECOMP(cr) && !(fr->nbv->bUseGPU))
            {
                type_ = SRCPU_verlet;
            }
            else
            {
                switch (inputrec->cutoff_scheme)
                {
                    case ecutsVERLET:
                        type_ = Verlet;
                        break;
                    case ecutsGROUP:
                        type_ = Group;
                        break;
                    default:
                        gmx_incons("Invalid cut-off scheme passed!");
                }
            }

            setSchedule();

        }

        void setScheduleType(ScheduleType type)
        {
            type_ = type;
            setSchedule();
        }

        void initSchedule(FILE *fplog, t_commrec *cr, const gmx_multisim_t *ms, t_inputrec *inputrec, t_nrnb *nrnb,
                          gmx_wallcycle *wcycle, gmx_localtop_t *top, gmx_groups_t *groups,
                          t_state *state, PaddedArrayRef<RVec> force, tensor *vir_force,
                          t_mdatoms *mdatoms, gmx_enerdata_t *enerd, t_fcdata *fcd, t_graph *graph, t_forcerec *fr,
                          gmx_vsite_t *vsite, rvec *mu_tot, gmx_edsam *ed, MDLoopSharedPrimitives *sp, Awh *awh, gmx_enfrot *enforcedRotation)
        {
            if (!root_)
            {
                selectSchedule(inputrec, cr, fr);
            }
            root_->init(fplog, cr, ms, inputrec, nrnb, wcycle, top, groups, state, force,
                        vir_force, mdatoms, enerd, fcd, graph, fr, vsite, mu_tot, ed, sp, awh, enforcedRotation);
        }

        std::shared_ptr<IForceSchedule> getSchedule()
        {
            return root_;
        }

};

}      // namespace gmx

#endif //GMX_MDTYPES_IFORCESCHEDULE_H
