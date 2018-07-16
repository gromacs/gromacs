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
 * \brief
 * Implements classes from iforceprovider.h.
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \ingroup module_mdtypes
 */

#include "iforceschedule.h"

using namespace gmx;

/*
 * Base class default constructor
 */
gmx::IForceSchedule::IForceSchedule() :
    log_(nullptr), cr_(nullptr), ms_(nullptr), inputrec_(nullptr), nrnb_(nullptr), wcycle_(nullptr),
    top_(nullptr), groups_(nullptr), state_(nullptr), force_(nullptr), vir_force_(nullptr),
    mdatoms_(nullptr), enerd_(nullptr), fcd_(nullptr), graph_(nullptr), fr_(nullptr),
    vsite_(nullptr), mu_tot_(nullptr), ed_(nullptr), sp_(nullptr), awh_(nullptr), enforcedRotation_(nullptr)
{
    // NULL Declaration
}

/*
 * Base class destructor
 */
gmx::IForceSchedule::~IForceSchedule()
{
    /*
     * Delete stuff owned by the schedule
     */
}

void gmx::IForceSchedule::init(FILE* &fplog, t_commrec* &cr, const gmx_multisim_t* &ms, t_inputrec* &inputrec, t_nrnb* &nrnb,
                               gmx_wallcycle* &wcycle, gmx_localtop_t* &top, gmx_groups_t* &groups,
                               t_state* &state, PaddedRVecVector* &force, tensor* &vir_force,
                               t_mdatoms* &mdatoms, gmx_enerdata_t* &enerd, t_fcdata* &fcd, t_graph* &graph, t_forcerec* &fr,
                               gmx_vsite_t* &vsite, rvec* &mu_tot, gmx_edsam* &ed, MDLoopSharedPrimitives* &sp, Awh* &awh, gmx_enfrot *enforcedRotation)

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
}
