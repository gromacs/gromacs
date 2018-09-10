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
/*! \libinternal \file
 *
 * \brief This file declares the state of an MD run
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_MDSTATE_H
#define GMX_MDLIB_MDSTATE_H

#include <memory>

#include "gromacs/awh/awh.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/imd/imd.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/mdsetup.h"
#include "gromacs/mdlib/repl_ex.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdrun/integrator.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/multireadsinglewritepointer.h"
#include "gromacs/utility/logger.h"

namespace gmx {
    class SimulationSetup {
    public:
        SimulationSetup(Integrator* integrator);

        bool bGStatEveryStep;
        bool bGStat;
        bool bCalcVir;
        bool bCalcEnerStep;
        bool bCalcEner;
        bool bNS;
        bool bNStList;
        bool bSimAnn;
        bool bStopCM;
        bool bFirstStep;
        bool bInitStep;
        bool bLastStep;
        bool bDoDHDL;
        bool bDoFEP;
        bool bDoExpanded;
        bool do_ene;
        bool do_log;
        bool do_verbose;
        bool bRerunWarnNoV;
        bool bForceUpdate;
        bool bMasterState;
        bool bSumEkinhOld;
        bool bDoReplEx;
        bool bExchanged;
        bool bNeedRepartition;
        bool bTemp;
        bool bPres;
        bool bTrotter;
        bool bPMETune;
        bool bPMETunePrinting;
        bool bIMDstep;
        bool bRerunMD;
        bool startingFromCheckpoint;
        bool useReplicaExchange;
        bool simulationsShareState;
        bool resetCountersIsLocal;

        int nstfep;
        int nstglobalcomm;
        int nstSignalComm;

        int64_t step;
        int64_t step_rel;
    };

    class MDState {
    public:
        // gmx_mdoutf       *outf = nullptr;
        multiReadSingleWritePointer<gmx_mdoutf> outf;
        // double            t, t0, lam0[efptNR];
        multiReadSingleWritePointer<double> t, t0;
        double lam0[efptNR];
        tensor force_vir, shake_vir, total_vir, pres;
        rvec mu_tot;
        multiReadSingleWritePointer<t_vcm> vcm;
        matrix parrinellorahmanMu, M;
        multiReadSingleWritePointer<gmx_repl_ex> repl_ex;
        multiReadSingleWritePointer<gmx_localtop_t> top;
        multiReadSingleWritePointer<t_mdebin> mdebin;
        multiReadSingleWritePointer<gmx_enerdata_t> enerd;
        PaddedRVecVector f{};
        multiReadSingleWritePointer<gmx_global_stat> gstat;
        multiReadSingleWritePointer<gmx_update_t> upd;
        multiReadSingleWritePointer<t_graph> graph;
        multiReadSingleWritePointer<gmx_groups_t> groups;
        multiReadSingleWritePointer<gmx_ekindata_t> ekind;
        multiReadSingleWritePointer<gmx_shellfc_t> shellfc;
        multiReadSingleWritePointer<t_extmass> MassQ;
        int **trotter_seq;

        /* PME load balancing data for GPU kernels */
        multiReadSingleWritePointer<pme_load_balancing_t> pme_loadbal;

        multiReadSingleWritePointer<EssentialDynamics> ed;

        std::unique_ptr<t_state> stateInstance;
        multiReadSingleWritePointer<t_state> state;

        multiReadSingleWritePointer<t_mdatoms> mdatoms;

        multiReadSingleWritePointer<Awh> awh;

        multiReadSingleWritePointer<SimulationSignals> signals;

        /* Domain decomposition could incorrectly miss a bonded
           interaction, but checking for that requires a global
           communication stage, which does not otherwise happen in DD
           code. So we do that alongside the first global energy reduction
           after a new DD is made. These variables handle whether the
           check happens, and the result it returns. */
        bool shouldCheckNumberOfBondedInteractions = false;
        int totalNumberOfBondedInteractions = -1;

        // rerun only
        t_trxframe rerun_fr;
        t_trxstatus *status;

        void init(
                Integrator *integrator,
                std::shared_ptr<SimulationSetup> setup);
    };

}  // namespace gmx

#endif //GMX_MDLIB_MDSTATE_H
