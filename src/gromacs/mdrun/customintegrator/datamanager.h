//
// Created by Pascal Merz on 4/21/18.
//

/*

   Extend state to be aware of time. Allow elements to tell state
   when next state is needed. Save this state if needed.
   Notifications? (Hey, state is at t, need it?)
   Handle whether state is constrained or not!

   Constraints need previous position (could take care of it by itself?)
   Trajectory writing needs specific state
   EDR writing needs specific state
   Checkpointing needs specific state
   Updates / force need whatever is available?
   MC needs valid state.

   Maybe step setup can announce?

   (EDR output and compute globals, loop element)

   Some logging module that gets input from elements? Writes when called?

 */

#ifndef GROMACS_DATAMANAGER_H
#define GROMACS_DATAMANAGER_H

#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/vcm.h"

#include "gromacs/mdlib/repl_ex.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdlib/mdrun.h"

#include "gromacs/mdrun/integrator.h"



namespace gmx
{
class DataManager
{
    public:
        DataManager(Integrator &integrator_ref);
        ~DataManager();

        void preStep();
        void postStep();

        Integrator       &integrator;

        gmx_mdoutf       *outf = nullptr;
        gmx_int64_t       step, step_rel;
        double            elapsed_time;
        double            t, t0, lam0[efptNR];
        gmx_bool          bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
        gmx_bool          bNS, bNStList, bSimAnn, bStopCM,
                          bFirstStep, bInitStep, bLastStep = FALSE;
        gmx_bool          bDoDHDL = FALSE, bDoFEP = FALSE, bDoExpanded = FALSE;
        gmx_bool          do_ene, do_log, do_verbose,
                          bForceUpdate = FALSE, bCPT;
        gmx_bool          bMasterState;
        int               force_flags, cglo_flags;
        tensor            force_vir, shake_vir, total_vir, tmp_vir, pres;
        int               i, m;
        rvec              mu_tot;
        t_vcm            *vcm;
        matrix            parrinellorahmanMu, M;
        gmx_repl_ex_t     repl_ex = nullptr;
        int               nchkpt  = 1;
        gmx_localtop_t   *top;
        t_mdebin         *mdebin   = nullptr;
        gmx_enerdata_t   *enerd;
        PaddedRVecVector  f {};
        gmx_global_stat_t gstat;
        gmx_update_t     *upd   = nullptr;
        t_graph          *graph = nullptr;
        gmx_groups_t     *groups;
        gmx_ekindata_t   *ekind;
        gmx_shellfc_t    *shellfc;
        gmx_bool          bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
        gmx_bool          bResetCountersHalfMaxH = FALSE;
        gmx_bool          bTemp, bPres, bTrotter;
        real              dvdl_constr;
        rvec             *cbuf        = nullptr;
        int               cbuf_nalloc = 0;
        matrix            lastbox;
        int               lamnew  = 0;
        /* for FEP */
        int               nstfep = 0;
        double            cycles;
        real              saved_conserved_quantity = 0;
        real              last_ekin                = 0;
        t_extmass         MassQ;
        int             **trotter_seq;
        char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
        int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition*/


        /* PME load balancing data for GPU kernels */
        pme_load_balancing_t *pme_loadbal      = nullptr;
        gmx_bool              bPMETune         = FALSE;
        gmx_bool              bPMETunePrinting = FALSE;

        /* Interactive MD */
        gmx_bool          bIMDstep = FALSE;

        /* Domain decomposition could incorrectly miss a bonded
           interaction, but checking for that requires a global
           communication stage, which does not otherwise happen in DD
           code. So we do that alongside the first global energy reduction
           after a new DD is made. These variables handle whether the
           check happens, and the result it returns. */
        bool              shouldCheckNumberOfBondedInteractions = false;
        int               totalNumberOfBondedInteractions       = -1;

        SimulationSignals signals;
        // Most global communication stages don't propagate mdrun
        // signals, and will use this object to achieve that.
        SimulationSignaller        nullSignaller = SimulationSignaller(nullptr, nullptr, nullptr, false, false);

        gmx_edsam                 *ed = nullptr;

        const gmx_bool             bRerunMD      = FALSE;
        int                        nstglobalcomm;

        std::unique_ptr<t_state>   stateInstance;
        t_state   *                state;
        t_mdatoms                 *mdatoms;

        const ContinuationOptions &continuationOptions;
        bool                       startingFromCheckpoint;
        real max_hours;
        int  nstSignalComm;

    private:
        static void reset_all_counters(FILE *fplog, const gmx::MDLogger &mdlog, t_commrec *cr,
                                       gmx_int64_t step,
                                       gmx_int64_t *step_rel, t_inputrec *ir,
                                       gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                                       gmx_walltime_accounting_t walltime_accounting,
                                       struct nonbonded_verlet_t *nbv,
                                       struct gmx_pme_t *pme);
};
}


#endif //GROMACS_DATAMANAGER_H
