#ifndef _statemanager_h
#define _statemanager_h

#include "gmxpre.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/logger.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/essentialdynamics/edsam.h"

#include <memory>
#include <map>
#include <string>

#include "programs/mdrun/repl_ex.h"
struct gmx_membed_t;

    /*
        The StateManager class holds all the data that will be needed for the
        integrator.

        *** !!! Temporary Solution: !!!***
        Currently this class holds all the data that the standard do_md()
        function in the GROMACS code holds. All elements in the integrator
        sequences (preRunSequence and runSequence) have full read and write
        access to all the data.

        A general discussion is needed to decide on the data structures,
        categorization and data access levels for the elements. Please see the
        accompanying document, section IV (Data reorganization) for more details.
    */
class StateManager
{
public:
    /*
       Variables which were originally passed as arguments in do_md()
     */
    FILE *fplog;
    t_commrec *cr;
    const gmx::MDLogger *mdlog;
    int nfile;
    const t_filenm *fnm;
    const gmx_output_env_t *oenv;
    gmx_bool bVerbose;
    int nstglobalcomm;
    gmx_vsite_t *vsite;
    gmx_constr_t constr;
    int stepout;
    t_inputrec *ir;
    gmx_mtop_t *top_global;
    t_fcdata *fcd;
    t_state *state_global;
    energyhistory_t *energyHistory;
    t_mdatoms *mdatoms;
    t_nrnb *nrnb;
    gmx_wallcycle_t wcycle;
    gmx_edsam_t ed;
    t_forcerec *fr;
    int repl_ex_nst;
    int repl_ex_nex;
    int repl_ex_seed;
    gmx_membed_t *membed;
    real cpt_period;
    real max_hours;
    int imdport;
    unsigned long Flags;
    gmx_walltime_accounting_t walltime_accounting;

    /*
        Variables which were originally declared as local variables in do_md()
     */
    gmx_mdoutf_t    outf;
    gmx_int64_t     step, step_rel;
    double          elapsed_time;
    double          t, t0, lam0[efptNR];
    gmx_bool        bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool        bNS, bNStList, bSimAnn, bStopCM,
                    bFirstStep, startingFromCheckpoint, bInitStep, bLastStep,
                    bBornRadii, bUsingEnsembleRestraints;
    gmx_bool          bDoDHDL, bDoFEP, bDoExpanded;
    gmx_bool          do_ene, do_log, do_verbose, bCPT;
    gmx_bool          bMasterState;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, tmp_vir, pres;
    int               i, m;
    rvec              mu_tot;
    t_vcm            *vcm;
    matrix            parrinellorahmanMu, M;
    gmx_repl_ex_t     repl_ex;
    int               nchkpt;
    gmx_localtop_t   *top;
    t_mdebin         *mdebin;
    gmx_enerdata_t   *enerd;
    PaddedRVecVector  f {};
    gmx_global_stat_t gstat;
    gmx_update_t     *upd;
    t_graph          *graph;
    gmx_groups_t     *groups;
    gmx_ekindata_t   *ekind;
    gmx_shellfc_t    *shellfc;
    gmx_bool          bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
    gmx_bool          bResetCountersHalfMaxH;
    gmx_bool          bTemp, bPres, bTrotter;
    real              dvdl_constr;
    matrix            lastbox;
    int               lamnew;
    /* for FEP */
    int               nstfep;
    double            cycles;
    real              saved_conserved_quantity;
    real              last_ekin;
    t_extmass         MassQ;
    int             **trotter_seq;
    char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
    /* compare to get_stop_condition*/
    int               handled_stop_condition = gmx_stop_cond_none;
    std::unique_ptr<t_state> stateInstance;
    t_state *                state;

    typedef struct t_simulationState
    {
        t_state state;
        gmx_enerdata_t enerd;
        // index of 
        t_commrec cr;
    }t_simulationState;
    
    std::map<std::string, t_simulationState> simstate_map;

    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t *pme_loadbal;
    gmx_bool              bPMETune;
    gmx_bool              bPMETunePrinting;

    /* Interactive MD */
    gmx_bool          bIMDstep;

#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif
    /* Domain decomposition could incorrectly miss a bonded
       interaction, but checking for that requires a global
       communication stage, which does not otherwise happen in DD
       code. So we do that alongside the first global energy reduction
       after a new DD is made. These variables handle whether the
       check happens, and the result it returns. */
    bool              shouldCheckNumberOfBondedInteractions;
    int               totalNumberOfBondedInteractions;

    gmx::SimulationSignals signals;
    // Most global communication stages don't propagate mdrun
    // signals, and will use this object to achieve that.
    gmx::SimulationSignaller nullSignaller;
    
    // pass the do_md input arguments to StateManager
    // Initialize various data
    StateManager(FILE *_fplog, t_commrec *_cr, const gmx::MDLogger &_mdlog,
                 int _nfile, const t_filenm _fnm[],
                 const gmx_output_env_t *_oenv, gmx_bool _bVerbose,
                 int _nstglobalcomm,
                 gmx_vsite_t *_vsite, gmx_constr_t _constr,
                 int _stepout, t_inputrec *_ir,
                 gmx_mtop_t *_top_global,
                 t_fcdata *_fcd,
                 t_state *_state_global,
                 energyhistory_t *_energyHistory,
                 t_mdatoms *_mdatoms,
                 t_nrnb *_nrnb, gmx_wallcycle_t _wcycle,
                 gmx_edsam_t _ed, t_forcerec *_fr,
                 int _repl_ex_nst, int _repl_ex_nex, int _repl_ex_seed,
                 gmx_membed_t *_membed,
                 real _cpt_period, real _max_hours,
                 int _imdport,
                 unsigned long _Flags,
                 gmx_walltime_accounting_t _walltime_accounting):
            outf(nullptr),
            bLastStep(FALSE),
            bDoDHDL(FALSE),
            bDoFEP(FALSE),
            bDoExpanded(FALSE),
            repl_ex(nullptr),
            nchkpt(1),
            mdebin(nullptr),
            upd(nullptr),
            graph(nullptr),
            bResetCountersHalfMaxH(FALSE),
            lamnew(0),
            nstfep(0),
            saved_conserved_quantity(0),
            last_ekin(0),
            pme_loadbal(nullptr),
            bPMETune(FALSE),
            bPMETunePrinting(FALSE),
            bIMDstep(FALSE),
            shouldCheckNumberOfBondedInteractions(FALSE),
            totalNumberOfBondedInteractions(-1),
            nullSignaller(constructNullSignaller())
    {

        // assign pointers any changes to these variables will be visible
        // outside this class
        fplog = _fplog;
        cr = _cr;
        fnm = _fnm;
        oenv = _oenv;
        vsite = _vsite;
        ir = _ir;
        top_global = _top_global;
        fcd = _fcd;
        state_global = _state_global;
        energyHistory = _energyHistory;
        mdatoms = _mdatoms;
        nrnb = _nrnb;
        fr = _fr;
        membed = _membed;
        mdlog = &_mdlog;
        
        nfile = _nfile;
        bVerbose = _bVerbose;
        nstglobalcomm = _nstglobalcomm;
        constr = _constr;
        stepout = _stepout;
        wcycle = _wcycle;
        ed = _ed;  
        repl_ex_nst = _repl_ex_nst;
        repl_ex_nex = _repl_ex_nex;
        repl_ex_seed = _repl_ex_seed;
        cpt_period = _cpt_period;
        max_hours = _max_hours;
        imdport = _imdport;
        Flags = _Flags;
        walltime_accounting = _walltime_accounting;
    };

    // Return a nullSignaller object, to be used to intialize the constructor
    gmx::SimulationSignaller constructNullSignaller()
    {
        gmx::SimulationSignaller nullSignaller(nullptr, nullptr, false, false);
        return nullSignaller;
    }

    // temporary functions to change the StateManager variables
    
    /* Essentially, the parts of do_mdrun that occur at the beginning of
    * the program, before the loop starts. Will be modularized more in
    * the future as data structures are cleaned, but for now, just
    * excised from the loop to make the logic clearer. */
    void loopSetup();

    /* code to set up for each iteration of the integrator (i.e. each
    * step). Again, mostly cutting and pasting from do_mdrun, but can
    * eventually be refactored to be more modular as data structures
    * change. */
    void stepSetup();

    /* code to end each iteration of the main integrator (i.e. each
    * step). Again, mostly cutting and pasting from do_mdrun, but can
    * eventually be refactored to be more modular as data structures
    * change. */
    void stepTeardown();

    /* Essentially, the parts of do_mdrun that occur at the end of
    * the program, before the loop starts. Will be modularized more in
    * the future as data structures are cleaned, but for now, just
    * excised from the loop to make the logic clearer. */
    void loopTeardown();
};
#endif



