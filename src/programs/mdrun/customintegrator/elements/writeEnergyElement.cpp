#include "writeEnergyElement.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/utility/fatalerror.h"

#include "programs/mdrun/customintegrator/statemanager.h"

void WriteEnergy::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void WriteEnergy::run()
{

    /* Output stuff */
    if (MASTER(_data->cr))
    {
        if (_data->fplog && _data->do_log && _data->bDoExpanded)
        {
            /* only needed if doing expanded ensemble */
            PrintFreeEnergyInfoToFile(_data->fplog, _data->ir->fepvals, _data->ir->expandedvals, _data->ir->bSimTemp ? _data->ir->simtempvals : nullptr,
                                      _data->state_global->dfhist, _data->state->fep_state, _data->ir->nstlog, _data->step);
        }
        if (_data->bCalcEner)
        {
            upd_mdebin(_data->mdebin, _data->bDoDHDL, _data->bCalcEnerStep,
                       _data->t, _data->mdatoms->tmass, _data->enerd, _data->state,
                       _data->ir->fepvals, _data->ir->expandedvals, _data->lastbox,
                       _data->shake_vir, _data->force_vir, _data->total_vir, _data->pres,
                       _data->ekind, _data->mu_tot, _data->constr);
        }
        else
        {
            upd_mdebin_step(_data->mdebin);
        }

        gmx_bool do_dr  = do_per_step(_data->step, _data->ir->nstdisreout);
        gmx_bool do_or  = do_per_step(_data->step, _data->ir->nstorireout);

        print_ebin(mdoutf_get_fp_ene(_data->outf), _data->do_ene, do_dr, do_or, _data->do_log ? _data->fplog : nullptr,
                   _data->step, _data->t,
                   eprNORMAL, _data->mdebin, _data->fcd, _data->groups, &(_data->ir->opts));

        if (_data->ir->bPull)
        {
            pull_print_output(_data->ir->pull_work, _data->step, _data->t);
        }

        if (do_per_step(_data->step, _data->ir->nstlog))
        {
            if (fflush(_data->fplog) != 0)
            {
                gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
            }
        }
    }
}