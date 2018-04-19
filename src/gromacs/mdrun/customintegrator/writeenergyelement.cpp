//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/utility/fatalerror.h"
#include "writeenergyelement.h"

void gmx::WriteEnergyElement::run()
{
    if (dataManager.bCalcEner)
    {
        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

        /* use the directly determined last velocity, not actually the averaged half steps */
        dataManager.enerd->term[F_ETOT] = dataManager.enerd->term[F_EPOT] + dataManager.enerd->term[F_EKIN];

        if (integratorHasConservedEnergyQuantity(integrator.inputrec))
        {
            {
                dataManager.enerd->term[F_ECONSERVED] = dataManager.enerd->term[F_ETOT] + NPT_energy(integrator.inputrec, dataManager.state, &dataManager.MassQ);
            }
        }
        /* #########  END PREPARING EDR OUTPUT  ###########  */
    }

    /* Output stuff */
    if (MASTER(integrator.cr))
    {
        if (dataManager.bCalcEner)
        {
            upd_mdebin(dataManager.mdebin, dataManager.bDoDHDL, dataManager.bCalcEnerStep,
                       dataManager.t, dataManager.mdatoms->tmass, dataManager.enerd, dataManager.state,
                       integrator.inputrec->fepvals, integrator.inputrec->expandedvals, dataManager.lastbox,
                       dataManager.shake_vir, dataManager.force_vir, dataManager.total_vir, dataManager.pres,
                       dataManager.ekind, dataManager.mu_tot, integrator.constr);
        }
        else
        {
            upd_mdebin_step(dataManager.mdebin);
        }

        gmx_bool do_dr  = do_per_step(dataManager.step, integrator.inputrec->nstdisreout);
        gmx_bool do_or  = do_per_step(dataManager.step, integrator.inputrec->nstorireout);

        print_ebin(mdoutf_get_fp_ene(dataManager.outf), dataManager.do_ene, do_dr, do_or, dataManager.do_log ? integrator.fplog : nullptr,
                   dataManager.step, dataManager.t,
                   eprNORMAL, dataManager.mdebin, integrator.fcd, dataManager.groups, &(integrator.inputrec->opts), integrator.inputrec->awh);

        if (integrator.inputrec->bPull)
        {
            pull_print_output(integrator.inputrec->pull_work, dataManager.step, dataManager.t);
        }

        // TODO: does that belong here or rather at end of loop?

        if (do_per_step(dataManager.step, integrator.inputrec->nstlog))
        {
            if (fflush(integrator.fplog) != 0)
            {
                gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
            }
        }
    }

}
