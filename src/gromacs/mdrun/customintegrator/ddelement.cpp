//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/inputrec.h"

#include "ddelement.h"

void gmx::DDElement::run()
{
    if (dataManager.bFirstStep && integrator.inputrec->bContinuation)
    {
        return;
    }
    gmx_bool bMasterState = FALSE;

    // TODO: Does that belong here?
    /* Correct the new box if it is too skewed */
    if (inputrecDynamicBox(integrator.inputrec))
    {
        if (correct_box(integrator.fplog, dataManager.step, dataManager.state->box, dataManager.graph))
        {
            bMasterState = TRUE;
        }
    }
    if (DOMAINDECOMP(integrator.cr) && bMasterState)
    {
        dd_collect_state(integrator.cr->dd, dataManager.state, integrator.state_global);
    }

    if (DOMAINDECOMP(integrator.cr))
    {
        /* Repartition the domain decomposition */
        dd_partition_system(integrator.fplog, dataManager.step,
                            integrator.cr, bMasterState,
                            dataManager.nstglobalcomm,
                            integrator.state_global,
                            integrator.top_global,
                            integrator.inputrec,
                            dataManager.state,
                            &dataManager.f,
                            integrator.mdAtoms,
                            dataManager.top,
                            integrator.fr,
                            integrator.vsite,
                            integrator.constr,
                            integrator.nrnb,
                            integrator.wcycle,
                            dataManager.do_verbose && !dataManager.bPMETunePrinting);
        dataManager.shouldCheckNumberOfBondedInteractions = true;
        update_realloc(dataManager.upd, dataManager.state->natoms);
    }
}
