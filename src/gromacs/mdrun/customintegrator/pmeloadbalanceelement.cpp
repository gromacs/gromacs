//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/mdrun.h"


#include "pmeloadbalanceelement.h"

void gmx::PMELoadBalanceElement::loopSetup()
{
    pme_loadbal_init(&dataManager.pme_loadbal,
                     integrator.cr,
                     integrator.mdlog,
                     integrator.inputrec,
                     dataManager.state->box,
                     integrator.fr->ic,
                     integrator.fr->nbv->listParams.get(),
                     integrator.fr->pmedata,
                     use_GPU(integrator.fr->nbv),
                     &dataManager.bPMETunePrinting);
}

void gmx::PMELoadBalanceElement::run()
{
    /* PME grid + cut-off optimization with GPUs or PME nodes */
    pme_loadbal_do(dataManager.pme_loadbal,
                   integrator.cr,
                   (integrator.mdrunOptions.verbose && MASTER(integrator.cr)) ? stderr : nullptr,
                   integrator.fplog, integrator.mdlog,
                   integrator.inputrec, integrator.fr,
                   dataManager.state, integrator.wcycle,
                   dataManager.step, dataManager.step_rel,
                   &dataManager.bPMETunePrinting);
}

void gmx::PMELoadBalanceElement::loopTeardown()
{
    pme_loadbal_done(
            dataManager.pme_loadbal,
            integrator.fplog, integrator.mdlog,
            use_GPU(integrator.fr->nbv)
            );
}
