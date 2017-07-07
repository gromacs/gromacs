#include "gmxpre.h"

#include "helperFunctionsTemporary.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdtypes/commrec.h"

void checkNumberOfBondedInteractions(FILE *fplog, t_commrec *cr, int totalNumberOfBondedInteractions,
                                            gmx_mtop_t *top_global, gmx_localtop_t *top_local, t_state *state,
                                            bool *shouldCheckNumberOfBondedInteractions)
{
    if (*shouldCheckNumberOfBondedInteractions)
    {
        if (totalNumberOfBondedInteractions != cr->dd->nbonded_global)
        {
            dd_print_missing_interactions(fplog, cr, totalNumberOfBondedInteractions, top_global, top_local, state); // Does not return
        }
        *shouldCheckNumberOfBondedInteractions = false;
    }
}