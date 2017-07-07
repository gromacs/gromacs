
#ifndef _helperFunctionsTemporary_h
#define _helperFunctionsTemporary_h

#include <stdio.h>

#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/group.h"


struct t_commrec;
struct t_inputrec;
struct t_state;
struct t_graph;
struct gmx_shellfc_t;

void checkNumberOfBondedInteractions(FILE *fplog, t_commrec *cr, int totalNumberOfBondedInteractions,
                                            gmx_mtop_t *top_global, gmx_localtop_t *top_local, t_state *state,
                                            bool *shouldCheckNumberOfBondedInteractions);

#endif