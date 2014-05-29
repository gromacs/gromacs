#ifndef GMX_REPLICAEXCHANGE_STATE_H
#define GMX_REPLICAEXCHANGE_STATE_H

#include "gromacs/legacyheaders/types/state.h"
#include "gromacs/legacyheaders/types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Exchange state between two replicas
 *
 * \param[in]    ms     Multi-simulation manager object
 * \param[in]    b      Index of exchange partner
 * \param[inout] state  State object (both source and sink)
 */
void
exchange_state(const gmx_multisim_t *ms, int b, t_state *state);

/*! \brief Exchange state between two replicas
 */
void
copy_state_nonatomdata(const t_state *source, t_state *sink);

#ifdef __cplusplus
}
#endif

#endif
