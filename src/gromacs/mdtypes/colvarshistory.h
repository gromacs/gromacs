#ifndef GMX_MDLIB_COLVARSHISTORY_H
#define GMX_MDLIB_COLVARSHISTORY_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

/* Helper structure to be able to make colvars group(s) whole
 *
 * To also make colvars group(s) whole, we save the last whole configuration
 * of the atoms in the checkpoint file.
 */
typedef struct colvarshistory_t
{
    gmx_bool         bFromCpt;     // Did we start from a checkpoint file?
    int              n_atoms; // Number of colvars atoms
    rvec*            xa_old_whole; // Last known whole positions of the colvars atoms
    rvec*            xa_old_whole_p; // Pointer to these positions
} colvarshistory_t;

#endif
