#ifndef GMX_EWALD_PME_SOLVE_H
#define GMX_EWALD_PME_SOLVE_H

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_pme;

void get_pme_ener_vir_q(struct gmx_pme *work, int nthread,
                   real *mesh_energy, matrix vir);

void get_pme_ener_vir_lj(struct gmx_pme *work, int nthread,
                         real *mesh_energy, matrix vir);

int solve_pme_yzx(struct gmx_pme *pme, t_complex *grid,
                  real ewaldcoeff, real vol,
                  gmx_bool bEnerVir,
                  int nthread, int thread);

int solve_pme_lj_yzx(struct gmx_pme *pme, t_complex **grid, gmx_bool bLB,
                     real ewaldcoeff, real vol,
                     gmx_bool bEnerVir, int nthread, int thread);

#ifdef __cplusplus
}
#endif

#endif
