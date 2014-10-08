/* TODO find out what this file should be called */
#ifndef GMX_EWALD_PME_GRID
#define GMX_EWALD_PME_GRID

#ifdef __cplusplus
extern "C" {
#endif
#include "config.h"

#include "pme-internal.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef GMX_MPI
void
gmx_sum_qgrid_dd(struct gmx_pme *pme, real *grid, int direction);
#endif

int
copy_pmegrid_to_fftgrid(struct gmx_pme *pme, real *pmegrid, real *fftgrid, int grid_index);

int
copy_fftgrid_to_pmegrid(struct gmx_pme *pme, const real *fftgrid, real *pmegrid, int grid_index,
                            int nthread, int thread);

void
wrap_periodic_pmegrid(struct gmx_pme *pme, real *pmegrid);

void
unwrap_periodic_pmegrid(struct gmx_pme *pme, real *pmegrid);

void
pmegrid_init(pmegrid_t *grid,
                  int cx, int cy, int cz,
                  int x0, int y0, int z0,
                  int x1, int y1, int z1,
                  gmx_bool set_alignment,
                  int pme_order,
                  real *ptr);

void
pmegrids_init(pmegrids_t *grids,
                   int nx, int ny, int nz, int nz_base,
                   int pme_order,
                   gmx_bool bUseThreads,
                   int nthread,
                   int overlap_x,
              int overlap_y);

void
pmegrids_destroy(pmegrids_t *grids);

void
make_gridindex5_to_localindex(int n, int local_start, int local_range,
                              int **global_to_local,
                              real **fraction_shift);

void
set_grid_alignment(int gmx_unused *pmegrid_nz, int gmx_unused pme_order);

void
reuse_pmegrids(const pmegrids_t *old, pmegrids_t *new);

#ifdef __cplusplus
}
#endif

#endif

