
#ifndef _localpressure_h_
#define _localpressure_h_

#include "types/simple.h"

typedef struct
{
    int        nx;
    int        ny;
    int        nz;
    matrix *   current_grid;   /* allocated length is nx*ny*nz */
    matrix *   longrange_grid;   /* allocated length is nx*ny*nz */
    matrix *   sum_grid;
    int        nframes;
    real       spacing;
    matrix     box;
    matrix     invbox;
    bool       bLR;
} gmx_localp_grid_t;


void 
grid_lookup(gmx_localp_grid_t *grid,rvec pt, int *i, int *j, int *k);


void
gmx_spread_local_virial_on_grid(gmx_localp_grid_t *    grid,
                                real                   ix,
                                real                   iy,
                                real                   iz,
                                real                   jx,
                                real                   jy,
                                real                   jz,
                                real                   fx,
                                real                   fy,
                                real                   fz);

void
gmx_spread_local_virial_on_grid_mat(gmx_localp_grid_t *    grid,
                                    real                   ix, 
                                    real                   iy,
                                    real                   iz,
                                    matrix		 pressure);

void
gmx_spread_local_gaspressure_on_grid(gmx_localp_grid_t * grid,
                                     rvec x[],
                                     rvec v[],
                                     real mass[],
                                     int natoms);
#endif
