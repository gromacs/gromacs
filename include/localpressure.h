
#ifndef _localpressure_h_
#define _localpressure_h_

#include "types/simple.h"

typedef struct
{
    int        nx;
    int        ny;
    int        nz;
    matrix *   current_grid;   /* allocated length is nx*ny*nz */
    matrix *   sum_grid;
    matrix *   sum_grid_total;
    matrix *   tmp_grid;
    int        nframes;
    int        nframes_total;
    real       spacing;
    matrix     box;
    matrix     invbox;
    gmx_bool   calc_localp;    /* true at nstcalcenergy steps */
    gmx_bool   CGlocalp;       /* Charge-group based local pressure */
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
