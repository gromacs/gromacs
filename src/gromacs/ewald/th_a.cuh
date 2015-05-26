#ifndef GMX_EWALD_TH_A_H
#define GMX_EWALD_TH_A_H

#include "gromacs/utility/real.h"

const int TH = 32;
enum th_id {
  TH_ID_THETA, TH_ID_DTHETA, TH_ID_FRACTX, TH_ID_COEFFICIENT,

  TH_ID_GRID,
  TH_ID_I0, TH_ID_J0, TH_ID_K0,
  TH_ID_THX, TH_ID_THY, TH_ID_THZ,

  TH_ID_END
};
enum th_loc {
  TH_LOC_HOST, TH_LOC_CUDA, TH_LOC_END
};

real *th_a(th_id id, int thread, int size, th_loc loc);
int *th_i(th_id id, int thread, int size, th_loc loc);

void th_cpy(void *dest, void *src, int size, th_loc dest_loc);

#endif
