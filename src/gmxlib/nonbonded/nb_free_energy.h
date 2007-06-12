#ifndef _nb_free_energy_h_
#define _nb_free_energy_h_

#include <typedefs.h>

void
gmx_nb_free_energy_kernel(int                  icoul,
                          int                  ivdw,
                          int                  nri,
                          int *                iinr,
                          int *                jindex,
                          int *                jjnr,
                          int *                shift,
                          real *               shiftvec,
                          real *               fshift,
                          int *                gid,
                          real *               x,
                          real *               f,
                          real *               chargeA,
                          real *               chargeB,
                          real                 facel,
                          real                 krf,
                          real                 crf,
                          real                 ewc,
                          real *               Vc,
                          int *                typeA,
                          int *                typeB,
                          int                  ntype,
                          real *               nbfp,
                          real *               Vvdw,
                          real                 tabscale,
                          real *               VFtab,
                          real *               lambda,
                          int                  nlambda,
                          real *               dvdlambda,
                          real *               deltaH,
                          real                 alpha,
                          int                  lam_power,
                          real                 def_sigma6,
                          int *                outeriter,
                          int *                inneriter);

#endif

