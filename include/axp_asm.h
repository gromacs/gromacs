#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Function definitions for alpha AXP assembly routines */

/* Double precision inverse square root */
void sqrtiv_(double *indata, double *outdata, int *n);

/* Double precision square root */
void sqrtv_(double *indata, double *outdata, int *n);

/* Single precision inverse square root */
void ssqrtiv_(float *indata, float *outdata, int *n);

/* Single precision square root */
void ssqrtv_(float *indata, float *outdata, int *n);

