#ifndef _correl_h
#define _correl_h

#include "typedefs.h"
#include "gmx_fft.h"


extern void correl(real data1[],real data2[],int n,real ans[]);


typedef struct {
  int n;
  gmx_fft_t  fft_setup;
  real *buf1,*buf2,*abuf;
} correl_t;

extern correl_t *init_correl(int n);
extern void done_correl(correl_t *c);

#endif
