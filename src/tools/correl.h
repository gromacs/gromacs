#ifndef _correl_h
#define _correl_h

#include "typedefs.h"
#include "fftw_wrapper.h"

typedef struct {
  int n;
  rfftw_plan p_fw,p_bw;
  real *buf1,*buf2,*abuf;
} correl_t;

extern void four1(real data[],int nn,int isign);
extern void correl(real data1[],real data2[],int n,real ans[]);
extern void correl_fftw(correl_t *c,real data1[],real data2[],real ans[]);
extern correl_t *init_correl(int n);
extern void done_correl(correl_t *c);

#endif
