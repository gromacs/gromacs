#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "nsgrid.h"

enum { eNL_VDWQQ, eNL_VDW, eNL_QQ, 
       eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE, 
       eNL_VDWQQ_WATER, eNL_QQ_WATER, 
       eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER, 
       eNL_NR };

#define MAX_CG 1024

typedef struct {
  int     ncg;
  int     nj;
  atom_id jcg[MAX_CG];
} t_ns_buf;

typedef unsigned long t_excl;

typedef struct {
  atom_id  *simple_aaj;
  t_grid   *grid;
  t_excl   *bexcl;
  bool     *bHaveVdW;
  t_ns_buf **ns_buf;
  bool     *bExcludeAlleg;
  int      nra_alloc;
  int      cg_alloc;
  atom_id  **nl_sr;
  int      *nsr;
  atom_id  **nl_lr_ljc;
  atom_id  **nl_lr_one;
  int      *nlr_ljc;
  int      *nlr_one;
} gmx_ns_t;
