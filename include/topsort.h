#ifndef _topsort_h
#define _topsort_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef CPLUSPLUS 
extern "C" {
#endif

#include "typedefs.h"

  /* Analyze the idef and set idef->ilsort */
  extern void gmx_analyze_ilist_fe(t_idef *idef,t_inputrec *ir);

  /* Sort all the bonded ilists in idef to have the perturbed ones at the end
   * and set nr_nr_nonperturbed in ilist.
   */
  extern void gmx_sort_ilist_fe(t_idef *idef);

#ifdef CPLUSPLUS
}
#endif

#endif
