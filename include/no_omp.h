#ifndef _NO_OMP_H_
#define _NO_OMP_H_

/* This header file is meant to enable compiling OpenMP-parallelized code 
   without OpenMP, for single-threaded use. It should be included as an 
   alternative for omp.h 
 */

#ifndef GMX_OMPENMP
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num()  { return 0; }
#endif /* GMX_OMPENMP */

#endif /* _NO_OMP_H_ */ 
