#include <pthread.h>

/* Since most fortran compilers dont support threads started
 * in a calling c program we call these wrapper syncronization
 * routines from the fortran innerloops
 */
#ifdef USE_FORTRAN
void FUNCTION(inlsync)(int *nri,int *nthreads,int *count,int *ii0,
		       int *ii1, pthread_mutex_t *mtx)
{
  int t0,t1;
  pthread_mutex_lock(mtx);
  t0=*count;
  t1=t0+(*nri-t0)/nthreads+1;
  *count=t1;
  pthread_mutex_unlock(mtx);
  *ii0=t0;
  *ii1=t1;
}
#endif

