#include <math.h>
#include "typedefs.h"
#include "callf77.h"
#include "fatal.h"

/* This file provides the interface to fortran routines in a machine
 * independent way.
 */

extern void FUNCTION(ffillbuf) (void);

void fillbuf(void)
{
#ifdef USEF77
#ifdef FINVSQRT
  FUNCTION(ffillbuf)();
#endif
#else
  fatal_error(0,"fillbuf called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

real FUNCTION(cerfc)(real x)
{
	return erfc(x);
}
