#include <math.h>
#include "typedefs.h"
#include "callf77.h"
#include "fatal.h"

/* This file provides the interface to fortran routines in a machine
 * independent way.
 */



real FUNCTION(cerfc)(real x)
{
	return erfc(x);
}

real FUNCTION(cpow)(real x,real y)
{
	return pow(x,y);
}
