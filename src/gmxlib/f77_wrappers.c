#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <typedefs.h>
#include <callf77.h>

/* These are wrapper routines for f77 to use the C
 * math library routines.
 */

#ifdef USE_FORTRAN
real F77_FUNC(cerfc,CERFC)(real *x)
{
	return erfc(*x);
}

real F77_FUNC(cpow,CPOW)(real *x,real *y)
{
	return pow(*x,*y);
}
#endif
