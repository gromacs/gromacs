#ifdef DFFTW
#include <dfftw.h>
#endif

#ifdef SFFTW
#include <sfftw.h>
#endif

#ifdef FFTW
#include <fftw.h>
#endif

int
main()
{
#ifdef DOUBLE
  /* Cause a compile-time error if fftw_real is not double */
  int _array_ [1 - 2 * !((sizeof(fftw_real)) == sizeof(double))];
#else
  /* Cause a compile-time error if fftw_real is not float */
  int _array_ [1 - 2 * !((sizeof(fftw_real)) == sizeof(float))]; 
#endif
  return 0;
}
