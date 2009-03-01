#include <mkl_dfti.h>

int
main()
{
  DFTI_DESCRIPTOR *desc;
  MKL_LONG nx = 10;

  DftiCreateDescriptor(&desc,DFTI_SINGLE,DFTI_COMPLEX,1,nx);
}
