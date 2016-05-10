#include <sys/types.h>

int main()
{
  /* Cause a compile-time error if off_t is smaller than 64 bits */
  int off_t_is_large[sizeof(off_t)-7];
  return off_t_is_large[0];
}

