#ifdef __CYGWIN__
    /* Pipes need POSIX things, not just std ones */
    #define _POSIX_C_SOURCE 200809L
#endif
#include <stdio.h>

int
main()
{
  FILE *fp;

  fp = popen("/tmp/xyz","r");
  return (fp==NULL);
}
