#include <sys/types.h>
#include <signal.h>

int
main()
{
  return *(signal (0, 0)) (0) == 1;
}    
