#include <stdio.h>

int
main()
{
  FILE *fp;

  fp = popen("/tmp/xyz","r");
  return (fp==NULL);
}
