#include <stdio.h>
#include "vec.h"
#include "lutab.h"

int main(int argc,char *argv[])
{
  real x,y,z,diff,av;
  int  i;

  init_lookup_table(stdout);
  printf("%12s  %12s  %12s  %12s  %12s\n","X","Y","Z","Dy","Dy/Z");
  for(i=1; (i<1000); i++) {
    x = i*1.0;
    y = invsqrt(x);
    z = 1.0/sqrt(x);
    diff = y-z;
    av   = 0.5*(y+z);
    printf("%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n",x,y,z,diff,diff/z);
  }
}
