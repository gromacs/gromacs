#include <stdio.h>
#include <stdlib.h>
#include "poldata.h"

int main(int argc,char*argv[])
{
  gmx_poldata_t pd;
  
  if (argc < 3) {
    fprintf(stderr,"Usage: %s infile outfile\n",argv[0]);
    exit(1);
  }
  pd = gmx_poldata_read(argv[1]);
  gmx_poldata_write(argv[2],pd);
  
  return 0;
}
