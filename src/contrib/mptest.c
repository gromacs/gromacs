#include <stdio.h>
#include <stdlib.h>
#include "molprop_xml.h"

int main(int argc,char*argv[])
{
  gmx_molprop_t *mpt;
  int nmolprop;
  
  if (argc < 3) {
    fprintf(stderr,"Usage: %s infile outfile\n",argv[0]);
    exit(1);
  }
  mpt = gmx_molprops_read(argv[1],&nmolprop);
  printf("Read %d molecules from %s\n",nmolprop,argv[1]);
  gmx_molprops_write(argv[2],nmolprop,mpt);
  
  
  return 0;
}
