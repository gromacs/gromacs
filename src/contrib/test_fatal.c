#include "gromacs/utility/fatalerror.h"

void my_func(char *msg)
{
  fprintf(stderr,"Welcome to my_func\n%s\n",msg);
  exit(1);
}

int main(int argc,char *argv[])
{
  int n = -3;
  int choice;
  char *end;
  
  /* set_gmx_error_handler(my_func);*/
  
  if (argc <= 1)
    gmx_fatal(FARGS,"Expected an integer argument to %s",argv[0]);
  choice = strtol(argv[1], &end, 10); 
  if (end!='\0')
    gmx_fatal(FARGS,"Expected an integer argument to %s",argv[0]);
  
  
  switch (choice) {
  case 1:
    gmx_error("pme","oooo");
    break;
  case 2:
    gmx_fatal(FARGS,"Passing string %s to you %f  %d %x","lll",8.3,34,34);
    break;
  case 3:
    range_check(n,1,5);
    break;
  default:
    range_check(choice,1,3);
  }
  return 0;
}
