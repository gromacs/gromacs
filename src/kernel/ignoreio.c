#include <stdio.h>
#include <strdb.h>

void read_ignore()
{
  int  nign;
  char **ign;
  char *db = "ignore.dat";
  
  nign = get_lines(db,&ign);
  
}
