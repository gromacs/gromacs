#include <stdio.h>
#include <stdlib.h>
#include "gmx_fatal.h"
#include "macros.h"
#include "gmx_qhop_db.h"

int main(int argc,char *argv[])
{
  gmx_qhop_db db;
  char *donors[] = { "H3O+", "ACE" };
  char *acceptors[] = { "H2O", "GLU" };
  t_qhop_parameters qp;
  int i,j;
  
  if ((db = gmx_qhop_db_read("ffoplsaa")) == NULL) 
    gmx_fatal(FARGS,"Can not read qhop database information");
  if (gmx_qhop_db_write("koe.dat",db) != 1)
    gmx_fatal(FARGS,"Can not write qhop database information");
  
  for(i=0; (i<asize(donors)); i++) {
    for(j=0; (j<asize(acceptors)); j++) {
      if (gmx_qhop_db_get_parameters(db,donors[i],acceptors[j],&qp) == 1) {
	printf("Found qhop parameters for donor %s and acceptor %s\n",
	       donors[i],acceptors[j]);
      }
      else {
	printf("Could not find qhop parameters for donor %s and acceptor %s\n",
	       donors[i],acceptors[j]);
      }
    }
  }
  
  if (gmx_qhop_db_done(db) != 1)
    gmx_fatal(FARGS,"Error destroying qhop data");
    
  return 0;
}
