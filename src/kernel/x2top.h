#ifndef _x2top_h
#define _x2top_h
	
#include <stdio.h>
	
typedef struct {
  char *elem,*type;
  int  nbonds;
  char **bond;
} t_nm2type;

extern t_nm2type *rd_nm2type(int *nnm);
/* Read the name 2 type database. nnm is the number of entries */

extern void dump_nm2type(FILE *fp,int nnm,t_nm2type nm2t[]);
/* Dump the database for debugging. Can be reread by the program */

extern char *nm2type(int nnm,t_nm2type nm2t[],char *elem,int nbonds);
/* Try to determine the atomtype (force field dependent) for an element */

#endif
