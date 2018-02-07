#ifndef GMXEXTERN_HPP
#define GMXEXTERN_HPP

#include "gromacs/topology/atoms.h"
#include "gromacs/gmxana/eigio.h"

void get_stx_coordnum(const char *infile,int *natoms);
void read_stx_conf(const char *infile,
                          t_symtab *symtab, char **name, t_atoms *atoms,
                          rvec x[], rvec *v, int *ePBC, matrix box);

#endif
