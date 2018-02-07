#ifndef GMXEXTERN_HPP
#define GMXEXTERN_HPP


#include "gromacs/topology/atoms.h"

extern "C" {

//void get_index(t_atoms *atoms, const char *fnm, int ngrps,
//                int isize[], int *index[], char *grpnames[]);
void read_eigenvectors(const char* file, int* natoms, gmx_bool* bFit,
                                  rvec** xref, gmx_bool* bDMR, rvec** xav, gmx_bool* bDMA, int* nvec,
                                  int** eignr, rvec*** eigvec, real** eigval);
/**
 * declare write_eigenvectors this function
 * eigio.h is not in the includes in recent gromacs versions, but its functions
 * are availabe from the gmxana libraries
 */
void write_eigenvectors(const char* trnname, const int natoms, const real mat[],
                                   const gmx_bool bReverse, const int begin, const int end,
                                   const int WriteXref, const rvec* xref, const gmx_bool bDMR,
                                   const rvec xav[], const gmx_bool bDMA, const real* eigval);

void get_stx_coordnum(const char *infile,int *natoms);
void read_stx_conf(const char *infile,char *title,t_atoms *atoms,
                   rvec x[],rvec *v,int *ePBC,matrix box);

enum { eWXR_NO,
       eWXR_YES,
       eWXR_NOFIT }; //< from eigio.h (see write_eigenvectors)
};

#endif
