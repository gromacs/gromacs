/*
 * $Id: molprop_util.h,v 1.18 2009/05/27 13:44:55 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */

#ifndef _molprop_util_hpp
#define _molprop_util_hpp

#include "grompp.h"
#include "atomprop.h"
#include "molprop.hpp"
#include "poldata.h"
#include "molselect.hpp"

extern void generate_formula(std::vector<alexandria::MolProp>& mp,gmx_atomprop_t ap);

alexandria::MolProp atoms_2_molprop(char *molname,int natoms,char **smnames,
                                    gmx_atomprop_t ap,gmx_poldata_t pd);

/* Return number of atoms, 0 means failure */
extern int molprop_2_topology(alexandria::MolProp mp,gmx_atomprop_t ap,
                              gmx_poldata_t pd,
                              t_symtab *tab,const char *lot,
                              t_topology **top,const char *q_algorithm,
                              rvec **x,t_params plist[F_NRE],
                              int nexcl,t_excls **excls);

extern void merge_doubles(std::vector<alexandria::MolProp> &mp,
                          char *doubles,gmx_bool bForceMerge);
				     				     
extern void merge_xml(int nfile,char **infiles,
                      std::vector<alexandria::MolProp> &mp,
                      char *outf,char *sorted,char *doubles,
                      gmx_atomprop_t ap,gmx_poldata_t pd,
                      gmx_bool bForceMerge,gmx_bool bForceGenComp,
                      double th_toler,double ph_toler);

typedef struct {
    int  n;
    int  *count;
    char **method,**basis,**type,**lot;
    int  nconf;
    char **conf;
} t_qmcount;

/* Check the available molprops to see what kind of calculations are stored in there */
extern t_qmcount *find_calculations(std::vector<alexandria::MolProp> mp,
                                    MolPropObservable mpo,const char *fc_str);

/*! \brief
 * Enumerated type for MolPropSort function
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum MolPropSortAlgorithm { 
    MPSA_MOLNAME, 
    MPSA_FORMULA, 
    MPSA_COMPOSITION, 
    MPSA_SELECTION,
    MPSA_NR
};

/*! \brief
 * Sorts a vector of molprops
 *
 * Function that uses the std::sort routine and can apply different sorting
 * keys.
 * 
 * \param[inout]  mp        The vector of MolProp
 * \param[in]     mpsa      The algorithm used for sorting
 * \param[in]     apt       Database of atom properties
 * \param[in]     mgs       Optional structure containing selection criteria
 * \ingroup module_alexandria
 */
extern void MolPropSort(std::vector<alexandria::MolProp> &mp,
                        MolPropSortAlgorithm mpsa,gmx_atomprop_t apt,
                        gmx_molselect_t gms);

#endif
