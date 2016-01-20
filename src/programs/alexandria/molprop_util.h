/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/topology/atomprop.h"

#include "molprop.h"
#include "molselect.h"
#include "poldata.h"

struct t_topology;

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

namespace alexandria
{

    void generate_composition(std::vector<MolProp> &mp,
                              const Poldata &pd);
                                 
    void generate_formula(std::vector<MolProp> &mp, 
                          gmx_atomprop_t ap);

    /* Return number of atoms, 0 means failure */
    bool molprop_2_topology2(MolProp mp, gmx_atomprop_t ap,
                             t_symtab *tab, const char *lot,
                             t_topology **top, const char *q_algorithm,
                             rvec **x, t_params plist[F_NRE],
                             int nexcl, t_excls **excls);

    int merge_doubles(std::vector<alexandria::MolProp> &mp,
                      char *doubles, bool bForceMerge);

    int merge_xml(int nfile, char **infiles,
                  std::vector<alexandria::MolProp> &mp,
                  char *outf, char *sorted, char *doubles,
                  gmx_atomprop_t ap,
                  const Poldata &pd,
                  bool bForceMerge);

    typedef struct {
        int    n;
        int   *count;
        char **method, **basis, **type, **lot;
        int    nconf;
        char **conf;
    } t_qmcount;
    
    /* Check the available molprops to see what kind of calculations are stored in there */
    t_qmcount *find_calculations(std::vector<alexandria::MolProp> &mp,
                                 MolPropObservable mpo, 
                                 const char *fc_str);
    
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
    void MolPropSort(std::vector<MolProp> &mp,
                     MolPropSortAlgorithm mpsa,
                     gmx_atomprop_t apt,
                     gmx_molselect *gms);
    
} // namespace alexandria

#endif
    
