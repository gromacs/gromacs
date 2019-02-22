/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_GROMPP_IMPL_H
#define GMX_GMXPREPROCESS_GROMPP_IMPL_H

#include <string>

#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#define MAXSLEN 32

/*! \libinternal \brief
 * Describes a single parameter in a force field.
 */
struct FFParameter
{
    //! The atom list (eg. bonds: particle, i = a[0] (ai), j = a[1] (aj))
    std::array<int, MAXATOMLIST>    atoms = {-1};
    //! Force parameters (eg. b0 = c[0])
    std::array<real, MAXFORCEPARAM> forceParm {NOTSET};
    //! A string (instead of parameters),read from the .rtp file in pdb2gmx
    std::string                     paramString;

    //! @{
    //! Access the individual elements set for the parameter.
    const int       &ai() const { return atoms[0]; }
    int             &ai() { return atoms[0]; }
    const int       &aj() const { return atoms[1]; }
    int             &aj() { return atoms[1]; }
    const int       &ak() const { return atoms[2]; }
    int             &ak() { return atoms[2]; }
    const int       &al() const { return atoms[3]; }
    int             &al() { return atoms[3]; }
    const int       &am() const { return atoms[4]; }
    int             &am() { return atoms[4]; }

    const real      &c0() const { return forceParm[0]; }
    real            &c0() { return forceParm[0]; }
    const real      &c1() const { return forceParm[1]; }
    real            &c1() { return forceParm[1]; }
    const real      &c2() const { return forceParm[2]; }
    real            &c2() { return forceParm[2]; }
    //! @}
};

/*! \libinternal \brief
 * All parameters for a system.
 */
struct SystemParameters
{                       // NOLINT (clang-analyzer-optin.performance.Padding)
    //! The different parameters in the system.
    std::vector<FFParameter> param;
    //! CMAP grid spacing.
    int                      cmakeGridSpacing = -1;
    //! Number of cmap angles.
    int                      cmapAngles = -1;
    //! CMAP grid data.
    std::vector<real>        cmap;
    //! The five atomtypes followed by a number that identifies the type.
    std::vector<int>         cmapAtomTypes;

    //! Number of parameters.
    int nr() const { return param.size(); }
    //! Elements in cmap grid data.
    int ncmap() const { return cmap.size(); }
    //! Number of elements in cmapAtomTypes.
    int nct() const { return cmapAtomTypes.size(); }
};

struct t_excls
{
    int            nr;      /* The number of exclusions             */
    int           *e;       /* The excluded atoms                   */
};


/*! \libinternal \brief
 * Holds the molecule information during preprocessing.
 */
struct MoleculeInformation
{
    //! Name of the molecule.
    char                              **name = nullptr;
    //!Number of exclusions per atom.
    int                                 nrexcl = 0;
    //! Have exclusions been generated?.
    bool                                excl_set = false;
    //! Has the mol been processed.
    bool                                bProcessed = false;
    //! Atoms in the moelcule.
    t_atoms                             atoms;
    //! Charge groups in the molecule
    t_block                             cgs;
    //! Molecules separated in datastructure.
    t_block                             mols;
    //! Exclusions in the molecule.
    t_blocka                            excls;
    //! Parameters in old style.
    std::array<SystemParameters, F_NRE> plist;

    /*! \brief
     * Initializer.
     *
     * This should be removed as soon as the underlying datastructures
     * have been cleaned up to use proper initialization and can be copy
     * constructed.
     */
    void initMolInfo();

    /*! \brief
     * Partial clean up function.
     *
     * Should be removed once this datastructure actually owns all its own memory and
     * elements of it are not stolen by other structures and properly copy constructed
     * or moved.
     * Cleans up the mols and plist datastructures but not cgs and excls.
     */
    void partialCleanUp();

    /*! \brief
     * Full clean up function.
     *
     * Should be removed once the destructor can always do this.
     */
    void fullCleanUp();
};

struct t_mols
{
    std::string name;
    int         nr;
};

bool is_int(double x);
/* Returns TRUE when x is integer */

#endif
