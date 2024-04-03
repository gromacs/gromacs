/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_GMXPREPROCESS_GROMPP_IMPL_H
#define GMX_GMXPREPROCESS_GROMPP_IMPL_H

#include <string>

#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"

namespace gmx
{
template<typename>
class ArrayRef;
}

/*! \libinternal \brief
 * Describes an interaction of a given type, plus its parameters.
 */
class InteractionOfType
{
public:
    //! Constructor that initializes vectors.
    InteractionOfType(gmx::ArrayRef<const int>  atoms,
                      gmx::ArrayRef<const real> params,
                      const std::string&        name = "");
    /*!@{*/
    //! Access the individual elements set for the parameter.
    const int& ai() const;
    const int& aj() const;
    const int& ak() const;
    const int& al() const;
    const int& am() const;

    const real& c0() const;
    const real& c1() const;
    const real& c2() const;

    const std::string& interactionTypeName() const;
    /*!@}*/

    /*! \brief Renumbers atom Ids.
     *
     *  Enforces that ai() is less than the opposite terminal atom index,
     *  with the number depending on the interaction type.
     */
    void sortAtomIds();

    //! Set single force field parameter.
    void setForceParameter(int pos, real value);

    //! View on all atoms numbers that are actually set.
    gmx::ArrayRef<int> atoms() { return atoms_; }
    //! Const view on all atoms numbers that are actually set.
    gmx::ArrayRef<const int> atoms() const { return atoms_; }
    //! View on all of the force field parameters
    gmx::ArrayRef<const real> forceParam() const { return forceParam_; }
    //! View on all of the force field parameters
    gmx::ArrayRef<real> forceParam() { return forceParam_; }

private:
    //! Return if we have a bond parameter, means two atoms right now.
    bool isBond() const { return atoms_.size() == 2; }
    //! Return if we have an angle parameter, means three atoms right now.
    bool isAngle() const { return atoms_.size() == 3; }
    //! Return if we have a dihedral parameter, means four atoms right now.
    bool isDihedral() const { return atoms_.size() == 4; }
    //! Return if we have a cmap parameter, means five atoms right now.
    bool isCmap() const { return atoms_.size() == 5; }
    //! Enforce that atom id ai() is less than aj().
    void sortBondAtomIds();
    //! Enforce that atom id ai() is less than ak(). Does not change aj().
    void sortAngleAtomIds();
    /*! \brief Enforce order of atoms in dihedral.
     *
     * Changes atom order if needed to enforce that ai() is less than al().
     * If ai() and al() are swapped, aj() and ak() are swapped as well,
     * independent of their previous order.
     */
    void sortDihedralAtomIds();
    //! The atom list (eg. bonds: particle, i = atoms[0] (ai), j = atoms[1] (aj))
    std::vector<int> atoms_;
    //! Force parameters (eg. b0 = forceParam[0])
    std::array<real, MAXFORCEPARAM> forceParam_;
    //! Used with forcefields whose .rtp files name the interaction types (e.g. GROMOS), rather than look them up from the atom names.
    std::string interactionTypeName_;
};

/*! \libinternal \brief
 * A set of interactions of a given type
 * (found in the enumeration in ifunc.h), complete with
 * atom indices and force field function parameters.
 *
 * This is used for containing the data obtained from the
 * lists of interactions of a given type in a [moleculetype]
 * topology file definition.
 */
struct InteractionsOfType
{ // NOLINT (clang-analyzer-optin.performance.Padding)
    //! The different parameters in the system.
    std::vector<InteractionOfType> interactionTypes;
    //! CMAP grid spacing.
    int cmapGridSpacing_ = -1;
    //! Number of CMAP dihedral angle pairs.
    int numCmaps_ = -1;
    //! CMAP grid data.
    std::vector<real> cmap;
    //! The five atomtypes followed by a number that identifies the type.
    std::vector<int> cmapAtomTypes;

    //! Number of parameters.
    size_t size() const { return interactionTypes.size(); }
    //! Elements in cmap grid data.
    std::size_t ncmap() const { return cmap.size(); }
    //! Number of elements in cmapAtomTypes.
    std::size_t nct() const { return cmapAtomTypes.size(); }
};

struct t_excls
{
    int  nr; /* The number of exclusions             */
    int* e;  /* The excluded atoms                   */
};


/*! \libinternal \brief
 * Holds the molecule information during preprocessing.
 */
struct MoleculeInformation
{
    //! Name of the molecule.
    char** name = nullptr;
    //! Number of exclusions per atom.
    int nrexcl = 0;
    //! Has the mol been processed.
    bool bProcessed = false;
    //! Atoms in the moelcule.
    t_atoms atoms;
    //! Molecules separated in datastructure.
    t_block mols;
    //! Exclusions in the molecule.
    gmx::ListOfLists<int> excls;
    //! Interactions of a defined type.
    std::array<InteractionsOfType, F_NRE> interactions;

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
