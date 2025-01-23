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

#ifndef GMX_GMXPREPROCESS_GEN_VSITE_H
#define GMX_GMXPREPROCESS_GEN_VSITE_H

#include <filesystem>
#include <optional>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

class PreprocessingAtomTypes;
struct t_atoms;
struct InteractionsOfType;
struct PreprocessResidue;
struct t_symtab;

namespace gmx
{
template<typename>
class ArrayRef;
}

/* stuff for pdb2gmx */

//! Struct for handling vsite information for vsite generation inside grompp
struct VsiteTypeAndSign
{
    //! The interaction type
    std::optional<int> ftype;
    //! Tells whether we should swap the sign of vsites that can have different orientations
    bool swapSign = false;
};

//! Turn all hydrogens that can be turned into virtual sites into virtual sites
void do_vsites(gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
               PreprocessingAtomTypes*                atype,
               t_atoms*                               at,
               t_symtab*                              symtab,
               std::vector<gmx::RVec>*                x,
               gmx::ArrayRef<InteractionsOfType>      plist,
               std::vector<VsiteTypeAndSign>*         vsiteTypeAndSign,
               real                                   mHmult,
               bool                                   bVSiteAromatics,
               const std::filesystem::path&           ffdir);

/*! \brief Optionally, change masses of hydrogens
 *
 * \param[in] psb        List of bond interations
 * \param[in] vsiteType  Array with vsites with their type and sign, only type is used
 * \param[in,out] at     Atom information, masses may be modified
 * \param[in] mHmult     Factor to multiply the masses of hydrogens with
 * \param[in] deuterate  When false, subtract the increase in hydrogen mass from the bonded heavt atom
 */
void do_h_mass(const InteractionsOfType&             psb,
               gmx::ArrayRef<const VsiteTypeAndSign> vsiteType,
               t_atoms*                              at,
               real                                  mHmult,
               bool                                  deuterate);

#endif
