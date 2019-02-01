/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018,2019, by the GROMACS development team, led by
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
/*! \file
 * \libinternal \brief
 * Methods to modify atoms during preprocessing.
 */
#ifndef GMX_GMXPREPROCESS_HACKBLOCK_H
#define GMX_GMXPREPROCESS_HACKBLOCK_H

#include <cstdio>

#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"

struct t_atom;

/*! \libinternal \brief
 * Used for reading .rtp/.tdb
 * ebtsBONDS must be the first, new types can be added to the end
 * these *MUST* correspond to the arrays in hackblock.cpp
 */
enum {
    ebtsBONDS, ebtsANGLES, ebtsPDIHS, ebtsIDIHS, ebtsEXCLS, ebtsCMAP, ebtsNR
};
//! Names for interaction type entries
extern const char *btsNames[ebtsNR];
//! Numbers for atoms in the interactions.
extern const int   btsNiatoms[ebtsNR];

/* if changing any of these structs, make sure that all of the
   free/clear/copy/merge_t_* functions stay updated */

/* BONDEDS */
struct t_rbonded
{
    char  *a[MAXATOMLIST]; /* atom names */
    char  *s;              /* optional define string which gets copied from
                              .rtp/.tdb to .top and will be parsed by cpp
                              during grompp */
    bool     match;        /* boolean to mark that the entry has been found */
    char*   &ai() { return a[0]; }
    char*   &aj() { return a[1]; }
    char*   &ak() { return a[2]; }
    char*   &al() { return a[3]; }
    char*   &am() { return a[4]; }
};

struct t_rbondeds
{
    int        type;     /* The type of bonded interaction */
    int        nb;       /* number of bondeds */
    t_rbonded *b;        /* bondeds */
};

/* RESIDUES (rtp) */
struct t_restp
{
    char         *resname;
    /* The base file name this rtp entry was read from */
    char         *filebase;
    /* atom data */
    int           natom;
    t_atom       *atom;
    char       ***atomname;
    int          *cgnr;
    /* Bonded interaction setup */
    bool          bKeepAllGeneratedDihedrals;
    int           nrexcl;
    bool          bGenerateHH14Interactions;
    bool          bRemoveDihedralIfWithImproper;
    /* list of bonded interactions to add */
    t_rbondeds    rb[ebtsNR];
};

enum class HackType
{
    Add,
    Delete,
    Replace
};

/*! \internal \brief
 * Block to modify individual residues
 */
struct HackBlock
{
    //! Number of new are deleted atoms. NOT always equal to atom.size()!
    int                        nr;
    //! Old name for entry.
    std::string                oname;
    //! New name for entry.
    std::string                nname;
    //! New atom data.
    std::vector<t_atom>        atom;
    //! Chargegroup number.
    int                        cgnr = NOTSET;
    //! Type of attachement.
    int                        tp = 0;
    //! Number of control atoms.
    int                        nctl = 0;
    //! Name of control atoms.
    std::array<std::string, 4> a;
    //! Is an atom to be hacked already present?
    bool                       bAlreadyPresent = false;
    //! Are coordinates for a new atom already set?
    bool                       bXSet = false;
    //! New position for hacked atom.
    rvec                       newx = {NOTSET};
    //! New atom index number after additions.
    int                        newi = -1;

    /*! \brief
     * Get type of hack.
     *
     * This depends on the setting of oname and nname
     * for legacy reasons. If oname is empty, we are adding,
     * if oname is set and nname is empty, an atom is deleted,
     * if both are set replacement is going on. If both are unset,
     * an error is thrown.
     */
    HackType type() const;

    //! Control atom i name.
    const char* ai() const { return a[0].c_str(); }
    //! Control atom j name.
    const char* aj() const { return a[1].c_str(); }
    //! Control atom k name.
    const char* ak() const { return a[2].c_str(); }
    //! Control atom l name.
    const char* al() const { return a[3].c_str(); }
};
/*!\libinternal \brief
 * A set of modifications to apply to atoms.
 */
struct AtomModificationBlock
{
    //! Name of block
    std::string                    name;
    //! File that entry was read from.
    std::string                    filebase;
    //! List of changes to atoms.
    std::vector<HackBlock>         hack;
    //! List of bonded interactions to add.
    std::array<t_rbondeds, ebtsNR> rb;
    //! Number of atoms to modify
    int                            nhack() const { return hack.size(); }
};

//! Free t_restp
void free_t_restp(int nrtp, t_restp **rtp);

/*!\brief
 * Clear up memory.
 *
 * \param[in] amb Datastructure to clean.
 * \todo Remove once the underlying data has been cleaned up.
 */
void freeModificationBlock(gmx::ArrayRef<AtomModificationBlock> amb);

/*! \internal \brief
 * Reset modification block.
 *
 * \param[inout] amb Block to reset.
 * \todo Remove once constructor/destructor takes care of all of this.
 */
void clearModificationBlock(AtomModificationBlock *amb);

/*!\brief
 * Add bond information in \p s to \p d.
 *
 * \param[in] s Source information to copy.
 * \param[inout] d Destination to copy to.
 * \param[in] bMin don't copy bondeds with atoms starting with '-'.
 * \param[in] bPlus don't copy bondeds with atoms starting with '+'.
 * \returns if bonds were removed at the termini.
 */
bool merge_t_bondeds(gmx::ArrayRef<const t_rbondeds> s, gmx::ArrayRef<t_rbondeds> d,
                     bool bMin, bool bPlus);

//! Copy t_restp.
void copy_t_restp(t_restp *s, t_restp *d);
/*!\brief
 * Copy all information from datastructure.
 *
 * \param[in] s Source information.
 * \param[inout] d Destination to copy to.
 */
void copyModificationBlocks(const AtomModificationBlock &s, AtomModificationBlock *d);

/*!\brief
 * Add the individual modifications in \p s to \p d.
 *
 * \param[in] s Source information.
 * \param[inout] d Destination to copy to.
 */
void mergeAtomModifications(const AtomModificationBlock &s, AtomModificationBlock *d);

//! \copydoc mergeAtomModifications
void mergeAtomAndBondModifications(const AtomModificationBlock &s, AtomModificationBlock *d);

#endif
