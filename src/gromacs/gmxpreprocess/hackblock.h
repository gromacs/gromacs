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

/* Block to hack residues */
struct t_hack
{
    int      nr;      /* Number of atoms to hack    */
    char    *oname;   /* Old name                   */
    char    *nname;   /* New name                   */
    /* the type of hack depends on the setting of oname and nname:
     * if oname==NULL                we're adding, must have tp>0 also!
     * if oname!=NULL && nname==NULL we're deleting
     * if oname!=NULL && nname!=NULL we're replacing
     */
    t_atom     *atom; /* New atom data              */
    int         cgnr; /* chargegroup number. if not read will be NOTSET */
    int         tp;   /* Type of attachment (1..11) */
    int         nctl; /* How many control atoms there are */
    char       *a[4]; /* Control atoms i,j,k,l	  */
    bool        bAlreadyPresent;
    bool        bXSet;
    rvec        newx; /* calculated new position    */
    int         newi; /* new atom index number (after additions) */
    const char* ai() const { return a[0]; }
    const char* aj() const { return a[1]; }
    const char* ak() const { return a[2]; }
    const char* al() const { return a[3]; }
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
    //! Number of atoms to modify
    int                            nhack = 0;
    //! Max number to modify.
    int                            maxhack = 0;
    //! List of changes to atoms.
    t_hack                        *hack = nullptr;
    //! List of bonded interactions to add.
    std::array<t_rbondeds, ebtsNR> rb;
};

//! Free t_restp
void free_t_restp(int nrtp, t_restp **rtp);
//! Free t_hack
void free_t_hack(int nh, t_hack **h);

/*!\brief
 * Clear up memory.
 *
 * \param[in] amb Datastructure to clean.
 * \todo Remove once the underlying data has been cleaned up.
 */
void freeModificationBlock(gmx::ArrayRef<AtomModificationBlock> amb);

/*!\brief
 * Reset the datastructure to initial state.
 * \param[inout] amb Datastructure to reset.
 * \todo Remove once underlying data has been cleaned up.
 */
void clearModificationBlock(AtomModificationBlock *amb);

//! Reset t_hack
void clear_t_hack(t_hack *hack);

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
//! Copy t_hack.
void copy_t_hack(const t_hack *s, t_hack *d);

/*!\brief
 * Copy all information from datastructure.
 *
 * \param[in] s Source information.
 * \param[inout] d Destination to copy to.
 */
void copyModificationBlocks(const AtomModificationBlock &s, AtomModificationBlock *d);

/*!\brief
 * Add information in \p s to \p d.
 *
 * \param[in] ns Number of source elements to add
 * \param[in] s Source information to add.
 * \param[in] nd Counter for destination elements.
 * \param[inout] d Destination to add information to.
 */
void merge_hacks_lo(int ns, const t_hack *s, int *nd, t_hack **d);

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
