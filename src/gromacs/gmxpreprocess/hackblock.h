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

#ifndef GMX_GMXPREPROCESS_HACKBLOCK_H
#define GMX_GMXPREPROCESS_HACKBLOCK_H

#include <cstdio>

#include <vector>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"

struct t_atom;

/* Used for reading .rtp/.tdb */
/* ebtsBONDS must be the first, new types can be added to the end */
/* these *MUST* correspond to the arrays in hackblock.c */
enum {
    ebtsBONDS, ebtsANGLES, ebtsPDIHS, ebtsIDIHS, ebtsEXCLS, ebtsCMAP, ebtsNR
};
extern const char *btsNames[ebtsNR];
extern const int   btsNiatoms[ebtsNR];

/* if changing any of these structs, make sure that all of the
   free/clear/copy/merge_t_* functions stay updated */

/* BONDEDS */
struct t_rbonded
{
    std::array<std::string ,MAXATOMLIST> a; /* atom names */
    std::string s;   /* optional define string which gets copied from
                                          .rtp/.tdb to .top and will be parsed by cpp
                                          during grompp */
    bool     match = false;            /* boolean to mark that the entry has been found */
    const char*   ai() const { return a[0].c_str(); }
    const char*   aj() const { return a[1].c_str(); }
    const char*   ak() const { return a[2].c_str(); }
    const char*   al() const { return a[3].c_str(); }
    const char*   am() const { return a[4].c_str(); }
};

struct t_rbondeds
{
    int                    type; /* The type of bonded interaction */
    std::vector<t_rbonded> b;    /* bondeds */
    int                    nb() const { return b.size(); }
};

/* RESIDUES (rtp) */
struct t_restp
{
    std::string resname;
    /* The base file name this rtp entry was read from */
    std::string filebase;
    /* atom data */
    std::vector<t_atom>  atom;
    std::vector<char **> atomname;
    std::vector<int>     cgnr;

    int                  natom() const { return atom.size(); }

    /* Bonded interaction setup */
    bool          bKeepAllGeneratedDihedrals    = false;
    int           nrexcl                        = -1;
    bool          bGenerateHH14Interactions     = false;
    bool          bRemoveDihedralIfWithImproper = false;
    /* list of bonded interactions to add */
    std::array<t_rbondeds, ebtsNR>    rb;
};

/* Block to hack residues */
struct t_hack
{
    std::string oname;                        /* Old name                   */
    std::string nname;                        /* New name                   */
    /* the type of hack depends on the setting of oname and nname:
     * if oname==NULL                we're adding, must have tp>0 also!
     * if oname!=NULL && nname==NULL we're deleting
     * if oname!=NULL && nname!=NULL we're replacing
     */
    std::vector<t_atom> atom;                        /* New atom data              */
    int                 cgnr            = -1;        /* chargegroup number. if not read will be NOTSET */
    int                 tp              = 0;         /* Type of attachment (1..11) */
    int                 nctl            = -1;        /* How many control atoms there are */
    std::array<std::string, 4> a; /* Control atoms i,j,k,l	  */
    bool                bAlreadyPresent = false;
    bool                bXSet           = false;
    bool                bIsNewOrDelete  = false; /* Hack to set nr of atoms to add or remove if not replaced */
    int                 nrAddRemove     = -1;
    rvec                newx;      /* calculated new position    */
    int                 newi = -1; /* new atom index number (after additions) */
    int                 nr() const { if (bIsNewOrDelete) return nrAddRemove; else return atom.size(); }
    const char*              ai() const { return a[0].c_str(); }
    const char*              aj() const { return a[1].c_str(); }
    const char*              ak() const { return a[2].c_str(); }
    const char*              al() const { return a[3].c_str(); }
};

struct t_hackblock
{
    std::string name; /* Name of hack block (residue or terminus) */
    std::string filebase; /* The base file name this entry was read from */
    std::vector<t_hack> hack;               /* Hack list                                */
    /* list of bonded interactions to add */
    std::array<t_rbondeds, ebtsNR> rb;

    int                 nhack() const { return hack.size(); }
};

void free_t_restp(gmx::ArrayRef<t_restp> rtp);
/* free the whole datastructure */

void clear_t_hackblock(t_hackblock *hb);
void clear_t_hack(t_hack *hack);
/* reset struct */

bool merge_t_bondeds(const std::array<t_rbondeds, ebtsNR> &s, std::array<t_rbondeds, ebtsNR> *d,
                     bool bMin, bool bPlus);
/* add s[].b[] to d[].b[]
 * If bMin==TRUE, don't copy bondeds with atoms starting with '-'
 * If bPlus==TRUE, don't copy bondeds with atoms starting with '+'
 * Returns if bonds were removed at the termini.
 */

void copy_t_restp(const t_restp *s, t_restp *d);
void copy_t_hack(const t_hack &s, t_hack *d);
void copy_t_hackblock(const t_hackblock &s, t_hackblock *d);
/* make copy of whole datastructure */

void merge_hacks_lo(gmx::ArrayRef<const t_hack> s, std::vector<t_hack> *d);
/* add s[] to *d[] */

void merge_hacks(const t_hackblock &s, t_hackblock *d);
/* add s->hacks[] to d->hacks[] */

void merge_t_hackblock(const t_hackblock &s, t_hackblock *d);
/* add s->hacks[] and s->rb[] to d*/

void dump_hb(FILE *out, int nres, t_hackblock hb[]);
/* print out whole datastructure */

#endif
