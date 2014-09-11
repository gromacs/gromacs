/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/topology/symtab.h"

#ifdef __cplusplus
extern "C" {
#endif

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
typedef struct {
    char  *a[MAXATOMLIST]; /* atom names */
    char  *s;              /* optional define string which gets copied from
                              .rtp/.tdb to .top and will be parsed by cpp
                              during grompp */
    gmx_bool match;        /* boolean to mark that the entry has been found */
} t_rbonded;

typedef struct {
    int        type;     /* The type of bonded interaction */
    int        nb;       /* number of bondeds */
    t_rbonded *b;        /* bondeds */
} t_rbondeds;

/* RESIDUES (rtp) */
typedef struct {
    char         *resname;
    /* The base file name this rtp entry was read from */
    char         *filebase;
    /* atom data */
    int           natom;
    t_atom       *atom;
    char       ***atomname;
    int          *cgnr;
    /* Bonded interaction setup */
    gmx_bool      bKeepAllGeneratedDihedrals;
    int           nrexcl;
    gmx_bool      bGenerateHH14Interactions;
    gmx_bool      bRemoveDihedralIfWithImproper;
    /* list of bonded interactions to add */
    t_rbondeds    rb[ebtsNR];
} t_restp;

/* Block to hack residues */
typedef struct {
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
    gmx_bool    bAlreadyPresent;
    gmx_bool    bXSet;
    rvec        newx; /* calculated new position    */
    atom_id     newi; /* new atom index number (after additions) */
} t_hack;

typedef struct {
    char      *name;     /* Name of hack block (residue or terminus) */
    char      *filebase; /* The base file name this entry was read from */
    int        nhack;    /* Number of atoms to hack                  */
    int        maxhack;  /* used for efficient srenew-ing            */
    t_hack    *hack;     /* Hack list                                */
    /* list of bonded interactions to add */
    t_rbondeds rb[ebtsNR];
} t_hackblock;

/* all libraries and other data to protonate a structure or trajectory */
typedef struct {
    gmx_bool        bInit; /* true after init; set false by init_t_protonate */
    /* force field name: */
    char            FF[10];
    /* libarary data: */
    int            *nab;
    t_hack        **ab;
    t_hackblock    *ah, *ntdb, *ctdb;
    t_hackblock   **sel_ntdb, **sel_ctdb;
    int             nah;
    t_symtab        tab;
    /* residue indices (not numbers!) of the N and C termini */
    int            *rN, *rC;
    gpp_atomtype_t  atype;
    /* protonated topology: */
    t_atoms        *patoms;
    /* unprotonated topology: */
    t_atoms        *upatoms;

} t_protonate;

typedef struct {
    char *res1, *res2;
    char *atom1, *atom2;
    char *newres1, *newres2;
    int   nbond1, nbond2;
    real  length;
} t_specbond;

t_specbond *get_specbonds(int *nspecbond);
void done_specbonds(int nsb, t_specbond sb[]);

void free_t_restp(int nrtp, t_restp **rtp);
void free_t_hack(int nh, t_hack **h);
void free_t_hackblock(int nhb, t_hackblock **hb);
/* free the whole datastructure */

void clear_t_hackblock(t_hackblock *hb);
void clear_t_hack(t_hack *hack);
/* reset struct */

gmx_bool merge_t_bondeds(t_rbondeds s[], t_rbondeds d[],
                         gmx_bool bMin, gmx_bool bPlus);
/* add s[].b[] to d[].b[]
 * If bMin==TRUE, don't copy bondeds with atoms starting with '-'
 * If bPlus==TRUE, don't copy bondeds with atoms starting with '+'
 * Returns if bonds were removed at the termini.
 */

void copy_t_restp(t_restp *s, t_restp *d);
void copy_t_hack(t_hack *s, t_hack *d);
void copy_t_hackblock(t_hackblock *s, t_hackblock *d);
/* make copy of whole datastructure */

void merge_hacks_lo(int ns, t_hack *s, int *nd, t_hack **d);
/* add s[] to *d[] */

void merge_hacks(t_hackblock *s, t_hackblock *d);
/* add s->hacks[] to d->hacks[] */

void merge_t_hackblock(t_hackblock *s, t_hackblock *d);
/* add s->hacks[] and s->rb[] to d*/

void dump_hb(FILE *out, int nres, t_hackblock hb[]);
/* print out whole datastructure */

void init_t_protonate(t_protonate *protonate);
/* initialize t_protein struct */

#ifdef __cplusplus
}
#endif

#endif
