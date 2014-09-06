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

#ifndef GMX_GMXPREPROCESS_PDB2TOP_H
#define GMX_GMXPREPROCESS_PDB2TOP_H

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* this *MUST* correspond to array in pdb2top.c */
enum {
    ehisA, ehisB, ehisH, ehis1, ehisNR
};
extern const char *hh[ehisNR];

typedef struct {
    int   res1, res2;
    char *a1, *a2;
} t_ssbond;

void choose_ff(const char *ffsel,
               char *forcefield, int ff_maxlen,
               char *ffdir, int ffdir_maxlen);
/* Find force fields in the current and libdirs and choose an ff.
 * If ffsel!=NULL: search for ffsel.
 * If ffsel==NULL: interactive selection.
 */

void choose_watermodel(const char *wmsel, const char *ffdir,
                       char **watermodel);
/* Choose, possibly interactively, which water model to include,
 * based on the wmsel command line option choice and watermodels.dat
 * in ffdir.
 */

void get_hackblocks_rtp(t_hackblock **hb, t_restp **restp,
                        int nrtp, t_restp rtp[],
                        int nres, t_resinfo *resinfo,
                        int nterpairs,
                        t_hackblock **ntdb, t_hackblock **ctdb,
                        int *rn, int *rc);
/* Get the database entries for the nres residues in resinfo
 * and store them in restp and hb.
 */

void match_atomnames_with_rtp(t_restp restp[], t_hackblock hb[],
                              t_atoms *pdba, rvec *x,
                              gmx_bool bVerbose);
/* Check if atom in pdba need to be deleted of renamed due to tdb or hdb.
 * If renaming involves atoms added wrt to the rtp database,
 * add these atoms to restp.
 */

void print_top_comment(FILE *out, const char *filename, const char *ffdir, gmx_bool bITP);

void print_top_header(FILE *out, const char *filename, gmx_bool bITP,
                      const char *ffdir, real mHmult);

void print_top_mols(FILE *out,
                    const char *title, const char *ffdir, const char *water,
                    int nincl, char **incls,
                    int nmol, t_mols *mols);

void write_top(FILE *out, char *pr, char *molname,
               t_atoms *at, gmx_bool bRTPresname,
               int bts[], t_params plist[], t_excls excls[],
               gpp_atomtype_t atype, int *cgnr, int nrexcl);
/* NOTE: nrexcl is not the size of *excl! */


void pdb2top(FILE *top_file, char *posre_fn, char *molname,
             t_atoms *atoms, rvec **x,
             gpp_atomtype_t atype, struct t_symtab *tab,
             int nrtp, t_restp rtp[],
             t_restp *restp, t_hackblock *hb,
             gmx_bool bAllowMissing,
             gmx_bool bVsites, gmx_bool bVsiteAromatics,
             const char *ffdir,
             real mHmult,
             int nssbonds, t_ssbond ssbonds[],
             real long_bond_dist, real short_bond_dist,
             gmx_bool bDeuterate, gmx_bool bChargeGroups, gmx_bool bCmap,
             gmx_bool bRenumRes, gmx_bool bRTPresname);
/* Create a topology ! */

void print_sums(t_atoms *atoms, gmx_bool bSystem);

#ifdef __cplusplus
}
#endif

#endif
