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

#ifndef GMX_GMXPREPROCESS_PDB2TOP_H
#define GMX_GMXPREPROCESS_PDB2TOP_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

class PreprocessingAtomTypes;
struct t_atoms;
struct t_excls;
struct MoleculePatchDatabase;
struct t_mols;
struct InteractionsOfType;
struct t_resinfo;
struct PreprocessResidue;
struct DisulfideBond;
struct t_symtab;

/* this *MUST* correspond to array in pdb2top.c */
enum
{
    ehisA,
    ehisB,
    ehisH,
    ehis1,
    ehisNR
};
extern const char* hh[ehisNR];

void choose_ff(const char* ffsel, char* forcefield, int ff_maxlen, char* ffdir, int ffdir_maxlen);
/* Find force fields in the current and libdirs and choose an ff.
 * If ffsel!=NULL: search for ffsel.
 * If ffsel==NULL: interactive selection.
 */

void choose_watermodel(const char* wmsel, const char* ffdir, char** watermodel);
/* Choose, possibly interactively, which water model to include,
 * based on the wmsel command line option choice and watermodels.dat
 * in ffdir.
 */

void get_hackblocks_rtp(std::vector<MoleculePatchDatabase>*    globalPatches,
                        std::vector<PreprocessResidue>*        usedPpResidues,
                        gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
                        int                                    nres,
                        t_resinfo*                             resinfo,
                        int                                    nterpairs,
                        t_symtab*                              symtab,
                        gmx::ArrayRef<MoleculePatchDatabase*>  ntdb,
                        gmx::ArrayRef<MoleculePatchDatabase*>  ctdb,
                        gmx::ArrayRef<const int>               rn,
                        gmx::ArrayRef<const int>               rc,
                        bool                                   bAllowMissing);
/* Get the database entries for the nres residues in resinfo
 * and store them in restp and hb.
 */

void match_atomnames_with_rtp(gmx::ArrayRef<PreprocessResidue>     usedPpResidues,
                              gmx::ArrayRef<MoleculePatchDatabase> globalPatches,
                              t_atoms*                             pdba,
                              t_symtab*                            symtab,
                              gmx::ArrayRef<gmx::RVec>             x,
                              bool                                 bVerbose);
/* Check if atom in pdba need to be deleted of renamed due to tdb or hdb.
 * If renaming involves atoms added wrt to the rtp database,
 * add these atoms to restp.
 */

void print_top_comment(FILE* out, const char* filename, const char* ffdir, bool bITP);

void print_top_header(FILE* out, const char* filename, bool bITP, const char* ffdir, real mHmult);

void print_top_mols(FILE*                            out,
                    const char*                      title,
                    const char*                      ffdir,
                    const char*                      water,
                    gmx::ArrayRef<const std::string> incls,
                    gmx::ArrayRef<const t_mols>      mols);

void write_top(FILE*                                   out,
               const char*                             pr,
               const char*                             molname,
               t_atoms*                                at,
               bool                                    bRTPresname,
               int                                     bts[],
               gmx::ArrayRef<const InteractionsOfType> plist,
               t_excls                                 excls[],
               PreprocessingAtomTypes*                 atype,
               int*                                    cgnr,
               int                                     nrexcl);
/* NOTE: nrexcl is not the size of *excl! */

void pdb2top(FILE*                                  top_file,
             const char*                            posre_fn,
             const char*                            molname,
             t_atoms*                               atoms,
             std::vector<gmx::RVec>*                x,
             PreprocessingAtomTypes*                atype,
             t_symtab*                              tab,
             gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
             gmx::ArrayRef<PreprocessResidue>       usedPpResidues,
             gmx::ArrayRef<MoleculePatchDatabase>   globalPatches,
             bool                                   bAllowMissing,
             bool                                   bVsites,
             bool                                   bVsiteAromatics,
             const char*                            ffdir,
             real                                   mHmult,
             gmx::ArrayRef<const DisulfideBond>     ssbonds,
             real                                   long_bond_dist,
             real                                   short_bond_dist,
             bool                                   bDeuterate,
             bool                                   bChargeGroups,
             bool                                   bCmap,
             bool                                   bRenumRes,
             bool                                   bRTPresname);
/* Create a topology ! */

void print_sums(const t_atoms* atoms, bool bSystem);

#endif
