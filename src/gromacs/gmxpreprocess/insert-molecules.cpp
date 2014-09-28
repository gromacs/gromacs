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
#include "gmxpre.h"

#include "insert-molecules.h"

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "gromacs/gmxpreprocess/read-conformation.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* enum for random rotations of inserted solutes */
enum {
    en_rot, en_rotXYZ, en_rotZ, en_rotNone, en_NR
};

static void generate_trial_conf(int atomCount, const rvec xin[],
                                const rvec offset, int enum_rot, gmx_rng_t rng,
                                rvec xout[])
{
    for (int i = 0; i < atomCount; i++)
    {
        copy_rvec(xin[i], xout[i]);
    }
    real alfa = 0.0, beta = 0.0, gamma = 0.0;
    switch (enum_rot)
    {
        case en_rotXYZ:
            alfa  = 2*M_PI * gmx_rng_uniform_real(rng);
            beta  = 2*M_PI * gmx_rng_uniform_real(rng);
            gamma = 2*M_PI * gmx_rng_uniform_real(rng);
            break;
        case en_rotZ:
            alfa  = beta = 0.;
            gamma = 2*M_PI * gmx_rng_uniform_real(rng);
            break;
        case en_rotNone:
            alfa = beta = gamma = 0.;
            break;
    }
    if (enum_rot == en_rotXYZ || (enum_rot == en_rotZ))
    {
        rotate_conf(atomCount, xout, NULL, alfa, beta, gamma);
    }
    for (int i = 0; i < atomCount; i++)
    {
        rvec_inc(xout[i], offset);
    }
}

static bool is_insertion_allowed(gmx::AnalysisNeighborhoodSearch *search,
                                 const real *exclusionDistances,
                                 int atomCount, const rvec *x,
                                 const real *exclusionDistances_insrt)
{
    gmx::AnalysisNeighborhoodPositions  pos(x, atomCount);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search->startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const real r1 = exclusionDistances[pair.refIndex()];
        const real r2 = exclusionDistances_insrt[pair.testIndex()];
        if (pair.distance2() < sqr(r1 + r2))
        {
            return false;
        }
    }
    return true;
}

static void merge_atoms_noalloc(t_atoms *atoms, const t_atoms *atoms_add)
{
    int resnr = 0;
    if (atoms->nr > 0)
    {
        resnr = atoms->resinfo[atoms->atom[atoms->nr-1].resind].nr;
    }
    int prevResInd = -1;
    for (int i = 0; i < atoms_add->nr; ++i)
    {
        if (atoms_add->atom[i].resind != prevResInd)
        {
            prevResInd = atoms_add->atom[i].resind;
            ++resnr;
            atoms->resinfo[atoms->nres]    = atoms_add->resinfo[prevResInd];
            atoms->resinfo[atoms->nres].nr = resnr;
            ++atoms->nres;
        }
        atoms->atom[atoms->nr]        = atoms_add->atom[i];
        atoms->atomname[atoms->nr]    = atoms_add->atomname[i];
        atoms->atom[atoms->nr].resind = atoms->nres-1;
        ++atoms->nr;
    }
}

static void insert_mols(int nmol_insrt, int ntry, int seed,
                        t_atoms *atoms, rvec **x, real **exclusionDistances,
                        t_atoms *atoms_insrt, rvec *x_insrt, real *exclusionDistances_insrt,
                        int ePBC, matrix box,
                        const char* posfn, const rvec deltaR, int enum_rot)
{
    t_pbc            pbc;
    rvec            *x_n;
    int              mol;
    int              trial;
    double         **rpos;

    const real       maxInsertRadius
        = *std::max_element(exclusionDistances_insrt,
                            exclusionDistances_insrt + atoms_insrt->nr);
    real             maxRadius = maxInsertRadius;
    if (atoms->nr > 0)
    {
        const real maxExistingRadius
            = *std::max_element(*exclusionDistances,
                                *exclusionDistances + atoms->nr);
        maxRadius = std::max(maxInsertRadius, maxExistingRadius);
    }

    // TODO: Make all of this exception-safe.
    gmx::AnalysisNeighborhood nb;
    nb.setCutoff(maxInsertRadius + maxRadius);

    gmx_rng_t        rng = gmx_rng_init(seed);
    set_pbc(&pbc, ePBC, box);

    snew(x_n, atoms_insrt->nr);

    /* With -ip, take nmol_insrt from file posfn */
    if (posfn != NULL)
    {
        int ncol;
        nmol_insrt = read_xvg(posfn, &rpos, &ncol);
        if (ncol != 3)
        {
            gmx_fatal(FARGS, "Expected 3 columns (x/y/z coordinates) in file %s\n", ncol, posfn);
        }
        fprintf(stderr, "Read %d positions from file %s\n\n", nmol_insrt, posfn);
    }

    {
        const int finalAtomCount    = atoms->nr + nmol_insrt * atoms_insrt->nr;
        const int finalResidueCount = atoms->nres + nmol_insrt * atoms_insrt->nres;
        srenew(atoms->resinfo,      finalResidueCount);
        srenew(atoms->atomname,     finalAtomCount);
        srenew(atoms->atom,         finalAtomCount);
        srenew(*x,                  finalAtomCount);
        srenew(*exclusionDistances, finalAtomCount);
    }

    trial = mol = 0;
    while ((mol < nmol_insrt) && (trial < ntry*nmol_insrt))
    {
        fprintf(stderr, "\rTry %d", trial++);
        rvec offset_x;
        if (posfn == NULL)
        {
            /* insert at random positions */
            offset_x[XX] = box[XX][XX] * gmx_rng_uniform_real(rng);
            offset_x[YY] = box[YY][YY] * gmx_rng_uniform_real(rng);
            offset_x[ZZ] = box[ZZ][ZZ] * gmx_rng_uniform_real(rng);
        }
        else
        {
            /* Insert at positions taken from option -ip file */
            offset_x[XX] = rpos[XX][mol] + deltaR[XX]*(2 * gmx_rng_uniform_real(rng)-1);
            offset_x[YY] = rpos[YY][mol] + deltaR[YY]*(2 * gmx_rng_uniform_real(rng)-1);
            offset_x[ZZ] = rpos[ZZ][mol] + deltaR[ZZ]*(2 * gmx_rng_uniform_real(rng)-1);
        }
        generate_trial_conf(atoms_insrt->nr, x_insrt, offset_x, enum_rot, rng,
                            x_n);
        gmx::AnalysisNeighborhoodPositions pos(*x, atoms->nr);
        gmx::AnalysisNeighborhoodSearch    search = nb.initSearch(&pbc, pos);
        if (is_insertion_allowed(&search, *exclusionDistances, atoms_insrt->nr,
                                 x_n, exclusionDistances_insrt))
        {
            const int firstIndex = atoms->nr;
            for (int i = 0; i < atoms_insrt->nr; ++i)
            {
                copy_rvec(x_n[i], (*x)[firstIndex + i]);
                (*exclusionDistances)[firstIndex + i] = exclusionDistances_insrt[i];
            }
            merge_atoms_noalloc(atoms, atoms_insrt);
            ++mol;
            fprintf(stderr, " success (now %d atoms)!\n", atoms->nr);
        }
    }
    gmx_rng_destroy(rng);
    srenew(atoms->resinfo,  atoms->nres);
    srenew(atoms->atomname, atoms->nr);
    srenew(atoms->atom,     atoms->nr);
    srenew(*x,              atoms->nr);
    srenew(*exclusionDistances, atoms->nr);

    fprintf(stderr, "\n");
    /* print number of molecules added */
    fprintf(stderr, "Added %d molecules (out of %d requested) of %s\n",
            mol, nmol_insrt, *atoms_insrt->resinfo[0].name);

    sfree(x_n);
}

int gmx_insert_molecules(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] inserts [TT]-nmol[tt] copies of the system specified in",
        "the [TT]-ci[tt] input file. The insertions take place either into",
        "vacant space in the solute conformation given with [TT]-f[tt], or",
        "into an empty box given by [TT]-box[tt]. Specifying both [TT]-f[tt]",
        "and [TT]-box[tt] behaves like [TT]-f[tt], but places a new box",
        "around the solute before insertions. Any velocities present are",
        "discarded.[PAR]",

        "By default, the insertion positions are random (with initial seed",
        "specified by [TT]-seed[tt]). The program iterates until [TT]-nmol[tt]",
        "molecules have been inserted in the box. Molecules are not inserted",
        "where the distance between any existing atom and any atom of the",
        "inserted molecule is less than the sum based on the van der Waals",
        "radii of both atoms. A database ([TT]vdwradii.dat[tt]) of van der",
        "Waals radii is read by the program, and the resulting radii scaled",
        "by [TT]-scale[tt]. If radii are not found in the database, those"
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].[PAR]",

        "A total of [TT]-nmol[tt] * [TT]-try[tt] insertion attempts are made",
        "before giving up. Increase [TT]-try[tt] if you have several small",
        "holes to fill. Option [TT]-rot[tt] specifies whether the insertion",
        "molecules are randomly oriented before insertion attempts.[PAR]",

        "Alternatively, the molecules can be inserted only at positions defined in",
        "positions.dat ([TT]-ip[tt]). That file should have 3 columns (x,y,z),",
        "that give the displacements compared to the input molecule position",
        "([TT]-ci[tt]). Hence, if that file should contain the absolute",
        "positions, the molecule must be centered on (0,0,0) before using",
        "[THISMODULE] (e.g. from [gmx-editconf] [TT]-center[tt]).",
        "Comments in that file starting with # are ignored. Option [TT]-dr[tt]",
        "defines the maximally allowed displacements during insertial trials.",
        "[TT]-try[tt] and [TT]-rot[tt] work as in the default mode (see above).",
        "[PAR]",
    };

    const char *bugs[] = {
        "Molecules must be whole in the initial configurations."
    };

    /* parameter data */
    gmx_bool       bProt, bBox;
    real          *exclusionDistances       = NULL;
    real          *exclusionDistances_insrt = NULL;
    gmx_atomprop_t aps;

    /* protein configuration data */
    char          *title = NULL;
    t_atoms       *atoms, *atoms_insrt;
    rvec          *x    = NULL, *x_insrt = NULL;
    int            ePBC = -1;
    matrix         box;

    t_filenm       fnm[] = {
        { efSTX, "-f", "protein", ffOPTRD },
        { efSTX, "-ci", "insert",  ffREAD},
        { efDAT, "-ip", "positions",  ffOPTRD},
        { efSTO, NULL,  NULL,      ffWRITE},
    };
#define NFILE asize(fnm)

    static int      nmol_ins               = 0, nmol_try = 10, seed = 1997, enum_rot;
    static real     defaultDistance        = 0.105, scaleFactor = 0.57;
    static rvec     new_box                = {0.0, 0.0, 0.0}, deltaR = {0.0, 0.0, 0.0};
    output_env_t    oenv;
    const char     *enum_rot_string[] = {NULL, "xyz", "z", "none", NULL};
    t_pargs         pa[]              = {
        { "-box",    FALSE, etRVEC, {new_box},
          "Box size (in nm)" },
        { "-nmol",   FALSE, etINT, {&nmol_ins},
          "Number of extra molecules to insert" },
        { "-try",    FALSE, etINT, {&nmol_try},
          "Try inserting [TT]-nmol[tt] times [TT]-try[tt] times" },
        { "-seed",   FALSE, etINT, {&seed},
          "Random generator seed"},
        { "-radius",   FALSE, etREAL, {&defaultDistance},
          "Default van der Waals distance"},
        { "-scale", FALSE, etREAL, {&scaleFactor},
          "Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 g/l for proteins in water." },
        { "-dr",    FALSE, etRVEC, {deltaR},
          "Allowed displacement in x/y/z from positions in [TT]-ip[tt] file" },
        { "-rot", FALSE,  etENUM, {enum_rot_string},
          "rotate inserted molecules randomly" }
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    bProt     = opt2bSet("-f", NFILE, fnm);
    bBox      = opt2parg_bSet("-box", asize(pa), pa);
    enum_rot  = nenum(enum_rot_string);

    /* check input */
    if (nmol_ins <= 0 && !opt2bSet("-ip", NFILE, fnm))
    {
        gmx_fatal(FARGS, "Either -nmol must be larger than 0, "
                  "or positions must be given with -ip");
    }
    if (!bProt && !bBox)
    {
        gmx_fatal(FARGS, "When no solute (-f) is specified, "
                  "a box size (-box) must be specified");
    }

    aps = gmx_atomprop_init();

    snew(atoms, 1);
    init_t_atoms(atoms, 0, FALSE);
    if (bProt)
    {
        /* Generate a solute configuration */
        const char *conf_prot = opt2fn("-f", NFILE, fnm);
        title                 = readConformation(conf_prot, atoms, &x, NULL,
                                                 &ePBC, box);
        exclusionDistances = makeExclusionDistances(atoms, aps, defaultDistance, scaleFactor);
        if (atoms->nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", conf_prot);
            sfree(title);
            title = NULL;
            bProt = FALSE;
        }
    }
    if (bBox)
    {
        ePBC = epbcXYZ;
        clear_mat(box);
        box[XX][XX] = new_box[XX];
        box[YY][YY] = new_box[YY];
        box[ZZ][ZZ] = new_box[ZZ];
    }
    if (det(box) == 0)
    {
        gmx_fatal(FARGS, "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }
    snew(atoms_insrt, 1);
    init_t_atoms(atoms_insrt, 0, FALSE);
    {
        int         ePBC_dummy;
        matrix      box_dummy;
        const char *conf_insrt = opt2fn("-ci", NFILE, fnm);
        char       *title_ins
            = readConformation(conf_insrt, atoms_insrt, &x_insrt, NULL,
                               &ePBC_dummy, box_dummy);
        if (atoms_insrt->nr == 0)
        {
            gmx_fatal(FARGS, "No molecule in %s, please check your input", conf_insrt);
        }
        sfree(title);
        title                    = title_ins;
        exclusionDistances_insrt = makeExclusionDistances(atoms_insrt, aps, defaultDistance, scaleFactor);
    }

    gmx_atomprop_destroy(aps);

    /* add nmol_ins molecules of atoms_ins
       in random orientation at random place */
    insert_mols(nmol_ins, nmol_try, seed,
                atoms, &x, &exclusionDistances,
                atoms_insrt, x_insrt, exclusionDistances_insrt,
                ePBC, box,
                opt2fn_null("-ip", NFILE, fnm), deltaR, enum_rot);

    /* write new configuration to file confout */
    const char *confout = ftp2fn(efSTO, NFILE, fnm);
    fprintf(stderr, "Writing generated configuration to %s\n", confout);
    write_sto_conf(confout, title, atoms, x, NULL, ePBC, box);

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n",
            atoms->nr, atoms->nres);

    sfree(exclusionDistances);
    sfree(exclusionDistances_insrt);
    sfree(x);
    done_atom(atoms);
    done_atom(atoms_insrt);
    sfree(atoms);
    sfree(atoms_insrt);
    sfree(title);

    return 0;
}
