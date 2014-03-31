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
#include "insert-molecules.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "sysstuff.h"
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/utilities.h"
#include "gromacs/fileio/confio.h"
#include "macros.h"
#include "gromacs/random/random.h"
#include "gromacs/fileio/futil.h"
#include "atomprop.h"
#include "names.h"
#include "vec.h"
#include "gmx_fatal.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "addconf.h"
#include "read-conformation.h"
#include "pbc.h"
#include "xvgr.h"

static gmx_bool in_box(t_pbc *pbc, rvec x)
{
    rvec box_center, dx;
    int  shift;

    /* pbc_dx_aiuc only works correctly with the rectangular box center */
    calc_box_center(ecenterRECT, pbc->box, box_center);

    shift = pbc_dx_aiuc(pbc, x, box_center, dx);

    return (shift == CENTRAL);
}

/* This is a (maybe) slow workaround to avoid the neighbor searching in addconf.c, which
 * leaks memory (May 2012). The function could be deleted as soon as the memory leaks
 * there are fixed.
 * However, when inserting a small molecule in a system containing not too many atoms,
 * allPairsDistOk is probably even faster than the other code.
 */
static gmx_bool
allPairsDistOk(t_atoms *atoms, rvec *x, real *exclusionDistances,
               int ePBC, matrix box,
               t_atoms *atoms_insrt, rvec *x_n, real *exclusionDistances_insrt)
{
    int   i, j;
    rvec  dx;
    real  n2, r2;
    t_pbc pbc;

    set_pbc(&pbc, ePBC, box);
    for (i = 0; i < atoms->nr; i++)
    {
        for (j = 0; j < atoms_insrt->nr; j++)
        {
            pbc_dx(&pbc, x[i], x_n[j], dx);
            n2 = norm2(dx);
            r2 = sqr(exclusionDistances[i]+exclusionDistances_insrt[j]);
            if (n2 < r2)
            {
                return FALSE;
            }
        }
    }
    return TRUE;
}

/* enum for random rotations of inserted solutes */
enum {
    en_rot, en_rotXYZ, en_rotZ, en_rotNone, en_NR
};

static char *insert_mols(const char *mol_insrt, int nmol_insrt, int ntry, int seed,
                         t_atoms *atoms, rvec **x, real **exclusionDistances, int ePBC, matrix box,
                         gmx_atomprop_t aps,
                         real defaultDistance, real scaleFactor, real rshell,
                         const output_env_t oenv,
                         const char* posfn, const rvec deltaR, int enum_rot,
                         gmx_bool bCheckAllPairDist)
{
    t_pbc            pbc;
    static  char    *title_insrt;
    t_atoms          atoms_insrt;
    rvec            *x_insrt, *x_n;
    real            *exclusionDistances_insrt;
    int              ePBC_insrt;
    matrix           box_insrt;
    int              i, mol, onr, ncol;
    real             alfa = 0., beta = 0., gamma = 0.;
    rvec             offset_x;
    int              trial;
    double         **rpos;
    gmx_rng_t        rng;

    rng = gmx_rng_init(seed);
    set_pbc(&pbc, ePBC, box);

    /* read number of atoms of insert molecules */
    {
        int natoms;
        get_stx_coordnum(mol_insrt, &natoms);
        if (natoms == 0)
        {
            gmx_fatal(FARGS, "No molecule in %s, please check your input\n", mol_insrt);
        }
        init_t_atoms(&atoms_insrt, natoms, FALSE);
    }
    /* allocate memory for atom coordinates of insert molecules */
    snew(x_insrt, atoms_insrt.nr);
    snew(atoms_insrt.resinfo, atoms_insrt.nr);
    snew(atoms_insrt.atomname, atoms_insrt.nr);
    snew(atoms_insrt.atom, atoms_insrt.nr);
    atoms_insrt.pdbinfo = NULL;
    snew(x_n, atoms_insrt.nr);
    snew(title_insrt, STRLEN);

    /* read residue number, residue names, atomnames, coordinates etc. */
    fprintf(stderr, "Reading molecule configuration \n");
    read_stx_conf(mol_insrt, title_insrt, &atoms_insrt, x_insrt, NULL,
                  &ePBC_insrt, box_insrt);
    fprintf(stderr, "%s\nContaining %d atoms in %d residue\n",
            title_insrt, atoms_insrt.nr, atoms_insrt.nres);
    srenew(atoms_insrt.resinfo, atoms_insrt.nres);

    /* initialise distance arrays for inserted molecules */
    exclusionDistances_insrt = makeExclusionDistances(&atoms_insrt, aps, defaultDistance, scaleFactor);

    /* With -ip, take nmol_insrt from file posfn */
    if (posfn != NULL)
    {
        nmol_insrt = read_xvg(posfn, &rpos, &ncol);
        if (ncol != 3)
        {
            gmx_fatal(FARGS, "Expected 3 columns (x/y/z coordinates) in file %s\n", ncol, posfn);
        }
        fprintf(stderr, "Read %d positions from file %s\n\n", nmol_insrt, posfn);
    }

    srenew(atoms->resinfo, (atoms->nres+nmol_insrt*atoms_insrt.nres));
    srenew(atoms->atomname, (atoms->nr+atoms_insrt.nr*nmol_insrt));
    srenew(atoms->atom, (atoms->nr+atoms_insrt.nr*nmol_insrt));
    srenew(*x, (atoms->nr+atoms_insrt.nr*nmol_insrt));
    srenew(*exclusionDistances, (atoms->nr+atoms_insrt.nr*nmol_insrt));

    trial = mol = 0;
    while ((mol < nmol_insrt) && (trial < ntry*nmol_insrt))
    {
        fprintf(stderr, "\rTry %d", trial++);
        for (i = 0; (i < atoms_insrt.nr); i++)
        {
            copy_rvec(x_insrt[i], x_n[i]);
        }
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
            rotate_conf(atoms_insrt.nr, x_n, NULL, alfa, beta, gamma);
        }
        if (posfn == NULL)
        {
            /* insert at random positions */
            offset_x[XX] = box[XX][XX] * gmx_rng_uniform_real(rng);
            offset_x[YY] = box[YY][YY] * gmx_rng_uniform_real(rng);
            offset_x[ZZ] = box[ZZ][ZZ] * gmx_rng_uniform_real(rng);
            make_new_box(atoms_insrt.nr, x_n, box_insrt, offset_x, TRUE);
            if (!in_box(&pbc, x_n[0]) || !in_box(&pbc, x_n[atoms_insrt.nr-1]))
            {
                continue;
            }
        }
        else
        {
            /* Insert at positions taken from option -ip file */
            offset_x[XX] = rpos[XX][mol] + deltaR[XX]*(2 * gmx_rng_uniform_real(rng)-1);
            offset_x[YY] = rpos[YY][mol] + deltaR[YY]*(2 * gmx_rng_uniform_real(rng)-1);
            offset_x[ZZ] = rpos[ZZ][mol] + deltaR[ZZ]*(2 * gmx_rng_uniform_real(rng)-1);
            for (i = 0; i < atoms_insrt.nr; i++)
            {
                rvec_inc(x_n[i], offset_x);
            }
        }

        onr = atoms->nr;

        /* This is a (maybe) slow workaround to avoid too many calls of add_conf, which
         * leaks memory (status May 2012). If the momory leaks in add_conf() are fixed,
         * this check could be removed. Note, however, that allPairsDistOk is probably
         * even faster than add_conf() when inserting a small molecule into a moderately
         * small system.
         */
        if (bCheckAllPairDist && !allPairsDistOk(atoms, *x, *exclusionDistances, ePBC, box, &atoms_insrt, x_n, exclusionDistances_insrt))
        {
            continue;
        }

        add_conf(atoms, x, NULL, exclusionDistances, FALSE, ePBC, box, TRUE,
                 &atoms_insrt, x_n, NULL, exclusionDistances_insrt, FALSE, rshell, 0, oenv);

        if (atoms->nr == (atoms_insrt.nr+onr))
        {
            mol++;
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
            mol, nmol_insrt, *atoms_insrt.resinfo[0].name);

    sfree(x_n);
    sfree(exclusionDistances_insrt);
    done_atom(&atoms_insrt);

    return title_insrt;
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
        "Molecules must be whole in the initial configurations.",
        "Many repeated neighbor searchings with -ci blows up the allocated memory. "
        "Option -allpair avoids this using all-to-all distance checks (slow for large systems)"
    };

    /* parameter data */
    gmx_bool       bProt, bBox;
    const char    *conf_prot, *confout;
    real          *exclusionDistances = NULL;
    char          *title_ins          = NULL;
    gmx_atomprop_t aps;

    /* protein configuration data */
    char          *title = NULL;
    t_atoms       *atoms;
    rvec          *x    = NULL;
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
    static gmx_bool bCheckAllPairDist      = FALSE;
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
          "rotate inserted molecules randomly" },
        { "-allpair",    FALSE, etBOOL, {&bCheckAllPairDist},
          "Avoid momory leaks during neighbor searching with option -ci. May be slow for large systems." },
    };

    if (!parse_common_args(&argc, argv, PCA_BE_NICE, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    bProt     = opt2bSet("-f", NFILE, fnm);
    bBox      = opt2parg_bSet("-box", asize(pa), pa);
    enum_rot  = nenum(enum_rot_string);

    /* check input */
    const char *insertionMoleculeFileName = opt2fn("-ci", NFILE, fnm);
    if (!gmx_fexist(insertionMoleculeFileName))
    {
        gmx_fatal(FARGS,
                  "A molecule conformation to insert is required in -ci. %s was not found!",
                  insertionMoleculeFileName);
    }
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
        conf_prot = opt2fn("-f", NFILE, fnm);
        title     = readConformation(conf_prot, atoms, &x, NULL,
                                     &ePBC, box);
        exclusionDistances = makeExclusionDistances(atoms, aps, defaultDistance, scaleFactor);
        if (atoms->nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", conf_prot);
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

    /* add nmol_ins molecules of atoms_ins
       in random orientation at random place */
    title_ins = insert_mols(insertionMoleculeFileName, nmol_ins, nmol_try, seed,
                            atoms, &x, &exclusionDistances, ePBC, box, aps,
                            defaultDistance, scaleFactor, 0,
                            oenv, opt2fn_null("-ip", NFILE, fnm), deltaR, enum_rot,
                            bCheckAllPairDist);

    /* write new configuration to file confout */
    confout = ftp2fn(efSTO, NFILE, fnm);
    fprintf(stderr, "Writing generated configuration to %s\n", confout);
    if (bProt)
    {
        write_sto_conf(confout, title, atoms, x, NULL, ePBC, box);
        /* print box sizes and box type to stderr */
        fprintf(stderr, "%s\n", title);
    }
    else
    {
        write_sto_conf(confout, title_ins, atoms, x, NULL, ePBC, box);
    }

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n",
            atoms->nr, atoms->nres);

    gmx_atomprop_destroy(aps);
    sfree(exclusionDistances);
    sfree(x);
    done_atom(atoms);
    sfree(atoms);

    return 0;
}
