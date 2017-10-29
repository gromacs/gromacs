/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "x2top.h"

#include <cmath>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/hackblock.h"
#include "gromacs/gmxpreprocess/nm2type.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

char atp[7] = "HCNOSX";
#define NATP (asize(atp)-1)

double blen[NATP][NATP] = {
    {  0.00,  0.108, 0.105, 0.10, 0.10, 0.10 },
    {  0.108, 0.15,  0.14,  0.14, 0.16, 0.14 },
    {  0.105, 0.14,  0.14,  0.14, 0.16, 0.14 },
    {  0.10,  0.14,  0.14,  0.14, 0.17, 0.14 },
    {  0.10,  0.16,  0.16,  0.17, 0.20, 0.17 },
    {  0.10,  0.14,  0.14,  0.14, 0.17, 0.17 }
};

#define MARGIN_FAC 1.1

static gmx_bool is_bond(int nnm, t_nm2type nmt[], char *ai, char *aj, real blen)
{
    int i, j;

    for (i = 0; (i < nnm); i++)
    {
        for (j = 0; (j < nmt[i].nbonds); j++)
        {
            if ((((gmx_strncasecmp(ai, nmt[i].elem, 1) == 0) &&
                  (gmx_strncasecmp(aj, nmt[i].bond[j], 1) == 0)) ||
                 ((gmx_strncasecmp(ai, nmt[i].bond[j], 1) == 0) &&
                  (gmx_strncasecmp(aj, nmt[i].elem, 1) == 0))) &&
                (fabs(blen-nmt[i].blen[j]) <= 0.1*nmt[i].blen[j]))
            {
                return TRUE;
            }
        }
    }
    return FALSE;
}

static void mk_bonds(int nnm, t_nm2type nmt[],
                     t_atoms *atoms, const rvec x[], t_params *bond, int nbond[],
                     gmx_bool bPBC, matrix box)
{
    t_param b;
    int     i, j;
    t_pbc   pbc;
    rvec    dx;
    real    dx2;

    for (i = 0; (i < MAXATOMLIST); i++)
    {
        b.a[i] = -1;
    }
    for (i = 0; (i < MAXFORCEPARAM); i++)
    {
        b.c[i] = 0.0;
    }

    if (bPBC)
    {
        set_pbc(&pbc, -1, box);
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        if ((i % 10) == 0)
        {
            fprintf(stderr, "\ratom %d", i);
            fflush(stderr);
        }
        for (j = i+1; (j < atoms->nr); j++)
        {
            if (bPBC)
            {
                pbc_dx(&pbc, x[i], x[j], dx);
            }
            else
            {
                rvec_sub(x[i], x[j], dx);
            }

            dx2 = iprod(dx, dx);
            if (is_bond(nnm, nmt, *atoms->atomname[i], *atoms->atomname[j],
                        std::sqrt(dx2)))
            {
                b.ai() = i;
                b.aj() = j;
                b.c0() = std::sqrt(dx2);
                add_param_to_list (bond, &b);
                nbond[i]++;
                nbond[j]++;
                if (debug)
                {
                    fprintf(debug, "Bonding atoms %s-%d and %s-%d\n",
                            *atoms->atomname[i], i+1, *atoms->atomname[j], j+1);
                }
            }
        }
    }
    fprintf(stderr, "\ratom %d\n", i);
    fflush(stderr);
}

static int *set_cgnr(t_atoms *atoms, gmx_bool bUsePDBcharge, real *qtot, real *mtot)
{
    int     i, n = 1;
    int    *cgnr;
    double  qt = 0;

    *qtot = *mtot = 0;
    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        if (atoms->pdbinfo && bUsePDBcharge)
        {
            atoms->atom[i].q = atoms->pdbinfo[i].bfac;
        }
        qt     += atoms->atom[i].q;
        *qtot  += atoms->atom[i].q;
        *mtot  += atoms->atom[i].m;
        cgnr[i] = n;
        if (is_int(qt))
        {
            n++;
            qt = 0;
        }
    }
    return cgnr;
}

static gpp_atomtype_t set_atom_type(t_symtab *tab, t_atoms *atoms, t_params *bonds,
                                    int *nbonds, int nnm, t_nm2type nm2t[])
{
    gpp_atomtype_t atype;
    int            nresolved;

    atype = init_atomtype();
    snew(atoms->atomtype, atoms->nr);
    nresolved = nm2type(nnm, nm2t, tab, atoms, atype, nbonds, bonds);
    if (nresolved != atoms->nr)
    {
        gmx_fatal(FARGS, "Could only find a forcefield type for %d out of %d atoms",
                  nresolved, atoms->nr);
    }

    fprintf(stderr, "There are %d different atom types in your sample\n",
            get_atomtype_ntypes(atype));

    return atype;
}

static void lo_set_force_const(t_params *plist, real c[], int nrfp, gmx_bool bRound,
                               gmx_bool bDih, gmx_bool bParam)
{
    int    i, j;
    double cc;
    char   buf[32];

    for (i = 0; (i < plist->nr); i++)
    {
        if (!bParam)
        {
            for (j = 0; j < nrfp; j++)
            {
                c[j] = NOTSET;
            }
        }
        else
        {
            if (bRound)
            {
                sprintf(buf, "%.2e", plist->param[i].c[0]);
                sscanf(buf, "%lf", &cc);
                c[0] = cc;
            }
            else
            {
                c[0] = plist->param[i].c[0];
            }
            if (bDih)
            {
                c[0] *= c[2];
                c[0]  = ((int)(c[0] + 3600)) % 360;
                if (c[0] > 180)
                {
                    c[0] -= 360;
                }
                /* To put the minimum at the current angle rather than the maximum */
                c[0] += 180;
            }
        }
        GMX_ASSERT(nrfp <= MAXFORCEPARAM/2, "Only 6 parameters may be used for an interaction");
        for (j = 0; (j < nrfp); j++)
        {
            plist->param[i].c[j]      = c[j];
            plist->param[i].c[nrfp+j] = c[j];
        }
        set_p_string(&(plist->param[i]), "");
    }
}

static void set_force_const(t_params plist[], real kb, real kt, real kp, gmx_bool bRound,
                            gmx_bool bParam)
{
    real c[MAXFORCEPARAM];

    c[0] = 0;
    c[1] = kb;
    lo_set_force_const(&plist[F_BONDS], c, 2, bRound, FALSE, bParam);
    c[1] = kt;
    lo_set_force_const(&plist[F_ANGLES], c, 2, bRound, FALSE, bParam);
    c[1] = kp;
    c[2] = 3;
    lo_set_force_const(&plist[F_PDIHS], c, 3, bRound, TRUE, bParam);
}

static void calc_angles_dihs(t_params *ang, t_params *dih, const rvec x[], gmx_bool bPBC,
                             matrix box)
{
    int    i, ai, aj, ak, al, t1, t2, t3;
    rvec   r_ij, r_kj, r_kl, m, n;
    real   th, costh, ph;
    t_pbc  pbc;

    if (bPBC)
    {
        set_pbc(&pbc, epbcXYZ, box);
    }
    if (debug)
    {
        pr_rvecs(debug, 0, "X2TOP", box, DIM);
    }
    for (i = 0; (i < ang->nr); i++)
    {
        ai = ang->param[i].ai();
        aj = ang->param[i].aj();
        ak = ang->param[i].ak();
        th = RAD2DEG*bond_angle(x[ai], x[aj], x[ak], bPBC ? &pbc : nullptr,
                                r_ij, r_kj, &costh, &t1, &t2);
        if (debug)
        {
            fprintf(debug, "X2TOP: ai=%3d aj=%3d ak=%3d r_ij=%8.3f r_kj=%8.3f th=%8.3f\n",
                    ai, aj, ak, norm(r_ij), norm(r_kj), th);
        }
        ang->param[i].c0() = th;
    }
    for (i = 0; (i < dih->nr); i++)
    {
        ai = dih->param[i].ai();
        aj = dih->param[i].aj();
        ak = dih->param[i].ak();
        al = dih->param[i].al();
        ph = RAD2DEG*dih_angle(x[ai], x[aj], x[ak], x[al], bPBC ? &pbc : nullptr,
                               r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);
        if (debug)
        {
            fprintf(debug, "X2TOP: ai=%3d aj=%3d ak=%3d al=%3d r_ij=%8.3f r_kj=%8.3f r_kl=%8.3f ph=%8.3f\n",
                    ai, aj, ak, al, norm(r_ij), norm(r_kj), norm(r_kl), ph);
        }
        dih->param[i].c0() = ph;
    }
}

static void dump_hybridization(FILE *fp, t_atoms *atoms, int nbonds[])
{
    int i;

    for (i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "Atom %5s has %1d bonds\n", *atoms->atomname[i], nbonds[i]);
    }
}

static void print_pl(FILE *fp, t_params plist[], int ftp, const char *name,
                     char ***atomname)
{
    int i, j, nral, nrfp;

    if (plist[ftp].nr > 0)
    {
        fprintf(fp, "\n");
        fprintf(fp, "[ %s ]\n", name);
        nral = interaction_function[ftp].nratoms;
        nrfp = interaction_function[ftp].nrfpA;
        for (i = 0; (i < plist[ftp].nr); i++)
        {
            for (j = 0; (j < nral); j++)
            {
                fprintf(fp, "  %5s", *atomname[plist[ftp].param[i].a[j]]);
            }
            for (j = 0; (j < nrfp); j++)
            {
                if (plist[ftp].param[i].c[j] != NOTSET)
                {
                    fprintf(fp, "  %10.3e", plist[ftp].param[i].c[j]);
                }
            }
            fprintf(fp, "\n");
        }
    }
}

static void print_rtp(const char *filenm, const char *title, t_atoms *atoms,
                      t_params plist[], gpp_atomtype_t atype, int cgnr[])
{
    FILE *fp;
    int   i, tp;
    char *tpnm;

    fp = gmx_fio_fopen(filenm, "w");
    fprintf(fp, "; %s\n", title);
    fprintf(fp, "\n");
    fprintf(fp, "[ %s ]\n", *atoms->resinfo[0].name);
    fprintf(fp, "\n");
    fprintf(fp, "[ atoms ]\n");
    for (i = 0; (i < atoms->nr); i++)
    {
        tp = atoms->atom[i].type;
        if ((tpnm = get_atomtype_name(tp, atype)) == nullptr)
        {
            gmx_fatal(FARGS, "tp = %d, i = %d in print_rtp", tp, i);
        }
        fprintf(fp, "%-8s  %12s  %8.4f  %5d\n",
                *atoms->atomname[i], tpnm,
                atoms->atom[i].q, cgnr[i]);
    }
    print_pl(fp, plist, F_BONDS, "bonds", atoms->atomname);
    print_pl(fp, plist, F_ANGLES, "angles", atoms->atomname);
    print_pl(fp, plist, F_PDIHS, "dihedrals", atoms->atomname);
    print_pl(fp, plist, F_IDIHS, "impropers", atoms->atomname);

    gmx_fio_fclose(fp);
}

int gmx_x2top(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] generates a primitive topology from a coordinate file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.",
        "The program can also make an [REF].rtp[ref] entry, which you can then add",
        "to the [REF].rtp[ref] database.[PAR]",
        "When [TT]-param[tt] is set, equilibrium distances and angles",
        "and force constants will be printed in the topology for all",
        "interactions. The equilibrium distances and angles are taken",
        "from the input coordinates, the force constant are set with",
        "command line options.",
        "The force fields somewhat supported currently are:[PAR]",
        "G53a5  GROMOS96 53a5 Forcefield (official distribution)[PAR]",
        "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
        "The corresponding data files can be found in the library directory",
        "with name [TT]atomname2type.n2t[tt]. Check Chapter 5 of the manual for more",
        "information about file formats. By default, the force field selection",
        "is interactive, but you can use the [TT]-ff[tt] option to specify",
        "one of the short names above on the command line instead. In that",
        "case [THISMODULE] just looks for the corresponding file.[PAR]",
    };
    const char        *bugs[] = {
        "The atom type selection is primitive. Virtually no chemical knowledge is used",
        "Periodic boundary conditions screw up the bonding",
        "No improper dihedrals are generated",
        "The atoms to atomtype translation table is incomplete ([TT]atomname2type.n2t[tt] file in the data directory). Please extend it and send the results back to the GROMACS crew."
    };
    FILE              *fp;
    t_params           plist[F_NRE];
    t_excls           *excls;
    gpp_atomtype_t     atype;
    t_nextnb           nnb;
    t_nm2type         *nm2t;
    t_mols             mymol;
    int                nnm;
    char               forcefield[32], ffdir[STRLEN];
    rvec              *x; /* coordinates? */
    int               *nbonds, *cgnr;
    int                bts[] = { 1, 1, 1, 2 };
    matrix             box;    /* box length matrix */
    int                natoms; /* number of atoms in one molecule  */
    int                epbc;
    gmx_bool           bRTP, bTOP, bOPLS;
    t_symtab           symtab;
    real               qtot, mtot;
    char               n2t[STRLEN];
    gmx_output_env_t  *oenv;

    t_filenm           fnm[] = {
        { efSTX, "-f", "conf", ffREAD  },
        { efTOP, "-o", "out",  ffOPTWR },
        { efRTP, "-r", "out",  ffOPTWR }
    };
#define NFILE asize(fnm)
    real               kb                            = 4e5, kt = 400, kp = 5;
    t_restp            rtp_header_settings           = { 0 };
    gmx_bool           bRemoveDihedralIfWithImproper = FALSE;
    gmx_bool           bGenerateHH14Interactions     = TRUE;
    gmx_bool           bKeepAllGeneratedDihedrals    = FALSE;
    int                nrexcl                        = 3;
    gmx_bool           bParam                        = TRUE, bRound = TRUE;
    gmx_bool           bPairs                        = TRUE, bPBC = TRUE;
    gmx_bool           bUsePDBcharge                 = FALSE, bVerbose = FALSE;
    const char        *molnm                         = "ICE";
    const char        *ff                            = "oplsaa";
    t_pargs            pa[]                          = {
        { "-ff",     FALSE, etSTR, {&ff},
          "Force field for your simulation. Type \"select\" for interactive selection." },
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose output in the top file." },
        { "-nexcl", FALSE, etINT,  {&nrexcl},
          "Number of exclusions" },
        { "-H14",    FALSE, etBOOL, {&bGenerateHH14Interactions},
          "Use 3rd neighbour interactions for hydrogen atoms" },
        { "-alldih", FALSE, etBOOL, {&bKeepAllGeneratedDihedrals},
          "Generate all proper dihedrals" },
        { "-remdih", FALSE, etBOOL, {&bRemoveDihedralIfWithImproper},
          "Remove dihedrals on the same bond as an improper" },
        { "-pairs",  FALSE, etBOOL, {&bPairs},
          "Output 1-4 interactions (pairs) in topology file" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-pbc",    FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions." },
        { "-pdbq",  FALSE, etBOOL, {&bUsePDBcharge},
          "Use the B-factor supplied in a [REF].pdb[ref] file for the atomic charges" },
        { "-param", FALSE, etBOOL, {&bParam},
          "Print parameters in the output" },
        { "-round",  FALSE, etBOOL, {&bRound},
          "Round off measured values" },
        { "-kb",    FALSE, etREAL, {&kb},
          "Bonded force constant (kJ/mol/nm^2)" },
        { "-kt",    FALSE, etREAL, {&kt},
          "Angle force constant (kJ/mol/rad^2)" },
        { "-kp",    FALSE, etREAL, {&kp},
          "Dihedral angle force constant (kJ/mol/rad^2)" }
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    bRTP = opt2bSet("-r", NFILE, fnm);
    bTOP = opt2bSet("-o", NFILE, fnm);
    /* C89 requirements mean that these struct members cannot be used in
     * the declaration of pa. So some temporary variables are needed. */
    rtp_header_settings.bRemoveDihedralIfWithImproper = bRemoveDihedralIfWithImproper;
    rtp_header_settings.bGenerateHH14Interactions     = bGenerateHH14Interactions;
    rtp_header_settings.bKeepAllGeneratedDihedrals    = bKeepAllGeneratedDihedrals;
    rtp_header_settings.nrexcl = nrexcl;

    if (!bRTP && !bTOP)
    {
        gmx_fatal(FARGS, "Specify at least one output file");
    }

    /* Force field selection, interactive or direct */
    choose_ff(strcmp(ff, "select") == 0 ? nullptr : ff,
              forcefield, sizeof(forcefield),
              ffdir, sizeof(ffdir));

    bOPLS = (strcmp(forcefield, "oplsaa") == 0);


    mymol.name = gmx_strdup(molnm);
    mymol.nr   = 1;

    /* Init parameter lists */
    init_plist(plist);

    /* Read coordinates */
    t_topology *top;
    snew(top, 1);
    read_tps_conf(opt2fn("-f", NFILE, fnm), top, &epbc, &x, nullptr, box, FALSE);
    t_atoms  *atoms = &top->atoms;
    natoms = atoms->nr;
    if (atoms->pdbinfo == nullptr)
    {
        snew(atoms->pdbinfo, natoms);
    }

    sprintf(n2t, "%s", ffdir);
    nm2t = rd_nm2type(n2t, &nnm);
    if (nnm == 0)
    {
        gmx_fatal(FARGS, "No or incorrect atomname2type.n2t file found (looking for %s)",
                  n2t);
    }
    else
    {
        printf("There are %d name to type translations in file %s\n", nnm, n2t);
    }
    if (debug)
    {
        dump_nm2type(debug, nnm, nm2t);
    }
    printf("Generating bonds from distances...\n");
    snew(nbonds, atoms->nr);
    mk_bonds(nnm, nm2t, atoms, x, &(plist[F_BONDS]), nbonds, bPBC, box);

    open_symtab(&symtab);
    atype = set_atom_type(&symtab, atoms, &(plist[F_BONDS]), nbonds, nnm, nm2t);

    /* Make Angles and Dihedrals */
    snew(excls, atoms->nr);
    printf("Generating angles and dihedrals from bonds...\n");
    init_nnb(&nnb, atoms->nr, 4);
    gen_nnb(&nnb, plist);
    print_nnb(&nnb, "NNB");
    gen_pad(&nnb, atoms, &rtp_header_settings, plist, excls, nullptr, TRUE);
    done_nnb(&nnb);

    if (!bPairs)
    {
        plist[F_LJ14].nr = 0;
    }
    fprintf(stderr,
            "There are %4d %s dihedrals, %4d impropers, %4d angles\n"
            "          %4d pairs,     %4d bonds and  %4d atoms\n",
            plist[F_PDIHS].nr,
            bOPLS ? "Ryckaert-Bellemans" : "proper",
            plist[F_IDIHS].nr, plist[F_ANGLES].nr,
            plist[F_LJ14].nr, plist[F_BONDS].nr, atoms->nr);

    calc_angles_dihs(&plist[F_ANGLES], &plist[F_PDIHS], x, bPBC, box);

    set_force_const(plist, kb, kt, kp, bRound, bParam);

    cgnr = set_cgnr(atoms, bUsePDBcharge, &qtot, &mtot);
    printf("Total charge is %g, total mass is %g\n", qtot, mtot);
    if (bOPLS)
    {
        bts[2] = 3;
        bts[3] = 1;
    }

    if (bTOP)
    {
        fp = ftp2FILE(efTOP, NFILE, fnm, "w");
        print_top_header(fp, ftp2fn(efTOP, NFILE, fnm), TRUE, ffdir, 1.0);

        write_top(fp, nullptr, mymol.name, atoms, FALSE, bts, plist, excls, atype,
                  cgnr, rtp_header_settings.nrexcl);
        print_top_mols(fp, mymol.name, ffdir, nullptr, 0, nullptr, 1, &mymol);

        gmx_ffclose(fp);
    }
    if (bRTP)
    {
        print_rtp(ftp2fn(efRTP, NFILE, fnm), "Generated by x2top",
                  atoms, plist, atype, cgnr);
    }

    if (debug)
    {
        dump_hybridization(debug, atoms, nbonds);
    }
    close_symtab(&symtab);
    sfree(mymol.name);

    printf("\nWARNING: topologies generated by %s can not be trusted at face value.\n",
           output_env_get_program_display_name(oenv));
    printf("         Please verify atomtypes and charges by comparison to other\n");
    printf("         topologies.\n");

    return 0;
}
