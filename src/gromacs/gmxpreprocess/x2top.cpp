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
#include "gmxpre.h"

#include "x2top.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/nm2type.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "hackblock.h"

struct gmx_output_env_t;

static bool is_bond(int nnm, t_nm2type nmt[], char* ai, char* aj, real blen)
{
    int i, j;

    for (i = 0; (i < nnm); i++)
    {
        for (j = 0; (j < nmt[i].nbonds); j++)
        {
            if ((((gmx::equalCaseInsensitive(ai, nmt[i].elem, 1))
                  && (gmx::equalCaseInsensitive(aj, nmt[i].bond[j], 1)))
                 || ((gmx::equalCaseInsensitive(ai, nmt[i].bond[j], 1))
                     && (gmx::equalCaseInsensitive(aj, nmt[i].elem, 1))))
                && (std::fabs(blen - nmt[i].blen[j]) <= 0.1 * nmt[i].blen[j]))
            {
                return TRUE;
            }
        }
    }
    return FALSE;
}

static void mk_bonds(int                 nnm,
                     t_nm2type           nmt[],
                     t_atoms*            atoms,
                     const rvec          x[],
                     InteractionsOfType* bond,
                     int                 nbond[],
                     bool                bPBC,
                     matrix              box)
{
    int   i, j;
    t_pbc pbc;
    rvec  dx;
    real  dx2;

    std::array<real, MAXFORCEPARAM> forceParam = { 0.0 };
    if (bPBC)
    {
        set_pbc(&pbc, PbcType::Unset, box);
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        for (j = i + 1; (j < atoms->nr); j++)
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
            if (is_bond(nnm, nmt, *atoms->atomname[i], *atoms->atomname[j], std::sqrt(dx2)))
            {
                forceParam[0]             = std::sqrt(dx2);
                std::vector<int> atomPair = { i, j };
                add_param_to_list(bond, InteractionOfType(atomPair, forceParam));
                nbond[i]++;
                nbond[j]++;
            }
        }
    }
}

static int* set_cgnr(t_atoms* atoms, bool bUsePDBcharge, real* qtot, real* mtot)
{
    int    i, n = 1;
    int*   cgnr;
    double qt = 0;

    *qtot = *mtot = 0;
    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        if (atoms->pdbinfo && bUsePDBcharge)
        {
            atoms->atom[i].q = atoms->pdbinfo[i].bfac;
        }
        qt += atoms->atom[i].q;
        *qtot += atoms->atom[i].q;
        *mtot += atoms->atom[i].m;
        cgnr[i] = n;
        if (is_int(qt))
        {
            n++;
            qt = 0;
        }
    }
    return cgnr;
}

static void set_atom_type(PreprocessingAtomTypes* atypes,
                          t_atoms*                atoms,
                          InteractionsOfType*     bonds,
                          int*                    nbonds,
                          int                     nnm,
                          t_nm2type               nm2t[],
                          const gmx::MDLogger&    logger)
{
    int nresolved;

    snew(atoms->atomtype, atoms->nr);
    nresolved = nm2type(nnm, nm2t, atoms, atypes, nbonds, bonds);
    if (nresolved != atoms->nr)
    {
        gmx_fatal(FARGS, "Could only find a forcefield type for %d out of %d atoms", nresolved, atoms->nr);
    }

    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted("There are %zu different atom types in your sample", atypes->size());
}

static void lo_set_force_const(InteractionsOfType* plist, real c[], int nrfp, bool bRound, bool bDih, bool bParam)
{
    double cc;
    char   buf[32];

    for (auto& param : plist->interactionTypes)
    {
        if (!bParam)
        {
            for (int j = 0; j < nrfp; j++)
            {
                c[j] = NOTSET;
            }
        }
        else
        {
            if (bRound)
            {
                sprintf(buf, "%.2e", param.c0());
                sscanf(buf, "%lf", &cc);
                c[0] = cc;
            }
            else
            {
                c[0] = param.c0();
            }
            if (bDih)
            {
                c[0] *= c[2];
                c[0] = (static_cast<int>(c[0] + 3600)) % 360;
                if (c[0] > 180)
                {
                    c[0] -= 360;
                }
                /* To put the minimum at the current angle rather than the maximum */
                c[0] += 180;
            }
        }
        GMX_ASSERT(nrfp <= MAXFORCEPARAM / 2, "Only 6 parameters may be used for an interaction");
        std::array<real, MAXFORCEPARAM> forceParam;
        for (int j = 0; (j < nrfp); j++)
        {
            forceParam[j]        = c[j];
            forceParam[nrfp + j] = c[j];
        }
        param = InteractionOfType(param.atoms(), forceParam);
    }
}

static void set_force_const(gmx::ArrayRef<InteractionsOfType> plist, real kb, real kt, real kp, bool bRound, bool bParam)
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

static void calc_angles_dihs(InteractionsOfType* ang, InteractionsOfType* dih, const rvec x[], bool bPBC, matrix box)
{
    int   t1, t2, t3;
    rvec  r_ij, r_kj, r_kl, m, n;
    real  costh;
    t_pbc pbc;

    if (bPBC)
    {
        set_pbc(&pbc, PbcType::Xyz, box);
    }
    for (auto& angle : ang->interactionTypes)
    {
        int  ai = angle.ai();
        int  aj = angle.aj();
        int  ak = angle.ak();
        real th = gmx::c_rad2Deg
                  * bond_angle(x[ai], x[aj], x[ak], bPBC ? &pbc : nullptr, r_ij, r_kj, &costh, &t1, &t2);
        angle.setForceParameter(0, th);
    }
    for (auto dihedral : dih->interactionTypes)
    {
        int  ai = dihedral.ai();
        int  aj = dihedral.aj();
        int  ak = dihedral.ak();
        int  al = dihedral.al();
        real ph =
                gmx::c_rad2Deg
                * dih_angle(
                        x[ai], x[aj], x[ak], x[al], bPBC ? &pbc : nullptr, r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);
        dihedral.setForceParameter(0, ph);
    }
}

static void dump_hybridization(FILE* fp, t_atoms* atoms, int nbonds[])
{
    for (int i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "Atom %5s has %1d bonds\n", *atoms->atomname[i], nbonds[i]);
    }
}

static void
print_pl(FILE* fp, gmx::ArrayRef<const InteractionsOfType> plist, int ftp, const char* name, char*** atomname)
{
    if (!plist[ftp].interactionTypes.empty())
    {
        fprintf(fp, "\n");
        fprintf(fp, "[ %s ]\n", name);
        int nrfp = interaction_function[ftp].nrfpA;
        for (const auto& param : plist[ftp].interactionTypes)
        {
            gmx::ArrayRef<const int>  atoms      = param.atoms();
            gmx::ArrayRef<const real> forceParam = param.forceParam();
            for (const auto& atom : atoms)
            {
                fprintf(fp, "  %5s", *atomname[atom]);
            }
            for (int j = 0; (j < nrfp); j++)
            {
                if (forceParam[j] != NOTSET)
                {
                    fprintf(fp, "  %10.3e", forceParam[j]);
                }
            }
            fprintf(fp, "\n");
        }
    }
}

static void print_rtp(const char*                             filenm,
                      const char*                             title,
                      t_atoms*                                atoms,
                      gmx::ArrayRef<const InteractionsOfType> plist,
                      PreprocessingAtomTypes*                 atypes,
                      int                                     cgnr[])
{
    FILE* fp;
    int   i, tp;

    fp = gmx_fio_fopen(filenm, "w");
    fprintf(fp, "; %s\n", title);
    fprintf(fp, "\n");
    fprintf(fp, "[ %s ]\n", *atoms->resinfo[0].name);
    fprintf(fp, "\n");
    fprintf(fp, "[ atoms ]\n");
    for (i = 0; (i < atoms->nr); i++)
    {
        tp        = atoms->atom[i].type;
        auto tpnm = atypes->atomNameFromAtomType(tp);
        if (!tpnm.has_value())
        {
            gmx_fatal(FARGS, "tp = %d, i = %d in print_rtp", tp, i);
        }
        fprintf(fp, "%-8s  %12s  %8.4f  %5d\n", *atoms->atomname[i], tpnm->c_str(), atoms->atom[i].q, cgnr[i]);
    }
    print_pl(fp, plist, F_BONDS, "bonds", atoms->atomname);
    print_pl(fp, plist, F_ANGLES, "angles", atoms->atomname);
    print_pl(fp, plist, F_PDIHS, "dihedrals", atoms->atomname);
    print_pl(fp, plist, F_IDIHS, "impropers", atoms->atomname);

    gmx_fio_fclose(fp);
}

int gmx_x2top(int argc, char* argv[])
{
    const char* desc[] = {
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
    const char* bugs[] = {
        "The atom type selection is primitive. Virtually no chemical knowledge is used",
        "Periodic boundary conditions screw up the bonding",
        "No improper dihedrals are generated",
        ("The atoms to atomtype translation table is incomplete ([TT]atomname2type.n2t[tt] file in "
         "the data directory). Please extend it and send the results back to the GROMACS crew.")
    };
    FILE*                                 fp;
    std::array<InteractionsOfType, F_NRE> plist;
    t_excls*                              excls;
    t_nm2type*                            nm2t;
    t_mols                                mymol;
    int                                   nnm;
    char                                  forcefield[32];
    rvec*                                 x; /* coordinates? */
    int *                                 nbonds, *cgnr;
    int                                   bts[] = { 1, 1, 1, 2 };
    matrix                                box;    /* box length matrix */
    int                                   natoms; /* number of atoms in one molecule  */
    PbcType                               pbcType;
    bool                                  bRTP, bTOP, bOPLS;
    real                                  qtot, mtot;
    gmx_output_env_t*                     oenv;

    t_filenm fnm[] = { { efSTX, "-f", "conf", ffREAD },
                       { efTOP, "-o", "out", ffOPTWR },
                       { efRTP, "-r", "out", ffOPTWR } };
#define NFILE asize(fnm)
    real              kb = 4e5, kt = 400, kp = 5;
    PreprocessResidue rtp_header_settings;
    bool              bRemoveDihedralIfWithImproper = FALSE;
    bool              bGenerateHH14Interactions     = TRUE;
    bool              bKeepAllGeneratedDihedrals    = FALSE;
    int               nrexcl                        = 3;
    bool              bParam = TRUE, bRound = TRUE;
    bool              bPairs = TRUE, bPBC = TRUE;
    bool              bUsePDBcharge = FALSE, bVerbose = FALSE;
    const char*       molnm = "ICE";
    const char*       ff    = "oplsaa";
    t_pargs           pa[]  = {
        { "-ff",
          FALSE,
          etSTR,
          { &ff },
          "Force field for your simulation. Type \"select\" for interactive selection." },
        { "-v", FALSE, etBOOL, { &bVerbose }, "Generate verbose output in the top file." },
        { "-nexcl", FALSE, etINT, { &nrexcl }, "Number of exclusions" },
        { "-H14",
          FALSE,
          etBOOL,
          { &bGenerateHH14Interactions },
          "Use 3rd neighbour interactions for hydrogen atoms" },
        { "-alldih",
          FALSE,
          etBOOL,
          { &bKeepAllGeneratedDihedrals },
          "Generate all proper dihedrals" },
        { "-remdih",
          FALSE,
          etBOOL,
          { &bRemoveDihedralIfWithImproper },
          "Remove dihedrals on the same bond as an improper" },
        { "-pairs", FALSE, etBOOL, { &bPairs }, "Output 1-4 interactions (pairs) in topology file" },
        { "-name", FALSE, etSTR, { &molnm }, "Name of your molecule" },
        { "-pbc", FALSE, etBOOL, { &bPBC }, "Use periodic boundary conditions." },
        { "-pdbq",
          FALSE,
          etBOOL,
          { &bUsePDBcharge },
          "Use the B-factor supplied in a [REF].pdb[ref] file for the atomic charges" },
        { "-param", FALSE, etBOOL, { &bParam }, "Print parameters in the output" },
        { "-round", FALSE, etBOOL, { &bRound }, "Round off measured values" },
        { "-kb", FALSE, etREAL, { &kb }, "Bonded force constant (kJ/mol/nm^2)" },
        { "-kt", FALSE, etREAL, { &kt }, "Angle force constant (kJ/mol/rad^2)" },
        { "-kp", FALSE, etREAL, { &kp }, "Dihedral angle force constant (kJ/mol/rad^2)" }
    };

    if (!parse_common_args(
                &argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
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
    rtp_header_settings.nrexcl                        = nrexcl;

    if (!bRTP && !bTOP)
    {
        gmx_fatal(FARGS, "Specify at least one output file");
    }

    gmx::LoggerBuilder builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Info, &gmx::TextOutputFile::standardOutput());
    builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &gmx::TextOutputFile::standardError());
    gmx::LoggerOwner logOwner(builder.build());
    gmx::MDLogger    logger(logOwner.logger());


    /* Force field selection, interactive or direct */
    auto ffdir = choose_ff(strcmp(ff, "select") == 0 ? nullptr : ff, forcefield, sizeof(forcefield), logger);

    bOPLS = (strcmp(forcefield, "oplsaa") == 0);


    mymol.name = gmx_strdup(molnm);
    mymol.nr   = 1;

    /* Read coordinates */
    t_topology* top;
    snew(top, 1);
    read_tps_conf(opt2fn("-f", NFILE, fnm), top, &pbcType, &x, nullptr, box, FALSE);
    t_atoms* atoms = &top->atoms;
    natoms         = atoms->nr;
    if (atoms->pdbinfo == nullptr)
    {
        snew(atoms->pdbinfo, natoms);
    }

    nm2t = rd_nm2type(ffdir, &nnm);
    if (nnm == 0)
    {
        gmx_fatal(FARGS,
                  "No or incorrect atomname2type.n2t file found (looking for %s)",
                  ffdir.string().c_str());
    }
    else
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "There are %d name to type translations in file %s", nnm, ffdir.string().c_str());
    }
    if (debug)
    {
        dump_nm2type(debug, nnm, nm2t);
    }
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("Generating bonds from distances...");
    snew(nbonds, atoms->nr);
    mk_bonds(nnm, nm2t, atoms, x, &(plist[F_BONDS]), nbonds, bPBC, box);

    PreprocessingAtomTypes atypes;
    set_atom_type(&atypes, atoms, &(plist[F_BONDS]), nbonds, nnm, nm2t, logger);

    /* Make Angles and Dihedrals */
    snew(excls, atoms->nr);
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted("Generating angles and dihedrals from bonds...");
    gen_pad(atoms, gmx::arrayRefFromArray(&rtp_header_settings, 1), plist, excls, {}, TRUE, {});

    if (!bPairs)
    {
        plist[F_LJ14].interactionTypes.clear();
    }
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "There are %4zu %s dihedrals, %4zu impropers, %4zu angles\n"
                    "          %4zu pairs,     %4zu bonds and  %4d atoms\n",
                    plist[F_PDIHS].size(),
                    bOPLS ? "Ryckaert-Bellemans" : "proper",
                    plist[F_IDIHS].size(),
                    plist[F_ANGLES].size(),
                    plist[F_LJ14].size(),
                    plist[F_BONDS].size(),
                    atoms->nr);

    calc_angles_dihs(&plist[F_ANGLES], &plist[F_PDIHS], x, bPBC, box);

    set_force_const(plist, kb, kt, kp, bRound, bParam);

    cgnr = set_cgnr(atoms, bUsePDBcharge, &qtot, &mtot);
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("Total charge is %g, total mass is %g", qtot, mtot);
    if (bOPLS)
    {
        bts[2] = 3;
        bts[3] = 1;
    }

    if (bTOP)
    {
        fp = ftp2FILE(efTOP, NFILE, fnm, "w");
        print_top_header(fp, ftp2fn(efTOP, NFILE, fnm), TRUE, ffdir, 1.0);

        write_top(fp,
                  {},
                  mymol.name.c_str(),
                  atoms,
                  FALSE,
                  bts,
                  plist,
                  excls,
                  &atypes,
                  cgnr,
                  rtp_header_settings.nrexcl);
        print_top_mols(fp, mymol.name.c_str(), ffdir, nullptr, {}, gmx::arrayRefFromArray(&mymol, 1));

        gmx_ffclose(fp);
    }
    if (bRTP)
    {
        print_rtp(ftp2fn(efRTP, NFILE, fnm), "Generated by x2top", atoms, plist, &atypes, cgnr);
    }

    if (debug)
    {
        dump_hybridization(debug, atoms, nbonds);
    }

    GMX_LOG(logger.warning)
            .asParagraph()
            .appendTextFormatted(
                    "Topologies generated by %s can not be trusted at face value. "
                    "Please verify atomtypes and charges by comparison to other topologies.",
                    output_env_get_program_display_name(oenv));

    return 0;
}
