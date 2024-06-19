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

#include "topio.h"

#include <cassert>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>

#include <sys/types.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/gmxcpp.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/topshake.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/vsite_parm.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct t_nbparam;

#define OPENDIR '['  /* starting sign for directive */
#define CLOSEDIR ']' /* ending sign for directive   */

static void gen_pairs(const InteractionsOfType& nbs, InteractionsOfType* pairs, real fudge, CombinationRule comb)
{
    real scaling;
    int  ntp = nbs.size();
    int  nnn = static_cast<int>(std::sqrt(static_cast<double>(ntp)));
    GMX_ASSERT(nnn * nnn == ntp,
               "Number of pairs of generated non-bonded parameters should be a perfect square");
    int nrfp  = NRFP(F_LJ);
    int nrfpA = interaction_function[F_LJ14].nrfpA;
    int nrfpB = interaction_function[F_LJ14].nrfpB;

    if ((nrfp != nrfpA) || (nrfpA != nrfpB))
    {
        gmx_incons("Number of force parameters in gen_pairs wrong");
    }

    fprintf(stderr, "Generating 1-4 interactions: fudge = %g\n", fudge);
    pairs->interactionTypes.clear();
    int                             i = 0;
    std::array<int, 2>              atomNumbers;
    std::array<real, MAXFORCEPARAM> forceParam = { NOTSET };
    for (const auto& type : nbs.interactionTypes)
    {
        /* Copy type.atoms */
        atomNumbers = { i / nnn, i % nnn };
        /* Copy normal and FEP parameters and multiply by fudge factor */
        gmx::ArrayRef<const real> existingParam = type.forceParam();
        GMX_RELEASE_ASSERT(2 * nrfp <= MAXFORCEPARAM,
                           "Can't have more parameters than half of maximum parameter number");
        for (int j = 0; j < nrfp; j++)
        {
            /* If we are using sigma/epsilon values, only the epsilon values
             * should be scaled, but not sigma.
             * The sigma values have even indices 0,2, etc.
             */
            if ((comb == CombinationRule::Arithmetic || comb == CombinationRule::GeomSigEps)
                && (j % 2 == 0))
            {
                scaling = 1.0;
            }
            else
            {
                scaling = fudge;
            }

            forceParam[j]        = scaling * existingParam[j];
            forceParam[nrfp + j] = scaling * existingParam[j];
        }
        pairs->interactionTypes.emplace_back(atomNumbers, forceParam);
        i++;
    }
}

double check_mol(const gmx_mtop_t* mtop, WarningHandler* wi)
{
    char   buf[256];
    int    i, ri;
    double q;
    real   m, mB;

    /* Check mass and charge */
    q = 0.0;

    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const t_atoms* atoms = &mtop->moltype[molb.type].atoms;
        for (i = 0; (i < atoms->nr); i++)
        {
            q += molb.nmol * atoms->atom[i].q;
            m               = atoms->atom[i].m;
            mB              = atoms->atom[i].mB;
            ParticleType pt = atoms->atom[i].ptype;
            /* If the particle is an atom or a nucleus it must have a mass,
             * else, if it is a shell, a vsite or a bondshell it can have mass zero
             */
            if (((m <= 0.0) || (mB <= 0.0)) && ((pt == ParticleType::Atom) || (pt == ParticleType::Nucleus)))
            {
                ri = atoms->atom[i].resind;
                sprintf(buf,
                        "atom %s (Res %s-%d) has mass %g (state A) / %g (state B)\n",
                        *(atoms->atomname[i]),
                        *(atoms->resinfo[ri].name),
                        atoms->resinfo[ri].nr,
                        m,
                        mB);
                wi->addError(buf);
            }
            else if (((m != 0) || (mB != 0)) && (pt == ParticleType::VSite))
            {
                ri = atoms->atom[i].resind;
                sprintf(buf,
                        "virtual site %s (Res %s-%d) has non-zero mass %g (state A) / %g (state "
                        "B)\n"
                        "     Check your topology.\n",
                        *(atoms->atomname[i]),
                        *(atoms->resinfo[ri].name),
                        atoms->resinfo[ri].nr,
                        m,
                        mB);
                wi->addError(buf);
                /* The following statements make LINCS break! */
                /* atoms->atom[i].m=0; */
            }
        }
    }
    return q;
}

/*! \brief Describe molecule and involved atoms for a dihedral interaction of a given type
 *
 * Searches in the dihedrals of type 3 (Ryckaert-Bellemans or Fourier) for an interaction
 * type matching the interactionType input parameters, returning for the first match
 * a string with the name of the  molecule and the 4 involved atoms.
 *
 * Precondition: the interaction should exist in the topology
 *
 */
static std::string describeAtomsForRBDihedralOfGivenType(const gmx_mtop_t& mtop, int interactionType)
{
    for (const auto& molt : mtop.moltype)
    {
        const int* ia = molt.ilist[F_RBDIHS].iatoms.data();
        for (int i = 0; (i < molt.ilist[F_RBDIHS].size());)
        {
            const int type  = ia[0];
            const int ftype = mtop.ffparams.functype[type];
            const int nra   = interaction_function[ftype].nratoms;
            if (type == interactionType)
            {
                return gmx::formatString(
                        "First such dihedral in molecule %s, involving atoms %d %d %d %d",
                        *(molt.name),
                        ia[1],
                        ia[2],
                        ia[3],
                        ia[4]);
            }
            ia += nra + 1;
            i += nra + 1;
        }
    }
    gmx_fatal(FARGS, "Precondition violation: could not find RB interaction of given type %d", interactionType);
}

void checkRBDihedralSum(const gmx_mtop_t& mtop, const t_inputrec& ir, WarningHandler* wi)
{
    /*
     * The sum of the RB dihedral coefficient being zero is relevant when:
     *     - Free energy computation is being performed (dHdl)
     *     - Comparing energies between force field ports and/or other MD codes
     *  because this affect the value of the potential energy.
     *
     *  We can clearly detect problems in the first situation: free energy is enabled, stateA
     *  and stateB have a different sum (=> potential energy "offset" at 0 degree), so we emit a
     *  warning. The second case is more subtle, since formally the potential is up to a constant,
     *  which does not affect the dynamics. We therefore only emit a note.
     */

    // Mistakes here are typically due to using the wrong formula to port dihedrals, not numerical
    // issues, so we use a relatively large tolerance
    const real         absoluteTolerance        = 0.01;
    int                numSumCoefficientNotZero = 0;
    std::optional<int> indexOfFirstSumCoefficientNotZero;
    int                numSumCoefficientDifferentInStateAStateB = 0;
    std::optional<int> indexOfSumCoefficientDifferentInStateAStateB;

    for (int j = 0; j < gmx::ssize(mtop.ffparams.functype); ++j)
    {
        int ftype = mtop.ffparams.functype[j];
        if (ftype == F_RBDIHS)
        {
            const t_iparams& params = mtop.ffparams.iparams[j];
            const real       sum_a =
                    std::accumulate(std::begin(params.rbdihs.rbcA), std::end(params.rbdihs.rbcA), 0.0);
            const real sum_b =
                    std::accumulate(std::begin(params.rbdihs.rbcB), std::end(params.rbdihs.rbcB), 0.0);

            if (std::abs(sum_a - sum_b) > absoluteTolerance)
            {
                numSumCoefficientDifferentInStateAStateB++;
                if (!indexOfSumCoefficientDifferentInStateAStateB.has_value())
                {
                    indexOfSumCoefficientDifferentInStateAStateB = j;
                }
            }

            if (std::abs(sum_a) > absoluteTolerance || std::abs(sum_b) > absoluteTolerance)
            {
                numSumCoefficientNotZero++;
                if (!indexOfFirstSumCoefficientNotZero.has_value())
                {
                    indexOfFirstSumCoefficientNotZero = j;
                }
            }
        }
    }

    // At this stage of grompp, we only have the interactions that are going to be used in the
    // simulation - this eliminates warning unrelated to the user's system, but also unfortunately
    // means that we cannot easily identify the file/line where the offending parameters were
    // originally defined. We can however go through the molecule types and print one instance of
    // the offending dihedral.

    auto generateMessage = [](int numDihedrals, const std::string& note, const std::string& involvedAtoms) {
        return gmx::formatString(
                "%d dihedrals with function type 3 (Ryckaert-Bellemans or Fourier) have "
                "coefficients %s"
                "\n%s",
                numDihedrals,
                note.c_str(),
                involvedAtoms.c_str());
    };

    if (numSumCoefficientNotZero > 0)
    {
        std::string involvedAtoms =
                describeAtomsForRBDihedralOfGivenType(mtop, indexOfFirstSumCoefficientNotZero.value());
        std::string note =
                "that do not sum to zero. This does not affect the simulation and can "
                "be ignored, unless you are comparing potential energy values with other force "
                "field ports and/or MD software.";
        std::string message = generateMessage(numSumCoefficientNotZero, note, involvedAtoms);
        wi->addNote(message);
    }

    if (numSumCoefficientDifferentInStateAStateB > 0 && ir.efep != FreeEnergyPerturbationType::No)
    {
        std::string involvedAtoms = describeAtomsForRBDihedralOfGivenType(
                mtop, indexOfSumCoefficientDifferentInStateAStateB.value());
        std::string note =
                "whose sums do not match in state A and B. This could introduce an "
                "undesired offset in dHdl values.";
        std::string message =
                generateMessage(numSumCoefficientDifferentInStateAStateB, note, involvedAtoms);
        wi->addWarning(message);
    }
}

/*! \brief Returns the rounded charge of a molecule, when close to integer, otherwise returns the original charge.
 *
 * The results of this routine are only used for checking and for
 * printing warning messages. Thus we can assume that charges of molecules
 * should be integer. If the user wanted non-integer molecular charge,
 * an undesired warning is printed and the user should use grompp -maxwarn 1.
 *
 * \param qMol     The total, unrounded, charge of the molecule
 * \param sumAbsQ  The sum of absolute values of the charges, used for determining the tolerance for the rounding.
 */
static double roundedMoleculeCharge(double qMol, double sumAbsQ)
{
    /* We use a tolerance of 1e-6 for inaccuracies beyond the 6th decimal
     * of the charges for ascii float truncation in the topology files.
     * Although the summation here uses double precision, the charges
     * are read and stored in single precision when real=float. This can
     * lead to rounding errors of half the least significant bit.
     * Note that, unfortunately, we can not assume addition of random
     * rounding errors. It is not entirely unlikely that many charges
     * have a near half-bit rounding error with the same sign.
     */
    double tolAbs = 1e-6;
    double tol    = std::max(tolAbs, 0.5 * GMX_REAL_EPS * sumAbsQ);
    double qRound = std::round(qMol);
    if (std::abs(qMol - qRound) <= tol)
    {
        return qRound;
    }
    else
    {
        return qMol;
    }
}

static void sum_q(const t_atoms* atoms, int numMols, double* qTotA, double* qTotB)
{
    /* sum charge */
    double qmolA    = 0;
    double qmolB    = 0;
    double sumAbsQA = 0;
    double sumAbsQB = 0;
    for (int i = 0; i < atoms->nr; i++)
    {
        qmolA += atoms->atom[i].q;
        qmolB += atoms->atom[i].qB;
        sumAbsQA += std::abs(atoms->atom[i].q);
        sumAbsQB += std::abs(atoms->atom[i].qB);
    }

    *qTotA += numMols * roundedMoleculeCharge(qmolA, sumAbsQA);
    *qTotB += numMols * roundedMoleculeCharge(qmolB, sumAbsQB);
}

static void get_nbparm(char* nb_str, char* comb_str, VanDerWaalsPotential* nb, CombinationRule* comb, WarningHandler* wi)
{
    *nb = VanDerWaalsPotential::Count;
    for (auto i : gmx::EnumerationArray<VanDerWaalsPotential, bool>::keys())
    {
        if (gmx_strcasecmp(nb_str, enumValueToString(i)) == 0)
        {
            *nb = i;
        }
    }
    if (*nb == VanDerWaalsPotential::Count)
    {
        int integerValue = strtol(nb_str, nullptr, 10);
        if ((integerValue < 1) || (integerValue >= static_cast<int>(VanDerWaalsPotential::Count)))
        {
            std::string message =
                    gmx::formatString("Invalid nonbond function selector '%s' using %s",
                                      nb_str,
                                      enumValueToString(VanDerWaalsPotential::LJ));
            wi->addError(message);
            *nb = VanDerWaalsPotential::LJ;
        }
        else
        {
            *nb = static_cast<VanDerWaalsPotential>(integerValue);
        }
    }
    *comb = CombinationRule::Count;
    for (auto i : gmx::EnumerationArray<CombinationRule, bool>::keys())
    {
        if (gmx_strcasecmp(comb_str, enumValueToString(i)) == 0)
        {
            *comb = i;
        }
    }
    if (*comb == CombinationRule::Count)
    {
        int integerValue = strtol(comb_str, nullptr, 10);
        if ((integerValue < 1) || (integerValue >= static_cast<int>(CombinationRule::Count)))
        {
            std::string message =
                    gmx::formatString("Invalid combination rule selector '%s' using %s",
                                      comb_str,
                                      enumValueToString(CombinationRule::Geometric));
            wi->addError(message);
            *comb = CombinationRule::Geometric;
        }
        else
        {
            *comb = static_cast<CombinationRule>(integerValue);
        }
    }
}

/*! \brief Parses define and include flags.
 *
 * Returns a vector of parsed include/define flags, with an extra nullptr entry at the back
 * for consumers that expect null-terminated char** structures.
 */
static std::vector<char*> cpp_opts(const char* define, const char* include, WarningHandler* wi)
{
    int         n, len;
    const char* cppadds[2];
    const char* option[2] = { "-D", "-I" };
    const char* nopt[2]   = { "define", "include" };
    const char* ptr;
    const char* rptr;
    char*       buf;
    char        warn_buf[STRLEN];

    cppadds[0] = define;
    cppadds[1] = include;
    std::vector<char*> cppOptions;
    for (n = 0; (n < 2); n++)
    {
        if (cppadds[n])
        {
            ptr = cppadds[n];
            while (*ptr != '\0')
            {
                while ((*ptr != '\0') && isspace(*ptr))
                {
                    ptr++;
                }
                rptr = ptr;
                while ((*rptr != '\0') && !isspace(*rptr))
                {
                    rptr++;
                }
                len = (rptr - ptr);
                if (len > 2)
                {
                    snew(buf, (len + 1));
                    strncpy(buf, ptr, len);
                    if (strstr(ptr, option[n]) != ptr)
                    {
                        wi->setFileAndLineNumber("mdp file", -1);
                        sprintf(warn_buf, "Malformed %s option %s", nopt[n], buf);
                        wi->addWarning(warn_buf);
                    }
                    else
                    {
                        cppOptions.emplace_back(gmx_strdup(buf));
                    }
                    sfree(buf);
                    ptr = rptr;
                }
            }
        }
    }
    // Users of cppOptions expect a null last element.
    cppOptions.emplace_back(nullptr);
    return cppOptions;
}


static void make_atoms_sys(gmx::ArrayRef<const gmx_molblock_t>      molblock,
                           gmx::ArrayRef<const MoleculeInformation> molinfo,
                           t_atoms*                                 atoms)
{
    atoms->nr   = 0;
    atoms->atom = nullptr;

    for (const gmx_molblock_t& molb : molblock)
    {
        const t_atoms& mol_atoms = molinfo[molb.type].atoms;

        srenew(atoms->atom, atoms->nr + molb.nmol * mol_atoms.nr);

        for (int m = 0; m < molb.nmol; m++)
        {
            for (int a = 0; a < mol_atoms.nr; a++)
            {
                atoms->atom[atoms->nr++] = mol_atoms.atom[a];
            }
        }
    }
}


static char** read_topol(const char*                                 infile,
                         const std::optional<std::filesystem::path>& outfile,
                         const char*                                 define,
                         const char*                                 include,
                         t_symtab*                                   symtab,
                         PreprocessingAtomTypes*                     atypes,
                         std::vector<MoleculeInformation>*           molinfo,
                         std::unique_ptr<MoleculeInformation>*       intermolecular_interactions,
                         gmx::ArrayRef<InteractionsOfType>           interactions,
                         CombinationRule*                            combination_rule,
                         double*                                     reppow,
                         t_gromppopts*                               opts,
                         real*                                       fudgeQQ,
                         std::vector<gmx_molblock_t>*                molblock,
                         bool*                ffParametrizedWithHBondConstraints,
                         bool                 bFEP,
                         bool                 bZero,
                         bool                 usingFullRangeElectrostatics,
                         WarningHandler*      wi,
                         const gmx::MDLogger& logger)
{
    FILE*                out;
    int                  sl;
    char *               pline = nullptr, **title = nullptr;
    char                 line[STRLEN], errbuf[256], comb_str[256], nb_str[256];
    char                 genpairs[32];
    char *               dirstr, *dummy2;
    int                  nrcopies, nscan, ncombs, ncopy;
    double               fLJ, fQQ, fPOW;
    MoleculeInformation* mi0 = nullptr;
    DirStack*            DS;
    Directive            d, newd;
    t_nbparam **         nbparam, **pair;
    real                 fudgeLJ = -1; /* Multiplication factor to generate 1-4 from LJ */
    bool                 bReadDefaults, bReadMolType, bGenPairs, bWarn_copy_A_B;
    double               qt = 0, qBt = 0; /* total charge */
    int                  dcatt = -1, nmol_couple;
    /* File handling variables */
    int         status;
    bool        done;
    gmx_cpp_t   handle;
    char*       tmp_line = nullptr;
    char        warn_buf[STRLEN];
    const char* floating_point_arithmetic_tip =
            "Total charge should normally be an integer. See\n"
            "https://manual.gromacs.org/current/user-guide/floating-point.html\n"
            "for discussion on how close it should be to an integer.\n";
    /* We need to open the output file before opening the input file,
     * because cpp_open_file can change the current working directory.
     */
    if (outfile)
    {
        out = gmx_fio_fopen(outfile.value(), "w");
    }
    else
    {
        out = nullptr;
    }

    /* open input file */
    auto cpp_opts_return = cpp_opts(define, include, wi);
    status               = cpp_open_file(infile, &handle, cpp_opts_return.data());
    if (status != 0)
    {
        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
    }

    /* some local variables */
    DS_Init(&DS);                   /* directive stack	*/
    d       = Directive::d_invalid; /* first thing should be a directive */
    nbparam = nullptr;              /* The temporary non-bonded matrix */
    pair    = nullptr;              /* The temporary pair interaction matrix */
    std::vector<std::vector<gmx::ExclusionBlock>> exclusionBlocks;
    VanDerWaalsPotential                          nb_funct = VanDerWaalsPotential::LJ;

    *reppow = 12.0; /* Default value for repulsion power     */

    /* Init the number of CMAP torsion angles  and grid spacing */
    interactions[F_CMAP].cmapGridSpacing_ = 0;
    interactions[F_CMAP].numCmaps_        = 0;

    bWarn_copy_A_B = bFEP;

    PreprocessingBondAtomType bondAtomType;
    /* parse the actual file */
    bReadDefaults = FALSE;
    bGenPairs     = FALSE;
    bReadMolType  = FALSE;
    nmol_couple   = 0;

    do
    {
        status = cpp_read_line(&handle, STRLEN, line);
        done   = (status == eCPP_EOF);
        if (!done)
        {
            if (status != eCPP_OK)
            {
                gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
            }
            else if (out)
            {
                fprintf(out, "%s\n", line);
            }

            wi->setFileAndLineNumber(cpp_cur_file(&handle), cpp_cur_linenr(&handle));

            pline = gmx_strdup(line);

            /* Strip trailing '\' from pline, if it exists */
            sl = strlen(pline);
            if ((sl > 0) && (pline[sl - 1] == CONTINUE))
            {
                pline[sl - 1] = ' ';
            }

            /* build one long line from several fragments - necessary for CMAP */
            while (continuing(line))
            {
                status = cpp_read_line(&handle, STRLEN, line);
                wi->setFileAndLineNumber(cpp_cur_file(&handle), cpp_cur_linenr(&handle));

                /* Since we depend on the '\' being present to continue to read, we copy line
                 * to a tmp string, strip the '\' from that string, and cat it to pline
                 */
                tmp_line = gmx_strdup(line);

                sl = strlen(tmp_line);
                if ((sl > 0) && (tmp_line[sl - 1] == CONTINUE))
                {
                    tmp_line[sl - 1] = ' ';
                }

                done = (status == eCPP_EOF);
                if (!done)
                {
                    if (status != eCPP_OK)
                    {
                        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
                    }
                    else if (out)
                    {
                        fprintf(out, "%s\n", line);
                    }
                }

                srenew(pline, strlen(pline) + strlen(tmp_line) + 1);
                strcat(pline, tmp_line);
                sfree(tmp_line);
            }

            /* skip trailing and leading spaces and comment text */
            strip_comment(pline);
            trim(pline);

            /* if there is something left... */
            if (static_cast<int>(strlen(pline)) > 0)
            {
                if (pline[0] == OPENDIR)
                {
                    /* A directive on this line: copy the directive
                     * without the brackets into dirstr, then
                     * skip spaces and tabs on either side of directive
                     */
                    dirstr = gmx_strdup((pline + 1));
                    if ((dummy2 = strchr(dirstr, CLOSEDIR)) != nullptr)
                    {
                        (*dummy2) = 0;
                    }
                    trim(dirstr);

                    if ((newd = str2dir(dirstr)) == Directive::d_invalid)
                    {
                        sprintf(errbuf, "Invalid directive %s", dirstr);
                        wi->addError(errbuf);
                    }
                    else
                    {
                        /* Directive found */
                        if (DS_Check_Order(DS, newd))
                        {
                            DS_Push(&DS, newd);
                            d = newd;
                        }
                        else
                        {
                            /* we should print here which directives should have
                               been present, and which actually are */
                            gmx_fatal(FARGS,
                                      "%s\nInvalid order for directive %s",
                                      cpp_error(&handle, eCPP_SYNTAX),
                                      enumValueToString(newd));
                            /* d = Directive::d_invalid; */
                        }

                        if (d == Directive::d_intermolecular_interactions)
                        {
                            if (*intermolecular_interactions == nullptr)
                            {
                                /* We (mis)use the moleculetype processing
                                 * to process the intermolecular interactions
                                 * by making a "molecule" of the size of the system.
                                 */
                                *intermolecular_interactions = std::make_unique<MoleculeInformation>();
                                mi0                          = intermolecular_interactions->get();
                                mi0->initMolInfo();
                                make_atoms_sys(*molblock, *molinfo, &mi0->atoms);
                            }
                        }
                    }
                    sfree(dirstr);
                }
                else if (d != Directive::d_invalid)
                {
                    /* Not a directive, just a plain string
                     * use a gigantic switch to decode,
                     * if there is a valid directive!
                     */
                    switch (d)
                    {
                        case Directive::d_defaults:
                            if (bReadDefaults)
                            {
                                gmx_fatal(FARGS,
                                          "%s\nFound a second defaults directive.\n",
                                          cpp_error(&handle, eCPP_SYNTAX));
                            }
                            bReadDefaults = TRUE;
                            nscan         = sscanf(
                                    pline, "%s%s%s%lf%lf%lf", nb_str, comb_str, genpairs, &fLJ, &fQQ, &fPOW);
                            if (nscan < 2)
                            {
                                too_few(wi);
                            }
                            else
                            {
                                bGenPairs = FALSE;
                                fudgeLJ   = 1.0;
                                *fudgeQQ  = 1.0;

                                get_nbparm(nb_str, comb_str, &nb_funct, combination_rule, wi);
                                if (nscan >= 3)
                                {
                                    bGenPairs = (gmx::equalCaseInsensitive(genpairs, "Y", 1));
                                    if (nb_funct != VanDerWaalsPotential::LJ && bGenPairs)
                                    {
                                        gmx_fatal(FARGS,
                                                  "Generating pair parameters is only supported "
                                                  "with LJ non-bonded interactions");
                                    }
                                }
                                if (nscan >= 4)
                                {
                                    fudgeLJ = fLJ;
                                }
                                if (nscan >= 5)
                                {
                                    *fudgeQQ = fQQ;
                                }
                                if (nscan >= 6)
                                {
                                    *reppow = fPOW;
                                }
                            }
                            nb_funct = static_cast<VanDerWaalsPotential>(ifunc_index(
                                    Directive::d_nonbond_params, static_cast<int>(nb_funct)));

                            break;
                        case Directive::d_atomtypes:
                            push_at(atypes,
                                    &bondAtomType,
                                    pline,
                                    static_cast<int>(nb_funct),
                                    &nbparam,
                                    bGenPairs ? &pair : nullptr,
                                    wi);
                            break;

                        case Directive::d_bondtypes: // Intended to fall through
                        case Directive::d_constrainttypes:
                            push_bt(d, interactions, 2, nullptr, &bondAtomType, pline, wi);
                            break;
                        case Directive::d_pairtypes:
                            if (bGenPairs)
                            {
                                push_nbt(d, pair, atypes, pline, F_LJ14, wi);
                            }
                            else
                            {
                                push_bt(d, interactions, 2, atypes, nullptr, pline, wi);
                            }
                            break;
                        case Directive::d_angletypes:
                            push_bt(d, interactions, 3, nullptr, &bondAtomType, pline, wi);
                            break;
                        case Directive::d_dihedraltypes:
                            /* Special routine that can read both 2 and 4 atom dihedral definitions. */
                            push_dihedraltype(d, interactions, &bondAtomType, pline, wi);
                            break;

                        case Directive::d_nonbond_params:
                            push_nbt(d, nbparam, atypes, pline, static_cast<int>(nb_funct), wi);
                            break;

                        case Directive::d_implicit_genborn_params: // NOLINT bugprone-branch-clone
                            // Skip this line, so old topologies with
                            // GB parameters can be read.
                            break;

                        case Directive::d_implicit_surface_params:
                            // Skip this line, so that any topologies
                            // with surface parameters can be read
                            // (even though these were never formally
                            // supported).
                            break;

                        case Directive::d_cmaptypes:
                            push_cmaptype(d, interactions, NRAL(F_CMAP), atypes, &bondAtomType, pline, wi);
                            break;

                        case Directive::d_moleculetype:
                        {
                            if (!bReadMolType)
                            {
                                int ntype;
                                if (opts->couple_moltype != nullptr
                                    && (opts->couple_lam0 == ecouplamNONE || opts->couple_lam0 == ecouplamQ
                                        || opts->couple_lam1 == ecouplamNONE
                                        || opts->couple_lam1 == ecouplamQ))
                                {
                                    dcatt = add_atomtype_decoupled(
                                            atypes, &nbparam, bGenPairs ? &pair : nullptr);
                                }
                                ntype  = atypes->size();
                                ncombs = (ntype * (ntype + 1)) / 2;
                                generate_nbparams(*combination_rule,
                                                  static_cast<int>(nb_funct),
                                                  &(interactions[static_cast<int>(nb_funct)]),
                                                  atypes,
                                                  wi);
                                ncopy = copy_nbparams(nbparam,
                                                      static_cast<int>(nb_funct),
                                                      &(interactions[static_cast<int>(nb_funct)]),
                                                      ntype);
                                GMX_LOG(logger.info)
                                        .asParagraph()
                                        .appendTextFormatted(
                                                "Generated %d of the %d non-bonded parameter "
                                                "combinations",
                                                ncombs - ncopy,
                                                ncombs);
                                free_nbparam(nbparam, ntype);
                                if (bGenPairs)
                                {
                                    gen_pairs((interactions[static_cast<int>(nb_funct)]),
                                              &(interactions[F_LJ14]),
                                              fudgeLJ,
                                              *combination_rule);
                                    ncopy = copy_nbparams(
                                            pair, static_cast<int>(nb_funct), &(interactions[F_LJ14]), ntype);
                                    GMX_LOG(logger.info)
                                            .asParagraph()
                                            .appendTextFormatted(
                                                    "Generated %d of the %d 1-4 parameter "
                                                    "combinations",
                                                    ncombs - ncopy,
                                                    ncombs);
                                    free_nbparam(pair, ntype);
                                }
                                /* Copy GBSA parameters to atomtype array? */

                                bReadMolType = TRUE;
                            }

                            push_molt(symtab, molinfo, pline, wi);
                            exclusionBlocks.emplace_back();
                            mi0                    = &molinfo->back();
                            mi0->atoms.haveMass    = TRUE;
                            mi0->atoms.haveCharge  = TRUE;
                            mi0->atoms.haveType    = TRUE;
                            mi0->atoms.haveBState  = TRUE;
                            mi0->atoms.havePdbInfo = FALSE;
                            break;
                        }
                        case Directive::d_atoms:
                            push_atom(symtab, &(mi0->atoms), atypes, pline, wi);
                            break;

                        case Directive::d_pairs:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_bond(d,
                                      interactions,
                                      mi0->interactions,
                                      &(mi0->atoms),
                                      atypes,
                                      pline,
                                      FALSE,
                                      bGenPairs,
                                      *fudgeQQ,
                                      bZero,
                                      &bWarn_copy_A_B,
                                      wi);
                            break;
                        case Directive::d_pairs_nb:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_bond(d,
                                      interactions,
                                      mi0->interactions,
                                      &(mi0->atoms),
                                      atypes,
                                      pline,
                                      FALSE,
                                      FALSE,
                                      1.0,
                                      bZero,
                                      &bWarn_copy_A_B,
                                      wi);
                            break;

                        case Directive::d_vsites1:
                        case Directive::d_vsites2:
                        case Directive::d_vsites3:
                        case Directive::d_vsites4:
                        case Directive::d_bonds:
                        case Directive::d_angles:
                        case Directive::d_constraints:
                        case Directive::d_settles:
                        case Directive::d_position_restraints:
                        case Directive::d_angle_restraints:
                        case Directive::d_angle_restraints_z:
                        case Directive::d_distance_restraints:
                        case Directive::d_orientation_restraints:
                        case Directive::d_dihedral_restraints:
                        case Directive::d_dihedrals:
                        case Directive::d_polarization:
                        case Directive::d_water_polarization:
                        case Directive::d_thole_polarization:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_bond(d,
                                      interactions,
                                      mi0->interactions,
                                      &(mi0->atoms),
                                      atypes,
                                      pline,
                                      TRUE,
                                      bGenPairs,
                                      *fudgeQQ,
                                      bZero,
                                      &bWarn_copy_A_B,
                                      wi);
                            break;
                        case Directive::d_cmap:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_cmap(d, interactions, mi0->interactions, &(mi0->atoms), atypes, pline, wi);
                            break;

                        case Directive::d_vsitesn:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_vsitesn(d, mi0->interactions, &(mi0->atoms), pline, wi);
                            break;
                        case Directive::d_exclusions:
                            GMX_ASSERT(!exclusionBlocks.empty(),
                                       "exclusionBlocks must always be allocated so exclusions can "
                                       "be processed");
                            if (exclusionBlocks.back().empty())
                            {
                                GMX_RELEASE_ASSERT(mi0,
                                                   "Need to have a valid MoleculeInformation "
                                                   "object to work on");
                                exclusionBlocks.back().resize(mi0->atoms.nr);
                            }
                            push_excl(pline, exclusionBlocks.back(), wi);
                            break;
                        case Directive::d_system:
                            trim(pline);
                            title = put_symtab(symtab, pline);
                            break;
                        case Directive::d_molecules:
                        {
                            int  whichmol;
                            bool bCouple;

                            push_mol(*molinfo, pline, &whichmol, &nrcopies, wi);
                            mi0 = &((*molinfo)[whichmol]);
                            molblock->resize(molblock->size() + 1);
                            molblock->back().type = whichmol;
                            molblock->back().nmol = nrcopies;

                            bCouple = (opts->couple_moltype != nullptr
                                       && (gmx_strcasecmp("system", opts->couple_moltype) == 0
                                           || strcmp(*(mi0->name), opts->couple_moltype) == 0));
                            if (bCouple)
                            {
                                nmol_couple += nrcopies;
                            }

                            if (mi0->atoms.nr == 0)
                            {
                                gmx_fatal(FARGS, "Molecule type '%s' contains no atoms", *mi0->name);
                            }
                            GMX_LOG(logger.info)
                                    .asParagraph()
                                    .appendTextFormatted(
                                            "Excluding %d bonded neighbours molecule type '%s'",
                                            mi0->nrexcl,
                                            *mi0->name);
                            sum_q(&mi0->atoms, nrcopies, &qt, &qBt);
                            if (!mi0->bProcessed)
                            {
                                generate_excl(mi0->nrexcl, mi0->atoms.nr, mi0->interactions, &(mi0->excls));
                                gmx::mergeExclusions(&(mi0->excls), exclusionBlocks[whichmol]);
                                make_shake(mi0->interactions, &mi0->atoms, opts->nshake, logger);

                                if (bCouple)
                                {
                                    convert_moltype_couple(mi0,
                                                           dcatt,
                                                           *fudgeQQ,
                                                           opts->couple_lam0,
                                                           opts->couple_lam1,
                                                           opts->bCoupleIntra,
                                                           static_cast<int>(nb_funct),
                                                           &(interactions[static_cast<int>(nb_funct)]),
                                                           wi);
                                }
                                stupid_fill_block(&mi0->mols, mi0->atoms.nr, TRUE);
                                mi0->bProcessed = TRUE;
                            }
                            break;
                        }
                        default:
                            GMX_LOG(logger.warning)
                                    .asParagraph()
                                    .appendTextFormatted("case: %d", static_cast<int>(d));
                            gmx_incons("unknown directive");
                    }
                }
            }
            sfree(pline);
            pline = nullptr;
        }
    } while (!done);

    // Check that all strings defined with -D were used when processing topology
    std::string unusedDefineWarning = checkAndWarnForUnusedDefines(*handle);
    if (!unusedDefineWarning.empty())
    {
        wi->addWarning(unusedDefineWarning);
    }

    for (char* element : cpp_opts_return)
    {
        sfree(element);
    }

    if (out)
    {
        gmx_fio_fclose(out);
    }

    /* List of GROMACS define names for force fields that have been
     * parametrized using constraints involving hydrogens only.
     *
     * We should avoid hardcoded names, but this is hopefully only
     * needed temparorily for discouraging use of constraints=all-bonds.
     */
    const std::array<std::string, 3> ffDefines = { "_FF_AMBER", "_FF_CHARMM", "_FF_OPLSAA" };
    *ffParametrizedWithHBondConstraints        = false;
    for (const std::string& ffDefine : ffDefines)
    {
        if (cpp_find_define(&handle, ffDefine))
        {
            *ffParametrizedWithHBondConstraints = true;
        }
    }

    if (cpp_find_define(&handle, "_FF_GROMOS96") != nullptr)
    {
        wi->addWarning(
                "The GROMOS force fields have been parametrized with a physically incorrect "
                "multiple-time-stepping scheme for a twin-range cut-off. When used with "
                "a single-range cut-off (or a correct Trotter multiple-time-stepping scheme), "
                "physical properties, such as the density, might differ from the intended values. "
                "Since there are researchers actively working on validating GROMOS with modern "
                "integrators we have not yet removed the GROMOS force fields, but you should be "
                "aware of these issues and check if molecules in your system are affected before "
                "proceeding. "
                "Further information is available at "
                "https://gitlab.com/gromacs/gromacs/-/issues/2884, "
                "and a longer explanation of our decision to remove physically incorrect "
                "algorithms "
                "can be found at https://doi.org/10.26434/chemrxiv.11474583.v1 .");
    }
    // TODO: Update URL for Issue #2884 in conjunction with updating grompp.warn in regressiontests.

    cpp_done(handle);

    if (opts->couple_moltype)
    {
        if (nmol_couple == 0)
        {
            gmx_fatal(FARGS, "Did not find any molecules of type '%s' for coupling", opts->couple_moltype);
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "Coupling %d copies of molecule type '%s'", nmol_couple, opts->couple_moltype);
    }

    /* this is not very clean, but fixes core dump on empty system name */
    if (!title)
    {
        title = put_symtab(symtab, "");
    }

    if (std::fabs(qt) > 1e-4)
    {
        sprintf(warn_buf, "System has non-zero total charge: %.6f\n%s\n", qt, floating_point_arithmetic_tip);
        wi->addNote(warn_buf);
    }
    if (std::fabs(qBt) > 1e-4 && !gmx_within_tol(qBt, qt, 1e-6))
    {
        sprintf(warn_buf, "State B has non-zero total charge: %.6f\n%s\n", qBt, floating_point_arithmetic_tip);
        wi->addNote(warn_buf);
    }
    if (usingFullRangeElectrostatics && (std::fabs(qt) > 1e-4 || std::fabs(qBt) > 1e-4))
    {
        wi->addWarning(
                "You are using Ewald electrostatics in a system with net charge. This can lead to "
                "severe artifacts, such as ions moving into regions with low dielectric, due to "
                "the uniform background charge. We suggest to neutralize your system with counter "
                "ions, possibly in combination with a physiological salt concentration.");
        please_cite(stdout, "Hub2014a");
    }

    DS_Done(&DS);

    if (*intermolecular_interactions != nullptr)
    {
        sfree(intermolecular_interactions->get()->atoms.atom);
    }

    return title;
}

char** do_top(bool                                        bVerbose,
              const char*                                 topfile,
              const std::optional<std::filesystem::path>& topppfile,
              t_gromppopts*                               opts,
              bool                                        bZero,
              t_symtab*                                   symtab,
              gmx::ArrayRef<InteractionsOfType>           interactions,
              CombinationRule*                            combination_rule,
              double*                                     repulsion_power,
              real*                                       fudgeQQ,
              PreprocessingAtomTypes*                     atypes,
              std::vector<MoleculeInformation>*           molinfo,
              std::unique_ptr<MoleculeInformation>*       intermolecular_interactions,
              const t_inputrec*                           ir,
              std::vector<gmx_molblock_t>*                molblock,
              bool*                                       ffParametrizedWithHBondConstraints,
              WarningHandler*                             wi,
              const gmx::MDLogger&                        logger)
{
    char** title;

    if (bVerbose)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("processing topology...");
    }
    title = read_topol(topfile,
                       topppfile,
                       opts->define,
                       opts->include,
                       symtab,
                       atypes,
                       molinfo,
                       intermolecular_interactions,
                       interactions,
                       combination_rule,
                       repulsion_power,
                       opts,
                       fudgeQQ,
                       molblock,
                       ffParametrizedWithHBondConstraints,
                       ir->efep != FreeEnergyPerturbationType::No,
                       bZero,
                       usingFullElectrostatics(ir->coulombtype),
                       wi,
                       logger);

    if ((*combination_rule != CombinationRule::Geometric) && (ir->vdwtype == VanDerWaalsType::User))
    {
        wi->addWarning(
                "Using sigma/epsilon based combination rules with"
                " user supplied potential function may produce unwanted"
                " results");
    }

    return title;
}

/*! \brief
 * Exclude molecular interactions for QM atoms handled by MiMic.
 *
 * Update the exclusion lists to include all QM atoms of this molecule,
 * replace bonds between QM atoms with CONNBOND and
 * set charges of QM atoms to 0.
 *
 * \param[in,out] molt molecule type with QM atoms
 * \param[in] grpnr group informatio
 * \param[in,out] ir input record
 * \param[in] logger Handle to logging interface.
 */
static void generate_qmexcl_moltype(gmx_moltype_t*       molt,
                                    const unsigned char* grpnr,
                                    t_inputrec*          ir,
                                    const gmx::MDLogger& logger)
{
    /* This routine expects molt->ilist to be of size F_NRE and ordered. */

    /* generates the exclusions between the individual QM atoms, as
     * these interactions should be handled by the QM subroutines and
     * not by the gromacs routines
     */
    int   qm_max = 0, qm_nr = 0, link_nr = 0;
    int * qm_arr = nullptr, *link_arr = nullptr;
    bool* bQMMM;

    /* First we search and select the QM atoms in an qm_arr array that
     * we use to create the exclusions.
     *
     * we take the possibility into account that a user has defined more
     * than one QM group:
     *
     * for that we also need to do this an ugly work-about just in case
     * the QM group contains the entire system...
     */

    /* we first search for all the QM atoms and put them in an array
     */
    for (int j = 0; j < ir->opts.ngQM; j++)
    {
        for (int i = 0; i < molt->atoms.nr; i++)
        {
            if (qm_nr >= qm_max)
            {
                qm_max += 100;
                srenew(qm_arr, qm_max);
            }
            if ((grpnr ? grpnr[i] : 0) == j)
            {
                qm_arr[qm_nr++]        = i;
                molt->atoms.atom[i].q  = 0.0;
                molt->atoms.atom[i].qB = 0.0;
            }
        }
    }
    /* bQMMM[..] is an array containin TRUE/FALSE for atoms that are
     * QM/not QM. We first set all elements to false. Afterwards we use
     * the qm_arr to change the elements corresponding to the QM atoms
     * to TRUE.
     */
    snew(bQMMM, molt->atoms.nr);
    for (int i = 0; i < molt->atoms.nr; i++)
    {
        bQMMM[i] = FALSE;
    }
    for (int i = 0; i < qm_nr; i++)
    {
        bQMMM[qm_arr[i]] = TRUE;
    }

    /* We remove all bonded interactions (i.e. bonds,
     * angles, dihedrals, 1-4's), involving the QM atoms. The way they
     * are removed is as follows: if the interaction invloves 2 atoms,
     * it is removed if both atoms are QMatoms. If it involves 3 atoms,
     * it is removed if at least two of the atoms are QM atoms, if the
     * interaction involves 4 atoms, it is removed if there are at least
     * 2 QM atoms.  Since this routine is called once before any forces
     * are computed, the top->idef.il[N].iatom[] array (see idef.h) can
     * be rewritten at this poitn without any problem. 25-9-2002 */

    /* first check whether we already have CONNBONDS.
     * Note that if we don't, we don't add a param entry and set ftype=0,
     * which is ok, since CONNBONDS does not use parameters.
     */
    int ftype_connbond = 0;
    int ind_connbond   = 0;
    if (!molt->ilist[F_CONNBONDS].empty())
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("nr. of CONNBONDS present already: %d",
                                     molt->ilist[F_CONNBONDS].size() / 3);
        ftype_connbond = molt->ilist[F_CONNBONDS].iatoms[0];
        ind_connbond   = molt->ilist[F_CONNBONDS].size();
    }
    /* now we delete all bonded interactions, except the ones describing
     * a chemical bond. These are converted to CONNBONDS
     */
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (!(interaction_function[ftype].flags & IF_BOND) || ftype == F_CONNBONDS)
        {
            continue;
        }
        int nratoms = interaction_function[ftype].nratoms;
        int j       = 0;
        while (j < molt->ilist[ftype].size())
        {
            bool bexcl;

            if (nratoms == 2)
            {
                /* Remove an interaction between two atoms when both are
                 * in the QM region. Note that we don't have to worry about
                 * link atoms here, as they won't have 2-atom interactions.
                 */
                int a1 = molt->ilist[ftype].iatoms[1 + j + 0];
                int a2 = molt->ilist[ftype].iatoms[1 + j + 1];
                bexcl  = (bQMMM[a1] && bQMMM[a2]);
                /* A chemical bond between two QM atoms will be copied to
                 * the F_CONNBONDS list, for reasons mentioned above.
                 */
                if (bexcl && IS_CHEMBOND(ftype))
                {
                    InteractionList& ilist = molt->ilist[F_CONNBONDS];
                    ilist.iatoms.resize(ind_connbond + 3);
                    ilist.iatoms[ind_connbond++] = ftype_connbond;
                    ilist.iatoms[ind_connbond++] = a1;
                    ilist.iatoms[ind_connbond++] = a2;
                }
            }
            else
            {
                /* MM interactions have to be excluded if they are included
                 * in the QM already. Because we use a link atom (H atom)
                 * when the QM/MM boundary runs through a chemical bond, this
                 * means that as long as one atom is MM, we still exclude,
                 * as the interaction is included in the QM via:
                 * QMatom1-QMatom2-QMatom-3-Linkatom.
                 */
                int numQmAtoms = 0;
                for (int jj = j + 1; jj < j + 1 + nratoms; jj++)
                {
                    if (bQMMM[molt->ilist[ftype].iatoms[jj]])
                    {
                        numQmAtoms++;
                    }
                }

                /* MiMiC treats link atoms as quantum atoms - therefore
                 * we do not need do additional exclusions here */
                bexcl = numQmAtoms == nratoms;

                if (bexcl && ftype == F_SETTLE)
                {
                    gmx_fatal(FARGS,
                              "Can not apply QM to molecules with SETTLE, replace the moleculetype "
                              "using QM and SETTLE by one without SETTLE");
                }
            }
            if (bexcl)
            {
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                InteractionList& ilist = molt->ilist[ftype];
                for (int k = j; k < ilist.size() - (nratoms + 1); k++)
                {
                    ilist.iatoms[k] = ilist.iatoms[k + (nratoms + 1)];
                }
                ilist.iatoms.resize(ilist.size() - (nratoms + 1));
            }
            else
            {
                j += nratoms + 1; /* the +1 is for the functype */
            }
        }
    }
    /* Now, we search for atoms bonded to a QM atom because we also want
     * to exclude their nonbonded interactions with the QM atoms. The
     * reason for this is that this interaction is accounted for in the
     * linkatoms interaction with the QMatoms and would be counted
     * twice.  */

    /* creating the exclusion block for the QM atoms. Each QM atom has
     * as excluded elements all the other QMatoms (and itself).
     */
    t_blocka qmexcl;
    qmexcl.nr  = molt->atoms.nr;
    qmexcl.nra = qm_nr * (qm_nr + link_nr) + link_nr * qm_nr;
    snew(qmexcl.index, qmexcl.nr + 1);
    snew(qmexcl.a, qmexcl.nra);
    int j = 0;
    for (int i = 0; i < qmexcl.nr; i++)
    {
        qmexcl.index[i] = j;
        if (bQMMM[i])
        {
            for (int k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k + j] = qm_arr[k];
            }
            for (int k = 0; k < link_nr; k++)
            {
                qmexcl.a[qm_nr + k + j] = link_arr[k];
            }
            j += (qm_nr + link_nr);
        }
    }
    qmexcl.index[qmexcl.nr] = j;

    /* and merging with the exclusions already present in sys.
     */

    std::vector<gmx::ExclusionBlock> qmexcl2(molt->atoms.nr);
    gmx::blockaToExclusionBlocks(&qmexcl, qmexcl2);
    gmx::mergeExclusions(&(molt->excls), qmexcl2);

    /* Finally, we also need to get rid of the pair interactions of the
     * classical atom bonded to the boundary QM atoms with the QMatoms,
     * as this interaction is already accounted for by the QM, so also
     * here we run the risk of double counting! We proceed in a similar
     * way as we did above for the other bonded interactions: */
    for (int i = F_LJ14; i < F_COUL14; i++)
    {
        int nratoms = interaction_function[i].nratoms;
        int j       = 0;
        while (j < molt->ilist[i].size())
        {
            int  a1    = molt->ilist[i].iatoms[j + 1];
            int  a2    = molt->ilist[i].iatoms[j + 2];
            bool bexcl = (bQMMM[a1] && bQMMM[a2]);
            if (bexcl)
            {
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                InteractionList& ilist = molt->ilist[i];
                for (int k = j; k < ilist.size() - (nratoms + 1); k++)
                {
                    ilist.iatoms[k] = ilist.iatoms[k + (nratoms + 1)];
                }
                ilist.iatoms.resize(ilist.size() - (nratoms + 1));
            }
            else
            {
                j += nratoms + 1; /* the +1 is for the functype */
            }
        }
    }

    free(qm_arr);
    free(bQMMM);
    free(link_arr);
} /* generate_qmexcl */

void generate_qmexcl(gmx_mtop_t* sys, t_inputrec* ir, const gmx::MDLogger& logger)
{
    /* This routine expects molt->molt[m].ilist to be of size F_NRE and ordered.
     */

    unsigned char*  grpnr;
    int             mol, nat_mol;
    gmx_molblock_t* molb;
    bool            bQMMM;

    grpnr = sys->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics].data();

    for (size_t mb = 0; mb < sys->molblock.size(); mb++)
    {
        molb    = &sys->molblock[mb];
        nat_mol = sys->moltype[molb->type].atoms.nr;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            bQMMM = FALSE;
            for (int i = 0; i < nat_mol; i++)
            {
                if ((grpnr ? grpnr[i] : 0) < (ir->opts.ngQM))
                {
                    bQMMM = TRUE;
                }
            }

            if (bQMMM)
            {
                if (molb->nmol > 1)
                {
                    /* We need to split this molblock */
                    if (mol > 0)
                    {
                        /* Split the molblock at this molecule */
                        auto pos = sys->molblock.begin() + mb + 1;
                        sys->molblock.insert(pos, sys->molblock[mb]);
                        sys->molblock[mb].nmol = mol;
                        sys->molblock[mb + 1].nmol -= mol;
                        mb++;
                        molb = &sys->molblock[mb];
                    }
                    if (molb->nmol > 1)
                    {
                        /* Split the molblock after this molecule */
                        auto pos = sys->molblock.begin() + mb + 1;
                        sys->molblock.insert(pos, sys->molblock[mb]);
                        molb                   = &sys->molblock[mb];
                        sys->molblock[mb].nmol = 1;
                        sys->molblock[mb + 1].nmol -= 1;
                    }

                    /* Create a copy of a moltype for a molecule
                     * containing QM atoms and append it in the end of the list
                     */
                    std::vector<gmx_moltype_t> temp(sys->moltype.size());
                    for (size_t i = 0; i < sys->moltype.size(); ++i)
                    {
                        copy_moltype(&sys->moltype[i], &temp[i]);
                    }
                    sys->moltype.resize(sys->moltype.size() + 1);
                    for (size_t i = 0; i < temp.size(); ++i)
                    {
                        copy_moltype(&temp[i], &sys->moltype[i]);
                    }
                    copy_moltype(&sys->moltype[molb->type], &sys->moltype.back());
                    /* Copy the exclusions to a new array, since this is the only
                     * thing that needs to be modified for QMMM.
                     */
                    sys->moltype.back().excls = sys->moltype[molb->type].excls;
                    /* Set the molecule type for the QMMM molblock */
                    molb->type = sys->moltype.size() - 1;
                }
                generate_qmexcl_moltype(&sys->moltype[molb->type], grpnr, ir, logger);
            }
            if (grpnr)
            {
                grpnr += nat_mol;
            }
        }
    }
}
