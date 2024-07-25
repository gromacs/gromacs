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

#include "toppush.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringtoenumvalueconverter.h"
#include "gromacs/utility/stringutil.h"

void generate_nbparams(CombinationRule         comb,
                       int                     ftype,
                       InteractionsOfType*     interactions,
                       PreprocessingAtomTypes* atypes,
                       WarningHandler*         wi)
{
    constexpr int c_nrfp2 = 2;

    int  nr, nrfp;
    real c, bi, bj, ci, cj, ci0, ci1, ci2, cj0, cj1, cj2;

    /* Lean mean shortcuts */
    nr   = atypes->size();
    nrfp = NRFP(ftype);
    interactions->interactionTypes.clear();

    std::array<real, MAXFORCEPARAM> forceParam = { NOTSET };
    /* Fill the matrix with force parameters */
    // Prefetch the parameters to improve cache hits and avoid dereference and call overhead
    std::vector<std::pair<real, real>> cPrefetch;
    cPrefetch.reserve(nr);
    for (int i = 0; i < nr; i++)
    {
        cPrefetch.emplace_back(*atypes->atomNonBondedParamFromAtomType(i, 0),
                               *atypes->atomNonBondedParamFromAtomType(i, 1));
    }
    interactions->interactionTypes.reserve(nr * nr);
    switch (ftype)
    {
        case F_LJ:
            switch (comb)
            {
                case CombinationRule::Geometric:
                    // Geometric combination rules, c6 and c12 are independent
                    GMX_RELEASE_ASSERT(nrfp == c_nrfp2, "nfrp should be 2");
                    for (int i = 0; (i < nr); i++)
                    {
                        for (int j = 0; (j < nr); j++)
                        {
                            for (int nf = 0; (nf < c_nrfp2); nf++)
                            {
                                ci = (nf == 0 ? cPrefetch[i].first : cPrefetch[i].second);
                                cj = (nf == 0 ? cPrefetch[j].first : cPrefetch[j].second);
                                c  = std::sqrt(ci * cj);
                                forceParam[nf] = c;
                            }
                            interactions->interactionTypes.emplace_back(gmx::ArrayRef<const int>{},
                                                                        forceParam);
                        }
                    }
                    break;

                case CombinationRule::Arithmetic:
                    /* c0 and c1 are sigma and epsilon */
                    for (int i = 0; (i < nr); i++)
                    {
                        for (int j = 0; (j < nr); j++)
                        {
                            ci0           = cPrefetch[i].first;
                            cj0           = cPrefetch[j].first;
                            ci1           = cPrefetch[i].second;
                            cj1           = cPrefetch[j].second;
                            forceParam[0] = (std::fabs(ci0) + std::fabs(cj0)) * 0.5;
                            /* Negative sigma signals that c6 should be set to zero later,
                             * so we need to propagate that through the combination rules.
                             */
                            if (ci0 < 0 || cj0 < 0)
                            {
                                forceParam[0] *= -1;
                            }
                            forceParam[1] = std::sqrt(ci1 * cj1);
                            interactions->interactionTypes.emplace_back(gmx::ArrayRef<const int>{},
                                                                        forceParam);
                        }
                    }

                    break;
                case CombinationRule::GeomSigEps:
                    /* c0 and c1 are sigma and epsilon */
                    for (int i = 0; (i < nr); i++)
                    {
                        for (int j = 0; (j < nr); j++)
                        {
                            ci0           = cPrefetch[i].first;
                            cj0           = cPrefetch[j].first;
                            ci1           = cPrefetch[i].second;
                            cj1           = cPrefetch[j].second;
                            forceParam[0] = std::sqrt(std::fabs(ci0 * cj0));
                            /* Negative sigma signals that c6 should be set to zero later,
                             * so we need to propagate that through the combination rules.
                             */
                            if (ci0 < 0 || cj0 < 0)
                            {
                                forceParam[0] *= -1;
                            }
                            forceParam[1] = std::sqrt(ci1 * cj1);
                            interactions->interactionTypes.emplace_back(gmx::ArrayRef<const int>{},
                                                                        forceParam);
                        }
                    }

                    break;
                default:
                    auto message =
                            gmx::formatString("No such combination rule %s", enumValueToString(comb));
                    warning_error_and_exit(wi, message, FARGS);
            }
            break;

        case F_BHAM:
            /* Buckingham rules */
            for (int i = 0; (i < nr); i++)
            {
                for (int j = 0; (j < nr); j++)
                {
                    ci0           = cPrefetch[i].first;
                    cj0           = cPrefetch[j].first;
                    ci2           = *atypes->atomNonBondedParamFromAtomType(i, 2);
                    cj2           = *atypes->atomNonBondedParamFromAtomType(j, 2);
                    bi            = cPrefetch[i].second;
                    bj            = cPrefetch[j].second;
                    forceParam[0] = std::sqrt(ci0 * cj0);
                    if ((bi == 0) || (bj == 0))
                    {
                        forceParam[1] = 0;
                    }
                    else
                    {
                        forceParam[1] = 2.0 / (1 / bi + 1 / bj);
                    }
                    forceParam[2] = std::sqrt(ci2 * cj2);
                    interactions->interactionTypes.emplace_back(gmx::ArrayRef<const int>{}, forceParam);
                }
            }

            break;
        default:
            auto message = gmx::formatString("Invalid nonbonded type %s",
                                             interaction_function[ftype].longname);
            wi->addError(message);
    }
}

/*! \brief Used to temporarily store the explicit non-bonded parameter
 * combinations, which will be copied to InteractionsOfType. */
struct t_nbparam
{
    //! Has this combination been set.
    bool bSet;
    //! The non-bonded parameters
    real c[4];
};

static void realloc_nb_params(PreprocessingAtomTypes* atypes, t_nbparam*** nbparam, t_nbparam*** pair)
{
    /* Add space in the non-bonded parameters matrix */
    int atnr = atypes->size();
    srenew(*nbparam, atnr);
    snew((*nbparam)[atnr - 1], atnr);
    if (pair)
    {
        srenew(*pair, atnr);
        snew((*pair)[atnr - 1], atnr);
    }
}

int copy_nbparams(t_nbparam** param, int ftype, InteractionsOfType* interactions, int nr)
{
    int nrfp, ncopy;

    nrfp = NRFP(ftype);

    ncopy = 0;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            GMX_RELEASE_ASSERT(param, "Must have valid parameters");
            if (param[i][j].bSet)
            {
                for (int f = 0; f < nrfp; f++)
                {
                    interactions->interactionTypes[nr * i + j].setForceParameter(f, param[i][j].c[f]);
                    interactions->interactionTypes[nr * j + i].setForceParameter(f, param[i][j].c[f]);
                }
                ncopy++;
            }
        }
    }

    return ncopy;
}

void free_nbparam(t_nbparam** param, int nr)
{
    int i;

    GMX_RELEASE_ASSERT(param, "Must have valid parameters");
    for (i = 0; i < nr; i++)
    {
        GMX_RELEASE_ASSERT(param[i], "Must have valid parameters");
        sfree(param[i]);
    }
    sfree(param);
}

static void copy_B_from_A(int ftype, double* c)
{
    int nrfpA, nrfpB, i;

    nrfpA = NRFPA(ftype);
    nrfpB = NRFPB(ftype);

    /* Copy the B parameters from the first nrfpB A parameters */
    for (i = 0; (i < nrfpB); i++)
    {
        c[nrfpA + i] = c[i];
    }
}

//! Local definition that supersedes the central one, as we only want the leading letter
static const char* enumValueToLetterAsString(ParticleType enumValue)
{
    static constexpr gmx::EnumerationArray<ParticleType, const char*> particleTypeLetters = {
        "A", "N", "S", "B", "V"
    };
    return particleTypeLetters[enumValue];
}

void push_at(PreprocessingAtomTypes*    at,
             PreprocessingBondAtomType* bondAtomType,
             char*                      line,
             int                        nb_funct,
             t_nbparam***               nbparam,
             t_nbparam***               pair,
             WarningHandler*            wi)
{
    int     nfields, nfp0 = -1;
    int     nread;
    char    type[STRLEN], btype[STRLEN], ptype[STRLEN];
    double  m, q;
    double  c[MAXFORCEPARAM];
    char    tmpfield[12][100]; /* Max 12 fields of width 100 */
    t_atom* atom;
    int     atomnr;
    bool    have_atomic_number;
    bool    have_bonded_type;

    snew(atom, 1);

    /* First assign input line to temporary array */
    nfields = sscanf(line,
                     "%s%s%s%s%s%s%s%s%s%s%s%s",
                     tmpfield[0],
                     tmpfield[1],
                     tmpfield[2],
                     tmpfield[3],
                     tmpfield[4],
                     tmpfield[5],
                     tmpfield[6],
                     tmpfield[7],
                     tmpfield[8],
                     tmpfield[9],
                     tmpfield[10],
                     tmpfield[11]);

    /* Comments on optional fields in the atomtypes section:
     *
     * The force field format is getting a bit old. For OPLS-AA we needed
     * to add a special bonded atomtype, and for Gerrit Groenhofs QM/MM stuff
     * we also needed the atomic numbers.
     * To avoid making all old or user-generated force fields unusable we
     * have introduced both these quantities as optional columns, and do some
     * acrobatics to check whether they are present or not.
     * This will all look much nicer when we switch to XML... sigh.
     *
     * Field 0 (mandatory) is the nonbonded type name. (string)
     * Field 1 (optional)  is the bonded type (string)
     * Field 2 (optional)  is the atomic number (int)
     * Field 3 (mandatory) is the mass (numerical)
     * Field 4 (mandatory) is the charge (numerical)
     * Field 5 (mandatory) is the particle type (single character)
     * This is followed by a number of nonbonded parameters.
     *
     * The safest way to identify the format is the particle type field.
     *
     * So, here is what we do:
     *
     * A. Read in the first six fields as strings
     * B. If field 3 (starting from 0) is a single char, we have neither
     *    bonded_type or atomic numbers.
     * C. If field 5 is a single char we have both.
     * D. If field 4 is a single char we check field 1. If this begins with
     *    an alphabetical character we have bonded types, otherwise atomic numbers.
     */

    if (nfields < 6)
    {
        too_few(wi);
        return;
    }

    if ((strlen(tmpfield[5]) == 1) && isalpha(tmpfield[5][0]))
    {
        have_bonded_type   = TRUE;
        have_atomic_number = TRUE;
    }
    else if ((strlen(tmpfield[3]) == 1) && isalpha(tmpfield[3][0]))
    {
        have_bonded_type   = FALSE;
        have_atomic_number = FALSE;
    }
    else
    {
        // Attempt parsing field 1 to integer. If successful, *end == '\0'
        char* end;
        strtol(tmpfield[1], &end, 10);

        // If conversion fails, we do not have an atomic number but a bonded type
        have_bonded_type   = (*end != 0);
        have_atomic_number = !have_bonded_type;
    }

    /* optional fields */
    atomnr = -1;

    switch (nb_funct)
    {

        case F_LJ:
            nfp0 = 2;

            if (have_atomic_number)
            {
                if (have_bonded_type)
                {
                    nread = sscanf(
                            line, "%s%s%d%lf%lf%s%lf%lf", type, btype, &atomnr, &m, &q, ptype, &c[0], &c[1]);
                    if (nread < 8)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%d%lf%lf%s%lf%lf", type, &atomnr, &m, &q, ptype, &c[0], &c[1]);
                    if (nread < 7)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }
            else
            {
                if (have_bonded_type)
                {
                    /* !have_atomic_number && have_bonded_type */
                    nread = sscanf(line, "%s%s%lf%lf%s%lf%lf", type, btype, &m, &q, ptype, &c[0], &c[1]);
                    if (nread < 7)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* !have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%lf%lf%s%lf%lf", type, &m, &q, ptype, &c[0], &c[1]);
                    if (nread < 6)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }

            if (!have_bonded_type)
            {
                strcpy(btype, type);
            }

            if (!have_atomic_number)
            {
                atomnr = -1;
            }

            break;

        case F_BHAM:
            nfp0 = 3;

            if (have_atomic_number)
            {
                if (have_bonded_type)
                {
                    nread = sscanf(
                            line, "%s%s%d%lf%lf%s%lf%lf%lf", type, btype, &atomnr, &m, &q, ptype, &c[0], &c[1], &c[2]);
                    if (nread < 9)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* have_atomic_number && !have_bonded_type */
                    nread = sscanf(
                            line, "%s%d%lf%lf%s%lf%lf%lf", type, &atomnr, &m, &q, ptype, &c[0], &c[1], &c[2]);
                    if (nread < 8)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }
            else
            {
                if (have_bonded_type)
                {
                    /* !have_atomic_number && have_bonded_type */
                    nread = sscanf(
                            line, "%s%s%lf%lf%s%lf%lf%lf", type, btype, &m, &q, ptype, &c[0], &c[1], &c[2]);
                    if (nread < 8)
                    {
                        too_few(wi);
                        return;
                    }
                }
                else
                {
                    /* !have_atomic_number && !have_bonded_type */
                    nread = sscanf(line, "%s%lf%lf%s%lf%lf%lf", type, &m, &q, ptype, &c[0], &c[1], &c[2]);
                    if (nread < 7)
                    {
                        too_few(wi);
                        return;
                    }
                }
            }

            if (!have_bonded_type)
            {
                strcpy(btype, type);
            }

            if (!have_atomic_number)
            {
                atomnr = -1;
            }

            break;

        default:
            auto message = gmx::formatString("Invalid function type %d in push_at", nb_funct);
            warning_error_and_exit(wi, message, FARGS);
    }
    for (int j = nfp0; (j < MAXFORCEPARAM); j++)
    {
        c[j] = 0.0;
    }
    std::array<real, MAXFORCEPARAM> forceParam;

    if (strlen(type) == 1 && isdigit(type[0]))
    {
        warning_error_and_exit(wi, "Atom type names can't be single digits.", FARGS);
    }

    if (strlen(btype) == 1 && isdigit(btype[0]))
    {
        warning_error_and_exit(wi, "Bond atom type names can't be single digits.", FARGS);
    }

    /* Hack to read old topologies */
    if (gmx_strcasecmp(ptype, "D") == 0)
    {
        sprintf(ptype, "V");
    }
    static const gmx::StringToEnumValueConverter<ParticleType, enumValueToLetterAsString, gmx::StringCompareType::CaseInsensitive, gmx::StripStrings::No>
                                s_stringToParticleType;
    std::optional<ParticleType> pt = s_stringToParticleType.valueFrom(ptype);
    if (!pt)
    {
        auto message = gmx::formatString("Invalid particle type %s", ptype);
        warning_error_and_exit(wi, message, FARGS);
    }

    atom->q     = q;
    atom->m     = m;
    atom->ptype = pt.value();
    for (int i = 0; i < MAXFORCEPARAM; i++)
    {
        forceParam[i] = c[i];
    }

    InteractionOfType interactionType({}, forceParam, "");

    auto batype_nr = bondAtomType->addBondAtomType(btype);

    auto atomType = at->atomTypeFromName(type);
    if (atomType.has_value())
    {
        auto message = gmx::formatString(
                "Atomtype %s was defined previously (e.g. in the forcefield files), "
                "and has now been defined again. This could happen e.g. if you would "
                "use a self-contained molecule .itp file that duplicates or replaces "
                "the contents of the standard force-field files. You should check "
                "the contents of your files and remove such repetition. If you know "
                "you should override the previous definition, then you could choose "
                "to suppress this warning with -maxwarn.",
                type);
        wi->addWarning(message);
        auto newAtomType = at->setType(*atomType, *atom, type, interactionType, batype_nr, atomnr);
        if (!newAtomType.has_value())
        {
            auto message = gmx::formatString("Replacing atomtype %s failed", type);
            warning_error_and_exit(wi, message, FARGS);
        }
    }
    else
    {
        at->addType(*atom, type, interactionType, batype_nr, atomnr);
        /* Add space in the non-bonded parameters matrix */
        realloc_nb_params(at, nbparam, pair);
    }
    sfree(atom);
}

//! Return whether the contents of \c a and \c b are the same, considering also reversed order.
template<typename T>
static bool equalEitherForwardOrBackward(gmx::ArrayRef<const T> a, gmx::ArrayRef<const T> b)
{
    return (std::equal(a.begin(), a.end(), b.begin()) || std::equal(a.begin(), a.end(), b.rbegin()));
}

static void push_bondtype(InteractionsOfType*      bt,
                          const InteractionOfType& b,
                          int                      nral,
                          int                      ftype,
                          bool                     bAllowRepeat,
                          const char*              line,
                          WarningHandler*          wi)
{
    int nr   = bt->size();
    int nrfp = NRFP(ftype);

    /* If bAllowRepeat is TRUE, we allow multiple entries as long as they
       are on directly _adjacent_ lines.
     */

    /* First check if our atomtypes are _identical_ (not reversed) to the previous
       entry. If they are not identical we search for earlier duplicates. If they are
       we can skip it, since we already searched for the first line
       in this group.
     */

    bool isContinuationOfBlock = false;
    if (bAllowRepeat && nr > 1)
    {
        isContinuationOfBlock               = true;
        gmx::ArrayRef<const int> newParAtom = b.atoms();
        gmx::ArrayRef<const int> sysParAtom = bt->interactionTypes[nr - 2].atoms();
        for (int j = 0; j < nral; j++)
        {
            if (newParAtom[j] != sysParAtom[j])
            {
                isContinuationOfBlock = false;
            }
        }
    }

    /* Search for earlier duplicates if this entry was not a continuation
       from the previous line.
     */
    bool addBondType = true;
    bool haveWarned  = false;
    bool haveErrored = false;
    for (int i = 0; (i < nr); i++)
    {
        gmx::ArrayRef<const int> bParams    = b.atoms();
        gmx::ArrayRef<const int> testParams = bt->interactionTypes[i].atoms();
        GMX_RELEASE_ASSERT(bParams.size() == testParams.size(),
                           "Number of atoms needs to be the same between parameters");
        if (equalEitherForwardOrBackward(bParams, testParams))
        {
            GMX_ASSERT(nrfp <= MAXFORCEPARAM,
                       "This is ensured in other places, but we need this assert to keep the clang "
                       "analyzer happy");
            const bool identicalParameters = std::equal(bt->interactionTypes[i].forceParam().begin(),
                                                        bt->interactionTypes[i].forceParam().begin() + nrfp,
                                                        b.forceParam().begin());

            if (!bAllowRepeat || identicalParameters)
            {
                addBondType = false;
            }

            if (!identicalParameters)
            {
                if (bAllowRepeat)
                {
                    /* With dihedral type 9 we only allow for repeating
                     * of the same parameters with blocks with 1 entry.
                     * Allowing overriding is too complex to check.
                     */
                    if (!isContinuationOfBlock && !haveErrored)
                    {
                        wi->addError(
                                "Encountered a second block of parameters for dihedral "
                                "type 9 for the same atoms, with either different parameters "
                                "and/or the first block has multiple lines. This is not "
                                "supported.");
                        haveErrored = true;
                    }
                }
                else if (!haveWarned)
                {
                    auto message = gmx::formatString(
                            "Bondtype %s was defined previously (e.g. in the forcefield files), "
                            "and has now been defined again. This could happen e.g. if you would "
                            "use a self-contained molecule .itp file that duplicates or replaces "
                            "the contents of the standard force-field files. You should check "
                            "the contents of your files and remove such repetition. If you know "
                            "you should override the previous definition, then you could choose "
                            "to suppress this warning with -maxwarn.%s",
                            interaction_function[ftype].longname,
                            (ftype == F_PDIHS) ? "\nUse dihedraltype 9 to allow several "
                                                 "multiplicity terms. Only consecutive "
                                                 "lines are combined. Non-consective lines "
                                                 "overwrite each other."
                                               : "");
                    wi->addWarning(message);

                    fprintf(stderr, "  old:                                         ");
                    gmx::ArrayRef<const real> forceParam = bt->interactionTypes[i].forceParam();
                    for (int j = 0; j < nrfp; j++)
                    {
                        fprintf(stderr, " %g", forceParam[j]);
                    }
                    fprintf(stderr, " \n  new: %s\n\n", line);

                    haveWarned = true;
                }
            }

            if (!identicalParameters && !bAllowRepeat)
            {
                /* Overwrite the parameters with the latest ones */
                // TODO considering improving the following code by replacing with:
                // std::copy(b->c, b->c + nrfp, bt->param[i].c);
                gmx::ArrayRef<const real> forceParam = b.forceParam();
                for (int j = 0; j < nrfp; j++)
                {
                    bt->interactionTypes[i].setForceParameter(j, forceParam[j]);
                }
            }
        }
    }

    if (addBondType)
    {
        /* fill the arrays up and down */
        bt->interactionTypes.emplace_back(b.atoms(), b.forceParam(), b.interactionTypeName());
        /* need to store force values because they might change below */
        std::vector<real> forceParam(b.forceParam().begin(), b.forceParam().end());

        /* The definitions of linear angles depend on the order of atoms,
         * that means that for atoms i-j-k, with certain parameter a, the
         * corresponding k-j-i angle will have parameter 1-a.
         */
        if (ftype == F_LINEAR_ANGLES)
        {
            forceParam[0] = 1 - forceParam[0];
            forceParam[2] = 1 - forceParam[2];
        }
        std::vector<int>         atoms;
        gmx::ArrayRef<const int> oldAtoms = b.atoms();
        for (auto oldAtom = oldAtoms.rbegin(); oldAtom != oldAtoms.rend(); oldAtom++)
        {
            atoms.emplace_back(*oldAtom);
        }
        bt->interactionTypes.emplace_back(atoms, forceParam, b.interactionTypeName());
    }
}

static std::vector<int> atomTypesFromAtomNames(const PreprocessingAtomTypes*    atomTypes,
                                               const PreprocessingBondAtomType* bondAtomTypes,
                                               gmx::ArrayRef<const char[20]>    atomNames,
                                               WarningHandler*                  wi)
{

    GMX_RELEASE_ASSERT(!(!atomNames.empty() && !atomTypes && !bondAtomTypes),
                       "Need to have either valid atomtypes or bondatomtypes object");

    std::vector<int> atomTypesFromAtomNames;
    for (const auto& name : atomNames)
    {
        if (atomTypes != nullptr)
        {
            auto atomType = atomTypes->atomTypeFromName(name);
            if (!atomType.has_value())
            {
                auto message = gmx::formatString("Unknown atomtype %s\n", name);
                warning_error_and_exit(wi, message, FARGS);
            }
            atomTypesFromAtomNames.emplace_back(*atomType);
        }
        else if (bondAtomTypes != nullptr)
        {
            auto bondAtomType = bondAtomTypes->bondAtomTypeFromName(name);
            if (!bondAtomType.has_value())
            {
                auto message = gmx::formatString("Unknown bond_atomtype %s\n", name);
                warning_error_and_exit(wi, message, FARGS);
            }
            atomTypesFromAtomNames.emplace_back(*bondAtomType);
        }
    }
    return atomTypesFromAtomNames;
}


void push_bt(Directive                         d,
             gmx::ArrayRef<InteractionsOfType> bt,
             int                               nral,
             PreprocessingAtomTypes*           at,
             PreprocessingBondAtomType*        bondAtomType,
             char*                             line,
             WarningHandler*                   wi)
{
    const char* formal[MAXATOMLIST + 1] = {
        "%s", "%s%s", "%s%s%s", "%s%s%s%s", "%s%s%s%s%s", "%s%s%s%s%s%s", "%s%s%s%s%s%s%s"
    };
    const char* formnl[MAXATOMLIST + 1] = { "%*s",
                                            "%*s%*s",
                                            "%*s%*s%*s",
                                            "%*s%*s%*s%*s",
                                            "%*s%*s%*s%*s%*s",
                                            "%*s%*s%*s%*s%*s%*s",
                                            "%*s%*s%*s%*s%*s%*s%*s" };
    const char* formlf                  = "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
    int         i, ft, ftype, nn, nrfp, nrfpA;
    char        f1[STRLEN];
    char        alc[MAXATOMLIST + 1][20];
    /* One force parameter more, so we can check if we read too many */
    double c[MAXFORCEPARAM + 1];

    if ((bondAtomType && at) || (!bondAtomType && !at))
    {
        gmx_incons("You should pass either bondAtomType or at to push_bt");
    }

    /* Make format string (nral ints+functype) */
    if ((nn = sscanf(line, formal[nral], alc[0], alc[1], alc[2], alc[3], alc[4], alc[5])) != nral + 1)
    {
        auto message = gmx::formatString("Not enough atomtypes (%d instead of %d)", nn - 1, nral);
        wi->addError(message);
        return;
    }

    ft    = strtol(alc[nral], nullptr, 10);
    ftype = ifunc_index(d, ft);
    nrfp  = NRFP(ftype);
    nrfpA = interaction_function[ftype].nrfpA;
    strcpy(f1, formnl[nral]);
    strcat(f1, formlf);
    if ((nn = sscanf(
                 line, f1, &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8], &c[9], &c[10], &c[11], &c[12]))
        != nrfp)
    {
        if (nn == nrfpA)
        {
            /* Copy the B-state from the A-state */
            copy_B_from_A(ftype, c);
        }
        else
        {
            if (nn < nrfpA)
            {
                wi->addError("Not enough parameters");
            }
            else if (nn > nrfpA && nn < nrfp)
            {
                wi->addError("Too many parameters or not enough parameters for topology B");
            }
            else if (nn > nrfp)
            {
                wi->addError("Too many parameters");
            }
            for (i = nn; (i < nrfp); i++)
            {
                c[i] = 0.0;
            }
        }
    }
    std::vector<int> atomTypes =
            atomTypesFromAtomNames(at, bondAtomType, gmx::arrayRefFromArray(alc, nral), wi);
    std::array<real, MAXFORCEPARAM> forceParam;
    for (int i = 0; (i < nrfp); i++)
    {
        forceParam[i] = c[i];
    }
    push_bondtype(&(bt[ftype]), InteractionOfType(atomTypes, forceParam), nral, ftype, FALSE, line, wi);
}


void push_dihedraltype(Directive                         d,
                       gmx::ArrayRef<InteractionsOfType> bt,
                       PreprocessingBondAtomType*        bondAtomType,
                       char*                             line,
                       WarningHandler*                   wi)
{
    const char* formal[MAXATOMLIST + 1] = {
        "%s", "%s%s", "%s%s%s", "%s%s%s%s", "%s%s%s%s%s", "%s%s%s%s%s%s", "%s%s%s%s%s%s%s"
    };
    const char* formnl[MAXATOMLIST + 1] = { "%*s",
                                            "%*s%*s",
                                            "%*s%*s%*s",
                                            "%*s%*s%*s%*s",
                                            "%*s%*s%*s%*s%*s",
                                            "%*s%*s%*s%*s%*s%*s",
                                            "%*s%*s%*s%*s%*s%*s%*s" };
    const char* formlf[MAXFORCEPARAM]   = {
        "%lf",
        "%lf%lf",
        "%lf%lf%lf",
        "%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
    };
    int    i, ft, ftype, nn, nrfp, nrfpA, nral;
    char   f1[STRLEN];
    char   alc[MAXATOMLIST + 1][20];
    double c[MAXFORCEPARAM];
    bool   bAllowRepeat;

    /* This routine accepts dihedraltypes defined from either 2 or 4 atoms.
     *
     * We first check for 2 atoms with the 3th column being an integer
     * defining the type. If this isn't the case, we try it with 4 atoms
     * and the 5th column defining the dihedral type.
     */
    nn = sscanf(line, formal[4], alc[0], alc[1], alc[2], alc[3], alc[4]);
    if (nn >= 3 && strlen(alc[2]) == 1 && isdigit(alc[2][0]))
    {
        nral = 2;
        ft   = strtol(alc[nral], nullptr, 10);
        /* Move atom types around a bit and use 'X' for wildcard atoms
         * to create a 4-atom dihedral definition with arbitrary atoms in
         * position 1 and 4.
         */
        if (alc[2][0] == '2')
        {
            /* improper - the two atomtypes are 1,4. Use wildcards for 2,3 */
            strcpy(alc[3], alc[1]);
            sprintf(alc[2], "X");
            sprintf(alc[1], "X");
            /* alc[0] stays put */
        }
        else
        {
            /* proper - the two atomtypes are 2,3. Use wildcards for 1,4 */
            sprintf(alc[3], "X");
            strcpy(alc[2], alc[1]);
            strcpy(alc[1], alc[0]);
            sprintf(alc[0], "X");
        }
    }
    else if (nn == 5 && strlen(alc[4]) == 1 && isdigit(alc[4][0]))
    {
        nral = 4;
        ft   = strtol(alc[nral], nullptr, 10);
    }
    else
    {
        auto message = gmx::formatString(
                "Incorrect number of atomtypes for dihedral (%d instead of 2 or 4)", nn - 1);
        wi->addError(message);
        return;
    }

    if (ft == 9)
    {
        /* Previously, we have always overwritten parameters if e.g. a torsion
           with the same atomtypes occurs on multiple lines. However, CHARMM and
           some other force fields specify multiple dihedrals over some bonds,
           including cosines with multiplicity 6 and somethimes even higher.
           Thus, they cannot be represented with Ryckaert-Bellemans terms.
           To add support for these force fields, Dihedral type 9 is identical to
           normal proper dihedrals, but repeated entries are allowed.
         */
        bAllowRepeat = TRUE;
        ft           = 1;
    }
    else
    {
        bAllowRepeat = FALSE;
    }


    ftype = ifunc_index(d, ft);
    nrfp  = NRFP(ftype);
    nrfpA = interaction_function[ftype].nrfpA;

    strcpy(f1, formnl[nral]);
    strcat(f1, formlf[nrfp - 1]);

    /* Check number of parameters given */
    if ((nn = sscanf(
                 line, f1, &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8], &c[9], &c[10], &c[11]))
        != nrfp)
    {
        if (nn == nrfpA)
        {
            /* Copy the B-state from the A-state */
            copy_B_from_A(ftype, c);
        }
        else
        {
            if (nn < nrfpA)
            {
                wi->addError("Not enough parameters");
            }
            else if (nn > nrfpA && nn < nrfp)
            {
                wi->addError("Too many parameters or not enough parameters for topology B");
            }
            else if (nn > nrfp)
            {
                wi->addError("Too many parameters");
            }
            for (i = nn; (i < nrfp); i++)
            {
                c[i] = 0.0;
            }
        }
    }

    std::vector<int>                atoms;
    std::array<real, MAXFORCEPARAM> forceParam;
    for (int i = 0; (i < 4); i++)
    {
        if (!strcmp(alc[i], "X"))
        {
            atoms.emplace_back(-1);
        }
        else
        {
            auto atomNumber = bondAtomType->bondAtomTypeFromName(alc[i]);
            if (!atomNumber.has_value())
            {
                auto message = gmx::formatString("Unknown bond_atomtype %s", alc[i]);
                warning_error_and_exit(wi, message, FARGS);
            }
            atoms.emplace_back(*atomNumber);
        }
    }
    for (int i = 0; (i < nrfp); i++)
    {
        forceParam[i] = c[i];
    }
    /* Always use 4 atoms here, since we created two wildcard atoms
     * if there wasn't of them 4 already.
     */
    push_bondtype(&(bt[ftype]), InteractionOfType(atoms, forceParam), 4, ftype, bAllowRepeat, line, wi);
}


void push_nbt(Directive d, t_nbparam** nbt, PreprocessingAtomTypes* atypes, char* pline, int nb_funct, WarningHandler* wi)
{
    /* swap the atoms */
    const char* form3 = "%*s%*s%*s%lf%lf%lf";
    const char* form4 = "%*s%*s%*s%lf%lf%lf%lf";
    const char* form5 = "%*s%*s%*s%lf%lf%lf%lf%lf";
    char        a0[80], a1[80];
    int         i, f, n, ftype, nrfp;
    double      c[4], dum;
    real        cr[4];
    t_nbparam*  nbp;
    bool        bId;

    if (sscanf(pline, "%s%s%d", a0, a1, &f) != 3)
    {
        too_few(wi);
        return;
    }

    ftype = ifunc_index(d, f);

    if (ftype != nb_funct)
    {
        auto message = gmx::formatString("Trying to add %s while the default nonbond type is %s",
                                         interaction_function[ftype].longname,
                                         interaction_function[nb_funct].longname);
        wi->addError(message);
        return;
    }

    /* Get the force parameters */
    nrfp = NRFP(ftype);
    if (ftype == F_LJ14)
    {
        n = sscanf(pline, form4, &c[0], &c[1], &c[2], &c[3]);
        if (n < 2)
        {
            too_few(wi);
            return;
        }
        /* When the B topology parameters are not set,
         * copy them from topology A
         */
        GMX_ASSERT(nrfp <= NRFP(F_LJ14), "LJ-14 cannot have more than 4 parameters");
        for (i = n; i < nrfp; i++)
        {
            c[i] = c[i - 2];
        }
    }
    else if (ftype == F_LJC14_Q)
    {
        n = sscanf(pline, form5, &c[0], &c[1], &c[2], &c[3], &dum);
        if (n != 4)
        {
            incorrect_n_param(wi);
            return;
        }
    }
    else if (nrfp == 2)
    {
        if (sscanf(pline, form3, &c[0], &c[1], &dum) != 2)
        {
            incorrect_n_param(wi);
            return;
        }
    }
    else if (nrfp == 3)
    {
        if (sscanf(pline, form4, &c[0], &c[1], &c[2], &dum) != 3)
        {
            incorrect_n_param(wi);
            return;
        }
    }
    else
    {
        auto message =
                gmx::formatString("Number of force parameters for nonbonded interactions is %d", nrfp);
        warning_error_and_exit(wi, message, FARGS);
    }
    for (i = 0; (i < nrfp); i++)
    {
        cr[i] = c[i];
    }

    /* Put the parameters in the matrix */
    auto ai = atypes->atomTypeFromName(a0);
    if (!ai.has_value())
    {
        auto message = gmx::formatString("Atomtype %s not found", a0);
        warning_error_and_exit(wi, message, FARGS);
    }
    auto aj = atypes->atomTypeFromName(a1);
    if (!aj.has_value())
    {
        auto message = gmx::formatString("Atomtype %s not found", a1);
        warning_error_and_exit(wi, message, FARGS);
    }
    nbp = &(nbt[std::max(*ai, *aj)][std::min(*ai, *aj)]);

    if (nbp->bSet)
    {
        bId = TRUE;
        for (i = 0; i < nrfp; i++)
        {
            bId = bId && (nbp->c[i] == cr[i]);
        }
        if (!bId)
        {
            auto message = gmx::formatString(
                    "Non-bonded parameters were defined previously (e.g. in the forcefield files), "
                    "and have now been defined again. This could happen e.g. if you would "
                    "use a self-contained molecule .itp file that duplicates or replaces "
                    "the contents of the standard force-field files. You should check "
                    "the contents of your files and remove such repetition. If you know "
                    "you should override the previous definitions, then you could choose "
                    "to suppress this warning with -maxwarn.");
            wi->addWarning(message);
            fprintf(stderr, "  old:");
            for (i = 0; i < nrfp; i++)
            {
                fprintf(stderr, " %g", nbp->c[i]);
            }
            fprintf(stderr, " new\n%s\n", pline);
        }
    }
    nbp->bSet = TRUE;
    for (i = 0; i < nrfp; i++)
    {
        nbp->c[i] = cr[i];
    }
}

void push_cmaptype(Directive                         d,
                   gmx::ArrayRef<InteractionsOfType> bt,
                   int                               nral,
                   PreprocessingAtomTypes*           atomtypes,
                   PreprocessingBondAtomType*        bondAtomType,
                   char*                             line,
                   WarningHandler*                   wi)
{
    GMX_ASSERT(nral == NRAL(F_CMAP), "CMAP requires 5 atoms per interaction");

    const char* formal = "%s%s%s%s%s%s%s%s%n";

    int  ft, ftype, nn, nrfp, nrfpA, nrfpB;
    int  start, nchar_consumed;
    int  nxcmap, nycmap, ncmap, read_cmap, sl, nct;
    char s[20], alc[MAXATOMLIST + 2][20];

    /* Keep the compiler happy */
    read_cmap = 0;
    start     = 0;

    /* Here we can only check for < 8 */
    if ((nn = sscanf(line, formal, alc[0], alc[1], alc[2], alc[3], alc[4], alc[nral], alc[nral + 1], alc[nral + 2], &nchar_consumed))
        < nral + 3)
    {
        auto message = gmx::formatString(
                "Incorrect number of atomtypes for cmap type (%d instead of %d)", nn - 3, nral);
        wi->addError(message);
        return;
    }
    start += nchar_consumed;

    ft = strtol(alc[nral], nullptr, 10);
    GMX_RELEASE_ASSERT(ft == 1, "Invalid function type for cmap type: must be 1");
    nxcmap = strtol(alc[nral + 1], nullptr, 10);
    nycmap = strtol(alc[nral + 2], nullptr, 10);

    /* Check for equal grid spacing in x and y dims */
    if (nxcmap != nycmap)
    {
        auto message = gmx::formatString(
                "Not the same grid spacing in x and y for cmap grid: x=%d, y=%d", nxcmap, nycmap);
        wi->addError(message);
    }

    ncmap = nxcmap * nycmap;
    ftype = ifunc_index(d, ft);
    nrfpA = strtol(alc[nral + 1], nullptr, 10) * strtol(alc[nral + 1], nullptr, 10);
    nrfpB = strtol(alc[nral + 2], nullptr, 10) * strtol(alc[nral + 2], nullptr, 10);
    nrfp  = nrfpA + nrfpB;

    /* Read in CMAP parameters */
    sl = 0;
    for (int i = 0; i < ncmap; i++)
    {
        while (isspace(*(line + start + sl)))
        {
            sl++;
        }
        nn = sscanf(line + start + sl, " %s ", s);
        sl += strlen(s);
        bt[F_CMAP].cmap.emplace_back(strtod(s, nullptr));

        if (nn == 1)
        {
            read_cmap++;
        }
        else
        {
            auto message = gmx::formatString(
                    "Error in reading cmap parameter for atomtypes %s %s %s %s %s: found %d, "
                    "expected %d",
                    alc[0],
                    alc[1],
                    alc[2],
                    alc[3],
                    alc[4],
                    read_cmap,
                    ncmap);
            wi->addError(message);
        }
    }
    if ((nn = sscanf(line + start + sl, " %s ", s)))
    {
        if (nn == 1)
        {
            auto message = gmx::formatString(
                    "One or more unread cmap parameters exist for atomtypes %s %s %s %s %s",
                    alc[0],
                    alc[1],
                    alc[2],
                    alc[3],
                    alc[4]);
            wi->addError(message);
        }
    }

    /* Check do that we got the number of parameters we expected */
    if (read_cmap == nrfpA)
    {
        for (int i = 0; i < ncmap; i++)
        {
            bt[F_CMAP].cmap.emplace_back(bt[F_CMAP].cmap[i]);
        }
    }
    else
    {
        if (read_cmap < nrfpA)
        {
            wi->addError("Not enough cmap parameters");
        }
        else if (read_cmap > nrfpA && read_cmap < nrfp)
        {
            wi->addError("Too many cmap parameters or not enough parameters for topology B");
        }
        else if (read_cmap > nrfp)
        {
            wi->addError("Too many cmap parameters");
        }
    }


    /* Set grid spacing and the number of grids (we assume these numbers to be the same for all
     * grids so we can safely assign them each time
     */
    bt[F_CMAP].cmapGridSpacing_ = nxcmap; /* Or nycmap, they need to be equal */

    for (int i = 0; (i < nral); i++)
    {
        /* Assign a grid number to each cmap_type */
        GMX_RELEASE_ASSERT(bondAtomType != nullptr, "Need valid PreprocessingBondAtomType object");
        auto cmapBondAtomType = bondAtomType->bondAtomTypeFromName(alc[i]);
        if (!cmapBondAtomType)
        {
            auto message = gmx::formatString(
                    "Unknown bond_atomtype for %s in cmap atomtypes %s %s %s %s %s",
                    alc[i],
                    alc[0],
                    alc[1],
                    alc[2],
                    alc[3],
                    alc[4]);
            wi->addError(message);
            continue;
        }
        bt[F_CMAP].cmapAtomTypes.emplace_back(*cmapBondAtomType);
    }

    /* Assign a type number to this cmap */
    bt[F_CMAP].cmapAtomTypes.emplace_back(bt[F_CMAP].numCmaps_);
    bt[F_CMAP].numCmaps_++;

    /* Check for the correct number of atoms (again) */
    nct = (nral + 1) * bt[F_CMAP].numCmaps_;
    if (bt[F_CMAP].nct() != static_cast<std::size_t>(nct))
    {
        auto message = gmx::formatString(
                "Incorrect number of atomtypes (%d) in cmap type %d\n", nct, bt[F_CMAP].numCmaps_);
        wi->addError(message);
    }
    std::vector<int> atomTypes =
            atomTypesFromAtomNames(atomtypes, bondAtomType, gmx::constArrayRefFromArray(alc, nral), wi);
    std::array<real, MAXFORCEPARAM> forceParam = { NOTSET };

    /* Push the bond to the bondlist */
    push_bondtype(&(bt[ftype]), InteractionOfType(atomTypes, forceParam), nral, ftype, FALSE, line, wi);
}


static void push_atom_now(t_symtab*       symtab,
                          t_atoms*        at,
                          int             atomnr,
                          int             atomicnumber,
                          int             type,
                          char*           ctype,
                          ParticleType    ptype,
                          char*           resnumberic,
                          char*           resname,
                          char*           name,
                          real            m0,
                          real            q0,
                          int             typeB,
                          char*           ctypeB,
                          real            mB,
                          real            qB,
                          WarningHandler* wi)
{
    int           j, resind = 0, resnr;
    unsigned char ric;
    int           nr = at->nr;

    if (((nr == 0) && (atomnr != 1)) || (nr && (atomnr != at->nr + 1)))
    {
        auto message = gmx::formatString(
                "Atoms in the .top are not numbered consecutively from 1 (rather, "
                "atomnr = %d, while at->nr = %d)",
                atomnr,
                at->nr);
        warning_error_and_exit(wi, message, FARGS);
    }

    j = strlen(resnumberic) - 1;
    if (isdigit(resnumberic[j]))
    {
        ric = ' ';
    }
    else
    {
        ric = resnumberic[j];
        if (j == 0 || !isdigit(resnumberic[j - 1]))
        {
            auto message =
                    gmx::formatString("Invalid residue number '%s' for atom %d", resnumberic, atomnr);
            warning_error_and_exit(wi, message, FARGS);
        }
    }
    resnr = strtol(resnumberic, nullptr, 10);

    if (nr > 0)
    {
        resind = at->atom[nr - 1].resind;
    }
    if (nr == 0 || strcmp(resname, *at->resinfo[resind].name) != 0
        || resnr != at->resinfo[resind].nr || ric != at->resinfo[resind].ic)
    {
        if (nr == 0)
        {
            resind = 0;
        }
        else
        {
            resind++;
        }
        at->nres = resind + 1;
        srenew(at->resinfo, at->nres);
        at->resinfo[resind].name = put_symtab(symtab, resname);
        at->resinfo[resind].nr   = resnr;
        at->resinfo[resind].ic   = ric;
    }
    else
    {
        resind = at->atom[at->nr - 1].resind;
    }

    /* New atom instance
     * get new space for arrays
     */
    srenew(at->atom, nr + 1);
    srenew(at->atomname, nr + 1);
    srenew(at->atomtype, nr + 1);
    srenew(at->atomtypeB, nr + 1);

    /* fill the list */
    at->atom[nr].type  = type;
    at->atom[nr].ptype = ptype;
    at->atom[nr].q     = q0;
    at->atom[nr].m     = m0;
    at->atom[nr].typeB = typeB;
    at->atom[nr].qB    = qB;
    at->atom[nr].mB    = mB;

    at->atom[nr].resind     = resind;
    at->atom[nr].atomnumber = atomicnumber;
    at->atomname[nr]        = put_symtab(symtab, name);
    at->atomtype[nr]        = put_symtab(symtab, ctype);
    at->atomtypeB[nr]       = put_symtab(symtab, ctypeB);
    at->nr++;
}

void push_atom(t_symtab* symtab, t_atoms* at, PreprocessingAtomTypes* atypes, char* line, WarningHandler* wi)
{
    int  cgnumber, atomnr, nscan;
    char id[STRLEN], ctype[STRLEN], ctypeB[STRLEN], resnumberic[STRLEN], resname[STRLEN],
            name[STRLEN], check[STRLEN];
    double m, q, mb, qb;
    real   m0, q0, mB, qB;

    /* Fixed parameters */
    if (sscanf(line, "%s%s%s%s%s%d", id, ctype, resnumberic, resname, name, &cgnumber) != 6)
    {
        too_few(wi);
        return;
    }
    sscanf(id, "%d", &atomnr);
    auto type = atypes->atomTypeFromName(ctype);
    if (!type.has_value())
    {
        auto message = gmx::formatString("Atomtype %s not found", ctype);
        warning_error_and_exit(wi, message, FARGS);
    }
    ParticleType ptype = *atypes->atomParticleTypeFromAtomType(*type);

    /* Set default from type */
    q0         = *atypes->atomChargeFromAtomType(*type);
    m0         = *atypes->atomMassFromAtomType(*type);
    auto typeB = type;
    qB         = q0;
    mB         = m0;

    /* Optional parameters */
    nscan = sscanf(line, "%*s%*s%*s%*s%*s%*s%lf%lf%s%lf%lf%s", &q, &m, ctypeB, &qb, &mb, check);

    /* Nasty switch that falls thru all the way down! */
    if (nscan > 0)
    {
        q0 = qB = q;
        if (nscan > 1)
        {
            m0 = mB = m;
            if (nscan > 2)
            {
                typeB = atypes->atomTypeFromName(ctypeB);
                if (!typeB.has_value())
                {
                    auto message = gmx::formatString("Atomtype %s not found", ctypeB);
                    warning_error_and_exit(wi, message, FARGS);
                }
                qB = *atypes->atomChargeFromAtomType(*typeB);
                mB = *atypes->atomMassFromAtomType(*typeB);
                if (nscan > 3)
                {
                    qB = qb;
                    if (nscan > 4)
                    {
                        mB = mb;
                        if (nscan > 5)
                        {
                            wi->addError("Too many parameters");
                        }
                    }
                }
            }
        }
    }

    push_atom_now(symtab,
                  at,
                  atomnr,
                  *atypes->atomNumberFromAtomType(*type),
                  *type,
                  ctype,
                  ptype,
                  resnumberic,
                  resname,
                  name,
                  m0,
                  q0,
                  *typeB,
                  typeB == type ? ctype : ctypeB,
                  mB,
                  qB,
                  wi);
}

void push_molt(t_symtab* symtab, std::vector<MoleculeInformation>* mol, char* line, WarningHandler* wi)
{
    char type[STRLEN];
    int  nrexcl;

    if ((sscanf(line, "%s%d", type, &nrexcl)) != 2)
    {
        wi->addError("Expected a molecule type name and nrexcl");
    }

    /* Test if this moleculetype overwrites another */
    const auto found = std::find_if(
            mol->begin(), mol->end(), [&type](const auto& m) { return strcmp(*(m.name), type) == 0; });
    if (found != mol->end())
    {
        auto message = gmx::formatString("moleculetype %s is redefined", type);
        warning_error_and_exit(wi, message, FARGS);
    }

    mol->emplace_back();
    mol->back().initMolInfo();

    /* Fill in the values */
    mol->back().name   = put_symtab(symtab, type);
    mol->back().nrexcl = nrexcl;
}

static bool findIfAllNBAtomsMatch(gmx::ArrayRef<const int> atomsFromParameterArray,
                                  gmx::ArrayRef<const int> atomsFromCurrentParameter,
                                  const t_atoms*           at,
                                  bool                     bB)
{
    if (atomsFromParameterArray.size() != atomsFromCurrentParameter.size())
    {
        return false;
    }
    else if (bB)
    {
        for (gmx::Index i = 0; i < atomsFromCurrentParameter.ssize(); i++)
        {
            if (at->atom[atomsFromCurrentParameter[i]].typeB != atomsFromParameterArray[i])
            {
                return false;
            }
        }
        return true;
    }
    else
    {
        for (gmx::Index i = 0; i < atomsFromCurrentParameter.ssize(); i++)
        {
            if (at->atom[atomsFromCurrentParameter[i]].type != atomsFromParameterArray[i])
            {
                return false;
            }
        }
        return true;
    }
}

static bool default_nb_params(int                               ftype,
                              gmx::ArrayRef<InteractionsOfType> bt,
                              t_atoms*                          at,
                              InteractionOfType*                p,
                              int                               c_start,
                              bool                              bB,
                              bool                              bGenPairs)
{
    int                ti, tj, ntype;
    bool               bFound;
    InteractionOfType* pi    = nullptr;
    int                nr    = bt[ftype].size();
    int                nral  = NRAL(ftype);
    int                nrfpA = interaction_function[ftype].nrfpA;
    int                nrfpB = interaction_function[ftype].nrfpB;

    if ((!bB && nrfpA == 0) || (bB && nrfpB == 0))
    {
        return TRUE;
    }

    bFound = FALSE;
    if (bGenPairs)
    {
        /* First test the generated-pair position to save
         * time when we have 1000*1000 entries for e.g. OPLS...
         */
        ntype = static_cast<int>(std::sqrt(static_cast<double>(nr)));
        GMX_ASSERT(ntype * ntype == nr,
                   "Number of pairs of generated non-bonded parameters should be a perfect square");
        if (bB)
        {
            ti = at->atom[p->ai()].typeB;
            tj = at->atom[p->aj()].typeB;
        }
        else
        {
            ti = at->atom[p->ai()].type;
            tj = at->atom[p->aj()].type;
        }
        pi = &(bt[ftype].interactionTypes[ntype * ti + tj]);
        if (pi->atoms().ssize() < nral)
        {
            /* not initialized yet with atom names */
            bFound = false;
        }
        else
        {
            bFound = ((ti == pi->ai()) && (tj == pi->aj()));
        }
    }

    gmx::ArrayRef<const int> paramAtoms = p->atoms();
    /* Search explicitly if we didnt find it */
    if (!bFound)
    {
        auto foundParameter =
                std::find_if(bt[ftype].interactionTypes.begin(),
                             bt[ftype].interactionTypes.end(),
                             [&paramAtoms, &at, &bB](const auto& param) {
                                 return findIfAllNBAtomsMatch(param.atoms(), paramAtoms, at, bB);
                             });
        if (foundParameter != bt[ftype].interactionTypes.end())
        {
            bFound = true;
            pi     = &(*foundParameter);
        }
    }

    if (bFound)
    {
        gmx::ArrayRef<const real> forceParam = pi->forceParam();
        if (bB)
        {
            if (nrfpA + nrfpB > MAXFORCEPARAM)
            {
                gmx_incons("Too many force parameters");
            }
            for (int j = c_start; j < nrfpB; j++)
            {
                p->setForceParameter(nrfpA + j, forceParam[j]);
            }
        }
        else
        {
            for (int j = c_start; j < nrfpA; j++)
            {
                p->setForceParameter(j, forceParam[j]);
            }
        }
    }
    else
    {
        for (int j = c_start; j < nrfpA; j++)
        {
            p->setForceParameter(j, 0.0);
        }
    }
    return bFound;
}

static bool default_cmap_params(gmx::ArrayRef<InteractionsOfType> bondtype,
                                t_atoms*                          at,
                                PreprocessingAtomTypes*           atypes,
                                InteractionOfType*                p,
                                bool                              bB,
                                int*                              cmap_type,
                                int*                              nparam_def,
                                WarningHandler*                   wi)
{
    int  nparam_found;
    int  ct;
    bool bFound = false;

    nparam_found = 0;
    ct           = 0;

    /* Match the current cmap angle against the list of cmap_types */
    for (std::size_t i = 0; i < bondtype[F_CMAP].nct() && !bFound; i += NRAL(F_CMAP) + 1)
    {
        if (bB) {}
        else
        {
            if ((atypes->bondAtomTypeFromAtomType(at->atom[p->ai()].type)
                 == bondtype[F_CMAP].cmapAtomTypes[i])
                && (atypes->bondAtomTypeFromAtomType(at->atom[p->aj()].type)
                    == bondtype[F_CMAP].cmapAtomTypes[i + 1])
                && (atypes->bondAtomTypeFromAtomType(at->atom[p->ak()].type)
                    == bondtype[F_CMAP].cmapAtomTypes[i + 2])
                && (atypes->bondAtomTypeFromAtomType(at->atom[p->al()].type)
                    == bondtype[F_CMAP].cmapAtomTypes[i + 3])
                && (atypes->bondAtomTypeFromAtomType(at->atom[p->am()].type)
                    == bondtype[F_CMAP].cmapAtomTypes[i + 4]))
            {
                /* Found cmap torsion */
                bFound       = true;
                ct           = bondtype[F_CMAP].cmapAtomTypes[i + NRAL(F_CMAP)];
                nparam_found = 1;
            }
        }
    }

    /* If we did not find a matching type for this cmap torsion */
    if (!bFound)
    {
        auto message = gmx::formatString("Unknown cmap torsion between atoms %d %d %d %d %d",
                                         p->ai() + 1,
                                         p->aj() + 1,
                                         p->ak() + 1,
                                         p->al() + 1,
                                         p->am() + 1);
        warning_error_and_exit(wi, message, FARGS);
    }

    *nparam_def = nparam_found;
    *cmap_type  = ct;

    return bFound;
}

/* Returns the number of exact atom type matches, i.e. non wild-card matches,
 * returns -1 when there are no matches at all.
 */
static int findNumberOfDihedralAtomMatches(const InteractionOfType&       bondType,
                                           const gmx::ArrayRef<const int> atomTypes)
{
    GMX_RELEASE_ASSERT(atomTypes.size() == 4, "Dihedrals have 4 atom types");
    const gmx::ArrayRef<const int> bondTypeAtomTypes = bondType.atoms();
    GMX_RELEASE_ASSERT(bondTypeAtomTypes.size() == 4, "Dihedral types have 4 atom types");
    int numExactMatches = 0;
    if (std::equal(bondTypeAtomTypes.begin(),
                   bondTypeAtomTypes.end(),
                   atomTypes.begin(),
                   [&numExactMatches](int bondTypeAtomType, int atomType) {
                       if (bondTypeAtomType == atomType)
                       {
                           // Found an exact atom type match
                           ++numExactMatches;
                           return true;
                       }
                       else if (bondTypeAtomType == -1)
                       {
                           // Found a wildcard atom type match
                           return true;
                       }
                       // Atom types do not match
                       return false;
                   }))
    {
        return numExactMatches;
    }
    return -1;
}

static std::vector<InteractionOfType>::iterator
defaultInteractionsOfType(int                               ftype,
                          gmx::ArrayRef<InteractionsOfType> bondType,
                          const gmx::ArrayRef<const int>    atomTypes,
                          int*                              nparam_def)
{
    int nparam_found = 0;

    if (ftype == F_PDIHS || ftype == F_RBDIHS || ftype == F_IDIHS || ftype == F_PIDIHS)
    {
        int nmatch_max = -1;

        /* For dihedrals we allow wildcards. We choose the first type
         * that has the most exact matches, i.e. non-wildcard matches.
         */
        auto prevPos = bondType[ftype].interactionTypes.end();
        auto pos     = bondType[ftype].interactionTypes.begin();
        while (pos != bondType[ftype].interactionTypes.end() && nmatch_max < 4)
        {
            pos = std::find_if(bondType[ftype].interactionTypes.begin(),
                               bondType[ftype].interactionTypes.end(),
                               [&atomTypes, &nmatch_max](const auto& elemBondType) {
                                   return (findNumberOfDihedralAtomMatches(elemBondType, atomTypes)
                                           > nmatch_max);
                               });
            if (pos != bondType[ftype].interactionTypes.end())
            {
                prevPos    = pos;
                nmatch_max = findNumberOfDihedralAtomMatches(*pos, atomTypes);
            }
        }

        if (prevPos != bondType[ftype].interactionTypes.end())
        {
            nparam_found++;

            /* Find additional matches for this dihedral - necessary
             * for ftype==9.
             * The rule in that case is that additional matches
             * HAVE to be on adjacent lines!
             */
            bool bSame = true;
            // Advance iterator (like std::advance) without incrementing past end (UB)
            const auto safeAdvance = [](auto& it, auto n, auto end) {
                it = end - it > n ? it + n : end;
            };
            /* Continue from current iterator position */
            auto       nextPos = prevPos;
            const auto endIter = bondType[ftype].interactionTypes.end();
            safeAdvance(nextPos, 2, endIter);
            for (; nextPos < endIter && bSame; safeAdvance(nextPos, 2, endIter))
            {
                bSame = (prevPos->ai() == nextPos->ai() && prevPos->aj() == nextPos->aj()
                         && prevPos->ak() == nextPos->ak() && prevPos->al() == nextPos->al());
                if (bSame)
                {
                    nparam_found++;
                }
                /* nparam_found will be increased as long as the numbers match */
            }
        }
        *nparam_def = nparam_found;
        return prevPos;
    }
    else /* Not a dihedral */
    {
        auto found = std::find_if(
                bondType[ftype].interactionTypes.begin(),
                bondType[ftype].interactionTypes.end(),
                [&atomTypes](const auto& param) {
                    return std::equal(param.atoms().begin(), param.atoms().end(), atomTypes.begin());
                });
        if (found != bondType[ftype].interactionTypes.end())
        {
            nparam_found = 1;
        }
        *nparam_def = nparam_found;
        return found;
    }
}


void push_bond(Directive                         d,
               gmx::ArrayRef<InteractionsOfType> bondtype,
               gmx::ArrayRef<InteractionsOfType> bond,
               t_atoms*                          at,
               PreprocessingAtomTypes*           atypes,
               char*                             line,
               bool                              bBonded,
               bool                              bGenPairs,
               real                              fudgeQQ,
               bool                              bZero,
               bool*                             bWarn_copy_A_B,
               WarningHandler*                   wi)
{
    const char* aaformat[MAXATOMLIST] = { "%d%d",       "%d%d%d",       "%d%d%d%d",
                                          "%d%d%d%d%d", "%d%d%d%d%d%d", "%d%d%d%d%d%d%d" };
    const char* asformat[MAXATOMLIST] = {
        "%*s%*s",          "%*s%*s%*s",          "%*s%*s%*s%*s",
        "%*s%*s%*s%*s%*s", "%*s%*s%*s%*s%*s%*s", "%*s%*s%*s%*s%*s%*s%*s"
    };
    const char* ccformat = "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
    int         nral, nral_fmt, nread, ftype;
    char        format[STRLEN];
    /* One force parameter more, so we can check if we read too many */
    double                           cc[MAXFORCEPARAM + 1];
    std::array<int, MAXATOMLIST + 1> aa;
    bool                             bFoundA = FALSE, bFoundB = FALSE, bDef, bSwapParity = FALSE;
    int                              nparam_defA, nparam_defB;

    nparam_defA = nparam_defB = 0;

    ftype = ifunc_index(d, 1);
    nral  = NRAL(ftype);
    for (int j = 0; j < nral; j++)
    {
        aa[j] = NOTSET;
    }
    bDef = (NRFP(ftype) > 0);

    if (ftype == F_SETTLE)
    {
        /* SETTLE acts on 3 atoms, but the topology format only specifies
         * the first atom (for historical reasons).
         */
        nral_fmt = 1;
    }
    else
    {
        nral_fmt = nral;
    }

    nread = sscanf(line, aaformat[nral_fmt - 1], &aa[0], &aa[1], &aa[2], &aa[3], &aa[4], &aa[5]);

    if (ftype == F_SETTLE)
    {
        aa[3] = aa[1];
        aa[1] = aa[0] + 1;
        aa[2] = aa[0] + 2;
    }

    if (nread < nral_fmt)
    {
        too_few(wi);
        return;
    }
    else if (nread > nral_fmt)
    {
        /* this is a hack to allow for virtual sites with swapped parity */
        bSwapParity = (aa[nral] < 0);
        if (bSwapParity)
        {
            aa[nral] = -aa[nral];
        }
        ftype = ifunc_index(d, aa[nral]);
        if (bSwapParity)
        {
            switch (ftype)
            {
                case F_VSITE3FAD:
                case F_VSITE3OUT: break;
                default:
                    auto message =
                            gmx::formatString("Negative function types only allowed for %s and %s",
                                              interaction_function[F_VSITE3FAD].longname,
                                              interaction_function[F_VSITE3OUT].longname);
                    warning_error_and_exit(wi, message, FARGS);
            }
        }
    }


    /* Check for double atoms and atoms out of bounds, then convert to 0-based indexing */
    for (int i = 0; (i < nral); i++)
    {
        if (aa[i] < 1 || aa[i] > at->nr)
        {
            auto message = gmx::formatString(
                    "Atom index (%d) in %s out of bounds (1-%d).\n"
                    "This probably means that you have inserted topology section \"%s\"\n"
                    "in a part belonging to a different molecule than you intended to.\n"
                    "In that case move the \"%s\" section to the right molecule.",
                    aa[i],
                    enumValueToString(d),
                    at->nr,
                    enumValueToString(d),
                    enumValueToString(d));
            warning_error_and_exit(wi, message, FARGS);
        }
        for (int j = i + 1; (j < nral); j++)
        {
            if (aa[i] == aa[j])
            {
                auto message = gmx::formatString(
                        "Duplicate atom index (%d) in %s", aa[i], enumValueToString(d));
                if (ftype == F_ANGRES)
                {
                    /* Since the angle restraints uses 2 pairs of atoms to
                     * defines an angle between vectors, it can be useful
                     * to use one atom twice, so we only issue a note here.
                     */
                    wi->addNote(message);
                }
                else
                {
                    wi->addError(message);
                }
            }
        }

        // Convert to 0-based indexing
        --aa[i];
    }

    // These are the atom indices for this interaction
    auto atomIndices = gmx::ArrayRef<const int>(aa).subArray(0, nral);

    // Look up the A-state atom types for this interaction
    std::vector<int> atomTypes(atomIndices.size());
    std::transform(atomIndices.begin(), atomIndices.end(), atomTypes.begin(), [at, atypes](const int atomIndex) {
        return atypes->bondAtomTypeFromAtomType(at->atom[atomIndex].type).value();
    });
    // Look up the B-state atom types for this interaction
    std::vector<int> atomTypesB(atomIndices.size());
    std::transform(atomIndices.begin(), atomIndices.end(), atomTypesB.begin(), [at, atypes](const int atomIndex) {
        return atypes->bondAtomTypeFromAtomType(at->atom[atomIndex].typeB).value();
    });

    /* default force parameters  */
    /* need to have an empty but initialized param array for some reason */
    std::array<real, MAXFORCEPARAM> forceParam = { 0.0 };

    /* Get force params for normal and free energy perturbation
     * studies, as determined by types!
     */
    InteractionOfType param(atomIndices, forceParam, "");

    std::vector<InteractionOfType>::iterator foundAParameter = bondtype[ftype].interactionTypes.end();
    std::vector<InteractionOfType>::iterator foundBParameter = bondtype[ftype].interactionTypes.end();
    if (bBonded)
    {
        if (NRFPA(ftype) == 0)
        {
            bFoundA = true;
        }
        else
        {
            foundAParameter = defaultInteractionsOfType(ftype, bondtype, atomTypes, &nparam_defA);
            if (foundAParameter != bondtype[ftype].interactionTypes.end())
            {
                /* Copy the A-state and B-state default parameters. */
                GMX_ASSERT(NRFPA(ftype) + NRFPB(ftype) <= MAXFORCEPARAM,
                           "Bonded interactions may have at most 12 parameters");
                gmx::ArrayRef<const real> defaultParam = foundAParameter->forceParam();
                for (int j = 0; (j < NRFPA(ftype) + NRFPB(ftype)); j++)
                {
                    param.setForceParameter(j, defaultParam[j]);
                }
                bFoundA = true;
            }
        }

        if (NRFPB(ftype) == 0)
        {
            bFoundB = true;
        }
        else
        {
            foundBParameter = defaultInteractionsOfType(ftype, bondtype, atomTypesB, &nparam_defB);
            if (foundBParameter != bondtype[ftype].interactionTypes.end())
            {
                /* Copy only the B-state default parameters */
                gmx::ArrayRef<const real> defaultParam = foundBParameter->forceParam();
                for (int j = NRFPA(ftype); (j < NRFP(ftype)); j++)
                {
                    param.setForceParameter(j, defaultParam[j]);
                }
                bFoundB = true;
            }
        }
    }
    else if (ftype == F_LJ14)
    {
        bFoundA = default_nb_params(ftype, bondtype, at, &param, 0, FALSE, bGenPairs);
        bFoundB = default_nb_params(ftype, bondtype, at, &param, 0, TRUE, bGenPairs);
    }
    else if (ftype == F_LJC14_Q)
    {
        /* Fill in the A-state charges as default parameters */
        param.setForceParameter(0, fudgeQQ);
        param.setForceParameter(1, at->atom[param.ai()].q);
        param.setForceParameter(2, at->atom[param.aj()].q);
        /* The default LJ parameters are the standard 1-4 parameters */
        bFoundA = default_nb_params(F_LJ14, bondtype, at, &param, 3, FALSE, bGenPairs);
        bFoundB = TRUE;
    }
    else if (ftype == F_LJC_PAIRS_NB)
    {
        /* Defaults are not supported here */
        bFoundA = FALSE;
        bFoundB = TRUE;
    }
    else
    {
        gmx_incons("Unknown function type in push_bond");
    }

    if (nread > nral_fmt)
    {
        /* Manually specified parameters - in this case we discard multiple torsion info! */

        strcpy(format, asformat[nral_fmt - 1]);
        strcat(format, ccformat);

        nread = sscanf(line,
                       format,
                       &cc[0],
                       &cc[1],
                       &cc[2],
                       &cc[3],
                       &cc[4],
                       &cc[5],
                       &cc[6],
                       &cc[7],
                       &cc[8],
                       &cc[9],
                       &cc[10],
                       &cc[11],
                       &cc[12]);

        if ((nread == NRFPA(ftype)) && (NRFPB(ftype) != 0))
        {
            /* We only have to issue a warning if these atoms are perturbed! */
            bool                     bPert      = false;
            gmx::ArrayRef<const int> paramAtoms = param.atoms();
            for (int j = 0; (j < nral); j++)
            {
                bPert = bPert || PERTURBED(at->atom[paramAtoms[j]]);
            }

            if (bPert && *bWarn_copy_A_B)
            {
                auto message = gmx::formatString(
                        "Some parameters for bonded interaction involving "
                        "perturbed atoms are specified explicitly in "
                        "state A, but not B - copying A to B");
                wi->addWarning(message);
                *bWarn_copy_A_B = FALSE;
            }

            /* If only the A parameters were specified, copy them to the B state */
            /* The B-state parameters correspond to the first nrfpB
             * A-state parameters.
             */
            for (int j = 0; (j < NRFPB(ftype)); j++)
            {
                cc[nread++] = cc[j];
            }
        }

        /* If nread was 0 or EOF, no parameters were read => use defaults.
         * If nread was nrfpA we copied above so nread=nrfp.
         * If nread was nrfp we are cool.
         * For F_LJC14_Q we allow supplying fudgeQQ only.
         * Anything else is an error!
         */
        if ((nread != 0) && (nread != EOF) && (nread != NRFP(ftype)) && !(ftype == F_LJC14_Q && nread == 1))
        {
            auto message = gmx::formatString(
                    "Incorrect number of parameters - found %d, expected %d "
                    "or %d for %s (after the function type).",
                    nread,
                    NRFPA(ftype),
                    NRFP(ftype),
                    interaction_function[ftype].longname);
            warning_error_and_exit(wi, message, FARGS);
        }

        for (int j = 0; (j < nread); j++)
        {
            param.setForceParameter(j, cc[j]);
        }
        /* Check whether we have to use the defaults */
        if (nread == NRFP(ftype))
        {
            bDef = FALSE;
        }
    }
    else
    {
        nread = 0;
    }
    /* nread now holds the number of force parameters read! */

    if (bDef)
    {
        /* Use defaults */
        /* When we have multiple terms it would be very dangerous to allow perturbations to a different atom type! */
        if (ftype == F_PDIHS)
        {
            if ((nparam_defA != nparam_defB)
                || ((nparam_defA > 1 || nparam_defB > 1) && (foundAParameter != foundBParameter)))
            {
                auto message = gmx::formatString(
                        "Cannot automatically perturb a torsion with multiple terms to different "
                        "form.\n"
                        "Please specify perturbed parameters manually for this torsion in your "
                        "topology!");
                wi->addError(message);
            }
        }

        if (nread > 0 && nread < NRFPA(ftype))
        {
            /* Issue an error, do not use defaults */
            auto message = gmx::formatString(
                    "Not enough parameters, there should be at least %d (or 0 for defaults)", NRFPA(ftype));
            wi->addError(message);
        }

        if (nread == 0 || nread == EOF)
        {
            if (!bFoundA)
            {
                if (interaction_function[ftype].flags & IF_VSITE)
                {
                    for (int j = 0; j < MAXFORCEPARAM; j++)
                    {
                        param.setForceParameter(j, NOTSET);
                    }
                    if (bSwapParity)
                    {
                        /* flag to swap parity of vsi  te construction */
                        param.setForceParameter(1, -1);
                    }
                }
                else
                {
                    if (bZero)
                    {
                        fprintf(stderr,
                                "NOTE: No default %s types, using zeroes\n",
                                interaction_function[ftype].longname);
                    }
                    else
                    {
                        auto message = gmx::formatString("No default %s types",
                                                         interaction_function[ftype].longname);
                        wi->addError(message);
                    }
                }
            }
            else
            {
                if (bSwapParity)
                {
                    switch (ftype)
                    {
                        case F_VSITE3FAD: param.setForceParameter(0, 360 - param.c0()); break;
                        case F_VSITE3OUT: param.setForceParameter(2, -param.c2()); break;
                    }
                }
            }
            if (!bFoundB)
            {
                /* We only have to issue a warning if these atoms are perturbed! */
                bool                     bPert      = false;
                gmx::ArrayRef<const int> paramAtoms = param.atoms();
                for (int j = 0; (j < nral); j++)
                {
                    bPert = bPert || PERTURBED(at->atom[paramAtoms[j]]);
                }

                if (bPert)
                {
                    auto message = gmx::formatString(
                            "No default %s types for perturbed atoms, "
                            "using normal values",
                            interaction_function[ftype].longname);
                    wi->addWarning(message);
                }
            }
        }
    }

    gmx::ArrayRef<const real> paramValue = param.forceParam();
    if ((ftype == F_PDIHS || ftype == F_ANGRES || ftype == F_ANGRESZ) && paramValue[5] != paramValue[2])
    {
        auto message = gmx::formatString("%s multiplicity can not be perturbed %f!=%f",
                                         interaction_function[ftype].longname,
                                         paramValue[2],
                                         paramValue[5]);
        warning_error_and_exit(wi, message, FARGS);
    }

    if (IS_TABULATED(ftype) && param.c0() != param.c2())
    {
        auto message = gmx::formatString("%s table number can not be perturbed %d!=%d",
                                         interaction_function[ftype].longname,
                                         gmx::roundToInt(param.c0()),
                                         gmx::roundToInt(param.c0()));
        warning_error_and_exit(wi, message, FARGS);
    }

    /* Dont add R-B dihedrals where all parameters are zero (no interaction) */
    if (ftype == F_RBDIHS)
    {

        int nr = 0;
        for (int i = 0; i < NRFP(ftype); i++)
        {
            if (paramValue[i] != 0.0)
            {
                nr++;
            }
        }
        if (nr == 0)
        {
            return;
        }
    }

    /* Put the values in the appropriate arrays */
    add_param_to_list(&bond[ftype], param);

    /* Push additional torsions from FF for ftype==9 if we have them.
     * We have already checked that the A/B states do not differ in this case,
     * so we do not have to double-check that again, or the vsite stuff.
     * In addition, those torsions cannot be automatically perturbed.
     */
    if (bDef && ftype == F_PDIHS)
    {
        for (int i = 1; i < nparam_defA; i++)
        {
            /* Advance pointer! */
            foundAParameter += 2;
            gmx::ArrayRef<const real> forceParam = foundAParameter->forceParam();
            for (int j = 0; j < (NRFPA(ftype) + NRFPB(ftype)); j++)
            {
                param.setForceParameter(j, forceParam[j]);
            }
            /* And push the next term for this torsion */
            add_param_to_list(&bond[ftype], param);
        }
    }
}

void push_cmap(Directive                         d,
               gmx::ArrayRef<InteractionsOfType> bondtype,
               gmx::ArrayRef<InteractionsOfType> bond,
               t_atoms*                          at,
               PreprocessingAtomTypes*           atypes,
               char*                             line,
               WarningHandler*                   wi)
{
    const char* aaformat[] = { "%d%d%d%d%d%d", "%d%d%d%d%d%d%d", "%d%d%d%d%d%d%d%d" };

    int  ftype, nral, nread, ncmap_params;
    int  cmap_type;
    int  aa[MAXATOMLIST + 1];
    bool bFound;

    ftype        = ifunc_index(d, 1);
    nral         = NRAL(ftype);
    ncmap_params = 0;

    nread = sscanf(line, aaformat[0], &aa[0], &aa[1], &aa[2], &aa[3], &aa[4], &aa[nral]);

    if (nread < nral)
    {
        too_few(wi);
        return;
    }
    else if (nread == nral)
    {
        ftype = ifunc_index(d, 1);
    }
    GMX_RELEASE_ASSERT(aa[nral] == 1, "Invalid function type for cmap torsion: must be 1");

    /* Check for double atoms and atoms out of bounds */
    for (int i = 0; i < nral; i++)
    {
        if (aa[i] < 1 || aa[i] > at->nr)
        {
            auto message = gmx::formatString(
                    "Atom index (%d) in %s out of bounds (1-%d).\n"
                    "This probably means that you have inserted topology section \"%s\"\n"
                    "in a part belonging to a different molecule than you intended to.\n"
                    "In that case move the \"%s\" section to the right molecule.",
                    aa[i],
                    enumValueToString(d),
                    at->nr,
                    enumValueToString(d),
                    enumValueToString(d));
            warning_error_and_exit(wi, message, FARGS);
        }

        for (int j = i + 1; (j < nral); j++)
        {
            if (aa[i] == aa[j])
            {
                auto message = gmx::formatString(
                        "Duplicate atom index (%d) in %s", aa[i], enumValueToString(d));
                wi->addError(message);
            }
        }
    }

    /* default force parameters  */
    std::vector<int> atoms;
    for (int j = 0; (j < nral); j++)
    {
        atoms.emplace_back(aa[j] - 1);
    }
    std::array<real, MAXFORCEPARAM> forceParam = { 0.0 };
    InteractionOfType               param(atoms, forceParam, "");
    /* Get the cmap type for this cmap angle */
    bFound = default_cmap_params(bondtype, at, atypes, &param, FALSE, &cmap_type, &ncmap_params, wi);

    /* We want exactly one parameter (the cmap type in state A (currently no state B) back */
    if (bFound && ncmap_params == 1)
    {
        /* Put the values in the appropriate arrays */
        param.setForceParameter(0, cmap_type);
        add_param_to_list(&bond[ftype], param);
    }
    else
    {
        /* This is essentially the same check as in default_cmap_params() done one more time */
        auto message =
                gmx::formatString("Unable to assign a cmap type to torsion %d %d %d %d and %d\n",
                                  param.ai() + 1,
                                  param.aj() + 1,
                                  param.ak() + 1,
                                  param.al() + 1,
                                  param.am() + 1);
        warning_error_and_exit(wi, message, FARGS);
    }
}


void push_vsitesn(Directive d, gmx::ArrayRef<InteractionsOfType> bond, t_atoms* at, char* line, WarningHandler* wi)
{
    char*   ptr;
    int     type, ftype, n, ret, nj, a;
    int*    atc    = nullptr;
    double *weight = nullptr, weight_tot;

    std::array<real, MAXFORCEPARAM> forceParam = { 0.0 };
    ptr                                        = line;
    ret                                        = sscanf(ptr, "%d%n", &a, &n);
    ptr += n;
    if (ret == 0)
    {
        auto message =
                gmx::formatString("Expected an atom index in section \"%s\"", enumValueToString(d));
        warning_error_and_exit(wi, message, FARGS);
    }

    sscanf(ptr, "%d%n", &type, &n);
    ptr += n;
    ftype         = ifunc_index(d, type);
    int firstAtom = a - 1;

    weight_tot = 0;
    nj         = 0;
    do
    {
        ret = sscanf(ptr, "%d%n", &a, &n);
        ptr += n;
        if (ret > 0)
        {
            if (nj % 20 == 0)
            {
                srenew(atc, nj + 20);
                srenew(weight, nj + 20);
            }
            atc[nj] = a - 1;
            switch (type)
            {
                case 1: weight[nj] = 1; break;
                case 2:
                    /* Here we use the A-state mass as a parameter.
                     * Note that the B-state mass has no influence.
                     */
                    weight[nj] = at->atom[atc[nj]].m;
                    break;
                case 3:
                    weight[nj] = -1;
                    ret        = sscanf(ptr, "%lf%n", &(weight[nj]), &n);
                    ptr += n;
                    if (weight[nj] < 0)
                    {
                        auto message = gmx::formatString(
                                "No weight or negative weight found for vsiten "
                                "constructing atom %d (atom index %d)",
                                nj + 1,
                                atc[nj] + 1);
                        warning_error_and_exit(wi, message, FARGS);
                    }
                    break;
                default:
                    auto message = gmx::formatString("Unknown vsiten type %d", type);
                    warning_error_and_exit(wi, message, FARGS);
            }
            weight_tot += weight[nj];
            nj++;
        }
    } while (ret > 0);

    if (nj == 0)
    {
        auto message = gmx::formatString("Expected more than one atom index in section \"%s\"",
                                         enumValueToString(d));
        warning_error_and_exit(wi, message, FARGS);
    }

    if (weight_tot == 0)
    {
        warning_error_and_exit(wi, "The total mass of the construting atoms is zero", FARGS);
    }

    for (int j = 0; j < nj; j++)
    {
        std::vector<int> atoms = { firstAtom, atc[j] };
        forceParam[0]          = nj;
        forceParam[1]          = weight[j] / weight_tot;
        /* Put the values in the appropriate arrays */
        add_param_to_list(&bond[ftype], InteractionOfType(atoms, forceParam));
    }

    sfree(atc);
    sfree(weight);
}

void push_mol(gmx::ArrayRef<MoleculeInformation> mols, char* pline, int* whichmol, int* nrcopies, WarningHandler* wi)
{
    char type[STRLEN];

    if (sscanf(pline, "%s%d", type, nrcopies) != 2)
    {
        too_few(wi);
        return;
    }

    /* Search moleculename.
     * Here we originally only did case insensitive matching. But because
     * some PDB files can have many chains and use case to generate more
     * chain-identifiers, which in turn end up in our moleculetype name,
     * we added support for case-sensitivity.
     */
    int nrcs    = 0;
    int nrci    = 0;
    int matchci = -1;
    int matchcs = -1;
    int i       = 0;
    for (const auto& mol : mols)
    {
        if (strcmp(type, *(mol.name)) == 0)
        {
            nrcs++;
            matchcs = i;
        }
        if (gmx_strcasecmp(type, *(mol.name)) == 0)
        {
            nrci++;
            matchci = i;
        }
        i++;
    }

    if (nrcs == 1)
    {
        // select the case sensitive match
        *whichmol = matchcs;
    }
    else
    {
        // avoid matching case-insensitive when we have multiple matches
        if (nrci > 1)
        {
            auto message = gmx::formatString(
                    "For moleculetype '%s' in [ system ] %d case insensitive "
                    "matches, but %d case sensitive matches were found. Check "
                    "the case of the characters in the moleculetypes.",
                    type,
                    nrci,
                    nrcs);
            warning_error_and_exit(wi, message, FARGS);
        }
        if (nrci == 1)
        {
            // select the unique case insensitive match
            *whichmol = matchci;
        }
        else
        {
            auto message = gmx::formatString("No such moleculetype %s", type);
            warning_error_and_exit(wi, message, FARGS);
        }
    }
}

void push_excl(char* line, gmx::ArrayRef<gmx::ExclusionBlock> b2, WarningHandler* wi)
{
    int  i, j;
    int  n;
    char base[STRLEN], format[STRLEN];

    if (sscanf(line, "%d", &i) == 0)
    {
        return;
    }

    if ((1 <= i) && (i <= b2.ssize()))
    {
        i--;
    }
    else
    {
        return;
    }
    strcpy(base, "%*d");
    do
    {
        strcpy(format, base);
        strcat(format, "%d");
        n = sscanf(line, format, &j);
        if (n == 1)
        {
            if ((1 <= j) && (j <= b2.ssize()))
            {
                j--;
                b2[i].atomNumber.push_back(j);
                /* also add the reverse exclusion! */
                b2[j].atomNumber.push_back(i);
                strcat(base, "%*d");
            }
            else
            {
                auto message = gmx::formatString("Invalid Atomnr j: %d, b2->nr: %zu\n", j, b2.size());
                warning_error_and_exit(wi, message, FARGS);
            }
        }
    } while (n == 1);
}

int add_atomtype_decoupled(PreprocessingAtomTypes* at, t_nbparam*** nbparam, t_nbparam*** pair)
{
    t_atom atom;
    int    nr;

    /* Define an atom type with all parameters set to zero (no interactions) */
    atom.q = 0.0;
    atom.m = 0.0;
    /* Type for decoupled atoms could be anything,
     * this should be changed automatically later when required.
     */
    atom.ptype = ParticleType::Atom;

    std::array<real, MAXFORCEPARAM> forceParam = { 0.0 };
    nr = at->addType(atom, "decoupled", InteractionOfType({}, forceParam, ""), -1, 0);

    /* Add space in the non-bonded parameters matrix */
    realloc_nb_params(at, nbparam, pair);

    return nr;
}

static void convert_pairs_to_pairsQ(gmx::ArrayRef<InteractionsOfType> interactions, real fudgeQQ, t_atoms* atoms)
{
    /* Add the pair list to the pairQ list */
    std::vector<InteractionOfType> paramnew;

    gmx::ArrayRef<const InteractionOfType> paramp1 = interactions[F_LJ14].interactionTypes;
    gmx::ArrayRef<const InteractionOfType> paramp2 = interactions[F_LJC14_Q].interactionTypes;

    /* Fill in the new F_LJC14_Q array with the old one. NOTE:
       it may be possible to just ADD the converted F_LJ14 array
       to the old F_LJC14_Q array, but since we have to create
       a new sized memory structure, better just to deep copy it all.
     */


    for (const auto& param : paramp2)
    {
        paramnew.emplace_back(param);
    }

    for (const auto& param : paramp1)
    {
        std::vector<real> forceParam = {
            fudgeQQ, atoms->atom[param.ai()].q, atoms->atom[param.aj()].q, param.c0(), param.c1()
        };
        paramnew.emplace_back(param.atoms(), forceParam, "");
    }

    /* now assign the new data to the F_LJC14_Q structure */
    interactions[F_LJC14_Q].interactionTypes = paramnew;

    /* Empty the LJ14 pairlist */
    interactions[F_LJ14].interactionTypes.clear();
}

static void generate_LJCpairsNB(MoleculeInformation* mol, int nb_funct, InteractionsOfType* nbp, WarningHandler* wi)
{
    int     n, ntype;
    t_atom* atom;

    n    = mol->atoms.nr;
    atom = mol->atoms.atom;

    ntype = static_cast<int>(std::sqrt(static_cast<double>(nbp->size())));
    GMX_ASSERT(ntype * ntype == gmx::ssize(*nbp),
               "Number of pairs of generated non-bonded parameters should be a perfect square");

    /* Add a pair interaction for all non-excluded atom pairs */
    const auto& excls = mol->excls;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            bool pairIsExcluded = false;
            for (const int atomK : excls[i])
            {
                if (atomK == j)
                {
                    pairIsExcluded = true;
                }
            }
            if (!pairIsExcluded)
            {
                if (nb_funct != F_LJ)
                {
                    auto message = gmx::formatString(
                            "Can only generate non-bonded pair interactions "
                            "for Van der Waals type Lennard-Jones");
                    warning_error_and_exit(wi, message, FARGS);
                }
                std::vector<int>  atoms      = { i, j };
                std::vector<real> forceParam = {
                    atom[i].q,
                    atom[j].q,
                    nbp->interactionTypes[ntype * atom[i].type + atom[j].type].c0(),
                    nbp->interactionTypes[ntype * atom[i].type + atom[j].type].c1()
                };
                add_param_to_list(&mol->interactions[F_LJC_PAIRS_NB], InteractionOfType(atoms, forceParam));
            }
        }
    }
}

static void set_excl_all(gmx::ListOfLists<int>* excl)
{
    /* Get rid of the current exclusions and exclude all atom pairs */
    const int        numAtoms = excl->ssize();
    std::vector<int> exclusionsForAtom(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        exclusionsForAtom[i] = i;
    }
    excl->clear();
    for (int i = 0; i < numAtoms; i++)
    {
        excl->pushBack(exclusionsForAtom);
    }
}

static void decouple_atoms(t_atoms*        atoms,
                           int             atomtype_decouple,
                           int             couple_lam0,
                           int             couple_lam1,
                           const char*     mol_name,
                           WarningHandler* wi)
{
    int i;

    for (i = 0; i < atoms->nr; i++)
    {
        t_atom* atom;

        atom = &atoms->atom[i];

        if (atom->qB != atom->q || atom->typeB != atom->type)
        {
            auto message = gmx::formatString(
                    "Atom %d in molecule type '%s' has different A and B state "
                    "charges and/or atom types set in the topology file as well "
                    "as through the mdp option '%s'. You can not use both "
                    "these methods simultaneously.",
                    i + 1,
                    mol_name,
                    "couple-moltype");
            warning_error_and_exit(wi, message, FARGS);
        }

        if (couple_lam0 == ecouplamNONE || couple_lam0 == ecouplamVDW)
        {
            atom->q = 0.0;
        }
        if (couple_lam0 == ecouplamNONE || couple_lam0 == ecouplamQ)
        {
            atom->type = atomtype_decouple;
        }
        if (couple_lam1 == ecouplamNONE || couple_lam1 == ecouplamVDW)
        {
            atom->qB = 0.0;
        }
        if (couple_lam1 == ecouplamNONE || couple_lam1 == ecouplamQ)
        {
            atom->typeB = atomtype_decouple;
        }
    }
}

void convert_moltype_couple(MoleculeInformation* mol,
                            int                  atomtype_decouple,
                            real                 fudgeQQ,
                            int                  couple_lam0,
                            int                  couple_lam1,
                            bool                 bCoupleIntra,
                            int                  nb_funct,
                            InteractionsOfType*  nbp,
                            WarningHandler*      wi)
{
    convert_pairs_to_pairsQ(mol->interactions, fudgeQQ, &mol->atoms);
    if (!bCoupleIntra)
    {
        generate_LJCpairsNB(mol, nb_funct, nbp, wi);
        set_excl_all(&mol->excls);
    }
    decouple_atoms(&mol->atoms, atomtype_decouple, couple_lam0, couple_lam1, *mol->name, wi);
}
