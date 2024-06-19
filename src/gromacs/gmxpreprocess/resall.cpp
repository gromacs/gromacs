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

#include "resall.h"

#include <cctype>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringtoenumvalueconverter.h"
#include "gromacs/utility/stringutil.h"

#include "hackblock.h"

PreprocessingAtomTypes read_atype(const std::filesystem::path& ffdir)
{
    FILE*  in;
    char   buf[STRLEN], name[STRLEN];
    double m;
    auto   atom = std::make_unique<t_atom>();

    auto                   files = fflib_search_file_end(ffdir, ".atp", TRUE);
    PreprocessingAtomTypes at;

    for (const auto& filename : files)
    {
        in = fflib_open(filename);
        while (!feof(in))
        {
            /* Skip blank or comment-only lines */
            do
            {
                if (fgets2(buf, STRLEN, in) != nullptr)
                {
                    strip_comment(buf);
                    trim(buf);
                }
            } while ((feof(in) == 0) && strlen(buf) == 0);

            if (sscanf(buf, "%s%lf", name, &m) == 2)
            {
                atom->m = m;
                at.addType(*atom, name, InteractionOfType({}, {}), 0, 0);
            }
            else
            {
                gmx_fatal(FARGS, "Invalid atomtype format: '%s'", buf);
            }
        }
        gmx_ffclose(in);
    }
    return at;
}

static void print_resatoms(FILE* out, const PreprocessingAtomTypes& atype, const PreprocessResidue& rtpDBEntry)
{
    fprintf(out, "[ %s ]\n", rtpDBEntry.resname.c_str());
    fprintf(out, " [ atoms ]\n");

    for (int j = 0; (j < rtpDBEntry.natom()); j++)
    {
        int  tp   = rtpDBEntry.atom[j].type;
        auto tpnm = atype.atomNameFromAtomType(tp);
        if (!tpnm.has_value())
        {
            gmx_fatal(FARGS, "Incorrect atomtype (%d)", tp);
        }
        fprintf(out,
                "%6s  %6s  %8.3f  %6d\n",
                *(rtpDBEntry.atomname[j]),
                tpnm->c_str(),
                rtpDBEntry.atom[j].q,
                rtpDBEntry.cgnr[j]);
    }
}

static bool read_atoms(FILE* in, char* line, PreprocessResidue* r0, t_symtab* tab, PreprocessingAtomTypes* atype)
{
    int    cg;
    char   buf[256], buf1[256];
    double q;

    /* Read Atoms */
    r0->atom.clear();
    r0->atomname.clear();
    r0->cgnr.clear();
    while (get_a_line(in, line, STRLEN) && (strchr(line, '[') == nullptr))
    {
        if (sscanf(line, "%s%s%lf%d", buf, buf1, &q, &cg) != 4)
        {
            return FALSE;
        }
        r0->atomname.push_back(put_symtab(tab, buf));
        r0->atom.emplace_back();
        r0->atom.back().q = q;
        r0->cgnr.push_back(cg);
        auto j = atype->atomTypeFromName(buf1);
        if (!j.has_value())
        {
            gmx_fatal(FARGS,
                      "Atom type %s (residue %s) not found in atomtype "
                      "database",
                      buf1,
                      r0->resname.c_str());
        }
        r0->atom.back().type = *j;
        r0->atom.back().m    = *atype->atomMassFromAtomType(*j);
    }

    return TRUE;
}

static bool read_bondeds(BondedTypes bt, FILE* in, char* line, PreprocessResidue* rtpDBEntry)
{
    char str[STRLEN];

    while (get_a_line(in, line, STRLEN) && (strchr(line, '[') == nullptr))
    {
        int n = 0;
        int ni;
        rtpDBEntry->rb[bt].b.emplace_back();
        BondedInteraction* newBond = &rtpDBEntry->rb[bt].b.back();
        for (int j = 0; j < enumValueToNumIAtoms(bt); j++)
        {
            if (sscanf(line + n, "%s%n", str, &ni) == 1)
            {
                newBond->a[j] = str;
            }
            else
            {
                return FALSE;
            }
            n += ni;
        }
        while (isspace(line[n]))
        {
            n++;
        }
        rtrim(line + n);
        newBond->s = line + n;
    }

    return TRUE;
}

static void print_resbondeds(FILE* out, BondedTypes bt, const PreprocessResidue& rtpDBEntry)
{
    if (!rtpDBEntry.rb[bt].b.empty())
    {
        fprintf(out, " [ %s ]\n", enumValueToString(bt));

        for (const auto& b : rtpDBEntry.rb[bt].b)
        {
            for (int j = 0; j < enumValueToNumIAtoms(bt); j++)
            {
                fprintf(out, "%6s ", b.a[j].c_str());
            }
            if (!b.s.empty())
            {
                fprintf(out, "    %s", b.s.c_str());
            }
            fprintf(out, "\n");
        }
    }
}

static void check_rtp(gmx::ArrayRef<const PreprocessResidue> rtpDBEntry,
                      const std::filesystem::path&           libfn,
                      const gmx::MDLogger&                   logger)
{
    /* check for double entries, assuming list is already sorted */
    for (auto it = rtpDBEntry.begin() + 1; it != rtpDBEntry.end(); it++)
    {
        auto prev = it - 1;
        if (gmx::equalCaseInsensitive(prev->resname, it->resname))
        {
            GMX_LOG(logger.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "Double entry %s in file %s", it->resname.c_str(), libfn.string().c_str());
        }
    }
}

static std::optional<BondedTypes> get_bt(char* header)
{
    gmx::StringToEnumValueConverter<BondedTypes, enumValueToString> converter;
    return converter.valueFrom(header);
}

/* print all the BondedTypes type numbers */
static void print_resall_header(FILE* out, gmx::ArrayRef<const PreprocessResidue> rtpDBEntry)
{
    fprintf(out, "[ bondedtypes ]\n");
    fprintf(out,
            "; bonds  angles  dihedrals  impropers all_dihedrals nr_exclusions  HH14  "
            "remove_dih\n");
    fprintf(out,
            " %5d  %6d  %9d  %9d  %14d  %14d %14d %14d\n\n",
            rtpDBEntry[0].rb[0].type,
            rtpDBEntry[0].rb[1].type,
            rtpDBEntry[0].rb[2].type,
            rtpDBEntry[0].rb[3].type,
            static_cast<int>(rtpDBEntry[0].bKeepAllGeneratedDihedrals),
            rtpDBEntry[0].nrexcl,
            static_cast<int>(rtpDBEntry[0].bGenerateHH14Interactions),
            static_cast<int>(rtpDBEntry[0].bRemoveDihedralIfWithImproper));
}


static void print_resall_log(const gmx::MDLogger& logger, gmx::ArrayRef<const PreprocessResidue> rtpDBEntry)
{
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("[ bondedtypes ]");
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "; bonds  angles  dihedrals  impropers all_dihedrals nr_exclusions  HH14  "
                    "remove_dih");
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(" %5d  %6d  %9d  %9d  %14d  %14d %14d %14d",
                                 rtpDBEntry[0].rb[0].type,
                                 rtpDBEntry[0].rb[1].type,
                                 rtpDBEntry[0].rb[2].type,
                                 rtpDBEntry[0].rb[3].type,
                                 static_cast<int>(rtpDBEntry[0].bKeepAllGeneratedDihedrals),
                                 rtpDBEntry[0].nrexcl,
                                 static_cast<int>(rtpDBEntry[0].bGenerateHH14Interactions),
                                 static_cast<int>(rtpDBEntry[0].bRemoveDihedralIfWithImproper));
}


void print_resall(FILE* out, gmx::ArrayRef<const PreprocessResidue> rtpDBEntry, const PreprocessingAtomTypes& atype)
{
    if (rtpDBEntry.empty())
    {
        return;
    }

    print_resall_header(out, rtpDBEntry);

    for (const auto& r : rtpDBEntry)
    {
        if (r.natom() > 0)
        {
            print_resatoms(out, atype, r);
            for (auto bt : gmx::EnumerationWrapper<BondedTypes>{})
            {
                print_resbondeds(out, bt, r);
            }
        }
    }
}

void readResidueDatabase(const std::filesystem::path&    rrdb,
                         std::vector<PreprocessResidue>* rtpDBEntry,
                         PreprocessingAtomTypes*         atype,
                         t_symtab*                       tab,
                         const gmx::MDLogger&            logger,
                         bool                            bAllowOverrideRTP)
{
    FILE* in;
    char  line[STRLEN], header[STRLEN];
    int   nparam;
    int   dum1, dum2, dum3;
    bool  bNextResidue, bError;

    auto filebase = fflib_filename_base(rrdb);

    in = fflib_open(rrdb);

    PreprocessResidue header_settings;

    /* these bonded parameters will overwritten be when  *
     * there is a [ bondedtypes ] entry in the .rtp file */
    header_settings.rb[BondedTypes::Bonds].type             = 1; /* normal bonds     */
    header_settings.rb[BondedTypes::Angles].type            = 1; /* normal angles    */
    header_settings.rb[BondedTypes::ProperDihedrals].type   = 1; /* normal dihedrals */
    header_settings.rb[BondedTypes::ImproperDihedrals].type = 2; /* normal impropers */
    header_settings.rb[BondedTypes::Exclusions].type        = 1; /* normal exclusions */
    header_settings.rb[BondedTypes::Cmap].type              = 1; /* normal cmap torsions */

    header_settings.bKeepAllGeneratedDihedrals    = FALSE;
    header_settings.nrexcl                        = 3;
    header_settings.bGenerateHH14Interactions     = TRUE;
    header_settings.bRemoveDihedralIfWithImproper = TRUE;

    /* Column 5 & 6 aren't really bonded types, but we include
     * them here to avoid introducing a new section:
     * Column 5 : This controls the generation of dihedrals from the bonding.
     *            All possible dihedrals are generated automatically. A value of
     *            1 here means that all these are retained. A value of
     *            0 here requires generated dihedrals be removed if
     *              * there are any dihedrals on the same central atoms
     *                specified in the residue topology, or
     *              * there are other identical generated dihedrals
     *                sharing the same central atoms, or
     *              * there are other generated dihedrals sharing the
     *                same central bond that have fewer hydrogen atoms
     * Column 6: Number of bonded neighbors to exclude.
     * Column 7: Generate 1,4 interactions between two hydrogen atoms
     * Column 8: Remove proper dihedrals if centered on the same bond
     *           as an improper dihedral
     */
    get_a_line(in, line, STRLEN);
    if (!get_header(line, header))
    {
        gmx_fatal(FARGS, "in .rtp file at line:\n%s\n", line);
    }
    if (gmx::equalCaseInsensitive("bondedtypes", header, 5))
    {
        get_a_line(in, line, STRLEN);
        if ((nparam = sscanf(line,
                             "%d %d %d %d %d %d %d %d",
                             &header_settings.rb[BondedTypes::Bonds].type,
                             &header_settings.rb[BondedTypes::Angles].type,
                             &header_settings.rb[BondedTypes::ProperDihedrals].type,
                             &header_settings.rb[BondedTypes::ImproperDihedrals].type,
                             &dum1,
                             &header_settings.nrexcl,
                             &dum2,
                             &dum3))
            < 4)
        {
            gmx_fatal(FARGS,
                      "need 4 to 8 parameters in the header of .rtp file %s at line:\n%s\n",
                      rrdb.string().c_str(),
                      line);
        }
        header_settings.bKeepAllGeneratedDihedrals    = (dum1 != 0);
        header_settings.bGenerateHH14Interactions     = (dum2 != 0);
        header_settings.bRemoveDihedralIfWithImproper = (dum3 != 0);
        get_a_line(in, line, STRLEN);
        if (nparam < 5)
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendTextFormatted("Using default: not generating all possible dihedrals");
            header_settings.bKeepAllGeneratedDihedrals = FALSE;
        }
        if (nparam < 6)
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendTextFormatted("Using default: excluding 3 bonded neighbors");
            header_settings.nrexcl = 3;
        }
        if (nparam < 7)
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendTextFormatted("Using default: generating 1,4 H--H interactions");
            header_settings.bGenerateHH14Interactions = TRUE;
        }
        if (nparam < 8)
        {
            GMX_LOG(logger.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "Using default: removing proper dihedrals found on the same bond as a "
                            "proper dihedral");
            header_settings.bRemoveDihedralIfWithImproper = TRUE;
        }
    }
    else
    {
        GMX_LOG(logger.warning)
                .asParagraph()
                .appendTextFormatted(
                        "Reading .rtp file without '[ bondedtypes ]' directive, "
                        "Will proceed as if the entry was:");
        print_resall_log(logger, gmx::arrayRefFromArray(&header_settings, 1));
    }
    /* We don't know the current size of rrtp, but simply realloc immediately */
    auto oldArrayEnd = rtpDBEntry->end();
    while (!feof(in))
    {
        /* Initialise rtp entry structure */
        rtpDBEntry->push_back(header_settings);
        PreprocessResidue* res = &rtpDBEntry->back();
        if (!get_header(line, header))
        {
            gmx_fatal(FARGS, "in .rtp file at line:\n%s\n", line);
        }
        res->resname  = header;
        res->filebase = filebase.string();

        get_a_line(in, line, STRLEN);
        bError       = FALSE;
        bNextResidue = FALSE;
        do
        {
            if (!get_header(line, header))
            {
                bError = TRUE;
            }
            else
            {
                auto bt = get_bt(header);
                if (bt.has_value())
                {
                    /* header is an bonded directive */
                    bError = !read_bondeds(*bt, in, line, res);
                }
                else if (gmx::equalCaseInsensitive("atoms", header, 5))
                {
                    /* header is the atoms directive */
                    bError = !read_atoms(in, line, res, tab, atype);
                }
                else
                {
                    /* else header must be a residue name */
                    bNextResidue = TRUE;
                }
            }
            if (bError)
            {
                gmx_fatal(FARGS, "in .rtp file in residue %s at line:\n%s\n", res->resname.c_str(), line);
            }
        } while ((feof(in) == 0) && !bNextResidue);

        if (res->natom() == 0)
        {
            gmx_fatal(FARGS, "No atoms found in .rtp file in residue %s\n", res->resname.c_str());
        }

        auto found = std::find_if(
                rtpDBEntry->begin(), rtpDBEntry->end() - 1, [&res](const PreprocessResidue& entry) {
                    return gmx::equalCaseInsensitive(entry.resname, res->resname);
                });

        if (found != rtpDBEntry->end() - 1)
        {
            if (found >= oldArrayEnd)
            {
                gmx_fatal(FARGS,
                          "Found a second entry for '%s' in '%s'",
                          res->resname.c_str(),
                          rrdb.string().c_str());
            }
            if (bAllowOverrideRTP)
            {
                GMX_LOG(logger.warning)
                        .asParagraph()
                        .appendTextFormatted(
                                "Found another rtp entry for '%s' in '%s',"
                                " ignoring this entry and keeping the one from '%s.rtp'",
                                res->resname.c_str(),
                                rrdb.string().c_str(),
                                found->filebase.c_str());
                /* We should free all the data for this entry.
                 * The current code gives a lot of dangling pointers.
                 */
                rtpDBEntry->erase(rtpDBEntry->end() - 1);
            }
            else
            {
                gmx_fatal(FARGS,
                          "Found rtp entries for '%s' in both '%s' and '%s'. If you want the first "
                          "definition to override the second one, set the -rtpo option of pdb2gmx.",
                          res->resname.c_str(),
                          found->filebase.c_str(),
                          rrdb.string().c_str());
            }
        }
    }
    gmx_ffclose(in);

    std::sort(rtpDBEntry->begin(), rtpDBEntry->end(), [](const PreprocessResidue& a, const PreprocessResidue& b) {
        return std::lexicographical_compare(
                a.resname.begin(),
                a.resname.end(),
                b.resname.begin(),
                b.resname.end(),
                [](const char& c1, const char& c2) { return std::toupper(c1) < std::toupper(c2); });
    });

    check_rtp(*rtpDBEntry, rrdb, logger);
}

/************************************************************
 *
 *                  SEARCH   ROUTINES
 *
 ***********************************************************/
static bool is_sign(char c)
{
    return (c == '+' || c == '-');
}

/* Compares if the strins match except for a sign at the end
 * of (only) one of the two.
 */
static int neq_str_sign(const char* a1, const char* a2)
{
    int l1, l2, lm;

    l1 = static_cast<int>(strlen(a1));
    l2 = static_cast<int>(strlen(a2));
    lm = std::min(l1, l2);

    if (lm >= 1 && ((l1 == l2 + 1 && is_sign(a1[l1 - 1])) || (l2 == l1 + 1 && is_sign(a2[l2 - 1])))
        && gmx::equalCaseInsensitive(a1, a2, lm))
    {
        return lm;
    }
    else
    {
        return 0;
    }
}

std::string searchResidueDatabase(const std::string&                     key,
                                  gmx::ArrayRef<const PreprocessResidue> rtpDBEntry,
                                  const gmx::MDLogger&                   logger)
{
    int         nbest, best, besti;
    std::string bestbuf;

    nbest = 0;
    besti = -1;
    /* We want to match at least one character */
    best = 1;
    for (auto it = rtpDBEntry.begin(); it != rtpDBEntry.end(); it++)
    {
        if (gmx::equalCaseInsensitive(key, it->resname))
        {
            besti = std::distance(rtpDBEntry.begin(), it);
            nbest = 1;
            break;
        }
        else
        {
            /* Allow a mismatch of at most a sign character (with warning) */
            int n = neq_str_sign(key.c_str(), it->resname.c_str());
            if (n >= best && n + 1 >= gmx::Index(key.length()) && n + 1 >= gmx::Index(it->resname.length()))
            {
                if (n == best)
                {
                    if (nbest == 1)
                    {
                        bestbuf = rtpDBEntry[besti].resname;
                    }
                    if (nbest >= 1)
                    {
                        bestbuf.append(" or ");
                        bestbuf.append(it->resname);
                    }
                }
                else
                {
                    nbest = 0;
                }
                besti = std::distance(rtpDBEntry.begin(), it);
                best  = n;
                nbest++;
            }
        }
    }
    if (nbest > 1)
    {
        gmx_fatal(FARGS,
                  "Residue '%s' not found in residue topology database, looks a bit like %s",
                  key.c_str(),
                  bestbuf.c_str());
    }
    else if (besti == -1)
    {
        gmx_fatal(FARGS, "Residue '%s' not found in residue topology database", key.c_str());
    }
    if (!gmx::equalCaseInsensitive(rtpDBEntry[besti].resname, key))
    {
        GMX_LOG(logger.warning)
                .asParagraph()
                .appendTextFormatted(
                        "'%s' not found in residue topology database, "
                        "trying to use '%s'",
                        key.c_str(),
                        rtpDBEntry[besti].resname.c_str());
    }

    return rtpDBEntry[besti].resname;
}

gmx::ArrayRef<const PreprocessResidue>::const_iterator
getDatabaseEntry(const std::string& rtpname, gmx::ArrayRef<const PreprocessResidue> rtpDBEntry)
{
    auto found = std::find_if(
            rtpDBEntry.begin(), rtpDBEntry.end(), [&rtpname](const PreprocessResidue& entry) {
                return gmx::equalCaseInsensitive(rtpname, entry.resname);
            });
    if (found == rtpDBEntry.end())
    {
        /* This should never happen, since searchResidueDatabase should have been called
         * before calling getDatabaseEntry.
         */
        gmx_fatal(FARGS, "Residue type '%s' not found in residue topology database", rtpname.c_str());
    }

    return found;
}
