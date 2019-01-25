/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "resall.h"

#include <cctype>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

void read_atype(const char *ffdir, t_symtab *tab, PreprocessingAtomType *at)
{
    FILE                    *in;
    char                     buf[STRLEN], name[STRLEN];
    double                   m;
    int                      nratt = 0;
    t_atom                  *a;

    std::vector<std::string> files = fflib_search_file_end(ffdir, ".atp", TRUE);
    snew(a, 1);
    t_param                  nb;

    for (const auto &filename : files)
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
            }
            while ((feof(in) == 0) && strlen(buf) == 0);

            if (sscanf(buf, "%s%lf", name, &m) == 2)
            {
                a->m = m;
                at->addType(tab, a, name, &nb, 0, 0);
                fprintf(stderr, "\rAtomtype %d", ++nratt);
                fflush(stderr);
            }
            else
            {
                fprintf(stderr, "\nInvalid format: %s\n", buf);
            }
        }
        gmx_ffclose(in);
    }
    fprintf(stderr, "\n");
}

static void print_resatoms(FILE *out, PreprocessingAtomType *atype, const t_restp &rtp)
{
    /* fprintf(out,"%5s\n",rtp->resname);
       fprintf(out,"%5d\n",rtp->natom); */
    fprintf(out, "[ %s ]\n", rtp.resname.c_str());
    fprintf(out, " [ atoms ]\n");

    for (int j = 0; (j < rtp.natom()); j++)
    {
        int         tp   = rtp.atom[j].type;
        const char* tpnm = atype->atomNameFromType(tp);
        if (tpnm == nullptr)
        {
            gmx_fatal(FARGS, "Incorrect atomtype (%d)", tp);
        }
        fprintf(out, "%6s  %6s  %8.3f  %6d\n",
                *(rtp.atomname[j]), tpnm, rtp.atom[j].q, rtp.cgnr[j]);
    }
}

static bool read_atoms(FILE *in, char *line,
                       t_restp *r0, t_symtab *tab, PreprocessingAtomType *atype)
{
    char   buf[256], buf1[256];

    /* Read Atoms */
    r0->atom.clear();
    r0->atomname.clear();
    r0->cgnr.clear();
    while (get_a_line(in, line, STRLEN) && (strchr(line, '[') == nullptr))
    {
        int    cg;
        double q;
        if (sscanf(line, "%s%s%lf%d", buf, buf1, &q, &cg) != 4)
        {
            return FALSE;
        }
        r0->atomname.push_back(put_symtab(tab, buf));
        r0->atom.push_back(t_atom());
        r0->atom.back().q   = q;
        r0->cgnr.push_back(cg);
        int j               = atype->atomTypeFromString(buf1);
        if (j == NOTSET)
        {
            gmx_fatal(FARGS, "Atom type %s (residue %s) not found in atomtype "
                      "database", buf1, r0->resname.c_str());
        }
        r0->atom.back().type = j;
        r0->atom.back().m    = atype->atomMassAFromType(j);
    }

    return TRUE;
}

static bool read_bondeds(int bt, FILE *in, char *line, t_restp *rtp)
{
    char str[STRLEN];

    while (get_a_line(in, line, STRLEN) && (strchr(line, '[') == nullptr))
    {
        int n = 0;
        int j;
        int ni;
        rtp->rb[bt].b.push_back(t_rbonded());
        for (j = 0; j < btsNiatoms[bt]; j++)
        {
            if (sscanf(line+n, "%s%n", str, &ni) == 1)
            {
                rtp->rb[bt].b.back().a[j] = str;
            }
            else
            {
                return FALSE;
            }
            n += ni;
        }
        for (; j < MAXATOMLIST; j++)
        {
            rtp->rb[bt].b.back().a[j].clear();
        }
        while (isspace(line[n]))
        {
            n++;
        }
        rtrim(line+n);
        rtp->rb[bt].b.back().s = line+n;
    }

    return TRUE;
}

static void print_resbondeds(FILE *out, int bt, const t_restp &rtp)
{
    if (rtp.rb[bt].nb())
    {
        fprintf(out, " [ %s ]\n", btsNames[bt]);

        for (int i = 0; i < rtp.rb[bt].nb(); i++)
        {
            for (int j = 0; j < btsNiatoms[bt]; j++)
            {
                fprintf(out, "%6s ", rtp.rb[bt].b[i].a[j].c_str());
            }
            if (rtp.rb[bt].b[i].s[0])
            {
                fprintf(out, "    %s", rtp.rb[bt].b[i].s.c_str());
            }
            fprintf(out, "\n");
        }
    }
}

static void check_rtp(gmx::ArrayRef<const t_restp> rtp, const char *libfn)
{
    /* check for double entries, assuming list is already sorted */
    for (int i = 1; (i < rtp.size()); i++)
    {
        if (gmx::equalCaseInsensitive(rtp[i-1].resname, rtp[i].resname))
        {
            fprintf(stderr, "WARNING double entry %s in file %s\n",
                    rtp[i].resname.c_str(), libfn);
        }
    }
}

static int get_bt(char* header)
{
    for (int i = 0; i < ebtsNR; i++)
    {
        if (gmx_strcasecmp(btsNames[i], header) == 0)
        {
            return i;
        }
    }
    return NOTSET;
}

static void clear_t_restp(t_restp *rrtp)
{
    rrtp->~t_restp();
}


/* print all the ebtsNR type numbers */
static void print_resall_header(FILE *out, gmx::ArrayRef<const t_restp> rtp)
{
    fprintf(out, "[ bondedtypes ]\n");
    fprintf(out, "; bonds  angles  dihedrals  impropers all_dihedrals nr_exclusions  HH14  remove_dih\n");
    fprintf(out, " %5d  %6d  %9d  %9d  %14d  %14d %14d %14d\n\n",
            rtp[0].rb[0].type,
            rtp[0].rb[1].type,
            rtp[0].rb[2].type,
            rtp[0].rb[3].type,
            static_cast<int>(rtp[0].bKeepAllGeneratedDihedrals),
            rtp[0].nrexcl,
            static_cast<int>(rtp[0].bGenerateHH14Interactions),
            static_cast<int>(rtp[0].bRemoveDihedralIfWithImproper));
}

void print_resall(FILE *out, gmx::ArrayRef<t_restp> rtp,
                  PreprocessingAtomType *atype)
{
    if (rtp.empty())
    {
        return;
    }

    print_resall_header(out, rtp);

    for (int i = 0; i < rtp.size(); i++)
    {
        if (rtp[i].natom() > 0)
        {
            print_resatoms(out, atype, rtp[i]);
            for (int bt = 0; bt < ebtsNR; bt++)
            {
                print_resbondeds(out, bt, rtp[i]);
            }
        }
    }
}

void read_resall(const char *rrdb, std::vector<t_restp> *rtp,
                 PreprocessingAtomType *atype, t_symtab *tab,
                 bool bAllowOverrideRTP)
{
    FILE         *in;
    char          filebase[STRLEN], line[STRLEN], header[STRLEN];
    int           bt, nparam;
    int           dum1, dum2, dum3;
    bool          bNextResidue, bError;

    fflib_filename_base(rrdb, filebase, STRLEN);

    in = fflib_open(rrdb);

    t_restp header_settings;

    /* these bonded parameters will overwritten be when  *
     * there is a [ bondedtypes ] entry in the .rtp file */
    header_settings.rb[ebtsBONDS].type  = 1; /* normal bonds     */
    header_settings.rb[ebtsANGLES].type = 1; /* normal angles    */
    header_settings.rb[ebtsPDIHS].type  = 1; /* normal dihedrals */
    header_settings.rb[ebtsIDIHS].type  = 2; /* normal impropers */
    header_settings.rb[ebtsEXCLS].type  = 1; /* normal exclusions */
    header_settings.rb[ebtsCMAP].type   = 1; /* normal cmap torsions */

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
    if (gmx_strncasecmp("bondedtypes", header, 5) == 0)
    {
        get_a_line(in, line, STRLEN);
        if ((nparam = sscanf(line, "%d %d %d %d %d %d %d %d",
                             &header_settings.rb[ebtsBONDS].type, &header_settings.rb[ebtsANGLES].type,
                             &header_settings.rb[ebtsPDIHS].type, &header_settings.rb[ebtsIDIHS].type,
                             &dum1, &header_settings.nrexcl, &dum2, &dum3)) < 4)
        {
            gmx_fatal(FARGS, "need 4 to 8 parameters in the header of .rtp file %s at line:\n%s\n", rrdb, line);
        }
        header_settings.bKeepAllGeneratedDihedrals    = (dum1 != 0);
        header_settings.bGenerateHH14Interactions     = (dum2 != 0);
        header_settings.bRemoveDihedralIfWithImproper = (dum3 != 0);
        get_a_line(in, line, STRLEN);
        if (nparam < 5)
        {
            fprintf(stderr, "Using default: not generating all possible dihedrals\n");
            header_settings.bKeepAllGeneratedDihedrals = FALSE;
        }
        if (nparam < 6)
        {
            fprintf(stderr, "Using default: excluding 3 bonded neighbors\n");
            header_settings.nrexcl = 3;
        }
        if (nparam < 7)
        {
            fprintf(stderr, "Using default: generating 1,4 H--H interactions\n");
            header_settings.bGenerateHH14Interactions = TRUE;
        }
        if (nparam < 8)
        {
            fprintf(stderr, "Using default: removing proper dihedrals found on the same bond as a proper dihedral\n");
            header_settings.bRemoveDihedralIfWithImproper = TRUE;
        }
    }
    else
    {
        fprintf(stderr,
                "Reading .rtp file without '[ bondedtypes ]' directive,\n"
                "Will proceed as if the entry was:\n");
        print_resall_header(stderr, gmx::arrayRefFromArray(&header_settings, 1));
    }
    /* We don't know the current size of rrtp, but simply realloc immediately */
    std::vector<t_restp> rrtp   = *rtp;
    while (!feof(in))
    {
        /* Initialise rtp entry structure */
        rrtp.push_back(header_settings);
        if (!get_header(line, header))
        {
            gmx_fatal(FARGS, "in .rtp file at line:\n%s\n", line);
        }
        rrtp.back().resname  = header;
        rrtp.back().filebase = filebase;

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
                bt = get_bt(header);
                if (bt != NOTSET)
                {
                    /* header is an bonded directive */
                    bError = !read_bondeds(bt, in, line, &rrtp.back());
                }
                else if (gmx_strncasecmp("atoms", header, 5) == 0)
                {
                    /* header is the atoms directive */
                    bError = !read_atoms(in, line, &(rrtp.back()), tab, atype);
                }
                else
                {
                    /* else header must be a residue name */
                    bNextResidue = TRUE;
                }
            }
            if (bError)
            {
                gmx_fatal(FARGS, "in .rtp file in residue %s at line:\n%s\n",
                          rrtp.back().resname.c_str(), line);
            }
        }
        while ((feof(in) == 0) && !bNextResidue);

        if (rrtp.back().atom.empty())
        {
            gmx_fatal(FARGS, "No atoms found in .rtp file in residue %s\n",
                      rrtp.back().resname.c_str());
        }

        int firstrtp = -1;
        for (int i = 0; i < gmx::index(rrtp.size() - 1); i++)
        {
            if (gmx::equalCaseInsensitive(rrtp[i].resname, rrtp.back().resname))
            {
                firstrtp = i;
            }
        }

        if (firstrtp == -1)
        {
            fprintf(stderr, "\rResidue %lu", rrtp.size());
            fflush(stderr);
        }
        else
        {
            if (firstrtp >= gmx::index(rtp->size()))
            {
                gmx_fatal(FARGS, "Found a second entry for '%s' in '%s'",
                          rrtp.back().resname.c_str(), rrdb);
            }
            if (bAllowOverrideRTP)
            {
                fprintf(stderr, "\n");
                fprintf(stderr, "Found another rtp entry for '%s' in '%s', ignoring this entry and keeping the one from '%s.rtp'\n",
                        rrtp.back().resname.c_str(), rrdb, rrtp[firstrtp].filebase.c_str());
                /* We should free all the data for this entry.
                 * The current code gives a lot of dangling pointers.
                 */
                clear_t_restp(&rrtp.back());
                rrtp.erase(rrtp.end() - 1);
            }
            else
            {
                gmx_fatal(FARGS, "Found rtp entries for '%s' in both '%s' and '%s'. If you want the first definition to override the second one, set the -rtpo option of pdb2gmx.", rrtp.back().resname.c_str(), rrtp[firstrtp].filebase.c_str(), rrdb);
            }
        }
    }
    gmx_ffclose(in);

    check_rtp(rrtp, rrdb);

    *rtp     = rrtp;
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
static int neq_str_sign(const char *a1, const char *a2)
{
    int l1, l2, lm;

    l1 = static_cast<int>(strlen(a1));
    l2 = static_cast<int>(strlen(a2));
    lm = std::min(l1, l2);

    if (lm >= 1 &&
        ((l1 == l2+1 && is_sign(a1[l1-1])) ||
         (l2 == l1+1 && is_sign(a2[l2-1]))) &&
        gmx_strncasecmp(a1, a2, lm) == 0)
    {
        return lm;
    }
    else
    {
        return 0;
    }
}

std::string search_rtp(const char *key, gmx::ArrayRef<t_restp> rtp)
{
    int  nbest, best, besti;
    char bestbuf[STRLEN];

    nbest =  0;
    besti = -1;
    /* We want to match at least one character */
    best  =  1;
    for (int i = 0; (i < rtp.size()); i++)
    {
        if (gmx::equalCaseInsensitive(key, rtp[i].resname))
        {
            besti = i;
            nbest = 1;
            break;
        }
        else
        {
            /* Allow a mismatch of at most a sign character (with warning) */
            int n = neq_str_sign(key, rtp[i].resname.c_str());
            if (n >= best &&
                n+1 >= static_cast<int>(strlen(key)) &&
                n+1 >= static_cast<int>(rtp[i].resname.length()))
            {
                if (n == best)
                {
                    if (nbest == 1)
                    {
                        strcpy(bestbuf, rtp[besti].resname.c_str());
                    }
                    if (nbest >= 1)
                    {
                        strcat(bestbuf, " or ");
                        strcat(bestbuf, rtp[i].resname.c_str());
                    }
                }
                else
                {
                    nbest = 0;
                }
                besti = i;
                best  = n;
                nbest++;
            }
        }
    }
    if (nbest > 1)
    {
        gmx_fatal(FARGS, "Residue '%s' not found in residue topology database, looks a bit like %s", key, bestbuf);
    }
    else if (besti == -1)
    {
        gmx_fatal(FARGS, "Residue '%s' not found in residue topology database", key);
    }
    if (!gmx::equalCaseInsensitive(rtp[besti].resname, key))
    {
        fprintf(stderr,
                "\nWARNING: '%s' not found in residue topology database, "
                "trying to use '%s'\n\n", key, rtp[besti].resname.c_str());
    }

    return rtp[besti].resname;
}

RestPIt get_restp(const char *rtpname, gmx::ArrayRef<t_restp> rtp)
{
    auto found = std::find_if(rtp.begin(), rtp.end(),
                              [&rtpname](const t_restp &rp)
                              { return gmx::equalCaseInsensitive(rtpname, rp.resname); });
    if (found == rtp.end())
    {

        /* This should never happen, since search_rtp should have been called
         * before calling get_restp.
         */
        gmx_fatal(FARGS, "Residue type '%s' not found in residue topology database", rtpname);
    }

    return found;
}
