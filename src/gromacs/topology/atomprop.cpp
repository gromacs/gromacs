/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "atomprop.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <memory>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

/* NOTFOUND should be smallest, others larger in increasing priority */
enum
{
    NOTFOUND = -4,
    WILDCARD,
    WILDPROT,
    PROTEIN
};

//! Basic entries in AtomProperty.
struct BaseEntry
{
    //! Default constructor.
    BaseEntry(const std::string& aName, const std::string& rName) :
        atomName(aName),
        residueName(rName),
        isAvailable(false),
        value(0.0)
    {
    }
    //! Name for atom.
    std::string atomName;
    //! Name for residue.
    std::string residueName;
    //! Is property available.
    bool isAvailable;
    //! Value set for property.
    real value;
};

//! Conglomeration of atom property entries.
struct AtomProperty
{
    //! Has property been set.
    bool isSet = false;
    //! Database the property is coming from.
    std::string db;
    //! Default value for property.
    double def = 0.0;
    //! Basic entries for properties.
    std::vector<BaseEntry> entry;
};

//! Implementation detail type for Atomproperties.
class AtomProperties::Impl
{
public:
    //! Should user be warned about error.
    bool bWarned = false;
    //! Should user be warned about vdW not found.
    bool bWarnVDW = false;
    //! The different atom properties.
    AtomProperty prop[epropNR];
    //! The residue types.
    ResidueType restype;
};

/*! \brief
 * Find number of matching characters in entry.
 *
 * If not all characters are matching, return NOTFOUND.
 * If the length of the database entry is different from the search,
 * also return NOTFOUND.
 *
 * \param[in] search Entry to compare to database.
 * \param[in] database Name of the database entry to compare to.
 * \returns Number of matching characters or NOTFOUND.
 */
static int compareToDatabase(const std::string& search, const std::string& database)
{
    if (database.length() > search.length())
    {
        return NOTFOUND;
    }
    size_t matches = 0;
    for (size_t i = 0; i < database.length(); i++)
    {
        if (search[i] == database[i])
        {
            matches++;
        }
    }
    if (matches == database.length())
    {
        return matches;
    }
    else
    {
        return NOTFOUND;
    }
}

/*! \brief
 * Finds the index for the property being searched.
 *
 * \param[in] ap Property to search for.
 * \param[in] restype Residuetypes in database.
 * \param[in] residueName The name of the residue to look for.
 * \param[in] atomName The name of the atom to look for.
 * \param[in] bExact Do we have the correct match.
 * \returns The index for the property.
 */
static int findPropertyIndex(AtomProperty*      ap,
                             ResidueType*       restype,
                             const std::string& residueName,
                             const std::string& atomName,
                             gmx_bool*          bExact)
{
    int j = NOTFOUND;

    bool bProtein  = restype->namedResidueHasType(residueName, "Protein");
    bool bProtWild = residueName == "AAA";
    int  malen     = NOTFOUND;
    int  mrlen     = NOTFOUND;
    for (size_t i = 0; (i < ap->entry.size()); i++)
    {
        int rlen = compareToDatabase(residueName, ap->entry[i].residueName);
        if (rlen == NOTFOUND)
        {
            if ((ap->entry[i].residueName == "*") || (ap->entry[i].residueName == "???"))
            {
                rlen = WILDCARD;
            }
            else if (ap->entry[i].residueName == "AAA")
            {
                rlen = WILDPROT;
            }
        }
        int alen = compareToDatabase(atomName, ap->entry[i].atomName);
        if ((alen > NOTFOUND) && (rlen > NOTFOUND))
        {
            if (((alen > malen) && (rlen >= mrlen)) || ((rlen > mrlen) && (alen >= malen)))
            {
                malen = alen;
                mrlen = rlen;
                j     = i;
            }
        }
    }

    *bExact = ((malen == static_cast<long int>(atomName.length()))
               && ((mrlen == static_cast<long int>(residueName.length())) || ((mrlen == WILDPROT) && bProtWild)
                   || ((mrlen == WILDCARD) && !bProtein && !bProtWild)));

    if (debug)
    {
        fprintf(debug, "searching residue: %4s atom: %4s\n", residueName.c_str(), atomName.c_str());
        if (j == NOTFOUND)
        {
            fprintf(debug, " not successful\n");
        }
        else
        {
            fprintf(debug, " match: %4s %4s\n", ap->entry[j].residueName.c_str(),
                    ap->entry[j].atomName.c_str());
        }
    }
    return j;
}

/*! \brief
 * Add new property to list.
 *
 * \param[in] ap Atomproperty to add.
 * \param[in] restype Residue type database to use.
 * \param[in] residueName Name of the residue.
 * \param[in] atomName Name of the atom.
 * \param[in] propValue Value of property.
 * \param[in] line Where to add property.
 */
static void addProperty(AtomProperty*      ap,
                        ResidueType*       restype,
                        const std::string& residueName,
                        const std::string& atomName,
                        real               propValue,
                        int                line)
{
    bool bExact;
    int  j = findPropertyIndex(ap, restype, residueName, atomName, &bExact);

    if (!bExact)
    {
        ap->entry.emplace_back(BaseEntry(atomName, residueName));

        j = ap->entry.size() - 1;
    }
    if (ap->entry[j].isAvailable)
    {
        if (ap->entry[j].value == propValue)
        {
            fprintf(stderr, "Warning double identical entries for %s %s %g on line %d in file %s\n",
                    residueName.c_str(), atomName.c_str(), propValue, line, ap->db.c_str());
        }
        else
        {
            fprintf(stderr,
                    "Warning double different entries %s %s %g and %g on line %d in file %s\n"
                    "Using last entry (%g)\n",
                    residueName.c_str(), atomName.c_str(), propValue, ap->entry[j].value, line,
                    ap->db.c_str(), propValue);
            ap->entry[j].value = propValue;
        }
    }
    else
    {
        ap->entry[j].isAvailable = TRUE;
        ap->entry[j].value       = propValue;
    }
}

/*! \brief
 * Read property value into structure.
 *
 * \param[in] ap Atomproperty to be read in.
 * \param[in] restype Library of residue types.
 * \param[in] factor Scaling factor for property.
 */
static void readProperty(AtomProperty* ap, ResidueType* restype, double factor)
{
    char line[STRLEN], resnm[32], atomnm[32];

    gmx::FilePtr fp      = gmx::openLibraryFile(ap->db);
    int          line_no = 0;
    while (get_a_line(fp.get(), line, STRLEN))
    {
        line_no++;
        double pp;
        if (sscanf(line, "%31s %31s %20lf", resnm, atomnm, &pp) == 3)
        {
            pp *= factor;
            addProperty(ap, restype, resnm, atomnm, pp, line_no);
        }
        else
        {
            fprintf(stderr, "WARNING: Error in file %s at line %d ignored\n", ap->db.c_str(), line_no);
        }
    }
    ap->isSet = TRUE;
}

/*! \brief
 * Set value for properties.
 *
 * \param[in] ap Atomproperty to set.
 * \param[in] restype Library of residue types.
 * \param[in] eprop Which property to set.
 * \param[in] haveBeenWarned If we already set a warning before
 * \returns True of warning should be printed.
 */
static bool setProperties(AtomProperty* ap, ResidueType* restype, int eprop, bool haveBeenWarned)
{
    const char* fns[epropNR] = { "atommass.dat", "vdwradii.dat", "dgsolv.dat", "electroneg.dat",
                                 "elements.dat" };
    double      fac[epropNR] = { 1.0, 1.0, 418.4, 1.0, 1.0 };
    double      def[epropNR] = { 12.011, 0.14, 0.0, 2.2, -1 };

    bool printWarning = false;
    if (!ap->isSet)
    {
        ap->db  = fns[eprop];
        ap->def = def[eprop];
        readProperty(ap, restype, fac[eprop]);

        if (debug)
        {
            fprintf(debug, "Entries in %s: %zu\n", ap->db.c_str(), ap->entry.size());
        }

        if ((!haveBeenWarned && (eprop == epropMass)) || (eprop == epropVDW))
        {
            printWarning = true;
        }
    }
    return printWarning;
}

AtomProperties::AtomProperties() : impl_(new Impl) {}

AtomProperties::~AtomProperties() {}

AtomProperty* AtomProperties::prop(int eprop)
{
    return &impl_->prop[eprop];
}

ResidueType* AtomProperties::restype()
{
    return &impl_->restype;
}

//! Print warning that vdW radii and masses are guessed.
static void printWarning()
{
    printf("\n"
           "WARNING: Masses and atomic (Van der Waals) radii will be guessed\n"
           "         based on residue and atom names, since they could not be\n"
           "         definitively assigned from the information in your input\n"
           "         files. These guessed numbers might deviate from the mass\n"
           "         and radius of the atom type. Please check the output\n"
           "         files if necessary.\n\n");
}

static void printvdwWarning(FILE* fp)
{
    if (nullptr != fp)
    {
        fprintf(fp, "NOTE: From version 5.0 %s uses the Van der Waals radii\n",
                gmx::getProgramContext().displayName());
        fprintf(fp, "from the source below. This means the results may be different\n");
        fprintf(fp, "compared to previous GROMACS versions.\n");
        please_cite(fp, "Bondi1964a");
    }
}

bool AtomProperties::setAtomProperty(int                eprop,
                                     const std::string& residueName,
                                     const std::string& atomName,
                                     real*              value)
{
    int         j;
    std::string tmpAtomName, tmpResidueName;
    gmx_bool    bExact;

    if (setProperties(prop(eprop), restype(), eprop, impl_->bWarned))
    {
        printWarning();
        impl_->bWarned = true;
    }
    if (isdigit(atomName[0]))
    {
        /* put digit after atomname */
        tmpAtomName.append(atomName.substr(1));
        tmpAtomName.append(1, atomName[0]);
    }
    else
    {
        tmpAtomName = atomName;
    }
    j = findPropertyIndex(&(impl_->prop[eprop]), &impl_->restype, residueName, tmpAtomName, &bExact);

    if (eprop == epropVDW && !impl_->bWarnVDW)
    {
        printvdwWarning(stdout);
        impl_->bWarnVDW = true;
    }
    if (j >= 0)
    {
        *value = impl_->prop[eprop].entry[j].value;
        return true;
    }
    else
    {
        *value = impl_->prop[eprop].def;
        return false;
    }
}


std::string AtomProperties::elementFromAtomNumber(int atomNumber)
{
    if (setProperties(prop(epropElement), restype(), epropElement, impl_->bWarned))
    {
        printWarning();
        impl_->bWarned = true;
    }
    for (const auto& e : prop(epropElement)->entry)
    {
        if (std::round(e.value) == atomNumber)
        {
            return e.atomName;
        }
    }
    return "";
}

int AtomProperties::atomNumberFromElement(const char* element)
{
    if (setProperties(prop(epropElement), restype(), epropElement, impl_->bWarned))
    {
        printWarning();
        impl_->bWarned = true;
    }
    for (const auto& e : prop(epropElement)->entry)
    {
        if (gmx_strcasecmp(e.atomName.c_str(), element) == 0)
        {
            return gmx::roundToInt(e.value);
        }
    }
    return -1;
}
