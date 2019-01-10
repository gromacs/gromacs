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

#include "gromacs/compat/make_unique.h"
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
enum {
    NOTFOUND = -4, WILDCARD, WILDPROT, PROTEIN
};

//! Conglomeration of atom property entries.
struct AtomProperty {
    //! Has property been set.
    bool   isSet = false;
    //! Number of properties.
    int    nprop = 0;
    //! Max number of properties.
    int    maxprop = 0;
    //! Database the property is coming from.
    char  *db = nullptr;
    //! Default value for property.
    double def = 0.0;
    //! Array of names for atoms.
    char **atomName = nullptr;
    //! Array of names for residues.
    char **residueName = nullptr;
    //! Array of flags if property is available.
    bool  *isAvailable = nullptr;
    //! Array of values for property.
    real  *value = nullptr;
};

//! Implementation detail type for Atomproperties.
class AtomProperties::Impl
{
    public:
        //! Should user be warned about error.
        bool               bWarned = false;
        //! Should user be warned about vdW not found.
        bool               bWarnVDW = false;
        //! The different atom properties.
        AtomProperty       prop[epropNR];
        //! The residue types.
        gmx_residuetype_t *restype = nullptr;
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
static int findNumberOfMatches(const char *search, const char *database)
{
    int i = 0;
    while (search[i] && database[i] && (search[i] == database[i]) )
    {
        i++;
    }

    if (database[i])
    {
        i = NOTFOUND;
    }
    return i;
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
static int findPropertyIndex(AtomProperty *ap, gmx_residuetype_t *restype,
                             char *residueName, char *atomName,
                             gmx_bool *bExact)
{
    int      j = NOTFOUND;

    bool     bProtein  = gmx_residuetype_is_protein(restype, residueName);
    bool     bProtWild = (strcmp(residueName, "AAA") == 0);
    int      malen     = NOTFOUND;
    int      mrlen     = NOTFOUND;
    for (int i = 0; (i < ap->nprop); i++)
    {
        int rlen = findNumberOfMatches(residueName, ap->residueName[i]);
        if (rlen == NOTFOUND)
        {
            if ( (strcmp(ap->residueName[i], "*") == 0) ||
                 (strcmp(ap->residueName[i], "???") == 0) )
            {
                rlen = WILDCARD;
            }
            else if (strcmp(ap->residueName[i], "AAA") == 0)
            {
                rlen = WILDPROT;
            }
        }
        int alen = findNumberOfMatches(atomName, ap->atomName[i]);
        if ( (alen > NOTFOUND) && (rlen > NOTFOUND))
        {
            if ( ( (alen > malen) && (rlen >= mrlen)) ||
                 ( (rlen > mrlen) && (alen >= malen) ) )
            {
                malen = alen;
                mrlen = rlen;
                j     = i;
            }
        }
    }

    *bExact = ((malen == static_cast<long int>(strlen(atomName))) &&
               ((mrlen == static_cast<long int>(strlen(residueName))) ||
                ((mrlen == WILDPROT) && bProtWild) ||
                ((mrlen == WILDCARD) && !bProtein && !bProtWild)));

    if (debug)
    {
        fprintf(debug, "searching residue: %4s atom: %4s\n", residueName, atomName);
        if (j == NOTFOUND)
        {
            fprintf(debug, " not successful\n");
        }
        else
        {
            fprintf(debug, " match: %4s %4s\n", ap->residueName[j], ap->atomName[j]);
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
static void addProperty(AtomProperty *ap, gmx_residuetype_t *restype,
                        char *residueName, char *atomName,
                        real propValue, int line)
{
    bool bExact;
    int  j = findPropertyIndex(ap, restype, residueName, atomName, &bExact);

    if (!bExact)
    {
        if (ap->nprop >= ap->maxprop)
        {
            ap->maxprop += 10;
            srenew(ap->residueName, ap->maxprop);
            srenew(ap->atomName, ap->maxprop);
            srenew(ap->value, ap->maxprop);
            srenew(ap->isAvailable, ap->maxprop);
            for (int i = ap->nprop; (i < ap->maxprop); i++)
            {
                ap->atomName[i]     = nullptr;
                ap->residueName[i]  = nullptr;
                ap->value[i]        = 0;
                ap->isAvailable[i]  = FALSE;
            }
        }
        ap->atomName[ap->nprop]     = gmx_strdup(atomName);
        ap->residueName[ap->nprop]  = gmx_strdup(residueName);
        j                           = ap->nprop;
        ap->nprop++;
    }
    if (ap->isAvailable[j])
    {
        if (ap->value[j] == propValue)
        {
            fprintf(stderr, "Warning double identical entries for %s %s %g on line %d in file %s\n",
                    residueName, atomName, propValue, line, ap->db);
        }
        else
        {
            fprintf(stderr, "Warning double different entries %s %s %g and %g on line %d in file %s\n"
                    "Using last entry (%g)\n",
                    residueName, atomName, propValue, ap->value[j], line, ap->db, propValue);
            ap->value[j] = propValue;
        }
    }
    else
    {
        ap->isAvailable[j] = TRUE;
        ap->value[j]       = propValue;
    }
}

/*! \brief
 * Read property value into structure.
 *
 * \param[in] ap Atomproperty to be read in.
 * \param[in] restype Library of residue types.
 * \param[in] factor Scaling factor for property.
 */
static void readProperty(AtomProperty *ap, gmx_residuetype_t *restype, double factor)
{
    char          line[STRLEN], resnm[32], atomnm[32];

    gmx::FilePtr  fp      = gmx::openLibraryFile(ap->db);
    int           line_no = 0;
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
            fprintf(stderr, "WARNING: Error in file %s at line %d ignored\n",
                    ap->db, line_no);
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
static bool setProperties(AtomProperty *ap, gmx_residuetype_t *restype, int eprop, bool haveBeenWarned)
{
    const char       *fns[epropNR]  = { "atommass.dat", "vdwradii.dat", "dgsolv.dat", "electroneg.dat", "elements.dat" };
    double            fac[epropNR]  = { 1.0,    1.0,  418.4, 1.0, 1.0 };
    double            def[epropNR]  = { 12.011, 0.14, 0.0, 2.2, -1 };

    bool              printWarning = false;
    if (!ap->isSet)
    {
        ap->db  = gmx_strdup(fns[eprop]);
        ap->def = def[eprop];
        readProperty(ap, restype, fac[eprop]);

        if (debug)
        {
            fprintf(debug, "Entries in %s: %d\n", ap->db, ap->nprop);
        }

        if ( (!haveBeenWarned && (eprop == epropMass) ) || (eprop == epropVDW))
        {
            printWarning = true;
        }

    }
    return printWarning;
}

AtomProperties::AtomProperties()
    : impl_(new Impl)
{
    gmx_residuetype_init(&impl_->restype);
}

static void destroy_prop(AtomProperty *ap)
{
    if (ap->isSet)
    {
        sfree(ap->db);

        for (int i = 0; i < ap->nprop; i++)
        {
            sfree(ap->atomName[i]);
            sfree(ap->residueName[i]);
        }
        sfree(ap->atomName);
        sfree(ap->residueName);
        sfree(ap->isAvailable);
        sfree(ap->value);
    }
}

AtomProperties::~AtomProperties()
{
    for (int p = 0; p < epropNR; p++)
    {
        destroy_prop(&impl_->prop[p]);
    }
    gmx_residuetype_destroy(impl_->restype);
}

AtomProperty *AtomProperties::prop(int eprop)
{
    return &impl_->prop[eprop];
}

gmx_residuetype_t *AtomProperties::restype()
{
    return impl_->restype;
}

void AtomProperties::printWarning()
{
    printf("\n"
           "WARNING: Masses and atomic (Van der Waals) radii will be guessed\n"
           "         based on residue and atom names, since they could not be\n"
           "         definitively assigned from the information in your input\n"
           "         files. These guessed numbers might deviate from the mass\n"
           "         and radius of the atom type. Please check the output\n"
           "         files if necessary.\n\n");
    impl_->bWarned = true;
}

static void printvdwWarning(FILE *fp)
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

bool AtomProperties::setAtomProperty(int         eprop,
                                     const char *residueName,
                                     const char *atomName,
                                     real       *value)
{
    int           j;
#define MAXQ 32
    char          atomname[MAXQ], resname[MAXQ];
    gmx_bool      bExact;

    if (setProperties(prop(eprop), restype(), eprop, impl_->bWarned))
    {
        printWarning();
    }
    if ((strlen(atomName) > MAXQ-1) || (strlen(residueName) > MAXQ-1))
    {
        if (debug)
        {
            fprintf(debug, "WARNING: will only compare first %d characters\n",
                    MAXQ-1);
        }
    }
    if (isdigit(atomName[0]))
    {
        int i;
        /* put digit after atomname */
        for (i = 1; i < MAXQ-1 && atomName[i] != '\0'; i++)
        {
            atomname[i-1] = atomName[i];
        }
        atomname[i-1] = atomName[0];
        atomname[i]   = '\0';
    }
    else
    {
        strncpy(atomname, atomName, MAXQ-1);
    }
    strncpy(resname, residueName, MAXQ-1);

    j = findPropertyIndex(&(impl_->prop[eprop]), impl_->restype, resname,
                          atomname, &bExact);

    if (eprop == epropVDW && !impl_->bWarnVDW)
    {
        printvdwWarning(stdout);
        impl_->bWarnVDW = TRUE;
    }
    if (j >= 0)
    {
        *value = impl_->prop[eprop].value[j];
        return TRUE;
    }
    else
    {
        *value = impl_->prop[eprop].def;
        return FALSE;
    }
}


const char *AtomProperties::elementFromAtomNumber(int atomNumber)
{
    if (setProperties(prop(epropElement), restype(), epropElement, impl_->bWarned))
    {
        printWarning();
    }
    for (int i = 0; (i < impl_->prop[epropElement].nprop); i++)
    {
        if (std::round(impl_->prop[epropElement].value[i]) == atomNumber)
        {
            return impl_->prop[epropElement].atomName[i];
        }
    }
    return nullptr;
}

int AtomProperties::atomNumberFromElement(const char *element)
{
    if (setProperties(prop(epropElement), restype(), epropElement, impl_->bWarned))
    {
        printWarning();
    }
    for (int i = 0; (i < impl_->prop[epropElement].nprop); i++)
    {
        if (gmx_strcasecmp(impl_->prop[epropElement].atomName[i], element) == 0)
        {
            return gmx::roundToInt(impl_->prop[epropElement].value[i]);
        }
    }
    return -1;
}
