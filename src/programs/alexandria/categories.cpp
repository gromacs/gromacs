/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "categories.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/futil.h"

#include "composition.h"
#include "molprop.h"
#include "molprop_tables.h"
#include "molprop_util.h"
#include "molselect.h"
#include "poldata.h"
#include "poldata_xml.h"

namespace alexandria
{

static bool CompareCategoryListElements(CategoryListElement ca,
                                        CategoryListElement cb)
{
    return (ca.getName().compare(cb.getName()) < 0);
}

static bool CompareStrings(std::string ca, std::string cb)
{
    return (ca.compare(cb) < 0);
}

void CategoryListElement::addMolecule(std::string molecule)
{
    for (std::vector<std::string>::iterator i = beginMolecules(); (i < endMolecules()); ++i)
    {
        if (i->compare(molecule) == 0)
        {
            return;
        }
    }
    molecule_.push_back(molecule);
}

void CategoryListElement::sortMolecules()
{
    std::sort(beginMolecules(), endMolecules(), CompareStrings);
}

bool CategoryListElement::hasMolecule(std::string molecule)
{
    return std::find(beginMolecules(), endMolecules(), molecule) != endMolecules();
}

void CategoryList::addCategory(std::string catname, std::string molecule)
{
    CategoryListElementIterator i;

    for (i = catli_.begin(); (i < catli_.end()); ++i)
    {
        if (i->getName().compare(catname) == 0)
        {
            i->addMolecule(molecule);
            break;
        }
    }
    if (i == catli_.end())
    {
        catli_.push_back(CategoryListElement(catname, molecule));
    }
}

void CategoryList::sortCategories()
{
    std::sort(beginCategories(), endCategories(), CompareCategoryListElements);
    for (CategoryListElementIterator i = beginCategories(); (i < endCategories()); ++i)
    {
        i->sortMolecules();
    }
}

void makeCategoryList(CategoryList         &cList,
                      std::vector<MolProp>  mp,
                      const MolSelect      &gms,
                      iMolSelect            ims)
{
    alexandria::CompositionSpecs  cs;
    const char                   *alex = cs.searchCS(alexandria::iCalexandria)->name();

    for (std::vector<alexandria::MolProp>::iterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if ((ims == gms.status(mpi->getIupac())) &&
            mpi->HasComposition(alex))
        {
            for (std::vector<std::string>::iterator si = mpi->BeginCategory(); (si < mpi->EndCategory()); si++)
            {
                cList.addCategory(*si, mpi->getIupac());
            }
        }
    }
    cList.sortCategories();
}

}
