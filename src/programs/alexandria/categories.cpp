/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
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
