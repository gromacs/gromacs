/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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


#ifndef CATEGORIES_H
#define CATEGORIES_H

#include <string.h>

#include <string>
#include <vector>

#include "gromacs/utility/real.h"

#include "molselect.h"

/*! \brief
 * Contains classes related to alexandria force field tools
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
namespace alexandria
{
class MolProp;

class CategoryListElement
{
    private:
        //! The name of the category
        std::string              cat_;
        //! Molecules in this category
        std::vector<std::string> molecule_;
    public:
        CategoryListElement(std::string cat, std::string molecule)
        {
            cat_   = cat;
            addMolecule(molecule);
        }

        void addMolecule(std::string molecule);

        int nMolecule() { return molecule_.size(); }

        bool hasMolecule(std::string molecule);

        std::string getName() { return cat_; }

        void sortMolecules();

        std::vector<std::string>::iterator beginMolecules() { return molecule_.begin(); }

        std::vector<std::string>::iterator endMolecules() { return molecule_.end(); }
};

typedef std::vector<CategoryListElement>::iterator CategoryListElementIterator;

class CategoryList
{
    private:
        std::vector<CategoryListElement> catli_;
    public:
        CategoryList() {};

        void addCategory(std::string catname, std::string molecule);

        void sortCategories();

        int nCategories() { return catli_.size(); }

        CategoryListElementIterator beginCategories() { return catli_.begin(); }

        CategoryListElementIterator endCategories() { return catli_.end(); }
};

void makeCategoryList(CategoryList         &cList,
                      std::vector<MolProp>  mp,
                      const MolSelect      &gms,
                      iMolSelect            ims);

}

#endif
