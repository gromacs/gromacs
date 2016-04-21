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
/*! \defgroup module_alexandria Processing of input files and force fields
 * \ingroup group_preprocessing
 * \brief
 * Provides tools for input processing based on the alexandria force field
 *
 * The tools in the alexandria directory are under heavy development still.
 * Once finished they will provide processing of small molecules files based
 * on quantum chemistry or other input source through the OpenBabel library.
 * Assigning of atom types, derivation of charges, and generation of
 * molecular topology files (or objects) is implemented.
 *
 * \author David van der Spoel <david.vanderspoel@gmail.com>
 * \inpublicapi
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
