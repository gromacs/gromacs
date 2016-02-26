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
 * \ingroup module_alexandria
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
        ~CategoryListElement() {};

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

        ~CategoryList() {};

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
