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
#include "gmxpre.h"

#include "composition.h"

namespace alexandria
{

CompositionSpecs::CompositionSpecs()
{
    cs_.push_back(CompositionSpec(iCalexandria, (const char *)"alexandria", (const char *)"Maaren2014a", (const char *)"AX"));
    cs_.push_back(CompositionSpec(iCbosque, (const char *)"bosque", (const char *)"Bosque2002a", (const char *)"BS"));
    cs_.push_back(CompositionSpec(iCmiller, (const char *)"miller", (const char *)"Miller1990a", (const char *)"MK"));
}

}
