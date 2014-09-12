/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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
