/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares gmx::SelectionCompiler.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_COMPILER_H
#define GMX_SELECTION_COMPILER_H

namespace gmx
{

class SelectionCollection;

/*! \internal
 * \brief
 * Implements selection compilation.
 *
 * This function is used to implement SelectionCollection::compile().
 * It prepares the selections in a selection collection for evaluation and
 * performs some optimizations.
 *
 * Before compilation, the selection collection should have been initialized
 * with gmx_ana_selcollection_parse_*().
 * The compiled selection collection can be passed to
 * gmx_ana_selcollection_evaluate() to evaluate the selection for a frame.
 * If an error occurs, \p coll is cleared.
 *
 * The covered fraction information in \p coll is initialized to
 * \ref CFRAC_NONE.
 *
 * See \ref page_module_selection_compiler.
 *
 * \param[in, out] coll Selection collection to work on.
 *
 * \ingroup module_selection
 */
void compileSelection(SelectionCollection* coll);

} // namespace gmx

#endif
