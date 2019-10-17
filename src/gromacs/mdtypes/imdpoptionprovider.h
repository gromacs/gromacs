/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::IMdpOptionProvider.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IMDPOPTIONPROVIDER_H
#define GMX_MDTYPES_IMDPOPTIONPROVIDER_H

namespace gmx
{

class IKeyValueTreeTransformRules;
class IOptionsContainerWithSections;
class KeyValueTreeObjectBuilder;

/*! \libinternal \brief
 * Interface for handling mdp/tpr input to a mdrun module.
 *
 * This interface provides a mechanism for modules to contribute data that
 * traditionally has been kept in t_inputrec.  This is essentially parameters
 * read from an mdp file and subsequently stored in a tpr file.
 *
 * The main method to implement is initMdpOptions().  This declares a set of
 * options in nested sections.  When declaring options, the module defines its
 * own internal variables where the values will be stored (see \ref
 * module_options for an overview).  The section structure for the options also
 * defines an internal data representation in the form of a generic key-value
 * tree.  This is used to store the values in the tpr file, and also for
 * broadcasting, gmx dump, and gmx check.
 *
 * Implementation of initMdpTransform() is required for populating the options
 * from the current flat mdp format.  It specifies how the mdp parameters map
 * into the options structure/internal key-value tree specified in
 * initMdpOptions().
 *
 * See \ref page_mdmodules for more details on how mdp parsing and related
 * functionality works with this interface.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IMdpOptionProvider
{
public:
    /*! \brief
     * Initializes a transform from mdp values to sectioned options.
     *
     * The transform is specified from a flat KeyValueTreeObject that
     * contains each mdp value as a property, to a structured
     * key-value tree that should match the options defined in
     * initMdpOptions().
     *
     * This method may be removed once the flat mdp file is replaced with a
     * more structure input file (that can be directly read into the
     * internal key-value tree), and there is no longer any need for
     * backward compatibility with old files.
     */
    virtual void initMdpTransform(IKeyValueTreeTransformRules* transform) = 0;
    /*! \brief
     * Initializes options that declare input (mdp) parameters for this
     * module.
     */
    virtual void initMdpOptions(IOptionsContainerWithSections* options) = 0;
    //! Prepares to write a flat key-value tree like an mdp file.
    virtual void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const = 0;

protected:
    virtual ~IMdpOptionProvider() {}
};

} // namespace gmx

#endif
