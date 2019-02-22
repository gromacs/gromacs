/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_GPP_ATOMTYPE_H
#define GMX_GMXPREPROCESS_GPP_ATOMTYPE_H

#include <cstdio>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_atom;
struct t_atomtypes;
class InteractionType;
struct InteractionTypeParameters;
struct t_symtab;

/*! \libinternal \brief
 * Storage of all atom types used during preprocessing of a simulation
 * input.
 */
class PreprocessingAtomTypes
{
    public:
        PreprocessingAtomTypes();
        //! Move constructor.
        PreprocessingAtomTypes(PreprocessingAtomTypes &&old) noexcept;
        //! Move assignment constructor.
        PreprocessingAtomTypes &operator=(PreprocessingAtomTypes &&old) noexcept;

        ~PreprocessingAtomTypes();

        /*! \brief
         *  Get atomtype from string input.
         *
         *  \param[in] str Input string to search type for.
         *  \returns Atomtype as integer.
         */
        int atomTypeFromString(const std::string &str) const;

        //! Get number of defined atom types.
        size_t size() const;

        /*! \brief
         * Get name of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The type name.
         */
        const char *atomNameFromType(int nt) const;

        /*! \brief
         * Get normal mass of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The mass for the atom in state A or NOTSET.
         */
        real atomMassAFromType(int nt) const;

        /*! \brief
         * Get mass for B state of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The mass for the atom in state B or NOTSET.
         */
        real atomMassBFromType(int nt) const;

        /*! \brief
         * Get normal charge of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The charge for the atom in state A or NOTSET.
         */
        real atomChargeAFromType(int nt) const;

        /*! \brief
         * Get charge for B state of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The charge for the atom in state B or NOTSET.
         */
        real atomChargeBFromType(int nt) const;

        /*! \brief
         * Get normal parameter type of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The parameter type or NOTSET.
         */
        int atomPTypeFromType(int nt) const;

        /*! \brief
         * Get bond atom parameter of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The bond atom parameter or NOTSET.
         */
        int bondAtomParameterFromType(int nt) const;

        /*! \brief
         * Get atomic number of atom from internal atom type number.
         *
         * \param[in] nt Internal number of atom type.
         * \returns The atomic number type or NOTSET.
         */
        int atomNumberFromType(int nt) const;

        /*! \brief
         * Get the value of \p param of type \p nt.
         *
         * \param[in] param The parameter value to find.
         * \param[in] nt The number of the type.
         * \returns The value of the parameter or NOTSET.
         */
        real atomNonBondedParamFromType(int nt, int param) const;

        /*! \brief
         * If a value is within the range of the current types or not.
         *
         * \param[in] nt Value to check.
         * \returns True if value is in range.
         */
        bool isSet(int nt) const;

        /*! \brief
         * Print data to file.
         *
         * \param[in] out File pointer.
         */
        void printTypes(FILE *out);

        /*! \brief
         * Set the values of an existing atom type \p nt.
         *
         * \param[in] nt Type that should be set.
         * \param[in] tab Symbol table.
         * \param[in] a Atom information.
         * \param[in] name Atom name.
         * \param[in] nb Nonbonded parameters.
         * \param[in] bondAtomType What kind of bonded interactions are there.
         * \param[in] atomNumber Atomic number of the entry.
         * \returns Number of the type set or NOTSET
         */
        int setType(int                    nt,
                    t_symtab              *tab,
                    const t_atom          &a,
                    const char            *name,
                    const InteractionType &nb,
                    int                    bondAtomType,
                    int                    atomNumber);

        /*! \brief
         * Add new atom type to database.
         *
         * \param[in] tab Symbol table.
         * \param[in] a Atom information.
         * \param[in] name Atom name.
         * \param[in] nb Nonbonded parameters.
         * \param[in] bondAtomType What kind of bonded interactions are there.
         * \param[in] atomNumber Atomic number of the entry.
         * \returns Number of entries in database.
         */
        int addType(t_symtab              *tab,
                    const t_atom          &a,
                    const char            *name,
                    const InteractionType &nb,
                    int                    bondAtomType,
                    int                    atomNumber);

        /*! \brief
         * Renumber existing atom type entries.
         *
         * \param[in] plist List of parameters.
         * \param[in] mtop Global topology.
         * \param[inout] wallAtomType Atom types of wall atoms, which may also be renumbered
         * \param[in] verbose If we want to print additional info.
         */
        void renumberTypes(gmx::ArrayRef<InteractionTypeParameters> plist,
                           gmx_mtop_t                              *mtop,
                           int                                     *wallAtomType,
                           bool                                     verbose);

        /*! \brief
         * Copy information to other structure.
         *
         * \param[inout] atypes Other datastructure to copy to.
         */
        void copyTot_atomtypes(t_atomtypes *atypes) const;
    private:
        class Impl;
        gmx::PrivateImplPointer<Impl> impl_;
};

#endif
