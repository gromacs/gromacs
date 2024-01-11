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
 * \brief Handling of selection parser symbol table.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SYMREC_H
#define GMX_SELECTION_SYMREC_H

#include <iterator>
#include <memory>
#include <string>

#include "external/boost/stl_interfaces/iterator_interface.hpp"

#include "selelem.h"

struct gmx_ana_selmethod_t;

namespace gmx
{

class SelectionParserSymbolTable;

/*! \internal
 * \brief
 * Single symbol for the selection parser.
 *
 * Public methods in this class do not throw.
 *
 * \ingroup module_selection
 */
class SelectionParserSymbol
{
public:
    //! Defines the type of the symbol.
    enum SymbolType
    {
        ReservedSymbol, //!< The symbol is a reserved keyword.
        VariableSymbol, //!< The symbol is a variable.
        MethodSymbol,   //!< The symbol is a selection method.
        PositionSymbol  //!< The symbol is a position keyword.
    };

    ~SelectionParserSymbol();

    //! Returns the name of the symbol.
    const std::string& name() const;
    //! Returns the type of the symbol.
    SymbolType type() const;

    /*! \brief
     * Returns the method associated with a \ref MethodSymbol symbol.
     *
     * \returns   The method associated with the symbol.
     *
     * Must only be called if type() returns \ref MethodSymbol.
     */
    gmx_ana_selmethod_t* methodValue() const;
    /*! \brief
     * Returns the selection tree associated with a \ref VariableSymbol symbol.
     *
     * \returns   The variable expression associated with the symbol.
     *
     * Must only be called if type() returns \ref VariableSymbol.
     */
    const SelectionTreeElementPointer& variableValue() const;

private:
    class Impl;

    /*! \brief
     * Initializes a new symbol with the given data.
     *
     * \param  impl  Implementation data.
     * \throws std::bad_alloc if out of memory.
     *
     * Only the parent symbol table creates symbol objects.
     */
    explicit SelectionParserSymbol(Impl* impl);

    std::unique_ptr<Impl> impl_;

    /*! \brief
     * Needed to call the constructor and for other initialization.
     */
    friend class SelectionParserSymbolTable;
};

/*! \internal
 * \brief
 * Forward iterator for iterating symbols of a given type.
 *
 * Behaves as standard C++ forward iterator.  To get an iterator, call
 * SelectionParserSymbolTable::beginIterator().  Each time the iterator is
 * incremented, it moves to the next symbol of the type given when the iterator
 * was created.  When there are no more symbols, the iterator will equal
 * SelectionParserSymbolTable::endIterator().  It is not allowed to dereference
 * or increment an iterator that has reached the end.
 *
 * Construction and assignment may throw std::bad_alloc if out of memory.
 * Other methods do not throw.
 *
 * \see SelectionParserSymbolTable::beginIterator()
 *
 * \ingroup module_selection
 */
class SelectionParserSymbolIterator :
    public gmx::boost::stl_interfaces::iterator_interface<SelectionParserSymbolIterator, std::forward_iterator_tag, const SelectionParserSymbol>
{
    using Base =
            gmx::boost::stl_interfaces::iterator_interface<SelectionParserSymbolIterator, std::forward_iterator_tag, const SelectionParserSymbol>;

public:
    //! Creates an independent copy of an iterator.
    SelectionParserSymbolIterator(const SelectionParserSymbolIterator& other);
    ~SelectionParserSymbolIterator();

    //! Creates an independent copy of an iterator.
    SelectionParserSymbolIterator& operator=(const SelectionParserSymbolIterator& other);

    //! Equality comparison for iterators.
    bool operator==(const SelectionParserSymbolIterator& other) const;
    //! Dereferences the iterator.
    reference operator*() const;
    //! Moves the iterator to the next symbol.
    SelectionParserSymbolIterator& operator++();
    using Base::                   operator++;

private:
    class Impl;

    /*! \brief
     * Initializes a new iterator with the given data.
     *
     * \param  impl  Implementation data.
     *
     * Only the parent symbol table can create non-default-constructed
     * iterators.
     */
    explicit SelectionParserSymbolIterator(Impl* impl);

    std::unique_ptr<Impl> impl_;

    /*! \brief
     * Needed to access the constructor.
     */
    friend class SelectionParserSymbolTable;
};

/*! \internal \brief
 * Symbol table for the selection parser.
 *
 * \ingroup module_selection
 */
class SelectionParserSymbolTable
{
public:
    /*! \brief
     * Creates a new symbol table.
     *
     * \throws std::bad_alloc if out of memory.
     *
     * The created table is initialized with reserved and position symbols.
     */
    SelectionParserSymbolTable();
    ~SelectionParserSymbolTable();

    /*! \brief
     * Finds a symbol by name.
     *
     * \param[in] name   Symbol name to find.
     * \returns   Pointer to the symbol with name \p name, or
     *      NULL if not found.
     *
     * Does not throw.
     */
    const SelectionParserSymbol* findSymbol(const std::string& name) const;

    /*! \brief
     * Returns the start iterator for iterating symbols of a given type.
     *
     * \param[in] type  Type of symbols to iterate over.
     * \returns   Iterator that points to the first symbol of type \p type.
     * \throws    std::bad_alloc if out of memory.
     *
     * \see SelectionParserSymbolIterator
     */
    SelectionParserSymbolIterator beginIterator(SelectionParserSymbol::SymbolType type) const;
    /*! \brief
     * Returns the end iterator for symbol iteration.
     *
     * \throws    std::bad_alloc if out of memory.
     *
     * Currently, the end value is the same for all symbol types.
     *
     * \see SelectionParserSymbolIterator
     */
    SelectionParserSymbolIterator endIterator() const;

    /*! \brief
     * Adds a new variable symbol.
     *
     * \param[in] name   Name of the new symbol.
     * \param[in] sel    Value of the variable.
     * \throws    std::bad_alloc if out of memory.
     * \throws    InvalidInputError if there was a symbol with the same
     *      name.
     */
    void addVariable(const char* name, const SelectionTreeElementPointer& sel);
    /*! \brief
     * Adds a new method symbol.
     *
     * \param[in] name   Name of the new symbol.
     * \param[in] method Method that this symbol represents.
     * \throws    std::bad_alloc if out of memory.
     * \throws    APIError if there was a symbol with the same name.
     */
    void addMethod(const char* name, gmx_ana_selmethod_t* method);

private:
    class Impl;

    std::unique_ptr<Impl> impl_;

    /*! \brief
     * Needed to access implementation types.
     */
    friend class SelectionParserSymbolIterator;
};

} // namespace gmx

#endif
