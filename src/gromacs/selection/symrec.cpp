/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements classes in symrec.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "symrec.h"

#include <map>
#include <string>
#include <utility>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/uniqueptr.h"

#include "poscalc.h"
#include "selelem.h"

namespace gmx
{

/********************************************************************
 * SelectionParserSymbol
 */

/*! \internal \brief
 * Private implementation class for SelectionParserSymbol.
 *
 * \ingroup module_selection
 */
class SelectionParserSymbol::Impl
{
    public:
        /*! \brief
         * Initializes a symbol.
         *
         * \param[in] type  Type for the symbol.
         * \param[in] name  Name for the symbol.
         *
         * The symbol table is responsible for initializing the \a meth_ and
         * \a var_ members as appropriate.
         */
        Impl(SymbolType type, const char *name)
            : name_(name), type_(type), meth_(NULL)
        {
        }

        //! Name of the symbol.
        std::string                     name_;
        //! Type of the symbol.
        SymbolType                      type_;
        //! Pointer to the method structure (\ref MethodSymbol).
        gmx_ana_selmethod_t            *meth_;
        //! Pointer to the variable value (\ref VariableSymbol).
        SelectionTreeElementPointer     var_;
};

SelectionParserSymbol::SelectionParserSymbol(Impl *impl)
    : impl_(impl)
{
}

SelectionParserSymbol::~SelectionParserSymbol()
{
}

const std::string &
SelectionParserSymbol::name() const
{
    return impl_->name_;
}

SelectionParserSymbol::SymbolType
SelectionParserSymbol::type() const
{
    return impl_->type_;
}

gmx_ana_selmethod_t *
SelectionParserSymbol::methodValue() const
{
    GMX_RELEASE_ASSERT(type() == MethodSymbol,
                       "Attempting to get method handle for a non-method symbol");
    return impl_->meth_;
}

const gmx::SelectionTreeElementPointer &
SelectionParserSymbol::variableValue() const
{
    GMX_RELEASE_ASSERT(type() == VariableSymbol,
                       "Attempting to get variable value for a non-variable symbol");
    return impl_->var_;
}

/********************************************************************
 * SelectionParserSymbolTable::Impl
 */

/*! \internal
 * \brief
 * Private implementation class for SelectionParserSymbolTable.
 *
 * All methods in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_selection
 */
class SelectionParserSymbolTable::Impl
{
    public:
        //! Smart pointer type for managing a SelectionParserSymbol.
        typedef gmx::gmx_unique_ptr<SelectionParserSymbol>::type
            SymbolPointer;
        //! Container type for the list of symbols.
        typedef std::map<std::string, SymbolPointer> SymbolMap;

        /*! \brief
         * Adds a symbol to the symbol list.
         *
         * \param[in] symbol  Symbol to add.
         */
        void addSymbol(SymbolPointer symbol);
        //! Adds the reserved symbols to this symbol table.
        void addReservedSymbols();
        //! Adds the position symbols to this symbol table.
        void addPositionSymbols();

        //! Symbols in this symbol table.
        SymbolMap               symbols_;
};

void
SelectionParserSymbolTable::Impl::addSymbol(SymbolPointer symbol)
{
    symbols_.insert(std::make_pair(symbol->name(), move(symbol)));
}

void
SelectionParserSymbolTable::Impl::addReservedSymbols()
{
    const char *const sym_reserved[] = {
        "group",
        "to",
        "not",
        "and",
        "or",
        "xor",
        "yes",
        "no",
        "on",
        "off"
    };

    for (size_t i = 0; i < asize(sym_reserved); ++i)
    {
        SymbolPointer sym(new SelectionParserSymbol(
                                  new SelectionParserSymbol::Impl(
                                          SelectionParserSymbol::ReservedSymbol, sym_reserved[i])));
        addSymbol(move(sym));
    }
}

void
SelectionParserSymbolTable::Impl::addPositionSymbols()
{
    const char *const *postypes
        = gmx::PositionCalculationCollection::typeEnumValues;
    for (int i = 0; postypes[i] != NULL; ++i)
    {
        SymbolPointer sym(new SelectionParserSymbol(
                                  new SelectionParserSymbol::Impl(
                                          SelectionParserSymbol::PositionSymbol, postypes[i])));
        addSymbol(move(sym));
    }
}

/********************************************************************
 * SelectionParserSymbolIterator
 */

/*! \internal \brief
 * Private implementation class for SelectionParserSymbolIterator.
 *
 * \ingroup module_selection
 */
class SelectionParserSymbolIterator::Impl
{
    public:
        //! Shorthand for the underlying iterator type.
        typedef SelectionParserSymbolTable::Impl::SymbolMap::const_iterator
            IteratorType;

        /*! \brief
         * Constructs an end iterator.
         *
         * \param[in] end  Iterator to the end of the iterated container.
         */
        explicit Impl(IteratorType end)
            : iter_(end), end_(end)
        {
        }
        /*! \brief
         * Constructs an iterator.
         *
         * \param[in] iter Iterator to the current symbol.
         * \param[in] end  Iterator to the end of the iterated container.
         */
        Impl(IteratorType iter, IteratorType end)
            : iter_(iter), end_(end)
        {
        }

        //! Underlying iterator to the symbol container.
        IteratorType            iter_;
        //! End of the symbol container being iterated.
        IteratorType            end_;
};

SelectionParserSymbolIterator::SelectionParserSymbolIterator(Impl *impl)
    : impl_(impl)
{
}

SelectionParserSymbolIterator::SelectionParserSymbolIterator(
        const SelectionParserSymbolIterator &other)
    : impl_(new Impl(*other.impl_))
{
}

SelectionParserSymbolIterator::~SelectionParserSymbolIterator()
{
}

SelectionParserSymbolIterator &SelectionParserSymbolIterator::operator=(
        const SelectionParserSymbolIterator &other)
{
    impl_.reset(new Impl(*other.impl_));
    return *this;
}

bool SelectionParserSymbolIterator::operator==(
        const SelectionParserSymbolIterator &other) const
{
    return impl_->iter_ == other.impl_->iter_;
}

const SelectionParserSymbol &SelectionParserSymbolIterator::operator*() const
{
    return *impl_->iter_->second;
}

SelectionParserSymbolIterator &SelectionParserSymbolIterator::operator++()
{
    SelectionParserSymbol::SymbolType type = impl_->iter_->second->type();
    do
    {
        ++impl_->iter_;
    }
    while (impl_->iter_ != impl_->end_ && impl_->iter_->second->type() != type);
    return *this;
}

/********************************************************************
 * SelectionParserSymbolTable
 */

SelectionParserSymbolTable::SelectionParserSymbolTable()
    : impl_(new Impl)
{
    impl_->addReservedSymbols();
    impl_->addPositionSymbols();
}

SelectionParserSymbolTable::~SelectionParserSymbolTable()
{
}

const SelectionParserSymbol *
SelectionParserSymbolTable::findSymbol(const std::string &name) const
{
    Impl::SymbolMap::const_iterator sym = impl_->symbols_.lower_bound(name);
    if (sym == impl_->symbols_.end())
    {
        return NULL;
    }
    if (sym->second->name() == name)
    {
        return sym->second.get();
    }
    return NULL;
}

SelectionParserSymbolIterator
SelectionParserSymbolTable::beginIterator(SelectionParserSymbol::SymbolType type) const
{
    Impl::SymbolMap::const_iterator sym;
    Impl::SymbolMap::const_iterator end = impl_->symbols_.end();
    for (sym = impl_->symbols_.begin(); sym != end; ++sym)
    {
        if (sym->second->type() == type)
        {
            return SelectionParserSymbolIterator(
                    new SelectionParserSymbolIterator::Impl(sym, end));
        }
    }
    return endIterator();
}

SelectionParserSymbolIterator
SelectionParserSymbolTable::endIterator() const
{
    return SelectionParserSymbolIterator(
            new SelectionParserSymbolIterator::Impl(impl_->symbols_.end()));
}

void
SelectionParserSymbolTable::addVariable(const char                             *name,
                                        const gmx::SelectionTreeElementPointer &sel)
{
    // In the current parser implementation, a syntax error is produced before
    // this point is reached, but the check is here for robustness.
    Impl::SymbolMap::const_iterator other = impl_->symbols_.find(name);
    if (other != impl_->symbols_.end())
    {
        if (other->second->type() == SelectionParserSymbol::VariableSymbol)
        {
            GMX_THROW(InvalidInputError(
                              formatString("Reassigning variable '%s' is not supported",
                                           name)));
        }
        else
        {
            GMX_THROW(InvalidInputError(
                              formatString("Variable name '%s' conflicts with a reserved keyword",
                                           name)));
        }
    }
    Impl::SymbolPointer sym(new SelectionParserSymbol(
                                    new SelectionParserSymbol::Impl(
                                            SelectionParserSymbol::VariableSymbol, name)));
    sym->impl_->var_ = sel;
    impl_->addSymbol(move(sym));
}

void
SelectionParserSymbolTable::addMethod(const char          *name,
                                      gmx_ana_selmethod_t *method)
{
    if (impl_->symbols_.find(name) != impl_->symbols_.end())
    {
        GMX_THROW(APIError(
                          formatString("Method name '%s' conflicts with another symbol",
                                       name)));
    }
    Impl::SymbolPointer sym(new SelectionParserSymbol(
                                    new SelectionParserSymbol::Impl(
                                            SelectionParserSymbol::MethodSymbol, name)));
    sym->impl_->meth_ = method;
    impl_->addSymbol(move(sym));
}

} // namespace gmx
