/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares utility classes for testing multi-line strings against reference data.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TEXTBLOCKMATCHERS_H
#define GMX_TESTUTILS_TEXTBLOCKMATCHERS_H

#include <memory>
#include <string>

namespace gmx
{

class TextInputStream;

namespace test
{

class TestReferenceChecker;

/*! \libinternal \brief
 * Represents a text matcher, matching text stream contents against reference
 * data.
 *
 * Typical pattern of declaring such matchers is to
 *  - Create a factory that implements ITextBlockMatcherSettings,
 *  - Make that factory provide any necessary parameters that the matcher needs,
 *    using a "named parameter" idiom (see XvgMatch for an example), and
 *  - Make the factory create and return an instance of an internal
 *    implementation class that implements ITextBlockMatcher and provides
 *    the actual matching logic.
 *
 * Any method that then wants to accept a matcher can accept a
 * ITextBlockMatcherSettings.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class ITextBlockMatcher
{
    public:
        virtual ~ITextBlockMatcher();

        /*! \brief
         * Matches contents of a stream against reference data.
         *
         * \param  stream   Stream to match.
         * \param  checker  Checker to use for matching.
         *
         * The method can change the state of the provided checker (e.g., by
         * changing the default tolerance).
         * The caller is responsible of providing a checker where such state
         * changes do not matter.
         */
        virtual void checkStream(TextInputStream      *stream,
                                 TestReferenceChecker *checker) = 0;
};

//! Smart pointer for managing a ITextBlockMatcher.
typedef std::unique_ptr<ITextBlockMatcher> TextBlockMatcherPointer;

/*! \libinternal \brief
 * Represents a factory for creating a text matcher.
 *
 * See derived classes for available matchers.  Each derived class represents
 * one type of matcher (see ITextBlockMatcher), and provides any methods
 * necessary to pass parameters to such a matcher.  Methods that accept a
 * matcher can then take in this interface, and call createMatcher() to use the
 * matcher that the caller of the method specifies.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class ITextBlockMatcherSettings
{
    public:
        //! Factory method that constructs the matcher after parameters are set.
        virtual TextBlockMatcherPointer createMatcher() const = 0;

    protected:
        virtual ~ITextBlockMatcherSettings();
};

/*! \libinternal \brief
 * Use an exact text match (the contents should be exactly equal).
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class ExactTextMatch : public ITextBlockMatcherSettings
{
    public:
        virtual TextBlockMatcherPointer createMatcher() const;
};

/*! \libinternal \brief
 * Do not match the text (the contents are ignored).
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class NoTextMatch : public ITextBlockMatcherSettings
{
    public:
        virtual TextBlockMatcherPointer createMatcher() const;
};

} // namespace test
} // namespace gmx

#endif
