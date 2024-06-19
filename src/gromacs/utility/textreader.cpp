/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Implements gmx::TextReader.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/textreader.h"

#include <cstddef>

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/nodelete.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

// static
std::string TextReader::readFileToString(const std::string& filename)
{
    return readFileToString(std::filesystem::path{ filename });
}

// static
std::string TextReader::readFileToString(const std::filesystem::path& filename)
{
    TextReader  reader(filename);
    std::string result(reader.readAll());
    reader.close();
    return result;
}

//! Implementation class
class TextReader::Impl
{
public:
    //! Constructor.
    explicit Impl(const TextInputStreamPointer& stream) :
        stream_(stream),
        trimLeadingWhiteSpace_(false),
        trimTrailingWhiteSpace_(false),
        trimTrailingComment_(false),
        commentChar_(0)
    {
    }

    //! Stream used by this reader.
    TextInputStreamPointer stream_;
    //! Whether leading whitespace should be removed.
    bool trimLeadingWhiteSpace_;
    //! Whether trailing whitespace should be removed.
    bool trimTrailingWhiteSpace_;
    //! Whether a trailing comment should be removed.
    bool trimTrailingComment_;
    /*! \brief Character that denotes the start of a comment on a line.
     *
     * Zero until TextReader::setTrimTrailingComment is called to
     * activate such trimming with a given character. */
    char commentChar_;
};

TextReader::TextReader(const std::filesystem::path& filename) :
    impl_(new Impl(TextInputStreamPointer(new TextInputFile(filename))))
{
}

TextReader::TextReader(TextInputStream* stream) :
    impl_(new Impl(TextInputStreamPointer(stream, no_delete<TextInputStream>())))
{
}

TextReader::TextReader(const TextInputStreamPointer& stream) : impl_(new Impl(stream)) {}

TextReader::~TextReader() {}

bool TextReader::readLine(std::string* linePtr)
{
    if (!impl_->stream_->readLine(linePtr))
    {
        return false;
    }
    auto&      line              = *linePtr;
    const char whiteSpaceChars[] = " \t\r\n";
    if (impl_->trimLeadingWhiteSpace_)
    {
        const size_t endPos = line.find_first_not_of(whiteSpaceChars);
        if (endPos == std::string::npos)
        {
            line.resize(0);
        }
        else
        {
            line = line.substr(endPos, std::string::npos);
        }
    }
    if (impl_->trimTrailingComment_)
    {
        auto commentPos = line.find(impl_->commentChar_);
        if (commentPos != std::string::npos)
        {
            line.resize(commentPos);
        }
    }
    if (impl_->trimTrailingWhiteSpace_)
    {
        const size_t endPos = line.find_last_not_of(whiteSpaceChars);
        if (endPos == std::string::npos)
        {
            line.resize(0);
        }
        else
        {
            line.resize(endPos + 1);
        }
    }
    return true;
}

void TextReader::setTrimLeadingWhiteSpace(bool doTrimming)
{
    impl_->trimLeadingWhiteSpace_ = doTrimming;
}

void TextReader::setTrimTrailingWhiteSpace(bool doTrimming)
{
    impl_->trimTrailingWhiteSpace_ = doTrimming;
}

void TextReader::setTrimTrailingComment(bool doTrimming, char commentChar)
{
    impl_->trimTrailingComment_ = doTrimming;
    if (impl_->trimTrailingComment_)
    {
        impl_->commentChar_ = commentChar;
    }
}

std::string TextReader::readAll()
{
    std::string result;
    std::string line;
    while (readLine(&line))
    {
        result.append(line);
    }
    return result;
}

void TextReader::close()
{
    impl_->stream_->close();
}

} // namespace gmx
