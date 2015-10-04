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
/*! \internal \file
 * \brief
 * Implements gmx::TextReader.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "textreader.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/nodelete.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

// static
std::string TextReader::readFileToString(const char *filename)
{
    TextReader  reader(filename);
    std::string result(reader.readAll());
    reader.close();
    return result;
}

// static
std::string TextReader::readFileToString(const std::string &filename)
{
    return readFileToString(filename.c_str());
}

class TextReader::Impl
{
    public:
        explicit Impl(const TextInputStreamPointer &stream)
            : stream_(stream), lineIndex_(0), currentLine_() {}

        bool readLine()
        {
            /*! As each line is read we prepare for an error message
             * that might be given later by makeError(). */
            lineIndex_++;
            return stream_->readLine(&currentLine_);
        }

        InvalidInputError makeError(const std::string &message) const
        {
            if (!currentLine_.empty())
            {
                return InvalidInputError(message + formatString("\n on line %d, which was\n '%s'",
                                                                lineIndex_, stripString(currentLine_).c_str()));
            }
            else
            {
                return InvalidInputError(message);
            }
        }

        TextInputStreamPointer stream_;
        int                    lineIndex_;
        std::string            currentLine_;
};

TextReader::TextReader(const std::string &filename)
    : impl_(new Impl(TextInputStreamPointer(new TextInputFile(filename))))
{
}

TextReader::TextReader(TextInputStream *stream)
    : impl_(new Impl(TextInputStreamPointer(stream, no_delete<TextInputStream>())))
{
}

TextReader::TextReader(const TextInputStreamPointer &stream)
    : impl_(new Impl(stream))
{
}

TextReader::~TextReader()
{
}

bool TextReader::readLine()
{
    return impl_->readLine();
}

const std::string &TextReader::currentLine()
{
    const size_t endPos = impl_->currentLine_.find_last_not_of(" \t\r\n");
    if (endPos != std::string::npos)
    {
        impl_->currentLine_.resize(endPos + 1);
    }
    return impl_->currentLine_;
}

std::string TextReader::readAll()
{
    std::string result;
    while (readLine())
    {
        result.append(currentLine());
    }
    return result;
}

InvalidInputError TextReader::makeError(const std::string &message) const
{
    return impl_->makeError(message);
}

void TextReader::close()
{
    impl_->stream_->close();
}

} // namespace gmx
