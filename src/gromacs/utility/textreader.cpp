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
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "textreader.h"

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/nodelete.h"
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
            : stream_(stream)
        {
        }

        TextInputStreamPointer stream_;
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

bool TextReader::readLine(std::string *line)
{
    return impl_->stream_->readLine(line);
}

bool TextReader::readLineTrimmed(std::string *line)
{
    if (!readLine(line))
    {
        return false;
    }
    const size_t endPos = line->find_last_not_of(" \t\r\n");
    if (endPos != std::string::npos)
    {
        line->resize(endPos + 1);
    }
    return true;
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
