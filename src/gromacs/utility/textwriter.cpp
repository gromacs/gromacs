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
 * Implements gmx::TextWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "textwriter.h"

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/nodelete.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

class TextWriter::Impl
{
    public:
        explicit Impl(const TextOutputStreamPointer &stream)
            : stream_(stream)
        {
        }

        TextOutputStreamPointer stream_;
};

TextWriter::TextWriter(const std::string &filename)
    : impl_(new Impl(TextOutputStreamPointer(new TextOutputFile(filename))))
{
}

TextWriter::TextWriter(TextOutputStream *stream)
    : impl_(new Impl(TextOutputStreamPointer(stream, no_delete<TextOutputStream>())))
{
}

TextWriter::TextWriter(const TextOutputStreamPointer &stream)
    : impl_(new Impl(stream))
{
}

TextWriter::~TextWriter()
{
}

TextOutputStream &TextWriter::stream()
{
    return *impl_->stream_;
}

void TextWriter::writeString(const char *str)
{
    impl_->stream_->write(str);
}

void TextWriter::writeString(const std::string &str)
{
    impl_->stream_->write(str.c_str());
}

void TextWriter::writeLine(const char *line)
{
    const size_t length = std::strlen(line);
    writeString(line);
    if (length == 0 || line[length-1] != '\n')
    {
        writeLine();
    }
}

void TextWriter::writeLine(const std::string &line)
{
    writeLine(line.c_str());
}

void TextWriter::writeLine()
{
    writeString("\n");
}

void TextWriter::close()
{
    impl_->stream_->close();
}

} // namespace gmx
