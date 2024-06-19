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
 * Implements gmx::TextWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/textwriter.h"

#include <cstdarg>
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/nodelete.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

class TextWriter::Impl
{
public:
    explicit Impl(const TextOutputStreamPointer& stream) :
        stream_(stream), newLineCount_(2), currentLineLength_(0), pendingNewLine_(false)
    {
        wrapper_.settings().setKeepFinalSpaces(true);
    }

    void writeRawString(const char* str)
    {
        if (pendingNewLine_ && str[0] != '\n')
        {
            stream_->write("\n");
        }
        pendingNewLine_         = false;
        const char* lastNewLine = std::strrchr(str, '\n');
        if (lastNewLine == nullptr)
        {
            newLineCount_ = 0;
            currentLineLength_ += std::strlen(str);
        }
        else if (lastNewLine[1] != '\0')
        {
            newLineCount_ = 0;
            currentLineLength_ += std::strlen(lastNewLine + 1);
        }
        else
        {
            currentLineLength_ = 0;
            int newLineCount   = 0;
            while (lastNewLine >= str && *lastNewLine == '\n')
            {
                ++newLineCount;
                --lastNewLine;
            }
            if (lastNewLine >= str)
            {
                newLineCount_ = 0;
            }
            newLineCount_ += newLineCount;
        }
        stream_->write(str);
    }
    void writeRawString(const std::string& str) { writeRawString(str.c_str()); }

    void writeWrappedString(const std::string& str)
    {
        if (newLineCount_ > 0)
        {
            writeRawString(wrapper_.wrapToString(str));
        }
        else
        {
            writeRawString(str);
        }
    }

    TextOutputStreamPointer stream_;
    TextLineWrapper         wrapper_;
    int                     newLineCount_;
    int                     currentLineLength_;
    bool                    pendingNewLine_;
};

// static
void TextWriter::writeFileFromString(const std::filesystem::path& filename, const std::string& text)
{
    TextWriter file(filename);
    file.writeString(text);
    file.close();
}

TextWriter::TextWriter(const std::filesystem::path& filename) :
    impl_(new Impl(TextOutputStreamPointer(new TextOutputFile(filename))))
{
}

TextWriter::TextWriter(FILE* fp) : impl_(new Impl(TextOutputStreamPointer(new TextOutputFile(fp))))
{
}

TextWriter::TextWriter(TextOutputStream* stream) :
    impl_(new Impl(TextOutputStreamPointer(stream, no_delete<TextOutputStream>())))
{
}

TextWriter::TextWriter(const TextOutputStreamPointer& stream) : impl_(new Impl(stream)) {}

TextWriter::~TextWriter() {}

TextLineWrapperSettings& TextWriter::wrapperSettings()
{
    return impl_->wrapper_.settings();
}

void TextWriter::writeString(const char* str)
{
    if (impl_->wrapper_.isTrivial())
    {
        impl_->writeRawString(str);
    }
    else
    {
        impl_->writeWrappedString(str);
    }
}

void TextWriter::writeString(const std::string& str)
{
    impl_->writeWrappedString(str);
}

void TextWriter::writeStringFormatted(const char* fmt, ...)
{
    std::va_list ap;

    va_start(ap, fmt);
    writeString(formatStringV(fmt, ap));
    va_end(ap);
}

void TextWriter::writeLine(const char* line)
{
    writeString(line);
    ensureLineBreak();
}

void TextWriter::writeLine(const std::string& line)
{
    writeString(line);
    ensureLineBreak();
}

void TextWriter::writeLineFormatted(const char* fmt, ...)
{
    std::va_list ap;

    va_start(ap, fmt);
    writeString(formatStringV(fmt, ap));
    va_end(ap);
    ensureLineBreak();
}

void TextWriter::writeLine()
{
    impl_->writeRawString("\n");
}

void TextWriter::ensureLineBreak()
{
    if (impl_->newLineCount_ == 0)
    {
        impl_->writeRawString("\n");
    }
}

void TextWriter::ensureEmptyLine()
{
    ensureLineBreak();
    if (impl_->newLineCount_ < 2)
    {
        impl_->pendingNewLine_ = true;
    }
}

void TextWriter::close()
{
    impl_->stream_->close();
}

} // namespace gmx
