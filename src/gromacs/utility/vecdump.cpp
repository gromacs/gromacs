/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/utility/vecdump.h"

#include <cstdio>
#include <cstdlib>

#include <type_traits>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/txtdump.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

void dumpIntArrayRef(TextWriter* writer, const char* description, ArrayRef<const int> values, bool showIndices)
{
    if (!checkIfAvailable(writer, description, values.data()))
    {
        return;
    }

    printTitleAndSize(writer, description, values.size());
    {
        ScopedIndenter indenter(writer);
        for (int i = 0; i != values.ssize(); ++i)
        {
            writer->writeStringFormatted("%s[%d]=%d\n", description, showIndices ? i : -1, values[i]);
        }
    }
}

namespace
{

//! Common implementation method for several dumping methods
template<typename T>
std::enable_if_t<std::is_floating_point_v<T>, void> dumpFloatingPointArrayRef(TextWriter* writer,
                                                                              const char* description,
                                                                              ArrayRef<const T> values,
                                                                              bool showIndices)
{
    if (!checkIfAvailable(writer, description, values.data()))
    {
        return;
    }

    printTitleAndSize(writer, description, values.size());
    {
        ScopedIndenter indenter(writer);
        for (int i = 0; i != values.ssize(); ++i)
        {
            writer->writeStringFormatted("%s[%d]=%12.5e\n", description, showIndices ? i : -1, values[i]);
        }
    }
}

} // namespace

void dumpFloatArrayRef(TextWriter* writer, const char* description, ArrayRef<const float> values, bool showIndices)
{
    dumpFloatingPointArrayRef(writer, description, values, showIndices);
}

void dumpDoubleArrayRef(TextWriter* writer, const char* description, ArrayRef<const double> values, bool showIndices)
{
    dumpFloatingPointArrayRef(writer, description, values, showIndices);
}

void dumpRealArrayRef(TextWriter* writer, const char* description, ArrayRef<const real> values, bool showIndices)
{
    dumpFloatingPointArrayRef(writer, description, values, showIndices);
}

void dumpRvecArrayRef(TextWriter* writer, const char* description, ArrayRef<const RVec> values)
{
    if (!checkIfAvailable(writer, description, values.data()))
    {
        return;
    }

    const bool useLongFormat       = (std::getenv("GMX_PRINT_LONGFORMAT") != nullptr);
    const int  realFormatWidth     = useLongFormat ? 15 : 12;
    const int  realFormatPrecision = useLongFormat ? 8 : 5;

    writer->writeStringFormatted("%s (%dx%d):\n", description, values.ssize(), DIM);
    {
        ScopedIndenter indenter(writer);
        for (int i = 0; i < values.ssize(); i++)
        {
            writer->writeStringFormatted("%s[%5d]={%*.*e, %*.*e, %*.*e}\n",
                                         description,
                                         i,
                                         realFormatWidth,
                                         realFormatPrecision,
                                         values[i][XX],
                                         realFormatWidth,
                                         realFormatPrecision,
                                         values[i][YY],
                                         realFormatWidth,
                                         realFormatPrecision,
                                         values[i][ZZ]);
        }
    }
}

} // namespace gmx

void pr_ivec(FILE* fp, int indent, const char* title, const int vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%d\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_ivec_block(FILE* fp, int indent, const char* title, const int vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        i      = 0;
        while (i < n)
        {
            j = i + 1;
            while (j < n && vec[j] == vec[j - 1] + 1)
            {
                j++;
            }
            /* Print consecutive groups of 3 or more as blocks */
            if (j - i < 3)
            {
                while (i < j)
                {
                    pr_indent(fp, indent);
                    fprintf(fp, "%s[%d]=%d\n", title, bShowNumbers ? i : -1, vec[i]);
                    i++;
                }
            }
            else
            {
                pr_indent(fp, indent);
                fprintf(fp,
                        "%s[%d,...,%d] = {%d,...,%d}\n",
                        title,
                        bShowNumbers ? i : -1,
                        bShowNumbers ? j - 1 : -1,
                        vec[i],
                        vec[j - 1]);
                i = j;
            }
        }
    }
}

void pr_ivecs(FILE* fp, int indent, const char* title, const ivec vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={", title, bShowNumbers ? i : -1);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, "%d", vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

template<typename T>
static void printRealVector(FILE* fp, int indent, const char* title, const T vec[], int n, gmx_bool bShowNumbers)
{
    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (int i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%12.5e\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_rvec(FILE* fp, int indent, const char* title, const real vec[], int n, gmx_bool bShowNumbers)
{
    printRealVector<real>(fp, indent, title, vec, n, bShowNumbers);
}

void pr_fvec(FILE* fp, int indent, const char* title, const float vec[], int n, gmx_bool bShowNumbers)
{
    printRealVector<float>(fp, indent, title, vec, n, bShowNumbers);
}

void pr_dvec(FILE* fp, int indent, const char* title, const double vec[], int n, gmx_bool bShowNumbers)
{
    printRealVector<double>(fp, indent, title, vec, n, bShowNumbers);
}


void pr_rvecs(FILE* fp, int indent, const char* title, const rvec vec[], int n)
{
    const char* fshort = "%12.5e";
    const char* flong  = "%15.8e";
    const char* format;
    int         i, j;

    if (std::getenv("GMX_PRINT_LONGFORMAT") != nullptr)
    {
        format = flong;
    }
    else
    {
        format = fshort;
    }

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

void prRVecs(FILE* fp, int indent, const char* title, gmx::ArrayRef<const gmx::RVec> vec)
{
    pr_rvecs(fp, indent, title, as_rvec_array(vec.data()), vec.size());
}
