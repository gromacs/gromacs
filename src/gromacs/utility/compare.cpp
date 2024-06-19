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
/* This file is completely threadsafe - keep it that way! */

#include "gmxpre.h"

#include "gromacs/utility/compare.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"

void cmp_int(FILE* fp, const char* s, int index, int i1, int i2)
{
    if (i1 != i2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%d - %d)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%d - %d)\n", s, i1, i2);
        }
    }
}

void cmp_int64(FILE* fp, const char* s, int64_t i1, int64_t i2)
{
    if (i1 != i2)
    {
        fprintf(fp, "%s (", s);
        fprintf(fp, "%" PRId64, i1);
        fprintf(fp, " - ");
        fprintf(fp, "%" PRId64, i2);
        fprintf(fp, ")\n");
    }
}

void cmp_us(FILE* fp, const char* s, int index, unsigned short i1, unsigned short i2)
{
    if (i1 != i2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%hu - %hu)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%hu - %hu)\n", s, i1, i2);
        }
    }
}

void cmp_uc(FILE* fp, const char* s, int index, unsigned char i1, unsigned char i2)
{
    if (i1 != i2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%d - %d)\n", s, index, int{ i1 }, int{ i2 });
        }
        else
        {
            fprintf(fp, "%s (%d - %d)\n", s, int{ i1 }, int{ i2 });
        }
    }
}

gmx_bool cmp_bool(FILE* fp, const char* s, int index, gmx_bool b1, gmx_bool b2)
{
    if (b1 != b2)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%s - %s)\n", s, index, gmx::boolToString(b1), gmx::boolToString(b2));
        }
        else
        {
            fprintf(fp, "%s (%s - %s)\n", s, gmx::boolToString(b1), gmx::boolToString(b2));
        }
    }
    return b1 && b2;
}

void cmp_str(FILE* fp, const char* s, int index, const char* s1, const char* s2)
{
    if (std::strcmp(s1, s2) != 0)
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%d] (%s - %s)\n", s, index, s1, s2);
        }
        else
        {
            fprintf(fp, "%s (%s - %s)\n", s, s1, s2);
        }
    }
}

gmx_bool equal_real(real i1, real i2, real ftol, real abstol)
{
    return ((2 * std::fabs(i1 - i2) <= (std::fabs(i1) + std::fabs(i2)) * ftol)
            || std::fabs(i1 - i2) <= abstol);
}

gmx_bool equal_float(float i1, float i2, float ftol, float abstol)
{
    return ((2 * std::fabs(i1 - i2) <= (std::fabs(i1) + std::fabs(i2)) * ftol)
            || std::fabs(i1 - i2) <= abstol);
}

gmx_bool equal_double(double i1, double i2, real ftol, real abstol)
{
    return ((2 * std::fabs(i1 - i2) <= (std::fabs(i1) + std::fabs(i2)) * ftol)
            || std::fabs(i1 - i2) <= abstol);
}

void cmp_real(FILE* fp, const char* s, int index, real i1, real i2, real ftol, real abstol)
{
    if (!equal_real(i1, i2, ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%2d] (%e - %e)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%e - %e)\n", s, i1, i2);
        }
    }
}

void cmp_float(FILE* fp, const char* s, int index, float i1, float i2, float ftol, float abstol)
{
    if (!equal_float(i1, i2, ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%2d] (%e - %e)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%e - %e)\n", s, i1, i2);
        }
    }
}

void cmp_double(FILE* fp, const char* s, int index, double i1, double i2, double ftol, double abstol)
{
    if (!equal_double(i1, i2, ftol, abstol))
    {
        if (index != -1)
        {
            fprintf(fp, "%s[%2d] (%16.9e - %16.9e)\n", s, index, i1, i2);
        }
        else
        {
            fprintf(fp, "%s (%16.9e - %16.9e)\n", s, i1, i2);
        }
    }
}
