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
/*! \internal \file
 * \brief
 * Implements simple math functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_math
 */

#include "gmxpre.h"

#include "gromacs/math/functions.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>

#include <array>
#include <limits>

#if GMX_NATIVE_WINDOWS
#    include <intrin.h> // _BitScanReverse, _BitScanReverse64
#endif

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

unsigned int log2I(std::uint32_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(0) is undefined");
#if HAVE_BUILTIN_CLZ
    // gcc, clang. xor with sign bit should be optimized out
    return __builtin_clz(n) ^ 31U;
#elif HAVE_BITSCANREVERSE
    // icc, MSVC
    {
        unsigned long res;
        _BitScanReverse(&res, static_cast<unsigned long>(n));
        return static_cast<unsigned int>(res);
    }
#elif HAVE_CNTLZ4
    return 31 - __cntlz4(n);
#else
    // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup

    static const std::array<char, 256> log2TableByte = {
        { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
          5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 }
    };

    unsigned int result;
    unsigned int tmp1, tmp2;

    if ((tmp1 = n >> 16) != 0)
    {
        result = ((tmp2 = tmp1 >> 8) != 0) ? 24 + log2TableByte[tmp2] : 16 + log2TableByte[tmp1];
    }
    else
    {
        result = ((tmp2 = n >> 8) != 0) ? 8 + log2TableByte[tmp2] : log2TableByte[n];
    }
    return result;
#endif
}


unsigned int log2I(std::uint64_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(0) is undefined");
#if HAVE_BUILTIN_CLZLL
    // gcc, icc, clang. xor with sign bit should be optimized out
    return __builtin_clzll(n) ^ 63U;
#elif HAVE_BITSCANREVERSE64
    unsigned long res;
    _BitScanReverse64(&res, static_cast<unsigned __int64>(n));
    return static_cast<unsigned int>(res);
#elif HAVE_CNTLZ8
    return 63 - __cntlz8(n);
#else

    // No 64-bit log2 instrinsic available. Solve it by calling our internal
    // 32-bit version (which in turn might defer to a software solution)

    std::uint32_t high32Bits = static_cast<std::uint32_t>(n >> 32);
    std::uint32_t result;

    if (high32Bits)
    {
        result = log2I(high32Bits) + 32;
    }
    else
    {
        result = log2I(static_cast<std::uint32_t>(n));
    }

    return result;
#endif
}

unsigned int log2I(std::int32_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(n) for n<=0 is undefined");
    return log2I(static_cast<std::uint32_t>(n));
}

unsigned int log2I(std::int64_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(n) for n<=0 is undefined");
    return log2I(static_cast<std::uint64_t>(n));
}

std::int64_t greatestCommonDivisor(std::int64_t p, std::int64_t q)
{
    while (q != 0)
    {
        std::int64_t tmp = q;
        q                = p % q;
        p                = tmp;
    }
    return p;
}

/* These inverse error function implementations are simplified versions
 * of the algorithms and polynomia developed for Boost by
 * John Maddock in 2006, and modified by Jeremy William Murphy in 2015.
 *
 * You might prefer the original version to avoid any bugs we have
 * introduced, but in the spirit of not unnecessarily restricting a
 * more liberal license, you are most welcome to also use and redistribute
 * the Gromacs versions of these functions under the
 * Boost Software License, Version 1.0. (http://www.boost.org/LICENSE_1_0.txt):
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

double erfinv(double arg)
{
    double x = std::abs(arg);

    if (x > 1.0)
    {
        return std::nan("");
    }
    else if (arg == 1.0)
    {
        return std::numeric_limits<double>::infinity();
    }
    else if (arg == -1.0)
    {
        return -std::numeric_limits<double>::infinity();
    }
    else if (x == 0.0)
    {
        return 0.0;
    }

    double y = 1.0 - x;
    double p, q, res;

    if (x <= 0.5)
    {
        // Rational approximation for |x| in range [0,0.5]
        p = -0.00538772965071242932965;
        p = p * x + 0.00822687874676915743155;
        p = p * x + 0.0219878681111168899165;
        p = p * x - 0.0365637971411762664006;
        p = p * x - 0.0126926147662974029034;
        p = p * x + 0.0334806625409744615033;
        p = p * x - 0.00836874819741736770379;
        p = p * x - 0.000508781949658280665617;

        q = 0.000886216390456424707504;
        q = q * x - 0.00233393759374190016776;
        q = q * x + 0.0795283687341571680018;
        q = q * x - 0.0527396382340099713954;
        q = q * x - 0.71228902341542847553;
        q = q * x + 0.662328840472002992063;
        q = q * x + 1.56221558398423026363;
        q = q * x - 1.56574558234175846809;
        q = q * x - 0.970005043303290640362;
        q = q * x + 1.0;

        double t = x * (x + 10.0);
        res      = t * 0.0891314744949340820313 + t * p / q;
    }
    else if (x <= 0.75)
    {
        // Rational approx for |x| in range ]0.5,0.75]
        double z = y - 0.25;
        p        = -3.67192254707729348546;
        p        = p * z + 21.1294655448340526258;
        p        = p * z + 17.445385985570866523;
        p        = p * z - 44.6382324441786960818;
        p        = p * z - 18.8510648058714251895;
        p        = p * z + 17.6447298408374015486;
        p        = p * z + 8.37050328343119927838;
        p        = p * z + 0.105264680699391713268;
        p        = p * z - 0.202433508355938759655;

        q = 1.72114765761200282724;
        q = q * z - 22.6436933413139721736;
        q = q * z + 10.8268667355460159008;
        q = q * z + 48.5609213108739935468;
        q = q * z - 20.1432634680485188801;
        q = q * z - 28.6608180499800029974;
        q = q * z + 3.9713437953343869095;
        q = q * z + 6.24264124854247537712;
        q = q * z + 1.0;

        double t = std::sqrt(-2.0 * std::log(y));
        res      = t / (2.249481201171875 + p / q);
    }
    else
    {
        // Branch for 0.75 < x < 1 (meaning 0 < y <= 0.25)
        double t = std::sqrt(-std::log(y));

        if (t < 3.0)
        {
            // |x| in range ]0.75,0.99987659019591335]
            double z = t - 1.125;
            p        = -0.681149956853776992068e-9;
            p        = p * z + 0.285225331782217055858e-7;
            p        = p * z - 0.679465575181126350155e-6;
            p        = p * z + 0.00214558995388805277169;
            p        = p * z + 0.0290157910005329060432;
            p        = p * z + 0.142869534408157156766;
            p        = p * z + 0.337785538912035898924;
            p        = p * z + 0.387079738972604337464;
            p        = p * z + 0.117030156341995252019;
            p        = p * z - 0.163794047193317060787;
            p        = p * z - 0.131102781679951906451;

            q = 0.01105924229346489121;
            q = q * z + 0.152264338295331783612;
            q = q * z + 0.848854343457902036425;
            q = q * z + 2.59301921623620271374;
            q = q * z + 4.77846592945843778382;
            q = q * z + 5.38168345707006855425;
            q = q * z + 3.46625407242567245975;
            q = q * z + 1.0;

            res = t * 0.807220458984375 + t * p / q;
        }
        else
        {
            // |x| in range ]0.99987659019591335,1[
            double z = t - 3.0;
            p        = 0.266339227425782031962e-11;
            p        = p * z - 0.230404776911882601748e-9;
            p        = p * z + 0.460469890584317994083e-5;
            p        = p * z + 0.000157544617424960554631;
            p        = p * z + 0.00187123492819559223345;
            p        = p * z + 0.00950804701325919603619;
            p        = p * z + 0.0185573306514231072324;
            p        = p * z - 0.00222426529213447927281;
            p        = p * z - 0.0350353787183177984712;

            q = 0.764675292302794483503e-4;
            q = q * z + 0.00263861676657015992959;
            q = q * z + 0.0341589143670947727934;
            q = q * z + 0.220091105764131249824;
            q = q * z + 0.762059164553623404043;
            q = q * z + 1.3653349817554063097;
            q = q * z + 1.0;

            res = t * 0.93995571136474609375 + t * p / q;
        }
    }
    return std::copysign(res, arg);
}

float erfinv(float arg)
{
    float x = std::abs(arg);

    if (x > 1.0F)
    {
        return std::nan("");
    }
    else if (arg == 1.0F)
    {
        return std::numeric_limits<float>::infinity();
    }
    else if (arg == -1.0F)
    {
        return -std::numeric_limits<float>::infinity();
    }
    else if (x == 0.0F)
    {
        return 0.0F;
    }

    float y = 1.0F - x;
    float p, q, res;

    // It is likely possible to use polynomia of slightly
    // lower order in single precision, but to really
    // optimize it would also require changing the intervals,
    // and adopting factors to exact fp32 representation in
    // IEEE-754. Given that we don't use erfinv() in any
    // tight loops it's not needed for now, so we leave it
    // as an exercise to the developer reading this note.
    if (x <= 0.5F)
    {
        // Rational approximation for |x| in range [0,0.5]
        p = -0.00538772965071242932965F;
        p = p * x + 0.00822687874676915743155F;
        p = p * x + 0.0219878681111168899165F;
        p = p * x - 0.0365637971411762664006F;
        p = p * x - 0.0126926147662974029034F;
        p = p * x + 0.0334806625409744615033F;
        p = p * x - 0.00836874819741736770379F;
        p = p * x - 0.000508781949658280665617F;

        q = 0.000886216390456424707504F;
        q = q * x - 0.00233393759374190016776F;
        q = q * x + 0.0795283687341571680018F;
        q = q * x - 0.0527396382340099713954F;
        q = q * x - 0.71228902341542847553F;
        q = q * x + 0.662328840472002992063F;
        q = q * x + 1.56221558398423026363F;
        q = q * x - 1.56574558234175846809F;
        q = q * x - 0.970005043303290640362F;
        q = q * x + 1.0F;

        float t = x * (x + 10.0F);
        res     = t * 0.0891314744949340820313F + t * p / q;
    }
    else if (x <= 0.75F)
    {
        // Rational approx for |x| in range ]0.5,0.75]
        float z = y - 0.25F;
        p       = -3.67192254707729348546F;
        p       = p * z + 21.1294655448340526258F;
        p       = p * z + 17.445385985570866523F;
        p       = p * z - 44.6382324441786960818F;
        p       = p * z - 18.8510648058714251895F;
        p       = p * z + 17.6447298408374015486F;
        p       = p * z + 8.37050328343119927838F;
        p       = p * z + 0.105264680699391713268F;
        p       = p * z - 0.202433508355938759655F;

        q = 1.72114765761200282724F;
        q = q * z - 22.6436933413139721736F;
        q = q * z + 10.8268667355460159008F;
        q = q * z + 48.5609213108739935468F;
        q = q * z - 20.1432634680485188801F;
        q = q * z - 28.6608180499800029974F;
        q = q * z + 3.9713437953343869095F;
        q = q * z + 6.24264124854247537712F;
        q = q * z + 1.0F;

        float t = std::sqrt(-2.0F * std::log(y));
        res     = t / (2.249481201171875F + p / q);
    }
    else
    {
        // Branch for 0.75 < x < 1 (meaning 0 < y <= 0.25)
        float t = std::sqrt(-std::log(y));

        if (t < 3.0F)
        {
            // |x| in range ]0.75,0.99987659019591335]
            float z = t - 1.125F;
            p       = -0.681149956853776992068e-9F;
            p       = p * z + 0.285225331782217055858e-7F;
            p       = p * z - 0.679465575181126350155e-6F;
            p       = p * z + 0.00214558995388805277169F;
            p       = p * z + 0.0290157910005329060432F;
            p       = p * z + 0.142869534408157156766F;
            p       = p * z + 0.337785538912035898924F;
            p       = p * z + 0.387079738972604337464F;
            p       = p * z + 0.117030156341995252019F;
            p       = p * z - 0.163794047193317060787F;
            p       = p * z - 0.131102781679951906451F;

            q = 0.01105924229346489121F;
            q = q * z + 0.152264338295331783612F;
            q = q * z + 0.848854343457902036425F;
            q = q * z + 2.59301921623620271374F;
            q = q * z + 4.77846592945843778382F;
            q = q * z + 5.38168345707006855425F;
            q = q * z + 3.46625407242567245975F;
            q = q * z + 1.0F;

            res = t * 0.807220458984375F + t * p / q;
        }
        else
        {
            // |x| in range ]0.99987659019591335,1[
            float z = t - 3.0F;
            p       = 0.266339227425782031962e-11F;
            p       = p * z - 0.230404776911882601748e-9F;
            p       = p * z + 0.460469890584317994083e-5F;
            p       = p * z + 0.000157544617424960554631F;
            p       = p * z + 0.00187123492819559223345F;
            p       = p * z + 0.00950804701325919603619F;
            p       = p * z + 0.0185573306514231072324F;
            p       = p * z - 0.00222426529213447927281F;
            p       = p * z - 0.0350353787183177984712F;

            q = 0.764675292302794483503e-4F;
            q = q * z + 0.00263861676657015992959F;
            q = q * z + 0.0341589143670947727934F;
            q = q * z + 0.220091105764131249824F;
            q = q * z + 0.762059164553623404043F;
            q = q * z + 1.3653349817554063097F;
            q = q * z + 1.0F;

            res = t * 0.93995571136474609375F + t * p / q;
        }
    }
    return std::copysign(res, arg);
}


} // namespace gmx
