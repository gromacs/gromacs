/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Declares and defines the XDRWrapper class.
 *
 * This wraps the xdr calls needed for checkpointing.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */

#ifndef GROMACS_XDRWRAPPER_H
#define GROMACS_XDRWRAPPER_H

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

template <typename T>
struct xdrType
{
};

template <>
struct xdrType<int>
{
    static const int value = xdr_datatype_int;
};

template <>
struct xdrType<int64_t >
{
    static const int value = xdr_datatype_int64;
};

template <>
struct xdrType<float>
{
    static const int value = xdr_datatype_float;
};

template <>
struct xdrType<double>
{
    static const int value = xdr_datatype_double;
};

class XDRWrapper
{
    public:
        XDRWrapper(
            t_fileio                *fp,
            std::vector<std::string> names = std::vector<std::string>()) :
            xd_(gmx_fio_getxdr(fp)), names_(std::move(names)) {}

        XDR* getXDR()
        {
            return xd_;
        }

        template<typename T>
        int doArrayRef(gmx::ArrayRef<T> array)
        {
            bool_t res = 0;

            int    numElements      = static_cast<int>(array.size());
            int    numElemInTheFile = numElements;
            /* Read/write the vector element count */
            res = xdr_int(xd_, &numElemInTheFile);
            if (res == 0)
            {
                return -1;
            }
            GMX_RELEASE_ASSERT(numElements == numElemInTheFile,
                               gmx::formatString("Vector length mismatch, code count is %d, file count is %d\n",
                                                 numElements, numElemInTheFile).c_str());

            /* Read/write the element data type */
            constexpr int xdrTypeInTheCode = xdrType<T>::value;
            int           xdrTypeInTheFile = xdrTypeInTheCode;
            res = xdr_int(xd_, &xdrTypeInTheFile);
            if (res == 0)
            {
                return -1;
            }

            bool typesMatch = (xdrTypeInTheFile == xdrTypeInTheCode);
            if (!typesMatch)
            {
                char buf[STRLEN];
                sprintf(buf, "Type mismatch, code precision is %s, file precision is %s",
                        xdr_datatype_names[xdrTypeInTheCode],
                        xdr_datatype_names[xdrTypeInTheFile]);

                /* Matching int and real should never occur, but check anyhow */
                if (xdrTypeInTheFile == xdr_datatype_int ||
                    xdrTypeInTheCode == xdr_datatype_int)
                {
                    gmx_fatal(FARGS, "Type %s: incompatible checkpoint formats or corrupted checkpoint file.", buf);
                }
            }

            T   * vp = array.data();

            char *vChar;
            if (typesMatch)
            {
                vChar = reinterpret_cast<char *>(vp);
            }
            else
            {
                snew(vChar, numElemInTheFile*sizeOfXdrType(xdrTypeInTheFile));
            }
            res = xdr_vector(xd_, vChar,
                             numElemInTheFile, sizeOfXdrType(xdrTypeInTheFile),
                             xdrProc(xdrTypeInTheFile));
            if (res == 0)
            {
                return -1;
            }

            if (!typesMatch)
            {
                /* In the old code float-double conversion came for free.
                 * In the new code we still support it, mainly because
                 * the tip4p_continue regression test makes use of this.
                 * It's an open question if we do or don't want to allow this.
                 */
                convertArrayRealPrecision(vChar, vp, numElemInTheFile);
                sfree(vChar);
            }

            return 0;
        }

        template<typename T>
        int doVector(std::vector<T> *vector)
        {
            return doArrayRef<T>(*vector);
        }

        void doCptStringErr(gmx::ArrayRef<char> s)
        {
            char *data = s.data();
            if (xdr_string(xd_, &data, s.size()) == 0)
            {
                cp_error();
            }
        }

        int doCptInt(int *i)
        {
            if (xdr_int(xd_, i) == 0)
            {
                return -1;
            }
            return 0;
        }

        int doCptUChars(int n, unsigned char *i)
        {
            bool_t res = 1;
            for (int j = 0; j < n && res; j++)
            {
                res &= xdr_u_char(xd_, &i[j]);
            }
            if (res == 0)
            {
                return -1;
            }

            return 0;
        }

        void doCptIntErr(int *i)
        {
            if (doCptInt(i) < 0)
            {
                cp_error();
            }
        }

        void doCptBoolErr(bool *b)
        {
            int i   = static_cast<int>(*b);

            if (doCptInt(&i) < 0)
            {
                cp_error();
            }

            *b = (i != 0);
        }

        void doCptInt64Err(int64_t *i)
        {
            if (xdr_int64(xd_, i) == 0)
            {
                cp_error();
            }
        }

        void doCptDoubleErr(double *f)
        {
            if (xdr_double(xd_, f) == 0)
            {
                cp_error();
            }
        }

        void doCptRealErr(real *f)
        {
#if GMX_DOUBLE
            bool_t res = xdr_double(xd_, f);
#else
            bool_t res = xdr_float(xd_, f);
#endif
            if (res == 0)
            {
                cp_error();
            }
        }

    private:
        void cp_warning(FILE *fp)
        {
            fprintf(fp, "\nWARNING: Checkpoint file is corrupted or truncated\n\n");
        }

        [[ noreturn ]] void cp_error()
        {
            gmx_fatal(FARGS, "Checkpoint file corrupted/truncated, or maybe you are out of disk space?");
        }

        //! \brief Convert a double array, typed char*, to float
        gmx_unused static void convertArrayRealPrecision(const char *c, float *v, int n)
        {
            auto *d = reinterpret_cast<const double *>(c);
            for (int i = 0; i < n; i++)
            {
                v[i] = static_cast<float>(d[i]);
            }
        }

        //! \brief Convert a float array, typed char*, to double
        static void convertArrayRealPrecision(const char *c, double *v, int n)
        {
            auto *f = reinterpret_cast<const float *>(c);
            for (int i = 0; i < n; i++)
            {
                v[i] = static_cast<double>(f[i]);
            }
        }

        //! \brief Generate an error for trying to convert to integer
        static void convertArrayRealPrecision(const char gmx_unused *c, int gmx_unused *v, int gmx_unused n)
        {
            GMX_RELEASE_ASSERT(false, "We only expect type mismatches between float and double, not integer");
        }

        //! \brief Generate an error for trying to convert to integer
        static void convertArrayRealPrecision(const char gmx_unused *c, int64_t gmx_unused *v, int gmx_unused n)
        {
            GMX_RELEASE_ASSERT(false, "We only expect type mismatches between float and double, not integer");
        }

        //! \brief Returns size in byte of an xdr_datatype
        inline unsigned int sizeOfXdrType(int xdrType)
        {
            switch (xdrType)
            {
                case xdr_datatype_int:
                    return sizeof(int);
                case xdr_datatype_float:
                    return sizeof(float);
                case xdr_datatype_double:
                    return sizeof(double);
                default: GMX_RELEASE_ASSERT(false, "XDR data type not implemented");
            }

            return 0;
        }

        //! \brief Returns the XDR process function for i/o of an XDR type
        inline xdrproc_t xdrProc(int xdrType)
        {
            switch (xdrType)
            {
                case xdr_datatype_int:
                    return reinterpret_cast<xdrproc_t>(xdr_int);
                case xdr_datatype_float:
                    return reinterpret_cast<xdrproc_t>(xdr_float);
                case xdr_datatype_double:
                    return reinterpret_cast<xdrproc_t>(xdr_double);
                default: GMX_RELEASE_ASSERT(false, "XDR data type not implemented");
            }

            return nullptr;
        }

        XDR *xd_;
        const std::vector<std::string> names_;
};

}      // namespace gmx

#endif //GROMACS_XDRWRAPPER_H
