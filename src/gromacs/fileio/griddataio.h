/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Reading and writing routines for volume data formats ccp4, mrc and imod.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_GRIDDATAIO_H_
#define GMX_FILEIO_GRIDDATAIO_H_
#include <array>
#include <memory>
#include <string>

#include "gromacs/math/griddata/griddata.h"

namespace gmx
{
struct MrcMetaData;
/*! \brief
 * Read and write real-space real-valued volume data files
 * according to the electron microscopy data bank (EMDB) standard.
 *
 * The formatting guraranties compliance with 3D EM maps described in
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 * However, other ccp4, mrc, imod and map formats might be compatible.
 *
 * Future implementations for reading crystallographic data or image stacks might
 * want to split MrcFile into an abstract griddata base class and respective
 * child implementatons, if demand exists.
 */
class MrcFile
{
    public:
        MrcFile();
        ~MrcFile();

        /*! \brief Write real-spaced, real-valued griddata to file
         * with default metadata for 3D electron microscopy data.
         *
         * \param[in] filename name of the file to write the griddata to, typically *.cpp4, *.mrc or *.map
         * \param[in] grid_data real-valued, real-space data on a grid
         */
        void write(std::string filename, const GridDataReal3D &grid_data);

        /*! \brief Write real-spaced, real-valued griddata to file with user-defined metadata.
         *
         * \param[in] filename name of the file to write the griddata to, typically *.cpp4, *.mrc or *.map
         * \param[in] grid_data real-valued, real-space data on a grid
         * \param[in] meta struct with own metadata
         * \param[in] bOwnGridStats calculate min, max, mean and rms self if true, otherwise copy from metadata
         */
        void write_with_own_meta(std::string filename, GridDataReal3D &grid_data, const MrcMetaData &meta, bool bOwnGridStats);

        /*! \brief Reads real-spaced, real-valued griddata from file.
         *
         * \param[in] filename name of the file from which to read the griddata, typically *.cpp4, *.mrc or *.map
         * \returns grid_data will be filled with real-valued, real-space data on a grid upon succesful reading
         */
        GridDataReal3D read(std::string filename);
        /*! \brief Reads real-spaced, real-valued voxel data and the map header / metadata.
         *
         * \param[in] filename name of the file from which to read the griddata, typically *.cpp4, *.mrc or *.map
         * \param[in] meta returns the metadata from reading; previous content will be overwritten
         * \returns grid_data will be filled with real-valued, real-space data on a grid upon succesful reading
         */
        GridDataReal3D read_with_meta(std::string filename, MrcMetaData *meta);

        /*! \brief Reads
         *
         * \param[in] filename name of the file from which to read the griddata, typically *.cpp4, *.mrc or *.map
         * \param[in] meta returns the metadata from reading; previous content will be overwritten
         */
        void read_meta(std::string filename, MrcMetaData *meta);


    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

class Df3File
{
    private:
        class SuccessfulDf3Write
        {
            public:
                void writePovray();
                SuccessfulDf3Write(std::string filename, const GridDataReal3D &grid_data);
            private:
                std::string           filename_;
                const GridDataReal3D &gridData_;
        };
    public:
        SuccessfulDf3Write write(std::string filename, const GridDataReal3D &grid_data);
};

}      // namespace gmx
#endif /* end of include guard: GMX_FILEIO_GRIDDATAIO_H_ */
