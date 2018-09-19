/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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

#include "xplor.h"

#include <string.h>

#include <cmath>

#include "gromacs/utility/exceptions.h"

namespace gmx
{


XplorData::XplorData(FILE *inputStream)
{
    if (inputStream == nullptr)
    {
        GMX_THROW(FileIOError("Cannot open xplor file for reading."));
    }
    // determine which line has the ZYX string
    char tmp_char[80]; // auxiliary variable to copy tmp header information of the
                       // map
    int  row_ZYX = 0;
    do
    {
        GMX_IGNORE_RETURN_VALUE(fgets(tmp_char, 80, inputStream));
        row_ZYX++;
    }
    while (strncmp(tmp_char, "ZYX", 3) != 0);
    rewind(inputStream);

    // loop over lines before string ZYX
    decltype(extend) usedExtend;
    for (int row = 1; row <= row_ZYX; row++)
    {
        if (row < row_ZYX - 2)
        {
            // === lines 1-4: empty space or header remarks
            if (fgets(tmp_char, 80, inputStream) == nullptr)
            {
                GMX_THROW(FileIOError("Error reading header remarks."));
            }
        }
        // === Line 5 : Grid points ===
        if (row == row_ZYX - 2)
        {
            for (int dimension = XX; dimension <= ZZ; ++dimension)
            {
                if (fscanf(inputStream, "%8d", &(extend[dimension])) != 1)
                {
                    GMX_THROW(FileIOError("Error reading map extend."));
                }
                if (fscanf(inputStream, "%8d", &(crs_start[dimension])) != 1)
                {
                    GMX_THROW(FileIOError("Error reading column row section start."));
                }
                if (fscanf(inputStream, "%8d", &(crs_end[dimension])) != 1)
                {
                    GMX_THROW(FileIOError("Error reading column row section end."));
                }
            }
        }
        // === Line 6 : box size (in Angstrom) and angles (degree)===
        if (row == row_ZYX - 1)
        {
            for (int dimension = XX; dimension <= ZZ; ++dimension)
            {
                if (fscanf(inputStream, "%12f", &(cell_length[dimension])) != 1)
                {
                    GMX_THROW(FileIOError("Error reading cell length."));
                }
            }

            for (int dimension = XX; dimension <= ZZ; ++dimension)
            {
                if (fscanf(inputStream, "%12f", &(cell_angles[dimension])) != 1)
                {
                    GMX_THROW(FileIOError("Error reading cell angles."));
                }
            }
            // having read N,MIN,MAX, box and angle, it is possible to compute bin

            for (int i = 0; i < DIM; i++)
            {
                // set box in nm
                cell_length[i] /= nmToXplorUnits;
                // number of bins (working lattice)
                usedExtend[i] = crs_end[i] - crs_start[i] + 1;
            }
        }
        // === Line 7 : ZYX string ===
        if (row == row_ZYX)
        {
            if (fscanf(inputStream, "%s", tmp_char) != 1)
            {
                GMX_THROW(FileIOError("Error reading ZYX string."));
            }
        }
    }

    // ==== Line 8 onwards: the density map values =====
    for (int section = 0; section < usedExtend[ZZ]; section++) // loop over z
    {                                                          // Print the section value

        int sectionFromFile;
        if (fscanf(inputStream, "%d", &sectionFromFile) != 1)
        {
            GMX_THROW(FileIOError("Error reading section from file."));
        }

        if (section != sectionFromFile)
        {
            GMX_THROW(InconsistentInputError(
                              "File consistency error: Map section index " +
                              std::to_string(sectionFromFile) +
                              " read from file does not match current expected section " +
                              std::to_string(section) + "."));
        }
        for (int j = 0; j < usedExtend[YY]; j++)
        {
            for (int i = 0; i < usedExtend[XX]; i++)
            {
                float readDataItem;
                GMX_IGNORE_RETURN_VALUE(fscanf(inputStream, "%12f", &readDataItem));
                data.push_back(readDataItem);
            }
        }
    }

    int cross_check_line;
    if (fscanf(inputStream, "%d", &cross_check_line) != 1)
    {
        GMX_THROW(FileIOError("Error reading magic check number."));
    }
    if (cross_check_line != magicCheckNumber)
    {
        GMX_THROW(InconsistentInputError(
                          "File consistency error: Magic checking number is " +
                          std::to_string(cross_check_line) + " , but " +
                          std::to_string(magicCheckNumber) + " is expected."));
    }

    // === Last line : average and stdeviations ===
    if (fscanf(inputStream, "%f", &mean_value) != 1)
    {
        GMX_THROW(FileIOError("Error reading meam value."));
    }
    if (fscanf(inputStream, "%f", &rms_value) != 1)
    {
        GMX_THROW(FileIOError("Error reading root mean square deviation number."));
    }


}

void XplorData::write(FILE *outputStream) const
{
    // === line 1 : Empty line ===
    fprintf(outputStream, "\n");

    // === Lines 2-4 : Title and remarks ===
    fprintf(outputStream, "       2\n");

    IVec nBins = {crs_end[XX] - crs_start[XX] + 1, crs_end[YY] - crs_start[YY] + 1, crs_end[ZZ] - crs_start[ZZ] + 1 };
    fprintf(outputStream, "REMARKS Lattice Nbins: %d %d %d\n", nBins[XX],
            nBins[YY], nBins[ZZ]);

    fprintf(outputStream, "REMARKS Lattice resolution (nm)");
    for (int i = 0; i < DIM; i++)
    {
        fprintf(outputStream, " %f", cell_length[i] / extend[i]);
    }
    fprintf(outputStream, "\n");

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        fprintf(outputStream, "%8d", extend[dimension]);
        fprintf(outputStream, "%8d", crs_start[dimension]);
        fprintf(outputStream, "%8d", crs_end[dimension]);
    }
    fprintf(outputStream, "\n");

    // === Line 6 : box size (Angstrom) and angles ===

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        fprintf(outputStream, "%12.5E", cell_length[dimension] * nmToXplorUnits);
    }
    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        fprintf(outputStream, "%12.5E", cell_angles[dimension]);
    }
    fprintf(outputStream, "\n");
    // === Line 7 : string "ZYX"===
    fprintf(outputStream, "ZYX\n");

    // write density data slice by slice
    int       ncols          = 0;
    const int sliceSize      = nBins[XX] * nBins[YY];
    const int numberOfSlices = data.size()/sliceSize;
    auto      currentData    = data.begin();
    for (int sliceIndex = 0; sliceIndex < numberOfSlices; ++sliceIndex)
    {
        fprintf(outputStream, "       %d\n", sliceIndex);
        for (int sliceIndex = 0; sliceIndex < sliceSize; ++sliceIndex)
        {
            fprintf(outputStream, "%12.5E", *currentData);
            ++currentData;
            ++ncols;
            if (ncols == 6 || (sliceIndex == sliceSize - 1))
            {
                fprintf(outputStream, "\n");
                ncols = 0;
            }
        }
    }

    //  === Footer lines
    fprintf(outputStream, "   ");
    fprintf(outputStream, "%d\n", magicCheckNumber);
    fprintf(outputStream, "%13.5E%13.5E\n", mean_value, rms_value);
};

const float XplorData::nmToXplorUnits   = 10;
const int   XplorData::magicCheckNumber = -9999;
}
