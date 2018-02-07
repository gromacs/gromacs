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
/*TODO: change to ! \file
 * \brief
 * Implements utility functions and classes for Fca.
 */
#include "gmxpre.h"

#include "utils.h"

#include <cstdio>

#include <vector>

#include "gromacs/fileio/confio.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

namespace fca_utils
{

/*TODO: change to ! \brief
 * Read coordinates out of STX file.
 */
int InputFunctions::read_conffile(const char* confin, char* title, rvec* x[])
{
    int     natoms;
    t_atoms confat;
    matrix  box;

    printf("read coordnumber from file %s\n", confin);
    get_stx_coordnum(confin, &natoms);
    printf("number of coordinates in file %d\n", natoms);
    /*  if (natoms != ncoords){
        fatal_error(0,"number of coordinates in coordinate file (%s, %d)\n"
        "does not match topology (= %d)",
        confin,natoms,ncoords);
        }
     */
    /* make space for coordinates and velocities */
    init_t_atoms(&confat, natoms, FALSE);
    printf("init_t\n");
    snew(*x, natoms);
    read_stx_conf(confin, nullptr, &title, &confat, *x, nullptr, nullptr /*epbcXYZ*/, box);
    return natoms;
}

/*TODO: change to ! \brief
 * Copies coordinates from x to edx which are given in index
 */
void InputFunctions::filter2x(rvec* x, int nindex, int index[], int ngro,
                              int igro[], rvec* xori, const char* structure)
{
    int pos, i;
    for (i = 0; i < nindex; i++)
    {
        for (pos = 0; pos < ngro - 1 && igro[pos] != index[i]; ++pos)
        {
        }       /*search element in igro*/
        if (igro[pos] != index[i])
        {
            fprintf(stdout, "Couldn't find atom with index %d in structure %s",
                    index[i], structure);
            exit(1);
        }
        copy_rvec(xori[pos], x[i]);
    }
}

bool InputFunctions::get_structure(t_atoms* atoms, const char* IndexFile,
                                   const char* StructureFile, rvec* x, int natoms,
                                   int index[])
{
    int     * igro; /*index corresponding to target or origin structure*/
    int       ngro;
    rvec    * xtar;
    char      title[STRLEN];
    char    * grpname;

    const int ntar = read_conffile(StructureFile, title, &xtar);
    get_index(atoms, IndexFile, 1, &ngro, &igro, &grpname);
    filter2x(x, natoms, index, ngro, igro, xtar, StructureFile);
    if (ngro != ntar)
    {
        return false;
    }
    return true;
}

int InputFunctions::count_number_columns(FILE* fp)
{
    constexpr size_t BUFFSIZE = 20000;
    char             buffer[BUFFSIZE];
    if (fgets(buffer, BUFFSIZE, fp) == nullptr)
    {
        gmx_fatal(FARGS, "Error reading buffer\n");
    }

    // skip first spaces
    int idxRead = 0;
    while (isspace(buffer[idxRead]) && buffer[idxRead])
    {
        ++idxRead;
    }       //load_frame_step spaces

    int cptColumns   = 0;
    while (buffer[idxRead])
    {
        // if it is a char count it as a column
        if (!isspace(buffer[idxRead]))
        {
            ++cptColumns;
        }
        // skip following char
        while (!isspace(buffer[idxRead]) && buffer[idxRead])
        {
            ++idxRead;
        }       //load_frame_step rest of number
        // skip following space
        while (isspace(buffer[idxRead]) && buffer[idxRead])
        {
            ++idxRead;
        }       //load_frame_step spaces
    }
    fprintf(stdout, "analyze line %s\n", buffer);
    fprintf(stdout, "found %d cols\n", cptColumns);
    return cptColumns;
}

std::vector< std::unique_ptr< real[] > > InputFunctions::read_fca_proj(const char pos_file[], int* dim)
{
    FILE* fp   = gmx_ffopen(pos_file, "r");
    (*dim) = count_number_columns(fp);
    rewind(fp);
    fprintf(stdout, "read vectors of dimensionality %d\n", *dim);

    std::vector< std::unique_ptr< real[] > > projx;
    while (true)
    {
        std::unique_ptr< real[] > nextFrame(new real[*dim]);

        int cptRead = 0;
        while (cptRead < (*dim) && fca_utils::fscanf_floating_point(fp, &(nextFrame[cptRead])) > 0)
        {
            cptRead += 1;
        }

        if (cptRead != (*dim))
        {
            break;
        }

        projx.emplace_back(std::move(nextFrame));
    }
    gmx_ffclose(fp);
    return projx;
}

}
