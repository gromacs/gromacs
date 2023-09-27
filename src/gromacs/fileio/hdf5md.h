/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
/* This file was inspired by ch5md by Pierre de Buyl (BSD license). */

#ifdef GMX_USE_HDF5

#ifndef GMX_FILEIO_HDF5MD_H
#define GMX_FILEIO_HDF5MD_H

#include <string>

#include <h5xx/h5xx.hpp>

class GmxHdf5MdElement
{
private:
    h5xx::attribute group;
    h5xx::attribute step;
    h5xx::attribute time;
    h5xx::attribute value;
    h5xx::datatype datatype;
    int isTime;
    int currentStep;
    GmxHdf5MdParticlesGroup *particlesGroup;
public:
    GmxHdf5MdElement();
    ~GmxHdf5MdElement();
    int append(void *data, int step, double time);
}

class GmxHdf5MdParticlesGroup
{
private:
    h5xx::attribute group;
    GmxHdf5MdElement position;
    h5xx::attribute box;
    GmxHdf5MdElement boxEdges;
    GmxHdf5MdElement image;
    GmxHdf5MdElement velocity;
    GmxHdf5MdElement force;
/*    GmxHdf5MdElement mass;
    GmxHdf5MdElement species;
    GmxHdf5MdElement id;
    GmxHdf5MdElement charge;*/
    int localSizeMax;
public:
    GmxHdf5MdParticlesGroup();
    ~GmxHdf5MdParticlesGroup();
    int createBox(int dim, char *boundary[], bool isTime, double value[], GmxHdf5MdElement *link);
}
class GmxHdf5MdFile
{
private:
    h5xx::attribute id;
    h5xx::attribute particles;
    h5xx::attribute observables;
    h5xx::attribute parameters;
public:
    GmxHdf5MdFile();
    GmxHdf5MdFile(std::string fileName);
    ~GmxHdf5MdFile();
    int createParticlesGroup(std::string name);
//     GmxHdf5MdElement createTimeData(hid_t loc, std::string name, int rank, int dims[], hid_t datatype, GmxHdf5MdElement *link);
//     GmxHdf5MdElement createFixedDataSimple(hid_t loc, std::string name, int rank, int dims[], hid_t datatype, void *data);
//     GmxHdf5MdElement createFixedDataScalar(hid_t loc, std::string name, hid_t datatype, void *data);
//     int writeStringAttribute(hid_t loc, std::string objectName, std::string attributeName, std::string value);
};



#endif // GMX_FILEIO_HDF5MD_H

#endif // GMX_USE_HDF5
