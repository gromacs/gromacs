/*
 * PDB.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "PDB.h"
#include "gromacs/utility/fatalerror.h"
#include "Index.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

namespace fda_analysis {

PDB::PDB(std::string const& pdbFilename)
{
    std::ifstream pdbFile(pdbFilename);
    if (!pdbFile) gmx_fatal(FARGS, "Error opening pdb file.");

    std::string line;
    while (getline(pdbFile, line))
	{
        if (line.substr(0,4) == "ATOM" or line.substr(0,6) == "HETATM")
        {
			indices_.push_back(stoi(line.substr(6,5)));
			Coordinate coord;
			coord[0] = stod(line.substr(30,8));
			coord[1] = stod(line.substr(38,8));
			coord[2] = stod(line.substr(46,8));
			coordinates_.push_back(coord);
		}
	}
}

PDB::PDB(std::string const& pdbFilename, std::vector<int> groupAtoms)
{
    std::ifstream pdbFile(pdbFilename);
    if (!pdbFile) gmx_fatal(FARGS, "Error opening pdb file.");

    std::sort(groupAtoms.begin(), groupAtoms.end());

    std::string line;
    std::vector<int>::const_iterator iterGroupAtoms = groupAtoms.begin();
    while (getline(pdbFile, line))
	{
        if (line.substr(0,4) == "ATOM" or line.substr(0,6) == "HETATM")
        {
        	if (stoi(line.substr(6,5)) == *iterGroupAtoms + 1)
        	{
        		indices_.push_back(stoi(line.substr(6,5)));
        		Coordinate coord;
        		coord[0] = stod(line.substr(30,8));
				coord[1] = stod(line.substr(38,8));
				coord[2] = stod(line.substr(46,8));
        		coordinates_.push_back(coord);
				++iterGroupAtoms;
        	}
		}
	}
}

void PDB::writePaths(std::string const& filename, std::vector< std::vector<int> > const& shortestPaths,
    std::vector<double> const& forceMatrix, bool append) const
{
    std::ofstream pdb;
    if (append) pdb.open(filename, std::ofstream::app);
    else pdb.open(filename);
    if (!pdb) gmx_fatal(FARGS, "Error opening pdb file.");

    std::stringstream connections;
    int i, j;
    int numAtom = 1;
    int numNetwork = 0;
    std::vector<int>::const_iterator i1, i2;
	int dim = std::sqrt(forceMatrix.size());

    real currentForce;
    bool valueToLargeForPDB = false;

    for (auto const& path : shortestPaths)
    {
		for (size_t n = 0; n < path.size() - 1; ++n)
		{
			i = path[n];
			j = path[n+1];

            currentForce = forceMatrix[i*dim+j];
            if (currentForce > 999.99) valueToLargeForPDB = true;

			writeAtomToPDB(pdb, numAtom, indices_[i], coordinates_[i], currentForce, numNetwork);
			writeAtomToPDB(pdb, numAtom + 1, indices_[j], coordinates_[j], currentForce, numNetwork);

			connections << "CONECT" << std::setw(5) << numAtom << std::setw(5) << numAtom + 1 << std::endl;

			numAtom += 2;
		}
		++numNetwork;
    }

    if (virginValueToLargeForPDB and valueToLargeForPDB) {
        gmx_warning("Force values larger than 999.99 are detected. Therefore, the general PDB format of the b-factor column of Real(6.2) is broken. "
                    "It is tested that it works for Pymol and VMD, but it is not guaranteed that it will work for other visualization programs.");
        virginValueToLargeForPDB = false;
    }

    pdb << connections.str()
        << "ENDMDL" << std::endl;
}

void PDB::updateCoordinates(const rvec x[])
{
	checkSanity();

	for (size_t i = 0; i != indices_.size(); ++i)
	{
		coordinates_[i][0] = x[indices_[i] - 1][0] * 10.0;
		coordinates_[i][1] = x[indices_[i] - 1][1] * 10.0;
		coordinates_[i][2] = x[indices_[i] - 1][2] * 10.0;
	}
}

void PDB::print() const
{
	checkSanity();

	std::cerr << "\nPDB status:" << std::endl;
	for (size_t i = 0; i != indices_.size(); ++i)
		std::cerr << indices_[i] << " " << coordinates_[i][0] << " " << coordinates_[i][1] << " " << coordinates_[i][2] << std::endl;
}

void PDB::writeAtomToPDB(std::ofstream& os, int num, int index, Coordinate const& coord, double force,
    int numNetwork) const
{
	os << "ATOM"
	   << std::setw(7) << num
	   << std::setw(15) << index
       << std::fixed << std::setprecision(3) << std::setw(12) << coord[0]
       << std::setw(8) << coord[1]
       << std::setw(8) << coord[2]
       << std::setprecision(2) << std::setw(6) << 1.0
       << std::setprecision(2) << std::setw(6) << force
       << std::setw(8) << colors[numNetwork % 32]
       << std::endl;
}

void PDB::checkSanity() const
{
	 if (indices_.size() != coordinates_.size()) gmx_fatal(FARGS, "PDB sanity check failed.");
}

std::vector<std::string> PDB::colors =
{
	"AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ",
    "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT",
    "AU", "AV", "AW", "AX", "AY", "AZ", "BA", "BB", "BC", "BD",
    "BE", "BF"
};

bool PDB::virginValueToLargeForPDB = true;

} // namespace fda_analysis
