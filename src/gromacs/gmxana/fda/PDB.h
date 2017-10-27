/*
 * PDB.h
 *
 *  Created on: Feb 23, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef PDB_H_
#define PDB_H_

#include <string>
#include <vector>
#include <boost/array.hpp> // back-compatibility to gcc-4.7.2
#include "gromacs/math/vectypes.h"

namespace fda_analysis {

/**
 * Simplified class for PDB file format for fda analytics.
 */
class PDB
{
public:

	// Plain C-like arrays can not be used within STL containers, since they are neither
	// copy-constructible nor assignable.
    typedef boost::array<real, 3> Coordinate;

	PDB(std::string const& pdbFilename);

	/// Only the atoms listed in @groupAtoms will be stored.
	PDB(std::string const& pdbFilename, std::vector<int> groupAtoms);

    void writePaths(std::string const& filename, std::vector< std::vector<int> > const& shortestPaths,
        std::vector<double> const& forceMatrix, bool append) const;

    /// Update with coordinates of trajectory file.
    /// Values will be converted from nm into Angstrom.
    void updateCoordinates(const rvec x[]);

    std::vector<Coordinate> getCoordinates() const { return coordinates_; }

    /// Print current state
    void print() const;

private:

	void writeAtomToPDB(std::ofstream& os, int num, int index, Coordinate const& coord, double force,
        int numNetwork) const;

	/// Abort if not sanity.
	void checkSanity() const;

	std::vector<uint32_t> indices_;

	std::vector<Coordinate> coordinates_;

    static std::vector<std::string> colors;

    /// Flag for printing warning only one time per run.
    static bool virginValueToLargeForPDB;

};

} // namespace fda_analysis

#endif /* PDB_H_ */
