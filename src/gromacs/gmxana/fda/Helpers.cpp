/*
 * Helpers.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "Helpers.h"
#include "gromacs/utility/fatalerror.h"
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>

namespace fda_analysis {

std::vector<double> parseScalarFileFormat(std::string const& filename,
	int nbParticles, int frame)
{
    int nbParticles2 = nbParticles * nbParticles;

    std::vector<double> array(nbParticles2, 0.0);
    std::ifstream pfFile(filename);
    if (!pfFile) gmx_fatal(FARGS, "Error opening file.");

    std::string line;
    bool foundFrame = false;
    int frameNumber;
    std::string tmp;
    int i, j;
    double value;

    while (getline(pfFile, line))
	{
		if (line.find("frame") != std::string::npos)
		{
			if (foundFrame) break;
			std::istringstream iss(line);
			iss >> tmp >> frameNumber;
			if (frameNumber != frame) continue;
	        foundFrame = true;
	        continue;
		}
		if (foundFrame) {
			std::istringstream iss(line);
            iss >> i >> j >> value;
            if (i >= nbParticles || j >= nbParticles)
            	gmx_fatal(FARGS, "Index is larger than dimension.");
            array[i*nbParticles+j] = value;
            array[j*nbParticles+i] = value;
		}
	}

    if (!foundFrame) gmx_fatal(FARGS, "Frame not found.");
    return array;
}

std::vector<double> getAveragedForcematrix(std::string const& filename,
    int nbParticles)
{
    std::ifstream pfFile(filename);
    if (!pfFile) gmx_fatal(FARGS, "Error opening file.");

    int nbParticles2 = nbParticles * nbParticles;
    std::vector<double> forcematrix(nbParticles2, 0.0);

    std::string line;
    int numberOfFrames = 0;
    std::string tmp;
    int i, j;
    double value;

    while (getline(pfFile, line))
	{
		if (line.find("frame") != std::string::npos) {
			++numberOfFrames;
			continue;
		}
		std::istringstream iss(line);
		iss >> i >> j >> value;
		if (i < 0 || i >= nbParticles || j < 0 || j >= nbParticles)
			gmx_fatal(FARGS, "Index error in getAveragedForcematrix.");
		forcematrix[i*nbParticles+j] = value;
	}

    if (!numberOfFrames) gmx_fatal(FARGS, "No frame found.");
    for (auto & elem : forcematrix) elem /= numberOfFrames;

    return forcematrix;
}

std::vector<real> readStress(std::string const& filename, int& nbFrames,
    int& nbParticles)
{
    std::ifstream file(filename);
    if (!file) gmx_fatal(FARGS, "Error opening file %s.", filename.c_str());

    std::string line;
    std::vector<real> stress;
    nbFrames = 0;
    nbParticles = 0;

    if (getline(file, line))
    {
		std::stringstream ss(line);
		real value;
		while(ss >> value) {
		   stress.push_back(value);
			++nbParticles;
		}
		++nbFrames;
    }

    while (getline(file, line))
    {
		std::stringstream ss(line);
		real value;
		while(ss >> value) {
		   stress.push_back(value);
		}
		++nbFrames;
    }

    return stress;
}

size_t getNumberOfFrames(std::string const& filename)
{
    std::ifstream pfFile(filename);
    if (!pfFile) gmx_fatal(FARGS, "Error opening file.");

    std::string line;
    int numberOfFrames = 0;

    while (getline(pfFile, line))
	{
		if (line.find("frame") != std::string::npos) {
			++numberOfFrames;
		}
	}

    return numberOfFrames;
}

size_t getMaxIndexSecondColumnFirstFrame(std::string const& filename)
{
    std::ifstream pfFile(filename);
    if (!pfFile) gmx_fatal(FARGS, "Error opening file.");

    std::string line;
    int i, j, maxIndex = 0;

    getline(pfFile, line); // skip first line

    while (getline(pfFile, line))
	{
		if (line.find("frame") != std::string::npos) break;

		std::stringstream ss(line);
		ss >> i >> j;

		if (j > maxIndex) maxIndex = j;
	}

    return maxIndex;
}

bool isInteger(std::string const& str)
{
	for (auto c : str) if (!std::isdigit(c)) return false;
	return true;
}

bool hasExtension(std::string const& str, std::string const& extension)
{
	return str.size() > extension.size() && str.compare(str.size() - extension.size(), extension.size(), extension) == 0;
}

FrameType getFrameTypeAndSkipValue(std::string const& frameString, int& value)
{
	FrameType frameType;
	value = 1;

	std::istringstream iss(frameString);
	std::istream_iterator<std::string> iterCur(iss), iterEnd;
	if (iterCur == iterEnd) gmx_fatal(FARGS, "Error in frame option.");

	std::string firstToken = *iterCur;
    ++iterCur;
	if (iterCur == iterEnd) {
		if (isInteger(firstToken)) {
			frameType = SINGLE;
			value = std::atoi(firstToken.c_str());
		}
		else frameType = EnumParser<FrameType>()(firstToken);
	} else {
	    frameType = EnumParser<FrameType>()(firstToken);
		value = std::atoi(iterCur->c_str());

		// Abort if more than two tokens are found
	    ++iterCur;
		if (iterCur != iterEnd) gmx_fatal(FARGS, "Error in frame option.");
	}

	if (frameType == ALL && value != 1) gmx_fatal(FARGS, "Frame type \"all\" does not expect a number.");

	return frameType;
}

} // namespace fda_analysis
