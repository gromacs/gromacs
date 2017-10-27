/*
 * TextSplitter.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <fstream>
#include "TextSplitter.h"
#include "gromacs/utility/fatalerror.h"

TextSplitter::TextSplitter(std::string const& filename, bool ignoreComments)
{
	size_t first, last;
	std::string line, token;
	std::ifstream is(filename);
	if (!is) gmx_fatal(FARGS, "Error reading file %s\n", filename.c_str());

    while (getline(is, line))
    {
    	// trim whitespaces
    	first = line.find_first_not_of(' ');
    	last  = line.find_last_not_of(' ');
    	line = line.substr(first, last - first + 1);

    	// ignore comments
		if (ignoreComments and line.compare(0, 2, "/*") == 0 and line.compare(line.size()-2, 2, "*/") == 0) continue;

		while (is >> token)
		{
			if (ignoreComments and token.compare(0, 2, "//") == 0) break;
			if (isType<int>(token)) integers_.push_back(std::stoi(token));
			else if (isType<double>(token)) reals_.push_back(std::stod(token));
			else strings_.push_back(token);
		}
    }
}
