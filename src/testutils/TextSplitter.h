/*
 * TextSplitter.h
 *
 *  Created on: Sep 1, 2014
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef TEXTSPLITTER_H_
#define TEXTSPLITTER_H_

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

/**
 * Splitting an ASCII file into arrays of integer, floating point numbers and strings
 * to compare if the two files are logically identical. Logically identical means that
 * whitespaces will be ignored, and the difference of two floating point numbers should
 * be within a given accuracy.
 *
 * The comments /\star ... \star/ and // ... will be ignore by default.
 */
class TextSplitter
{
public:

	TextSplitter(std::string const& filename, bool ignoreComments = true);

private:

	template <class Comparer>
	friend bool equal(TextSplitter const& t1, TextSplitter const& t2, Comparer const& comparer);

	std::vector<int> integers_;
	std::vector<std::string> strings_;
	std::vector<double> reals_;

};

template <typename T>
bool isType(std::string const& s)
{
    std::istringstream iss(s);
    T x;
    return iss >> x && !iss.ignore();
}

template <class Comparer>
bool equal(TextSplitter const& t1, TextSplitter const& t2, Comparer const& comparer)
{
	if (t1.integers_.size() != t2.integers_.size()) {
		std::cout << "TextSplitter: wrong integers size" << std::endl;
		return false;
	}
	if (t1.strings_.size() != t2.strings_.size()) {
		std::cout << "TextSplitter: wrong strings size" << std::endl;
		return false;
	}
	if (t1.reals_.size() != t2.reals_.size()) {
		std::cout << "TextSplitter: wrong reals size" << std::endl;
		return false;
	}

	for (size_t i(0); i != t1.integers_.size(); ++i) if (t1.integers_[i] != t2.integers_[i]) {
		std::cout << "TextSplitter: wrong integer " << t1.integers_[i] << " " << t2.integers_[i] << std::endl;
		return false;
	}
	for (size_t i(0); i != t1.strings_.size(); ++i) if (t1.strings_[i] != t2.strings_[i]) {
		std::cout << "TextSplitter: wrong string " << t1.strings_[i] << " " << t2.strings_[i] << std::endl;
		return false;
	}
	for (size_t i(0); i != t1.reals_.size(); ++i) {
		if (!comparer(t1.reals_[i], t2.reals_[i])) {
		    std::cout << std::scientific << std::setprecision(20) << "TextSplitter: wrong real " << t1.reals_[i] << " " << t2.reals_[i] << std::endl;
		    return false;
	    }
	}

    return true;
}

#endif /* TEXTSPLITTER_H_ */
