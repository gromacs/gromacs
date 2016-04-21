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

#ifndef STS_RANGE_H
#define STS_RANGE_H

class Ratio {
public:
    Ratio(int n, int d=1) : nom_(n), denom_(d) {} ;
    operator int() { return nom_/denom_; } 
private:
    int nom_, denom_;
    friend Ratio operator*(Ratio, int);
};

template<class T>
class Range {
public:
    Range(T s, T e) : start(s), end(e) {} 
    explicit Range(T e) : start(0), end(e) {}
    template<class R>
    operator Range<R>() { return Range<R>(start, end); }
    Range<T> subset(Range<Ratio> p) {
        Range<T> r = p * (end - start);
        return r + start;
    }
    T start, end;
};

template<class T>
Range<T> operator*(Range<T> r, int n) { return Range<T>(r.start*n, r.end*n); }
template<class T>
Range<T> operator*(int n, Range<T> r) { return r*n; }

template<class T>
Range<T> operator+(Range<T> r, int n) { return Range<T>(r.start+n, r.end+n); }
template<class T>
Range<T> operator+(int n, Range<T> r) { return r+n; }

#endif // STS_RANGE_H
