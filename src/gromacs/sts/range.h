#ifndef STS_RANGE_H
#define STS_RANGE_H

#include <cassert>
#include <string>

// TODO: Support negative ratios in all operations
class Ratio {
public:
    Ratio()               : nom_(0), denom_(1) {}
    Ratio(int n, int d=1) : nom_(n), denom_(d) {}
    operator int() const { return nom_/denom_; }
    Ratio reduce() {
        int gcd = Ratio::gcd(nom_, denom_);
        nom_ /= gcd;
        denom_ /= gcd;
        return *this;
    }
    std::string toString() const {
        return std::to_string(nom_) + "/" + std::to_string(denom_);
    }
    static int gcd(int a, int b) {
        assert(a >= 0 && b >= 0);
        if (a == 0) return b;
        if (b == 0) return a;
        return gcd(b,a%b);
    }
private:
    int nom_, denom_;
    friend Ratio operator*(Ratio, int);
    friend Ratio operator*(Ratio, Ratio);
    friend Ratio operator+(Ratio, Ratio);
    friend Ratio operator+=(Ratio &, Ratio);
    friend Ratio operator-(Ratio, Ratio);
};

template<class T>
class Range {
public:
    Range() {}
    Range(T s, T e) : start(s), end(e) {} 
    explicit Range(T e) : start(0), end(e) {}
    template<class R>
    operator Range<R>() const { return Range<R>(start, end); }
    Range<T> subset(Range<Ratio> p) const {
        Range<T> r = p * (end - start);
        return r + start;
    }
    std::string toString() const {
        return start.toString() + " " + end.toString();
    }
    T len() const {
        return end - start;
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
