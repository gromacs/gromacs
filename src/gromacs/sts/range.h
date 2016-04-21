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
