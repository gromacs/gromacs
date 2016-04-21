#include "range.h"

Ratio operator*(Ratio r, int64_t n) { return Ratio(r.nom_*n,r.denom_); }
Ratio operator*(int64_t n, Ratio r) { return r*n; }
Ratio operator*(Ratio r1, Ratio r2) {
   return Ratio(r1.nom_ * r2.nom_ , r1.denom_ * r2.denom_).reduce();
}
Ratio operator+=(Ratio &r1, Ratio r2) {
    r1.nom_ = r1.nom_ * r2.denom_ + r2.nom_ * r1.denom_;
    r1.denom_ *= r2.denom_;
    return r1.reduce();
}
Ratio operator+(Ratio r1, Ratio r2) {
    return r1 += r2;
}
Ratio operator-(Ratio r1, Ratio r2) {
    Ratio r3 = { r1.nom_ * r2.denom_ - r2.nom_ * r1.denom_ , r1.denom_ * r2.denom_ };
    return r3.reduce();
}
