#ifndef SZ3_BIOMD_DECOMPOSITION_HPP
#define SZ3_BIOMD_DECOMPOSITION_HPP

// ALGO_BIOMD: molecular-dynamics coordinates {frames, atoms, 3} (nm, absolute bound) as streams of integer symbols,
// one after the other as [count, symbols]. Every coordinate goes on a lattice,
//   q = round(x * inv_step), |x - q step| <= eb,
// and all prediction is integer arithmetic on that lattice, which the decoder repeats exactly.
//  * rigid water (O, H1, H2): H1 on the sphere |H1 - O| = r, H2 on the circle that r and the H-H distance leave
//  * other atoms: from their bond partner, one of the previous MAX_BOND_OFFSET atoms (intra a delta, else on the sphere
//    of the previous frame's bond length), else a delta
//  * frame 0 intra; for the frames after it water O's from the previous frame, and per other atom group intra or from
//    the previous frame, whichever costs less on a sample of frame 1
// The layout (which atoms are water, the bond of each other atom) is found on a chunk's first frame and stored with it.
// In a chunk of one frame without water and with under a quarter of its atoms bonded (coarse-grained runs), unbonded
// atoms are their points in their box, like intra water O's.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {
namespace biomd {

// A bond partner is searched among the previous MAX_BOND_OFFSET atoms: topologies list an atom shortly after the atom
// it is bonded to.
constexpr int MAX_BOND_OFFSET = 4;
// Deltas from here on (lattice units) are an escape symbol and a raw byte, which keeps the Huffman alphabet small.
constexpr uint32_t ESCAPE = 2048;
// |q| <= 2^28: displacements stay within 2^29, and every symbol within 32 bits
constexpr double MAX_LATTICE = double(1 << 28);
// The lattice spans |x| < B, the smallest power of two above 2^LATTICE_SPAN_BITS eb (for float 2048 nm at 5e-4 nm, 128
// nm at 5e-5 nm). A larger span costs ratio, as the step is 2 eb less the float ulp at B (choose_step()); double's ulp
// is small enough for B to reach MAX_LATTICE.
template <class T>
constexpr int LATTICE_SPAN_BITS = sizeof(T) == 4 ? 21 : 27;

// atom groups, each with its own predictor mode for the frames after the first (water O's: the previous frame)
enum { G_WATER_O, G_WATER_H, G_BONDED, G_UNBONDED, NUM_GROUPS };
// symbol streams: the layout (gaps between water O's, the bond of each other atom); water O from the previous frame;
// water H1 on the sphere around O (face, kept coordinates, radial residual); water H2 on the circle (coordinate and
// side, residuals); bonded atoms intra as deltas from their partner, else on the sphere around it (kept coordinates,
// radial residual); unbonded atoms as deltas
enum {
    S_WATER_GAP,
    S_BOND_REF,
    S_WATER_O,
    S_WATER_H1_FACE,
    S_WATER_H1_KEPT,
    S_WATER_H1_RADIAL,
    S_WATER_H2_AXIS,
    S_WATER_H2_RESIDUAL,
    S_BOND_DELTA,
    S_BOND_KEPT,
    S_BOND_RADIAL,
    S_UNBONDED,
    NUM_STREAMS
};
// the group whose atoms write the stream (-1: the layout streams, one symbol per atom at most), and the most symbols
// one of them (water: one molecule) writes per frame; the two size the streams' buffers
inline constexpr int STREAM_GROUP[NUM_STREAMS] = {-1,        -1,        G_WATER_O, G_WATER_H, G_WATER_H, G_WATER_H,
                                                  G_WATER_H, G_WATER_H, G_BONDED,  G_BONDED,  G_BONDED,  G_UNBONDED};
inline constexpr int MAX_SYMBOLS_PER_UNIT[NUM_STREAMS] = {1, 1, 3, 1, 2, 1, 1, 2, 3, 2, 1, 3};

// Values BIOMD does not code (NaN or Inf outside trailing fill, coordinates beyond the lattice); SZ_compress_bioMD
// stores the chunk losslessly.
struct Fallback : std::runtime_error {
    using std::runtime_error::runtime_error;
};

inline int64_t round_half_away(double y) { return int64_t(y + std::copysign(0.5, y)); }
inline int64_t isqrt_round(int64_t n) { return n <= 0 ? 0 : int64_t(std::sqrt(double(n)) + 0.5); }
// d = atom i - atom j
inline void displacement(const int32_t *q, size_t i, size_t j, int64_t d[3]) {
    for (int c = 0; c < 3; c++) d[c] = int64_t(q[3 * i + c]) - q[3 * j + c];
}

// The largest axis of p, dropped by the sphere code, and the other two, kept.
ALWAYS_INLINE void split_axes(const int64_t p[3], int &drop, int &keep0, int &keep1) {
    drop = std::llabs(p[1]) > std::llabs(p[0]) ? 1 : 0;
    if (std::llabs(p[2]) > std::llabs(p[drop])) drop = 2;
    keep0 = (drop + 1) % 3;
    keep1 = (drop + 2) % 3;
}

// ------------------------------------------------------------------------------------------------ layout
enum : uint8_t { K_WATER_O, K_WATER_H, K_OTHER };
struct Layout {
    std::vector<uint8_t> kind;          // per atom: K_WATER_O (then its H at i + 1, i + 2), K_WATER_H, K_OTHER
    std::vector<uint8_t> bond;          // other atoms: 0 = none, else the offset of the partner
    double water_oh = 0, water_hh = 0;  // O-H and H-H distances of the water (nm)
};

// Rigid water on frame x: kind, water_oh, water_hh. Without water, water_oh stays 0 and detect_layout() clears kind.
template <class T>
void detect_water(const T *x, size_t atoms, Layout &layout) {
    layout.kind.assign(atoms, K_OTHER);
    auto dist2 = [x](size_t a, size_t b) {
        float dx = float(x[3 * a] - x[3 * b]), dy = float(x[3 * a + 1] - x[3 * b + 1]),
              dz = float(x[3 * a + 2] - x[3 * b + 2]);
        return nofma(dx * dx) + nofma(dy * dy) + nofma(dz * dz);
    };
    // Candidates: O-H in [0.08, 0.125] nm, H-H in [0.13, 0.2] nm, probed at PROBES places across the frame. Rigid
    // water gives sharp O-H / H-H peaks, flexible CH2/NH2 groups broad ones.
    constexpr float OH2_MIN = 0.0064f, OH2_MAX = 0.015625f, HH2_MIN = 0.0169f, HH2_MAX = 0.04f;  // nm^2
    constexpr size_t PROBES = 300, MIN_COUNT = 16;  // fewer candidates, or waters, than MIN_COUNT: no water
    std::vector<float> cand_oh, cand_hh;
    const size_t probe = std::max<size_t>(1, atoms / PROBES);
    for (size_t i = 0; i + 2 < atoms; i += probe)
        // five starts: a probe can land on any site of a 3- to 5-site water, and one of them is then its O
        for (size_t o = 0; o < 5 && i + o + 2 < atoms; o++) {
            const size_t k = i + o;
            const float a = dist2(k, k + 1), b = dist2(k, k + 2), c = dist2(k + 1, k + 2);
            if (a > OH2_MIN && a < OH2_MAX && b > OH2_MIN && b < OH2_MAX && c > HH2_MIN && c < HH2_MAX) {
                cand_hh.push_back(std::sqrt(c));
                cand_oh.push_back(std::sqrt(a));
                cand_oh.push_back(std::sqrt(b));
                break;
            }
        }
    if (cand_hh.size() < MIN_COUNT) return;
    // tolerances (nm), the same at every bound: rigid water keeps its geometry far tighter than `tight`, and `loose`
    // stays below the O-H / C-H difference (> 0.01 nm) that separates water from CH2/NH2
    constexpr double tight = 0.001, tol = 0.002, loose = 0.01;
    auto peak = [](std::vector<float> v) {     // the median of the heaviest window of PEAK_WINDOW
        constexpr float PEAK_WINDOW = 0.004f;  // nm
        std::sort(v.begin(), v.end());
        size_t best = 0, n = 0;
        for (size_t lo = 0, hi = 0; hi < v.size(); hi++) {
            while (v[hi] - v[lo] > PEAK_WINDOW) lo++;
            if (hi - lo + 1 > n) {
                n = hi - lo + 1;
                best = (lo + hi) / 2;
            }
        }
        return double(v[best]);
    };
    const double hh = peak(cand_hh);
    std::vector<float> oh_of_hh;
    for (size_t c = 0; c < cand_hh.size(); c++)
        if (std::fabs(cand_hh[c] - hh) < tol)
            oh_of_hh.push_back(cand_oh[2 * c]), oh_of_hh.push_back(cand_oh[2 * c + 1]);
    const double oh = peak(oh_of_hh);
    // rigidity test: rigid water puts nearly every nearby candidate within `tight` of both peaks, flexible water and
    // CH2/NH2 groups under half of them
    size_t n_tight = 0, n_loose = 0;
    for (size_t c = 0; c < cand_hh.size(); c++) {
        const double dh = std::fabs(cand_hh[c] - hh), da = std::fabs(cand_oh[2 * c] - oh),
                     db = std::fabs(cand_oh[2 * c + 1] - oh);
        if (dh < loose && da < loose && db < loose) {
            n_loose++;
            n_tight += dh < tight && da < tight && db < tight;
        }
    }
    if (n_loose < MIN_COUNT || n_tight * 10 < n_loose * 7) return;
    const float oh2_lo = float((oh - tol) * (oh - tol)), oh2_hi = float((oh + tol) * (oh + tol));
    const float hh2_lo = float((hh - tol) * (hh - tol)), hh2_hi = float((hh + tol) * (hh + tol));
    size_t waters = 0;
    for (size_t i = 0; i + 2 < atoms;) {
        float a = dist2(i, i + 1), b, c;
        if (a > oh2_lo && a < oh2_hi && (b = dist2(i, i + 2)) > oh2_lo && b < oh2_hi &&
            (c = dist2(i + 1, i + 2)) > hh2_lo && c < hh2_hi) {
            layout.kind[i] = K_WATER_O;
            layout.kind[i + 1] = layout.kind[i + 2] = K_WATER_H;
            waters++;
            i += 3;
        } else {
            i++;
        }
    }
    if (waters < MIN_COUNT) return;
    layout.water_oh = oh;
    layout.water_hh = hh;
}

// Bonds of the other atoms on frame x: each takes the nearest of its previous MAX_BOND_OFFSET atoms if that is within
// 0.01 .. 0.2 nm (covalent bonds, and the virtual site of 4-site water at 0.015 nm from its O; coarse-grained bonds are
// longer). Under one bonded atom in 16 among other atoms probed across the frame, the S_BOND_REF symbols (one per
// atom) would cost more than the bonds save, and none is kept.
template <class T>
void detect_bonds(const T *x, size_t atoms, Layout &layout) {
    auto partner = [&](size_t i) {  // the offset of the bond partner of other atom i, or 0
        float best = 1e30f;
        uint8_t offset = 0;  // of the nearest previous atom
        for (unsigned o = 1; o <= unsigned(MAX_BOND_OFFSET) && o <= i; o++) {
            const float dx = float(x[3 * i] - x[3 * (i - o)]), dy = float(x[3 * i + 1] - x[3 * (i - o) + 1]),
                        dz = float(x[3 * i + 2] - x[3 * (i - o) + 2]),
                        v = nofma(dx * dx) + nofma(dy * dy) + nofma(dz * dz);
            offset = v < best ? uint8_t(o) : offset;
            best = std::min(best, v);
        }
        return best > 0.0001f && best < 0.04f ? offset : uint8_t(0);
    };
    layout.bond.assign(atoms, 0);
    size_t probed = 0, bonded = 0;
    for (size_t i = 1; i < atoms; i += std::max<size_t>(1, atoms / 256))
        if (layout.kind[i] == K_OTHER) probed++, bonded += partner(i) != 0;
    if (bonded * 16 < probed) return;
    for (size_t i = 1; i < atoms; i++)
        if (layout.kind[i] == K_OTHER) layout.bond[i] = partner(i);
}

// ------------------------------------------------------------------------------------------------ symbols
// the position in each stream, and raw bits (LSB first)
struct SymbolWriter {
    int *cursor[NUM_STREAMS];
    BitAppender *raw;
    ALWAYS_INLINE void put(int s, uint32_t v) { *cursor[s]++ = int(v); }
    ALWAYS_INLINE void put_bits(uint64_t v, int bits) { raw->put(v, bits); }  // bits <= 56
};
// bits of a sample, for the predictor choice
struct CostEstimate {
    double bits = 0;
    inline void put(int, uint32_t v) { bits += std::log2(double(v) + 1.0) + 1; }
    inline void put_bits(uint64_t, int n) { bits += n; }
};
// a cursor per stream, and one on the raw bits
struct SymbolReader {
    const int *cursor[NUM_STREAMS], *end[NUM_STREAMS];
    BitConsumer raw{nullptr, nullptr};
    ALWAYS_INLINE uint32_t get(int s) {
        if (cursor[s] == end[s]) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        return uint32_t(*cursor[s]++);
    }
    ALWAYS_INLINE uint64_t get_bits(int bits) { return raw.get(bits); }  // bits <= 56
};

// v as one symbol below ESCAPE, else as ESCAPE + v / 256 and its low byte as raw bits
template <class Sink>
ALWAYS_INLINE void put_escaped(Sink &out, int s, uint32_t v) {
    if (v < ESCAPE) return out.put(s, v);
    out.put(s, ESCAPE + (v >> 8));
    out.put_bits(v & 255, 8);
}
ALWAYS_INLINE uint32_t get_escaped(SymbolReader &in, int s) {
    const uint32_t v = in.get(s);
    return v < ESCAPE ? v : (v - ESCAPE) << 8 | uint32_t(in.get_bits(8));
}

struct FrameContext {
    int32_t *cur, *prev;  // this frame and the one before, on the lattice
    const Layout *layout;
    int64_t water_oh2, water_hh2;  // squared O-H and H-H distances of the water (lattice units)
    // the box of the water O's: corner, sides, and the bits of a point in it (x + Rx (y + Ry z)), or -1 if that
    // takes more than 56: then each coordinate in the bits of its side
    int32_t box_min[3];
    uint64_t box_size[3];
    int box_bits, side_bits[3];
    bool unbonded_box;  // unbonded atoms as their points in the box
};
inline void set_water_box(FrameContext &frame, const uint32_t *size) {  // sides of at most 2^30
    uint64_t *R = frame.box_size;
    for (int c = 0; c < 3; c++) R[c] = size[c];
    // bits of the largest value below n
    auto bits_below = [](uint64_t n) { return int(vector_bit_width(std::vector<uint64_t>{n - 1})); };
    frame.box_bits = R[0] * R[1] < (uint64_t(1) << 56) / R[2] ? bits_below(R[0] * R[1] * R[2]) : -1;
    for (int c = 0; c < 3; c++) frame.side_bits[c] = bits_below(R[c]);
}

// a water O of an intra frame: its point in the box
template <class Sink>
inline void put_box_point(const FrameContext &frame, const int32_t *p, Sink &out) {
    const uint64_t *R = frame.box_size;
    uint64_t u[3];
    for (int c = 0; c < 3; c++) u[c] = uint64_t(int64_t(p[c]) - frame.box_min[c]);
    if (frame.box_bits >= 0) return out.put_bits(u[0] + R[0] * (u[1] + R[1] * u[2]), frame.box_bits);
    for (int c = 0; c < 3; c++) out.put_bits(u[c], frame.side_bits[c]);
}
inline void get_box_point(const FrameContext &frame, int32_t *p, SymbolReader &in) {
    const uint64_t *R = frame.box_size;
    uint64_t u[3];
    if (frame.box_bits >= 0) {
        const uint64_t v = in.get_bits(frame.box_bits), yz = v / R[0];
        u[0] = v - yz * R[0];
        u[1] = yz % R[1];
        u[2] = yz / R[1];
    } else {
        for (int c = 0; c < 3; c++) u[c] = in.get_bits(frame.side_bits[c]);
    }
    for (int c = 0; c < 3; c++) p[c] = int32_t(frame.box_min[c] + int64_t(u[c]));
}

// coordinate c of unbonded atom i: from the atom before (intra) or from the previous frame
ALWAYS_INLINE int64_t predict_unbonded(const FrameContext &frame, int mode, size_t i, int c) {
    return mode ? frame.prev[3 * i + c] : (i ? frame.cur[3 * (i - 1) + c] : 0);
}

// Atom i on the sphere |d|^2 = r2 around atom j (r2 = 0: |p|^2 up to 2^28, mode 1 alone): the largest axis of d is
// dropped and rebuilt from the other two, with a radial residual. Intra (mode 0) the face (2 axis + sign) is sent with
// the residual and the two kept coordinates as they are; else the vector p of the previous frame fixes the axis and the
// sign, and the kept coordinates are residuals against p.
template <class Sink>
ALWAYS_INLINE void put_sphere(const FrameContext &frame, int mode, size_t i, size_t j, int64_t r2, int s_kept,
                              int s_radial, Sink &out) {
    int64_t d[3], p[3] = {0, 0, 0};
    int drop, k0, k1;
    displacement(frame.cur, i, j, d);
    if (mode) displacement(frame.prev, i, j, p);
    if (!r2) r2 = std::min(p[0] * p[0] + p[1] * p[1] + p[2] * p[2], int64_t(1) << 28);
    split_axes(mode ? p : d, drop, k0, k1);
    const bool neg = (mode ? p[drop] : d[drop]) < 0;
    const uint32_t radial = zigzag((neg ? -d[drop] : d[drop]) - isqrt_round(r2 - d[k0] * d[k0] - d[k1] * d[k1]));
    if (mode == 0) out.put(S_WATER_H1_FACE, uint32_t(2 * drop + neg) + 6 * std::min(radial, 15u));
    if (mode != 0 || radial >= 15) out.put(s_radial, mode == 0 ? radial - 15 : radial);
    out.put(s_kept, zigzag(d[k0] - p[k0]));
    out.put(s_kept, zigzag(d[k1] - p[k1]));
}
ALWAYS_INLINE void get_sphere(const FrameContext &frame, int mode, size_t i, size_t j, int64_t r2, int s_kept,
                              int s_radial, SymbolReader &in) {
    int64_t d[3], p[3] = {0, 0, 0};
    int drop, k0, k1;
    bool neg;
    uint32_t radial;
    if (mode == 0) {
        const uint32_t v = in.get(S_WATER_H1_FACE);
        drop = int(v % 6 / 2);
        k0 = (drop + 1) % 3;
        k1 = (drop + 2) % 3;
        neg = v & 1;
        radial = v / 6;
        if (radial == 15) radial += in.get(s_radial);
    } else {
        displacement(frame.prev, i, j, p);
        if (!r2) r2 = std::min(p[0] * p[0] + p[1] * p[1] + p[2] * p[2], int64_t(1) << 28);
        split_axes(p, drop, k0, k1);
        neg = p[drop] < 0;
        radial = in.get(s_radial);
    }
    d[k0] = p[k0] + unzigzag(in.get(s_kept));
    d[k1] = p[k1] + unzigzag(in.get(s_kept));
    // unsigned: a corrupt stream's squares may pass 2^63, a valid one's stay below 2^60
    const uint64_t rest = uint64_t(r2) - uint64_t(d[k0]) * uint64_t(d[k0]) - uint64_t(d[k1]) * uint64_t(d[k1]);
    const int64_t r = isqrt_round(int64_t(rest)) + unzigzag(radial);
    d[drop] = neg ? -r : r;
    for (int c = 0; c < 3; c++) frame.cur[3 * i + c] = int32_t(frame.cur[3 * j + c] + d[c]);
}

// H2 - O = d on the circle |d|^2 = oh2, 2 d.h1 = oh2 + |h1|^2 - hh2, h1 = H1 - O. The coordinate of d on the axis
// where h1 is smallest leaves two points of the circle, one on each side of the plane through h1 and that axis; the
// side of d (intra: sent, else as a change from the side of the prediction p) picks one, and two residuals correct it.
// Integer products and a square root and a division of exact doubles, so the decoder finds the same point. Beyond the
// lattice sizes of a water (a molecule split by the boundary) the point is {coord on the axis, 0, 0}.
ALWAYS_INLINE int circle_axis(const int64_t h1[3]) {
    int axis = std::llabs(h1[1]) < std::llabs(h1[0]) ? 1 : 0;
    return std::llabs(h1[2]) < std::llabs(h1[axis]) ? 2 : axis;
}
// the side of v: 1 when (-h1[j], h1[i]) . (v[i], v[j]) < 0, for the other two axes i, j
ALWAYS_INLINE int circle_side(const int64_t h1[3], int axis, const int64_t v[3]) {
    const int i = (axis + 1) % 3, j = (axis + 2) % 3;
    return h1[i] * v[j] - h1[j] * v[i] < 0;
}
ALWAYS_INLINE void circle_point(const int64_t h1[3], int64_t oh2, int64_t hh2, int axis, int64_t coord, int side,
                                int64_t point[3]) {
    const int i = (axis + 1) % 3, j = (axis + 2) % 3;
    const int64_t u = h1[i], w = h1[j], lim = 16384;
    const bool ok = std::llabs(u) < lim && std::llabs(w) < lim && std::llabs(h1[axis]) < lim &&
                    std::llabs(coord) < lim && u * u + w * w > 0;
    const int64_t g2 = ok ? u * u + w * w : 1;
    const int64_t M = ok ? oh2 + g2 + h1[axis] * h1[axis] - hh2 - 2 * h1[axis] * coord : 0;  // 2 (u P_i + w P_j)
    const int64_t s = ok ? (side ? -1 : 1) * isqrt_round(4 * g2 * (oh2 - coord * coord) - M * M) : 0;
    point[axis] = coord;
    point[i] = ok ? round_half_away(double(M * u - s * w) / double(2 * g2)) : 0;
    point[j] = ok ? round_half_away(double(M * w + s * u) / double(2 * g2)) : 0;
}
// H2 - O of the previous frame, or 0 beyond the lattice sizes of a water, so that 2 zigzag(d - p) + 1 fits 32 bits
ALWAYS_INLINE void previous_h2(const FrameContext &frame, size_t o, int64_t p[3]) {
    displacement(frame.prev, o + 2, o, p);
    if (std::llabs(p[0]) >= 16384 || std::llabs(p[1]) >= 16384 || std::llabs(p[2]) >= 16384) p[0] = p[1] = p[2] = 0;
}
template <class Sink>
ALWAYS_INLINE void put_circle(const FrameContext &frame, int mode, size_t o, Sink &out) {
    int64_t h1[3], d[3], p[3], point[3];
    displacement(frame.cur, o + 1, o, h1);
    displacement(frame.cur, o + 2, o, d);
    const int axis = circle_axis(h1), side = circle_side(h1, axis, d);
    if (mode) {
        previous_h2(frame, o, p);
        out.put(S_WATER_H2_AXIS, zigzag(d[axis] - p[axis]) * 2 + uint32_t(side != circle_side(h1, axis, p)));
    } else {
        out.put(S_WATER_H2_AXIS, zigzag(d[axis]) * 2 + uint32_t(side));
    }
    circle_point(h1, frame.water_oh2, frame.water_hh2, axis, d[axis], side, point);
    for (int c = 0; c < 3; c++)
        if (c != axis) out.put(S_WATER_H2_RESIDUAL, zigzag(d[c] - point[c]));
}
ALWAYS_INLINE void get_circle(const FrameContext &frame, int mode, size_t o, SymbolReader &in) {
    int64_t h1[3], p[3] = {0, 0, 0}, point[3];
    displacement(frame.cur, o + 1, o, h1);
    const int axis = circle_axis(h1);
    if (mode) previous_h2(frame, o, p);
    const uint32_t v = in.get(S_WATER_H2_AXIS);
    const int side = int(v & 1) ^ (mode ? circle_side(h1, axis, p) : 0);
    circle_point(h1, frame.water_oh2, frame.water_hh2, axis, p[axis] + unzigzag(v >> 1), side, point);
    for (int c = 0; c < 3; c++)
        frame.cur[3 * (o + 2) + c] =
            int32_t(frame.cur[3 * o + c] + (c == axis ? point[c] : point[c] + unzigzag(in.get(S_WATER_H2_RESIDUAL))));
}

// A water O: from the previous frame, or (intra) its point in the box.
template <class Sink>
ALWAYS_INLINE void put_water_o(const FrameContext &frame, int mode, size_t i, Sink &out) {
    if (mode)
        for (int c = 0; c < 3; c++) out.put(S_WATER_O, zigzag(frame.cur[3 * i + c] - frame.prev[3 * i + c]));
    else
        put_box_point(frame, frame.cur + 3 * i, out);
}

// One atom: a water (O and the rest of its molecule), or another atom. Water H come with their O.
template <class Sink>
ALWAYS_INLINE void put_atom(const FrameContext &frame, const int *mode, size_t i, Sink &out) {
    const Layout &layout = *frame.layout;
    const int32_t *q = frame.cur;
    if (layout.kind[i] == K_WATER_O) {
        put_water_o(frame, mode[G_WATER_O], i, out);
        put_sphere(frame, mode[G_WATER_H], i + 1, i, frame.water_oh2, S_WATER_H1_KEPT, S_WATER_H1_RADIAL, out);
        put_circle(frame, mode[G_WATER_H], i, out);
    } else if (layout.kind[i] == K_OTHER) {
        const size_t j = i - layout.bond[i];
        if (j != i && mode[G_BONDED])
            put_sphere(frame, 1, i, j, 0, S_BOND_KEPT, S_BOND_RADIAL, out);
        else if (j != i)
            for (int c = 0; c < 3; c++) put_escaped(out, S_BOND_DELTA, zigzag(q[3 * i + c] - q[3 * j + c]));
        else if (frame.unbonded_box)
            put_box_point(frame, q + 3 * i, out);
        else
            for (int c = 0; c < 3; c++)
                put_escaped(out, S_UNBONDED, zigzag(q[3 * i + c] - predict_unbonded(frame, mode[G_UNBONDED], i, c)));
    }
}
ALWAYS_INLINE void get_atom(const FrameContext &frame, const int *mode, size_t i, SymbolReader &in) {
    const Layout &layout = *frame.layout;
    int32_t *q = frame.cur;
    if (layout.kind[i] == K_WATER_O) {
        if (mode[G_WATER_O])
            for (int c = 0; c < 3; c++) q[3 * i + c] = int32_t(frame.prev[3 * i + c] + unzigzag(in.get(S_WATER_O)));
        else
            get_box_point(frame, q + 3 * i, in);
        get_sphere(frame, mode[G_WATER_H], i + 1, i, frame.water_oh2, S_WATER_H1_KEPT, S_WATER_H1_RADIAL, in);
        get_circle(frame, mode[G_WATER_H], i, in);
    } else if (layout.kind[i] == K_OTHER) {
        const size_t j = i - layout.bond[i];
        if (j != i && mode[G_BONDED])
            get_sphere(frame, 1, i, j, 0, S_BOND_KEPT, S_BOND_RADIAL, in);
        else if (j != i)
            for (int c = 0; c < 3; c++) q[3 * i + c] = int32_t(q[3 * j + c] + unzigzag(get_escaped(in, S_BOND_DELTA)));
        else if (frame.unbonded_box)
            get_box_point(frame, q + 3 * i, in);
        else
            for (int c = 0; c < 3; c++)
                q[3 * i + c] =
                    int32_t(predict_unbonded(frame, mode[G_UNBONDED], i, c) + unzigzag(get_escaped(in, S_UNBONDED)));
    }
}

// q = round(x * inv_step), inv_step the reciprocal of the step
template <class T>
inline void quantize(const T *x, size_t n, double inv_step, int32_t *q) {
    for (size_t i = 0; i < n; i++) {
        const double y = double(x[i]) * inv_step;
        q[i] = int32_t(y + std::copysign(0.5, y));
    }
}

// the box of the water O's of frame q
inline void water_box(const int32_t *q, const std::vector<uint32_t> &waters, int32_t lo[3], int32_t hi[3]) {
    for (int c = 0; c < 3; c++) lo[c] = hi[c] = waters.empty() ? 0 : q[3 * waters[0] + c];
    for (uint32_t i : waters)
        for (int c = 0; c < 3; c++) {
            lo[c] = std::min(lo[c], q[3 * i + c]);
            hi[c] = std::max(hi[c], q[3 * i + c]);
        }
}

// Per group but the water O's, the mode (0 intra, 1 previous frame) that costs fewer bits on a sample of about 256 of
// its atoms.
inline void choose_modes(const FrameContext &frame, const std::vector<uint32_t> *const group_atoms[NUM_GROUPS],
                         int mode[NUM_GROUPS]) {
    for (int g = G_WATER_H; g < NUM_GROUPS; g++) {
        const auto &sample = *group_atoms[g];
        double best = 1e300;
        for (int m = 0; m < 2 && !sample.empty(); m++) {
            CostEstimate cost;
            int trial[NUM_GROUPS] = {mode[0], mode[1], mode[2], mode[3]};
            trial[g] = m;
            const size_t stride = std::max<size_t>(1, sample.size() / 256);
            for (size_t k = 0; k < sample.size(); k += stride) put_atom(frame, trial, sample[k], cost);
            if (cost.bits < best) best = cost.bits, mode[g] = m;
        }
    }
}

}  // namespace biomd

template <class T, uint N>
class SZBioMDDecomposition : public concepts::DecompositionInterface<T, int, N> {
    static_assert(std::is_floating_point<T>::value, "SZ3 BioMD: data must be float or double");

   public:
    explicit SZBioMDDecomposition(const Config &conf)
        : frames_(N == 3 ? conf.dims[0] : 1), atoms_(N >= 2 ? conf.dims[N - 2] : 1), error_bound_(conf.absErrorBound) {
        if (N > 3 || conf.dims[N - 1] != 3) throw std::invalid_argument("SZ3 BioMD: data must be {frames, atoms, 3}");
    }

    // The symbol streams of the chunk, one after the other as [count, symbols]; the rest goes to save().
    std::vector<int> compress(const Config & /*conf*/, T *data) override {
        using namespace biomd;
        const size_t frame_values = atoms_ * 3;
        find_fill_frames(data);
        choose_step(data);
        Layout layout;
        std::vector<uint32_t> waters, bonded, unbonded;  // water O's, other atoms
        detect_layout(data, layout, waters, bonded, unbonded);
        // one frame without water and mostly unbonded (coarse-grained): box points cost about the bits of deltas, and
        // less time
        unbonded_box_ = coded_frames_ == 1 && waters.empty() && bonded.size() * 4 < atoms_;

        // a buffer per stream, as large as the symbols the chunk can put in it (every value is written before it is
        // read)
        const std::vector<uint32_t> *group_atoms[NUM_GROUPS] = {&waters, &waters, &bonded, &unbonded};
        std::unique_ptr<int[]> buffers[NUM_STREAMS];
        SymbolWriter writer;
        for (int s = 0; s < NUM_STREAMS; s++) {
            const size_t n = STREAM_GROUP[s] < 0 ? atoms_ : coded_frames_ * group_atoms[STREAM_GROUP[s]]->size();
            buffers[s].reset(new int[MAX_SYMBOLS_PER_UNIT[s] * n + 1]);
            writer.cursor[s] = buffers[s].get();
        }
        raw_bits_.clear();
        BitAppender raw(raw_bits_);
        writer.raw = &raw;
        for (size_t k = 0; k < waters.size(); k++) writer.put(S_WATER_GAP, waters[k] - (k ? waters[k - 1] : 0));
        for (size_t i = 0; i < atoms_; i++)
            if (layout.kind[i] == K_OTHER) writer.put(S_BOND_REF, layout.bond[i]);

        std::unique_ptr<int32_t[]> lattice(new int32_t[2 * frame_values]);  // written before read
        FrameContext frame = frame_context(lattice.get(), layout);
        const int intra[NUM_GROUPS] = {0, 0, 0, 0};
        std::fill(modes_ + G_WATER_H, modes_ + NUM_GROUPS, 0);
        for (size_t t = 0; t < coded_frames_; t++) {
            quantize(data + t * frame_values, frame_values, 1.0 / step_, frame.cur);
            if (t == 0) {  // sides of 1 without water O's or unbonded box points
                int32_t box_max[3];
                water_box(frame.cur, unbonded_box_ ? unbonded : waters, frame.box_min, box_max);
                for (int c = 0; c < 3; c++) box_min_[c] = frame.box_min[c], box_size_[c] = box_max[c] - box_min_[c] + 1;
                set_water_box(frame, box_size_);
            }
            // per group, the predictor that is cheapest on a sample of frame 1, for every frame after the first
            if (t == 1) choose_modes(frame, group_atoms, modes_);
            for (size_t i = 0; i < atoms_; i++) put_atom(frame, t ? modes_ : intra, i, writer);
            std::swap(frame.cur, frame.prev);
        }
        raw.flush();

        std::vector<int> bins;
        size_t n = NUM_STREAMS;
        for (int s = 0; s < NUM_STREAMS; s++) {
            if (size_t(writer.cursor[s] - buffers[s].get()) > size_t(std::numeric_limits<int>::max()))
                throw std::invalid_argument("SZ3 BioMD: a stream of more symbols than an int counts");
            n += size_t(writer.cursor[s] - buffers[s].get());
        }
        bins.reserve(n);
        for (int s = 0; s < NUM_STREAMS; s++) {
            bins.push_back(int(writer.cursor[s] - buffers[s].get()));
            bins.insert(bins.end(), buffers[s].get(), writer.cursor[s]);
        }
        return bins;
    }

    // The chunk from its symbol streams, each as [count, symbols], after load().
    T *decompress(const Config & /*conf*/, std::vector<int> &quant_inds, T *dec_data) override {
        using namespace biomd;
        const size_t frame_values = atoms_ * 3;
        SymbolReader reader;
        for (size_t k = 0, s = 0; s < NUM_STREAMS; s++) {
            const size_t n = k < quant_inds.size() ? size_t(quant_inds[k]) : 0;
            if (k >= quant_inds.size() || n > quant_inds.size() - k - 1)
                throw std::runtime_error("SZ3 BioMD: corrupt stream");
            reader.cursor[s] = quant_inds.data() + k + 1;
            reader.end[s] = reader.cursor[s] + n;
            k += n + 1;
        }
        reader.raw = BitConsumer(raw_bits_.data(), raw_bits_.data() + raw_bits_.size());
        Layout layout;
        read_layout(reader, layout);
        std::unique_ptr<int32_t[]> lattice(new int32_t[2 * frame_values]);  // written before read
        FrameContext frame = frame_context(lattice.get(), layout);
        const int intra[NUM_GROUPS] = {0, 0, 0, 0};
        std::copy(box_min_, box_min_ + 3, frame.box_min);
        set_water_box(frame, box_size_);
        for (size_t t = 0; t < coded_frames_; t++) {
            for (size_t i = 0; i < atoms_; i += layout.kind[i] == K_WATER_O ? 3 : 1)  // a water's H come with its O
                get_atom(frame, t ? modes_ : intra, i, reader);
            T *x = dec_data + t * frame_values;
            for (size_t i = 0; i < frame_values; i++) x[i] = T(double(frame.cur[i]) * step_);
            std::swap(frame.cur, frame.prev);
        }
        // a stream that holds more than the frames took is corrupt
        for (int s = 0; s < NUM_STREAMS; s++)
            if (reader.cursor[s] != reader.end[s]) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        if (!reader.raw.at_end()) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        std::fill(dec_data + coded_frames_ * frame_values, dec_data + frames_ * frame_values, fill_value_);
        return dec_data;
    }

    void save(uchar *&c) override {
        write(uint64_t(coded_frames_), c);
        write(fill_value_, c);
        write(step_, c);
        write(water_oh2_, c);
        write(water_hh2_, c);
        write(uint64_t(num_waters_), c);
        for (int m : modes_) write(uint8_t(m), c);
        write(uint8_t(unbonded_box_), c);
        if (has_box()) {
            write(box_min_, 3, c);
            write(box_size_, 3, c);
        }
        write(uint64_t(raw_bits_.size()), c);
        if (!raw_bits_.empty()) write(raw_bits_.data(), raw_bits_.size(), c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uint64_t waters, coded, raw_size = 0;
        read(coded, c, remaining_length);
        if (coded > frames_) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        coded_frames_ = size_t(coded);
        read(fill_value_, c, remaining_length);
        read(step_, c, remaining_length);
        read(water_oh2_, c, remaining_length);
        read(water_hh2_, c, remaining_length);
        read(waters, c, remaining_length);
        num_waters_ = size_t(std::min<uint64_t>(waters, atoms_));
        for (int &m : modes_) {
            uint8_t v;
            read(v, c, remaining_length);
            if (v > 1) throw std::runtime_error("SZ3 BioMD: corrupt stream");
            m = v;
        }
        uint8_t unbonded_box;
        read(unbonded_box, c, remaining_length);
        if (unbonded_box > 1) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        unbonded_box_ = unbonded_box && coded_frames_ == 1;
        if (has_box()) {
            read(box_min_, 3, c, remaining_length);
            read(box_size_, 3, c, remaining_length);
        }
        read(raw_size, c, remaining_length);
        if (raw_size > remaining_length) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        raw_bits_.resize(size_t(raw_size));
        if (raw_size) read(raw_bits_.data(), raw_bits_.size(), c, remaining_length);
        bool ok = std::isfinite(step_) && step_ > 0 && waters <= atoms_ && water_oh2_ >= 0 &&
                  water_oh2_ < (int64_t(1) << 28) && water_hh2_ >= 0 && water_hh2_ < (int64_t(1) << 30);
        for (uint32_t side : box_size_) ok = ok && side >= 1 && side <= (uint32_t(1) << 30);
        if (!ok) throw std::runtime_error("SZ3 BioMD: corrupt stream");
    }

    // a bound on what save() writes
    size_t size_est() override { return 96 + raw_bits_.size(); }

    // every bin, a count or a symbol, is in [0, INT_MAX]
    std::pair<int, int> get_out_range() override { return {0, std::numeric_limits<int>::max()}; }

   private:
    // frame 0 stores the box of its water O's, or of its unbonded atoms
    bool has_box() const { return unbonded_box_ || num_waters_; }

    // Trailing frames all of one value (the unwritten rest of a chunk) are stored as that value; a chunk may be all
    // fill (an unwritten chunk of NaN).
    void find_fill_frames(const T *data) {
        const size_t frame_values = atoms_ * 3;
        fill_value_ = data[(frames_ - 1) * frame_values];
        auto is_fill = [&](size_t t) {
            for (size_t i = t * frame_values; i < (t + 1) * frame_values; i++)
                if (memcmp(&data[i], &fill_value_, sizeof(T)) != 0) return false;
            return true;
        };
        for (coded_frames_ = frames_; coded_frames_ > 0 && is_fill(coded_frames_ - 1);) coded_frames_--;
    }

    // Lattice step: 2 eb - m ulp for the binade below B, so that q step rounded to T stays within eb for every
    // |q step| < B. Rounding q step to float adds ulp / 2 to the step / 2 of the lattice, and the double arithmetic
    // (x * inv_step rounded to q is off by up to 2 double ulps of x, q step by 1 more) under 2^-20 float ulps: m = 1 +
    // 2^-10 for float. For double those ulps are its own: m = 8. B, the power of two above
    // 2^LATTICE_SPAN_BITS eb, depends on the bound alone, so a decompressed chunk compressed again lands on the same
    // lattice points. A chunk with NaN or Inf goes to lossless compression; one with coordinates beyond B, or under a
    // bound for which m ulp reaches eb or B passes the largest value of T, also to lossless compression.
    void choose_step(const T *data) {
        if (!(error_bound_ > 0) || !std::isfinite(error_bound_))
            throw std::invalid_argument("SZ3 BioMD: the error bound must be positive and finite");
        using Bits = typename std::conditional<sizeof(T) == 4, uint32_t, uint64_t>::type;
        Bits max_bits = 0;  // the largest |x| as bits: NaN and inf are above every finite value
        for (size_t i = 0; i < coded_frames_ * atoms_ * 3; i++) {
            Bits b;
            memcpy(&b, &data[i], sizeof(T));
            b &= Bits(~Bits(0)) >> 1;
            max_bits = b > max_bits ? b : max_bits;
        }
        // compared as bits: under -ffinite-math-only a floating-point test for NaN may be compiled away
        const Bits infinity_bits = Bits(sizeof(T) == 4 ? 0x7f800000ull : 0x7ff0000000000000ull);
        if (max_bits >= infinity_bits)
            throw biomd::Fallback("SZ3 BioMD: NaN or Inf in a frame that is not trailing fill");
        T max_abs;
        memcpy(&max_abs, &max_bits, sizeof(T));
        int exponent;
        std::frexp(std::ldexp(error_bound_, biomd::LATTICE_SPAN_BITS<T>), &exponent);
        const double span = std::ldexp(1.0, exponent);
        const double margin = (sizeof(T) == 4 ? 1 + 1.0 / 1024 : 8) *
                              std::max(double(std::numeric_limits<T>::denorm_min()),
                                       std::ldexp(1.0, exponent - std::numeric_limits<T>::digits));
        step_ = 2.0 * error_bound_ - margin;
        // every lattice point |q| step below B, with q rounded as quantize() does (after max|x| < B + step, which
        // keeps q within MAX_LATTICE): decompressed values land on the same points, and compress the same again
        if (!(margin < error_bound_) || !std::isfinite(1.0 / step_) ||
            !(span <= double(std::numeric_limits<T>::max())) || !(span / step_ < biomd::MAX_LATTICE) ||
            !(double(max_abs) < span + step_) ||
            !(nofma(std::fabs(double(biomd::round_half_away(nofma(double(max_abs) * (1.0 / step_))))) * step_) < span))
            throw biomd::Fallback("SZ3 BioMD: coordinates beyond the lattice the error bound allows");
    }

    // The layout of the first frame, its geometry on the lattice, and the atoms of each group.
    void detect_layout(const T *data, biomd::Layout &layout, std::vector<uint32_t> &waters,
                       std::vector<uint32_t> &bonded, std::vector<uint32_t> &unbonded) {
        using namespace biomd;
        detect_water(data, atoms_, layout);
        const double oh = layout.water_oh / step_, hh = layout.water_hh / step_;
        water_oh2_ = oh < 16384 ? round_half_away(nofma(oh * oh)) : 0;  // squares that load() accepts
        water_hh2_ = hh < 32768 ? round_half_away(nofma(hh * hh)) : 0;
        if (!(water_oh2_ > 0 && water_oh2_ < (int64_t(1) << 28) && water_hh2_ > 0 && water_hh2_ < (int64_t(1) << 30))) {
            water_oh2_ = water_hh2_ = 0;  // no water, or too many lattice steps across one for exact products
            std::fill(layout.kind.begin(), layout.kind.end(), K_OTHER);
        }
        detect_bonds(data, atoms_, layout);
        for (size_t i = 0; i < atoms_; i++) {
            if (layout.kind[i] == K_WATER_O) waters.push_back(uint32_t(i));
            if (layout.kind[i] == K_OTHER) (layout.bond[i] ? bonded : unbonded).push_back(uint32_t(i));
        }
        num_waters_ = waters.size();
    }

    // The layout from its streams: the water O's by their gaps, then the bond of each other atom.
    void read_layout(biomd::SymbolReader &reader, biomd::Layout &layout) const {
        using namespace biomd;
        layout.kind.assign(atoms_, K_OTHER);
        for (size_t k = 0, o = 0; k < num_waters_; k++) {
            const uint32_t gap = reader.get(S_WATER_GAP);
            if ((k && gap < 3) || gap > atoms_ || o + gap + 2 >= atoms_)
                throw std::runtime_error("SZ3 BioMD: corrupt stream");
            o += gap;
            layout.kind[o] = K_WATER_O;
            layout.kind[o + 1] = layout.kind[o + 2] = K_WATER_H;
        }
        layout.bond.assign(atoms_, 0);
        for (size_t i = 0; i < atoms_; i++)
            if (layout.kind[i] == K_OTHER) {
                const uint32_t bond = reader.get(S_BOND_REF);
                if (bond > uint32_t(MAX_BOND_OFFSET) || bond > i) throw std::runtime_error("SZ3 BioMD: corrupt stream");
                layout.bond[i] = uint8_t(bond);
            }
    }

    // the lattice holds this frame and the one before
    biomd::FrameContext frame_context(int32_t *lattice, const biomd::Layout &layout) const {
        return {lattice,   lattice + atoms_ * 3, &layout, water_oh2_, water_hh2_, {0, 0, 0}, {1, 1, 1}, 0,
                {0, 0, 0}, unbonded_box_};
    }

    size_t frames_, atoms_, coded_frames_ = 1;  // coded: the frames before the fill
    T fill_value_ = 0;
    double error_bound_, step_ = 0;
    int64_t water_oh2_ = 0, water_hh2_ = 0;  // squared O-H and H-H distances of the water (lattice units)
    size_t num_waters_ = 0;
    int modes_[biomd::NUM_GROUPS] = {1, 0, 0, 0};  // per group: the predictor of the frames after the first
    bool unbonded_box_ = false;                    // one frame, its unbonded atoms as points in their box
    int32_t box_min_[3] = {0, 0, 0};               // the box of frame 0: corner
    uint32_t box_size_[3] = {1, 1, 1};             // and sides
    std::vector<uchar> raw_bits_;                  // intra water O's and the low bytes of escaped values
};

template <class T, uint N>
SZBioMDDecomposition<T, N> make_decomposition_biomd(const Config &conf) {
    return SZBioMDDecomposition<T, N>(conf);
}

}  // namespace SZ3
#endif
