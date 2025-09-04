// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARTYPES_H
#define COLVARTYPES_H

#include <sstream> // TODO specialize templates and replace this with iosfwd
#include <vector>

#ifdef COLVARS_LAMMPS
// Use open-source Jacobi implementation
#include "math_eigen_impl.h"
#endif

#include "colvarmodule.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

// ----------------------------------------------------------------------
/// Linear algebra functions and data types used in the collective
/// variables implemented so far
// ----------------------------------------------------------------------


/// \brief Arbitrary size array (one dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T> class colvarmodule::vector1d
{
protected:

  std::vector<T> data;

public:

  /// Default constructor
  inline vector1d(size_t const n = 0)
  {
    data.resize(n);
    reset();
  }

  /// Constructor from C array
  inline vector1d(size_t const n, T const *t)
  {
    data.resize(n);
    reset();
    size_t i;
    for (i = 0; i < size(); i++) {
      data[i] = t[i];
    }
  }

  /// Explicit Copy constructor
  inline vector1d(const vector1d&) = default;

  /// Explicit Copy assignement
  inline vector1d& operator=(const vector1d&) = default;

  /// Return a pointer to the data location
  inline T * c_array()
  {
    if (data.size() > 0) {
      return &(data[0]);
    } else {
      return NULL;
    }
  }

  /// Return a reference to the data
  inline std::vector<T> &data_array()
  {
    return data;
  }

  /// Return a reference to the data
  inline std::vector<T> const &data_array() const
  {
    return data;
  }

  inline ~vector1d()
  {
    data.clear();
  }

  /// Set all elements to zero
  inline void reset()
  {
    data.assign(data.size(), T(0.0));
  }

  inline size_t size() const
  {
    return data.size();
  }

  inline void resize(size_t const n)
  {
    data.resize(n);
  }

  inline void clear()
  {
    data.clear();
  }

  inline T & operator [] (size_t const i) {
    return data[i];
  }

  inline T const & operator [] (size_t const i) const {
    return data[i];
  }

  inline static void check_sizes(vector1d<T> const &v1, vector1d<T> const &v2)
  {
    if (v1.size() != v2.size()) {
      cvm::error("Error: trying to perform an operation between vectors of different sizes, "+
                 cvm::to_str(v1.size())+" and "+cvm::to_str(v2.size())+".\n");
    }
  }

  inline void operator += (vector1d<T> const &v)
  {
    check_sizes(*this, v);
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] += v[i];
    }
  }

  inline void operator -= (vector1d<T> const &v)
  {
    check_sizes(*this, v);
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] -= v[i];
    }
  }

  inline void operator *= (cvm::real a)
  {
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] *= a;
    }
  }

  inline void operator /= (cvm::real a)
  {
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] /= a;
    }
  }

  inline friend vector1d<T> operator + (vector1d<T> const &v1,
                                        vector1d<T> const &v2)
  {
    check_sizes(v1.size(), v2.size());
    vector1d<T> result(v1.size());
    size_t i;
    for (i = 0; i < v1.size(); i++) {
      result[i] = v1[i] + v2[i];
    }
    return result;
  }

  inline friend vector1d<T> operator - (vector1d<T> const &v1,
                                        vector1d<T> const &v2)
  {
    check_sizes(v1.size(), v2.size());
    vector1d<T> result(v1.size());
    size_t i;
    for (i = 0; i < v1.size(); i++) {
      result[i] = v1[i] - v2[i];
    }
    return result;
  }

  inline friend vector1d<T> operator * (vector1d<T> const &v, cvm::real a)
  {
    vector1d<T> result(v.size());
    size_t i;
    for (i = 0; i < v.size(); i++) {
      result[i] = v[i] * a;
    }
    return result;
  }

  inline friend vector1d<T> operator * (cvm::real a, vector1d<T> const &v)
  {
    return v * a;
  }

  inline friend vector1d<T> operator / (vector1d<T> const &v, cvm::real a)
  {
    vector1d<T> result(v.size());
    size_t i;
    for (i = 0; i < v.size(); i++) {
      result[i] = v[i] / a;
    }
    return result;
  }

  /// Inner product
  inline friend T operator * (vector1d<T> const &v1, vector1d<T> const &v2)
  {
    check_sizes(v1.size(), v2.size());
    T prod(0.0);
    size_t i;
    for (i = 0; i < v1.size(); i++) {
      prod += v1[i] * v2[i];
    }
    return prod;
  }

  /// Squared norm
  inline cvm::real norm2() const
  {
    cvm::real result = 0.0;
    size_t i;
    for (i = 0; i < this->size(); i++) {
      result += (*this)[i] * (*this)[i];
    }
    return result;
  }

  inline cvm::real norm() const
  {
    return cvm::sqrt(this->norm2());
  }

  inline cvm::real sum() const
  {
    cvm::real result = 0.0;
    size_t i;
    for (i = 0; i < this->size(); i++) {
      result += (*this)[i];
    }
    return result;
  }

  /// Slicing
  inline vector1d<T> const slice(size_t const i1, size_t const i2) const
  {
    if ((i2 < i1) || (i2 >= this->size())) {
      cvm::error("Error: trying to slice a vector using incorrect boundaries.\n");
    }
    vector1d<T> result(i2 - i1);
    size_t i;
    for (i = 0; i < (i2 - i1); i++) {
      result[i] = (*this)[i1+i];
    }
    return result;
  }

  /// Assign a vector to a slice of this vector
  inline void sliceassign(size_t const i1, size_t const i2,
                          vector1d<T> const &v)
  {
    if ((i2 < i1) || (i2 >= this->size())) {
      cvm::error("Error: trying to slice a vector using incorrect boundaries.\n");
    }
    size_t i;
    for (i = 0; i < (i2 - i1); i++) {
      (*this)[i1+i] = v[i];
    }
  }

  /// Formatted output

  inline size_t output_width(size_t real_width) const
  {
    return real_width*(this->size()) + 3*(this->size()-1) + 4;
  }

  inline friend std::istream & operator >> (std::istream &is,
                                            cvm::vector1d<T> &v)
  {
    if (v.size() == 0) return is;
    std::streampos const start_pos = is.tellg();
    char sep;
    if ( !(is >> sep) || !(sep == '(') ) {
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    }
    size_t count = 0;
    while ( (is >> v[count]) &&
            (count < (v.size()-1) ? ((is >> sep) && (sep == ',')) : true) ) {
      if (++count == v.size()) break;
    }
    if (count < v.size()) {
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
    }
    return is;
  }

  inline friend std::ostream & operator << (std::ostream &os,
                                            cvm::vector1d<T> const &v)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    os.width(2);
    os << "( ";
    size_t i;
    for (i = 0; i < v.size()-1; i++) {
      os.width(w); os.precision(p);
      os << v[i] << " , ";
    }
    os.width(w); os.precision(p);
    os << v[v.size()-1] << " )";
    return os;
  }

  inline std::string to_simple_string() const
  {
    if (this->size() == 0) return std::string("");
    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(cvm::cv_prec);
    os << (*this)[0];
    size_t i;
    for (i = 1; i < this->size(); i++) {
      os << " " << (*this)[i];
    }
    return os.str();
  }

  inline int from_simple_string(std::string const &s)
  {
    std::stringstream stream(s);
    size_t i = 0;
    if (this->size()) {
      while ((stream >> (*this)[i]) && (i < this->size())) {
        i++;
      }
      if (i < this->size()) {
        return COLVARS_ERROR;
      }
    } else {
      T input;
      while (stream >> input) {
        if ((i % 100) == 0) {
          data.reserve(data.size()+100);
        }
        data.resize(data.size()+1);
        data[i] = input;
        i++;
      }
    }
    return COLVARS_OK;
  }

};


/// \brief Arbitrary size array (two dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T> class colvarmodule::matrix2d
{
public:

  friend class row;
  size_t outer_length;
  size_t inner_length;

protected:

  class row {
  public:
    T * data;
    size_t length;
    inline row(T * const row_data, size_t const inner_length)
      : data(row_data), length(inner_length)
    {}
    inline T & operator [] (size_t const j) {
      return *(data+j);
    }
    inline T const & operator [] (size_t const j) const {
      return *(data+j);
    }
    inline operator vector1d<T>() const
    {
      return vector1d<T>(length, data);
    }
    inline int set(cvm::vector1d<T> const &v) const
    {
      if (v.size() != length) {
        return cvm::error("Error: setting a matrix row from a vector of "
                          "incompatible size.\n", COLVARS_BUG_ERROR);
      }
      for (size_t i = 0; i < length; i++) data[i] = v[i];
      return COLVARS_OK;
    }
  };

  std::vector<T> data;
  std::vector<row> rows;
  std::vector<T *> pointers;

public:

  /// Allocation routine, used by all constructors
  inline void resize(size_t const ol, size_t const il)
  {
    if ((ol > 0) && (il > 0)) {

      if (data.size() > 0) {
        // copy previous data
        size_t i, j;
        std::vector<T> new_data(ol * il);
        for (i = 0; i < outer_length; i++) {
          for (j = 0; j < inner_length; j++) {
            new_data[il*i+j] = data[inner_length*i+j];
          }
        }
        data.resize(ol * il);
        // copy them back
        data = new_data;
      } else {
        data.resize(ol * il);
      }

      outer_length = ol;
      inner_length = il;

      if (data.size() > 0) {
        // rebuild rows
        size_t i;
        rows.clear();
        rows.reserve(outer_length);
        pointers.clear();
        pointers.reserve(outer_length);
        for (i = 0; i < outer_length; i++) {
          rows.push_back(row(&(data[0])+inner_length*i, inner_length));
          pointers.push_back(&(data[0])+inner_length*i);
        }
     }
    } else {
      // zero size
      data.clear();
      rows.clear();
    }
  }

  /// Deallocation routine
  inline void clear() {
    rows.clear();
    data.clear();
  }

  /// Set all elements to zero
  inline void reset()
  {
    data.assign(data.size(), T(0.0));
  }

  inline size_t size() const
  {
    return data.size();
  }

  /// Default constructor
  inline matrix2d()
    : outer_length(0), inner_length(0)
  {
    this->resize(0, 0);
  }

  inline matrix2d(size_t const ol, size_t const il)
    : outer_length(ol), inner_length(il)
  {
    this->resize(outer_length, inner_length);
    reset();
  }

  /// Copy constructor
  inline matrix2d(matrix2d<T> const &m)
    : outer_length(m.outer_length), inner_length(m.inner_length)
  {
    // reinitialize data and rows arrays
    this->resize(outer_length, inner_length);
    // copy data
    data = m.data;
  }

  /// Destructor
  inline ~matrix2d() {
    this->clear();
  }

  /// Return a reference to the data
  inline std::vector<T> &data_array()
  {
    return data;
  }

  /// Return a reference to the data
  inline std::vector<T> const &data_array() const
  {
    return data;
  }

  inline row & operator [] (size_t const i)
  {
    return rows[i];
  }
  inline row const & operator [] (size_t const i) const
  {
    return rows[i];
  }

  /// Assignment
  inline matrix2d<T> & operator = (matrix2d<T> const &m)
  {
    if ((outer_length != m.outer_length) || (inner_length != m.inner_length)){
      this->clear();
      outer_length = m.outer_length;
      inner_length = m.inner_length;
      this->resize(outer_length, inner_length);
    }
    data = m.data;
    return *this;
  }

  /// Return the 2-d C array
  inline T ** c_array() {
    if (rows.size() > 0) {
      return &(pointers[0]);
    } else {
      return NULL;
    }
  }

  inline static void check_sizes(matrix2d<T> const &m1, matrix2d<T> const &m2)
  {
    if ((m1.outer_length != m2.outer_length) ||
        (m1.inner_length != m2.inner_length)) {
      cvm::error("Error: trying to perform an operation between "
                 "matrices of different sizes, "+
                 cvm::to_str(m1.outer_length)+"x"+
                 cvm::to_str(m1.inner_length)+" and "+
                 cvm::to_str(m2.outer_length)+"x"+
                 cvm::to_str(m2.inner_length)+".\n");
    }
  }

  inline void operator += (matrix2d<T> const &m)
  {
    check_sizes(*this, m);
    size_t i;
    for (i = 0; i < data.size(); i++) {
      data[i] += m.data[i];
    }
  }

  inline void operator -= (matrix2d<T> const &m)
  {
    check_sizes(*this, m);
    size_t i;
    for (i = 0; i < data.size(); i++) {
      data[i] -= m.data[i];
    }
  }

  inline void operator *= (cvm::real a)
  {
    size_t i;
    for (i = 0; i < data.size(); i++) {
      data[i] *= a;
    }
  }

  inline void operator /= (cvm::real a)
  {
    size_t i;
    for (i = 0; i < data.size(); i++) {
      data[i] /= a;
    }
  }

  inline friend matrix2d<T> operator + (matrix2d<T> const &m1,
                                        matrix2d<T> const &m2)
  {
    check_sizes(m1, m2);
    matrix2d<T> result(m1.outer_length, m1.inner_length);
    size_t i;
    for (i = 0; i < m1.data.size(); i++) {
      result.data[i] = m1.data[i] + m2.data[i];
    }
    return result;
  }

  inline friend matrix2d<T> operator - (matrix2d<T> const &m1,
                                        matrix2d<T> const &m2)
  {
    check_sizes(m1, m2);
    matrix2d<T> result(m1.outer_length, m1.inner_length);
    size_t i;
    for (i = 0; i < m1.data.size(); i++) {
      result.data[i] = m1.data[i] - m1.data[i];
    }
    return result;
  }

  inline friend matrix2d<T> operator * (matrix2d<T> const &m, cvm::real a)
  {
    matrix2d<T> result(m.outer_length, m.inner_length);
    size_t i;
    for (i = 0; i < m.data.size(); i++) {
      result.data[i] = m.data[i] * a;
    }
    return result;
  }

  inline friend matrix2d<T> operator * (cvm::real a, matrix2d<T> const &m)
  {
    return m * a;
  }

  inline friend matrix2d<T> operator / (matrix2d<T> const &m, cvm::real a)
  {
    matrix2d<T> result(m.outer_length, m.inner_length);
    size_t i;
    for (i = 0; i < m.data.size(); i++) {
      result.data[i] = m.data[i] * a;
    }
    return result;
  }

  /// vector-matrix multiplication
  inline friend vector1d<T> operator * (vector1d<T> const &v,
                                        matrix2d<T> const &m)
  {
    vector1d<T> result(m.inner_length);
    if (m.outer_length != v.size()) {
      cvm::error("Error: trying to multiply a vector and a matrix "
                 "of incompatible sizes, "+
                  cvm::to_str(v.size()) + " and " +
                 cvm::to_str(m.outer_length)+"x"+cvm::to_str(m.inner_length) +
                 ".\n");
    } else {
      size_t i, k;
      for (i = 0; i < m.inner_length; i++) {
        for (k = 0; k < m.outer_length; k++) {
          result[i] += m[k][i] * v[k];
        }
      }
    }
    return result;
  }

  /// Formatted output
  friend std::ostream & operator << (std::ostream &os,
                                     matrix2d<T> const &m)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    os.width(2);
    os << "( ";
    size_t i;
    for (i = 0; i < m.outer_length; i++) {
      os << " ( ";
      size_t j;
      for (j = 0; j < m.inner_length-1; j++) {
        os.width(w);
        os.precision(p);
        os << m[i][j] << " , ";
      }
      os.width(w);
      os.precision(p);
      os << m[i][m.inner_length-1] << " )";
    }

    os << " )";
    return os;
  }

  inline std::string to_simple_string() const
  {
    if (this->size() == 0) return std::string("");
    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(cvm::cv_prec);
    os << (*this)[0];
    size_t i;
    for (i = 1; i < data.size(); i++) {
      os << " " << data[i];
    }
    return os.str();
  }

  inline int from_simple_string(std::string const &s)
  {
    std::stringstream stream(s);
    size_t i = 0;
    while ((i < data.size()) && (stream >> data[i])) {
      i++;
    }
    if (i < data.size()) {
      return COLVARS_ERROR;
    }
    return COLVARS_OK;
  }

};


/// vector of real numbers with three components
class colvarmodule::rvector {

public:

  cvm::real x, y, z;

  inline rvector()
  {
    reset();
  }

  /// \brief Set all components to zero
  inline void reset()
  {
    set(0.0);
  }

  inline rvector(cvm::real x_i, cvm::real y_i, cvm::real z_i)
  {
    set(x_i, y_i, z_i);
  }

  inline rvector(cvm::vector1d<cvm::real> const &v)
  {
    set(v[0], v[1], v[2]);
  }

  inline rvector(cvm::real t)
  {
    set(t);
  }

  /// \brief Set all components to a scalar
  inline void set(cvm::real value)
  {
    x = y = z = value;
  }

  /// \brief Assign all components
  inline void set(cvm::real x_i, cvm::real y_i, cvm::real z_i)
  {
    x = x_i;
    y = y_i;
    z = z_i;
  }

  /// \brief Access cartesian components by index
  inline cvm::real & operator [] (int i) {
    return (i == 0) ? x : (i == 1) ? y : (i == 2) ? z : x;
  }

  /// \brief Access cartesian components by index
  inline cvm::real  operator [] (int i) const {
    return (i == 0) ? x : (i == 1) ? y : (i == 2) ? z : x;
  }

  inline cvm::vector1d<cvm::real> const as_vector() const
  {
    cvm::vector1d<cvm::real> result(3);
    result[0] = this->x;
    result[1] = this->y;
    result[2] = this->z;
    return result;
  }

  inline void operator += (cvm::rvector const &v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
  }

  inline void operator -= (cvm::rvector const &v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }

  inline void operator *= (cvm::real v)
  {
    x *= v;
    y *= v;
    z *= v;
  }

  inline void operator /= (cvm::real const& v)
  {
    x /= v;
    y /= v;
    z /= v;
  }

  inline cvm::real norm2() const
  {
    return (x*x + y*y + z*z);
  }

  inline cvm::real norm() const
  {
    return cvm::sqrt(this->norm2());
  }

  inline cvm::rvector unit() const
  {
    const cvm::real n = this->norm();
    return (n > 0. ? cvm::rvector(x, y, z)/n : cvm::rvector(1., 0., 0.));
  }

  static inline size_t output_width(size_t real_width)
  {
    return 3*real_width + 10;
  }


  static inline cvm::rvector outer(cvm::rvector const &v1,
                                   cvm::rvector const &v2)
  {
    return cvm::rvector( v1.y*v2.z - v2.y*v1.z,
                         -v1.x*v2.z + v2.x*v1.z,
                         v1.x*v2.y - v2.x*v1.y);
  }

  friend inline cvm::rvector operator - (cvm::rvector const &v)
  {
    return cvm::rvector(-v.x, -v.y, -v.z);
  }

  friend inline cvm::rvector operator + (cvm::rvector const &v1,
                                         cvm::rvector const &v2)
  {
    return cvm::rvector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
  }
  friend inline cvm::rvector operator - (cvm::rvector const &v1,
                                         cvm::rvector const &v2)
  {
    return cvm::rvector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
  }

  /// Inner (dot) product
  friend inline cvm::real operator * (cvm::rvector const &v1,
                                      cvm::rvector const &v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  friend inline cvm::rvector operator * (cvm::real a, cvm::rvector const &v)
  {
    return cvm::rvector(a*v.x, a*v.y, a*v.z);
  }

  friend inline cvm::rvector operator * (cvm::rvector const &v, cvm::real a)
  {
    return cvm::rvector(a*v.x, a*v.y, a*v.z);
  }

  friend inline cvm::rvector operator / (cvm::rvector const &v, cvm::real a)
  {
    return cvm::rvector(v.x/a, v.y/a, v.z/a);
  }

  std::string to_simple_string() const;
  int from_simple_string(std::string const &s);
};


/// \brief 2-dimensional array of real numbers with three components
/// along each dimension (works with colvarmodule::rvector)
class colvarmodule::rmatrix {

public:

  cvm::real xx, xy, xz, yx, yy, yz, zx, zy, zz;

  /// Default constructor
  inline rmatrix()
  {
    reset();
  }

  /// Constructor component by component
  inline rmatrix(cvm::real xxi, cvm::real xyi, cvm::real xzi,
                 cvm::real yxi, cvm::real yyi, cvm::real yzi,
                 cvm::real zxi, cvm::real zyi, cvm::real zzi)
  {
    xx = xxi;
    xy = xyi;
    xz = xzi;
    yx = yxi;
    yy = yyi;
    yz = yzi;
    zx = zxi;
    zy = zyi;
    zz = zzi;
  }


  inline void reset()
  {
    xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;
  }

  /// Return the determinant
  inline cvm::real determinant() const
  {
    return
      (  xx * (yy*zz - zy*yz))
      - (yx * (xy*zz - zy*xz))
      + (zx * (xy*yz - yy*xz));
  }

  inline cvm::rmatrix transpose() const
  {
    return cvm::rmatrix(xx, yx, zx,
                        xy, yy, zy,
                        xz, yz, zz);
  }

  inline friend cvm::rvector operator * (cvm::rmatrix const &m,
                                         cvm::rvector const &r)
  {
    return cvm::rvector(m.xx*r.x + m.xy*r.y + m.xz*r.z,
                        m.yx*r.x + m.yy*r.y + m.yz*r.z,
                        m.zx*r.x + m.zy*r.y + m.zz*r.z);
  }
};



/// \brief 1-dimensional vector of real numbers with four components and
/// a quaternion algebra
class colvarmodule::quaternion {

public:

  cvm::real q0, q1, q2, q3;

  /// Constructor component by component
  inline quaternion(cvm::real const qv[4])
    : q0(qv[0]), q1(qv[1]), q2(qv[2]), q3(qv[3])
  {}

  /// Constructor component by component
  inline quaternion(cvm::real q0i,
                    cvm::real q1i,
                    cvm::real q2i,
                    cvm::real q3i)
    : q0(q0i), q1(q1i), q2(q2i), q3(q3i)
  {}

  inline quaternion(cvm::vector1d<cvm::real> const &v)
    : q0(v[0]), q1(v[1]), q2(v[2]), q3(v[3])
  {}

  /// \brief Default constructor
  inline quaternion()
  {
    reset();
  }

  /// \brief Set all components to zero (null quaternion)
  inline void reset()
  {
    q0 = q1 = q2 = q3 = 0.0;
  }

  /// Tell the number of characters required to print a quaternion, given that of a real number
  static inline size_t output_width(size_t real_width)
  {
    return 4*real_width + 13;
  }

  std::string to_simple_string() const;
  int from_simple_string(std::string const &s);

  /// \brief Formatted output operator
  friend std::ostream & operator << (std::ostream &os, cvm::quaternion const &q);
  /// \brief Formatted input operator
  friend std::istream & operator >> (std::istream &is, cvm::quaternion &q);

  /// Access the quaternion as a 4-d array (return a reference)
  inline cvm::real & operator [] (int i) {
    switch (i) {
    case 0:
      return this->q0;
    case 1:
      return this->q1;
    case 2:
      return this->q2;
    case 3:
      return this->q3;
    default:
      cvm::error("Error: incorrect quaternion component.\n");
      return q0;
    }
  }

  /// Access the quaternion as a 4-d array (return a value)
  inline cvm::real operator [] (int i) const {
    switch (i) {
    case 0:
      return this->q0;
    case 1:
      return this->q1;
    case 2:
      return this->q2;
    case 3:
      return this->q3;
    default:
      cvm::error("Error: trying to access a quaternion "
                 "component which is not between 0 and 3.\n");
      return 0.0;
    }
  }

  inline cvm::vector1d<cvm::real> const as_vector() const
  {
    cvm::vector1d<cvm::real> result(4);
    result[0] = q0;
    result[1] = q1;
    result[2] = q2;
    result[3] = q3;
    return result;
  }

  /// Square norm of the quaternion
  inline cvm::real norm2() const
  {
    return q0*q0 + q1*q1 + q2*q2 + q3*q3;
  }

  /// Norm of the quaternion
  inline cvm::real norm() const
  {
    return cvm::sqrt(this->norm2());
  }

  /// Return the conjugate quaternion
  inline cvm::quaternion conjugate() const
  {
    return cvm::quaternion(q0, -q1, -q2, -q3);
  }

  inline void operator *= (cvm::real a)
  {
    q0 *= a; q1 *= a; q2 *= a; q3 *= a;
  }

  inline void operator /= (cvm::real a)
  {
    q0 /= a; q1 /= a; q2 /= a; q3 /= a;
  }

  inline void set_positive()
  {
    if (q0 > 0.0) return;
    q0 = -q0;
    q1 = -q1;
    q2 = -q2;
    q3 = -q3;
  }

  inline void operator += (cvm::quaternion const &h)
  {
    q0+=h.q0; q1+=h.q1; q2+=h.q2; q3+=h.q3;
  }
  inline void operator -= (cvm::quaternion const &h)
  {
    q0-=h.q0; q1-=h.q1; q2-=h.q2; q3-=h.q3;
  }

  /// Return the vector component
  inline cvm::rvector get_vector() const
  {
    return cvm::rvector(q1, q2, q3);
  }


  friend inline cvm::quaternion operator + (cvm::quaternion const &h,
                                            cvm::quaternion const &q)
  {
    return cvm::quaternion(h.q0+q.q0, h.q1+q.q1, h.q2+q.q2, h.q3+q.q3);
  }

  friend inline cvm::quaternion operator - (cvm::quaternion const &h,
                                            cvm::quaternion const &q)
  {
    return cvm::quaternion(h.q0-q.q0, h.q1-q.q1, h.q2-q.q2, h.q3-q.q3);
  }

  /// \brief Provides the quaternion product.  \b NOTE: for the inner
  /// product use: `h.inner (q);`
  friend inline cvm::quaternion operator * (cvm::quaternion const &h,
                                            cvm::quaternion const &q)
  {
    return cvm::quaternion(h.q0*q.q0 - h.q1*q.q1 - h.q2*q.q2 - h.q3*q.q3,
                           h.q0*q.q1 + h.q1*q.q0 + h.q2*q.q3 - h.q3*q.q2,
                           h.q0*q.q2 + h.q2*q.q0 + h.q3*q.q1 - h.q1*q.q3,
                           h.q0*q.q3 + h.q3*q.q0 + h.q1*q.q2 - h.q2*q.q1);
  }

  friend inline cvm::quaternion operator * (cvm::real c,
                                            cvm::quaternion const &q)
  {
    return cvm::quaternion(c*q.q0, c*q.q1, c*q.q2, c*q.q3);
  }
  friend inline cvm::quaternion operator * (cvm::quaternion const &q,
                                            cvm::real c)
  {
    return cvm::quaternion(q.q0*c, q.q1*c, q.q2*c, q.q3*c);
  }
  friend inline cvm::quaternion operator / (cvm::quaternion const &q,
                                            cvm::real c)
  {
    return cvm::quaternion(q.q0/c, q.q1/c, q.q2/c, q.q3/c);
  }


  /// \brief Rotate v through this quaternion (put it in the rotated
  /// reference frame)
  inline cvm::rvector rotate(cvm::rvector const &v) const
  {
    return ( (*this) * cvm::quaternion(0.0, v.x, v.y, v.z) *
             this->conjugate() ).get_vector();
  }

  /// \brief Rotate Q2 through this quaternion (put it in the rotated
  /// reference frame)
  inline cvm::quaternion rotate(cvm::quaternion const &Q2) const
  {
    cvm::rvector const vq_rot = this->rotate(Q2.get_vector());
    return cvm::quaternion(Q2.q0, vq_rot.x, vq_rot.y, vq_rot.z);
  }

  /// Return the 3x3 matrix associated to this quaternion
  inline cvm::rmatrix rotation_matrix() const
  {
    cvm::rmatrix R;

    R.xx = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R.yy = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    R.zz = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    R.xy = 2.0 * (q1*q2 - q0*q3);
    R.xz = 2.0 * (q0*q2 + q1*q3);

    R.yx = 2.0 * (q0*q3 + q1*q2);
    R.yz = 2.0 * (q2*q3 - q0*q1);

    R.zx = 2.0 * (q1*q3 - q0*q2);
    R.zy = 2.0 * (q0*q1 + q2*q3);

    return R;
  }


  /// \brief Multiply the given vector by the derivative of the given
  /// (rotated) position with respect to the quaternion
  /// \param pos The position \f$\mathbf{x}\f$.
  /// \param vec The vector \f$\mathbf{v}\f$.
  /// \return A quaternion (see the detailed documentation below).
  ///
  /// This function is mainly used for projecting the gradients or forces on
  /// the rotated atoms to the forces on quaternion. Assume this rotation can
  /// be represented as \f$R(\mathbf{q})\f$,
  /// where \f$\mathbf{q} := (q_0, q_1, q_2, q_3)\f$
  /// is the current quaternion, the function returns the following new
  /// quaternion:
  /// \f[
  /// \left(\mathbf{v}^\mathrm{T}\frac{\partial R(\mathbf{q})}{\partial q_0}\mathbf{x},
  ///       \mathbf{v}^\mathrm{T}\frac{\partial R(\mathbf{q})}{\partial q_1}\mathbf{x},
  ///       \mathbf{v}^\mathrm{T}\frac{\partial R(\mathbf{q})}{\partial q_2}\mathbf{x},
  ///       \mathbf{v}^\mathrm{T}\frac{\partial R(\mathbf{q})}{\partial q_3}\mathbf{x}\right)
  /// \f]
  /// where \f$\mathbf{v}\f$ is usually the gradient of \f$\xi\f$ with respect to
  /// the rotated frame \f$\tilde{\mathbf{X}}\f$,
  /// \f$\partial \xi / \partial \tilde{\mathbf{X}}\f$, or the force acting on it
  /// (\f$\mathbf{F}_{\tilde{\mathbf{X}}}\f$).
  /// By using the following loop in pseudo C++ code,
  /// either \f$\partial \xi / \partial \tilde{\mathbf{X}}\f$
  /// or \f$\mathbf{F}_{\tilde{\mathbf{X}}}\f$, can be projected to
  /// \f$\partial \xi / \partial \mathbf{q}\f$ or \f$\mathbf{F}_q\f$ into `sum_dxdq`:
  /// @code
  /// cvm::real sum_dxdq[4] = {0, 0, 0, 0};
  /// for (size_t i = 0; i < main_group_size(); ++i) {
  ///   const cvm::rvector v = grad_or_force_on_rotated_main_group(i);
  ///   const cvm::rvector x = unrotated_main_group_positions(i);
  ///   cvm::quaternion const dxdq = position_derivative_inner(x, v);
  ///   sum_dxdq[0] += dxdq[0];
  ///   sum_dxdq[1] += dxdq[1];
  ///   sum_dxdq[2] += dxdq[2];
  ///   sum_dxdq[3] += dxdq[3];
  /// }
  /// @endcode
  inline cvm::quaternion position_derivative_inner(cvm::rvector const &pos,
                                            cvm::rvector const &vec) const {
    return cvm::quaternion(2.0 * (vec.x * ( q0 * pos.x - q3 * pos.y + q2 * pos.z) +
                                  vec.y * ( q3 * pos.x + q0 * pos.y - q1 * pos.z) +
                                  vec.z * (-q2 * pos.x + q1 * pos.y + q0 * pos.z)),
                           2.0 * (vec.x * ( q1 * pos.x + q2 * pos.y + q3 * pos.z) +
                                  vec.y * ( q2 * pos.x - q1 * pos.y - q0 * pos.z) +
                                  vec.z * ( q3 * pos.x + q0 * pos.y - q1 * pos.z)),
                           2.0 * (vec.x * (-q2 * pos.x + q1 * pos.y + q0 * pos.z) +
                                  vec.y * ( q1 * pos.x + q2 * pos.y + q3 * pos.z) +
                                  vec.z * (-q0 * pos.x + q3 * pos.y - q2 * pos.z)),
                           2.0 * (vec.x * (-q3 * pos.x - q0 * pos.y + q1 * pos.z) +
                                  vec.y * ( q0 * pos.x - q3 * pos.y + q2 * pos.z) +
                                  vec.z * ( q1 * pos.x + q2 * pos.y + q3 * pos.z)));
  }


  /// \brief Return the cosine between the orientation frame
  /// associated to this quaternion and another
  inline cvm::real cosine(cvm::quaternion const &q) const
  {
    cvm::real const iprod = this->inner(q);
    return 2.0*iprod*iprod - 1.0;
  }

  /// \brief Square distance from another quaternion on the
  /// 4-dimensional unit sphere: returns the square of the angle along
  /// the shorter of the two geodesics
  inline cvm::real dist2(cvm::quaternion const &Q2) const
  {
    cvm::real const cos_omega = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;

    cvm::real const omega = cvm::acos( (cos_omega > 1.0) ? 1.0 :
                                       ( (cos_omega < -1.0) ? -1.0 : cos_omega) );

    // get the minimum distance: x and -x are the same quaternion
    if (cos_omega > 0.0)
      return omega * omega;
    else
      return (PI-omega) * (PI-omega);
  }

  /// Gradient of the square distance: returns a 4-vector equivalent
  /// to that provided by slerp
  inline cvm::quaternion dist2_grad(cvm::quaternion const &Q2) const
  {
    cvm::real const cos_omega = this->q0*Q2.q0 + this->q1*Q2.q1 + this->q2*Q2.q2 + this->q3*Q2.q3;
    cvm::real const omega = cvm::acos( (cos_omega > 1.0) ? 1.0 :
                                       ( (cos_omega < -1.0) ? -1.0 : cos_omega) );
    cvm::real const sin_omega = cvm::sin(omega);

    if (cvm::fabs(sin_omega) < 1.0E-14) {
      // return a null 4d vector
      return cvm::quaternion(0.0, 0.0, 0.0, 0.0);
    }

    cvm::quaternion const
      grad1((-1.0)*sin_omega*Q2.q0 + cos_omega*(this->q0-cos_omega*Q2.q0)/sin_omega,
            (-1.0)*sin_omega*Q2.q1 + cos_omega*(this->q1-cos_omega*Q2.q1)/sin_omega,
            (-1.0)*sin_omega*Q2.q2 + cos_omega*(this->q2-cos_omega*Q2.q2)/sin_omega,
            (-1.0)*sin_omega*Q2.q3 + cos_omega*(this->q3-cos_omega*Q2.q3)/sin_omega);

    if (cos_omega > 0.0) {
      return 2.0*omega*grad1;
    } else {
      return -2.0*(PI-omega)*grad1;
    }
  }

  /// \brief Choose the closest between Q2 and -Q2 and save it back.
  /// Not required for dist2() and dist2_grad()
  inline void match(cvm::quaternion &Q2) const
  {
    cvm::real const cos_omega = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;
    if (cos_omega < 0.0) Q2 *= -1.0;
  }

  /// \brief Inner product (as a 4-d vector) with Q2; requires match()
  /// if the largest overlap is looked for
  inline cvm::real inner(cvm::quaternion const &Q2) const
  {
    cvm::real const prod = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;
    return prod;
  }


};

#ifndef COLVARS_LAMMPS
namespace NR {
void diagonalize_matrix(cvm::real m[4][4],
                        cvm::real eigval[4],
                        cvm::real eigvec[4][4]);
}
#endif


/// \brief A rotation between two sets of coordinates (for the moment
/// a wrapper for colvarmodule::quaternion)
class colvarmodule::rotation
{
private:
  /// Correlation matrix C (3, 3)
  cvm::rmatrix C;

  /// Overlap matrix S (4, 4)
  cvm::real S[4][4];

  /// Eigenvalues of S
  cvm::real S_eigval[4];

  /// Eigenvectors of S
  cvm::real S_eigvec[4][4];

  /// Used for debugging gradients
  cvm::real S_backup[4][4];

public:
  /// \brief Perform gradient tests
  bool b_debug_gradients;

  /// \brief The rotation itself (implemented as a quaternion)
  cvm::quaternion q{1.0, 0.0, 0.0, 0.0};

  template <typename T1, typename T2>
  friend struct rotation_derivative;

  template<typename T1, typename T2>
  friend void debug_gradients(
    cvm::rotation &rot,
    const std::vector<T1> &pos1,
    const std::vector<T2> &pos2);

  /// \brief Calculate the optimal rotation and store the
  /// corresponding eigenvalue and eigenvector in the arguments l0 and
  /// q0; if the gradients have been previously requested, calculate
  /// them as well
  ///
  /// The method to derive the optimal rotation is defined in:
  /// Coutsias EA, Seok C, Dill KA.
  /// Using quaternions to calculate RMSD.
  /// J Comput Chem. 25(15):1849-57 (2004)
  /// DOI: 10.1002/jcc.20110  PubMed: 15376254
  void calc_optimal_rotation(std::vector<atom_pos> const &pos1,
                             std::vector<atom_pos> const &pos2);
  void calc_optimal_rotation(std::vector<cvm::atom> const &pos1,
                             std::vector<atom_pos> const &pos2);

  /// Initialize member data
  int init();

  /// Default constructor
  rotation();

  /// Constructor after a quaternion
  rotation(cvm::quaternion const &qi);

  /// Constructor after an axis of rotation and an angle (in radians)
  rotation(cvm::real angle, cvm::rvector const &axis);

  /// Destructor
  ~rotation();

  /// Return the rotated vector
  inline cvm::rvector rotate(cvm::rvector const &v) const
  {
    return q.rotate(v);
  }

  /// Return the inverse of this rotation
  inline cvm::rotation inverse() const
  {
    return cvm::rotation(this->q.conjugate());
  }

  /// Return the associated 3x3 matrix
  inline cvm::rmatrix matrix() const
  {
    return q.rotation_matrix();
  }

  /// \brief Return the spin angle (in degrees) with respect to the
  /// provided axis (which MUST be normalized)
  inline cvm::real spin_angle(cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real alpha = (180.0/PI) * 2.0 * cvm::atan2(axis * q_vec, q.q0);
    while (alpha >  180.0) alpha -= 360;
    while (alpha < -180.0) alpha += 360;
    return alpha;
  }

  /// \brief Return the derivative of the spin angle with respect to
  /// the quaternion
  inline cvm::quaternion dspin_angle_dq(cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real const iprod = axis * q_vec;

    if (q.q0 != 0.0) {

      cvm::real const dspindx =
        (180.0/PI) * 2.0 * (1.0 / (1.0 + (iprod*iprod)/(q.q0*q.q0)));

      return cvm::quaternion( dspindx * (iprod * (-1.0) / (q.q0*q.q0)),
                              dspindx * ((1.0/q.q0) * axis.x),
                              dspindx * ((1.0/q.q0) * axis.y),
                              dspindx * ((1.0/q.q0) * axis.z));
    } else {
      // (1/(1+x^2)) ~ (1/x)^2
      // The documentation of spinAngle discourages its use when q_vec and
      // axis are not close
      return cvm::quaternion((180.0/PI) * 2.0 * ((-1.0)/iprod), 0.0, 0.0, 0.0);
    }
  }

  /// \brief Return the projection of the orientation vector onto a
  /// predefined axis
  inline cvm::real cos_theta(cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real const alpha =
      (180.0/PI) * 2.0 * cvm::atan2(axis * q_vec, q.q0);

    cvm::real const cos_spin_2 = cvm::cos(alpha * (PI/180.0) * 0.5);
    cvm::real const cos_theta_2 = ( (cos_spin_2 != 0.0) ?
                                    (q.q0 / cos_spin_2) :
                                    (0.0) );
    // cos(2t) = 2*cos(t)^2 - 1
    return 2.0 * (cos_theta_2*cos_theta_2) - 1.0;
  }

  /// Return the derivative of the tilt wrt the quaternion
  inline cvm::quaternion dcos_theta_dq(cvm::rvector const &axis) const
  {
    cvm::rvector const q_vec = q.get_vector();
    cvm::real const iprod = axis * q_vec;

    cvm::real const cos_spin_2 = cvm::cos(cvm::atan2(iprod, q.q0));

    if (q.q0 != 0.0)  {

      cvm::real const d_cos_theta_dq0 =
        (4.0 * q.q0 / (cos_spin_2*cos_spin_2)) *
        (1.0 - (iprod*iprod)/(q.q0*q.q0) / (1.0 + (iprod*iprod)/(q.q0*q.q0)));

      cvm::real const d_cos_theta_dqn =
        (4.0 * q.q0 / (cos_spin_2*cos_spin_2) *
         (iprod/q.q0) / (1.0 + (iprod*iprod)/(q.q0*q.q0)));

      return cvm::quaternion(d_cos_theta_dq0,
                             d_cos_theta_dqn * axis.x,
                             d_cos_theta_dqn * axis.y,
                             d_cos_theta_dqn * axis.z);
    } else {

      cvm::real const d_cos_theta_dqn =
        (4.0 / (cos_spin_2*cos_spin_2 * iprod));

      return cvm::quaternion(0.0,
                             d_cos_theta_dqn * axis.x,
                             d_cos_theta_dqn * axis.y,
                             d_cos_theta_dqn * axis.z);
    }
  }

  /// \brief Whether to test for eigenvalue crossing
  static bool monitor_crossings;
  /// \brief Threshold for the eigenvalue crossing test
  static cvm::real crossing_threshold;

protected:

  /// \brief Previous value of the rotation (used to warn the user
  /// when the structure changes too much, and there may be an
  /// eigenvalue crossing)
  cvm::quaternion q_old;

  /// Build the correlation matrix C (used by calc_optimal_rotation())
  void build_correlation_matrix(std::vector<cvm::atom_pos> const &pos1,
                                std::vector<cvm::atom_pos> const &pos2);
  void build_correlation_matrix(std::vector<cvm::atom> const &pos1,
                                std::vector<cvm::atom_pos> const &pos2);

  /// \brief Actual implementation of `calc_optimal_rotation` (and called by it)
  void calc_optimal_rotation_impl();

  /// Compute the overlap matrix S (used by calc_optimal_rotation())
  void compute_overlap_matrix();

  /// Pointer to instance of Jacobi solver
  void *jacobi;
};


#endif
