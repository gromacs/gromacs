//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef MAPPED_VECTOR_H
#define MAPPED_VECTOR_H


#include <vector>
#include <map>


/**
** Container of objects that maintans their order (as vector does), but
** provides for fast searching. Objects must be unique, i.e., for any two
** objects in the container operator==() must yield false.
*/
template <typename T, typename StringCompareT = std::less<T> >
class mapped_vector
{
  private:
    typedef std::map<T, unsigned int, StringCompareT> tIndex;

    tIndex _index;

    std::vector<T> _vector;

    mutable std::pair<T, unsigned int> _current;

    unsigned int get_index(const T& inT) const;

    bool is_equal(const T& firstT, const T& secondT) const;

  public:
    mapped_vector();
    mapped_vector(const StringCompareT& cmp);
    mapped_vector(const mapped_vector& inMappedVector);
    ~mapped_vector();

    void push_back(const T& inT);
    unsigned int size() const;
    bool empty() const;

    void operator=(const mapped_vector& inMappedVector);
    void operator=(const std::vector<T>& inVector);
    bool operator==(const mapped_vector& inMappedVector);
    bool operator!=(const mapped_vector& inMappedVector);

    const T& operator[](unsigned int index) const;
    const std::vector<T>& get_vector() const;
    std::vector<T>& get_vector();

    void erase(const T& inT);
    void insert(const unsigned int index, const T& inT);
    void index_it();

    void clear();

    /// When not found, returns size()
    unsigned int find(const T& inT) const;
};


#endif
