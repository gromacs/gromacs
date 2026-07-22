//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/**
** \file mapped_ptr_vector.h
*/


#ifndef MAPPED_PTR_VECTOR_H
#define MAPPED_PTR_VECTOR_H


#include <string>
#include <vector>
#include <map>


/**
** This is a container of pointers to objects. The container maintains the
** order of the inserted elements (as vector does), but it provides for
** efficient element access, search, serialization and deserialization.
** Object names must be unique, i.e., for any two object names in the container
** operator==() must yield false.
*/
template <typename T, typename StringCompareT = std::less<std::string> >
class mapped_ptr_vector
{
  private:
    // The first integer is the index in _vector
    // The second integer is the index in the file
    typedef std::map<std::string, std::pair<unsigned int, unsigned int>,
      StringCompareT > tIndex;

    tIndex _index;

    std::vector<T*> _vector;

    std::string _currentName;
    std::pair<unsigned int, unsigned int> _currentIndices;

    bool is_equal(const std::string& first, const std::string& second,
      const typename tIndex::key_compare& keyComp) const;

  public:
    mapped_ptr_vector();
    mapped_ptr_vector(const StringCompareT& cmp);
    mapped_ptr_vector(const mapped_ptr_vector& inMappedVector);
    ~mapped_ptr_vector();

    void operator=(const mapped_ptr_vector& inMappedVector);

    unsigned int size() const;
    bool empty() const;
    void clear();

    bool operator==(const mapped_ptr_vector& inMappedVector);

    void push_back(T* inP, const unsigned int fileIndex = 0);
    void push_back(const std::string& name, const unsigned int fileIndex = 0);
    void push_back(const std::vector<std::string>& names,
      const std::vector<unsigned int>& fileIndices);
    void push_back(const std::vector<std::string>& names);

    /// Associate the object pointer to already entered object name
    void set(T* inP);

    T& operator[](unsigned int index);
    T& operator[](const std::string& name);

    // When not found, returns size()
    unsigned int find(const std::string& name);

    void rename(const std::string& oldName, const std::string& newName);

    /// Removes object pointer
    void erase(const std::string& name);

    /// Is object de-serialized
    bool is_read(const std::string& name);

    /// De-serialize the object
    void read(const std::string& name);

    /// Serialize the object
    unsigned int write(const std::string& name);

    std::pair<unsigned int, unsigned int> get_indices(const std::string& name);
    std::string get_name(const unsigned int index);
    void get_sorted_indices(std::vector<unsigned int>& sortedIndices);
};


#endif
