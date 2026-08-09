//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef MAPPED_VECTOR_C
#define MAPPED_VECTOR_C


#include <stdexcept>
#include <vector>

#include "mapped_vector.h"


using std::out_of_range;
using std::vector;


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::mapped_vector()
{

    _current.first.clear();
    _current.second = 0;

}


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::mapped_vector(const StringCompareT& cmp)
  : _index(cmp)
{

    _current.first.clear();
    _current.second = 0;

}


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::mapped_vector(
  const mapped_vector& inMappedVector)
{

    _index = inMappedVector._index;
    _vector = inMappedVector._vector;

}


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::~mapped_vector()
{

    clear();

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::push_back(const T& inT)
{

    _vector.push_back(inT);

    typename tIndex::value_type valuePair(inT, _vector.size() - 1);

    _index.insert(valuePair);

    _current.first = inT;
    _current.second = _vector.size() - 1;

}


template <typename T, typename StringCompareT>
unsigned int mapped_vector<T, StringCompareT>::size() const
{

    return(_vector.size());

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::empty() const
{

    return(_vector.empty());

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::operator=(
  const mapped_vector& inMappedVector)
{

    _index = inMappedVector._index;
    _vector = inMappedVector._vector;

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::operator=(const vector<T>& inVector)
{

    clear();

    for (unsigned int index = 0; index < inVector.size(); ++index)
    {
        push_back(inVector[index]);
    }

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::operator==(
  const mapped_vector& inMappedVector)
{

    return(_vector == inMappedVector._vector);

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::operator!=(
  const mapped_vector& inMappedVector)
{

    return(!operator==(inMappedVector));

}


template <typename T, typename StringCompareT>
const T& mapped_vector<T, StringCompareT>::operator[](unsigned int index) const
{

    if (index >= size())
    {
        throw out_of_range("Invalid index in mapped_vector::operator[]");
    }

    return(_vector[index]);

}


template <typename T, typename StringCompareT>
const vector<T>& mapped_vector<T, StringCompareT>::get_vector() const
{

    return(_vector);

}


template <typename T, typename StringCompareT>
vector<T>& mapped_vector<T, StringCompareT>::get_vector()
{

    return(_vector);

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::erase(const T& inT)
{

    unsigned int index = get_index(inT);
    if (index >= size())
    {
        throw out_of_range("Element not found in mapped_vector::erase");
    }

    if (is_equal(_current.first, _vector[index]))
    {
        _current.first.clear();
        _current.second = 0;
    }

    _vector.erase(_vector.begin() + index);

    _index.erase(inT);

    for (typename tIndex::iterator pos = _index.begin();
      pos != _index.end(); ++pos)
    {
        if (pos->second >= index)
        {
            --(pos->second);
        }
    }

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::insert(const unsigned int index,
  const T& inT)
{

    unsigned int existingIndex = get_index(inT);
    if (existingIndex != size())
    {
        throw out_of_range("Element exists in mapped_vector::insert");
    }

    _current.first = inT;
    _current.second = index;

    _vector.insert(_vector.begin() + index, inT);

    for (typename tIndex::iterator pos = _index.begin(); pos != _index.end();
      ++pos)
    {
        if (pos->second >= index)
        {
            ++(pos->second);
        }
    }

    typename tIndex::value_type valuePair(inT, index);

    _index.insert(valuePair);

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::index_it()
{

    for (unsigned int index = 0; index < _vector.size(); ++index)
    {
        typename tIndex::value_type valuePair(_vector[index], index);

        _index.insert(valuePair);
    }

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::clear()
{

    _index.clear();
    _vector.clear();

    _current.first.clear();
    _current.second = 0;

}


template <typename T, typename StringCompareT>
unsigned int mapped_vector<T, StringCompareT>::find(const T& inT) const
{

    return(get_index(inT));

}


template <typename T, typename StringCompareT>
unsigned int mapped_vector<T, StringCompareT>::get_index(const T& inT) const
{

    if (is_equal(_current.first, inT))
    {
        return(_current.second);
    }

    // Return index of found value or invalid index
    typename tIndex::const_iterator pos = _index.find(inT);
    if (pos != _index.end())
    {
        // Found
        _current.first = inT;
        _current.second = pos->second;

        return(pos->second);
    }
    else
    {
        // Not found. Return invalid index.
        return(_vector.size());
    }

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::is_equal(const T& firstT,
  const T& secondT) const
{

    typename tIndex::key_compare keyComp = _index.key_comp();

    return(!(keyComp(firstT, secondT) || keyComp(secondT, firstT)));

}


#endif

