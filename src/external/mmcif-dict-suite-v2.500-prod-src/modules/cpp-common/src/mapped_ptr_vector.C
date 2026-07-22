//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef MAPPED_PTR_VECTOR_C
#define MAPPED_PTR_VECTOR_C


#include <stdexcept>
#include <string>
#include <vector>

#include "Exceptions.h"
#include "mapped_ptr_vector.h"


using std::out_of_range;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;


template <typename T, typename StringCompareT>
mapped_ptr_vector<T, StringCompareT>::mapped_ptr_vector()
{


}


template <typename T, typename StringCompareT>
mapped_ptr_vector<T, StringCompareT>::mapped_ptr_vector(
  const StringCompareT& cmp) : _index(cmp)
{


}


template <typename T, typename StringCompareT>
mapped_ptr_vector<T, StringCompareT>::mapped_ptr_vector(
  const mapped_ptr_vector& inMappedPtrVector)
{

    _vector = inMappedPtrVector._vector;
    _index = inMappedPtrVector._index;
    _currentName = inMappedPtrVector._currentName;
    _currentIndices = inMappedPtrVector._currentIndices;

}


template <typename T, typename StringCompareT>
mapped_ptr_vector<T, StringCompareT>::~mapped_ptr_vector()
{

    clear();

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::operator=(const mapped_ptr_vector& inMappedPtrVector)
{

    _vector = inMappedPtrVector._vector;
    _index = inMappedPtrVector._index;
    _currentName = inMappedPtrVector._currentName;
    _currentIndices = inMappedPtrVector._currentIndices;

}


template <typename T, typename StringCompareT>
unsigned int mapped_ptr_vector<T, StringCompareT>::size() const
{

    return(_vector.size());

}


template <typename T, typename StringCompareT>
bool mapped_ptr_vector<T, StringCompareT>::empty() const
{

    return(_vector.empty());

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::clear()
{

    _vector.clear();

    _index.clear();

    _currentName.clear();

}


template <typename T, typename StringCompareT>
bool mapped_ptr_vector<T, StringCompareT>::operator==(
  const mapped_ptr_vector& inMappedPtrVector)
{

    return(_vector == inMappedPtrVector._vector);

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::push_back(T* inP,
  const unsigned int fileIndex)
{

    if (inP == NULL)
    {
        throw EmptyValueException("NULL vector",
          "mapped_ptr_vector::push_back");
    }

    _vector.push_back(inP);

    typename tIndex::value_type valuePair(inP->GetName(),
      make_pair(_vector.size() - 1, fileIndex));

    _index.insert(valuePair);

    _currentName = inP->GetName();
    _currentIndices = make_pair(_vector.size() - 1, fileIndex);

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::push_back(const string& name,
  const unsigned int fileIndex)
{

    _vector.push_back(NULL);

    typename tIndex::value_type valuePair(name,
      make_pair(_vector.size() - 1, fileIndex));

    _index.insert(valuePair);

    _currentName = name;
    _currentIndices = make_pair(_vector.size() - 1, fileIndex);
}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::push_back(
  const vector<string>& names, const vector<unsigned int>& fileIndices)
{

    if (names.size() != fileIndices.size())
    {
        throw out_of_range("Different sizes of names and fileIndices in"\
          " mapped_ptr_vector::push_back");
    }

    for (unsigned int nameI = 0; nameI < names.size(); ++nameI)
    {
        push_back(names[nameI], fileIndices[nameI]);
    }

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::push_back(const vector<string>& names)
{

    for (unsigned int nameI = 0; nameI < names.size(); ++nameI)
    {
        push_back(names[nameI], 0);
    }

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::set(T* inP)
{

    if (inP == NULL)
    {
        throw EmptyValueException("NULL vector",
          "mapped_ptr_vector::set");
    }

    pair<unsigned int, unsigned int> indices = get_indices(inP->GetName());

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::set");
    }

    _vector[indices.first] = inP;

}


template <typename T, typename StringCompareT>
T& mapped_ptr_vector<T, StringCompareT>::operator[](unsigned int index)
{

    if (index >= _vector.size())
    {
        throw out_of_range("Invalid index in"\
          " mapped_ptr_vector::operator[]");
    }

    return((T&)(*(_vector[index])));

}


template <typename T, typename StringCompareT>
T& mapped_ptr_vector<T, StringCompareT>::operator[](const string& name)
{

    pair<unsigned int, unsigned int> indices = get_indices(name);

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::operator[]");
    }

    return((T&)(*(_vector[indices.first])));

}


template <typename T, typename StringCompareT>
unsigned int mapped_ptr_vector<T, StringCompareT>::find(const string& name)
{

    pair<unsigned int, unsigned int> indices = get_indices(name);

    return(indices.first);

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::rename(const string& oldName,
  const string& newName)
{

    pair<unsigned int, unsigned int> indices = get_indices(oldName);

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::rename");
    }

    // Erase it from the map as it is about to change
    _index.erase(oldName);

    typename tIndex::value_type valuePair(newName, indices);

    _index.insert(valuePair);

    _vector[indices.first]->SetName(newName);

    typename tIndex::key_compare keyComp = _index.key_comp();

    if (is_equal(_currentName, oldName, keyComp))
    {
        _currentName = newName;
    }

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::erase(const string& name)
{

    pair<unsigned int, unsigned int> indices = get_indices(name);

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::erase");
    }

    _vector.erase(_vector.begin() + indices.first);

    _index.erase(name);
 
    // VLAD - PERFORMANCE - CAN THIS BE DONE MORE EFFICIENTLY USING FIND
    // AND NOT GOING THROUGH THE COMPLETE CONTAINER

    // Reduce by one all table indices greater than the table index of the
    // deleted table.
    for (typename tIndex::iterator pos = _index.begin(); pos != _index.end();
      ++pos)
    {
        if (pos->second.first > indices.first)
        {
            --(pos->second.first);
        }
    }

    typename tIndex::key_compare keyComp = _index.key_comp();

    if (is_equal(_currentName, name, keyComp))
    {
        _currentName.clear();
        _currentIndices = make_pair(_vector.size(), (unsigned int)0);
    }

}


template <typename T, typename StringCompareT>
bool mapped_ptr_vector<T, StringCompareT>::is_read(const string& name)
{

    pair<unsigned int, unsigned int> indices = get_indices(name);

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::is_read");
    }

    if (_vector[indices.first] != NULL)
    {
        return(true);
    }
    else
    {
        return(false);
    }

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::read(const string& name)
{

    pair<unsigned int, unsigned int> indices = get_indices(name);

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::read");
    }

    _vector[indices.first]->Read(indices.second);

}


template <typename T, typename StringCompareT>
unsigned int mapped_ptr_vector<T, StringCompareT>::write(const string& name)
{

    // VLAD TROUBLESHOOT POINT
    // VLAD PERFORMANCE

    pair<unsigned int, unsigned int> indices = get_indices(name);

    if (indices.first == _vector.size())
    {
        throw NotFoundException("Object not found",
          "mapped_ptr_vector::write");
    }

    // Erase it from the map as it is about to change
    _index.erase(name);

    indices.second = _vector[indices.first]->Write();

    typename tIndex::value_type valuePair(name, indices);

    _index.insert(valuePair);

    _currentName = name;
    _currentIndices = indices;

    return(indices.second);
}


template <typename T, typename StringCompareT>
pair<unsigned int, unsigned int> mapped_ptr_vector<T, StringCompareT>::get_indices(const string& name)
{

    if (_vector.empty())
    {
        // Empty container. Return invalid index.
        return(make_pair(_vector.size(), (unsigned int)0));
    }

    typename tIndex::key_compare keyComp = _index.key_comp();

    if (is_equal(name, _currentName, keyComp))
    {
        return(_currentIndices);
    }
    else
    {
        // Return index of found value or invalid index
        typename tIndex::iterator pos = _index.find(name);
        if (pos != _index.end())
        {
            // Update cache
            _currentName = name;
            _currentIndices = pos->second;
            // Found
            return(pos->second);
        }
        else
        {
            // Not found. Return invalid index.
            return(make_pair(_vector.size(), (unsigned int)0));
        } 
    }
    
}


template <typename T, typename StringCompareT>
string mapped_ptr_vector<T, StringCompareT>::get_name(const unsigned int index)
{

    if (index >= _vector.size())
    {
        throw out_of_range("Invalid index in"\
          " mapped_ptr_vector::get_name");
    }

    string ret;

    // Return index of found value or invalid index
    for (typename tIndex::iterator pos = _index.begin(); pos != _index.end();
      ++pos)
    if (pos->second.first == index)
    {
        ret = pos->first;
    }

    return(ret);

}


template <typename T, typename StringCompareT>
void mapped_ptr_vector<T, StringCompareT>::get_sorted_indices(
  vector<unsigned int>& sortedIndices)
{

    sortedIndices.clear();

    // Return index of found value or invalid index
    for (typename tIndex::iterator pos = _index.begin(); pos != _index.end();
      ++pos)
    {
        sortedIndices.push_back(pos->second.first);
    }

}


template <typename T, typename StringCompareT>
bool  mapped_ptr_vector<T, StringCompareT>::is_equal(const string& first,
  const string& second, const typename tIndex::key_compare& keyComp) const
{

    return(!(keyComp(first, second) || keyComp(second, first)));

}


#endif

