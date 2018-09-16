#include <iterator>

namespace gmx
{

namespace detail
{
//! Base class for forward iterator adapters to reduce boilerplate
template<typename T>
class ForwardIteratorAdapterBase
{
public:
    // IteratorBase traits
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using value_type      = typename std::iterator_traits<T>::value_type;
    using pointer         = typename std::iterator_traits<T>::pointer;
    using reference       = typename std::iterator_traits<T>::reference;
    using iterator_category = std::forward_iterator_tag;
    bool operator==(const ForwardIteratorAdapterBase &other) const {return it_ == other.it_;}
    bool operator!=(const ForwardIteratorAdapterBase &other) const {return !(*this == other);}
    reference operator*() {return *it_;}
protected:
    ForwardIteratorAdapterBase(T it) : it_(it) {}
    T it_;
};

} // namespace detail
   
//! Filter iterator. Takes start end of range and predicacte.
// Can only be Foward iterator because we need to store end and have to avoid
// skipping past the end. Predicate is stored as base close for empty base class
// optimizatin.
template<typename It, typename F>
class FilterIterator : public detail::ForwardIteratorAdapterBase<It>, private F
{
    using Base = detail::ForwardIteratorAdapterBase<It>;
public:
    FilterIterator(It it, It end, F f) : Base(it), F(f), end_(end) {}
    FilterIterator& operator++()
    {
        do { ++Base::it_; } while (Base::it_!=end_ && !(*this)(*Base::it_));
        return *this;
    }
    FilterIterator operator++(int) {FilterIterator ret = *this; ++(*this); return ret;}
private:
    It end_;
};

template<typename ItBegin, typename ItEnd>
class Range
{
public:
    Range(ItBegin begin, ItEnd end) : begin_(begin), end_(end) {}
    ItBegin begin() { return begin_; }
    ItEnd end() { return end_; }
private:
    ItBegin begin_;
    ItEnd end_;
};

template<typename R, typename F,
         typename It = FilterIterator<typename R::const_iterator, F>>
Range<It, It> makeFilteredView(const R &r, F f)
{
    It begin = It(r.begin(), r.end(), f);
    if (!f(*begin)) {++begin;}
    //Future: End iterator should be sentinel but requires C++17 to be
    //compatible with range loop.
    return {begin, It(r.end(), r.end(), f)};
}

} // namespace gmx
