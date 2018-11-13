//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

namespace std {
namespace experimental {
inline namespace fundamentals_v3 {

// [mdspan.accessor.basic]
template<class ElementType>
class accessor_basic;


template<class ElementType>
class accessor_basic {
public:
  using element_type  = ElementType;
  using pointer       = ElementType*;
  using offset_policy = accessor_basic;
  using reference     = ElementType&;

  constexpr typename offset_policy::pointer
    offset( pointer p , ptrdiff_t i ) const noexcept
      { return typename offset_policy::pointer(p+i); }

  constexpr reference access( pointer p , ptrdiff_t i ) const noexcept
    { return p[i]; }

  constexpr ElementType* decay( pointer p ) const noexcept
    { return p; }
};

}}} // std::experimental::fundamentals_v3
