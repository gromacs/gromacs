
#include <type_traits>

template<class T>
class atomic_ref {};

template<class T>
struct accessor_atomic {

  static_assert( std::is_trivially_copyable<T>::value );

	using element_type = T;
	using reference    = atomic_ref<T> ;
	using pointer      = T*;
	using offset       = accessor_atomic ;

  constexpr accessor_atomic() noexcept = default ;
  constexpr accessor_atomic( accessor_atomic && ) noexcept = default ;
  constexpr accessor_atomic( const accessor_atomic & ) noexcept = default ;
  accessor_atomic & operator =( accessor_atomic && ) noexcept = default ;
  accessor_atomic & operator =( const accessor_atomic & ) noexcept = default ;

  explicit constexpr accessor_atomic( pointer other ) noexcept
	  : ptr(other)
		{ assert( 0 == reinterpret_cast<uintptr_t>(ptr) % reference::required_alignment ); };

  constexpr reference operator[]( size_t i ) const noexcept
	  { return reference( ptr[i] ); }

  constexpr offset operator+( size_t i ) const noexcept
	  { return offset(ptr+i); }

  constexpr operator pointer() const
	  { assert(false /* cannot access raw data outside of atomic */); }
};

