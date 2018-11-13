


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

namespace std {
namespace experimental {
inline namespace fundamentals_v3 {



// [mdspan.basic]
template<class ElementType,
         class Extents,
         class LayoutPolicy = layout_right,
         class AccessorPolicy = accessor_basic<ElementType> >
class basic_mdspan ;

// [msspan.subspan]

class all_type { public: constexpr explicit all_type() = default; };

/* inline */ constexpr all_type all ;

}}} // experimental::fundamentals_v3




//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

namespace std {
namespace experimental {
inline namespace fundamentals_v3 {

template<class ElementType, class Extents, class LayoutPolicy, class AccessorPolicy>
class basic_mdspan {
public:

  // Domain and codomain types

  using extents_type     = Extents ;
  using layout_type      = LayoutPolicy ;
  using accessor_type    = AccessorPolicy ;
  using mapping_type     = typename layout_type::template mapping<extents_type> ;
  using element_type     = typename accessor_type::element_type ;
  using value_type       = typename remove_cv<element_type>::type ;
  using index_type       = ptrdiff_t ;
  using difference_type  = ptrdiff_t ;
  using pointer          = typename accessor_type::pointer;
  using reference        = typename accessor_type::reference;

  // [mdspan.basic.cons]

  constexpr basic_mdspan() noexcept : acc_(), map_(), ptr_() {}

  constexpr basic_mdspan(basic_mdspan&& other) noexcept = default;

  constexpr basic_mdspan(const basic_mdspan & other) noexcept = default;

  basic_mdspan& operator=(const basic_mdspan & other) noexcept = default;

  basic_mdspan& operator=(basic_mdspan&& other) noexcept = default;

  template<class OtherElementType,
           class OtherExtents,
           class OtherLayoutPolicy,
           class OtherAccessor>
  constexpr basic_mdspan(
    const basic_mdspan<OtherElementType,
                       OtherExtents,
                       OtherLayoutPolicy,
                       OtherAccessor> & rhs ) noexcept
    : acc_( rhs.acc_ )
    , map_( rhs.map_ )
    , ptr_( rhs.ptr_ )
    {}

  template<class OtherElementType,
           class OtherExtents,
           class OtherLayoutPolicy,
           class OtherAccessor>
  basic_mdspan & operator = (
    const basic_mdspan<OtherElementType,
                       OtherExtents,
                       OtherLayoutPolicy,
                       OtherAccessor> & rhs ) noexcept
    { acc_ = rhs.acc_ ; map_ = rhs.map_ ; ptr_ = rhs.ptr_ ; return *this ; }

  template<class... IndexType >
  explicit constexpr basic_mdspan
    ( pointer ptr , IndexType ... DynamicExtents ) noexcept
    : acc_(accessor_type()), map_( extents_type(DynamicExtents...) ), ptr_(ptr) {}

  constexpr basic_mdspan( pointer ptr , const array<ptrdiff_t,extents_type::rank_dynamic()> dynamic_extents)
    : acc_(accessor_type()), map_( extents_type(dynamic_extents)), ptr_(ptr) {}

  constexpr basic_mdspan( pointer ptr , const mapping_type m ) noexcept
    : acc_(accessor_type()), map_( m ), ptr_(ptr) {}
  
  constexpr basic_mdspan( pointer ptr , const mapping_type m , const accessor_type a ) noexcept
    : acc_(a), map_( m ), ptr_(ptr) {}

  // [mdspan.basic.mapping]

  // Enforce rank() <= sizeof...(IndexType)
  template<class... IndexType >
  constexpr typename enable_if<sizeof...(IndexType)==extents_type::rank(),reference>::type
  operator()( IndexType... indices) const noexcept
    { return acc_.access( ptr_ , map_( indices... ) ); }

  // Enforce rank() == 1
  template<class IndexType>
  constexpr typename enable_if<is_integral<IndexType>::value && 1==extents_type::rank(),reference>::type
  operator[]( const IndexType i ) const noexcept
    { return acc_.access( ptr_ , map_(i) ); }

  // [mdspan.basic.domobs]

  static constexpr int rank() noexcept
    { return extents_type::rank(); }

  static constexpr int rank_dynamic() noexcept
    { return extents_type::rank_dynamic(); }

  constexpr index_type static_extent( size_t k ) const noexcept
    { return map_.extents().static_extent( k ); }

  constexpr index_type extent( int k ) const noexcept
    { return map_.extents().extent( k ); }

  constexpr const extents_type & extents() const noexcept
    { return map_.extents(); }

  // [mdspan.basic.codomain]

  // ------------------------------

//  constexpr fundamentals_v3::span<element_type> span() const noexcept
//    { return fundamentals_v3::span<element_type>(acc_.decay(ptr_),map_.required_span_size()); }

  // ------------------------------

  // [mdspan.basic.obs]

  static constexpr bool is_always_unique()     noexcept { return mapping_type::is_always_unique(); }
  static constexpr bool is_always_strided()    noexcept { return mapping_type::is_always_strided(); }
  static constexpr bool is_always_contiguous() noexcept { return mapping_type::is_always_contiguous(); }

  constexpr bool is_unique() const noexcept  { return map_.is_unique(); }
  constexpr bool is_strided() const noexcept { return map_.is_strided(); }
  constexpr bool is_contiguous() const noexcept {return map_.is_contiguous();}

  constexpr index_type stride( size_t r ) const noexcept
    { return map_.stride(r); }

  constexpr mapping_type mapping() const noexcept { return map_ ; }

  constexpr accessor_type accessor() const noexcept { return acc_ ; } 

  constexpr pointer data() const noexcept { return ptr_ ; } 
private:

  accessor_type acc_ ;
  mapping_type map_ ;
  pointer  ptr_ ;
};


template<class T, ptrdiff_t... Indices>
using mdspan = basic_mdspan<T,extents<Indices...>,layout_right,accessor_basic<T> > ;

}}} // experimental::fundamentals_v3



