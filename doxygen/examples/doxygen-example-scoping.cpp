/*! \libinternal \file
 * \brief
 * ...
 *
 * The examples below assume that the file is documented like this:
 * with an \\libinternal definition at the beginning, with an intent to not
 * expose anything from the file in the public API.  Things work similarly for
 * the full documentation if you replace \\libinternal with \\internal
 * everywhere in the example.
 *
 * \ingroup module_example
 */


/*! \brief
 * Brief description for a free function.
 *
 * A free function is not extracted into the documentation unless the enclosing
 * scope (in this case, the file) is.  So a \\libinternal is not necessary.
 */
void gmx_function();

// Assume that the module_example group is defined in the public API.

//! \addtogroup module_example
//! \{

//! \cond libapi
/*! \brief
 * Brief description for a free function within \\addtogroup.
 *
 * In this case, the enclosing scope is actually the module_example group,
 * which is documented, so the function needs to be explicitly excluded.
 * \\libinternal does not work, since it would produce warnings about an
 * undocumented function, so the whole declaration is hidden from Doxygen.
 */
void gmx_function();
//! \endcond

//! \}

/*! \libinternal \brief
 * Brief description for a struct.
 *
 * Documented structs and classes are always extracted into the documentation,
 * so \\libinternal is necessary to exclude it.  Currently, undocumented
 * structs/classes do not produce warnings, so \\cond is not necessary.
 */
struct t_example
{
    int  member1; //!< Each non-private member should be documented.
    bool member2; //!< Otherwise, Doxygen will produce warnings.
};

// This namespace is documented in the public API.
namespace gmx
{

/*! \libinternal \brief
 * Brief description for a free function within a documented namespace.
 *
 * In this case, the enclosing scope is the documented namespace,
 * so a \\libinternal is necessary.
 */
void gmx_function();

} // namespace gmx
