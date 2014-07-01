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

// For modules that are only declared in the library API, \addtogroup
// cannot be used without an enclosing \cond.  Otherwise, it will create
// a dummy module with the identifier as the name...

//! \cond libapi
//! \addtogroup module_libmodule
//! \{

/*! \brief
 * Brief description.
 *
 * No \\libinternal is necessary here because of the enclosing \\cond.
 */
void gmx_function();

//! \}
//! \endcond

// An alternative to the above is use this, if the enclosing scope is only
// documented in the library API:

//! \libinternal \addtogroup module_libmodule
//! \{

//! Brief description.
void gmx_function()

//! \}

/*! \libinternal \brief
 * Brief description for a struct.
 *
 * Documented structs and classes from headers are always extracted into the
 * documentation, so \\libinternal is necessary to exclude it.
 * Currently, undocumented structs/classes do not produce warnings, so \\cond
 * is not necessary.
 */
struct t_example
{
    int  member1; //!< Each non-private member should be documented.
    bool member2; //!< Otherwise, Doxygen will produce warnings.
};

// This namespace is documented in the public API.
namespace gmx
{

//! \cond libapi
/*! \brief
 * Brief description for a free function within a documented namespace.
 *
 * In this case, the enclosing scope is the documented namespace,
 * so a \\cond is necessary to avoid warnings.
 */
void gmx_function();
//! \endcond

/*! \brief
 * Class meant for subclassing only within the module, but the subclasses will
 * be public.
 *
 * This base class still provides public methods that are visible through the
 * subclasses, so it should appear in the public documentation.
 * But it is not marked with \\inpublicapi.
 */
class BaseClass
{
    public:
        /*! \brief
         * A public method.
         *
         * This method also appears in the documentation of each subclass in
         * the public and library API docs.
         */
        void method();

    protected:
        // The \cond is necessary to exlude this documentation from the public
        // API, since the public API does not support subclassing.
        //! \cond internal
        //! A method that only subclasses inside the module see.
        void methodForSubclassToCall();

        //! A method that needs to be implemented by subclasses.
        virtual void virtualMethodToImplement() = 0;
        //! \endcond
};

} // namespace gmx
