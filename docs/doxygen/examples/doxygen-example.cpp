/*! \libinternal \file
 * \brief
 * Declares gmx::MyClass.
 *
 * More details.  The documentation is still extracted for the class even if
 * this whole comment block is missing.
 *
 * \author Example Author <example@author.com>
 * \inlibraryapi
 * \ingroup module_mymodule
 */

namespace gmx
{

/*! \libinternal
 * \brief
 * Brief description for the class.
 *
 * More details.  The \\libinternal tag is required for classes, since they are
 * extracted into the documentation even in the absence of documentation for
 * the enclosing scope.
 * The \\libinternal tag is on a separate line because of a bug in Doxygen
 * 1.8.5 (only affects \\internal, but for clarity it is also worked around
 * here).
 *
 * \inlibraryapi
 * \ingroup module_mymodule
 */
class MyClass
{
    public:
        // Trivial constructors or destructors do not require documentation.
        // But if a constructor takes parameters, it should be documented like
        // methods below.
        MyClass();
        ~MyClass();

        /*! \brief
         * Brief description for the method.
         *
         * \param[in] param1  Description of the first parameter.
         * \param[in] param2  Description of the second parameter.
         * \returns   Description of the return value.
         * \throws    std::bad_alloc if out of memory.
         *
         * More details describing the method.  It is not an error to put this
         * above the parameter block, but most existing code has it here.
         */
        int myMethod(int param1, const char *param2) const;

        //! Brief description for the accessor.
        int simpleAccessor() const { return var_; }
        /*! \brief
         * Alternative, more verbose way of specifying a brief description.
         */
        int anotherAccessor() const;
        /*! \brief
         * Brief description for another accessor that is so long that it does
         * not conveniently fit on a single line cannot be specified with //!.
         */
        int secondAccessor() const;

    private:
        // Private members (whether methods or variables) are currently ignored
        // by Doxygen, so they don't need to be documented.  Documentation
        // doesn't hurt, though.
        int var_;
};

} // namespace gmx
