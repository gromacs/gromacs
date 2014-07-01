// The struct itself is not documented; other comments within the declaration
// are ignored.

struct t_struct {

    // The comment tries to document both members at once, but it only
    // applies to the first.  The second produces warnings about missing
    // documentation (if the enclosing struct was documented).

    //! Angle parameters.
    double alpha, beta;
};


// This does not produce any brief documentation.
// An explicit \brief is required, or //! (C++) or /** */ (C) should be used.

/*! Brief comment. */
int gmx_function();


// This does not produce any documentation at all, since a ! is missing at
// the beginning.

/* \brief
 * Brief description.
 *
 * More details.
 */
int gmx_function();


// This puts the whole paragraph into the brief description.
// A short description is preferable, separated by an empty line from the rest
// of the text.

/*! \brief
 * Brief description. The description continues with all kinds of details about
 * what the function does and how it should be called.
 */
int gmx_function();


// This may be a Doxygen bug, but this does not produce any brief description.

/** \internal Brief description. */
int gmx_function();


// If the first declaration below appears in a header, and the second in a
// source file, then Doxygen does not associate them correctly and complains
// about missing documentation for the latter.  The solution is to explicitly
// add a namespace prefix also in the source file, even though the compiler
// does not require it.

// Header file
//! Example function with a namespace-qualified parameter type.
int gmx_function(const gmx::SomeClass &param);

// Source file
using gmx::SomeClass;

int gmx_function(const SomeClass &param);


// This puts the namespace into the mentioned module, instead of the contents
// of the namespace.  \addtogroup should go within the innermost scope.

//! \addtogroup module_example
//! \{

namespace gmx
{

//! Function intended to be part of module_example.
int gmx_function();

}
