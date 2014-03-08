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
