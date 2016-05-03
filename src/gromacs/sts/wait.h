#ifndef STS_WAIT_H
#define STS_WAIT_H

/*
 * Functions for spin waiting
 *
 * This file provides a small set of functions intended to handle all
 * thread waiting so that waiting is easier to optimize.
 */

// TODO: Add gmx_pause.
// TODO: Consider wait_until_not variant that doesn't return a value.
/*! \brief
 * Wait until atomic variable a is set to value v
 *
 * \param[in] a   atomic variable
 * \param[in] v   value
 */
template <typename A, typename T>
void wait_until(const A &a, T v) {
    while (a.load() != v);
}
/*! \brief
 * Wait until atomic variable a is not set to value v
 *
 * \param[in] a   atomic variable
 * \param[in] v   value
 * \returns       new value of a
 */
template <typename A, typename T>
T wait_until_not(const A &a, T v) {
    T v2;
    do {
        v2 = a.load();
    } while(v2 == v);
    return v2;
}

#endif // STS_WAIT_H
