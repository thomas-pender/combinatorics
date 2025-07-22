/**
 * @file diophantine.h
 * @brief Public interface for library routines used to find solutions to systems
 * of linear diophantine equations.
 *
 * Public interface for library routines used to find solutions to systems of
 * linear diophantine equations.
 *
 * @author Thomas Pender
 * @date 2025-07
 */
# ifndef INCLUDED_DIOPHANTINE_H
# define INCLUDED_DIOPHANTINE_H

# include <stddef.h>
# include <stdlib.h>
# include <stdint.h>

# include <containers/arrays.h>

/**
 * @brief Return all solutions to integer program.
 *
 * All solutions x to the integer program determined by Ax = b subject to
 * 0 <= l <= x <= u are calculated. The solutions are stored in an
 * <tt>arrays_t</tt> object which is returned by the function. The user may
 * optionally pass a function to apply to every encountered solution x. This
 * function returns an integer value. If the value returned is negative, the
 * solution x is not stored in the <tt>arrays_t</tt> object. If the value is
 * nonnegative, the solution x is stored in the <tt>arrays_t</tt> object.
 *
 * @param[in] A Coefficient matrix.
 * @param[in] b Constraint vector.
 * @param[in] l Nonegative lower bound for solution vector.
 * @param[in] u Nonegative upper bound for solution vector.
 * @param[in] nrows Number identities in the program, i.e., the number of rows of
 * A.
 * @param[in] ncols Number of decision variables, i.e., the number of columns of
 * A.
 * @param[in] apply User defined function to apply to each encountered solution
 * that returns an integer. If the returned value is negative, the solution is
 * appended to the solutions array; otherwise, the solution passed by.
 *
 * @return An <tt>arrays_t</tt> object containing the admissible solutions to the
 * system.
 */
extern arrays_t diophantine(const int32_t *restrict A[static 1],
                            const int32_t b[restrict static 1],
                            const int32_t l[restrict static 1],
                            const int32_t u[restrict static 1],
                            size_t nrows, size_t ncols, int apply(int32_t*));

/**
 * @brief Return all solutions to integer program.
 *
 * Reentrant version of <tt>diophantine</tt>.
 *
 * @param[in] A Coefficient matrix.
 * @param[in] b Constraint vector.
 * @param[in] l Nonegative lower bound for solution vector.
 * @param[in] u Nonegative upper bound for solution vector.
 * @param[in] nrows Number identities in the program, i.e., the number of rows of
 * A.
 * @param[in] ncols Number of decision variables, i.e., the number of columns of
 * A.
 * @param[in] apply User defined function to apply to each encountered solution
 * that returns an integer. If the returned value is negative, the solution is
 * appended to the solutions array; otherwise, the solution passed by.
 * @param[in] y Argument to be passed to <tt>apply</tt> parameter.
 *
 * @return An <tt>arrays_t</tt> object containing the admissible solutions to the
 * system.
 */
extern arrays_t diophantine_r(const int32_t *restrict A[static 1],
                              const int32_t b[restrict static 1],
                              const int32_t l[restrict static 1],
                              const int32_t u[restrict static 1],
                              size_t nrows, size_t ncols,
                              int apply(int32_t*, void*), void *y);

# endif
