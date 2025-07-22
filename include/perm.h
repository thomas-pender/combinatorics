/**
 * @file perm.h
 * @brief Public interface of permutation methods.
 *
 * Public interface of permutation methods. Permutations are assumed to act on the
 * right so that the composition ab of the permutations a and b is interpreted as
 * 'first apply, then apply b.'
 *
 * Permutations are assumed to act on the set {0, 1, ..., n - 1}. As such, they are
 * typically stored as arrays. The ith entry of the array is the image i^a. In
 * orther cases, it is better to store a permutation as an array of transversal
 * elements.
 *
 * @author Thomas Pender
 * @date 2025-07
 * @copyright GNU Public License
 */
# ifndef INCLUDED_PERM_H
# define INCLUDED_PERM_H

# include <stdio.h>
# include <stdint.h>
# include <stddef.h>
# include <stdlib.h>

/**
 * @brief Composition of two permutations.
 *
 * Composition of two permutations.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt> or <tt>b</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] a Permutation multiplicand.
 * @param[in] b Permutation multiplier.
 * @param[in] n The degree of the permutations.
 * @return The composition ab.
 */
static inline
uint32_t *perm_comp(const uint32_t a[static 1],
                    const uint32_t b[static 1], size_t n)
{
  uint32_t *c = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( size_t i = 0; i < n; i++ ) c[i] = b[a[i]];
  return c;
}

/**
 * @brief Composition of three permutations.
 *
 * Composition of three permutations. This is intended to be used mainly for
 * conjugation of b by c so that a = c^-1, but it will work for any three
 * permutations.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt>, <tt>b</tt> or <tt>c</tt> is not equal to
 * <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] a Permutation.
 * @param[in] b Permutation.
 * @param[in] c Permutation.
 * @param[in] n The degree of the permutations.
 * @return The composition abc.
 */
static inline
uint32_t *perm_ccomp(const uint32_t a[static 1],
                     const uint32_t b[static 1],
                     const uint32_t c[static 1], size_t n)
{
  uint32_t *p = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( size_t i = 0; i < n; i++ ) p[i] = c[b[a[i]]];
  return p;
}

/**
 * @brief Perform permutation composition in-place.
 *
 * Perform permutation composition in-place. The multiplicand is replaced by the
 * product.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt> or <tt>b</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in,out] a Permutation multiplicand.
 * @param[in] b Permutation multiplier.
 * @param[in] n The degree of the permutations.
 */
static inline
void perm_comp_inplace(uint32_t a[static 1],
                       const uint32_t b[static 1], size_t n)
{
  for ( size_t i = 0; i < n; i++ ) a[i] = b[a[i]];
}

/**
 * @brief Return inverse of permutation.
 *
 * Return inverse of permutation.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] a Permutation being inverted.
 * @param[in] n Degree of permutation.
 * @return Inverse of permutation <tt>a</tt>.
 */
static inline
uint32_t *perm_inv(const uint32_t a[restrict static 1], size_t n)
{
  uint32_t *c = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( size_t i = 0; i < n; i++ ) c[a[i]] = i;
  return c;
}

/**
 * @brief Return identity permutation.
 *
 * Return identity permutation. The identity permutation fixes every element of the
 * point set.
 *
 * @param[in] n Degree of permutation.
 * @return Identity permutation.
 */
static inline
uint32_t *perm_id(size_t n)
{
  uint32_t *p = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( size_t i = 0; i < n; i++ ) p[i] = i;
  return p;
}

/**
 * @brief Overwrite permutation with the identity permutation.
 *
 * Overwrite permutation with the identity permutation.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in,out] a Permutation being overwritten.
 * @param[in] n Degree of permutation.
 */
static inline
void perm_id_inplace(uint32_t a[static 1], size_t n)
{
  for ( size_t i = 0; i < n; i++ ) a[i] = i;
}

/**
 * @brief Check if permutation is the identity permutation.
 *
 * Check if permutation is the identity permutation.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>p</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] p Permutation being checked.
 * @param[in] n Degree of permuation.
 * @return 1 if <tt>p</tt> is the identity. -1 otherwise.
 */
static inline
int perm_id_cmp(const uint32_t p[static 1], size_t n)
{
  for ( size_t i = 0; i < n; i++ ) if ( p[i] != i ) return -1;
  return 1;
}

/**
 * @brief Return copy of permutation.
 *
 * Return copy of permutation.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>p</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] p Permutation being copied.
 * @param[in] n Degree of permutation.
 * @return Copy of permutation.
 */
static inline
uint32_t *perm_cpy(const uint32_t p[restrict static 1], size_t n)
{
  uint32_t *q = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( size_t i = 0; i < n; i++ ) q[i] = p[i];
  return q;
}

/**
 * @brief Copy permutation in place.
 *
 * Copy permutation in place. The contents of <tt>a</tt> are overwritten by those
 * of <tt>b</tt>.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt> or <tt>b</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] a Permutation being overwritten.
 * @param[in] b Permutation being copied.
 * @param[in] n Degree of permutations.
 */
static inline
void perm_cpy_inplace(uint32_t a[restrict static 1],
                      const uint32_t b[restrict static 1], size_t n)
{
  for ( size_t i = 0; i < n; i++ ) a[i] = b[i];
}

/**
 * @brief Print permutation in word form.
 *
 * Print permutation in word form.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of <tt>a</tt> is not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] a Permutation being printed.
 */
static inline
void perm_pr(const uint32_t a[static 1], size_t n)
{
  for ( size_t i = 0; i < n; i++ ) printf("%u ", a[i]);
}

/**
 * @brief Find the image of a point under the action of a permutation in array
 * form.
 *
 * Find the image of a point under the action of a permutation in array form.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd><tt>elem</tt> is larger than the degree of permutations.</dd>
 * <dd>The length of the array of permutations is not equal to <tt>l</tt>/</dd>
 * </dl>
 *
 * @param[in] elem Point being acted upon.
 * @param[in] u Array of permutations.
 * @param[in] l Length of array of permutations.
 * @return Image of <tt>elem</tt> under the action of the permutation array.
 */
static inline
uint32_t apply_elem(uint32_t elem, const uint32_t *u[static 1], size_t l)
{
  if ( l == 0 ) return elem;
  for ( long int i = l - 1; i >= 0; i-- ) elem = u[i][elem];
  return elem;
}

/**
 * @brief Find the image of a point under the action of the inverse of a
 * permutation in array form.
 *
 * Find the image of a point under the action of the inverse of a permutation in
 * array form.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd><tt>elem</tt> is larger than the degree of permutations.</dd>
 * <dd>The length of the array of permutations is not equal to <tt>l</tt>/</dd>
 * </dl>
 *
 * @param[in] elem Point being acted upon.
 * @param[in] u Array of permutations.
 * @param[in] l Length of array of permutations.
 * @return Image of <tt>elem</tt> under the action of the inverse of permutation
 * in array form.
 */
static inline
uint32_t apply_elem_inv(uint32_t elem, const uint32_t *u[static 1], size_t l)
{
  if ( l == 0 ) return elem;
  for ( long int i = 0; i < (long int)l; i++ ) elem = u[i][elem];
  return elem;
}

/**
 * @brief Convert permutation in array form to one in standard word form.
 *
 * Convert permutation in array form to one in standard word form.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of the array of permutations is not equal to <tt>l</tt>/</dd>
 * <dd>The degree of the permutations are not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] u Array of permutations.
 * @param[in] l Length of arrays of permutations.
 * @param[in] n Degree of permutations.
 * @return Permutation in standard word-form.
 */
static inline
uint32_t *arraytoelem(const uint32_t *restrict u[static 1], size_t l, size_t n)
{
  long int i, j;
  uint32_t *p = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( i = 0; i < (long int)n; i++ ) {
    p[i] = i;
    for ( j = l - 1; j >= 0; j-- ) p[i] = u[j][p[i]];
  }
  return p;
}

/**
 * @brief Convert permutation in array form to the inverse of the corresponding
 * permutation in standard word form.
 *
 * Convert permutation in array form to the inverse of the corresponding
 * permutation in standard word form.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of the array of permutations is not equal to <tt>l</tt>/</dd>
 * <dd>The degree of the permutations are not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] u Array of permutations.
 * @param[in] l Length of arrays of permutations.
 * @param[in] n Degree of permutations.
 * @return Permutation in standard word-form.
 */
static inline
uint32_t *arraytoeleminv(const uint32_t *restrict u[static 1], size_t l, size_t n)
{
  long int i, j;
  uint32_t *p = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( i = 0; i < (long int)n; i++ ) {
    p[i] = i;
    for ( j = 0; j < (long int)l; j++ ) p[i] = u[j][p[i]];
  }
  return p;
}

/**
 * @brief Print word-form of permutation in array form.
 *
 * Print word-form of permutation in array form.
 *
 * <dl>
 * <dt><strong>Unchecked Runtime Errors</strong></dt>
 * <dd>The length of the array of permutations is not equal to <tt>l</tt>/</dd>
 * <dd>The degree of the permutations are not equal to <tt>n</tt>.</dd>
 * </dl>
 *
 * @param[in] u Permutation in array form.
 * @param[in] l Length of array of permutations.
 * @param[in] n Degree of permutations.
 */
static inline
void array_print(const uint32_t *u[static 1], size_t l, size_t n)
{
  uint32_t *p = arraytoelem((const uint32_t *restrict*)u, l, n);
  for ( size_t i = 0; i < n; i++ ) printf("%u ", p[i]);
  free(p);
}

/**
 * @brief Lexicographic comparison of permutations.
 *
 * Lexicographic comparison of permutations.
 *
 * @param[in] _a First permutation.
 * @param[in] _b Second permutation.
 * @return -1 if _a < _b. 0 if _a = _b. 1 if _a > _b.
 */
static
int perm_cmp(const void *_a, const void *_b, void *_n)
{
  size_t n = *(size_t*)_n;
  uint32_t *a = (uint32_t*)_a;
  uint32_t *b = (uint32_t*)_b;
  for ( size_t i = 0; i < n; i++ ) {
    if ( a[i] < b[i] ) return -1;
    if ( a[i] > b[i] ) return 1;
  }
  return 0;
}

/**
 * @brief Reverse lexicographic comparison of permutations.
 *
 * Reverse lexicographic comparison of permutations.
 *
 * @param[in] _a First permutation.
 * @param[in] _b Second permutation.
 * @return -1 if _a < _b. 0 if _a = _b. 1 if _a > _b.
 */
static
int perm_cmprv(const void *_a, const void *_b, void *_n)
{
  size_t n = *(size_t*)_n;
  uint32_t *a = (uint32_t*)_a;
  uint32_t *b = (uint32_t*)_b;
  for ( size_t i = 0; i < n; i++ ) {
    if ( a[i] < b[i] ) return 1;
    if ( a[i] > b[i] ) return -1;
  }
  return 0;
}

# endif
