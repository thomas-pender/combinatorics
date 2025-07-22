/**
 * @file diophantine.c
 * @brief Implementation of diophantine equation solver.
 *
 * Implementation of diophantine equation solver.
 *
 * @author Thomas Pender
 * @date 2025-07
 */
# include <diophantine.h>

# define INIT_NUM 100

static
void diophantine_in(size_t j,
                    const int32_t *restrict A[static 1],
                    const int32_t b[restrict static 1],
                    const int32_t l[restrict static 1],
                    const int32_t u[restrict static 1],
                    size_t nrows, size_t ncols,
                    int apply(int32_t*),
                    const int32_t m[restrict static 1],
                    int32_t s[restrict static 1],
                    int32_t x[restrict static 1],
                    arrays_t sols)
{
  if ( j == ncols ) {
    if ( apply != NULL && apply(x) < 0 ) return;
    arrays_dynpush(sols, x);
    return;
  }

  int flag;
  int32_t _x, i;

  for ( _x = l[j]; _x <= u[j]; _x++ ) {
    x[j] = _x;
    flag = 1;
    int32_t *s2 = (int32_t*)malloc(nrows * sizeof(int32_t));

    for ( i = 0; i < (int32_t)nrows; i++ ) {
      s2[i] = s[i] + (A[i][j] * x[j]);
      if ( s2[i] > b[i] ) {
        free(s2);
        return;
      }
      if ( (int32_t)j == m[i] && s2[i] != b[i] ) {
        flag = -1;
        break;
      }
    }

    if ( flag == 1 )
      diophantine_in(j + 1, A, b, l, u, nrows, ncols, apply, m, s2, x, sols);

    free(s2);
  }
}

arrays_t diophantine(const int32_t *restrict A[static 1],
                     const int32_t b[restrict static 1],
                     const int32_t l[restrict static 1],
                     const int32_t u[restrict static 1],
                     size_t nrows, size_t ncols, int apply(int32_t*))
{
  int32_t i, j;
  int32_t *m = (int32_t*)calloc(nrows, sizeof(int32_t));
  int32_t *s = (int32_t*)calloc(nrows, sizeof(int32_t));
  int32_t *x = (int32_t*)calloc(ncols, sizeof(int32_t));
  arrays_t sols = arrays_new(ncols * sizeof(l[0]), INIT_NUM);

  for ( i = 0; i < (int32_t)nrows; i++ ) {
    for ( j = (int32_t)ncols - 1; j >= 0; j-- )
      if ( A[i][j] != 0 ) {
        m[i] = j;
        break;
      }
  }

  diophantine_in(0, A, b, l, u, nrows, ncols, apply, m, s, x, sols);

  free(x); free(s); free(m);

  return sols;
}

static
void diophantine_in_r(size_t j,
                      const int32_t *restrict A[static 1],
                      const int32_t b[restrict static 1],
                      const int32_t l[restrict static 1],
                      const int32_t u[restrict static 1],
                      size_t nrows, size_t ncols,
                      int apply(int32_t*, void*), void *y,
                      const int32_t m[restrict static 1],
                      int32_t s[restrict static 1],
                      int32_t x[restrict static 1],
                      arrays_t sols)
{
  if ( j == ncols ) {
    if ( apply != NULL && apply(x, y) < 0 ) return;
    arrays_dynpush(sols, x);
    return;
  }

  int flag;
  int32_t _x, i;

  for ( _x = l[j]; _x <= u[j]; _x++ ) {
    x[j] = _x;
    flag = 1;
    int32_t *s2 = (int32_t*)malloc(nrows * sizeof(int32_t));

    for ( i = 0; i < (int32_t)nrows; i++ ) {
      s2[i] = s[i] + (A[i][j] * x[j]);
      if ( s2[i] > b[i] ) {
        free(s2);
        return;
      }
      if ( (int32_t)j == m[i] && s2[i] != b[i] ) {
        flag = -1;
        break;
      }
    }

    if ( flag == 1 )
      diophantine_in_r(j + 1, A, b, l, u, nrows, ncols, apply, y, m, s2, x, sols);

    free(s2);
  }
}

arrays_t diophantine_r(const int32_t *restrict A[static 1],
                       const int32_t b[restrict static 1],
                       const int32_t l[restrict static 1],
                       const int32_t u[restrict static 1],
                       size_t nrows, size_t ncols,
                       int apply(int32_t*, void*), void *y)
{
  int32_t i, j;
  int32_t *m = (int32_t*)calloc(nrows, sizeof(int32_t));
  int32_t *s = (int32_t*)calloc(nrows, sizeof(int32_t));
  int32_t *x = (int32_t*)calloc(ncols, sizeof(int32_t));
  arrays_t sols = arrays_new(ncols * sizeof(l[0]), INIT_NUM);

  for ( i = 0; i < (int32_t)nrows; i++ ) {
    for ( j = (int32_t)ncols - 1; j >= 0; j-- )
      if ( A[i][j] != 0 ) {
        m[i] = j;
        break;
      }
  }

  diophantine_in_r(0, A, b, l, u, nrows, ncols, apply, y, m, s, x, sols);

  free(x); free(s); free(m);

  return sols;
}
