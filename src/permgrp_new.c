/**
 * @file permgrp_new.c
 * @brief Implementatino of <tt>permgrp_new</tt> of <tt>permgrp_t</tt> class.
 *
 * Implementatino of <tt>permgrp_new</tt> of <tt>permgrp_t</tt> class.
 *
 * @author Thomas Pender
 * @date 2025-07
 */
# include <config.h>

# include <permgrp.h>
# include <permgrp_struct_p.h>
# include <perm.h>

extern void orbit(tree_t*, uint32_t, size_t);
extern uint32_t *sift(permgrp_t, const uint32_t*restrict, size_t*);
extern void hash_insert(hashtabs_t*, const uint32_t*, size_t);
extern void swap(void*, void*, size_t);

# ifdef __DEBUG__
extern void print_grp_info(permgrp_t);
# endif

/* ---------------------------------------------------------------------------
 * auxiliary routines for initialization
 * ---------------------------------------------------------------------------*/

static inline
int fixed(const uint32_t p[restrict static 1],
          const uint32_t base[restrict static 1], size_t b)
{
  for ( size_t i = 0; i < b; i++ ) if ( p[base[i]] != base[i] ) return -1;
  return 1;
}

static inline
uint32_t moved(const uint32_t p[restrict static 1],
               const uint32_t base[restrict static 1],
               size_t n, size_t b)
{
  int flag;
  size_t i, j;
  for ( i = 0; i < n; i++ ) {
    flag = 1;
    for ( j = 0; j < b; j++ )
      if ( base[j] == i ) {
        flag = -1;
        break;
      }
    if ( flag == 1 && p[i] != i ) return i;
  }
  error(1, errno, "ERROR -- identity generator");
  return n;
}

static inline
size_t moved_index(permgrp_t G, const uint32_t p[static 1])
{
  for ( size_t i = G->base_size; i < G->degree; i++ )
    if ( p[G->base[i]] != G->base[i] ) return i;
  return G->degree;
}

static inline
void ordering(permgrp_t G)
{
  for ( size_t i = 0; i < G->degree; i++ ) G->order[G->base[i]] = i;
}

/* ---------------------------------------------------------------------------
 * Schreier-Sims algorithm
 * ---------------------------------------------------------------------------*/

typedef struct {
  permgrp_t G;
  size_t a;
  int *flag, *i;
} schreier_sims_t;

static inline
tree_t *tree_new(size_t n)
{
  tree_t *T = (tree_t*)malloc(sizeof(tree_t));
  T->T = (uint32_t**)calloc(n, sizeof(uint32_t*));
  T->T_inv = (uint32_t**)calloc(n, sizeof(uint32_t*));
  T->S = hashtabs_new(NULL, perm_cmp, NULL, hashtabs_stdhash, 0);
  T->size = 0;
  return T;
}

static
void schreier_sims_init(permgrp_t G)
{
  size_t i, j, ngens = arrays_nmem(G->gens);
  G->trees[0] = tree_new(G->degree);
  for ( i = 0; i < ngens; i++ )
    hash_insert(&G->trees[0]->S, arrays_at(G->gens, i), G->degree);
  orbit(G->trees[0], G->base[0], G->degree);

  for ( i = 1; i < G->base_size; i++ ) {
    G->trees[i] = tree_new(G->degree);
    for ( j = 0; j < ngens; j++ )
      if ( fixed(arrays_at(G->gens, j), G->base, i) > 0 )
        hash_insert(&G->trees[i]->S, arrays_at(G->gens, j), G->degree);

    if ( hashtabs_size(G->trees[i]->S) == 0 )
      hash_insert(&G->trees[i]->S, G->id, G->degree);

    orbit(G->trees[i], G->base[i], G->degree);
  }
}

static
int _schreier_sims_in(void **_x, void *_args)
{
  uint32_t *x = *(uint32_t**)_x;
  schreier_sims_t *args = (schreier_sims_t*)_args;

  int *flag = args->flag, *i = args->i;
  size_t a = args->a, j, k;
  permgrp_t G = args->G;

  uint32_t *ux = perm_comp(G->trees[*i]->T[a], x, G->degree);
  uint32_t elem = ux[G->base[*i]];
  uint32_t *g = perm_comp(ux, G->trees[*i]->T_inv[elem], G->degree);
  if ( perm_id_cmp(g, G->degree) < 0 ) {
    int y = 1;
    uint32_t *s = sift(G, g, &j);
    if ( j < G->base_size ) y = -1;
    else if ( perm_id_cmp(s, G->degree) < 0 ) {
      y = -1;
      size_t new_base_index = moved_index(G, s);
      swap(&G->base[G->base_size], &G->base[new_base_index], sizeof(G->base[0]));
      G->trees[G->base_size] = tree_new(G->degree);
      G->base_size++;
    }
    if ( y < 0 ) {
      if ( arrays_push(G->gens, s) < 0 )
        error(1, errno, "%s:%d array at capacity", __func__, __LINE__);

      for ( k = *i + 1; k <= j; k++ ) {
        hash_insert(&G->trees[k]->S,
                    arrays_at(G->gens, arrays_nmem(G->gens) - 1), G->degree);
        orbit(G->trees[k], G->base[k], G->degree);
      }
      *i = j;
      *flag = -1;
    }
    free(s);
  }
  free(g); free(ux);

  return *flag;
}

static
void schreier_sims(permgrp_t G)
{
  size_t a;
  int i = G->base_size - 1, flag;
  while ( i >= 0 ) {
    flag = 1;
    for ( a = 0; a < G->degree && flag == 1; a++ )
      if ( G->trees[i]->T[a] != NULL ) {
        size_t _i = i;
        schreier_sims_t *args = (schreier_sims_t*)malloc(sizeof(schreier_sims_t));
        args->G = G;
        args->a = a;
        args->i = &i;
        args->flag = &flag;
        hashtabs_map_r(G->trees[_i]->S, _schreier_sims_in, args);
        free(args);
      }
    if ( flag == 1 ) i--;
  }
}

/* ---------------------------------------------------------------------------
 * create new permgrp_t instance
 * ---------------------------------------------------------------------------*/

static
void permgrp_init(permgrp_t G, const uint32_t *restrict gens[static 1],
                  size_t ngens, size_t n, const uint32_t *restrict base, size_t b)
{
  G->base = (uint32_t*)malloc(n * sizeof(uint32_t));
  G->order = (uint32_t*)malloc(n * sizeof(uint32_t));
  G->base_size = b;
  G->degree = n;
  G->gens = arrays_new(n * sizeof(gens[0][0]), 100);
  G->id = perm_id(n);
  G->trees = (tree_t**)calloc(n, sizeof(tree_t*));

  /* deep copy generators */
  size_t i;
  for ( i = 0; i < ngens; i++ )
    if ( arrays_push(G->gens, gens[i]) < 0 )
      error(1, errno, "%s:%d generator array at capacity", __func__, __LINE__);

  /* initialize base prefix */
  if ( base != NULL && b > 0 ) for ( i = 0; i < b; i++ ) G->base[i] = base[i];

  /* extend base prefix to a proper initial base if necessary */
  for ( i = 0; i < ngens; i++ ) {
    if ( fixed(gens[i], G->base, G->base_size) > 0 ) {
      G->base[G->base_size] = moved(gens[i], G->base, n, G->base_size);
      G->base_size++;
    }
  }

  /* place remaining points for later use */
  int flag;
  size_t j, index = G->base_size;
  for ( i = 0; i < n; i++ ) {
    flag = 1;
    for ( j = 0; j < index; j++ )
      if ( G->base[j] == i ) {
        flag = -1;
        break;
      }
    if ( flag == 1 ) G->base[index++] = i;
  }
}

permgrp_t permgrp_new(const uint32_t *restrict gens[static 1], size_t ngens,
                      size_t n, const uint32_t *restrict base, size_t b)
{
  permgrp_t G;
  G = (permgrp_t)malloc(sizeof(*G));
  permgrp_init(G, gens, ngens, n, base, b);
  schreier_sims_init(G);
  schreier_sims(G);

  for ( size_t i = G->base_size; i < G->degree; i++ ) {
    if ( G->trees[i] == NULL ) {
      G->trees[i] = tree_new(G->degree);
      G->trees[i]->size = 1;
      G->trees[i]->T[G->base[i]] = perm_id(G->degree);
      G->trees[i]->T_inv[G->base[i]] = perm_id(G->degree);
      hash_insert(&G->trees[i]->S, G->id, G->degree);
    }
    else if ( hashtabs_size(G->trees[i]->S) == 0 ) {
      G->trees[i]->size = 1;
      hash_insert(&G->trees[i]->S, G->id, G->degree);
    }
    G->base_size++;
  }

  ordering(G);

  /* print_grp_info(G); */

  return G;
}

/* ---------------------------------------------------------------------------
 * update permgrp_t with new generator
 * ---------------------------------------------------------------------------*/

static
int permgrp_update_in(permgrp_t G, const uint32_t g[static 1])
{
  size_t j;
  int y = 1;
  uint32_t *s = sift(G, g, &j);
  if ( j < G->base_size ) y = -1;
  else if ( perm_id_cmp(s, G->degree) < 0 ) {
    y = -1;
    size_t new_base_index = moved_index(G, s);
    swap(&G->base[G->base_size], &G->base[new_base_index], sizeof(G->base[0]));
    G->trees[G->base_size] = tree_new(G->degree);
    G->base_size++;
  }

  if ( y >= 0 ) {
    free(s);
    return -1;
  }

  if ( arrays_push(G->gens, s) < 0 )
    error(1, errno, "%s:%d array at capacity", __func__, __LINE__);

  for ( size_t k = 0; k <= j; k++ ) {
    hash_insert(&G->trees[k]->S,
                arrays_at(G->gens, arrays_nmem(G->gens) - 1), G->degree);
    orbit(G->trees[k], G->base[k], G->degree);
  }
  free(s);

  return (int)j;
}

void permgrp_update(permgrp_t G, const uint32_t g[static 1])
{
  int i;
  if ( (i = permgrp_update_in(G, g)) < 0 ) return;

  size_t a;
  int flag;
  while ( i >= 0 ) {
    flag = 1;
    for ( a = 0; a < G->degree && flag == 1; a++ )
      if ( G->trees[i]->T[a] != NULL ) {
        size_t _i = i;
        schreier_sims_t *args = (schreier_sims_t*)malloc(sizeof(schreier_sims_t));
        args->G = G;
        args->a = a;
        args->i = &i;
        args->flag = &flag;
        hashtabs_map_r(G->trees[_i]->S, _schreier_sims_in, args);
        free(args);
      }
    if ( flag == 1 ) i--;
  }

  ordering(G);
}
