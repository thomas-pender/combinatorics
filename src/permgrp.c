/**
 * @file permgrp.c
 * @brief Main implementation file for <tt>permgrp_t</tt> class. NB: Some functions
 * are defined in their own files when their definitions are suitably complex.
 *
 * Main implementation file for <tt>permgrp_t</tt> class. NB: Some functions
 * are defined in their own files when their definitions are suitably complex.
 *
 * @author Thomas Pender
 * @date 2025-07
 */
# include <permgrp.h>
# include <permgrp_struct_p.h>
# include <perm.h>

# include <assert.h>

extern void swap(void*, void*, size_t);

size_t permgrp_degree(permgrp_t G)
{
  return G->degree;
}

size_t permgrp_baselen(permgrp_t G)
{
  return G->base_size;
}

void permgrp_free(permgrp_t *G)
{
  if ( G == NULL || *G == NULL ) return;
  size_t i, j;
  for ( i = 0; i < (*G)->degree; i++ )
    if ( (*G)->trees[i] != NULL ) {
      for ( j = 0; j < (*G)->degree; j++ )
        if ( (*G)->trees[i]->T[j] != NULL ) {
          free((*G)->trees[i]->T[j]);
          free((*G)->trees[i]->T_inv[j]);
        }
      free((*G)->trees[i]->T);
      free((*G)->trees[i]->T_inv);
      hashtabs_free(&(*G)->trees[i]->S);
      free((*G)->trees[i]);
    }
  free((*G)->trees); free((*G)->base); free((*G)->order); free((*G)->id);
  arrays_free(&(*G)->gens);
  free(*G);
  *G = NULL;
}

size_t permgrp_order(permgrp_t G)
{
  if ( G == NULL ) return 1;
  size_t num = 1;
  for ( size_t i = 0; i < G->base_size; i++ ) num *= G->trees[i]->size;
  return num;
}

size_t permgrp_sbgrporder(permgrp_t G, size_t i)
{
  size_t num = 1;
  for ( size_t j = i; j < G->base_size; j++ ) num *= G->trees[j]->size;
  return num;
}

void permgrp_new_base(permgrp_t G[static 1],
                      const uint32_t *base, size_t b)
{
  uint32_t **gens = (uint32_t**)malloc(arrays_nmem((*G)->gens) *
                                       sizeof(uint32_t*));
  for ( size_t i = 0; i < arrays_nmem((*G)->gens); i++ )
    gens[i] = arrays_at((*G)->gens, i);

  permgrp_t H = permgrp_new((const uint32_t**)gens,
                            arrays_nmem((*G)->gens), (*G)->degree, base, b);
  permgrp_swap(G, &H);
  permgrp_free(&H);

  free(gens);
}

permgrp_t permgrp_sym(size_t n)
{
  if ( n <= 1 ) return NULL;
  if ( n == 2 ) {
    uint32_t *gen = (uint32_t*)malloc(n * sizeof(uint32_t));
    gen[0] = 1; gen[1] = 0;
    permgrp_t Sym = permgrp_new((const uint32_t *restrict*)&gen, 1, n, NULL, 0);
    free(gen);
    return Sym;
  }

  uint32_t **gens = (uint32_t**)malloc(2 * sizeof(uint32_t*));
  gens[0] = (uint32_t*)malloc(n * sizeof(uint32_t));
  gens[1] = (uint32_t*)malloc(n * sizeof(uint32_t));

  size_t i;
  for ( i = 0; i < n; i++ ) {
    gens[0][i] = i;
    gens[1][i] = (i + 1) % n;
  }
  swap(&gens[0][0], &gens[0][1], sizeof(gens[0][0]));

  uint32_t *base = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( i = 0; i < n; i++ ) base[i] = i;

  permgrp_t Sym = permgrp_new((const uint32_t *restrict*)gens, 2, n, base, n);

  free(base);
  free(gens[0]); free(gens[1]); free(gens);

  return Sym;
}

uint32_t *permgrp_fundorb(permgrp_t G, size_t i)
{
  if ( G == NULL || i > G->base_size ) return NULL;
  uint32_t *orb = (uint32_t*)malloc(G->trees[i]->size * sizeof(uint32_t));
  for ( uint32_t j = 0, index = 0;
        j < G->degree && index < G->trees[i]->size; j++ )
    if ( G->trees[i]->T[j] != NULL ) orb[index++] = j;
  return orb;
}

size_t permgrp_fundorb_size(permgrp_t G, size_t i)
{
  return G->trees[i]->size;
}

uint32_t *permgrp_ordering(permgrp_t G)
{
  return G->order;
}

uint32_t *permgrp_rep(permgrp_t G, size_t i, uint32_t gamma)
{
  return G->trees[i]->T[gamma];
}

uint32_t *permgrp_repinv(permgrp_t G, size_t i, uint32_t gamma)
{
  return G->trees[i]->T_inv[gamma];
}


uint32_t *permgrp_base(permgrp_t G)
{
  return G->base;
}


size_t permgrp_fixed_base(permgrp_t G, const uint32_t p[static 1])
{
  for ( size_t i = 0; i < G->base_size; i++ )
    if ( G->base[i] != p[G->base[i]] ) return i;
  return G->base_size;
}

static
int gens_apply(void **_x, void *_gens)
{
  arrays_t gens = (arrays_t)_gens;
  uint32_t *x = *(uint32_t**)_x;
  arrays_push(gens, x);
  return 1;
}

arrays_t permgrp_generators(permgrp_t G)
{
  if ( G == NULL ) return NULL;
  arrays_t gens = arrays_new(arrays_size(G->gens), hashtabs_size(G->trees[0]->S));
  hashtabs_map_r(G->trees[0]->S, gens_apply, gens);
  return gens;
}

arrays_t permgrp_generatorssbgrp(permgrp_t G, size_t i)
{
  if ( G == NULL ) return NULL;
  arrays_t gens = arrays_new(arrays_size(G->gens), hashtabs_size(G->trees[i]->S));
  hashtabs_map_r(G->trees[i]->S, gens_apply, gens);
  return gens;
}

size_t permgrp_effective_len(permgrp_t G)
{
  size_t i = G->base_size;
  while ( i > 0 && G->trees[i - 1]->size <= 1 ) i--;
  return i;
}

arrays_t permgrp_sgs(permgrp_t G)
{
  return G->gens;
}

/* ---------------------------------------------------------------------------
 * Routines to update fundamental orbits. In many instances, one needs the
 * partitions of the point set under the actions of the groups in the point
 * stabilizer chain of the permutation group.
 * ---------------------------------------------------------------------------*/

typedef struct {
  uint32_t base, *a, *points, *min;
  sstacks_t *tmp;
  permgrp_t G;
} fundorb_args_t;

static
int hash_fundorb(void **_x, void *_args)
{
  uint32_t *x = *(uint32_t**)_x;
  fundorb_args_t *args = (fundorb_args_t*)_args;

  permgrp_t G = args->G;
  uint32_t *a = args->a, *points = args->points, *min = args->min;
  sstacks_t *tmp = args->tmp;

  uint32_t elem = x[args->base];
  if ( a[elem] == 0 ) {
    points[elem] = 0;
    a[elem] = 1;
    if ( G->order[elem] < G->order[*min] ) *min = elem;
    sstacks_push(tmp, &elem);
  }

  return 1;
}

static inline
uint32_t orb_next_point(const uint32_t points[static 1], size_t n)
{
  for ( uint32_t i = 0; i < n; i++ )
    if ( points[i] == 1 ) return i;
  return n;
}

void permgrp_orbrep_update(permgrp_t G, uint32_t orbrep[static 1], size_t i)
{
  size_t j;
  if ( G == NULL || G->trees[i]->S == NULL ) {
    for ( j = 0; j < G->degree; j++ ) orbrep[j] = 1;
    return;
  }

  uint32_t points[G->degree];
  for ( j = 0; j < G->degree; j++ ) {
    points[j] = 1;
    orbrep[j] = 0;
  }

  uint32_t base;
  while ( (base = orb_next_point(points, G->degree)) != G->degree ) {
    uint32_t min = base;
    sstacks_t *tmp = NULL;
    sstacks_new(tmp, sizeof(min), G->degree);
    sstacks_push(tmp, &base);

    uint32_t a[G->degree];
    for ( j = 0; j < G->degree; j++ ) a[j] = 0;
    a[base] = 1;
    points[base] = 0;

    while ( !sstacks_empty(tmp) ) {
      fundorb_args_t *args = (fundorb_args_t*)malloc(sizeof(fundorb_args_t));

      args->base = *sstacks_pop(tmp);
      args->a = &a[0];
      args->points = &points[0];
      args->min = &min;
      args->tmp = tmp;
      args->G = G;

      hashtabs_map_r(G->trees[i]->S, hash_fundorb, args);

      free(args);
    }

    sstacks_free(tmp);
    orbrep[min] = 1;
  }
}

