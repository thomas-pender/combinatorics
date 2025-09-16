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
# include <config.h>

# include <permgrp.h>
# include <permgrp_struct_p.h>
# include <perm.h>

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

size_t permgrp_orderfast(permgrp_t G)
{
  if ( G == NULL ) return 1;
  size_t num = 1;
  for ( size_t i = 0; i < G->base_size; i++ ) num *= G->trees[i]->size;
  return num;
}

size_t permgrp_sbgrporderfast(permgrp_t G, size_t i)
{
  size_t num = 1;
  for ( size_t j = i; j < G->base_size; j++ ) num *= G->trees[j]->size;
  return num;
}

# ifdef HAVE_GMP_H
void permgrp_order(mpz_t n, permgrp_t G)
{
  mpz_set_ui(n, 1UL);
  if ( G == NULL ) return;
  for ( size_t i = 0; i < G->base_size; i++ )
    mpz_mul_ui(n, n, G->trees[i]->size);
}

void permgrp_sbgrporder(mpz_t n, permgrp_t G, size_t i)
{
  mpz_set_ui(n, 1UL);
  if ( G == NULL ) return;
  for ( size_t j = i; j < G->base_size; j++ )
    mpz_mul_ui(n, n, G->trees[j]->size);
}
# endif

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

  size_t i;
  uint32_t *base = (uint32_t*)malloc(n * sizeof(uint32_t));
  for ( i = 0; i < n; i++ ) base[i] = i;

  if ( n == 2 ) {
    uint32_t *gen = (uint32_t*)malloc(n * sizeof(uint32_t));
    gen[0] = 1; gen[1] = 0;
    permgrp_t Sym = permgrp_new((const uint32_t *restrict*)&gen, 1, n, base, n);
    free(gen);
    free(base);
    return Sym;
  }

  uint32_t **gens = (uint32_t**)malloc(2 * sizeof(uint32_t*));
  gens[0] = (uint32_t*)malloc(n * sizeof(uint32_t));
  gens[1] = (uint32_t*)malloc(n * sizeof(uint32_t));

  for ( i = 0; i < n; i++ ) {
    gens[0][i] = i;
    gens[1][i] = (i + 1) % n;
  }
  swap(&gens[0][0], &gens[0][1], sizeof(gens[0][0]));

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

      args->base = *(uint32_t*)sstacks_pop(tmp);
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

/* ---------------------------------------------------------------------------
 * Setwize stabilizer routines.
 * ---------------------------------------------------------------------------*/

static inline
void _update_fund_orbs(permgrp_t K, size_t j, sets_t **fund_orbs, size_t m)
{
  size_t a, b;
  for ( a = 0; a < j; a++ ) {
    _EMPTYSET(fund_orbs[a], m);
    for ( b = 0; b < permgrp_degree(K); b++ )
      if ( permgrp_rep(K, a, b) != NULL )
        _ADDELEMENT(fund_orbs[a], b);
  }
}

static inline
uint32_t _orb_next_point(const uint32_t points[static 1], size_t deg)
{
  for ( uint32_t i = 0; i < deg; i++ ) if ( points[i] == 1 ) return i;
  return deg;
}

static inline
int _test(permgrp_t G, const sets_t *fund_orbs[static 1],
          const uint32_t *u[static 1], uint32_t gamma, size_t l)
{
  uint32_t elem;
  uint32_t *base = permgrp_base(G), *order = permgrp_ordering(G);
  for ( size_t i = 0; i < l; i++ ) {
    if ( _ISELEMENT(fund_orbs[i], base[l]) ) {
      elem = apply_elem(base[i], u, i + 1);
      if ( order[gamma] <= order[elem] ) return -1;
    }
  }
  return 1;
}

static inline
size_t _fixed_base(permgrp_t G, const uint32_t p[static 1], size_t j)
{
  uint32_t *base = permgrp_base(G);
  for ( size_t i = 0; i < j; i++ )
    if ( base[i] != p[base[i]] ) return i;
  return j;
}

static inline
void _orbit_apply(uint32_t orbit[static 1],
                  const uint32_t *restrict u[static 1],
                  size_t l, size_t orblen)
{
  if ( l == 0 ) return;
  for ( size_t i = 0; i < orblen; i++ )
    orbit[i] = apply_elem(orbit[i], (const uint32_t**)u, l);
}

static
int _orb_cmp(const void *_a, const void *_b, void *_c)
{
  uint32_t a = *(uint32_t*)_a;
  uint32_t b = *(uint32_t*)_b;
  uint32_t *order = (uint32_t*)_c;
  if ( order[a] < order[b] ) return -1;
  if ( order[a] > order[b] ) return 1;
  return 0;
}

static
void _setstab_in(size_t l, size_t j, size_t *skipto,
                 permgrp_t G, permgrp_t *K,
                 uint32_t *restrict u[restrict static 1],
                 uint32_t *restrict u_inv[restrict static 1],
                 sets_t *restrict fund_orbs[static 1],
                 uint32_t *restrict orb_decomp[restrict static 1],
                 const sets_t s[restrict static 1],
                 size_t m, uint32_t Kbase[restrict static 1])
{
  if ( l == j ) {
    uint32_t *p = arraytoelem((const uint32_t *restrict*)u, l, permgrp_degree(G));
    if ( perm_id_cmp(p, permgrp_degree(G)) < 0 ) {
      if ( *K == NULL )
        *K = permgrp_new((const uint32_t *restrict*)&p, 1, permgrp_degree(G),
                         NULL, 0);
      else permgrp_update(*K, p);
      permgrp_new_base(K, permgrp_base(G), permgrp_baselen(G));
      _update_fund_orbs(*K, j, (sets_t**)fund_orbs, m);
      permgrp_new_base(K, Kbase, l);
      *skipto = _fixed_base(G, p, j);
      for ( size_t i = 0; i <= *skipto; i++ )
        permgrp_orbrep_update(*K, orb_decomp[i], i);
    }
    free(p);
    return;
  }

  size_t i;
  uint32_t gamma;
  uint32_t *orbit = permgrp_fundorb(G, l);
  _orbit_apply(orbit, (const uint32_t *restrict*)u, l, permgrp_fundorb_size(G, l));
  qsort_r(orbit, permgrp_fundorb_size(G, l), sizeof(orbit[0]),
          _orb_cmp, permgrp_ordering(G));

  *skipto = l;
  if ( *K != NULL ) permgrp_orbrep_update(*K, orb_decomp[l], l);

  for ( i = 0; i < permgrp_fundorb_size(G, l) && l <= *skipto; i++ ) {
    if ( *K != NULL && orb_decomp[l][orbit[i]] != 1 ) continue;
    if ( !_ISELEMENT(s, orbit[i]) ) continue;
    gamma = apply_elem_inv(orbit[i], (const uint32_t**)u_inv, l);
    perm_cpy_inplace(u[l], permgrp_rep(G, l, gamma), permgrp_degree(G));
    perm_cpy_inplace(u_inv[l], permgrp_repinv(G, l, gamma), permgrp_degree(G));
    gamma = orbit[i];
    Kbase[l] = gamma;

    if ( *K != NULL ) {
      uint32_t *order = permgrp_ordering(G);
      if ( bit_setsize(fund_orbs[l], m) > 1 &&
           order[orbit[permgrp_fundorb_size(G, l) -
                       bit_setsize(fund_orbs[l], m) + 1]] <= order[gamma])
        continue;
      order = NULL;
      if ( _test(G, (const sets_t**)fund_orbs, (const uint32_t**)u, gamma, l) < 0 )
        continue;
      permgrp_new_base(K, Kbase, l + 1);
    }

    _setstab_in(l + 1, j, skipto, G, K, u, u_inv, fund_orbs, orb_decomp,
                s, m, Kbase);
  }

  free(orbit);
}

permgrp_t permgrp_setstab(permgrp_t G[static 1], sets_t s[static 1], size_t m)
{
  if ( G == NULL || *G == NULL ) return NULL;

  size_t i, j;

  /* rebase acting group */
  size_t set_size = bit_setsize(s, m);
  uint32_t *base = (uint32_t*)malloc(set_size * sizeof(uint32_t));
  uint32_t *old_base = (uint32_t*)malloc(permgrp_degree(*G) * sizeof(uint32_t));

  {
    size_t index = 0;
    for ( int z = -1; (z = bit_nextelement(s, m, z)) >= 0; ) base[index++] = z;
    uint32_t *Gbase = permgrp_base(*G);
    for ( i = 0; i < permgrp_degree(*G); i++ ) old_base[i] = Gbase[i];
    Gbase = NULL;
  }

  permgrp_new_base(G, base, set_size);

  /* get pointwise stabilizer of set */
  permgrp_t K = NULL;
  {
    arrays_t Agens = permgrp_generatorssbgrp(*G, set_size);
    if ( arrays_nmem(Agens) > 0 ) {
      uint32_t **gens = (uint32_t**)malloc(arrays_nmem(Agens) * sizeof(uint32_t*));
      size_t ngens = 0;
      for ( i = 0; i < arrays_nmem(Agens); i++ )
        if ( perm_id_cmp(arrays_at(Agens, i), permgrp_degree(*G)) < 0 )
          gens[ngens++] = arrays_at(Agens, i);

      if ( ngens > 0 )
        K = permgrp_new((const uint32_t *restrict*)gens, ngens,
                        permgrp_degree(*G), (const uint32_t *restrict)base,
                        set_size);

      free(gens);
    }
    arrays_free(&Agens);
  }

  uint32_t **u = (uint32_t**)malloc(set_size * sizeof(uint32_t*));
  uint32_t **u_inv = (uint32_t**)malloc(set_size * sizeof(uint32_t*));
  for ( i = 0; i < set_size; i++ ) {
    u[i] = (uint32_t*)malloc(permgrp_degree(*G) * sizeof(uint32_t));
    u_inv[i] = (uint32_t*)malloc(permgrp_degree(*G) * sizeof(uint32_t));
  }

  sets_t **fund_orbs = (sets_t**)malloc(set_size * sizeof(sets_t*));
  for ( i = 0; i < set_size; i++ )
    fund_orbs[i] = (sets_t*)malloc(m * sizeof(sets_t));
  if ( K != NULL ) _update_fund_orbs(K, set_size, fund_orbs, m);
  else {
    for ( i = 0; i < set_size; i++ ) {
      _EMPTYSET(fund_orbs[i], m);
      _ADDELEMENT(fund_orbs[i], base[i]);
    }
  }

  uint32_t **orb_decomp = (uint32_t**)malloc(set_size * sizeof(uint32_t*));
  for ( i = 0; i < set_size; i++ ) {
    orb_decomp[i] = (uint32_t*)malloc(permgrp_degree(*G) * sizeof(uint32_t));
    for ( j = 0; j < permgrp_degree(*G); j++ ) orb_decomp[i][j] = 1;
  }

  size_t skipto = set_size;
  uint32_t Kbase[set_size];

  _setstab_in(0, set_size, &skipto, *G, &K, u, u_inv, fund_orbs, orb_decomp,
              s, m, Kbase);

  permgrp_new_base(G, old_base, permgrp_degree(*G));

  for ( i = 0; i < set_size; i++ ) {
    free(orb_decomp[i]);
    free(fund_orbs[i]);
    free(u[i]);
    free(u_inv[i]);
  }
  free(orb_decomp);
  free(fund_orbs);
  free(u);
  free(u_inv);
  free(base);
  free(old_base);

  return K;
}
