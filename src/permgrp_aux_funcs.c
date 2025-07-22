/**
 * @file permgrp_aux_funcs.c
 * @brief Auxiliary functions used in the implementation of the <tt>permgrp_t</tt>
 * class.
 *
 * Auxiliary functions used in the implementation of the <tt>permgrp_t</tt> class.
 *
 * @author Thomas Pender
 * @date 2025-07
 */
# include <permgrp.h>
# include <permgrp_struct_p.h>
# include <perm.h>

# include <stdio.h>

/* ---------------------------------------------------------------------------
 * fundamental orbit generation
 * ---------------------------------------------------------------------------*/

typedef struct {
  uint32_t elem, *p;
} orb_pair_t;

typedef struct {
  size_t n;
  tree_t *T;
  orb_pair_t *pair;
  stacks_t tmp;
} hash_arg_t;

static
int _orbit_hash_apply(void **_x, void *_args)
{
  uint32_t *x = *(uint32_t**)_x;
  hash_arg_t *args = (hash_arg_t*)_args;

  uint32_t elem = x[args->pair->elem];
  if ( args->T->T[elem] == NULL ) {
    args->T->T[elem] = perm_comp(args->pair->p, x, args->n);
    args->T->T_inv[elem] = perm_inv(args->T->T[elem], args->n);
    args->T->size++;

    orb_pair_t *new_pair = (orb_pair_t*)malloc(sizeof(orb_pair_t));
    new_pair->p = args->T->T[elem];
    new_pair->elem = elem;

    stacks_push(args->tmp, new_pair);
  }

  return 1;
}

static
int _orbit_stack_free(void **x)
{
  free(*x);
  *x = NULL;
  return 1;
}

void orbit(tree_t T[static 1], uint32_t base, size_t n)
{
  stacks_t tmp = stacks_new();
  if ( T->T[base] == NULL ) {
    T->T[base] = perm_id(n);
    T->T_inv[base] = perm_id(n);
    T->size = 1;
    orb_pair_t *pair = (orb_pair_t*)malloc(sizeof(orb_pair_t));
    pair->elem = base;
    pair->p = T->T[base];
    stacks_push(tmp, pair);
  }
  else {
    for ( uint32_t a = 0, count = 0; a < n && count < T->size; a++ )
      if ( T->T[a] != NULL ) {
        count ++;
        orb_pair_t *pair = (orb_pair_t*)malloc(sizeof(orb_pair_t));
        pair->elem = a;
        pair->p = T->T[a];
        stacks_push(tmp, pair);
      }
  }

  while ( !stacks_empty(tmp) ) {
    hash_arg_t *args = (hash_arg_t*)malloc(sizeof(hash_arg_t));
    args->pair = (orb_pair_t*)stacks_pop(tmp);
    args->n = n;
    args->T = T;
    args->tmp = tmp;

    hashtabs_map_r(T->S, _orbit_hash_apply, args);

    free(args->pair);
    free(args);
  }

  stacks_map(tmp, _orbit_stack_free);
  stacks_free(&tmp);
}

/* ---------------------------------------------------------------------------
 * hashing/rehashing
 * ---------------------------------------------------------------------------*/

void hash_insert(hashtabs_t *t, const uint32_t p[static 1], size_t n)
{
  if ( hashtabs_loadfactor(*t) > ((3* hashtabs_size(*t)) >> 2UL) ) {
    hashtabs_t new_t = hashtabs_rehash_r(*t, &n, &n);
    hashtabs_swap(&new_t, t);
    hashtabs_free(&new_t);
  }
  hashtabs_insert_r(*t, p, &n, &n);
}

/* ---------------------------------------------------------------------------
 * auxiliary functions
 * ---------------------------------------------------------------------------*/

uint32_t *sift(permgrp_t G, const uint32_t p[restrict static 1], size_t *j)
{
  uint32_t *q = perm_cpy(p, G->degree);
  for ( size_t i = 0; i < G->base_size; i++ ) {
    uint32_t elem = q[G->base[i]];
    if ( G->trees[i]->T[elem] == NULL ) {
      *j = i;
      return q;
    }
    perm_comp_inplace(q, G->trees[i]->T_inv[elem], G->degree);
  }
  *j = G->base_size;
  return q;
}

void swap(void *a, void *b, size_t nbytes)
{
  char tmp[nbytes];
  memcpy(tmp, a, nbytes);
  memcpy(a, b, nbytes);
  memcpy(b, tmp, nbytes);
}

static
int print_hash_apply(void **_x, void *_y)
{
  size_t y = *(size_t*)_y;
  uint32_t *x = *(uint32_t**)_x;
  for ( size_t i = 0; i < y; i++ ) printf("%u ", x[i]);
  printf("\n");
  return 1;
}

void print_grp_info(permgrp_t G)
{
  size_t i, j;
  printf("grp order = %zu\n", permgrp_order(G));
  printf("degree = %zu\nbase length = %zu\n",
         permgrp_degree(G), permgrp_baselen(G));
  printf("base = ");
  for ( i = 0; i < G->base_size; i++ ) printf("%u ", G->base[i]);
  printf("\nordering = ");
  for ( i = 0; i < G->degree; i++ ) printf("%u ", G->order[i]);
  printf("\n\n");
  for ( i = 0; i < G->degree; i++ )
    if ( G->trees[i] != NULL ) {
      printf("tree %zu:\nsize = %zu\nbase elem = %u\n",
             i, G->trees[i]->size, G->base[i]);
      for ( j = 0; j < G->degree; j++ )
        if ( G->trees[i]->T[j] != NULL ) {
          printf("%2zu: ", j);
          perm_pr(G->trees[i]->T[j], G->degree);
          printf("\n    ");
          perm_pr(G->trees[i]->T_inv[j], G->degree);
          printf("\n");
        }
      printf("tree generators:\n");
      hashtabs_map_r(G->trees[i]->S, print_hash_apply, &G->degree);
      printf("\n");
    }
}
