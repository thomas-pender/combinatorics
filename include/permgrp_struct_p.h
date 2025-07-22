/**
 * @file permgrp_struct_p.h
 * @brief Struct definition for <tt>permgrp_t</tt> class.
 *
 * Struct definition for <tt>permgrp_t</tt> class.
 *
 * @author Thomas Pender
 * @date 2025-07
 * @copyright GNU Public License
 */
# ifndef INCLUDED_PERMGRP_STRUCT_P_H
# define INCLUDED_PERMGRP_STRUCT_P_H

# include <stdint.h>
# include <stddef.h>

# include <containers/containers.h>

/**
 * @brief Schreier tree data structure.
 *
 * Schreier tree data structure.
 */
typedef struct {
  uint32_t **T;     //< transversal elements
  uint32_t **T_inv; ///< inverses of transversal elements
  size_t size;      ///< size of fundamental orbit
  hashtabs_t S;     ///< generators of stabilizer
} tree_t;

/**
 * @brief Permutation group data structure.
 *
 * Permutation group data structure.
 */
struct permgrp_t {
  uint32_t *order;  ///< point ordering induced by base
  uint32_t *base;   ///< base of permutation group
  uint32_t *id;     ///< group identity
  arrays_t gens;    ///< strong generating set
  tree_t **trees;   ///< list of Schreier trees
  size_t degree;    ///< degree of permutation group
  size_t base_size; ///< length of base
};

# endif
