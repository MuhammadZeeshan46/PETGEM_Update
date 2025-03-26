/*
 * Filename: constats.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-12
 *
 * Description:
 * This file contains a collection of constants that are used throughout the PETGEM project.
 *
 * List of definitions:
 *
 * Usage:
 * Include this file in your source code to utilize the constants. 
 * For example:
 * #include "constants.h"
 * 
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <petsc.h>

// =============================================================================
// Declaration of constants
// =============================================================================

#define NUM_VERTICES_PER_ELEMENT  4
#define NUM_EDGES_PER_ELEMENT     6
#define NUM_FACES_PER_ELEMENT     4
#define NUM_EDGES_PER_FACE        3
#define NUM_VERTICES_PER_EDGE     2
#define NUM_VERTICES_PER_FACE     3
#define NUM_DIMENSIONS 			  3

#endif
