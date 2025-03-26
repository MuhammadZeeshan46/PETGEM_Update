/*
 * Filename: grid.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-04
 *
 * Description:
 * This file contains a collection of definitions for grid functions that are used
 * throughout the PETGEM project. These functions are based on DMPLex provided by PETSc.
 *
 * List of definitions:
 * 
 *
 * Usage:
 * Include this file in your source code to utilize the grid functions. 
 * For example:
 * #include "solver.h"
 * 
*/

#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

/* PETGEM functions */ 

#include <petsc.h>

// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode computeFields(DM dm, Vec x, userParams *params, petgemGrid *grid);


#endif
