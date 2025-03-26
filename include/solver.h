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

#ifndef SOLVER_H
#define SOLVER_H

/* PETGEM functions */ 

#include <petsc.h>

// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode solveSystem(DM dm, Mat A, Vec b, Vec x);


#endif
