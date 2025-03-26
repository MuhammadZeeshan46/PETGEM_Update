/*
 * Filename: assembly.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-06
 *
 * Description:
 * This file contains a collection of definitions for assembly functions that are used
 * throughout the PETGEM project.
 *
 * List of definitions:
 *
 * Usage:
 * Include this file in your source code to utilize the assembly functions. 
 * For example:
 * #include "assembly.h"
 * 
*/

#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <petsc.h>

// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode assemblySystem(DM dm, Mat A, Vec b, userParams *params, petgemGrid *grid, petgemSource *source, Vec *resistivity_x, Vec *resistivity_y, Vec *resistivity_z);




#endif

