/*
 * Filename: inputs.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-05-28
 *
 * Description:
 * This file contains a collection of definitions for user input functions that are used
 * throughout the PETGEM project. These functions include operations for printing and timers.
 *
 * List of definitions:
 * - void readUserParams(): Prints PETGEM Header. 
 *
 * Usage:
 * Include this file in your source code to utilize the common functions. 
 * For example:
 * #include "inputs.h"
 * 
*/

#ifndef INPUTS_H
#define INPUTS_H

// =============================================================================
// Declaration of structures
// =============================================================================
typedef struct {
    char meshFile[PETSC_MAX_PATH_LEN];
    char resistivityFile_x[PETSC_MAX_PATH_LEN];
    char resistivityFile_y[PETSC_MAX_PATH_LEN];
    char resistivityFile_z[PETSC_MAX_PATH_LEN];
    char receiversFile[PETSC_MAX_PATH_LEN];
    char outputDirectory[PETSC_MAX_PATH_LEN];
    PetscInt nord;
    PetscInt mode;
    PetscMPIInt numMPITasks;
} userParams;


// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode readUserParams(userParams *params, PetscMPIInt size);

#endif

