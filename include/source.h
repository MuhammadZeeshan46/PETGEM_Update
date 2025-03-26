/*
 * Filename: source.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-01
 *
 * Description:
 * This file contains a collection of definitions for source functions that are used
 * throughout the PETGEM project.
 *
 * List of definitions:
 * -  
 *
 * Usage:
 * Include this file in your source code to utilize the common functions. 
 * For example:
 * #include "source.h"
 * 
*/

#ifndef SOURCE_H
#define SOURCE_H

#include <petsc.h>

// =============================================================================
// Declaration of structures
// =============================================================================
typedef struct {
    PetscReal position[3]; // Source position (x, y, z)
    PetscReal freq;        // Frequency 
    PetscReal current;     // Electric current
    PetscReal length;      // Dipole length
    PetscReal dip;         // Dip
    PetscReal azimuth;     // Azimuth
    char polarization[PETSC_MAX_PATH_LEN];
} petgemSource;

// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode setupSource(petgemSource *source, int mode);

#endif 





