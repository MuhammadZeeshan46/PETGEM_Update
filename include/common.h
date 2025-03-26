/*
 * Filename: common.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-05-28
 *
 * Description:
 * This file contains a collection of definitions for common utility functions that are used
 * throughout the PETGEM project. These functions include operations for printing and timers.
 *
 * List of definitions:
 * - void printHeader(): Prints PETGEM Header. 
 * - void printFooter(): Prints PETGEM Footer. 
 * - void createDirectory(): Creates a directory. 
 *
 * Usage:
 * Include this file in your source code to utilize the common functions. 
 * For example:
 * #include "common.h"
 * 
*/

#ifndef COMMON_H
#define COMMON_H

// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode printHeader();
PetscErrorCode printFooter();
PetscErrorCode createDirectory(const char *path);

#endif

