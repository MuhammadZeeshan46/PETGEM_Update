/*
 * Filename: common.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-05-28
 *
 * Description:
 * This file contains a collection of functions for common utility functions that are used
 * throughout the PETGEM project. These functions include operations for printing and timers.
 *
 * List of functions:
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

/* C libraries */ 
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>

/* PETSc libraries */
#include <petscsys.h>

/* PETGEM funcions*/
#include "common.h"
#include "version.h"

// =============================================================================
// Function: printHeader
// =============================================================================

PetscErrorCode printHeader(){
    /*
    * Function: printHeader
    * ----------------------------
    *   Prints a header with PETGEM project information and the current year.
    *
    *   This function prints a formatted header to the console, containing
    *   information about the PETGEM project, including its name, purpose,
    *   GitHub repository, website, and the names and affiliations of the
    *   developers. It also includes the current year.
    *
    *   The header is printed using PETSc's parallel printing functions,
    *   ensuring that the output is consistent across all processes in the
    *   PETSc communicator.
    *
    *   Parameters:
    *     None
    *
    *   Returns:
    *     PetscErrorCode - Error code returned by PETSc functions. Returns
    *                      PETSC_SUCCESS on successful completion.
    *
    */

    PetscFunctionBeginUser;

    time_t now;
    time(&now);
    struct tm *local = localtime(&now);
    int year = local->tm_year + 1900;
    
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------------------------\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                                                                          -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                                   PETGEM                                 -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-          Parallel Edge-based Tool for Electromagnetic Modelling          -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                                                                          -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-          GitHub Repository: github.com/ocastilloreyes/petgem             -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                      Website: https://petgem.bsc.es/                     -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                                                                          -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------------------------\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                                                                          -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                         Octavio Castillo-Reyes                           -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-            Universitat Politècnica de Catalunya (UPC) - %d             -\n", year));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-               Barcelona Supercomputing Center (BSC) - %d               -\n", year));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "-                                                                          -\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------------------------\n"));

    PetscFunctionReturn(PETSC_SUCCESS);

}

// =============================================================================
// Function: printFooter
// =============================================================================

PetscErrorCode printFooter(){
    /*
    * Function: printFooter
    * ----------------------------
    *   Prints the finalization time of the PETGEM simulation in the format:
    *   "Finished: YYYY-MM-DD HH:MM:SS".
    *
    *   This function retrieves the current local time and prints it to the
    *   console in a formatted string. The printing is done using PETSc's
    *   parallel printing functions to ensure consistency across all processes
    *   in the PETSc communicator.
    *
    *   Parameters:
    *     None
    *
    *   Returns:
    *     PetscErrorCode - Error code returned by PETSc functions. Returns
    *                      PETSC_SUCCESS on successful completion.
    *
    */

    PetscFunctionBeginUser;

    time_t now;
    time(&now);
    struct tm *local = localtime(&now);
    
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Finished: %04d-%02d-%02d %02d:%02d:%02d", local->tm_year + 1900, local->tm_mon + 1, local->tm_mday, local->tm_hour, local->tm_min, local->tm_sec)); 
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n PETGEM version: %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH)); 
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------------------------\n"));

    PetscFunctionReturn(PETSC_SUCCESS);

}

// =============================================================================
// Function: createDirectory
// =============================================================================

PetscErrorCode createDirectory(const char *path) {
    /*
    * Function: createDirectory
    * ----------------------------
    *   Creates a directory if it does not already exist.
    *
    *   This function checks if the specified directory exists. If the directory
    *   does not exist, it attempts to create it with the specified permissions.
    *   The function uses PETSc error handling to report any issues encountered
    *   during the creation of the directory.
    *
    *   Parameters:
    *     path - A constant character pointer to the path of the directory to be created.
    *
    *   Returns:
    *     PetscErrorCode - Error code returned by PETSc functions. Returns
    *                      PETSC_SUCCESS on successful completion. If the
    *                      directory creation fails, it returns an appropriate
    *                      PETSc error code.
    */
    
     
    PetscFunctionBeginUser;
    
    /* Verify if the directory exists */
    struct stat st;
    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) { // Directory exists
            PetscFunctionReturn(PETSC_SUCCESS);
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, " Path exists but is not a directory: %s", path);
        }
    } else {
        /* Directory doesn't exist, create it */
        if (mkdir(path, 0755) && errno != EEXIST) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, " Error when creating output directory: %s", path);
        }
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}
