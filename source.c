/*
 * Filename: source.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-03
 *
 * Description:
 * This file contains functions for source (CSEM or MT) in a PETGEM simulation. 
 * It includes functions for parsing source data based on user-provided parameters. 
 * The functions in this file facilitate the setup and configuration of the PETGEM code.
 *
 * List of Functions:
 * 
*/

/* C libraries */ 


/* PETSc libraries */
#include <petscsys.h>

/* PETGEM functions */ 
#include "source.h"

// =============================================================================
// Function: setupSource
// =============================================================================
PetscErrorCode setupSource(petgemSource *source, int mode) {
    PetscFunctionBeginUser;
    
    /* CSEM source */
    if (mode == 0) { 
        PetscScalar  dataSource[8];
        PetscInt     numSourceEntries=8;
        PetscBool    csemFlg;

        PetscCall(PetscOptionsGetScalarArray(NULL, NULL, "-csem_source", dataSource, &numSourceEntries, &csemFlg));  

        // Get source parameters
        if (csemFlg && numSourceEntries == 8) {
            source->position[0] = dataSource[0];    // x-position
            source->position[1] = dataSource[1];    // y-position
            source->position[2] = dataSource[2];    // z-position
            source->freq    = dataSource[3];        // Frequency    
            source->current = dataSource[4];        // Electric current
            source->length  = dataSource[5];        // Dipole length
            source->dip     = dataSource[6];        // Dipole dip
            source->azimuth = dataSource[7];        // Dipole azimuth
            source->polarization[0] = '\0';         // Not applicable for CSEM
        } else {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: The number of entries for -csem_source must be 8.\n"));
            PetscFunctionReturn(PETSC_ERR_ARG_OUTOFRANGE);
        }    
    }
    /* MT source */
    else if (mode == 1) { 
        char        *dataSource[2];
        PetscInt    numSourceEntries=2;
        PetscBool   mtFlg;

        PetscCall(PetscOptionsGetStringArray(NULL, NULL, "-mt_source", dataSource, &numSourceEntries, &mtFlg));

        // Get source parameters
        if (mtFlg && numSourceEntries == 2) {
            source->freq            = atof(dataSource[0]);  // Frequency (Convert string to PETSc scalar)
            PetscCall(PetscStrncpy(source->polarization, dataSource[1], sizeof(source->polarization)));
            PetscCall(PetscFree(dataSource[0]));
            PetscCall(PetscFree(dataSource[1]));
        } else {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: The number of entries for -mt_source must be 2.\n"));
            PetscFunctionReturn(PETSC_ERR_ARG_OUTOFRANGE);
        }


    } else {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Source type not supported.\n"));
        PetscFunctionReturn(PETSC_ERR_SUP);
    }

    /* Print source data */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Source data:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Type: %s\n", mode == 0 ? "CSEM" : "MT"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Freq (Hz)       = %f\n", source->freq));
    if (mode == 0) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Position (xyz)  = [%f, %f, %f]\n", source->position[0], source->position[1], source->position[2]));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Current         = %f\n", source->current));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Length          = %f\n", source->length));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Dip             = %f\n", source->dip));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Azimuth         = %f\n", source->azimuth));
    } else {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Polarization    = %s\n", source->polarization));
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}
   
