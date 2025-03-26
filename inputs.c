/*
 * Filename: inputs.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-05-28
 *
 * Description:
 * This file contains functions for handling input data and user parameters in a PETGEM simulation. 
 * It includes functions for parsing input data, and processing user-provided parameters. 
 * The functions in this file facilitate the setup and configuration of the PETGEM code.
 *
 * List of Functions:
 * - void readUserParams(): Read user parameters for PETGEM simulation. 
*/

/* C libraries */ 
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>

/* PETSc libraries */
#include <petscsys.h>

/* PETGEM functions */ 
#include "common.h"
#include "inputs.h"  

// =============================================================================
// Function: readUserParams
// =============================================================================
PetscErrorCode readUserParams(userParams *params, PetscMPIInt size){
    /*
     * Function: readUserParams
     * ----------------------------
     *   Reads user parameters from command line options and sets the values in the provided Params struct.
     *
     *   This function retrieves the user-specified parameters for `nord`, `csem_source`, and `mt_source`
     *   from the command line options. It validates the values and ensures only one of `csem_source` 
     *   or `mt_source` is specified. If the `nord` parameter is outside the valid range (1 to 6), it 
     *   defaults to 1. The `mode` is set based on the presence of `csem_source` or `mt_source`.
     *
     *   Parameters:
     *     params - A pointer to a Params struct where the read values will be stored.
     *
     *   Returns:
     *     PetscErrorCode - Error code returned by PETSc functions. Returns 0 (PETSC_SUCCESS) on successful completion.
     *
     *   Preconditions:
     *     The PETSc library must be initialized before calling this function.
     *
     *   Postconditions:
     *     The `params` struct will contain the validated and processed user parameters.
     *
     *   Example usage:
     */
    
    PetscFunctionBeginUser;
    /* Declarations */
    char          meshFilename[PETSC_MAX_PATH_LEN];
    PetscBool     meshFilenameIsPresent;
    char          resistivityFilename[PETSC_MAX_PATH_LEN];
    PetscBool     resistivityFilenameIsPresent;
    char          receiversFilename[PETSC_MAX_PATH_LEN];
    PetscBool     receiversFilenameIsPresent;  
    char          outputDir[PETSC_MAX_PATH_LEN];
    PetscBool     outputDirIsPresent;  
    PetscInt      nord; // Basis order = 1, 2, 3, 4, 5, 6
    PetscBool     nordIsPresent;
    
    /* Read user parameters */
    PetscCall(PetscOptionsGetString(NULL, NULL, "-mesh_dm_plex_filename", meshFilename, sizeof(meshFilename), &meshFilenameIsPresent));  
    // Check user input
    if (!meshFilenameIsPresent) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Mesh file missing. Mandatory parameter required for simulation.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_NULL);
    }
    PetscCall(PetscStrncpy(params->meshFile, meshFilename, sizeof(params->meshFile)));

    PetscCall(PetscOptionsGetString(NULL, NULL, "-resistivity_filename_x", resistivityFilename, sizeof(resistivityFilename), &resistivityFilenameIsPresent));  
    // Check user input
    if (!resistivityFilenameIsPresent) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Resistivity file for x-component missing. Mandatory parameter required for simulation.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_NULL);
    }
    PetscCall(PetscStrncpy(params->resistivityFile_x, resistivityFilename, sizeof(params->resistivityFile_x)));

    PetscCall(PetscOptionsGetString(NULL, NULL, "-resistivity_filename_y", resistivityFilename, sizeof(resistivityFilename), &resistivityFilenameIsPresent));  
    // Check user input
    if (!resistivityFilenameIsPresent) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Resistivity file for y-component missing. Mandatory parameter required for simulation.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_NULL);
    }
    PetscCall(PetscStrncpy(params->resistivityFile_y, resistivityFilename, sizeof(params->resistivityFile_y)));

    PetscCall(PetscOptionsGetString(NULL, NULL, "-resistivity_filename_z", resistivityFilename, sizeof(resistivityFilename), &resistivityFilenameIsPresent));  
    // Check user input
    if (!resistivityFilenameIsPresent) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Resistivity file for z-component missing. Mandatory parameter required for simulation.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_NULL);
    }
    PetscCall(PetscStrncpy(params->resistivityFile_z, resistivityFilename, sizeof(params->resistivityFile_z)));

    PetscCall(PetscOptionsGetString(NULL, NULL, "-receivers_filename", receiversFilename, sizeof(receiversFilename), &receiversFilenameIsPresent));  
    // Check user input
    if (!receiversFilenameIsPresent) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Receivers file missing. Mandatory parameter required for simulation.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_NULL);
    }
    PetscCall(PetscStrncpy(params->receiversFile, receiversFilename, sizeof(params->receiversFile)));
    
    PetscCall(PetscOptionsGetString(NULL, NULL, "-output_dir", outputDir, sizeof(outputDir), &outputDirIsPresent));  
    // Check user input
    if (!outputDirIsPresent) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: Output directory missing. Mandatory parameter required for simulation.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_NULL);
    }
    PetscCall(PetscStrncpy(params->outputDirectory, outputDir, sizeof(params->outputDirectory)));
    
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-nord", &nord, &nordIsPresent));
    // Check user input
    if (nord < 1 || nord > 6) {
        nord = 1;
    }
    params->nord = nord;

    // Check user input
    PetscBool    csemSourceIsPresent, mtSourceIsPresent;
    PetscCall(PetscOptionsHasName(NULL, NULL, "-csem_source", &csemSourceIsPresent));
    PetscCall(PetscOptionsHasName(NULL, NULL, "-mt_source", &mtSourceIsPresent));
    
    if ((csemSourceIsPresent && mtSourceIsPresent) || !(csemSourceIsPresent || mtSourceIsPresent)) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "At least one (but not both) of the parameters csem_source or mt_source must be specified.\n"));
        PetscFunctionReturn(PETSC_ERR_ARG_IDN);       
    }

    /* csem mode = 0; mt mode = 1*/
    params->mode = csemSourceIsPresent ? 0 : 1;

    /* Number of MPI tasks */ 
    params->numMPITasks = size;
   
    /* Create output directory */
    createDirectory(outputDir);

    /* Print user params */ 
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Model data:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Mesh filename                     = %s\n", params->meshFile));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Resistivity x-component filename  = %s\n", params->resistivityFile_x));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Resistivity y-component filename  = %s\n", params->resistivityFile_y));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Resistivity z-component filename  = %s\n", params->resistivityFile_z));

    PetscFunctionReturn(PETSC_SUCCESS);

}
