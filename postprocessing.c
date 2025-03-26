/*
 * Filename: receivers.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-04
 *
 * Description:
 * This file contains functions for receivers (interpolation points) in a PETGEM simulation. 
 * It includes functions for parsing receivers data based on user-provided parameters. 
 * The functions in this file facilitate the setup and configuration of the PETGEM code.
 *
 * List of Functions:
 * 
*/

/* C libraries */ 


/* PETSc libraries */
#include <petscsys.h>

/* PETGEM functions */ 
#include "constants.h"
#include "inputs.h"
#include "grid.h"
#include "hvfem.h"
#include "postprocessing.h"

// =============================================================================
// Function: computeFields
// =============================================================================
PetscErrorCode computeFields(DM dm, Vec x, userParams *params, petgemGrid *grid){
    PetscFunctionBeginUser;
    
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Compute electric and magnetic fields:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Receivers filename        = %s\n", params->receiversFile));
        
    /* Variable declarations */
    PetscReal   jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS];
    PetscReal   **basisFunctions, **curlBasisFunctions, *XiEtaZeta;
    PetscReal   realReceiverCoords[NUM_DIMENSIONS];
    PetscScalar *cellCoords = NULL, *scalarReceiverCoords, *closureReceiver;
    PetscScalar tmpFields[6], coords[NUM_DIMENSIONS];
    PetscInt    numCoords, cell, globalSizeReceivers, numGlobalReceivers; //transitiveClosureCellSize, 
    PetscInt    numReceiversFoundGlobal, numReceiversFoundLocal;
    //PetscInt    *transitiveClosureCellPoints = NULL;
    PetscInt    closureSizeReceiver = grid->numDofInCell;
    PetscInt    cellOrientation[10], indexes[NUM_DIMENSIONS];
    PetscBool   isDG;
    PetscSF     receiverGlobalSF = NULL, receiverLocalSF = NULL;    
    PetscMPIInt rank;
    Vec         xLocal, receivers, receiverCoords, Ex, Ey, Ez, Hx, Hy, Hz;
    PetscViewer viewer; 
    
    char  fieldFileName[PETSC_MAX_PATH_LEN], outFileName[PETSC_MAX_PATH_LEN];
    const PetscSFNode   *receiverInCell;
    const PetscInt      *receiverFound;   
    const PetscScalar *arrayCoords;


    /* Load receiver data (sequential) */     
    PetscCall(VecCreate(PETSC_COMM_SELF, &receivers));
    PetscCall(VecSetFromOptions(receivers));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, params->receiversFile, FILE_MODE_READ, &viewer));
    PetscCall(VecLoad(receivers, viewer));
    PetscCall(VecSetBlockSize(receivers, NUM_DIMENSIONS));
    PetscCall(VecGetSize(receivers, &globalSizeReceivers));
    PetscCall(PetscViewerDestroy(&viewer));

    /* Verify receivers vector consistency */
    if (globalSizeReceivers % 3 != 0) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Error: Global size of the receivers vector (%d) is not divisible by 3, which is required for 3D points.\n", globalSizeReceivers););
        PetscFunctionReturn(PETSC_ERR_ARG_SIZ);       
    }
    else{
        numGlobalReceivers = globalSizeReceivers/3;
    }

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Number of receivers       = %d\n", numGlobalReceivers));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Postprocessing status     = started\n"));

    /* Check if all the receivers are within the computational domain */
    PetscCall(DMLocatePoints(dm, receivers, DM_POINTLOCATION_NONE, &receiverGlobalSF));
    PetscCall(PetscSFGetGraph(receiverGlobalSF, NULL, &numReceiversFoundLocal, NULL, NULL));
    MPI_Reduce(&numReceiversFoundLocal, &numReceiversFoundGlobal, 1, MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD);
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

    if (rank==0){
        if (numReceiversFoundGlobal <= 0){
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Error: Receiver coordinates are either not found or located outside the computational domain.\n"));
            PetscFunctionReturn(PETSC_ERR_ARG_WRONG);    
        }
        else if (numReceiversFoundGlobal != numGlobalReceivers){
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Warning: Some receiver coordinates are either not found or located outside the computational domain. Fields will not be computed for these receivers.\n"));
        }
    }

    /* Prepare receiver data for point location within loop */ 
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, NUM_DIMENSIONS, &receiverCoords));
    PetscCall(VecSetBlockSize(receiverCoords, NUM_DIMENSIONS));

    /* Create vectors for electric and magnetic fields */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &Ex));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &Ey));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &Ez));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &Hx));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &Hy));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &Hz));

    PetscCall(VecSetSizes(Ex, PETSC_DECIDE, numGlobalReceivers));
    PetscCall(VecSetSizes(Ey, PETSC_DECIDE, numGlobalReceivers));
    PetscCall(VecSetSizes(Ez, PETSC_DECIDE, numGlobalReceivers));
    PetscCall(VecSetSizes(Hx, PETSC_DECIDE, numGlobalReceivers));
    PetscCall(VecSetSizes(Hy, PETSC_DECIDE, numGlobalReceivers));
    PetscCall(VecSetSizes(Hz, PETSC_DECIDE, numGlobalReceivers));

    PetscCall(VecSetFromOptions(Ex));
    PetscCall(VecSetFromOptions(Ey));
    PetscCall(VecSetFromOptions(Ez));
    PetscCall(VecSetFromOptions(Hx));
    PetscCall(VecSetFromOptions(Hy));
    PetscCall(VecSetFromOptions(Hz));

    PetscCall(VecSetUp(Ex));
    PetscCall(VecSetUp(Ey));
    PetscCall(VecSetUp(Ez));
    PetscCall(VecSetUp(Hx));
    PetscCall(VecSetUp(Hy));
    PetscCall(VecSetUp(Hz));
    
    /* Allocate memory */
    PetscCall(PetscMalloc1(grid->numDofInCell, &closureReceiver));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &XiEtaZeta));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &basisFunctions));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &curlBasisFunctions));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(grid->numDofInCell, &basisFunctions[i]));
        PetscCall(PetscCalloc1(grid->numDofInCell, &curlBasisFunctions[i]));
    }

    /* Get local vector */
    PetscCall(DMGetLocalVector(dm, &xLocal));
    PetscCall(DMGlobalToLocal(dm, x, INSERT_VALUES, xLocal));

    PetscSection section;
    PetscCall(DMGetLocalSection(dm, &section));


    /* Compute electric/magnetic fields for each receiver */
    for (PetscInt i = 0; i < numGlobalReceivers; ++i){
        /* Get array data */
        PetscCall(VecGetArrayWrite(receiverCoords, &scalarReceiverCoords));
        
        /* Compute indexes  */
        indexes[0] = NUM_DIMENSIONS*i;      // x-coordinate index
        indexes[1] = NUM_DIMENSIONS*i + 1;  // y-coordinate index
        indexes[2] = NUM_DIMENSIONS*i + 2;  // z-coordinate index

        /* Get and setup receiver coordinates */
        PetscCall(VecGetValues(receivers, NUM_DIMENSIONS, indexes, coords));
        scalarReceiverCoords[0] = coords[0];
        scalarReceiverCoords[1] = coords[1];
        scalarReceiverCoords[2] = coords[2];
        realReceiverCoords[0] = PetscRealPart(coords[0]);
        realReceiverCoords[1] = PetscRealPart(coords[1]);
        realReceiverCoords[2] = PetscRealPart(coords[2]);

        /* Restore array data */ 
        PetscCall(VecRestoreArrayWrite(receiverCoords, &scalarReceiverCoords));

        /* Locate receiver within computational domain */
        PetscCall(DMLocatePoints(dm, receiverCoords, DM_POINTLOCATION_NONE, &receiverLocalSF));
        PetscCall(PetscSFGetGraph(receiverLocalSF, NULL, &numReceiversFoundLocal, &receiverFound, &receiverInCell));
        PetscCall(PetscSFDestroy(&receiverLocalSF));

        /* If receiver is located, compute fields */
        for (PetscInt j=0; j<numReceiversFoundLocal; j++){
            /* Get cell index in which receiver belongs */
            cell = receiverInCell[j].index;
            
            /* Get cell coordinates */ 
            PetscCall(DMPlexGetCellCoordinates(dm, cell, &isDG, &numCoords, &arrayCoords, &cellCoords));

            /* Compute jacobian and its inverse for receiverInCell */
            PetscCall(computeJacobian(cellCoords, jacobian, invJacobian));
            
            /* Transform xyz receiver position to XiEtaZeta coordinates (reference tetrahedral element) */
            PetscCall(tetrahedronXYZToXiEtaZeta(cellCoords, realReceiverCoords, XiEtaZeta));
            
            /* Restore cell coordinates */ 
            PetscCall(DMPlexRestoreCellCoordinates(dm, cell, &isDG, &numCoords, &arrayCoords, &cellCoords));

            /* Get orientations for faces, edges and vertices of receiverInCell */
            /* Orden convention:
                - Faces indices start on position 2
                - Edges indices start on position 2 + NUM_FACES_PER_ELEMENT * 2
                - Vertices indices start on position 2 + NUM_FACES_PER_ELEMENT * 2 + NUM_EDGES_PER_ELEMENT * 2
            */
            /* Compute cell orientation */
            for(PetscInt j = 0; j < 10; j++){
                cellOrientation[j] = 0;
            }

            PetscCall(computeCellOrientation(dm, cell, cellOrientation));

            /* Compute basis functions for receiverInCell */ 
            PetscCall(computeBasisFunctions(params->nord, cellOrientation, jacobian, invJacobian, XiEtaZeta, basisFunctions, curlBasisFunctions));

            /* Get clousure for receiverInCell */
            PetscCall(DMPlexVecGetClosure(dm, section, xLocal, cell, &closureSizeReceiver, &closureReceiver));

            /* Reset variables to zero */
            for (PetscInt k=0; k<grid->numDofInCell; k++){
                tmpFields[k] = 0.0 + 0.0*PETSC_i;
            }

            /* Interpolate fields at receiver i 
                tmpFields[0] = Ex
                tmpFields[1] = Ey
                tmpFields[2] = Ez
                tmpFields[3] = Hx
                tmpFields[4] = Hy
                tmpFields[5] = Hz
            */
            for (PetscInt k=0; k<grid->numDofInCell; k++){
                tmpFields[0] += (basisFunctions[0][k] * closureReceiver[k]);
                tmpFields[1] += (basisFunctions[1][k] * closureReceiver[k]);
                tmpFields[2] += (basisFunctions[2][k] * closureReceiver[k]);
                tmpFields[3] += (curlBasisFunctions[0][k] * closureReceiver[k]);
                tmpFields[4] += (curlBasisFunctions[1][k] * closureReceiver[k]);
                tmpFields[5] += (curlBasisFunctions[2][k] * closureReceiver[k]);
            }
       
            /* Set values to output vectors */
            PetscCall(VecSetValue(Ex, i, tmpFields[0], INSERT_VALUES));
            PetscCall(VecSetValue(Ey, i, tmpFields[1], INSERT_VALUES));
            PetscCall(VecSetValue(Ez, i, tmpFields[2], INSERT_VALUES));
            PetscCall(VecSetValue(Hx, i, tmpFields[3], INSERT_VALUES));
            PetscCall(VecSetValue(Hy, i, tmpFields[4], INSERT_VALUES));
            PetscCall(VecSetValue(Hz, i, tmpFields[5], INSERT_VALUES));            
        }
    }

    /* Restore local and global vector */ 
    PetscCall(DMRestoreLocalVector(dm, &xLocal));

    /* Perform global assembly */ 
    PetscCall(VecAssemblyBegin(Ex));
    PetscCall(VecAssemblyBegin(Ey));
    PetscCall(VecAssemblyBegin(Ez));
    PetscCall(VecAssemblyBegin(Hx));
    PetscCall(VecAssemblyBegin(Hy));
    PetscCall(VecAssemblyBegin(Hz));

    PetscCall(VecAssemblyEnd(Ex));
    PetscCall(VecAssemblyEnd(Ey));
    PetscCall(VecAssemblyEnd(Ez));
    PetscCall(VecAssemblyEnd(Hx));
    PetscCall(VecAssemblyEnd(Hy));
    PetscCall(VecAssemblyEnd(Hz));

    /* Write output vectors */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Electric field components:\n"));
    PetscCall(PetscStrncpy(fieldFileName, "/Ex.dat", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, params->outputDirectory, sizeof(outFileName)));
    PetscCall(PetscStrlcat(outFileName, fieldFileName, sizeof(outFileName)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "      X-component (Ex)       = %s\n", outFileName));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, outFileName, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(Ex, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscStrncpy(fieldFileName, " ", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, " ", sizeof(outFileName)));

    PetscCall(PetscStrncpy(fieldFileName, "/Ey.dat", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, params->outputDirectory, sizeof(outFileName)));
    PetscCall(PetscStrlcat(outFileName, fieldFileName, sizeof(outFileName)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "      Y-component (Ey)       = %s\n", outFileName));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, outFileName, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(Ey, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscStrncpy(fieldFileName, " ", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, " ", sizeof(outFileName)));

    PetscCall(PetscStrncpy(fieldFileName, "/Ez.dat", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, params->outputDirectory, sizeof(outFileName)));
    PetscCall(PetscStrlcat(outFileName, fieldFileName, sizeof(outFileName)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "      Z-component (Ez)       = %s\n", outFileName));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, outFileName, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(Ez, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscStrncpy(fieldFileName, " ", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, " ", sizeof(outFileName)));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Magnetic field components:\n"));
    PetscCall(PetscStrncpy(fieldFileName, "/Hx.dat", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, params->outputDirectory, sizeof(outFileName)));
    PetscCall(PetscStrlcat(outFileName, fieldFileName, sizeof(outFileName)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "      X-component (Hx)       = %s\n", outFileName));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, outFileName, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(Hx, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscStrncpy(fieldFileName, " ", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, " ", sizeof(outFileName)));

    PetscCall(PetscStrncpy(fieldFileName, "/Hy.dat", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, params->outputDirectory, sizeof(outFileName)));
    PetscCall(PetscStrlcat(outFileName, fieldFileName, sizeof(outFileName)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "      Y-component (Hy)       = %s\n", outFileName));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, outFileName, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(Hy, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscStrncpy(fieldFileName, " ", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, " ", sizeof(outFileName)));

    PetscCall(PetscStrncpy(fieldFileName, "/Hz.dat", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, params->outputDirectory, sizeof(outFileName)));
    PetscCall(PetscStrlcat(outFileName, fieldFileName, sizeof(outFileName)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "      Z-component (Hz)       = %s\n", outFileName));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, outFileName, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(Hz, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscStrncpy(fieldFileName, " ", sizeof(fieldFileName)));
    PetscCall(PetscStrncpy(outFileName, " ", sizeof(outFileName)));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Postprocessing status     = completed\n"));
    
    /* Free memory */
    PetscCall(PetscSFDestroy(&receiverGlobalSF));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(curlBasisFunctions[i]));
        PetscCall(PetscFree(basisFunctions[i]));        
    }
    PetscCall(PetscFree(curlBasisFunctions));
    PetscCall(PetscFree(basisFunctions));    
    PetscCall(PetscFree(closureReceiver));    
    PetscCall(PetscFree(XiEtaZeta));
    PetscCall(VecDestroy(&receiverCoords));
    PetscCall(VecDestroy(&receivers));    
    PetscCall(VecDestroy(&Ex));
    PetscCall(VecDestroy(&Ey));
    PetscCall(VecDestroy(&Ez));
    PetscCall(VecDestroy(&Hx));
    PetscCall(VecDestroy(&Hy));
    PetscCall(VecDestroy(&Hz));
    
    
    //PetscCall(VecDestroy(&xLocal));
    
      
    PetscFunctionReturn(PETSC_SUCCESS);
}
    
