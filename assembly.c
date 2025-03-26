/*
 * Filename: assembly.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-06
 *
 * Description:
 * This file contains functions for assembly linear system (CSEM or MT) in a PETGEM simulation. 
 *
 * List of Functions:
 * 
*/

/* C libraries */ 


/* PETSc libraries */
#include <petscsys.h>

/* PETGEM functions */ 
#include "inputs.h"
#include "grid.h"
#include "source.h"
#include "assembly.h"
#include "hvfem.h"
#include "constants.h"

// =============================================================================
// Function: setupSource
// =============================================================================
PetscErrorCode assemblySystem(DM dm, Mat A, Vec b, userParams *params, petgemGrid *grid, petgemSource *source, Vec *resistivity_x, Vec *resistivity_y, Vec *resistivity_z) {
    PetscFunctionBeginUser;
    
    /* Variables declaration */   
    PetscInt  m, n, numGaussPoints, numCoords;//, transitiveClosureCellSize;
    //PetscInt  *transitiveClosureCellPoints = NULL;
    PetscInt  cellOrientation[10];
    
    PetscReal   cellResistivity[NUM_DIMENSIONS];
    PetscReal   jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS];
    PetscReal   **gaussPoints, *weigths, **basisFunctions, **curlBasisFunctions, *XiEtaZeta;    
    PetscReal   **Me, **Ke;        
    PetscReal   omega;
    PetscReal   mu;

    const PetscScalar   *arrayCoords;
    PetscScalar constFactor;
    PetscScalar *localCellsResistivity_x, *localCellsResistivity_y, *localCellsResistivity_z;
    PetscScalar *cellCoords = NULL;
    PetscScalar *closureRHS;
    PetscScalar *closureLHS;

    PetscBool   isDG;
    
    PetscInt numDofIndices, *dofIndices;
    PetscSection section;
    
    /* Get linear system dimensions */ 
    PetscCall(MatGetSize(A, &m, &n));

    /* Print HEFEM statistics */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n HEFEM data:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Basis order             = %d\n", params->nord));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of dofs per vertex  = %d\n", grid->numDofInVertex));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of dofs per edge    = %d\n", grid->numDofInEdge));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of dofs per Face    = %d\n", grid->numDofInFace));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of dofs per volume  = %d\n", grid->numDofInVolume));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of dofs per cell    = %d\n", grid->numDofInCell));
    
    /* Print linear system statistics */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Assembly linear system:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of MPI tasks    = %d\n", params->numMPITasks));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Vector size         = %d\n", m));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Matrix size         = %d x %d\n", m, n));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Assembly status     = started\n"));

    /* Compute constants */
    omega = source->freq * 2.0 * PETSC_PI;
    mu    = 4.0 * PETSC_PI * 1.0e-7;
    constFactor = (0.0 + 1.0*PETSC_i) * omega * mu;

    /* Get the local values of the resistivity components */
    PetscCall(VecGetArray(*resistivity_x, &localCellsResistivity_x));
    PetscCall(VecGetArray(*resistivity_y, &localCellsResistivity_y));
    PetscCall(VecGetArray(*resistivity_z, &localCellsResistivity_z));

    /* Compute number of gauss points */
    PetscCall(computeNumGaussPoints3D(params->nord, &numGaussPoints));

    /* Allocate memory for gauss points */
    PetscCall(PetscCalloc1(numGaussPoints, &gaussPoints));
    for (PetscInt i = 0; i < numGaussPoints; i++) {
        PetscCall(PetscCalloc1(NUM_DIMENSIONS, &gaussPoints[i]));
    }
    PetscCall(PetscCalloc1(numGaussPoints, &weigths));

    /* Compute gauss points and its weigths */
    PetscCall(computeGaussPoints3D(numGaussPoints, gaussPoints, weigths));


    /* Allocate memory for RHS */
    PetscCall(PetscCalloc1(grid->numDofInCell, &closureRHS));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &XiEtaZeta));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &basisFunctions));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &curlBasisFunctions));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(grid->numDofInCell, &basisFunctions[i]));
        PetscCall(PetscCalloc1(grid->numDofInCell, &curlBasisFunctions[i]));
    }

    /* Get DM section */
    PetscCall(DMGetLocalSection(dm, &section));

    /* Perform finite element assembly for RHS */ 
    if (params->mode==0){ // CSEM mode
        /* Variables for CSEM source */ 
        PetscReal sourceRotationVector[NUM_DIMENSIONS], sourceVector[NUM_DIMENSIONS];
        PetscReal Dx[NUM_DIMENSIONS] = {0.0};
        PetscReal Dy[NUM_DIMENSIONS] = {0.1};
        PetscReal Dz[NUM_DIMENSIONS] = {0.2};

        /* Define dipole for total electric field formulation */
        Dx[0] = source->current * source->length;   // x-directed dipole
        Dy[1] = source->current * source->length;   // y-directed dipole
        Dz[2] = source->current * source->length;   // z-directed dipole

        /* Compute matrices for source rotation */
        PetscCall(vectorRotation(source->azimuth, source->dip, sourceRotationVector));

        /* Rotate source and setup electric field */
        sourceVector[0] = sourceRotationVector[0]*Dx[0] + sourceRotationVector[1]*Dy[0] + sourceRotationVector[2]*Dz[0];
        sourceVector[1] = sourceRotationVector[0]*Dx[1] + sourceRotationVector[1]*Dy[1] + sourceRotationVector[2]*Dz[1];
        sourceVector[2] = sourceRotationVector[0]*Dx[2] + sourceRotationVector[1]*Dy[2] + sourceRotationVector[2]*Dz[2];

        /* Locate source within computational domain */
        PetscInt sourceInCell = -1; 
        PetscCall(locatePoint(dm, source->position, &sourceInCell));

        /* Insert CSEM source */
        if (sourceInCell>=0){

            /* Get cell coordinates */ 
            PetscCall(DMPlexGetCellCoordinates(dm, sourceInCell, &isDG, &numCoords, &arrayCoords, &cellCoords));

            /* Compute jacobian and its inverse for sourceInCell */
            PetscCall(computeJacobian(cellCoords, jacobian, invJacobian));
            
            /* Transform xyz source position to XiEtaZeta coordinates (reference tetrahedral element) */
            PetscCall(tetrahedronXYZToXiEtaZeta(cellCoords, source->position, XiEtaZeta));

            /* Restore cell coordinates */ 
            PetscCall(DMPlexRestoreCellCoordinates(dm, sourceInCell, &isDG, &numCoords, &arrayCoords, &cellCoords));

            /* Get orientations for faces, edges and vertices of sourceInCell */
            /* Orden convention:
            - Faces indices start on position 2
            - Edges indices start on position 2 + NUM_FACES_PER_ELEMENT * 2
            - Vertices indices start on position 2 + NUM_FACES_PER_ELEMENT * 2 + NUM_EDGES_PER_ELEMENT * 2
            */
            /* Compute cell orientation */
            PetscCall(computeCellOrientation(dm, sourceInCell, cellOrientation));

            /* Compute basis functions for sourceInCell */ 
            PetscCall(computeBasisFunctions(params->nord, cellOrientation, jacobian, invJacobian, XiEtaZeta, basisFunctions, curlBasisFunctions));

            /* Get clousure indices for sourceInCell */
            PetscCall(DMPlexGetClosureIndices(dm, section, section, sourceInCell, PETSC_TRUE, &numDofIndices, &dofIndices, NULL, NULL));

            /* Compute contribution for clousure */
            for (PetscInt i = 0; i < grid->numDofInCell; i++){
                for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                    closureRHS[i] += (basisFunctions[j][i] * sourceVector[j]);
                }
                closureRHS[i] *= constFactor;
            }

            /* Add clousure to vector */
            PetscCall(VecSetValuesLocal(b, numDofIndices, dofIndices, closureRHS, INSERT_VALUES));

            /* Restore clousure indices for sourceInCell */
            PetscCall(DMPlexRestoreClosureIndices(dm, section, section, sourceInCell, PETSC_TRUE, &numDofIndices, &dofIndices, NULL, NULL));
        }
    } else if (params->mode==1){ // MT mode
        // TODO
    }
    else{
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Exiting: Source type not supported.\n"));
        PetscFunctionReturn(PETSC_ERR_SUP);
    }

    /* Perform global assembly for RHS */
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));
    
    /* Allocate memory for LHS */
    PetscCall(PetscCalloc1(grid->numDofInCell*grid->numDofInCell, &closureLHS)); 
    PetscCall(PetscCalloc1(grid->numDofInCell, &Me));
    PetscCall(PetscCalloc1(grid->numDofInCell, &Ke));
    for (PetscInt i = 0; i < grid->numDofInCell; i++){
        PetscCall(PetscCalloc1(grid->numDofInCell, &Me[i]));
        PetscCall(PetscCalloc1(grid->numDofInCell, &Ke[i]));
    }

    /* Perform finite element assembly for LHS */
    for (PetscInt i = grid->cellStart; i < grid->cellEnd; ++i) {
        /* Get spatial coordinates for cell i */ 
        PetscCall(DMPlexGetCellCoordinates(dm, i, &isDG, &numCoords, &arrayCoords, &cellCoords));

        /* Compute jacobian and its inverse for cell i */
        PetscCall(computeJacobian(cellCoords, jacobian, invJacobian));

        /* Restore coordinates for cell i */ 
        PetscCall(DMPlexRestoreCellCoordinates(dm, i, &isDG, &numCoords, &arrayCoords, &cellCoords));

        /* Get orientations for faces, edges and vertices of cell i */
        /* Orden convention:
        - Faces indices start on position 2
        - Edges indices start on position 2 + NUM_FACES_PER_ELEMENT * 2
        - Vertices indices start on position 2 + NUM_FACES_PER_ELEMENT * 2 + NUM_EDGES_PER_ELEMENT * 2
        */

        /* Compute cell orientation */
        for(PetscInt j = 0; j < 10; j++){
            cellOrientation[j] = 0;
        }
        PetscCall(computeCellOrientation(dm, i, cellOrientation));

        /* Get resistivity for cell i */ 
        cellResistivity[0] = PetscRealPart(localCellsResistivity_x[i]);
        cellResistivity[1] = PetscRealPart(localCellsResistivity_y[i]);
        cellResistivity[2] = PetscRealPart(localCellsResistivity_z[i]);

        /* Compute mass and stifness matrices for cell i */
        PetscCall(computeElementalMatrix(params->nord, cellOrientation, jacobian, invJacobian, numGaussPoints, gaussPoints, weigths, cellResistivity, Me, Ke));

        /* Get clousure indices for sourceInCell */
        PetscCall(DMPlexGetClosureIndices(dm, section, section, i, PETSC_TRUE, &numDofIndices, &dofIndices, NULL, NULL));

        /* Reset clousure */
        for(PetscInt j = 0; j < grid->numDofInCell; j++){
            for(PetscInt k = 0; k < grid->numDofInCell; k++){
                closureLHS[j * grid->numDofInCell + k] = 0.0 + 0.0*PETSC_i;
            }
        }    

        /* Compute elemental matrix for cell i */
        for(PetscInt j = 0; j < grid->numDofInCell; j++){
            for(PetscInt k = 0; k < grid->numDofInCell; k++){
                closureLHS[j * grid->numDofInCell + k] = Ke[j][k] - (constFactor * Me[j][k]);                
            }
        }    

        /* Add clousure to vector */
        PetscCall(MatSetValuesLocal(A, numDofIndices, dofIndices, numDofIndices, dofIndices, closureLHS, ADD_VALUES));

        /* Restore clousure indices for sourceInCell */
        PetscCall(DMPlexRestoreClosureIndices(dm, section, section, i, PETSC_TRUE, &numDofIndices, &dofIndices, NULL, NULL));
    }

    /* Perform global assembly for LHS */
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    /* End of assembly */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Assembly status     = completed\n"));

    /* Free memory */
    PetscCall(PetscFree(weigths));
    for (PetscInt i = 0; i < numGaussPoints; i++) {
        PetscCall(PetscFree(gaussPoints[i]));
    }
    PetscCall(PetscFree(gaussPoints));

    for (PetscInt i = 0; i < grid->numDofInCell; i++){
        PetscCall(PetscFree(Me[i]));
    PetscCall(PetscFree(Ke[i]));
    }
    PetscCall(PetscFree(Me));    
    PetscCall(PetscFree(Ke));
    PetscCall(PetscFree(closureLHS));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(basisFunctions[i]));
        PetscCall(PetscFree(curlBasisFunctions[i]));
    }
    PetscCall(PetscFree(basisFunctions));    
    PetscCall(PetscFree(curlBasisFunctions));
    PetscCall(PetscFree(XiEtaZeta));
    PetscCall(PetscFree(closureRHS));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}
    
