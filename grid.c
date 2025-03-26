/*
 * Filename: grid.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-04
 *
 * Description:
 * This file contains a collection of functions for grid functions that are used
 * throughout the PETGEM project. These functions are based on DMPlex provided by PETSc.
 *
 * List of functions:
 * 
 * 
 * Usage:
 * Include this file in your source code to utilize the common functions. 
 * For example:
 * #include "grid.h"
 * 
*/

/* C libraries */ 

/* PETSc libraries */
#include <petscdmplex.h>

/* PETGEM funcions*/
#include "constants.h"
#include "grid.h"
#include "inputs.h"

// =============================================================================
// Function: setupGrid
// =============================================================================
PetscErrorCode setupGrid(DM *dm, petgemGrid *grid, PetscInt nord) {

	PetscFunctionBeginUser;

    /* Create DM object */
    PetscCall(DMCreate(PETSC_COMM_WORLD, dm));
    PetscCall(DMSetOptionsPrefix(*dm, "mesh_"));
    PetscCall(DMSetType(*dm, DMPLEX));
    PetscCall(DMSetNumFields(*dm, 1));
    PetscCall(DMSetFromOptions(*dm));
    PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));

    /*  Create label for dirichlet boundary conditions
        Values for boundaries:
        - boundaries = 100  */
    DMLabel     labelBoundary;
    PetscCall(DMCreateLabel(*dm, "Boundary"));
    PetscCall(DMGetLabel(*dm, "Boundary", &labelBoundary));
    PetscCall(DMPlexMarkBoundaryFaces(*dm, 100, labelBoundary));
    PetscCall(DMPlexLabelComplete(*dm, labelBoundary));

    /*  Create PetscSection. The PETSc convention in 3 dimensions is to number
        first cells, then vertices, then faces, and then edges.
        The above statement is not always true and we should not rely on that.
        It may be true for meshes read from GMSH files, but not for others.
        We are only guaranteed that points at different depths (or different heights, depending
        from where we start looking at the DAG) are numbered contiguously, this is why we
        can get start and end. Also, note that this numbering is purely local, because it is
        used to perform local mesh traversals   */ 
    
    PetscSection  section;
    PetscInt      numComp[] = {1};
    IS            boundaryIS;
    PetscInt      numBC = 1;
    PetscInt      bcField[1] = {0};
    PetscInt      numDofInVertex, numDofInEdge, numDofInFace, numDofInVolume, numDofInCell;
    
    /* Compute DOFs for PETGEM basis functions at vertex, edges, faces, and volume */
    numDofInVertex      = 0;
    numDofInEdge        = nord;
    numDofInFace        = nord*(nord-1);
    numDofInVolume      = nord*(nord-1)*(nord-2)/2;
    numDofInCell        = nord*(nord+2)*(nord+3)/2;
    PetscInt numDof[4]  = {numDofInVertex, numDofInEdge, numDofInFace, numDofInVolume};

    PetscCall(DMLabelGetStratumIS(labelBoundary, 100, &boundaryIS)); // Get the IS for boundaries
    PetscCall(DMPlexCreateSection(*dm, NULL, numComp, numDof, numBC, bcField, NULL, &boundaryIS, NULL, &section));
    PetscCall(DMSetLocalSection(*dm, section));

    /* Compute mesh statistics (number of vertices, edges, faces, elements) */
    PetscInt    numCellsLocal=0, numCellsGlobal=0, numFacesLocal=0, numFacesGlobal=0;
    PetscInt    numEdgesLocal=0, numEdgesGlobal=0, numVerticesLocal=0, numVerticesGlobal=0;
    PetscInt    pStart, cellStart, cellEnd, faceStart, faceEnd, edgeStart, edgeEnd, vertexStart, vertexEnd;
    IS          globalPointNumbering;
    const PetscInt *gidxs; 
  
    /* Create point numbering */
    PetscCall(DMPlexCreatePointNumbering(*dm, &globalPointNumbering));
    PetscCall(ISGetIndices(globalPointNumbering, &gidxs));
    // pStart is almost always 0, but we support nonzero too 
    PetscCall(DMPlexGetChart(*dm, &pStart, NULL)); 

    /* Get numbering for cells (depth 3) */
    PetscCall(DMPlexGetDepthStratum(*dm, 3, &cellStart, &cellEnd)); 
    
    /* Get number of local cells */
    numCellsLocal = cellEnd - cellStart;
    
    /* Get total num of cells by MPI reduction */
    PetscCall(MPI_Allreduce(&numCellsLocal, &numCellsGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    /* Get numbering for faces (depth 2) */
    PetscCall(DMPlexGetDepthStratum(*dm, 2, &faceStart, &faceEnd)); 
    
    /* Compute local number of faces */
    for (PetscInt f=faceStart; f<faceEnd; f++){
        /* This is the global index of face f, using the 
           convention that if it is negative it is not owned in parallel */    
        if  (gidxs[f - pStart] >= 0){
            numFacesLocal += 1; 
        }
    }    

    /* Get total num of faces by MPI reduction */
    PetscCall(MPI_Allreduce(&numFacesLocal, &numFacesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    /* Get numbering for edges (depth 1) */
    PetscCall(DMPlexGetDepthStratum(*dm, 1, &edgeStart, &edgeEnd)); 

    /* Compute local number of edges */
    for (PetscInt e=edgeStart; e<edgeEnd; e++){
        /* This is the global index of edge e, using the 
           convention that if it is negative it is not owned in parallel */    
        if  (gidxs[e - pStart] >= 0){
            numEdgesLocal += 1; 
        }
    }    

    /* Get total num of edges by MPI reduction */
    PetscCall(MPI_Allreduce(&numEdgesLocal, &numEdgesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    /* Get numbering for vertices (depth 0) */
    PetscCall(DMPlexGetDepthStratum(*dm, 0, &vertexStart, &vertexEnd)); 

    /* Compute local number of vertices */
    for (PetscInt v=vertexStart; v<vertexEnd; v++){
        /* This is the global index of vertex v, using the convention
           that if it is negative it is not owned in parallel */
        if  (gidxs[v - pStart] >= 0){
            numVerticesLocal += 1; 
        }
    }    

    /* Get total num of vertices by MPI reduction */
    PetscCall(MPI_Allreduce(&numVerticesLocal, &numVerticesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    /* Get number of dimensions */
    PetscInt    dim; 
    PetscCall(DMGetDimension(*dm, &dim));

    /* Setup petgemGrid */  
    grid->numCellsLocal     = numCellsLocal;     
    grid->numCellsGlobal    = numCellsGlobal;    
    grid->numFacesLocal     = numFacesLocal;     
    grid->numFacesGlobal    = numFacesGlobal;    
    grid->numEdgesLocal     = numEdgesLocal;     
    grid->numEdgesGlobal    = numEdgesGlobal;    
    grid->numVerticesLocal  = numVerticesLocal;  
    grid->numVerticesGlobal = numVerticesGlobal; 
    grid->numDofInVertex    = numDofInVertex;    
    grid->numDofInEdge      = numDofInEdge;      
    grid->numDofInFace      = numDofInFace;      
    grid->numDofInVolume    = numDofInVolume;
    grid->numDofInCell      = numDofInCell;    
    grid->cellStart         = cellStart; 
    grid->cellEnd           = cellEnd;           
    grid->faceStart         = faceStart;         
    grid->faceEnd           = faceEnd;           
    grid->edgeStart         = edgeStart;         
    grid->edgeEnd           = edgeEnd;           
    grid->vertexStart       = vertexStart;       
    grid->vertexEnd         = vertexEnd;         
    grid->dim               = dim;         

    /* Print petgemGrid data */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Mesh data:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of cells     = %d\n", grid->numCellsGlobal));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of faces     = %d\n", grid->numFacesGlobal));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of edges     = %d\n", grid->numEdgesGlobal));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Num of vertices  = %d\n", grid->numVerticesGlobal));

    /* Restore global numbering and free memory */
    PetscCall(ISRestoreIndices(globalPointNumbering, &gidxs));
    PetscCall(ISDestroy(&globalPointNumbering));
    PetscCall(PetscSectionDestroy(&section));
    PetscCall(ISDestroy(&boundaryIS));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


// =============================================================================
// Function: readResistivityModel
// =============================================================================
PetscErrorCode readResistivityModel(DM dm, Vec *resistivity, char *resistivityFile) {

    PetscFunctionBeginUser;

    PetscInt    dim;
    PetscInt    numCompResistivity[] = {1};
    PetscInt    numDofResistivity[4] = {0};
    DM  dmAux = NULL;
    PetscSection    resistivitySection;
    PetscViewer viewer;
    Vec resistivityGlobal, resistivityLocal;

    /* Create auxiliary DMs (they share the mesh, the field information is different) */
    PetscCall(DMGetDimension(dm, &dim));
    numDofResistivity[dim] = 1; // one dof at the cell center
    PetscCall(DMClone(dm, &dmAux));
    PetscCall(DMSetNumFields(dmAux, 1));
    PetscCall(DMPlexCreateSection(dmAux, NULL, numCompResistivity, numDofResistivity, 0, NULL, NULL, NULL, NULL, &resistivitySection));
    PetscCall(DMSetLocalSection(dmAux, resistivitySection));
    PetscCall(DMGetGlobalVector(dmAux, &resistivityGlobal));
    PetscCall(DMGetLocalVector(dmAux, &resistivityLocal));
        
    /* Read resistivity binary file */
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, resistivityFile, FILE_MODE_READ, &viewer));
    PetscCall(VecLoad(resistivityGlobal, viewer));
    PetscCall(DMGlobalToLocal(dmAux, resistivityGlobal, INSERT_VALUES, resistivityLocal));

    /* Copy resistivity to output VECs */
    PetscCall(VecDuplicate(resistivityLocal, resistivity));
    PetscCall(VecCopy(resistivityLocal, *resistivity));
   
    /* Return PETSc vecs */
    PetscCall(DMRestoreGlobalVector(dmAux, &resistivityGlobal)); 
    PetscCall(DMRestoreLocalVector(dmAux, &resistivityLocal)); 

    /* Free memory */
    PetscCall(DMDestroy(&dmAux));
    PetscCall(VecDestroy(&resistivityGlobal));
    PetscCall(VecDestroy(&resistivityLocal));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscSectionDestroy(&resistivitySection));
                
    PetscFunctionReturn(PETSC_SUCCESS);
}

// =============================================================================
// Function: locateCSEMSource
// =============================================================================
PetscErrorCode locatePoint(DM dm, PetscReal *position, PetscInt *pointInCell) {

    PetscFunctionBeginUser;

    /* Declarations */
    Vec         pointCoordinates;
    PetscSF     pointSF = NULL;
    PetscInt    numPointsFound;
    PetscScalar *inputPointCoordinates;
    const PetscSFNode   *pointCell;
    const PetscInt      *pointFound;

    /* Prepare PETSc vector with point coordinates */
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, NUM_DIMENSIONS, &pointCoordinates));
    PetscCall(VecSetBlockSize(pointCoordinates, NUM_DIMENSIONS));
    PetscCall(VecGetArrayWrite(pointCoordinates, &inputPointCoordinates));
    inputPointCoordinates[0] = position[0];
    inputPointCoordinates[1] = position[1];
    inputPointCoordinates[2] = position[2];
    PetscCall(VecRestoreArrayWrite(pointCoordinates, &inputPointCoordinates));
    
    /* Search point within computational domain */
    PetscCall(DMLocatePoints(dm, pointCoordinates, DM_POINTLOCATION_NONE, &pointSF));
    PetscCall(PetscSFGetGraph(pointSF, NULL, &numPointsFound, &pointFound, &pointCell));

    for (PetscInt i = 0; i < numPointsFound; i++) {
        *pointInCell = pointCell[i].index; 
    }
    
    /* Perform validation (at least one MPI task must found the point) */         
    PetscMPIInt rank;
    PetscBool pointFoundGlobal = PETSC_FALSE;
    PetscBool pointFoundLocal;
    
    /* Each process sets its local flag */
    pointFoundLocal = (*pointInCell < 0) ? PETSC_FALSE : PETSC_TRUE;

    /* Gather all local flags to the master process */
    MPI_Reduce(&pointFoundLocal, &pointFoundGlobal, 1, MPI_INT, MPI_LOR, 0, PETSC_COMM_WORLD);
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0) {
        if (pointFoundGlobal == PETSC_FALSE){
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Exiting: CSEM source position not located. Verify source parameters or improve mesh quality.\n"));
            PetscFunctionReturn(PETSC_ERR_ARG_WRONG);            
        }    
    }
    
    /* Free memory */
    PetscCall(VecDestroy(&pointCoordinates));
    PetscCall(PetscSFDestroy(&pointSF));

    PetscFunctionReturn(PETSC_SUCCESS);
}
