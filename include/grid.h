/*
 * Filename: grid.h
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-04
 *
 * Description:
 * This file contains a collection of definitions for grid functions that are used
 * throughout the PETGEM project. These functions are based on DMPLex provided by PETSc.
 *
 * List of definitions:
 * 
 *
 * Usage:
 * Include this file in your source code to utilize the grid functions. 
 * For example:
 * #include "grid.h"
 * 
*/

#ifndef GRID_H
#define GRID_H

/* PETGEM functions */ 
#include "inputs.h"

#include <petsc.h>

// =============================================================================
// Declaration of structures
// =============================================================================
typedef struct {
    PetscInt numCellsLocal;     // Number of local cells
    PetscInt numCellsGlobal;    // Number of global cells
    PetscInt numFacesLocal;     // Number of local faces
    PetscInt numFacesGlobal;    // Number of global cells
    PetscInt numEdgesLocal;     // Number of local edges
    PetscInt numEdgesGlobal;    // Number of global edges
    PetscInt numVerticesLocal;  // Number of local vertices
    PetscInt numVerticesGlobal; // Number of global vertices
    PetscInt numDofInVertex;    // Number of dofs per vertex
    PetscInt numDofInEdge;      // Number of dofs per edge
    PetscInt numDofInFace;      // Number of dofs per vertex
    PetscInt numDofInVolume;    // Number of dofs per volume
    PetscInt numDofInCell;      // Number of dofs per cell
    PetscInt cellStart;         // Index of global cell start
    PetscInt cellEnd;           // Index of global cell end
    PetscInt faceStart;         // Index of global face start
    PetscInt faceEnd;           // Index of global face end
    PetscInt edgeStart;         // Index of global edge start
    PetscInt edgeEnd;           // Index of global edge end
    PetscInt vertexStart;       // Index of global vertex start
    PetscInt vertexEnd;         // Index of global vertex end
    PetscInt dim;               // Number of dimensions
} petgemGrid;


// =============================================================================
// Declaration of functions
// =============================================================================
PetscErrorCode setupGrid(DM *dm, petgemGrid *grid, PetscInt nord);

PetscErrorCode readResistivityModel(DM dm, Vec *resistivityLocal, char *resistivityFile);

PetscErrorCode locatePoint(DM dm, PetscReal *position, PetscInt *pointInCell);

#endif
