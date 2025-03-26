static char help[] = "PETGEM kernel for 3D CSEM modeling using high-order vector finite elements.\n\
  -nord                      <n>               : Basis order for finite element computations\n\
  -mesh_dm_plex_filename     <filename>        : Mesh file (Gmsh format)\n\
  -pc_type                   <type>            : Preconditioner type (lu)\n\
  -pc_factor_mat_solver_type <type>            : Solver type (mumps)\n";


/* C libraries */ 
#include <stdio.h>
#include <stdlib.h>

/* PETSc functions */   
#include <petscdmplex.h>

/* PETGEM functions */ 
#include "version.h"
#include "constants.h"
#include "common.h"
#include "inputs.h"  
#include "source.h" 
#include "grid.h" 
#include "assembly.h"
#include "solver.h"  
#include "postprocessing.h"    

#include "extrae_user_events.h"


int main(int argc, char **argv)
{

    /* Check if the --version option is provided */
    if (argc > 1 && strcmp(argv[1], "--version") == 0) {
        printf("PETGEM version %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
        return 0;
    }

    /* Variables declaration */
    PetscMPIInt     rank, size;
    DM              dm;
    Vec             resistivity_x, resistivity_y, resistivity_z, b, x;
    Mat             A;
    userParams      params;  
    petgemGrid      grid; 
    petgemSource    source;

    /* PETSC initialization */
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

    /* Print PETGEM header */
    PetscCall(printHeader());
   
    /* Parse user parameters */
    Extrae_event (1000, 1);
    PetscCall(readUserParams(&params, size));
    Extrae_event (1000, 0);
  
    /* Import mesh */    
Extrae_event (1000, 2);

    PetscCall(setupGrid(&dm, &grid, params.nord));

    Extrae_event (1000, 0);

    /* Create and setup source */
    Extrae_event (1000, 3);

    PetscCall(setupSource(&source, params.mode));
    Extrae_event (1000, 0);

    
    /* Import resistivity model */
    Extrae_event (1000, 4);
    PetscCall(readResistivityModel(dm, &resistivity_x, params.resistivityFile_x));
    PetscCall(readResistivityModel(dm, &resistivity_y, params.resistivityFile_y));
    PetscCall(readResistivityModel(dm, &resistivity_z, params.resistivityFile_z));
Extrae_event (1000, 0);
    /* Create and setup vectors and matrices for linear system */
Extrae_event (1000, 5);
    ISLocalToGlobalMapping mapping;
    PetscCall(DMGetLocalToGlobalMapping(dm, &mapping));
    PetscCall(DMCreateGlobalVector(dm, &b));
    PetscCall(DMCreateGlobalVector(dm, &x));    
    PetscCall(VecSetLocalToGlobalMapping(b, mapping));
    PetscCall(VecSetLocalToGlobalMapping(x, mapping));
    PetscCall(VecSetOption(b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));  
    PetscCall(VecSetFromOptions(b));
    PetscCall(VecSetFromOptions(x));
    PetscCall(DMSetAdjacency(dm, 0, PETSC_FALSE, PETSC_TRUE));
    PetscCall(DMCreateMatrix(dm, &A));
    PetscCall(MatSetLocalToGlobalMapping(A, mapping, mapping));
    PetscCall(MatSetFromOptions(A));
Extrae_event (1000, 0);
    /* Assembly linear system */
Extrae_event (1000, 6);
    PetscCall(assemblySystem(dm, A, b, &params, &grid, &source, &resistivity_x, &resistivity_y, &resistivity_z));
Extrae_event (1000, 0);
    /* Solve linear system */
Extrae_event (1000, 7);
    PetscCall(solveSystem(dm, A, b, x));
Extrae_event (1000, 0);
    /* Postprocessing solution */
Extrae_event (1000, 8);
    PetscCall(computeFields(dm, x, &params, &grid)); 
Extrae_event (1000, 0);
    /* Print PETGEM footer */
Extrae_event (1000, 9);
    PetscCall(printFooter());
Extrae_event (1000, 0);
    // /* Free memory */
    // Extrae_event (1000, 10);
    PetscCall(DMDestroy(&dm));
    PetscCall(VecDestroy(&resistivity_x));
    PetscCall(VecDestroy(&resistivity_y));
    PetscCall(VecDestroy(&resistivity_z));
    PetscCall(MatDestroy(&A));
    PetscCall(VecDestroy(&b));
    PetscCall(VecDestroy(&x));
    Extrae_event (1000, 0);
    
    /* PETSc finalize*/
    Extrae_event (1000, 11);
    PetscCall(PetscFinalize());
    Extrae_event (1000, 0);
    return 0;
}

